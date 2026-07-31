"""Resumable, parallel runner for per-instance solver sweeps.

Every corpus battery in `scripts/` has the same shape: pull rows from
duckdb, solve each, append a CSV row, be resumable across interruptions.
They ran strictly serially — the committed artifacts are ~70 CPU-hours
of sequential solving — and the one script that did parallelise
(`merit_sweep.py`) striped rows `i % nshards == shard` across manually
launched processes.

Static striping is the wrong shape for this workload. Solve times are
tail-dominated: the merit sweep's slowest 1% of rows account for 35% of
its total time, and the d=18 battery's slowest single row is 6,924 s
against a 7 s median. A stripe that happens to collect two hard rows
runs long after its siblings have gone idle. This module hands out work
dynamically instead — one queue, N workers, next row to whoever is free
— so the wall clock is bounded by the single hardest row rather than by
the unluckiest stripe.

Sizing: `default_jobs()` returns the *performance*-core count on Apple
Silicon rather than `os.cpu_count()`. These solves run roughly 2-3x
slower on an efficiency core, so a worker that lands on one becomes the
straggler that defines the wall clock; the E-cores are better left to
the OS. On this class of machine that is 4 workers, not 10.

Usage:

    from bb_lab.sweep import bb_distance_task, run_sweep

    n = run_sweep(
        items, bb_distance_task,
        out=Path("results.csv"), fieldnames=[...], key_field="instance_id",
        key=lambda it: it["instance_id"], jobs=6,
    )
"""

from __future__ import annotations

import csv
import multiprocessing
import os
import subprocess
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Any, Callable, Sequence

_WORK_DIR: Path | None = None


def default_jobs() -> int:
    """Worker count: performance cores on Apple Silicon, else cpu_count-1.

    `hw.perflevel0` is the performance cluster on Apple Silicon; the call
    fails on Intel Macs and Linux, where we fall back to leaving one core
    for the parent and the OS.
    """
    try:
        proc = subprocess.run(
            ["sysctl", "-n", "hw.perflevel0.logicalcpu"],
            capture_output=True, text=True, timeout=5,
        )
        if proc.returncode == 0 and proc.stdout.strip().isdigit():
            return max(1, int(proc.stdout.strip()))
    except (OSError, subprocess.SubprocessError):
        pass
    return max(1, (os.cpu_count() or 2) - 1)


def _as_tuple(x) -> tuple[str, ...]:
    return (x,) if isinstance(x, str) else tuple(x)


def completed_keys(
    out: Path, key_field: str | Sequence[str]
) -> set[tuple[str, ...]]:
    """Keys already present in `out`, for resume. Empty if absent.

    `key_field` may name several columns, for sweeps whose rows have no
    single id — `cell_hunt` draws fresh codes and identifies them by the
    canonical (A_poly, B_poly) pair already in its CSV, so resume costs
    it no schema change.
    """
    fields = _as_tuple(key_field)
    if not out.exists():
        return set()
    with out.open(newline="") as f:
        reader = csv.DictReader(f)
        if reader.fieldnames is None or any(
            fl not in reader.fieldnames for fl in fields
        ):
            return set()
        return {
            tuple(row[fl] for fl in fields)
            for row in reader
            if all(row.get(fl) for fl in fields)
        }


def worker_work_dir() -> Path:
    """Scratch directory owned by the calling worker.

    Instance-unique WCNF names already make a shared directory safe; a
    directory per worker additionally bounds how much scratch any one
    process can leave behind after a crash.
    """
    return _WORK_DIR if _WORK_DIR is not None else Path(".")


def _init_worker(work_root: str) -> None:
    global _WORK_DIR
    _WORK_DIR = Path(work_root) / f"w{os.getpid()}"
    _WORK_DIR.mkdir(parents=True, exist_ok=True)


def run_sweep(
    items: Sequence[Any],
    task: Callable[[Any], dict],
    *,
    out: Path,
    fieldnames: Sequence[str],
    key_field: str | Sequence[str],
    key: Callable[[Any], Any],
    jobs: int | None = None,
    work_root: Path | str = "sweep_work",
    report: Callable[[dict], str | None] | None = None,
    deadline: float | None = None,
) -> int:
    """Run `task` over `items` on a dynamic queue, appending rows to `out`.

    `task` must be a module-level function (it is pickled to the workers)
    taking one item and returning a dict over `fieldnames`. It is
    responsible for catching its own solver errors and encoding them in
    the row — an exception that escapes is recorded as an error row and
    the sweep continues.

    `key_field` names the CSV column(s) identifying a row, and `key(item)`
    must return a matching string (or tuple of strings, for several
    columns). Resume skips items whose key is already in `out`.

    Rows land in completion order, not submission order; resume is by
    `key_field`, so ordering carries no meaning. Returns the number of
    rows written.

    `deadline` (seconds) bounds the run: once it passes, queued rows that
    have not started are cancelled. Rows already in flight are allowed to
    finish, so the call can overrun by up to one solve — matching the
    serial scripts' "stop starting new solves" rather than killing work
    that is nearly done.
    """
    if multiprocessing.parent_process() is not None:
        raise RuntimeError(
            "run_sweep() was called from inside a worker process. The "
            "calling script is missing its `if __name__ == '__main__':` "
            "guard, so the `spawn` start method re-ran the module top "
            "level in every child — each of which would spawn its own "
            "pool. Wrap the call in a main guard."
        )
    out = Path(out)
    key_cols = _as_tuple(key_field)
    done = completed_keys(out, key_field)
    pending = [it for it in items if _as_tuple(key(it)) not in done]
    jobs = jobs or default_jobs()
    work_root = Path(work_root)
    work_root.mkdir(parents=True, exist_ok=True)

    print(
        f"sweep: {len(pending)} rows to run "
        f"({len(done)} already in {out.name}), {jobs} workers",
        flush=True,
    )
    if not pending:
        return 0

    write_header = not out.exists() or out.stat().st_size == 0
    written = 0
    timed_out = 0
    t0 = time.perf_counter()

    with out.open("a", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(fieldnames))
        if write_header:
            writer.writeheader()
            f.flush()

        ex = ProcessPoolExecutor(
            max_workers=jobs,
            initializer=_init_worker,
            initargs=(str(work_root),),
        )
        try:
            futures = {ex.submit(task, it): it for it in pending}
            stopped = False
            for fut in as_completed(futures):
                item = futures[fut]
                if fut.cancelled():
                    continue
                try:
                    row = fut.result()
                except Exception as e:  # worker died or task re-raised
                    row = {fn: None for fn in fieldnames}
                    row.update(zip(key_cols, _as_tuple(key(item))))
                    if "status" in row:
                        row["status"] = f"error:{type(e).__name__}"
                writer.writerow({fn: row.get(fn) for fn in fieldnames})
                f.flush()
                written += 1
                if row.get("status") == "timeout":
                    timed_out += 1
                line = report(row) if report else None
                if line:
                    print(line, flush=True)
                if written % 200 == 0:
                    rate = written / (time.perf_counter() - t0)
                    print(
                        f"  … {written}/{len(pending)} rows "
                        f"({rate:.1f}/s)", flush=True,
                    )
                if (deadline and not stopped
                        and time.perf_counter() - t0 > deadline):
                    stopped = True
                    dropped = sum(fu.cancel() for fu in futures)
                    print(
                        f"deadline reached — cancelled {dropped} rows that "
                        f"had not started; rerun to resume", flush=True,
                    )
        except KeyboardInterrupt:
            print("\ninterrupted — finished rows are already on disk; "
                  "rerun the same command to resume", flush=True)
            ex.shutdown(wait=False, cancel_futures=True)
            raise
        finally:
            ex.shutdown(wait=True)

    dt = time.perf_counter() - t0
    print(f"sweep: {written} rows in {dt / 60:.1f} min", flush=True)
    if timed_out and jobs > 1:
        # Measured on this workload: 4 concurrent solves inflate each
        # solve's own runtime ~20% (shared last-level cache and memory
        # bandwidth, not scheduling). A row that took 530 s serially hit
        # a 600 s cap under 4 workers — so parallelising can silently
        # *lose* exactly the slow, interesting rows unless the per-row
        # timeout is raised to match.
        print(
            f"WARNING: {timed_out} row(s) timed out at jobs={jobs}. "
            f"Concurrent solves run ~20% slower each; a timeout tuned on "
            f"a serial run will drop rows here. Raise the per-row timeout "
            f"or lower --jobs, then rerun to resume.", flush=True,
        )
    return written


def bb_distance_task(payload: dict) -> dict:
    """Solve one BB code's exact distance. The task every battery shares.

    `payload` carries the instance (`ell`, `m`, `A_poly`, `B_poly`), the
    run config (`binary`, `mode`, `timeout`), and a `passthrough` dict of
    columns copied verbatim into the result row.

    The `-cost-step=2` premise (every H_X row even AND every class even,
    so weight parity is a coset invariant) is verified here, per code,
    exactly as the serial scripts did — the flag is a soundness
    obligation on the caller, so it stays next to the solve.
    """
    # Imported inside the worker: these pull in numpy/pysat, and under
    # spawn the parent's copies are not inherited anyway.
    from .checks import bb_check_matrices
    from .group import ZmZn
    from .linalg import nullspace_f2, quotient_complement_basis
    from .maxsat_distance import maxsat_distance
    from .poly import Poly

    row = dict(payload.get("passthrough", {}))
    G = ZmZn(int(payload["ell"]), int(payload["m"]))
    checks = bb_check_matrices(
        Poly.from_string(payload["A_poly"], G),
        Poly.from_string(payload["B_poly"], G),
    )
    hx_even = not any(int(r.sum()) % 2 for r in checks.H_X)
    V = quotient_complement_basis(checks.H_X, nullspace_f2(checks.H_Z))
    step = hx_even and not any(int(v.sum()) % 2 for v in V)

    t0 = time.perf_counter()
    try:
        res = maxsat_distance(
            checks, payload["binary"],
            mode=payload.get("mode", "naive"),
            work_dir=worker_work_dir(),
            timeout=payload.get("timeout"),
            extra_args=("-cost-step=2",) if step else (),
        )
        d, secs, status = res.distance, res.solver_seconds, "ok"
    except subprocess.TimeoutExpired:
        d, secs, status = None, time.perf_counter() - t0, "timeout"
    except Exception as e:
        d, secs, status = None, time.perf_counter() - t0, \
            f"error:{type(e).__name__}"

    row.update(
        d=d, seconds=round(secs, 2), cost_step=int(step), status=status,
    )
    return row


def duckdb_rows(db: str, sql: str) -> list[dict]:
    """Fetch `sql` as dicts and close the connection.

    Sweeps must not hold a duckdb handle across `fork`/`spawn`; reading
    everything up front keeps the workers free of database state.
    """
    import duckdb

    con = duckdb.connect(db, read_only=True)
    try:
        cur = con.execute(sql)
        cols = [c[0] for c in cur.description]
        return [dict(zip(cols, r)) for r in cur.fetchall()]
    finally:
        con.close()
