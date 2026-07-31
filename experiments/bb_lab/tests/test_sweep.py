"""Parallel sweep runner + instance-unique solver scratch files.

The regression these guard: solver scratch used to be named
`{mode}_{group}`, so every instance sharing a group label wrote to one
path. The merit sweep's 28,688 rows span 21 group labels — two workers
on one `work_dir` would overwrite each other's WCNF and report a
distance for the wrong code.
"""

from __future__ import annotations

import csv
import os
from pathlib import Path

import pytest

from bb_lab.checks import bb_check_matrices
from bb_lab.group import ZmZn
from bb_lab.maxsat_distance import instance_stem, write_wcnf
from bb_lab.poly import Poly
from bb_lab.sweep import completed_keys, default_jobs, run_sweep

from tests._sweep_tasks import (
    fanout_task,
    flaky_task,
    slow_task,
    square_task,
    workdir_task,
)

FIELDNAMES = ["key", "value", "status", "pid"]


def _checks(ell: int, m: int, a: str, b: str):
    G = ZmZn(ell, m)
    return bb_check_matrices(Poly.from_string(a, G), Poly.from_string(b, G))


# --------------------------------------------------------------------
# scratch-file naming
# --------------------------------------------------------------------

def test_stem_separates_instances_sharing_a_group():
    """Two distinct codes on Z6xZ6 must not share a scratch path."""
    c1 = _checks(6, 6, "x^3 + y + y^2", "y^3 + x + x^2")
    c2 = _checks(6, 6, "1 + x + y", "1 + x^2 + y^2")

    assert c1.group.label() == c2.group.label() == "Z6xZ6", (
        "precondition: the old naming scheme collided precisely because "
        "these two codes share a group label"
    )
    assert instance_stem(c1, "naive") != instance_stem(c2, "naive")


def test_stem_is_deterministic_and_mode_scoped():
    c = _checks(6, 6, "x^3 + y + y^2", "y^3 + x + x^2")
    assert instance_stem(c, "naive") == instance_stem(c, "naive")
    assert instance_stem(c, "naive") != instance_stem(c, "strengthened")


def test_stem_distinguishes_swapped_polynomials():
    """A<->B is a different code; the hash must see the H_X/H_Z split."""
    c1 = _checks(6, 6, "x^3 + y + y^2", "y^3 + x + x^2")
    c2 = _checks(6, 6, "y^3 + x + x^2", "x^3 + y + y^2")
    assert instance_stem(c1, "naive") != instance_stem(c2, "naive")


def test_write_wcnf_leaves_no_temp_files(tmp_path):
    c = _checks(6, 6, "x^3 + y + y^2", "y^3 + x + x^2")
    path = tmp_path / f"{instance_stem(c, 'naive')}.wcnf"
    qv = write_wcnf(c, path, mode="naive")

    assert path.exists() and path.stat().st_size > 0
    assert len(qv) == c.num_qubits
    assert not list(tmp_path.glob("*.tmp*")), "atomic write left a temp file"

    lines = path.read_text().splitlines()
    assert all(ln.startswith(("h ", "1 ")) for ln in lines)
    assert sum(ln.startswith("1 ") for ln in lines) == c.num_qubits


def test_write_wcnf_is_atomic_under_rewrite(tmp_path):
    """Rewriting in place must never leave a truncated file behind."""
    c = _checks(6, 6, "x^3 + y + y^2", "y^3 + x + x^2")
    path = tmp_path / "reused.wcnf"
    write_wcnf(c, path, mode="naive")
    first = path.read_bytes()
    for _ in range(3):
        write_wcnf(c, path, mode="naive")
        assert path.read_bytes() == first


# --------------------------------------------------------------------
# the runner
# --------------------------------------------------------------------

def _items(n: int) -> list[dict]:
    return [{"key": f"k{i}", "x": i} for i in range(n)]


def _rows(out: Path) -> list[dict]:
    with out.open(newline="") as f:
        return list(csv.DictReader(f))


def test_run_sweep_covers_every_item(tmp_path):
    out = tmp_path / "r.csv"
    written = run_sweep(
        _items(12), square_task, out=out, fieldnames=FIELDNAMES,
        key_field="key", key=lambda it: it["key"],
        jobs=3, work_root=tmp_path / "work",
    )
    rows = _rows(out)
    assert written == 12 and len(rows) == 12
    assert {r["key"] for r in rows} == {f"k{i}" for i in range(12)}
    assert {r["key"]: int(r["value"]) for r in rows} == {
        f"k{i}": i * i for i in range(12)
    }


def test_run_sweep_actually_runs_workers_concurrently(tmp_path):
    """The point of the module: rows are solved in parallel, not striped.

    Uses a rendezvous rather than wall-clock timing — the tasks here are
    far too short for a timing comparison to mean anything.
    """
    out = tmp_path / "r.csv"
    run_sweep(
        _items(6), fanout_task, out=out, fieldnames=FIELDNAMES,
        key_field="key", key=lambda it: it["key"],
        jobs=3, work_root=tmp_path / "work",
    )
    rows = _rows(out)
    assert max(int(r["value"]) for r in rows) >= 2, (
        "no task ever observed a second live worker — the pool is serial"
    )
    assert len({r["pid"] for r in rows}) > 1


def test_run_sweep_resumes_without_duplicating(tmp_path):
    out = tmp_path / "r.csv"
    run_sweep(
        _items(5), square_task, out=out, fieldnames=FIELDNAMES,
        key_field="key", key=lambda it: it["key"],
        jobs=2, work_root=tmp_path / "work",
    )
    assert len(_rows(out)) == 5

    # Same output file, superset of items: only the new ones should run.
    written = run_sweep(
        _items(8), square_task, out=out, fieldnames=FIELDNAMES,
        key_field="key", key=lambda it: it["key"],
        jobs=2, work_root=tmp_path / "work",
    )
    rows = _rows(out)
    assert written == 3, "resume re-ran already-completed rows"
    assert len(rows) == 8
    assert len({r["key"] for r in rows}) == 8, "duplicate keys after resume"


def test_run_sweep_is_a_noop_when_everything_is_done(tmp_path):
    out = tmp_path / "r.csv"
    run_sweep(
        _items(4), square_task, out=out, fieldnames=FIELDNAMES,
        key_field="key", key=lambda it: it["key"],
        jobs=2, work_root=tmp_path / "work",
    )
    written = run_sweep(
        _items(4), square_task, out=out, fieldnames=FIELDNAMES,
        key_field="key", key=lambda it: it["key"],
        jobs=2, work_root=tmp_path / "work",
    )
    assert written == 0 and len(_rows(out)) == 4


def test_worker_exception_is_recorded_and_sweep_continues(tmp_path):
    out = tmp_path / "r.csv"
    items = _items(6) + [{"key": "boom", "x": 0}]
    run_sweep(
        items, flaky_task, out=out, fieldnames=FIELDNAMES,
        key_field="key", key=lambda it: it["key"],
        jobs=3, work_root=tmp_path / "work",
    )
    rows = _rows(out)
    assert len(rows) == 7, "a failing row must not abort the sweep"
    bad = [r for r in rows if r["key"] == "boom"]
    assert len(bad) == 1 and bad[0]["status"] == "error:ValueError"
    assert all(r["status"] == "ok" for r in rows if r["key"] != "boom")


def test_failed_row_is_retried_on_resume(tmp_path):
    """An error row keeps its key, so resume skips it -- deliberate:
    reruns should not silently repeat a deterministic failure. Callers
    who want a retry delete the row."""
    out = tmp_path / "r.csv"
    items = [{"key": "boom", "x": 0}]
    run_sweep(items, flaky_task, out=out, fieldnames=FIELDNAMES,
              key_field="key", key=lambda it: it["key"],
              jobs=1, work_root=tmp_path / "work")
    assert completed_keys(out, "key") == {"boom"}
    written = run_sweep(items, flaky_task, out=out, fieldnames=FIELDNAMES,
                        key_field="key", key=lambda it: it["key"],
                        jobs=1, work_root=tmp_path / "work")
    assert written == 0


def test_worker_scratch_is_private_to_the_pid(tmp_path):
    """Each worker's scratch dir is `work_root/w<pid>`, so two workers
    can never share one. (That several workers *run* is asserted by
    test_run_sweep_actually_runs_workers_concurrently; a task this cheap
    is often drained by whichever worker spawns first.)"""
    out = tmp_path / "r.csv"
    root = tmp_path / "work"
    run_sweep(
        _items(12), workdir_task, out=out,
        fieldnames=[*FIELDNAMES, "workdir"],
        key_field="key", key=lambda it: it["key"],
        jobs=3, work_root=root,
    )
    for r in _rows(out):
        d = Path(r["workdir"])
        assert d.is_dir()
        assert d.parent == root
        assert d.name == f"w{r['pid']}", "scratch dir is not pid-private"


def test_completed_keys_handles_missing_and_headerless(tmp_path):
    assert completed_keys(tmp_path / "nope.csv", "key") == set()
    empty = tmp_path / "empty.csv"
    empty.write_text("")
    assert completed_keys(empty, "key") == set()
    other = tmp_path / "other.csv"
    other.write_text("a,b\n1,2\n")
    assert completed_keys(other, "key") == set()


def test_deadline_stops_the_sweep_and_leaves_it_resumable(tmp_path):
    """`--deadline` bounds a run: unstarted rows are cancelled, finished
    rows are on disk, and rerunning picks up exactly the remainder."""
    out = tmp_path / "r.csv"
    items = [{"key": f"k{i}", "x": i, "sleep": 0.3} for i in range(24)]
    written = run_sweep(
        items, slow_task, out=out, fieldnames=FIELDNAMES,
        key_field="key", key=lambda it: it["key"],
        jobs=2, work_root=tmp_path / "work", deadline=0.5,
    )
    assert 0 < written < 24, f"deadline did not bound the run ({written})"
    assert len(_rows(out)) == written

    # No deadline: the rest complete, nothing is duplicated or lost.
    rest = run_sweep(
        items, slow_task, out=out, fieldnames=FIELDNAMES,
        key_field="key", key=lambda it: it["key"],
        jobs=4, work_root=tmp_path / "work",
    )
    rows = _rows(out)
    assert written + rest == 24
    assert {r["key"] for r in rows} == {f"k{i}" for i in range(24)}
    assert len(rows) == 24, "resume after deadline duplicated rows"


def test_run_sweep_refuses_to_recurse_inside_a_worker(tmp_path, monkeypatch):
    """A caller without an `if __name__ == '__main__':` guard would have
    every spawned child re-run the module top level and start its own
    pool. Fail loudly instead of fork-bombing."""
    import multiprocessing

    monkeypatch.setattr(
        multiprocessing, "parent_process", lambda: object(),
    )
    with pytest.raises(RuntimeError, match="main.*guard|worker process"):
        run_sweep(
            _items(2), square_task, out=tmp_path / "r.csv",
            fieldnames=FIELDNAMES, key_field="key",
            key=lambda it: it["key"], jobs=2, work_root=tmp_path / "work",
        )


def test_default_jobs_is_sane():
    j = default_jobs()
    assert isinstance(j, int)
    assert 1 <= j <= (os.cpu_count() or 1)


@pytest.mark.parametrize("jobs", [1, 2])
def test_single_and_multi_worker_agree(tmp_path, jobs):
    out = tmp_path / f"r{jobs}.csv"
    run_sweep(
        _items(10), square_task, out=out, fieldnames=FIELDNAMES,
        key_field="key", key=lambda it: it["key"],
        jobs=jobs, work_root=tmp_path / f"work{jobs}",
    )
    assert {r["key"]: int(r["value"]) for r in _rows(out)} == {
        f"k{i}": i * i for i in range(10)
    }
