"""Targeted Tandem battery for a projected distance rung (d = 16, 18, ...).

The obvious candidate pool -- rows carrying `d_ub = D` -- is worthless:
of the 1,167 solved corpus rows that carried d_ub = 16, ZERO turned out
to have d_exact = 16 (the mode was 10).  The stored bound is tight only
24% of the time overall.

What works instead is a *demonstrated-rung* projection: take the (n, k)
cells that have already produced d = D-2 among settled rows, and keep
the rows there whose d_ub still admits D.  On the 2026-07-29 run that
turned 0/1167 into 31/244 at D = 16 (12.7%), concentrated in
[[168,4]] at 25%.

Timing calibration from that run -- solve cost is a sharp function of
the answer, and the classes do not overlap:

    d <= 14   <= 64 s          (median 19 s at d=14)
    d = 16    357 - 833 s      (median 686 s)

i.e. a code still running past ~2 min is already evidence of d >= 16.
The d=16 codes were 12.7% of the battery and 92% of its compute, so
size the window by the expected number of *hits*, not by pool size.

Usage (from experiments/bb_lab):
  uv run python scripts/distance_battery.py --binary <tandem> \\
      --pool 154,6,100 --pool 168,6,50 --pool 168,4,50 --min-dub 16 \\
      --shard 0 --nshards 10 --out d16_0.csv --work-dir w0
"""

from __future__ import annotations

import argparse
import csv
import glob
import subprocess
import time
from pathlib import Path

import duckdb

from bb_lab.checks import bb_check_matrices
from bb_lab.group import ZmZn
from bb_lab.linalg import nullspace_f2, quotient_complement_basis
from bb_lab.maxsat_distance import maxsat_distance
from bb_lab.poly import Poly

DB = "/Users/stavanjain/Code/qec-lab/experiments/bb_lab/data/bb_instances.duckdb"

SQL = """
SELECT instance_id, group_struct, ell, m, n, k, A_poly, B_poly, d_ub
FROM bb_instances
WHERE d_exact IS NULL AND n=? AND k=? AND (d_ub >= ? OR d_ub IS NULL)
ORDER BY instance_id
"""

SQL_N = """
SELECT instance_id, group_struct, ell, m, n, k, A_poly, B_poly, d_ub
FROM bb_instances
WHERE d_exact IS NULL AND n=? AND (d_ub >= ? OR d_ub IS NULL)
ORDER BY k DESC, instance_id
"""


def _pool(spec: str) -> tuple[int, int, int]:
    n, k, take = (int(x) for x in spec.split(","))
    return n, k, take


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--binary", required=True)
    ap.add_argument("--pool", action="append", default=[], metavar="n,k,take",
                    help="cell to draw from; repeatable")
    ap.add_argument("--whole-n", type=int, default=0,
                    help="instead of --pool, take every open row at this n "
                         "(k desc). Use for the n=288 targets.")
    ap.add_argument("--min-dub", type=int, default=16,
                    help="only rows whose d_ub still admits this distance")
    ap.add_argument("--shard", type=int, default=0)
    ap.add_argument("--nshards", type=int, default=1)
    ap.add_argument("--per-code-timeout", type=float, default=1800.0)
    ap.add_argument("--deadline", type=float, default=0.0,
                    help="stop starting new solves after this many seconds")
    ap.add_argument("--exclude-glob", default="",
                    help="comma-separated CSV globs of already-tested rows")
    ap.add_argument("--work-dir", default=".")
    ap.add_argument("--out", default="distance_battery_results.csv")
    args = ap.parse_args()

    done: set[str] = set()
    for g in filter(None, args.exclude_glob.split(",")):
        for f in glob.glob(g):
            for r in csv.DictReader(open(f)):
                if r.get("d"):
                    done.add(r["instance_id"])

    con = duckdb.connect(DB, read_only=True)
    rows: list[tuple] = []
    if args.whole_n:
        rows = [r for r in con.execute(
            SQL_N, [args.whole_n, args.min_dub]).fetchall()
            if r[0] not in done]
    else:
        for n, k, take in map(_pool, args.pool):
            cand = [r for r in con.execute(SQL, [n, k, args.min_dub]).fetchall()
                    if r[0] not in done]
            rows += cand[:take]
    con.close()
    mine = [r for i, r in enumerate(rows) if i % args.nshards == args.shard]

    out = Path(args.out)
    new = not out.exists()
    work = Path(args.work_dir)
    work.mkdir(parents=True, exist_ok=True)
    t_start = time.perf_counter()
    print(f"shard {args.shard}: {len(mine)} of {len(rows)} codes", flush=True)

    with out.open("a", newline="") as f:
        w = csv.writer(f)
        if new:
            w.writerow(["instance_id", "group", "n", "k", "d", "d_ub_known",
                        "seconds", "cost_step", "status", "A_poly", "B_poly"])
        for iid, gs, ell, m, n, k, ap_, bp_, dub in mine:
            if args.deadline and time.perf_counter() - t_start > args.deadline:
                print("deadline reached — not starting further solves",
                      flush=True)
                break
            G = ZmZn(int(ell), int(m))
            ch = bb_check_matrices(Poly.from_string(ap_, G),
                                   Poly.from_string(bp_, G))
            hx_even = not any(int(r.sum()) % 2 for r in ch.H_X)
            V = quotient_complement_basis(ch.H_X, nullspace_f2(ch.H_Z))
            step = hx_even and not any(int(v.sum()) % 2 for v in V)
            t0 = time.perf_counter()
            try:
                r = maxsat_distance(
                    ch, args.binary, mode="naive", work_dir=work,
                    timeout=args.per_code_timeout,
                    extra_args=("-cost-step=2",) if step else (),
                )
                d, secs, status = r.distance, r.solver_seconds, "ok"
            except subprocess.TimeoutExpired:
                d, secs, status = None, time.perf_counter() - t0, "timeout"
            except Exception as e:
                d, secs, status = None, time.perf_counter() - t0, \
                    f"error:{type(e).__name__}"
            w.writerow([iid, gs, n, k, d, dub, round(secs, 2), int(step),
                        status, ap_, bp_])
            f.flush()
            if d and d >= args.min_dub:
                print(f"  *** d={d} [[{n},{k},{d}]] {gs} A={ap_} B={bp_} "
                      f"({secs:.0f}s)", flush=True)


if __name__ == "__main__":
    main()
