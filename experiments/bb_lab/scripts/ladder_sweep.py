"""Intermediate-ladder sweep: MaxCDCL-naive across corpus rungs.

Bridges the gross (n=144, ~3 s) → bb_288 (open) hardness gap with
measured points: calibration rungs (d_exact known — validates the
solver and fits the scaling curve) and frontier rungs (d_ub-only,
unscreened — every settled row is a new exact distance).

Usage (from experiments/bb_lab):
  uv run python scripts/ladder_sweep.py --mode calibration
  uv run python scripts/ladder_sweep.py --mode frontier --per-code-timeout 900
"""

from __future__ import annotations

import argparse
import csv
import subprocess
import time
from pathlib import Path

import duckdb

LAB_ROOT = Path(__file__).resolve().parent.parent
DB = "/Users/stavanjain/Code/qec-lab/experiments/bb_lab/data/bb_instances.duckdb"

from bb_lab.checks import bb_check_matrices  # noqa: E402
from bb_lab.group import ZmZn  # noqa: E402
from bb_lab.maxsat_distance import maxsat_distance  # noqa: E402
from bb_lab.poly import Poly  # noqa: E402

CALIBRATION_SQL = """
(SELECT * FROM bb_instances WHERE n=150 AND k=8 AND d_exact=12 LIMIT 6)
UNION ALL
(SELECT * FROM bb_instances WHERE n=150 AND k=8 AND d_exact=10 LIMIT 3)
UNION ALL
(SELECT * FROM bb_instances WHERE n=180 AND d_exact=10 ORDER BY k DESC LIMIT 4)
UNION ALL
(SELECT * FROM bb_instances WHERE n=210 AND d_exact IS NOT NULL
 ORDER BY d_exact DESC LIMIT 3)
"""

FRONTIER_SQL = """
SELECT * FROM bb_instances
WHERE n BETWEEN 150 AND 280 AND d_exact IS NULL
  AND d_ub BETWEEN {min_dub} AND {max_dub}
ORDER BY n DESC, d_ub ASC   -- larger n first; modest d_ub = likelier true
LIMIT {limit}
"""


def _already_swept(out: Path) -> set[str]:
    if not out.exists():
        return set()
    with out.open() as f:
        return {row.split(",")[0] for row in f.read().splitlines()[1:]}


def run(args) -> None:
    con = duckdb.connect(DB, read_only=True)
    sql = (
        CALIBRATION_SQL if args.mode == "calibration"
        else FRONTIER_SQL.format(
            limit=args.limit * 3,  # headroom for already-swept exclusion
            min_dub=args.min_dub, max_dub=args.max_dub,
        )
    )
    cols = [d[0] for d in con.execute("SELECT * FROM bb_instances LIMIT 0").description]
    rows = [dict(zip(cols, r)) for r in con.execute(sql).fetchall()]
    out = Path(args.out)
    done = _already_swept(out)
    rows = [r for r in rows if r["instance_id"] not in done][: args.limit]
    print(f"{args.mode}: {len(rows)} rungs", flush=True)
    new = not out.exists()
    with out.open("a", newline="") as f:
        w = csv.writer(f)
        if new:
            w.writerow([
                "instance_id", "n", "k", "d_exact_known", "d_ub_known",
                "d_found", "solver_seconds", "status",
            ])
        for row in rows:
            G = ZmZn(int(row["ell"]), int(row["m"]))
            checks = bb_check_matrices(
                Poly.from_string(row["A_poly"], G),
                Poly.from_string(row["B_poly"], G),
            )
            t0 = time.perf_counter()
            try:
                r = maxsat_distance(
                    checks, args.binary, mode="naive",
                    work_dir=args.work_dir, timeout=args.per_code_timeout,
                )
                status = "ok"
                d_found, secs = r.distance, r.solver_seconds
                if row["d_exact"] is not None and d_found != row["d_exact"]:
                    status = "MISMATCH"
            except subprocess.TimeoutExpired:
                status, d_found, secs = "timeout", None, time.perf_counter() - t0
            except Exception as e:  # keep sweeping; record the failure
                status, d_found, secs = f"error:{type(e).__name__}", None, 0.0
            w.writerow([
                row["instance_id"], row["n"], row["k"],
                row["d_exact"], row["d_ub"], d_found, round(secs, 2), status,
            ])
            f.flush()
            print(
                f"  [{row['instance_id']}] n={row['n']} k={row['k']} "
                f"d_known={row['d_exact']}/{row['d_ub']} → {d_found} "
                f"({secs:.1f}s, {status})",
                flush=True,
            )


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--mode", choices=["calibration", "frontier"],
                    default="calibration")
    ap.add_argument("--limit", type=int, default=20)
    ap.add_argument("--min-dub", type=int, default=12)
    ap.add_argument("--max-dub", type=int, default=99)
    ap.add_argument("--per-code-timeout", type=float, default=600.0)
    ap.add_argument("--binary", default=None, required=True)
    ap.add_argument("--work-dir", default=".")
    ap.add_argument("--out", default="ladder_results.csv")
    run(ap.parse_args())


if __name__ == "__main__":
    main()
