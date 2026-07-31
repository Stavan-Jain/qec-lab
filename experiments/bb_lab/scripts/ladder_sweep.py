"""Intermediate-ladder sweep: MaxCDCL-naive across corpus rungs.

Bridges the gross (n=144, ~3 s) → bb_288 (open) hardness gap with
measured points: calibration rungs (d_exact known — validates the
solver and fits the scaling curve) and frontier rungs (d_ub-only,
unscreened — every settled row is a new exact distance).

Rungs run on `bb_lab.sweep`'s dynamic queue; `--jobs` defaults to the
machine's performance-core count.

Usage (from experiments/bb_lab):
  uv run python scripts/ladder_sweep.py --mode calibration --binary <tandem>
  uv run python scripts/ladder_sweep.py --mode frontier --per-code-timeout 1800
"""

from __future__ import annotations

import argparse
from pathlib import Path

LAB_ROOT = Path(__file__).resolve().parent.parent
DB = "/Users/stavanjain/Code/qec-lab/experiments/bb_lab/data/bb_instances.duckdb"

from bb_lab.sweep import (  # noqa: E402
    bb_distance_task, completed_keys, default_jobs, duckdb_rows, run_sweep,
)

FIELDNAMES = [
    "instance_id", "n", "k", "d_exact_known", "d_ub_known",
    "d_found", "solver_seconds", "status",
]


def ladder_task(payload: dict) -> dict:
    """Solve, then flag any disagreement with a known exact distance.

    Column names differ from the shared task's (`d_found`,
    `solver_seconds`), so remap rather than change the CSV schema and
    strand the committed ladder results.
    """
    row = bb_distance_task(payload)
    row["d_found"] = row.pop("d")
    row["solver_seconds"] = row.pop("seconds")
    known = row.get("d_exact_known")
    if (row["status"] == "ok" and known is not None
            and row["d_found"] != known):
        row["status"] = "MISMATCH"
    return row

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
WHERE n BETWEEN {n_min} AND {n_max} AND d_exact IS NULL
  AND d_ub BETWEEN {min_dub} AND {max_dub}
ORDER BY n DESC, d_ub ASC   -- larger n first; modest d_ub = likelier true
LIMIT {limit}
"""


def run(args) -> None:
    sql = (
        CALIBRATION_SQL if args.mode == "calibration"
        else FRONTIER_SQL.format(
            limit=args.limit * 3,  # headroom for already-swept exclusion
            min_dub=args.min_dub, max_dub=args.max_dub,
            n_min=args.n_min, n_max=args.n_max,
        )
    )
    rows = duckdb_rows(DB, sql)
    out = Path(args.out)
    # Applied before --limit as well as inside run_sweep, so a resumed
    # run fills the limit with fresh rungs instead of re-counting done ones.
    done = {k[0] for k in completed_keys(out, "instance_id")}
    rows = [r for r in rows if r["instance_id"] not in done][: args.limit]
    print(f"{args.mode}: {len(rows)} rungs, "
          f"jobs={args.jobs or default_jobs()}", flush=True)

    items = [
        {
            "ell": r["ell"], "m": r["m"],
            "A_poly": r["A_poly"], "B_poly": r["B_poly"],
            "binary": args.binary, "mode": "naive",
            "timeout": args.per_code_timeout,
            "passthrough": {
                "instance_id": r["instance_id"], "n": r["n"], "k": r["k"],
                "d_exact_known": r["d_exact"], "d_ub_known": r["d_ub"],
            },
        }
        for r in rows
    ]

    def report(row: dict) -> str | None:
        tag = ""
        if row["status"] == "MISMATCH":
            tag = "  <<<<<< MISMATCH"
        elif row["status"] not in ("ok",):
            tag = f"  ({row['status']})"
        return (
            f"  [{row['instance_id']}] n={row['n']} k={row['k']} "
            f"d_known={row['d_exact_known']}/{row['d_ub_known']} → "
            f"{row['d_found']} ({row['solver_seconds']:.1f}s){tag}"
        )

    run_sweep(
        items, ladder_task,
        out=out, fieldnames=FIELDNAMES, key_field="instance_id",
        key=lambda it: it["passthrough"]["instance_id"],
        jobs=args.jobs, work_root=Path(args.work_dir), report=report,
    )


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--mode", choices=["calibration", "frontier"],
                    default="calibration")
    ap.add_argument("--limit", type=int, default=20)
    ap.add_argument("--min-dub", type=int, default=12)
    ap.add_argument("--max-dub", type=int, default=99)
    ap.add_argument("--n-min", type=int, default=150)
    ap.add_argument("--n-max", type=int, default=280)
    ap.add_argument("--tandem-step", action="store_true",
                    help="NO-OP, kept so existing invocations still run: "
                         "the shared task now verifies the coset-parity "
                         "premise per code and passes -cost-step=2 "
                         "whenever it holds, as the other batteries do")
    ap.add_argument("--jobs", type=int, default=None,
                    help="workers (default: performance-core count)")
    # Raised from 600 s alongside the move to a queue: concurrent solves
    # run ~20% slower each, and a cap tuned serially drops rows that
    # would otherwise settle.
    ap.add_argument("--per-code-timeout", type=float, default=1800.0)
    ap.add_argument("--binary", default=None, required=True)
    ap.add_argument("--work-dir", default="ladder_work")
    ap.add_argument("--out", default="ladder_results.csv")
    run(ap.parse_args())


if __name__ == "__main__":
    main()
