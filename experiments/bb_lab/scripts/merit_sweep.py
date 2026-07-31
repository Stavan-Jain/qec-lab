"""Merit sweep: settle every corpus row that could rival the gross code.

The comparison quantity is the BB figure of merit

    q = k * d^2 / n            gross [[144,12,12]] has q = 12.

Every corpus row has weight-6 checks (|A| = |B| = 3 on all 58,021
rows), so (n, k, d) are the only free parameters.  Under a distance
cap d <= D, a row can only reach q >= QMIN if k * D^2 / n >= QMIN --
rows failing that are excluded by arithmetic and never queried, which
is what makes a whole-corpus sweep affordable.

Per code: naive WCNF + Tandem, with `-cost-step=2` passed only when
the coset weight-parity premise is verified to hold for that code
(same gate as ladder_sweep.py).

Measured on the 2026-07-29 run: 28,003 rows at qmin=8, D=14, zero
timeouts, ~0.32 s/row mean over 6 manually launched shards.  Those
shards striped rows `i % nshards`, which handles this workload badly --
its slowest 1% of rows are 35% of its total time, so a stripe that
draws two hard rows runs long after the others idle.  `bb_lab.sweep`
now hands rows out dynamically instead; `--jobs` defaults to the
machine's performance-core count.

Usage (from experiments/bb_lab):
  uv run python scripts/merit_sweep.py --binary <tandem> --qmin 8

`--shard/--nshards` still partition the row set, for splitting a sweep
across several machines; within one machine leave them alone and let
`--jobs` do the work.
"""

from __future__ import annotations

import argparse
from pathlib import Path

from bb_lab.sweep import bb_distance_task, default_jobs, duckdb_rows, run_sweep

DB = "/Users/stavanjain/Code/qec-lab/experiments/bb_lab/data/bb_instances.duckdb"

SQL = """
SELECT instance_id, group_struct, ell, m, n, k, A_poly, B_poly, d_ub
FROM bb_instances
WHERE {settled} k > 0 AND k * {dcap} * {dcap} / n >= {qmin}
ORDER BY k * {dcap} * {dcap} / n DESC, n ASC, instance_id
"""

FIELDNAMES = [
    "instance_id", "group", "n", "k", "d", "q", "d_ub_known",
    "seconds", "cost_step", "status", "A_poly", "B_poly",
]


def merit_task(payload: dict) -> dict:
    """Distance solve plus the figure of merit it feeds."""
    row = bb_distance_task(payload)
    d = row.get("d")
    row["q"] = round(row["k"] * d * d / row["n"], 3) if d else None
    return row


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--binary", required=True)
    ap.add_argument("--jobs", type=int, default=None,
                    help="workers (default: performance-core count)")
    ap.add_argument("--shard", type=int, default=0)
    ap.add_argument("--nshards", type=int, default=1)
    ap.add_argument("--qmin", type=float, default=8.0,
                    help="keep rows whose best-case q at d=DCAP is >= QMIN")
    ap.add_argument("--dcap", type=float, default=14.0,
                    help="distance cap defining 'best case' in the screen")
    ap.add_argument("--per-code-timeout", type=float, default=600.0)
    ap.add_argument("--limit", type=int, default=None,
                    help="cap rows after sharding (smoke tests, benchmarks)")
    ap.add_argument("--include-settled", action="store_true",
                    help="also re-run rows that already have d_exact; the "
                         "screen is exhausted, so this is how you get a "
                         "row set to benchmark against")
    ap.add_argument("--work-dir", default="merit_work")
    ap.add_argument("--out", default="merit_sweep_results.csv")
    args = ap.parse_args()

    rows = duckdb_rows(DB, SQL.format(
        qmin=args.qmin, dcap=args.dcap,
        settled="" if args.include_settled else "d_exact IS NULL AND",
    ))
    rows = [r for i, r in enumerate(rows) if i % args.nshards == args.shard]
    if args.limit is not None:
        rows = rows[: args.limit]

    items = [
        {
            "ell": r["ell"], "m": r["m"],
            "A_poly": r["A_poly"], "B_poly": r["B_poly"],
            "binary": args.binary, "mode": "naive",
            "timeout": args.per_code_timeout,
            "passthrough": {
                "instance_id": r["instance_id"], "group": r["group_struct"],
                "n": r["n"], "k": r["k"], "d_ub_known": r["d_ub"],
                "A_poly": r["A_poly"], "B_poly": r["B_poly"],
            },
        }
        for r in rows
    ]

    def report(row: dict) -> str | None:
        if row.get("d") and row.get("q") and row["q"] >= args.qmin:
            return (f"  ** [[{row['n']},{row['k']},{row['d']}]] "
                    f"q={row['q']}  {row['group']}  "
                    f"A={row['A_poly']} B={row['B_poly']}")
        if row.get("status") not in (None, "ok"):
            return f"  !! [{row['instance_id'][:16]}] {row['status']}"
        return None

    print(f"merit sweep: shard {args.shard}/{args.nshards}, "
          f"{len(items)} rows, jobs={args.jobs or default_jobs()}", flush=True)
    run_sweep(
        items, merit_task,
        out=Path(args.out), fieldnames=FIELDNAMES,
        key_field="instance_id",
        key=lambda it: it["passthrough"]["instance_id"],
        jobs=args.jobs, work_root=Path(args.work_dir), report=report,
    )


if __name__ == "__main__":
    main()
