"""Write merit-sweep / distance-battery settlements back into the corpus DB.

Reads result CSVs (any of the sweep or battery formats -- all carry
`instance_id`, `d`, `cost_step`, `status`), sets `d_exact` +
`d_method` + `updated_at` on rows whose stored d_exact is still NULL,
and leaves everything else -- including the historically loose `d_ub`
-- untouched, matching scripts/ladder_writeback.py.

Idempotent: the UPDATE is guarded on `d_exact IS NULL`, so re-running
after a partial write is safe.  A row whose stored d_exact disagrees
with the CSV is reported as a CONFLICT and aborts the write rather
than being silently overwritten.

Provenance records which solver path settled each row:
  - 'maxsat-tandem@mse23+step2' -- the fork with the parity-verified
    -cost-step=2 (the premise held for that code);
  - 'maxsat-tandem@mse23'       -- the fork with no flags, i.e. the
    behaviourally-stock path, taken where the premise failed.
The *witness* is re-verified at solve time by maxsat_distance; the
optimality half is the solver's claim, which is what d_method exists
to record (the certified tier stays 'sat' / Lean).

Validate the solve path before trusting a large write: run the same
harness over rows already settled by a different engine (d_method
'sat%') and require exact agreement.  The 2026-07-29 write did this
over 36 stratified rows -- 36/36 -- before committing 28,003 updates.

Usage (from experiments/bb_lab):
  uv run python scripts/merit_writeback.py --csv 'results/shard*.csv' --dry-run
  uv run python scripts/merit_writeback.py --csv 'results/shard*.csv'
"""

from __future__ import annotations

import argparse
import csv
import datetime
import glob
import shutil
from pathlib import Path

import duckdb

DB = Path("/Users/stavanjain/Code/qec-lab/experiments/bb_lab/data/bb_instances.duckdb")

METHOD = {1: "maxsat-tandem@mse23+step2", 0: "maxsat-tandem@mse23"}


def load(patterns: str) -> list[dict]:
    rows: list[dict] = []
    for pat in patterns.split(","):
        for f in sorted(glob.glob(pat)):
            rows += [r for r in csv.DictReader(open(f))
                     if r["status"] == "ok" and r["d"]]
    return rows


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--csv", required=True,
                    help="comma-separated globs of result CSVs")
    ap.add_argument("--db", default=str(DB))
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    rows = load(args.csv)
    ids = [r["instance_id"] for r in rows]
    print(f"{len(rows)} settled rows, {len(set(ids))} distinct instance_ids")
    by_method: dict[str, int] = {}
    for r in rows:
        m = METHOD[int(r["cost_step"])]
        by_method[m] = by_method.get(m, 0) + 1
    print("  provenance split:", by_method)

    con = duckdb.connect(args.db, read_only=args.dry_run)
    n_new = n_already = n_missing = n_conflict = 0
    for r in rows:
        got = con.execute(
            "SELECT d_exact FROM bb_instances WHERE instance_id=?",
            [r["instance_id"]],
        ).fetchone()
        if got is None:
            n_missing += 1
        elif got[0] is None:
            n_new += 1
        elif got[0] != int(r["d"]):
            n_conflict += 1
            print(f"  CONFLICT {r['instance_id']}: stored {got[0]}, "
                  f"swept {r['d']}")
        else:
            n_already += 1
    print(f"  would update {n_new}; {n_already} already exact and agree; "
          f"{n_missing} ids missing; {n_conflict} CONFLICTS")
    if n_conflict:
        raise SystemExit("conflicts present — not writing")
    if args.dry_run:
        return

    stamp = datetime.datetime.now().strftime("%Y%m%d-%H%M%S")
    bak = Path(args.db).with_suffix(f".duckdb.bak-{stamp}")
    con.close()
    shutil.copy2(args.db, bak)
    print(f"backup: {bak}")

    con = duckdb.connect(args.db)
    n_upd = 0
    for r in rows:
        n_upd += len(con.execute(
            "UPDATE bb_instances SET d_exact=?, d_method=?, updated_at=now()"
            " WHERE instance_id=? AND d_exact IS NULL RETURNING instance_id",
            [int(r["d"]), METHOD[int(r["cost_step"])], r["instance_id"]],
        ).fetchall())
    con.commit()
    print(f"updated: {n_upd} rows")

    print("post-write d_method census:")
    for m, c, lo, hi in con.execute(
        "SELECT d_method, COUNT(*), MIN(d_exact), MAX(d_exact)"
        " FROM bb_instances WHERE d_exact IS NOT NULL GROUP BY 1 ORDER BY 2 DESC"
    ).fetchall():
        print(f"  {m:<32} {c:>6} rows, d in [{lo},{hi}]")
    print("corpus coverage:", con.execute(
        "SELECT COUNT(d_exact), COUNT(*) FROM bb_instances").fetchone())
    con.close()


if __name__ == "__main__":
    main()
