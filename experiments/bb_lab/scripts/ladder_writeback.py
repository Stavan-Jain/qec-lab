"""Write today's ladder settlements back into the corpus DB.

Reads `ladder_results.csv` (rows with status 'ok' and no prior
d_exact), sets `d_exact` + `d_method` + `updated_at` on the matching
corpus rows, and leaves everything else (including the historically
loose `d_ub`) untouched — matching the convention of the existing
'sat' rows. Idempotent: only rows whose stored d_exact is still NULL
are updated, so re-running after a partial write is safe.

Provenance strings record which engine settled each row:
  - the n=168 batches ran before the fork existed → stock MSE-2023
    MaxCDCL ('maxsat-maxcdcl@mse23');
  - the n=154 batches ran under Tandem with the parity-verified
    -cost-step=2 ('maxsat-tandem@mse23+step2').
In both cases the *witness* was independently re-verified at solve
time; the optimality half is the solver's claim — exactly what
d_method exists to record (the certified tier stays 'sat'/Lean).

Also writes the Wang–Mueller [[154,6,16]] result (Tandem, 761 s) if
that exact polynomial pair exists as a corpus row.

Usage (from experiments/bb_lab):
  uv run python scripts/ladder_writeback.py --dry-run
  uv run python scripts/ladder_writeback.py
"""

from __future__ import annotations

import argparse
import csv
import datetime
import shutil
from pathlib import Path

import duckdb

LAB_ROOT = Path(__file__).resolve().parent.parent
DB = Path("/Users/stavanjain/Code/qec-lab/experiments/bb_lab/data/bb_instances.duckdb")

METHOD_BY_N = {
    168: "maxsat-maxcdcl@mse23",
    154: "maxsat-tandem@mse23+step2",
}

WM154 = {  # published Quantum version of arXiv:2408.10001 (v4 row is a typo)
    "ell": 7, "m": 11,
    "A_poly": "1 + x*y + x^3*y^9",
    "B_poly": "1 + x^5*y^8 + x^4*y^9",
    "d": 16, "method": "maxsat-tandem@mse23+step2",
}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--csv", default=str(LAB_ROOT / "ladder_results.csv"))
    ap.add_argument("--db", default=str(DB))
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    rows = list(csv.DictReader(open(args.csv)))
    todo = [
        r for r in rows
        if r["status"] == "ok" and r["d_exact_known"] in ("", "None")
        and int(r["n"]) in METHOD_BY_N
    ]
    print(f"csv: {len(rows)} rows, {len(todo)} settlement candidates")

    if args.dry_run:
        con = duckdb.connect(args.db, read_only=True)
        n_new = n_already = n_missing = 0
        for r in todo:
            row = con.execute(
                "SELECT d_exact FROM bb_instances WHERE instance_id=?",
                [r["instance_id"]],
            ).fetchone()
            if row is None:
                n_missing += 1
            elif row[0] is not None:
                n_already += 1
            else:
                n_new += 1
        wm = con.execute(
            "SELECT instance_id, d_exact FROM bb_instances WHERE ell=? AND m=?"
            " AND A_poly=? AND B_poly=?",
            [WM154["ell"], WM154["m"], WM154["A_poly"], WM154["B_poly"]],
        ).fetchall()
        print(f"dry-run: {n_new} would update, {n_already} already exact, "
              f"{n_missing} ids missing; WM154 corpus match: {wm}")
        return

    stamp = datetime.datetime.now().strftime("%Y%m%d-%H%M%S")
    bak = Path(args.db).with_suffix(f".duckdb.bak-{stamp}")
    shutil.copy2(args.db, bak)
    print(f"backup: {bak}")

    con = duckdb.connect(args.db)
    n_upd = 0
    for r in todo:
        method = METHOD_BY_N[int(r["n"])]
        updated = con.execute(
            "UPDATE bb_instances SET d_exact=?, d_method=?, updated_at=now()"
            " WHERE instance_id=? AND d_exact IS NULL RETURNING instance_id",
            [int(r["d_found"]), method, r["instance_id"]],
        ).fetchall()
        n_upd += len(updated)
    wm_updated = con.execute(
        "UPDATE bb_instances SET d_exact=?, d_method=?, updated_at=now()"
        " WHERE ell=? AND m=? AND A_poly=? AND B_poly=? AND d_exact IS NULL"
        " RETURNING instance_id",
        [WM154["d"], WM154["method"], WM154["ell"], WM154["m"],
         WM154["A_poly"], WM154["B_poly"]],
    ).fetchall()
    con.commit()
    print(f"updated: {n_upd} ladder rows; WM154: {wm_updated}")

    # verification pass
    got = con.execute(
        "SELECT d_method, COUNT(*), MIN(d_exact), MAX(d_exact)"
        " FROM bb_instances WHERE d_method LIKE 'maxsat%' GROUP BY 1"
    ).fetchall()
    print("post-write maxsat provenance rows:", got)
    con.close()


if __name__ == "__main__":
    main()
