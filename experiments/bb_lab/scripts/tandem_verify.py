"""Tandem validation battery: ~100 codes with independently-known d.

Ground truth: corpus rows whose `d_exact` was settled by the repo's
pre-Tandem SAT pipelines (`d_method` 'sat' / 'sat-cadical@…'), plus the
published Bravyi-table codes — i.e. distances from a different solver
AND a different encoding than anything Tandem touches. Stratified
toward large n and d (the informative region), capped per (d, n)
bucket for diversity.

Per code: the coset-parity premise is verified and `-cost-step=2` is
passed only when it holds (so both fork paths get exercised); the
returned distance must equal the stored ground truth. Mismatches are
recorded and the sweep continues — a single mismatch is a red-alert
finding either way.

Usage (from experiments/bb_lab):
  uv run python scripts/tandem_verify.py --binary <tandem> [--limit 96]
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

LAB_ROOT = Path(__file__).resolve().parent.parent
DB = "/Users/stavanjain/Code/qec-lab/experiments/bb_lab/data/bb_instances.duckdb"

from bb_lab.sweep import (  # noqa: E402
    bb_distance_task, default_jobs, duckdb_rows, run_sweep,
)

FIELDNAMES = ["instance_id", "n", "k", "d_known", "d_method", "d_tandem",
              "parity_flag", "seconds", "status"]


def verify_task(payload: dict) -> dict:
    """Solve, then compare against the independently-known distance.

    A MISMATCH is the whole point of this battery, so it is recorded as
    a row rather than raised — one bad row must not stop the other 99.
    """
    row = bb_distance_task(payload)
    row["d_tandem"] = row.pop("d")
    row["parity_flag"] = bool(row.pop("cost_step"))
    if row["status"] == "ok":
        row["status"] = (
            "PASS" if row["d_tandem"] == row["d_known"] else "MISMATCH"
        )
    return row

# Stratified sample: per (d_exact, n) bucket, up to `cap` rows, largest
# n first inside each d; d weighted by the caps below.
CAPS = {12: 99, 10: 5, 8: 4, 6: 4, 4: 3, 2: 2}  # per-(d,n)-bucket caps

SAMPLE_SQL = """
WITH ranked AS (
  SELECT *, ROW_NUMBER() OVER (
      PARTITION BY d_exact, n ORDER BY instance_id
  ) AS rk
  FROM bb_instances
  WHERE d_exact IS NOT NULL AND d_exact BETWEEN 2 AND 12
    AND d_method IS NOT NULL
)
SELECT * FROM ranked
WHERE rk <= CASE d_exact
    WHEN 12 THEN {c12} WHEN 10 THEN {c10} WHEN 8 THEN {c8}
    WHEN 6 THEN {c6} WHEN 4 THEN {c4} ELSE {c2} END
ORDER BY d_exact DESC, n DESC
LIMIT {limit}
"""

BRAVYI = [
    ("bb_72_12_6", 6, 6, "x^3 + y + y^2", "y^3 + x + x^2", 6),
    ("bb_90_8_10", 15, 3, "x^9 + y + y^2", "1 + x^2 + x^7", 10),
    ("bb_108_8_10", 9, 6, "x^3 + y + y^2", "y^3 + x + x^2", 10),
    ("gross", 12, 6, "x^3 + y + y^2", "y^3 + x + x^2", 12),
]


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--binary", required=True)
    ap.add_argument("--limit", type=int, default=96)
    ap.add_argument("--jobs", type=int, default=None,
                    help="workers (default: performance-core count)")
    ap.add_argument("--per-code-timeout", type=float, default=1800.0)
    ap.add_argument("--resume", action="store_true",
                    help="keep an existing --out and solve only what is "
                         "missing; by default the battery starts fresh, "
                         "since a validation run wants every code checked "
                         "against the binary in hand")
    ap.add_argument("--work-dir", default="verify_work")
    ap.add_argument("--out", default="tandem_verify_results.csv")
    args = ap.parse_args()

    from bb_lab.checks import bb_check_matrices
    from bb_lab.codeparams import code_params
    from bb_lab.group import ZmZn
    from bb_lab.poly import Poly

    rows = duckdb_rows(DB, SAMPLE_SQL.format(
        c12=CAPS[12], c10=CAPS[10], c8=CAPS[8], c6=CAPS[6],
        c4=CAPS[4], c2=CAPS[2], limit=args.limit,
    ))
    codes = [
        (r["instance_id"], int(r["ell"]), int(r["m"]), r["A_poly"],
         r["B_poly"], int(r["d_exact"]), r["d_method"],
         int(r["n"]), int(r["k"]))
        for r in rows
    ]
    for cid, ell, m, a, b, d in BRAVYI:  # published rows carry no (n, k)
        G = ZmZn(ell, m)
        params = code_params(bb_check_matrices(
            Poly.from_string(a, G), Poly.from_string(b, G)))
        codes.append((cid, ell, m, a, b, d, "published", params.n, params.k))

    items = [
        {
            "ell": ell, "m": m, "A_poly": a, "B_poly": b,
            "binary": args.binary, "mode": "naive",
            "timeout": args.per_code_timeout,
            "passthrough": {"instance_id": cid, "n": n, "k": k,
                            "d_known": d_known, "d_method": method},
        }
        for cid, ell, m, a, b, d_known, method, n, k in codes
    ]

    out = Path(args.out)
    if not args.resume:
        out.unlink(missing_ok=True)
    print(f"battery: {len(items)} codes, "
          f"jobs={args.jobs or default_jobs()}", flush=True)

    def report(row: dict) -> str | None:
        marker = "" if row["status"] == "PASS" else "  <<<<<< ATTENTION"
        return (f"  [{row['instance_id'][:16]:<16}] n={row['n']:>3} "
                f"d={row['d_known']:>2} → {row['d_tandem']}  "
                f"({row['seconds']:6.2f}s) {row['status']}{marker}")

    run_sweep(
        items, verify_task,
        out=out, fieldnames=FIELDNAMES, key_field="instance_id",
        key=lambda it: it["passthrough"]["instance_id"],
        jobs=args.jobs, work_root=Path(args.work_dir), report=report,
    )

    with out.open(newline="") as f:
        done = list(csv.DictReader(f))
    n_pass = sum(r["status"] == "PASS" for r in done)
    n_fail = sum(r["status"] == "MISMATCH" for r in done)
    n_err = len(done) - n_pass - n_fail
    times = sorted(float(r["seconds"]) for r in done if r["status"] == "PASS")
    print(f"\nSUMMARY: {n_pass} pass / {n_fail} MISMATCH / {n_err} error",
          end="", flush=True)
    if times:
        print(f"   median solver {times[len(times) // 2]:.2f}s, "
              f"max {times[-1]:.2f}s", flush=True)
    else:
        print(flush=True)
    if n_fail:
        raise SystemExit(f"{n_fail} MISMATCH — a wrong distance is a "
                         f"red-alert finding either way")


if __name__ == "__main__":
    main()
