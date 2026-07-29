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
import time
from pathlib import Path

import duckdb

LAB_ROOT = Path(__file__).resolve().parent.parent
DB = "/Users/stavanjain/Code/qec-lab/experiments/bb_lab/data/bb_instances.duckdb"

from bb_lab.checks import bb_check_matrices  # noqa: E402
from bb_lab.group import ZmZn  # noqa: E402
from bb_lab.linalg import nullspace_f2, quotient_complement_basis  # noqa: E402
from bb_lab.maxsat_distance import maxsat_distance  # noqa: E402
from bb_lab.poly import Poly  # noqa: E402

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
    ap.add_argument("--per-code-timeout", type=float, default=600.0)
    ap.add_argument("--work-dir", default=".")
    ap.add_argument("--out", default="tandem_verify_results.csv")
    args = ap.parse_args()

    con = duckdb.connect(DB, read_only=True)
    cols = [d[0] for d in con.execute(
        "SELECT * FROM bb_instances LIMIT 0").description]
    sql = SAMPLE_SQL.format(
        c12=CAPS[12], c10=CAPS[10], c8=CAPS[8], c6=CAPS[6],
        c4=CAPS[4], c2=CAPS[2], limit=args.limit,
    )
    rows = [dict(zip(cols, r)) for r in con.execute(sql).fetchall()]

    codes = [
        (r["instance_id"], int(r["ell"]), int(r["m"]),
         r["A_poly"], r["B_poly"], int(r["d_exact"]), r["d_method"])
        for r in rows
    ]
    codes += [(cid, l, m, a, b, d, "published") for cid, l, m, a, b, d in BRAVYI]
    print(f"battery: {len(codes)} codes", flush=True)

    n_pass = n_fail = n_err = 0
    times: list[float] = []
    with Path(args.out).open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["instance_id", "n", "k", "d_known", "d_method",
                    "d_tandem", "parity_flag", "seconds", "status"])
        for cid, ell, m, a, b, d_known, method in codes:
            G = ZmZn(ell, m)
            checks = bb_check_matrices(
                Poly.from_string(a, G), Poly.from_string(b, G))
            V = quotient_complement_basis(
                checks.H_X, nullspace_f2(checks.H_Z))
            all_even = (
                not any(int(r_.sum()) % 2 for r_ in checks.H_X)
                and not any(int(v.sum()) % 2 for v in V)
            )
            extra = ("-cost-step=2",) if all_even else ()
            t0 = time.perf_counter()
            try:
                r = maxsat_distance(
                    checks, args.binary, mode="naive",
                    work_dir=args.work_dir,
                    timeout=args.per_code_timeout, extra_args=extra,
                )
                dt = time.perf_counter() - t0
                ok = r.distance == d_known
                status = "PASS" if ok else "MISMATCH"
                if ok:
                    n_pass += 1
                    times.append(r.solver_seconds)
                else:
                    n_fail += 1
                d_out = r.distance
            except Exception as e:
                dt = time.perf_counter() - t0
                status, d_out = f"error:{type(e).__name__}", None
                n_err += 1
            w.writerow([cid, checks.num_qubits, V.shape[0], d_known,
                        method, d_out, bool(extra), round(dt, 2), status])
            f.flush()
            marker = "" if status == "PASS" else "  <<<<<< ATTENTION"
            print(f"  [{cid[:16]:<16}] n={checks.num_qubits:>3} "
                  f"d={d_known:>2} → {d_out}  ({dt:6.2f}s) {status}{marker}",
                  flush=True)

    times.sort()
    if times:
        med = times[len(times) // 2]
        print(f"\nSUMMARY: {n_pass} pass / {n_fail} MISMATCH / {n_err} error"
              f"   median solver {med:.2f}s, max {times[-1]:.2f}s",
              flush=True)


if __name__ == "__main__":
    main()
