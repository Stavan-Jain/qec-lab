"""A/B the qec-lab MaxCDCL fork against stock on the session's d≤14 suite.

Configs (naive WCNF encoding throughout — the fork moves the analytic
content into *search mechanics*, not constraints):
  stock        — unpatched MSE-2023 MaxCDCL
  step         — fork, -cost-step=2   (coset weight parity as pruning)
  step+prime   — fork, + -prime-vars=<L0,R0 anchors>
  step+seed    — fork, -cost-step=2, incumbent seeded with known d

The parity flag is only passed after this script VERIFIES the instance
property (every H_X row even AND every logical-basis parity even ⟹ all
coset weights even). Every result is checked against the expected d.

Usage: uv run python scripts/fork_ab.py --stock BIN --fork BIN [--out CSV]
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
from bb_lab.maxsat_distance import maxsat_distance  # noqa: E402
from bb_lab.poly import Poly  # noqa: E402
from bb_lab.sat_distance import find_logical_z  # noqa: E402
from bb_lab.linalg import nullspace_f2, quotient_complement_basis  # noqa: E402

BRAVYI = [
    ("bb_72", 6, 6, "x^3 + y + y^2", "y^3 + x + x^2", 6),
    ("bb_90", 15, 3, "x^9 + y + y^2", "1 + x^2 + x^7", 10),
    ("bb_108", 9, 6, "x^3 + y + y^2", "y^3 + x + x^2", 10),
    ("gross", 12, 6, "x^3 + y + y^2", "y^3 + x + x^2", 12),
]
CORPUS = [  # (label, instance_id, expected d) — settled earlier today
    ("l150_d12", "cf2664aac573dc9e", 12),
    ("l168_d12", "628f9bfe885c9003", 12),
    ("l168_d14a", "24413d6313c49539", 14),
    ("l168_d14b", "2e81450272ba7b92", 14),
    ("l168_d14c", "8ad0d721d845b1c1", 14),
    ("l168_d14d", "0929c2bd0d18d231", 14),
]


def _all_even(checks) -> bool:
    if any(int(r.sum()) % 2 for r in checks.H_X):
        return False
    L_Z = find_logical_z(checks)
    V = quotient_complement_basis(checks.H_X, nullspace_f2(checks.H_Z))
    assert V.shape[0] == L_Z.shape[0]
    return not any(int(v.sum()) % 2 for v in V)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--stock", required=True)
    ap.add_argument("--fork", required=True)
    ap.add_argument("--work-dir", default=".")
    ap.add_argument("--out", default="fork_ab_results.csv")
    args = ap.parse_args()

    codes = []
    for label, ell, m, a, b, d in BRAVYI:
        G = ZmZn(ell, m)
        codes.append((label, bb_check_matrices(
            Poly.from_string(a, G), Poly.from_string(b, G)), d))
    con = duckdb.connect(DB, read_only=True)
    for label, iid, d in CORPUS:
        ell, m, a, b = con.execute(
            "SELECT ell,m,A_poly,B_poly FROM bb_instances WHERE instance_id=?",
            [iid],
        ).fetchone()
        G = ZmZn(int(ell), int(m))
        codes.append((label, bb_check_matrices(
            Poly.from_string(a, G), Poly.from_string(b, G)), d))

    rows = []
    for label, checks, d_exp in codes:
        assert _all_even(checks), f"{label}: parity premise fails"
        n = checks.num_qubits
        anchors = f"1,{n // 2 + 1}"  # qubit vars are ids 1..n; L0 and R0
        configs = [
            ("stock", args.stock, (), None),
            ("step", args.fork, ("-cost-step=2",), None),
            ("step+prime", args.fork,
             ("-cost-step=2", f"-prime-vars={anchors}"), None),
            ("step+seed", args.fork, ("-cost-step=2",), d_exp),
        ]
        times = {}
        for name, binary, extra, seed in configs:
            t0 = time.perf_counter()
            r = maxsat_distance(
                checks, binary, mode="naive", work_dir=args.work_dir,
                extra_args=extra, seed_ub=seed,
            )
            wall = time.perf_counter() - t0
            assert r.distance == d_exp, (label, name, r.distance, d_exp)
            times[name] = r.solver_seconds
            print(f"  {label:>10} {name:<11} d={r.distance} OK "
                  f"{r.solver_seconds:8.2f}s", flush=True)
        rows.append({"code": label, "n": n, "d": d_exp, **times})

    with Path(args.out).open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)

    print("\ncode          n   d    stock     step  step+prime  step+seed")
    for r in rows:
        print(f"{r['code']:<12}{r['n']:>4}{r['d']:>4}"
              f"{r['stock']:>9.2f}{r['step']:>9.2f}"
              f"{r['step+prime']:>12.2f}{r['step+seed']:>11.2f}")


if __name__ == "__main__":
    main()
