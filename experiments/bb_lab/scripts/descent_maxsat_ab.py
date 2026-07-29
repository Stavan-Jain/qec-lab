"""A/B: Tandem-naive vs the descent sector decomposition (experiment 1).

Configs per code (all under the Tandem fork binary, -cost-step=2
auto-gated per instance/side):
  naive          — the measured-best solve-mode baseline
  descent-min    — d = min(2·opt_a, opt_b); (b) = naive + sector clause
  descent-rows   — (b) + implied base parity rows
  descent-full   — (b) + rows + Λ-transport link

Every distance is asserted against the expected value; every witness
is re-verified independently of the solver.

Usage (from experiments/bb_lab):
  uv run python scripts/descent_maxsat_ab.py --binary <tandem> \
      [--codes gross l168_d14a] [--reps 3]
"""

from __future__ import annotations

import argparse
import time
from pathlib import Path

import duckdb
import yaml

LAB_ROOT = Path(__file__).resolve().parent.parent
DB = "/Users/stavanjain/Code/qec-lab/experiments/bb_lab/data/bb_instances.duckdb"

from bb_lab.checks import bb_check_matrices  # noqa: E402
from bb_lab.descent_sat import axis_decks  # noqa: E402
from bb_lab.group import ZmZn  # noqa: E402
from bb_lab.linalg import nullspace_f2, quotient_complement_basis  # noqa: E402
from bb_lab.maxsat_distance import (  # noqa: E402
    maxsat_distance,
    maxsat_distance_descent,
)
from bb_lab.poly import Poly  # noqa: E402
from bb_lab.sat_distance import find_logical_z  # noqa: E402

CORPUS_IDS = {  # settled earlier today (ladder sweeps)
    "l150_d12": ("cf2664aac573dc9e", 12),
    "l168_d12": ("628f9bfe885c9003", 12),
    "l168_d14a": ("24413d6313c49539", 14),
}


def load_code(name: str):
    if name in CORPUS_IDS:
        iid, d = CORPUS_IDS[name]
        con = duckdb.connect(DB, read_only=True)
        row = con.execute(
            "SELECT ell, m, A_poly, B_poly FROM bb_instances "
            "WHERE instance_id = ?", [iid],
        ).fetchone()
        con.close()
        ell, m, a_s, b_s = row
        G = ZmZn(int(ell), int(m))
        return G, Poly.from_string(a_s, G), Poly.from_string(b_s, G), d
    table = {
        r["code_id"]: r
        for r in yaml.safe_load(
            (LAB_ROOT / "instances" / "bravyi_table.yaml").read_text()
        )["instances"]
    }
    row = table[name]
    G = ZmZn(row["group"]["ell"], row["group"]["m"])
    return (
        G,
        Poly.from_string(row["polynomials"]["A"], G),
        Poly.from_string(row["polynomials"]["B"], G),
        row["parameters"]["d"],
    )


def _all_even(checks) -> bool:
    if any(int(r.sum()) % 2 for r in checks.H_X):
        return False
    V = quotient_complement_basis(checks.H_X, nullspace_f2(checks.H_Z))
    L_Z = find_logical_z(checks)
    assert V.shape[0] == L_Z.shape[0]
    return not any(int(v.sum()) % 2 for v in V)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--binary", required=True)
    ap.add_argument(
        "--codes", nargs="*", default=["bb_108_8_10", "gross", "l168_d14a"]
    )
    ap.add_argument("--reps", type=int, default=1)
    ap.add_argument("--work", default=None)
    args = ap.parse_args()
    work = Path(args.work) if args.work else LAB_ROOT / "wcnf_tmp"

    rows = []
    for name in args.codes:
        G, A, B, d_exp = load_code(name)
        checks = bb_check_matrices(A, B)
        decks = axis_decks(checks)
        assert decks, f"{name}: no axis deck"
        sigma = decks[-1]  # prefer the second axis (the empirical winner)
        for rep in range(args.reps):
            step = ("-cost-step=2",) if _all_even(checks) else ()
            t0 = time.perf_counter()
            rn = maxsat_distance(
                checks, args.binary, mode="naive",
                work_dir=work, extra_args=step,
            )
            assert rn.distance == d_exp, (name, "naive", rn.distance)
            naive_s = rn.solver_seconds
            print(
                f"{name:<12} rep{rep} naive        d={rn.distance}  "
                f"{naive_s:7.2f}s",
                flush=True,
            )
            rows.append((name, "naive", rep, naive_s, None, None))
            for variant in ("min", "rows", "full"):
                rd = maxsat_distance_descent(
                    checks, A, B, args.binary, sigma=sigma,
                    variant=variant, work_dir=work,
                )
                assert rd.distance == d_exp, (name, variant, rd.distance)
                tot = rd.solver_seconds
                print(
                    f"{name:<12} rep{rep} descent-{variant:<5}"
                    f"d={rd.distance}  {tot:7.2f}s"
                    f"  (a: {rd.seconds_a:5.2f}s opt={rd.opt_a}, "
                    f"b: {rd.seconds_b:6.2f}s opt={rd.opt_b}, "
                    f"sector={rd.sector}, vs naive "
                    f"{naive_s / max(tot, 1e-9):.2f}×)",
                    flush=True,
                )
                rows.append(
                    (name, f"descent-{variant}", rep, tot,
                     rd.seconds_a, rd.seconds_b)
                )

    print("\n=== medians per config ===")
    import statistics
    for name in args.codes:
        for cfg in ("naive", "descent-min", "descent-rows", "descent-full"):
            ts = [r[3] for r in rows if r[0] == name and r[1] == cfg]
            if ts:
                print(f"{name:<12} {cfg:<14} {statistics.median(ts):7.2f}s")


if __name__ == "__main__":
    main()
