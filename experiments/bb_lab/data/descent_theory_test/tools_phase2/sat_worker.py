"""SAT distance ladder worker (one code), run under a parent watchdog.

Writes one line per rung to --progress (w,{SAT|UNSAT},secs) so the
parent can extract the honest UNSAT floor when it kills us on timeout.
Prints 'DISTANCE <d> WITNESS_W <w>' on completion.

Usage: python sat_worker.py --ell L --m M --A 'poly' --B 'poly' \
         --wmax W --progress path
"""

from __future__ import annotations

import argparse
import os
import sys


def main() -> None:
    try:
        os.nice(10)
    except PermissionError:
        pass  # already at max niceness under a nice'd parent
    ap = argparse.ArgumentParser()
    ap.add_argument("--ell", type=int, required=True)
    ap.add_argument("--m", type=int, required=True)
    ap.add_argument("--A", required=True)
    ap.add_argument("--B", required=True)
    ap.add_argument("--wmax", type=int, default=None)
    ap.add_argument("--wmin", type=int, default=1)
    ap.add_argument("--progress", required=True)
    args = ap.parse_args()

    from bb_lab.group import AbelianGroup
    from bb_lab.poly import Poly
    from bb_lab.checks import bb_check_matrices
    from bb_lab.sat_distance import x_distance

    G = AbelianGroup((args.ell, args.m))
    A = Poly.from_string(args.A, G)
    B = Poly.from_string(args.B, G)
    checks = bb_check_matrices(A, B)

    pf = open(args.progress, "a", buffering=1)

    def prog(w: int, sat: bool, secs: float) -> None:
        pf.write(f"{w},{'SAT' if sat else 'UNSAT'},{secs:.2f}\n")
        pf.flush()

    res = x_distance(
        checks,
        weight_lower_bound=args.wmin,
        weight_upper_bound=args.wmax,
        progress=prog,
    )
    print(f"DISTANCE {res.distance} WITNESS_W {int(res.witness.sum())}", flush=True)


if __name__ == "__main__":
    main()
