"""A19 near-flat gate: the single query that closes d(C) = 24 if UNSAT.

Claim (the near-flat reduction, 2026-07-22): any nontrivial X-logical v of
C = (30,6) with |v| <= 22 must have at most 5 doubly-occupied y-fibers.
Proof: m := #{fibers with both sheets in v}; |p_y(v)| = |v| - 2m.
  [p_y(v)] != 0  =>  |p_y(v)| >= d(BY) = 12  =>  m <= 5.
  [p_y(v)] = 0, b := p_y(v) != 0  =>  b is a stabilizer with |b| = |v|-2m
      <= 22; |b| <= 20 contradicts the certified (M)@24 floors
      (bands <= 20 complete, all UNSAT); |b| = 22 forces m = 0 <= 5.
  b = 0  =>  v = tau(u), |u| <= 11 nonzero-class BY cycle — contradicts
      d(BY) = 12 (and the A16 floor). qed

Hence: UNSAT of [v in ker H_Z(C), v nontrivial, |v| <= 22, #doubles <= 5]
together with the banked certificates gives d(C) >= 24; with the tau-lift
witness, d(C) = 24 EXACT at solver grade.

Usage: uv run scripts/a19_nearflat_gate.py [--max-doubles 5] [--max-weight 22]
       (--max-doubles 0 = the flat stratum, expected fastest)
"""
import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np

LAB = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(LAB / "src"))

import pycryptosat
from pysat.card import CardEnc, EncType
from pysat.formula import IDPool

from bb_lab.checks import bb_check_matrices
from bb_lab.group import AbelianGroup
from bb_lab.poly import Poly
from bb_lab.sat_distance import find_logical_z

COVER = {"frame": (30, 6), "A": "x^9 + y + y^2", "B": "y^3 + x^25 + x^26"}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--max-doubles", type=int, default=5)
    ap.add_argument("--max-weight", type=int, default=22)
    args = ap.parse_args()

    Gc = AbelianGroup(COVER["frame"])
    cc = bb_check_matrices(Poly.from_string(COVER["A"], Gc),
                           Poly.from_string(COVER["B"], Gc))
    L_Z = find_logical_z(cc) % 2
    HZc = cc.H_Z % 2
    nc = cc.num_qubits
    Nc = 180
    fiber = []
    for blk in (0, 1):
        for a in range(30):
            for b in range(3):
                fiber.append((blk * Nc + a * 6 + b, blk * Nc + a * 6 + b + 3))
    assert len(fiber) == 180

    pool = IDPool()
    v = [pool.id() for _ in range(nc)]
    solver = pycryptosat.Solver()
    for row in HZc:
        idx = np.flatnonzero(row)
        if idx.size:
            solver.add_xor_clause([v[i] for i in idx], False)
    sig_vars = []
    for L in L_Z:
        idx = np.flatnonzero(L)
        s = pool.id()
        solver.add_xor_clause([s] + [v[i] for i in idx], False)
        sig_vars.append(s)
    solver.add_clause(sig_vars)                      # nontrivial

    if args.max_doubles == 0:
        for q0, q1 in fiber:                         # flat: no doubles
            solver.add_clause([-v[q0], -v[q1]])
    else:
        a_vars = []
        for q0, q1 in fiber:                         # a <-> v_q0 & v_q1
            a = pool.id()
            solver.add_clause([-a, v[q0]])
            solver.add_clause([-a, v[q1]])
            solver.add_clause([a, -v[q0], -v[q1]])
            a_vars.append(a)
        card = CardEnc.atmost(lits=a_vars, bound=args.max_doubles,
                              vpool=pool, encoding=EncType.seqcounter)
        for cl in card.clauses:
            solver.add_clause(cl)

    card = CardEnc.atmost(lits=v, bound=args.max_weight, vpool=pool,
                          encoding=EncType.seqcounter)
    for cl in card.clauses:
        solver.add_clause(cl)

    print(f"near-flat gate: |v| <= {args.max_weight}, doubles <= "
          f"{args.max_doubles} — solving...", flush=True)
    t0 = time.time()
    sat, model = solver.solve()
    dt = time.time() - t0
    if not sat:
        print(f"UNSAT ({dt:.0f}s) — gate CLOSED for doubles <= "
              f"{args.max_doubles}", flush=True)
        print(json.dumps({"gate": "UNSAT", "max_doubles": args.max_doubles,
                          "max_weight": args.max_weight,
                          "secs": round(dt, 1)}), flush=True)
    else:
        vv = np.array([1 if model[x] else 0 for x in v], dtype=np.uint8)
        assert not (HZc @ vv % 2).any()
        w = int(vv.sum())
        out = LAB / "data" / "a19" / f"nearflat_SAT_wit_w{w}.npy"
        np.save(out, vv)
        print(f"SAT ({dt:.0f}s): NONTRIVIAL LOGICAL OF WEIGHT {w} — "
              f"d(C) <= {w}!  witness -> {out}", flush=True)
        print(json.dumps({"gate": "SAT", "witness_weight": w,
                          "secs": round(dt, 1)}), flush=True)


if __name__ == "__main__":
    main()
