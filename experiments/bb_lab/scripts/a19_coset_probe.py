"""A19 doubly-old feasibility probe: element counts of one BY-logical coset.

Decides engineering-vs-moonshot for the doubly-old sector (A19 §4(iii)): the
lift-aware route (A20 pattern) needs, per BY-logical class, the census of
coset elements c with |c| <= 23, then per-element fiber-pinned lift queries.
This probe measures the element-count growth curve on ONE representative
coset (bands 12-16, raw counts, capped) to extrapolate total lift-query
volume across the <= 35 orbit classes.

Method: ISD finds a weight-12 nontrivial BY logical c0; SAT-enumerate the
affine coset c0 + rowspace(H_X) by bands (XOR system b = H_X^T f + c0,
atmost-W, single-model blocking, cap per band).
"""
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
from bb_lab.linalg import nullspace_f2, rref_f2, rank_f2
from bb_lab.poly import Poly

CAP = 20000


def isd_logical(ch, budget_s=90, seed=23):
    K = nullspace_f2(ch.H_Z)
    RX, _ = rref_f2(ch.H_X.copy())
    RX = RX[RX.any(axis=1)]
    rX = RX.shape[0]
    rng = np.random.default_rng(seed)
    best_w, best_v = 10**9, None
    t0 = time.time()
    while time.time() - t0 < budget_s:
        perm = rng.permutation(K.shape[1])
        R, _ = rref_f2(K[:, perm].copy())
        R = R[R.any(axis=1)]
        wts = R.sum(axis=1)
        for j in np.argsort(wts):
            w = int(wts[j])
            if w >= best_w:
                break
            v = np.zeros(K.shape[1], dtype=np.uint8)
            v[perm] = R[j]
            if rank_f2(np.vstack([RX, v[None, :]])) > rX:
                best_w, best_v = w, v
        if best_w == 12:
            break
    return best_w, best_v


def main():
    G = AbelianGroup((30, 3))
    ch = bb_check_matrices(Poly.from_string("x^9 + y + y^2", G),
                           Poly.from_string("1 + x^25 + x^26", G))
    HX = ch.H_X % 2
    n2 = ch.num_qubits
    N = 90

    w0, c0 = isd_logical(ch)
    assert w0 == 12, f"ISD found {w0}, expected 12"
    print(f"anchor: weight-12 BY logical found; censusing its coset", flush=True)

    results = {}
    for W in (12, 14, 16):
        t0 = time.time()
        pool = IDPool()
        b_vars = [pool.id() for _ in range(n2)]
        f_vars = [pool.id() for _ in range(N)]
        solver = pycryptosat.Solver()
        for j in range(n2):
            idx = np.flatnonzero(HX[:, j])
            solver.add_xor_clause([b_vars[j]] + [f_vars[i] for i in idx],
                                  bool(c0[j]))
        card = CardEnc.atmost(lits=b_vars, bound=W, vpool=pool,
                              encoding=EncType.seqcounter)
        for cl in card.clauses:
            solver.add_clause(cl)
        count = 0
        while True:
            sat, model = solver.solve()
            if not sat:
                results[W] = ("complete", count, round(time.time() - t0, 1))
                print(f"band <= {W}: COMPLETE, {count} elements "
                      f"({time.time()-t0:.0f}s)", flush=True)
                break
            b = np.array([1 if model[v] else 0 for v in b_vars],
                         dtype=np.uint8)
            assert int(b.sum()) <= W
            count += 1
            solver.add_clause([-b_vars[j] if b[j] else b_vars[j]
                               for j in range(n2)])
            if count >= CAP:
                results[W] = ("CAP", count, round(time.time() - t0, 1))
                print(f"band <= {W}: CAP {CAP} hit ({time.time()-t0:.0f}s)",
                      flush=True)
                break
            if count % 2000 == 0:
                print(f"  <= {W}: {count} elements...", flush=True)

    print(json.dumps({"coset_probe": {str(k): v for k, v in results.items()}}),
          flush=True)
    # extrapolation guide: raw stabilizer-band counts for comparison are
    # ~90x the orbit-reduced census ({12: 42, 14: 54, 16: 487} classes).


if __name__ == "__main__":
    main()
