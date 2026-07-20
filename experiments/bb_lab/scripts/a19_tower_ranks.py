"""Gross-template step 1 for Bravyi [[360,12,<=24]] = (30,6), A=x^9+y+y^2, B=y^3+x^25+x^26.

Tower:
    GB (15,3)  --x-->  BY (30,3)
      |                  |
      y                  y
      v                  v
    BX (15,6)  --x-->  C  (30,6)

Computes: k at every level (rung k-equality <=> deck condition (R), per A12),
CSS sanity, and complete weight-6 Z/X-logical searches (meet-in-the-middle,
exact for weight 6) at GB, BX, BY. PAR (odd polynomial weights) => distances even.
"""
import sys, itertools, time
from pathlib import Path

LAB = Path.home() / "Code/qec-lab/experiments/bb_lab"
sys.path.insert(0, str(LAB / "src"))

import numpy as np
from bb_lab.group import ZmZn
from bb_lab.poly import Poly
from bb_lab.checks import bb_check_matrices, assert_css_commutation
from bb_lab.codeparams import code_params
from bb_lab.linalg import rref_f2, rank_f2

CODES = {
    "C  (30,6)": (30, 6, "x^9 + y + y^2", "y^3 + x^25 + x^26"),
    "BX (15,6)": (15, 6, "x^9 + y + y^2", "y^3 + x^10 + x^11"),
    "BY (30,3)": (30, 3, "x^9 + y + y^2", "1 + x^25 + x^26"),
    "GB (15,3)": (15, 3, "x^9 + y + y^2", "1 + x^10 + x^11"),
}

built = {}
for name, (ell, m, As, Bs) in CODES.items():
    G = ZmZn(ell, m)
    A = Poly.from_string(As, G)
    B = Poly.from_string(Bs, G)
    ch = bb_check_matrices(A, B)
    assert_css_commutation(ch)
    p = code_params(ch)
    built[name] = (G, A, B, ch, p)
    print(f"{name}: n={p.n} k={p.k} rank(HX)={p.rank_HX} rank(HZ)={p.rank_HZ}")

kC = built["C  (30,6)"][4].k
for base, rung in [("BX (15,6)", "x-double -> C"), ("BY (30,3)", "y-double -> C")]:
    kb = built[base][4].k
    print(f"(R) rung {base} {rung}: k(base)={kb} k(cover)={kC} -> "
          f"{'R HOLDS (k-preserving)' if kb == kC else 'R FAILS'}")
kGB = built["GB (15,3)"][4].k
for cover in ["BX (15,6)", "BY (30,3)"]:
    kc = built[cover][4].k
    print(f"(R) rung GB (15,3) -> {cover}: k(base)={kGB} k(cover)={kc} -> "
          f"{'R HOLDS (k-preserving)' if kGB == kc else 'R FAILS'}")


def min_logicals_w6(name, kerside, quotside, label):
    """Complete search for weight<=6 vectors in ker(kerside) \\ rowspace(quotside).

    kerside: check matrix H (kernel = candidate operators), quotside: stabilizer rows.
    Meet-in-the-middle over weight-3 halves (PAR: odd weights impossible; w in {2,4,6}).
    Returns list of (weight, vector-support) for nontrivial logicals found.
    """
    H = kerside.astype(np.uint8)
    S = quotside.astype(np.uint8)
    ncols = H.shape[1]
    R, piv = rref_f2(S.copy())
    R = R[~np.all(R == 0, axis=1)]
    rankS = R.shape[0]

    def is_trivial(vec):
        M = np.vstack([R, vec[None, :]])
        return rank_f2(M) == rankS

    cols = H.T.copy()  # cols[i] = syndrome of e_i
    found = []
    # weight 2 and 4 first (cheap), then 6 via MITM on weight-3 halves
    t0 = time.time()
    syn1 = {}
    for i in range(ncols):
        syn1.setdefault(cols[i].tobytes(), []).append((i,))
    # w=2: two equal columns
    for w, half in [(2, 1), (4, 2), (6, 3)]:
        halves = {}
        for comb in itertools.combinations(range(ncols), half):
            s = np.bitwise_xor.reduce(cols[list(comb)], axis=0) if half > 1 else cols[comb[0]]
            halves.setdefault(s.tobytes(), []).append(comb)
        hits = []
        for s, combos in halves.items():
            if len(combos) > 1:
                for c1, c2 in itertools.combinations(combos, 2):
                    supp = set(c1) ^ set(c2)
                    if len(supp) == w:
                        hits.append(tuple(sorted(supp)))
        for supp in set(hits):
            vec = np.zeros(ncols, dtype=np.uint8)
            vec[list(supp)] = 1
            assert not np.any(H @ vec % 2), "not in kernel?!"
            if not is_trivial(vec):
                found.append((w, supp))
        if found:
            break  # minimum weight found; no need to go higher
    dt = time.time() - t0
    if found:
        wmin = min(f[0] for f in found)
        print(f"{name} {label}: MIN LOGICAL WEIGHT = {wmin} "
              f"({sum(1 for f in found if f[0]==wmin)} at min, complete through w=6) [{dt:.1f}s]")
        ex = next(f for f in found if f[0] == wmin)
        print(f"    example support (col indices): {ex[1]}")
    else:
        print(f"{name} {label}: no nontrivial logical of weight <= 6 "
              f"(complete; PAR => d >= 8) [{dt:.1f}s]")
    return found


for name in ["GB (15,3)", "BX (15,6)", "BY (30,3)"]:
    G, A, B, ch, p = built[name]
    min_logicals_w6(name, ch.H_X, ch.H_Z, "Z-logicals")
    min_logicals_w6(name, ch.H_Z, ch.H_X, "X-logicals")
