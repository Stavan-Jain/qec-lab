"""A23 phase 11b: exact small-strata verification, sigma side (a <= 7 case).

Derivation (sigma = the order-2 automorphism x -> x y^6, y -> y^4, with
B = x y^6 sigma(A), sigma(e0) = e0, P_1 := sigma(P) satisfies
P_1 * sigma(A) = 1 + e0):

  a-side case: |A*f + e0| = a <= 7.  Writing u'' for the residual and
  changing variables by sigma, the obligation becomes

      forall w in V with |w + e0| = a:   |Psi * w| >= 16 - a,

  same Psi = sigma(A)*P and V = im(A*.) as the b-side.  Enumeration:
  w = e0 + u'' with |u''| = a and GF(16) syndrome S(u'') = S(e0)
  (the coset layer; S(e0) = 1 on the block).

This script checks a <= 5 exactly (coset layers ~ N_a-sized).
"""

from __future__ import annotations

import sys
import time
from itertools import combinations, product
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from bb_lab.linalg import nullspace_f2, rref_f2  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
from a23_seam_calibration import A_SUPP, B_SUPP, conv_matrix  # noqa: E402


def main() -> None:
    MA = conv_matrix(A_SUPP)
    MB = conv_matrix(B_SUPP)
    D2 = np.vstack([MA, MB])
    K = nullspace_f2(D2)
    e0 = None
    for coeffs in product((0, 1), repeat=4):
        if not any(coeffs):
            continue
        e = np.zeros(75, dtype=np.uint8)
        for c, row in zip(coeffs, K):
            if c:
                e ^= row
        Me = conv_matrix([(j // 15, j % 15) for j in np.flatnonzero(e)])
        if all(np.array_equal((Me @ K[i]) % 2, K[i]) for i in range(4)):
            e0 = e
            break
    assert e0 is not None
    target = np.zeros(75, dtype=np.uint8)
    target[0] = 1
    target ^= e0
    aug = np.hstack([MA, target[:, None]])
    R, piv = rref_f2(aug)
    P = np.zeros(75, dtype=np.uint8)
    for r, p in enumerate(piv):
        P[p] = R[r, 75]
    MSA = conv_matrix([(0, 0), (0, 4), (1, 6)])
    Psi = (MSA @ P) % 2
    MPsi = conv_matrix([(j // 15, j % 15) for j in np.flatnonzero(Psi)])

    # sanity of the sigma-side derivation on random h:
    # |A h + e0| = a with w := sigma-change ... verify the equivalent claim
    # directly: for random h, set aa = |A h + e0|, bb = |sigma(A) h ... |
    # Instead verify the identity used:  for w in V (w = A g):
    #   |Psi w| = |sigma(A) g|  up to the killed block?  Psi*A = sigma(A)(1+e0):
    rng = np.random.default_rng(11)
    for _ in range(50):
        g = rng.integers(0, 2, 75).astype(np.uint8)
        w = (MA @ g) % 2
        lhs = int(((MPsi @ w) % 2).sum())
        rhs = int(((MSA @ ((g + Me @ g) % 2)) % 2).sum())  # sigma(A)(1+e0)g
        assert lhs == rhs, (lhs, rhs)
    print("identity Psi*A = sigma(A)*(1+e0) verified on random g")

    # syndromes
    kerAT = nullspace_f2(MA.T)
    syn = [int(sum(int(kerAT[bi, j]) << bi for bi in range(4))) for j in range(75)]

    # S(e0): e0 has syndrome components <e0, kerAT_i>
    s_e0 = 0
    for bi in range(4):
        s_e0 |= int((e0 & kerAT[bi]).sum() % 2) << bi
    print(f"S(e0) = {s_e0} (nonzero expected)")
    assert s_e0 != 0

    def pack(v: np.ndarray) -> int:
        return int.from_bytes(np.packbits(v).tobytes(), "big")

    colPsi = [pack(MPsi[:, g]) for g in range(75)]
    psi_e0 = pack((MPsi @ e0) % 2)

    t0 = time.time()
    for a in range(1, 6):
        mn = 10**9
        cnt = 0
        viol = 0
        for comb in combinations(range(75), a):
            s = 0
            for g in comb:
                s ^= syn[g]
            if s != s_e0:
                continue
            cnt += 1
            # w = e0 + u'' ; |Psi w| = |Psi e0 + sum cols|
            acc = psi_e0
            for g in comb:
                acc ^= colPsi[g]
            b = bin(acc).count("1")
            if b < mn:
                mn = b
            if a + b < 16:
                viol += 1
        print(f"sigma-side stratum a={a}: {cnt} elements, min b = "
              f"{mn if cnt else '-'}, min total = {mn + a if cnt else '-'}, "
              f"violations(<16): {viol} [{time.time()-t0:.0f}s]")


if __name__ == "__main__":
    main()
