"""A23 phase 11: exact small-strata verification of the final-form inequality.

For u in the weight-b layer of V = im(A*.) (b <= 5 here), check
    |Psi*u + e0| >= 16 - b,      Psi = sigma(A)*P  (so Psi*A*h = sigma(A)*h
                                  up to the killed block),
and record the exact per-stratum minimum of a = |Psi*u + e0|.

Layer enumeration: u in V <=> the GF(16) chi0-syndromes of its cells XOR
to 0.  We enumerate supports directly with a syndrome MITM (b <= 5:
1.16M elements, exact).

Validation: min over checked strata of (a + b) must be >= 16, with
equality expected at (a,b) = (10,6)-adjacent strata only (b = 6 not
checked here; b <= 5 should show slack >= 2).
"""

from __future__ import annotations

import sys
import time
from itertools import combinations
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from bb_lab.linalg import nullspace_f2, rref_f2  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
from a23_seam_calibration import (  # noqa: E402
    A_SUPP,
    B_SUPP,
    conv_matrix,
)


def main() -> None:
    MA = conv_matrix(A_SUPP)
    MB = conv_matrix(B_SUPP)
    D2 = np.vstack([MA, MB])
    K = nullspace_f2(D2)

    # e0
    from itertools import product

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

    # P with P*A = 1 + e0, Psi = sigma(A) * P where sigma(A) = 1 + y^4 + x y^6
    target = np.zeros(75, dtype=np.uint8)
    target[0] = 1
    target ^= e0
    aug = np.hstack([MA, target[:, None]])
    R, piv = rref_f2(aug)
    P = np.zeros(75, dtype=np.uint8)
    for r, p in enumerate(piv):
        P[p] = R[r, 75]
    SA_SUPP = [(0, 0), (0, 4), (1, 6)]  # sigma(A) = 1 + y^4 + x y^6
    MSA = conv_matrix(SA_SUPP)
    Psi = (MSA @ P) % 2
    MPsi = conv_matrix([(j // 15, j % 15) for j in np.flatnonzero(Psi)])
    print(f"|Psi| = {int(Psi.sum())}")

    # sanity: for random h, |sigma(A) h + e0| == |Psi (A h) + e0|
    rng = np.random.default_rng(5)
    for _ in range(50):
        h = rng.integers(0, 2, 75).astype(np.uint8)
        lhs = int(((MSA @ h + e0) % 2).sum())
        rhs = int(((MPsi @ ((MA @ h) % 2) + e0) % 2).sum())
        assert lhs == rhs, (lhs, rhs)
    print("identity |sigma(A) h + e0| = |Psi (A h) + e0| verified")

    # chi0 syndrome per cell: 4-bit GF(16) value; u in V <=> XOR of cell
    # syndromes = 0.  Get syndromes from the kernel basis pairing:
    # syndrome components = <u, k_i> for the 4 dual kernel elements
    # (ker A^T basis).
    kerAT = nullspace_f2(MA.T)
    assert kerAT.shape[0] == 4
    syn = np.zeros(75, dtype=np.uint8)
    for j in range(75):
        s = 0
        for bi in range(4):
            s |= int(kerAT[bi, j]) << bi
        syn[j] = s
    # columns of Psi + e0 offset, packed as 75-bit ints for speed
    def pack(v: np.ndarray) -> int:
        return int.from_bytes(np.packbits(v).tobytes(), "big")

    colPsi = [pack(MPsi[:, g]) for g in range(75)]
    e0p = pack(e0)

    t0 = time.time()
    for b in range(2, 6):
        mn = 10**9
        cnt = 0
        viol = 0
        # enumerate supports with syndrome XOR 0
        for comb in combinations(range(75), b):
            s = 0
            for g in comb:
                s ^= int(syn[g])
            if s:
                continue
            cnt += 1
            acc = e0p
            for g in comb:
                acc ^= colPsi[g]
            a = bin(acc).count("1")
            if a < mn:
                mn = a
            if a + b < 16:
                viol += 1
        print(f"stratum b={b}: {cnt} elements, min a = {mn}, "
              f"min total = {mn + b}, violations(<16): {viol} "
              f"[{time.time()-t0:.0f}s]")


if __name__ == "__main__":
    main()
