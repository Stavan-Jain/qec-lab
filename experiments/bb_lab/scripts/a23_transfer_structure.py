"""A23 phase 9: the transfer-operator reduction.

Conjectured structure (from character counting): ker(A*) = ker(B*) = ker d2
(all dim 4, the single F16 character block).  Then B*f is a linear function
of A*f, and the seam coset has a quantifier-free characterization:

  w in coset(zeta)  <=>  w_R = Q * w_L + t'(zeta)   and   e0-pin on w_L,

where P is a quasi-inverse (P*A = 1 + e0, e0 = block idempotent), Q = B*P,
t' = sC_R + Q*sC_L.  Verify all of it numerically and compute the sizes.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from bb_lab.linalg import nullspace_f2, rank_f2, rref_f2  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
from a23_seam_calibration import (  # noqa: E402
    A_SUPP,
    B_SUPP,
    conv_matrix,
    seam_maps,
)


def main() -> None:
    rng = np.random.default_rng(7)
    MA = conv_matrix(A_SUPP)
    MB = conv_matrix(B_SUPP)
    D2 = np.vstack([MA, MB])
    K = nullspace_f2(D2)

    # 1. kernel identity
    kerA = nullspace_f2(MA)
    kerB = nullspace_f2(MB)
    print(f"dim ker A = {kerA.shape[0]}, dim ker B = {kerB.shape[0]}, "
          f"dim ker d2 = {K.shape[0]}")
    rkA = rank_f2(np.vstack([kerA, K]))
    rkB = rank_f2(np.vstack([kerB, K]))
    print(f"ker A == ker d2: {rkA == K.shape[0] == kerA.shape[0]}; "
          f"ker B == ker d2: {rkB == K.shape[0] == kerB.shape[0]}")

    # 2. block idempotent e0: projection onto the F16 block along the rest.
    #    e0 = the unique idempotent with e0*R = ker A.  Compute as the
    #    group-algebra element: solve e0 in ker(A-side)| e0^2 = e0...
    #    Easier: e0 = sum over the block characters; numerically get the
    #    projector matrix Pi0 with im = ker A, ker = im(A-op on ...).
    #    In the group algebra, Pi0 = conv by eps0 where eps0 spans the
    #    1-dim space ker(MA) /\ "translation-eigen"... use: eps0 = the
    #    idempotent generator of ker A as an ideal.  ker A is an ideal
    #    (dim 4); its idempotent e satisfies e*v = v for v in ker A and
    #    e*u = 0 for u in the complementary ideal.
    #    Find e0 in ker A with MA-column-space test: e0 * k = k for k in K.
    #    Solve linear system over the 4-dim ker A: e0 in kerA-span with
    #    conv(e0, K[i]) = K[i].
    from itertools import product

    sols = []
    for coeffs in product((0, 1), repeat=kerA.shape[0]):
        if not any(coeffs):
            continue
        e = np.zeros(75, dtype=np.uint8)
        for c, row in zip(coeffs, kerA):
            if c:
                e ^= row
        # conv e with K[0]: e acts as identity on the block?
        conv_mat_e = conv_matrix(
            [(j // 15, j % 15) for j in np.flatnonzero(e)]
        )
        if all(
            np.array_equal((conv_mat_e @ K[i]) % 2, K[i]) for i in range(K.shape[0])
        ):
            sols.append(e)
    print(f"idempotent candidates in ker A acting as id on block: {len(sols)}")
    e0 = sols[0]
    Me0 = conv_matrix([(j // 15, j % 15) for j in np.flatnonzero(e0)])
    print(f"|e0| = {int(e0.sum())}; e0*e0 == e0: "
          f"{np.array_equal((Me0 @ e0) % 2, e0)}")

    # 3. quasi-inverse P: P*A = 1 + e0  (delta_0 + e0 as vectors)
    target = np.zeros(75, dtype=np.uint8)
    target[0] = 1
    target ^= e0
    # solve MA^T? conv is commutative: P*A = A*P: solve MA @ P = target
    aug = np.hstack([MA, target[:, None]])
    R, piv = rref_f2(aug)
    assert 75 not in piv, "P*A = 1+e0 unsolvable?!"
    P = np.zeros(75, dtype=np.uint8)
    for r, p in enumerate(piv):
        P[p] = R[r, 75]
    assert np.array_equal((MA @ P) % 2, target)
    print(f"quasi-inverse P found, |P| = {int(P.sum())}")

    MP = conv_matrix([(j // 15, j % 15) for j in np.flatnonzero(P)])
    Q = (MB @ P) % 2
    MQ = conv_matrix([(j // 15, j % 15) for j in np.flatnonzero(Q)])
    print(f"|Q| = |B*P| = {int(Q.sum())}")

    # 4. verify the coset characterization on random elements
    zeta = K[0]
    _, sC = seam_maps(zeta)
    sL, sR = sC[:75], sC[75:]
    tprime = (sR ^ (MQ @ sL) % 2) % 2
    print(f"|t'| = |sC_R + Q*sC_L| = {int(tprime.sum())}")
    ok = True
    for _ in range(200):
        f = rng.integers(0, 2, 75).astype(np.uint8)
        wL = (sL + MA @ f) % 2
        wR = (sR + MB @ f) % 2
        # check wR == Q*wL + t'
        pred = (MQ @ wL + tprime) % 2
        if not np.array_equal(pred, wR):
            ok = False
            break
    print(f"identity  w_R = Q*w_L + t'  on 200 random coset elements: {ok}")

    # e0-pin: e0*w_L constant = e0*sC_L over the coset; and nonzero
    pin = (Me0 @ sL) % 2
    ok2 = True
    for _ in range(50):
        f = rng.integers(0, 2, 75).astype(np.uint8)
        wL = (sL + MA @ f) % 2
        if not np.array_equal((Me0 @ wL) % 2, pin):
            ok2 = False
    print(f"e0-pin  e0*w_L = e0*sC_L  on coset: {ok2}; |pin| = {int(pin.sum())}, "
          f"pin nonzero: {pin.any()}")

    # converse: does (pin + relation) => membership?  dimension count:
    # solutions of [wR = Q wL + t', e0 wL = pin] : wL free with e0-pin (dim 71
    # affine) x wR determined => dim 71 = coset dim (im d2 = 71). match?
    print(f"dim im d2 = {rank_f2(D2)} (expect 71 = affine solution dim)")

    # 5. structure of Q for the record: is Q sigma-related to something?
    supQ = sorted((j // 15, j % 15) for j in np.flatnonzero(Q))
    print(f"P support: {sorted((j // 15, j % 15) for j in np.flatnonzero(P))}")
    print(f"Q support: {supQ}")


if __name__ == "__main__":
    main()
