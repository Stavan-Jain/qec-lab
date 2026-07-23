"""A23 session 2, phase 1: the A22-fibering site sweep for the final-form
inequality  |A*f + e0| + |B*f| >= 16.

The reduction (to be verified end-to-end here, then formalized):

  * z-fibers: G150 = Z5(x) x Z15(y) partitions into 15 fibers of 5 points
    F_(i,c) = {(i, c+3k) : k in 0..4},  sites (i,c) in Z5 x Z3.
  * delta-type of a fiber vector p in F2^5 BY WEIGHT ALONE:
    |p| in {0,5} -> O,  {1,4} -> M,  {2,3} -> D   (CRT F2^5 = F2 x GF16).
  * parity link: parity(fiber_s(A*f)) = parity(fiber_{s+xbar}(B*f)),
    xbar acting on sites as (i,c) -> (i+1,c); parity(fiber_s(e0)) = 0.
  * per-site cost bound: wp + wq >= c(type wp, type wq) whenever
    wp = wq mod 2, with c(O,O)=0, c(O,D)=2, c(M,M)=2, c(O,M)=4, c(M,D)=4,
    c(D,D)=4.  Every non-(O,O) site costs >= 2.
  * hence |A*f+e0| + |B*f| >= sum_s c(...); if >= 8 active sites, >= 16.
  * else active set inside a 7-site S: f satisfies the 64 affine F2
    conditions "off-S fibers are eps-only"; sweep all C(15,7) = 6435
    systems, enumerate solutions modulo K0 = span{N_s} + ker d2 (dim 19,
    cost-invariant), check cost >= 16 on the 2^k quotient reps.

Verdicts printed as V0..V6; the sweep must find ZERO violations and its
tight set must be exactly the SAT-witness cost pattern.
"""

from __future__ import annotations

import itertools
import json
import sys
import time
from pathlib import Path

import numpy as np

LAB = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(LAB / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from a23_seam_calibration import (  # noqa: E402
    A_SUPP,
    B_SUPP,
    conv_matrix,
)
from bb_lab.linalg import nullspace_f2, rref_f2  # noqa: E402

LX, LY = 5, 15
N = 75

MA = conv_matrix(A_SUPP)
MB = conv_matrix(B_SUPP)
D2 = np.vstack([MA, MB])

E0_SUPP = json.loads((LAB / "data" / "a23" / "final_form.json").read_text())[
    "e0_support"
]
E0 = np.zeros(N, dtype=np.uint8)
E0[E0_SUPP] = 1

# ---------------------------------------------------------------- sites
SITES = [(i, c) for i in range(5) for c in range(3)]
SIDX = {s: k for k, s in enumerate(SITES)}


def fiber_cells(i: int, c: int) -> list[int]:
    return [i * LY + ((c + 3 * k) % LY) for k in range(5)]


FIBERS = np.array([fiber_cells(i, c) for (i, c) in SITES])  # (15,5)

# every cell in exactly one fiber
assert sorted(FIBERS.reshape(-1).tolist()) == list(range(N))


def fiber_wts(v: np.ndarray) -> np.ndarray:
    """15-vector of fiber weights of v (length-75 F2 vector)."""
    return v[FIBERS].sum(axis=1)


def type_of_wt(w: int) -> str:
    return "O" if w in (0, 5) else ("M" if w in (1, 4) else "D")


COST = {
    ("O", "O"): 0, ("O", "M"): 4, ("O", "D"): 2,
    ("M", "O"): 4, ("M", "M"): 2, ("M", "D"): 4,
    ("D", "O"): 2, ("D", "M"): 4, ("D", "D"): 4,
}

XBAR = [SIDX[((i + 1) % 5, c)] for (i, c) in SITES]  # site s -> s + xbar


def cost_of(f: np.ndarray) -> int:
    """Type-based cost sum of the pair (A*f + e0, B*f)."""
    wu = fiber_wts((MA @ f + E0) % 2)
    wv = fiber_wts((MB @ f) % 2)
    return sum(
        COST[(type_of_wt(int(wu[s])), type_of_wt(int(wv[XBAR[s]])))]
        for s in range(15)
    )


def true_total(f: np.ndarray) -> int:
    return int(((MA @ f + E0) % 2).sum() + ((MB @ f) % 2).sum())


def main() -> None:  # noqa: PLR0915
    rng = np.random.default_rng(23230722)
    t0 = time.time()

    # --- V0: weight partition (trivial but pins the convention) --------
    for _ in range(20):
        v = rng.integers(0, 2, N).astype(np.uint8)
        assert int(v.sum()) == int(fiber_wts(v).sum())
    print("V0 fiber partition: PASS")

    # --- V1: e0 fiber profile -----------------------------------------
    we0 = fiber_wts(E0)
    prof = sorted(int(x) for x in we0)
    assert all(w in (2, 4) for w in we0), prof
    print(f"V1 e0 fiber weights all in {{2,4}} (profile {prof}): PASS")

    # --- V2: parity link ----------------------------------------------
    ok = True
    for _ in range(200):
        f = rng.integers(0, 2, N).astype(np.uint8)
        pu = fiber_wts((MA @ f) % 2) % 2
        pv = fiber_wts((MB @ f) % 2) % 2
        for s in range(15):
            if pu[s] != pv[XBAR[s]]:
                ok = False
    assert ok
    # ... and on the delta-basis (the exact Lean obligation shape)
    for g in range(N):
        f = np.zeros(N, dtype=np.uint8)
        f[g] = 1
        pu = fiber_wts((MA @ f) % 2) % 2
        pv = fiber_wts((MB @ f) % 2) % 2
        assert all(pu[s] == pv[XBAR[s]] for s in range(15))
    print("V2 parity link (200 random + 75 basis): PASS")

    # --- V3: per-site cost bound (all weight pairs) --------------------
    for wp in range(6):
        for wq in range(6):
            if (wp - wq) % 2 == 0:
                c = COST[(type_of_wt(wp), type_of_wt(wq))]
                assert wp + wq >= c, (wp, wq)
    for k, v in COST.items():
        if k != ("O", "O"):
            assert v >= 2
    print("V3 per-site cost bound + active >= 2: PASS")

    # --- V4: cost-sum is the eps-min of the true total -----------------
    # K0 eps-part: adding N_s to f changes A*f, B*f by A*N_s, B*N_s which
    # are eps-only; realizing every parity vector h needs Abar invertible.
    NS = np.zeros((15, N), dtype=np.uint8)
    for s in range(15):
        NS[s, FIBERS[s]] = 1
    # check A*N_s, B*N_s are eps-only (all fiber weights in {0,5})
    for s in range(15):
        for M in (MA, MB):
            w = fiber_wts((M @ NS[s]) % 2)
            assert all(int(x) in (0, 5) for x in w)
    mins_ok = 0
    for _ in range(8):
        f = rng.integers(0, 2, N).astype(np.uint8)
        c = cost_of(f)
        best = 10**9
        for mask in range(1 << 15):
            g = f.copy()
            for s in range(15):
                if (mask >> s) & 1:
                    g ^= NS[s]
            best = min(best, true_total(g))
        assert best >= c
        if best == c:
            mins_ok += 1
        assert c <= true_total(f)
    print(f"V4 cost = eps-min of true total: {mins_ok}/8 exact "
          f"(all >= cost): {'PASS' if mins_ok == 8 else 'PARTIAL'}")

    # --- V5: K0 and cost invariance ------------------------------------
    K = nullspace_f2(D2)
    assert len(K) == 4
    K0 = np.vstack([NS, np.array(K, dtype=np.uint8)])
    # rank 19
    R, piv = rref_f2(K0.copy())
    assert len(piv) == 19
    for _ in range(50):
        f = rng.integers(0, 2, N).astype(np.uint8)
        coeff = rng.integers(0, 2, 19).astype(np.uint8)
        kap = (coeff @ K0) % 2
        assert cost_of(f) == cost_of((f + kap) % 2)
    print("V5 K0 (dim 19) cost-invariance: PASS")

    # --- V6: THE SWEEP -------------------------------------------------
    # conditions for off-site g: fiber_g(A*f+e0) eps-only AND
    # fiber_{g+xbar}(B*f) eps-only:  4+4 affine F2 conditions each.
    # rows as (mask75, const).
    def site_conditions(g: int) -> tuple[np.ndarray, np.ndarray]:
        rows, consts = [], []
        cells_u = FIBERS[g]
        for k in range(4):
            rows.append((MA[cells_u[k]] ^ MA[cells_u[k + 1]]))
            consts.append(int(E0[cells_u[k]] ^ E0[cells_u[k + 1]]))
        cells_v = FIBERS[XBAR[g]]
        for k in range(4):
            rows.append((MB[cells_v[k]] ^ MB[cells_v[k + 1]]))
            consts.append(0)
        return np.array(rows, dtype=np.uint8), np.array(consts, dtype=np.uint8)

    SC = [site_conditions(g) for g in range(15)]

    n_cons = 0
    n_incons = 0
    kdims: dict[int, int] = {}
    total_reps = 0
    min_cost = 10**9
    tight = []
    violations = []
    t1 = time.time()
    for S in itertools.combinations(range(15), 7):
        off = [g for g in range(15) if g not in S]
        rows = np.vstack([SC[g][0] for g in off])
        consts = np.concatenate([SC[g][1] for g in off])
        aug = np.hstack([rows, consts[:, None]])
        R, piv = rref_f2(aug.copy())
        if N in piv:
            n_incons += 1
            continue
        n_cons += 1
        # particular solution
        p = np.zeros(N, dtype=np.uint8)
        for r, pv_ in enumerate(piv):
            p[pv_] = R[r, N]
        assert np.array_equal((rows @ p) % 2, consts)
        # homogeneous kernel
        NSol = nullspace_f2(rows)
        nk = len(NSol)
        k_extra = nk - 19
        assert k_extra >= 0, (S, nk)
        kdims[k_extra] = kdims.get(k_extra, 0) + 1
        # quotient reps: rows of NSol modulo K0 span
        # rref of [K0; NSol] -> pivots beyond K0's 19 give extras
        stack = np.vstack([K0, np.array(NSol, dtype=np.uint8)])
        Rs, ps = rref_f2(stack.copy())
        assert len(ps) == 19 + k_extra, (S, len(ps), k_extra)
        extras = Rs[19:19 + k_extra] if k_extra else np.zeros((0, N), np.uint8)
        # sanity: extras satisfy homogeneous system
        if k_extra:
            assert not ((rows @ extras.T) % 2).any()
        for t in range(1 << k_extra):
            f = p.copy()
            for j in range(k_extra):
                if (t >> j) & 1:
                    f ^= extras[j]
            c = cost_of(f)
            total_reps += 1
            if c < min_cost:
                min_cost = c
            if c < 16:
                violations.append((S, t, c))
            elif c == 16:
                tight.append((S, t))
    t2 = time.time()

    print(f"V6 sweep: {n_cons} consistent / {n_incons} inconsistent of 6435")
    print(f"   extra-kernel dims: {kdims}")
    print(f"   total quotient reps checked: {total_reps}")
    print(f"   min cost = {min_cost};  tight (=16): {len(tight)};  "
          f"VIOLATIONS: {len(violations)}")
    if violations:
        print("   !!!", violations[:5])
    print(f"   sweep time {t2 - t1:.1f}s, total {t2 - t0:.1f}s")

    # tight witness cross-check: SAT split (10,6)
    if tight and not violations:
        S, t = tight[0]
        print(f"   sample tight cell: S={S}")


if __name__ == "__main__":
    main()
