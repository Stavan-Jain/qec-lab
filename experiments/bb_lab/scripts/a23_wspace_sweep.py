"""A23 session 2, phase 2: the w-space (delta-coordinate) sweep.

Reformulates the site sweep of a23_site_sweep.py in the 120-dim
delta-data space, the exact shape the Lean certificates consume:

  w(f) : WIdx -> F2,   WIdx = site x Fin 8,  widx(s,r) = 8*siteIdx(s)+r,
    r < 4:  dcoord_r fiber_s(A*f + e0)   (dcoord_i p = p_i + p_4)
    r >= 4: dcoord_{r-4} fiber_{s+xbar}(B*f)

  * realizable w's satisfy 64 affine relations lam_j . w = c_j
    (lam = null space of L^T, rank L = 56, c = lam . w0);
  * site s active <=> its 8-bit block nonzero; cost from 4-bit nibbles
    (type O = 0000, M in {0001,0010,0100,1000,1111}, D = rest);
  * for every sitemask m with popcount <= 7 (16,384 masks): solve the
    64-row restricted affine system; inconsistent -> row-combo cert;
    consistent -> RREF pivot certs + particular + kernel extras
    (delta-normalized), enumerate 2^k reps, check cost >= 16.

Cross-checks the |S| = 7 layer against the f-space sweep and asserts
ZERO violations globally.  This is the validation pass; the emitter
(a23_gen_seam_sweep.py) re-runs this and writes SeamSweepData.lean.
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

from a23_seam_calibration import A_SUPP, B_SUPP, conv_matrix  # noqa: E402
from bb_lab.linalg import nullspace_f2, rref_f2  # noqa: E402

LX, LY, N = 5, 15, 75

MA = conv_matrix(A_SUPP)
MB = conv_matrix(B_SUPP)

E0_SUPP = json.loads((LAB / "data" / "a23" / "final_form.json").read_text())[
    "e0_support"
]
E0 = np.zeros(N, dtype=np.uint8)
E0[E0_SUPP] = 1

SITES = [(i, c) for i in range(5) for c in range(3)]  # siteIdx = 3*i + c
SIDX = {s: k for k, s in enumerate(SITES)}


def fpt(i: int, c: int, k: int) -> int:
    return i * LY + (c + 3 * k) % LY


def xbar(k: int) -> int:
    i, c = SITES[k]
    return SIDX[((i + 1) % 5, c)]


MSET = {0b0001, 0b0010, 0b0100, 0b1000, 0b1111}


def nibble_type(nib: int) -> int:
    return 0 if nib == 0 else (1 if nib in MSET else 2)


CTAB = [[0, 4, 2], [4, 2, 4], [2, 4, 4]]


def build_L_w0() -> tuple[np.ndarray, np.ndarray]:
    L = np.zeros((120, N), dtype=np.uint8)
    w0 = np.zeros(120, dtype=np.uint8)
    for sk, (i, c) in enumerate(SITES):
        for r in range(4):
            L[8 * sk + r] = MA[fpt(i, c, r)] ^ MA[fpt(i, c, 4)]
            w0[8 * sk + r] = E0[fpt(i, c, r)] ^ E0[fpt(i, c, 4)]
        i2, c2 = SITES[xbar(sk)]
        for r in range(4):
            L[8 * sk + 4 + r] = MB[fpt(i2, c2, r)] ^ MB[fpt(i2, c2, 4)]
    return L, w0


def cost_of_w(w: np.ndarray) -> int:
    tot = 0
    for sk in range(15):
        nu = int(sum(int(w[8 * sk + r]) << r for r in range(4)))
        nv = int(sum(int(w[8 * sk + 4 + r]) << r for r in range(4)))
        tot += CTAB[nibble_type(nu)][nibble_type(nv)]
    return tot


def rref_aug_identity(mat: np.ndarray) -> tuple[np.ndarray, np.ndarray, list[int]]:
    """RREF of mat with recorded row transformation: returns (R, T, pivots)
    with R = T @ mat (mod 2), R in RREF."""
    m, n = mat.shape
    R = mat.copy()
    T = np.eye(m, dtype=np.uint8)
    pivots: list[int] = []
    r = 0
    for cidx in range(n):
        piv = None
        for rr in range(r, m):
            if R[rr, cidx]:
                piv = rr
                break
        if piv is None:
            continue
        R[[r, piv]] = R[[piv, r]]
        T[[r, piv]] = T[[piv, r]]
        for rr in range(m):
            if rr != r and R[rr, cidx]:
                R[rr] ^= R[r]
                T[rr] ^= T[r]
        pivots.append(cidx)
        r += 1
        if r == m:
            break
    return R, T, pivots


def sweep() -> dict:
    L, w0 = build_L_w0()
    rankL = len(rref_f2(L.copy())[1])
    assert rankL == 56, rankL
    lam = np.array(nullspace_f2(L.T.copy()), dtype=np.uint8)  # rows lam_j
    assert lam.shape == (64, 120), lam.shape
    assert not ((lam @ L) % 2).any()
    consts = (lam @ w0) % 2

    # spot-verify relations on random f through the real Delta map
    rng = np.random.default_rng(7)
    for _ in range(50):
        f = rng.integers(0, 2, N).astype(np.uint8)
        w = (L @ f + w0) % 2
        assert np.array_equal((lam @ w) % 2, consts)

    stats = {"incons": 0, "cons": 0, "kdims": {}, "reps": 0,
             "min_cost": 10**9, "tight": 0, "violations": 0}
    certs: dict[int, dict] = {}
    t1 = time.time()
    for size in range(8):
        for S in itertools.combinations(range(15), size):
            m = sum(1 << s for s in S)
            on = [8 * s + r for s in S for r in range(8)]
            Lon = lam[:, on] if on else np.zeros((64, 0), dtype=np.uint8)
            # augmented system [Lon | consts]
            aug = np.hstack([Lon, consts[:, None]])
            R, T, piv = rref_aug_identity(aug)
            ncols = len(on)
            if ncols in piv:
                # inconsistent: the pivot row in the const column
                ridx = piv.index(ncols)
                combo = T[ridx]
                assert not ((combo @ Lon) % 2).any() if ncols else True
                assert (combo @ consts) % 2 == 1
                # full-mask support check: combo of full lam rows is
                # supported off the on-coords
                full = (combo @ lam) % 2
                assert not full[on].any()
                certs[m] = {"kind": 0, "rowsel": combo}
                stats["incons"] += 1
                continue
            stats["cons"] += 1
            # particular solution on on-coords
            p_on = np.zeros(ncols, dtype=np.uint8)
            for r_i, c_i in enumerate(piv):
                p_on[c_i] = R[r_i, ncols]
            free = [c for c in range(ncols) if c not in piv]
            k = len(free)
            stats["kdims"][k] = stats["kdims"].get(k, 0) + 1
            # kernel basis, delta-normalized on free cols
            basis = []
            for fc in free:
                vec = np.zeros(ncols, dtype=np.uint8)
                vec[fc] = 1
                for r_i, c_i in enumerate(piv):
                    vec[c_i] = R[r_i, fc]
                basis.append(vec)
            # pivot certs: row r_i of (R,T): combo T[r_i], on-restricted row
            # = e_{piv[r_i]} + free-support
            pivcerts = []
            for r_i, c_i in enumerate(piv):
                combo = T[r_i]
                row_on = (combo @ Lon) % 2
                assert row_on[c_i] == 1
                supp = set(np.flatnonzero(row_on))
                assert supp <= ({c_i} | set(free)), (m, c_i)
                # no const residue: combo pairs to 0 against consts?
                # (not needed for homogeneous use; not asserted)
                pivcerts.append((c_i, combo))
                for b_i, fc in enumerate(free):
                    assert basis[b_i][c_i] == row_on[fc], "duality"
            # embed to 120
            def emb(v_on: np.ndarray) -> np.ndarray:
                w = np.zeros(120, dtype=np.uint8)
                for ii, o in enumerate(on):
                    w[o] = v_on[ii]
                return w
            part = emb(p_on)
            extras = [emb(b) for b in basis]
            frees = [on[fc] for fc in free]
            # verify against full relations
            assert np.array_equal((lam @ part) % 2, consts)
            for e in extras:
                assert not ((lam @ e) % 2).any()
            # reps + costs
            costs = []
            for t in range(1 << k):
                w = part.copy()
                for j in range(k):
                    if (t >> j) & 1:
                        w ^= extras[j]
                c = cost_of_w(w)
                costs.append(c)
                stats["reps"] += 1
                stats["min_cost"] = min(stats["min_cost"], c)
                if c < 16:
                    stats["violations"] += 1
                elif c == 16:
                    stats["tight"] += 1
            certs[m] = {
                "kind": 1, "part": part, "extras": extras, "frees": frees,
                "piv": [(on[c_i], combo) for (c_i, combo) in pivcerts],
                "costs": costs,
            }
    stats["sweep_s"] = round(time.time() - t1, 1)
    return {"L": L, "w0": w0, "lam": lam, "consts": consts,
            "certs": certs, "stats": stats}


def main() -> None:
    out = sweep()
    st = out["stats"]
    print(f"w-sweep over {st['incons'] + st['cons']} masks (<=7): "
          f"{st['cons']} consistent / {st['incons']} inconsistent")
    print(f"  kernel dims: {st['kdims']}; reps {st['reps']}; "
          f"min cost {st['min_cost']}; tight {st['tight']}; "
          f"VIOLATIONS {st['violations']}  [{st['sweep_s']}s]")
    assert st["violations"] == 0
    assert st["min_cost"] == 16

    # cross-check |S|=7 layer against the f-space sweep counts
    n7c = sum(1 for m, c in out["certs"].items()
              if bin(m).count("1") == 7 and c["kind"] == 1)
    k7 = {}
    for m, c in out["certs"].items():
        if bin(m).count("1") == 7 and c["kind"] == 1:
            k = len(c["extras"])
            k7[k] = k7.get(k, 0) + 1
    print(f"  |S|=7 layer: {n7c} consistent, dims {k7} "
          f"(f-space run: 300, {{0: 120, 4: 180}})")
    assert n7c == 300 and k7 == {0: 120, 4: 180}
    print("CROSS-CHECK vs f-space sweep: PASS")


if __name__ == "__main__":
    main()
