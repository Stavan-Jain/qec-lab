"""A23 calibration: seam-coset structure of f2a6f17e:y ([[150,8,8]] -> [[300,8,16]]).

Everything in the REPO convention (QECLean BBChainComplex/BBCover):

  conv A f (h) = sum_x A(x) f(h-x)          (so matrix M_A[h,g] = A(h-g))
  d2 f  = (A*f | B*f)                       (block 0 = A, block 1 = B)
  d1 c  = B*c_L + A*c_R
  sec (a,b) = (a,b) with y in [0,15); deckS = (0,15)
  seamC zeta (h,j) = (P_j^30 * lift zeta)(h + (0,15)),  P_0 = A, P_1 = B
  seamN zeta (h,j) = (P_j^30 * lift zeta)(h)

Instance: base A = 1 + y + x, B = x y^6 + x y^10 + x^2 y^12 over Z5 x Z15;
cover = literal lift over Z5 x Z30.

Outputs (A23 calibration numbers):
  1. ker d2(base): dim, the 15 nonzero elements, weights, G-orbit structure,
     stabilizers.
  2. raw seam weights |seamC zeta| + parities  (the f=0 coset elements).
  3. delta2-injectivity: seamC zeta not in im d2 for zeta != 0.
  4. attack-line-1 audit: the boundary-part weight range 14 + |seamC zeta|
     vs the LightClassification cutoff 14.
  5. dual-detector calibration: dual cycle space dim, dual homology dim,
     pairing profile of each seam class against dual homology.
"""

from __future__ import annotations

import itertools
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from bb_lab.linalg import nullspace_f2, rank_f2, rref_f2  # noqa: E402

LX, LY = 5, 15
CLY = 30  # cover y

A_SUPP = [(0, 0), (0, 1), (1, 0)]  # 1 + y + x   as (x-exp, y-exp)
B_SUPP = [(1, 6), (1, 10), (2, 12)]  # x y^6 + x y^10 + x^2 y^12


def flat(a: int, b: int, ly: int = LY) -> int:
    return (a % LX) * ly + (b % ly)


def conv_matrix(supp: list[tuple[int, int]], lx: int = LX, ly: int = LY) -> np.ndarray:
    """M[h,g] = P(h-g): column g is the support translated by g."""
    n = lx * ly
    M = np.zeros((n, n), dtype=np.uint8)
    for ga in range(lx):
        for gb in range(ly):
            g = ga * ly + gb
            for (sa, sb) in supp:
                M[((ga + sa) % lx) * ly + ((gb + sb) % ly), g] ^= 1
    return M


def translate2(v75: np.ndarray, g: tuple[int, int]) -> np.ndarray:
    """(g.v)(h) = v(h-g) on C2 = (5,15) flattened."""
    arr = v75.reshape(LX, LY)
    return np.roll(np.roll(arr, g[0], axis=0), g[1], axis=1).reshape(-1)


def translate1(v150: np.ndarray, g: tuple[int, int]) -> np.ndarray:
    """Block-preserving translation on C1."""
    out = np.empty_like(v150)
    for blk in range(2):
        arr = v150[blk * 75 : (blk + 1) * 75].reshape(LX, LY)
        out[blk * 75 : (blk + 1) * 75] = np.roll(
            np.roll(arr, g[0], axis=0), g[1], axis=1
        ).reshape(-1)
    return out


def conv_cover(supp: list[tuple[int, int]], f: np.ndarray) -> np.ndarray:
    """(P * f)(h) = sum_x P(x) f(h-x) on the cover group (5,30)."""
    out = np.zeros_like(f)
    for (sa, sb) in supp:
        out ^= np.roll(np.roll(f, sa, axis=0), sb, axis=1)
    return out


def seam_maps(zeta75: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Return (seamN, seamC) as C1(base) vectors (length 150)."""
    lift = np.zeros((LX, CLY), dtype=np.uint8)
    lift[:, :LY] = zeta75.reshape(LX, LY)
    convA = conv_cover(A_SUPP, lift)
    convB = conv_cover(B_SUPP, lift)
    seamN = np.concatenate([convA[:, :LY].reshape(-1), convB[:, :LY].reshape(-1)])
    seamC = np.concatenate([convA[:, LY:].reshape(-1), convB[:, LY:].reshape(-1)])
    return seamN, seamC


def main() -> None:
    rng = np.random.default_rng(23)

    MA = conv_matrix(A_SUPP)
    MB = conv_matrix(B_SUPP)
    D2 = np.vstack([MA, MB])  # (150, 75): d2 f = D2 @ f
    # d1 c = B*c_L + A*c_R : (75, 150)
    D1 = np.hstack([MB, MA])

    assert not ((D1 @ D2) % 2).any(), "d1 . d2 != 0 -- convention bug"

    # --- 1. kernel of d2 ---------------------------------------------------
    K = nullspace_f2(D2)  # rows = basis
    dimk = K.shape[0]
    print(f"[1] dim ker d2(base) = {dimk}  ({2**dimk - 1} nonzero elements)")

    kernel_elts = []
    for mask in range(1, 2**dimk):
        v = np.zeros(75, dtype=np.uint8)
        for i in range(dimk):
            if (mask >> i) & 1:
                v ^= K[i]
        kernel_elts.append((mask, v))
    wts = sorted(int(v.sum()) for _, v in kernel_elts)
    print(f"    kernel element weights: {wts}")

    # orbit structure under G-translation
    key = {mask: bytes(v.tobytes()) for mask, v in kernel_elts}
    byte_to_mask = {v: m for m, v in key.items()}
    seen: set[int] = set()
    orbits: list[list[int]] = []
    stabs: list[list[tuple[int, int]]] = []
    for mask, v in kernel_elts:
        if mask in seen:
            continue
        orb = set()
        stab = []
        for ga in range(LX):
            for gb in range(LY):
                tv = translate2(v, (ga, gb))
                tb = tv.tobytes()
                assert tb in byte_to_mask, "kernel not G-stable?!"
                orb.add(byte_to_mask[tb])
                if tb == v.tobytes():
                    stab.append((ga, gb))
        orbits.append(sorted(orb))
        stabs.append(stab)
        seen |= orb
    print(f"    G-orbits on ker\\0: {len(orbits)} orbits, sizes "
          f"{[len(o) for o in orbits]}, stab orders {[len(s) for s in stabs]}")
    for i, (o, s) in enumerate(zip(orbits, stabs)):
        print(f"      orbit {i}: size {len(o)}, rep mask {o[0]}, "
              f"stab generated by {s[:4]}{'...' if len(s) > 4 else ''}")

    # --- 2. seam weights ---------------------------------------------------
    print("[2] raw seam weights (f = 0 coset elements):")
    seamCs = {}
    for mask, v in kernel_elts:
        sN, sC = seam_maps(v)
        # kernel element => seamN + seamC = d2 zeta = 0 => seamN == seamC
        assert np.array_equal(sN, sC), f"seamN != seamC for kernel elt {mask}"
        seamCs[mask] = sC
    wt_by_mask = {m: int(s.sum()) for m, s in seamCs.items()}
    spectrum: dict[int, int] = {}
    for m, w in wt_by_mask.items():
        spectrum[w] = spectrum.get(w, 0) + 1
    print(f"    |seamC| spectrum {{weight: count}}: {dict(sorted(spectrum.items()))}")
    odd = [m for m, w in wt_by_mask.items() if w % 2 == 1]
    print(f"    parity: {'ALL EVEN' if not odd else f'ODD at masks {odd}'}")
    print(f"    min |seamC| = {min(wt_by_mask.values())}, "
          f"max = {max(wt_by_mask.values())}")

    # per orbit
    for i, o in enumerate(orbits):
        ws = sorted(wt_by_mask[m] for m in o)
        print(f"      orbit {i} seam weights: {ws}")

    # --- 3. delta2 injectivity --------------------------------------------
    rk_D2 = rank_f2(D2)
    print(f"[3] rank d2 = {rk_D2} (boundaries dim), cycles dim = "
          f"{150 - rank_f2(D1)}, k = {150 - rank_f2(D1) - rk_D2}")
    # basis of boundaries: image of D2 = columns; rref rows of D2^T
    R, piv = rref_f2(D2.T)
    bnd_basis = R[: len(piv)]  # rows span im d2 (as 150-vectors)
    inj = True
    for mask, sC in seamCs.items():
        aug = np.vstack([bnd_basis, sC])
        if rank_f2(aug) == len(piv):  # seamC in boundaries
            inj = False
            print(f"    !! seamC of kernel elt {mask} IS a boundary")
    print(f"    delta2 injective on ker\\0: {inj}")

    # --- 4. attack-line-1 audit -------------------------------------------
    smin = min(wt_by_mask.values())
    print(f"[4] attack-line-1 audit: a |w|<=14 coset element gives boundary "
          f"part |b| <= 14 + |seamC| in [{14 + smin}, {14 + max(wt_by_mask.values())}]")
    print("    LightClassification cutoff = 14  ==> "
          f"gap of {14 + smin - 14} = |seamC|_min; reduction as stated "
          f"{'CLOSES' if smin == 0 else 'DOES NOT CLOSE'}")

    # --- 5. dual-side calibration -----------------------------------------
    # dual cycles: z with z . d2 f = 0 for all f  <=>  D2^T z = 0
    Zdual = nullspace_f2(D2.T)  # rows = basis, dim expect 79
    # dual boundaries: im D1^T (columns of D1^T = rows of D1)
    Rd, pivd = rref_f2(D1)
    dual_bnd = Rd[: len(pivd)]
    print(f"[5] dual cycles dim = {Zdual.shape[0]}, dual boundaries dim = "
          f"{len(pivd)}, dual homology dim = {Zdual.shape[0] - len(pivd)}")

    # pairing profile: for each orbit rep, which dual homology classes pair odd
    # dual homology basis: extend dual_bnd inside Zdual
    # build complement basis greedily
    cur = dual_bnd.copy()
    dual_log_basis = []
    for row in Zdual:
        aug = np.vstack([cur, row])
        if rank_f2(aug) > rank_f2(cur):
            cur = aug
            dual_log_basis.append(row)
    dual_log_basis = np.array(dual_log_basis, dtype=np.uint8)
    print(f"    dual logical basis rows: {dual_log_basis.shape[0]}")
    for i, o in enumerate(orbits):
        rep = o[0]
        sC = seamCs[rep]
        prof = [(int((sC & z).sum() % 2)) for z in dual_log_basis]
        print(f"    orbit {i} rep {rep}: pairing profile vs dual-log basis = {prof}")
        # count of odd-pairing dual classes among 2^8
        # pairing with class sum_{i in S} z_i = sum of profile entries
        n_odd = sum(
            1
            for msk in range(1, 2 ** len(prof))
            if sum(prof[j] for j in range(len(prof)) if (msk >> j) & 1) % 2 == 1
        )
        print(f"      odd-pairing dual classes: {n_odd}/255")

    # sanity: base witness from the lean data json, if present
    import json

    data_p = Path(__file__).resolve().parents[1] / "data/a17/f2a6_z5z30_lean_data.json"
    if data_p.exists():
        data = json.loads(data_p.read_text())
        w = np.zeros(150, dtype=np.uint8)
        for (a, b) in data["ustar_left"]:
            w[flat(a, b)] = 1
        for (a, b) in data["ustar_right"]:
            w[75 + flat(a, b)] = 1
        cyc = not ((D1 @ w) % 2).any()
        in_bnd = rank_f2(np.vstack([bnd_basis, w])) == len(piv)
        print(f"[6] ustar from lean data: weight {int(w.sum())}, cycle={cyc}, "
              f"in boundaries={in_bnd} (expect weight 8, True, False)")
        # reflect-swap into a dual cycle and verify
        def reflect(v75: np.ndarray) -> np.ndarray:
            arr = v75.reshape(LX, LY)
            out = np.zeros_like(arr)
            for a in range(LX):
                for b in range(LY):
                    out[(-a) % LX, (-b) % LY] = arr[a, b]
            return out.reshape(-1)

        z = np.concatenate([reflect(w[75:]), reflect(w[:75])])
        ok = not ((D2.T @ z) % 2).any()
        print(f"    reflect-swap(ustar) is dual cycle: {ok} (weight {int(z.sum())})")
        zb = rank_f2(np.vstack([dual_bnd, z])) == len(pivd)
        print(f"    ... and is a dual boundary: {zb} (expect False = logical-carrying)")


if __name__ == "__main__":
    main()
