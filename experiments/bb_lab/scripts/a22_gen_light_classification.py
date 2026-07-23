"""A22 session-2 generator: certificate data for the Lean `lightClassification`.

Emits `LightCertData.lean` + `LightTTData.lean` into the QECLean worktree
(QECLEAN_ROOT), certifying the analytic light-boundary classification of
`bb_neigh_z5z15_f2a6f17e` ([[150,8,8]], A = 1+y+x, B = xy^6+xy^10+x^2y^12
over Z5 x Z15) in the F2-ified (eps,delta) fibering convention:

  * site s in [0,15) <-> (i, b) = (s // 3, s % 3): x-coord i in Z5,
    w-coord b in Z3 (w = y^5); fiber coordinate a in Z5 (z = y^3);
    the cell of (s, a) on the u-block is (i, (3a+5b) % 15, blk=0), and
    the site's PAIRED v-fiber is ((i+1) % 5, (3a+5b) % 15, blk=1)
    (the xbar-pairing: v's site s+(1,0) shares u's eps-coordinate).
  * delta-nibble of a 5-bit fiber t: d_n = t_n ^ t_4 (reduction mod
    q(z) = 1+z+z^2+z^3+z^4); eps = parity(t).  (eps,delta) <-> t is a
    bijection with explicit inverse t_4 = eps ^ parity(d), t_n = d_n ^ t_4.
  * alpha-beta data of a 1-chain pair: position p = 8s + 4t + n,
    t = 0: alpha (u-side delta), t = 1: beta' (paired v-side delta).

Certificate content (all numpy-hard-asserted before writing):
  * COLPACK: transposed delta-data matrix of d2 (120 x 75-bit) — the
    orthogonality oracle for pivot functionals.
  * REP7: the 429 canonical size-7 site subsets (Z15-translation orbit
    reps of C(15,7)); SHIFT_ANS: for every <=7-popcount site mask the
    (tau, j) with shift(M, tau) subset-of REP7[j].
  * Per rep: row-combination pivot certificates (positions + 120-bit
    W-functionals from the left-null-space of the delta-data map,
    RREF-restricted to the rep's 56 byte-columns), delta-normalized
    kernel generators (120-bit) + d2-preimages (75-bit), and the
    survivor list SURV of coefficient vectors with mincost <= 14.
  * TT/TT_I/TT_C: the sorted 8,475-row translate table of the 113
    `ClassData.REPS` class representatives.

Global completeness assert: the swept survivor reconstructions (weight
<= 14, nonzero) are EXACTLY the TT mask set — i.e. the analytic sweep
re-derives the SAT-enumerated classification on the emitted data.

Regen: cd experiments/bb_lab && uv run python scripts/a22_gen_light_classification.py --force
"""

from __future__ import annotations

import argparse
import os
import re
import sys
from itertools import combinations
from pathlib import Path

import numpy as np

LAB_ROOT = Path(__file__).resolve().parent.parent
REPO_ROOT = LAB_ROOT.parent.parent
QECLEAN_ROOT = Path(os.environ.get("QECLEAN_ROOT", REPO_ROOT / "QECLean"))

# ------------------------------------------------------------ conventions

N5, N15 = 5, 15
NSITE = 15          # sites (i, b): i in Z5, b in Z3, s = 3*i + b
NPOS = 120          # alpha-beta positions: p = 8*s + 4*t + n
NCELL2 = 75         # 2-chain cells g = (x, y): idx = 15*x + y
NCELL1 = 150        # 1-chain cells (g, blk): idx = 75*blk + 15*x + y

A_SUPP = [(0, 0), (0, 1), (1, 0)]      # A = 1 + y + x
B_SUPP = [(1, 6), (1, 10), (2, 12)]    # B = xy^6 + xy^10 + x^2y^12


def cell2(x: int, y: int) -> int:
    return 15 * (x % 5) + (y % 15)


def cell1(x: int, y: int, blk: int) -> int:
    return 75 * blk + 15 * (x % 5) + (y % 15)


def site_iy(s: int) -> tuple[int, int]:
    return s // 3, s % 3


def ucell(s: int, a: int) -> int:
    i, b = site_iy(s)
    return cell1(i, 3 * a + 5 * b, 0)


def vcell(s: int, a: int) -> int:
    i, b = site_iy(s)
    return cell1(i + 1, 3 * a + 5 * b, 1)


def site_add(s: int, tau: int) -> int:
    i, b = site_iy(s)
    ti, tb = site_iy(tau)
    return 3 * ((i + ti) % 5) + ((b + tb) % 3)


def site_sub(s: int, tau: int) -> int:
    i, b = site_iy(s)
    ti, tb = site_iy(tau)
    return 3 * ((i - ti) % 5) + ((b - tb) % 3)


# ------------------------------------------------------------ boundary map

def boundary2(f: np.ndarray) -> np.ndarray:
    """d2 f = (A*f | B*f) as a 150-vector; conv (P*f)(g) = sum_h P(h) f(g-h)."""
    out = np.zeros(NCELL1, dtype=np.int8)
    for x in range(5):
        for y in range(15):
            au = 0
            bv = 0
            for (px, py) in A_SUPP:
                au ^= int(f[cell2(x - px, y - py)])
            for (px, py) in B_SUPP:
                bv ^= int(f[cell2(x - px, y - py)])
            out[cell1(x, y, 0)] = au
            out[cell1(x, y, 1)] = bv
    return out


def fiber_u(bvec: np.ndarray, s: int) -> list[int]:
    return [int(bvec[ucell(s, a)]) for a in range(5)]


def fiber_v(bvec: np.ndarray, s: int) -> list[int]:
    return [int(bvec[vcell(s, a)]) for a in range(5)]


def dnib(t: list[int]) -> int:
    return sum(((t[n] ^ t[4]) << n) for n in range(4))


def eps(t: list[int]) -> int:
    return sum(t) & 1


def abdata(bvec: np.ndarray) -> np.ndarray:
    """The 120-bit alpha-beta data of a 1-chain pair."""
    out = np.zeros(NPOS, dtype=np.int8)
    for s in range(NSITE):
        du = dnib(fiber_u(bvec, s))
        dv = dnib(fiber_v(bvec, s))
        for n in range(4):
            out[8 * s + n] = (du >> n) & 1
            out[8 * s + 4 + n] = (dv >> n) & 1
    return out


def sigma_fiber(h: int, d: int) -> list[int]:
    """The unique 5-bit fiber with eps = h, delta-nibble = d."""
    pard = bin(d).count("1") & 1
    t4 = h ^ pard
    t = [((d >> n) & 1) ^ t4 for n in range(4)] + [t4]
    return t


# W5T[h][d] = weight of sigma_fiber(h, d)
W5T = [[sum(sigma_fiber(h, d)) for d in range(16)] for h in range(2)]


def site_mincost(anib: int, bnib: int) -> int:
    return min(W5T[0][anib] + W5T[0][bnib], W5T[1][anib] + W5T[1][bnib])


# ------------------------------------------------------------ F2 linear algebra

def rref_with_transform(M: np.ndarray) -> tuple[np.ndarray, np.ndarray, list[int]]:
    """RREF over F2. Returns (R, T, pivcols) with R = T @ M (mod 2)."""
    R = M.copy().astype(np.int8)
    rows, cols = R.shape
    T = np.eye(rows, dtype=np.int8)
    piv: list[int] = []
    r = 0
    for c in range(cols):
        pr = None
        for rr in range(r, rows):
            if R[rr, c]:
                pr = rr
                break
        if pr is None:
            continue
        if pr != r:
            R[[r, pr]] = R[[pr, r]]
            T[[r, pr]] = T[[pr, r]]
        for rr in range(rows):
            if rr != r and R[rr, c]:
                R[rr] ^= R[r]
                T[rr] ^= T[r]
        piv.append(c)
        r += 1
        if r == rows:
            break
    return R, T, piv


def nullspace(M: np.ndarray) -> list[np.ndarray]:
    """Basis of {x : Mx = 0} over F2."""
    R, _, piv = rref_with_transform(M)
    cols = M.shape[1]
    free = [c for c in range(cols) if c not in piv]
    basis = []
    for fc in free:
        v = np.zeros(cols, dtype=np.int8)
        v[fc] = 1
        for ri, pc in enumerate(piv):
            v[pc] = R[ri, fc]
        basis.append(v)
    return basis


def solve_f2(M: np.ndarray, b: np.ndarray) -> np.ndarray:
    """One solution of Mx = b over F2 (asserts consistency)."""
    aug = np.concatenate([M, b.reshape(-1, 1)], axis=1).astype(np.int8)
    R, _, piv = rref_with_transform(aug)
    cols = M.shape[1]
    assert cols not in piv, "inconsistent system"
    x = np.zeros(cols, dtype=np.int8)
    for ri, pc in enumerate(piv):
        x[pc] = R[ri, cols]
    assert np.array_equal((M @ x) % 2, b % 2)
    return x


def pack_bits(v: np.ndarray) -> int:
    m = 0
    for p, bit in enumerate(v):
        if bit:
            m |= 1 << p
    return m


# ------------------------------------------------------------ main pipeline

def build_all() -> dict:
    print("== Phi (delta-data of d2) on the 75-point basis ==")
    PHI = np.zeros((NPOS, NCELL2), dtype=np.int8)   # columns = basis images
    DGEN_masks: list[int] = []
    for g in range(NCELL2):
        f = np.zeros(NCELL2, dtype=np.int8)
        f[g] = 1
        y = abdata(boundary2(f))
        PHI[:, g] = y
        DGEN_masks.append(pack_bits(y))
    rank_phi = len(rref_with_transform(PHI)[2])
    assert rank_phi == 56, f"rank Phi = {rank_phi} != 56"
    # global kernel of the delta-data map: dim 75 - 56 = 19??  No:
    # ker(abdata . d2) includes the eps-only directions.  The IMAGE dim is 56.
    print(f"   rank(Phi) = {rank_phi}  (image V'_delta dim = 56)  OK")

    COLPACK = [pack_bits(PHI[p, :]) for p in range(NPOS)]

    print("== left null space of Phi (the W-functional space, dim 64) ==")
    WBASIS = nullspace(PHI.T)   # {w : Phi^T w = 0} = {w : w . every column = 0}
    assert len(WBASIS) == 64, f"dim V'_delta-perp = {len(WBASIS)} != 64"
    WB = np.array(WBASIS, dtype=np.int8)           # 64 x 120
    assert np.all((WB @ PHI) % 2 == 0)

    print("== eps-relation sanity: eps(v-fiber at s+(1,0)) = eps(u-fiber at s) ==")
    rng = np.random.default_rng(220)
    for _ in range(25):
        f = rng.integers(0, 2, NCELL2).astype(np.int8)
        bv = boundary2(f)
        for s in range(NSITE):
            assert eps(fiber_v(bv, s)) == eps(fiber_u(bv, s)), "eps pairing broken"
    print("   25 random f x 15 sites OK")

    print("== weight formula sanity ==")
    for _ in range(25):
        f = rng.integers(0, 2, NCELL2).astype(np.int8)
        bv = boundary2(f)
        w_direct = int(bv.sum())
        w_sites = 0
        for s in range(NSITE):
            tu, tv = fiber_u(bv, s), fiber_v(bv, s)
            hu = eps(tu)
            assert tu == sigma_fiber(hu, dnib(tu))
            assert tv == sigma_fiber(hu, dnib(tv))
            w_sites += W5T[hu][dnib(tu)] + W5T[hu][dnib(tv)]
        assert w_direct == w_sites, "site weight formula broken"
    print("   25 random f OK (sigma-reconstruction + site weights exact)")

    print("== orbit reps of size-7 site subsets ==")
    def shift_mask(m: int, tau: int) -> int:
        out = 0
        for s in range(NSITE):
            if (m >> s) & 1:
                out |= 1 << site_sub(s, tau)
        return out

    def canon(m: int) -> tuple[int, int]:
        best, btau = None, 0
        for tau in range(NSITE):
            c = shift_mask(m, tau)
            if best is None or c < best:
                best, btau = c, tau
        return best, btau

    all7 = [sum(1 << s for s in c) for c in combinations(range(NSITE), 7)]
    reps_set: dict[int, None] = {}
    for m in all7:
        reps_set.setdefault(canon(m)[0], None)
    REP7 = sorted(reps_set.keys())
    assert len(REP7) == 429, f"orbit count {len(REP7)} != 429"
    rep_idx = {m: j for j, m in enumerate(REP7)}
    print(f"   429 orbit reps OK")

    print("== SHIFT_ANS for all popcount<=7 masks ==")
    SHIFT_ANS = np.zeros(1 << NSITE, dtype=np.int64)
    n_small = 0
    for m in range(1 << NSITE):
        if bin(m).count("1") > 7:
            continue
        n_small += 1
        # pad to exactly 7 sites, canonicalize, look up
        mm = m
        s = 0
        while bin(mm).count("1") < 7:
            if not (mm >> s) & 1:
                mm |= 1 << s
            s += 1
        cm, tau = canon(mm)
        j = rep_idx[cm]
        # verify: shift(m, tau) subset of REP7[j]
        assert shift_mask(m, tau) & ~REP7[j] == 0
        SHIFT_ANS[m] = tau + 15 * j
    assert n_small == 16384
    print(f"   {n_small} masks OK")

    print("== per-rep certificates ==")
    piv_pos_flat: list[int] = []
    piv_w_flat: list[int] = []
    sub_piv_off = [0]
    gen_free_flat: list[int] = []
    gen_mask_flat: list[int] = []
    gen_pre_flat: list[int] = []
    sub_gen_off = [0]
    surv_flat: list[int] = []
    sub_surv_off = [0]
    dim_hist: dict[int, int] = {}
    all_recon_masks: set[int] = set()
    n_wtpass = 0

    for j, repm in enumerate(REP7):
        S = [s for s in range(NSITE) if (repm >> s) & 1]
        bytepos = [8 * s + t4n for s in S for t4n in range(8)]
        bytemask_arr = np.zeros(NPOS, dtype=np.int8)
        for p in bytepos:
            bytemask_arr[p] = 1
        outpos = [p for p in range(NPOS) if not bytemask_arr[p]]
        assert len(bytepos) == 56 and len(outpos) == 64

        # solution space Y_S = {y in im Phi : supp(y) in bytepos}
        # via kernel K = {f : Phi f = 0 outside}
        Mout = PHI[outpos, :]                      # 64 x 75
        K = nullspace(Mout)
        img = np.array([(PHI @ v) % 2 for v in K], dtype=np.int8) if K else \
            np.zeros((0, NPOS), dtype=np.int8)
        # reduce image to a basis
        if len(img):
            Rimg, _, pivimg = rref_with_transform(img)
            Ybasis = Rimg[: len(pivimg)]
        else:
            Ybasis = np.zeros((0, NPOS), dtype=np.int8)
        dimY = Ybasis.shape[0]
        assert dimY % 4 == 0 and dimY <= 16, f"rep {j}: dimY = {dimY}"
        dim_hist[dimY] = dim_hist.get(dimY, 0) + 1

        # pivot certificates: RREF of WB restricted to byte columns
        WBS = WB[:, bytepos]                       # 64 x 56
        R, T, pivc = rref_with_transform(WBS)
        rank = len(pivc)
        assert rank == 56 - dimY, f"rep {j}: rank {rank} vs dimY {dimY}"
        piv_positions = [bytepos[c] for c in pivc]
        free_positions = [bytepos[c] for c in range(56) if c not in pivc]
        assert len(free_positions) == dimY
        Wfull = (T @ WB) % 2                       # combinations on full rows
        for ri in range(rank):
            w = Wfull[ri]
            assert np.all((w @ PHI) % 2 == 0)      # still in the left null space
            for ri2, p2 in enumerate(piv_positions):
                assert int(w[p2]) == (1 if ri2 == ri else 0)
            piv_pos_flat.append(piv_positions[ri])
            piv_w_flat.append(pack_bits(w))
        sub_piv_off.append(len(piv_pos_flat))

        # delta-normalized generators of Y_S on the free positions
        gens: list[np.ndarray] = []
        if dimY:
            # solve for gen with gen[free_i] = delta_ij within span(Ybasis)
            F = Ybasis[:, free_positions]          # dimY x dimY, invertible
            Finv_cols = []
            for i in range(dimY):
                e = np.zeros(dimY, dtype=np.int8)
                e[i] = 1
                Finv_cols.append(solve_f2(F.T.astype(np.int8), e))
            for i in range(dimY):
                gen = (Finv_cols[i] @ Ybasis) % 2
                for i2 in range(dimY):
                    assert int(gen[free_positions[i2]]) == (1 if i2 == i else 0)
                assert np.all(gen[outpos] == 0)
                gens.append(gen.astype(np.int8))
        for i, gen in enumerate(gens):
            pre = solve_f2(PHI, gen)               # Phi pre = gen
            gen_free_flat.append(free_positions[i])
            gen_mask_flat.append(pack_bits(gen))
            gen_pre_flat.append(pack_bits(pre))
        sub_gen_off.append(len(gen_free_flat))

        # survivor sweep for this rep
        genmasks = [pack_bits(g) for g in gens]
        survs = []
        for e in range(1, 1 << dimY):
            ab = 0
            for i in range(dimY):
                if (e >> i) & 1:
                    ab ^= genmasks[i]
            mc = 0
            for s in range(NSITE):
                anib = (ab >> (8 * s)) & 15
                bnib = (ab >> (8 * s + 4)) & 15
                if anib or bnib:
                    mc += site_mincost(anib, bnib)
            if mc <= 14:
                survs.append(e)
        surv_flat.extend(survs)
        sub_surv_off.append(len(surv_flat))

        # full (e, hS, out) sweep for the completeness assert
        for e in [0] + survs:
            ab = 0
            for i in range(dimY):
                if (e >> i) & 1:
                    ab ^= genmasks[i]
            nibs = [((ab >> (8 * s)) & 15, (ab >> (8 * s + 4)) & 15)
                    for s in range(NSITE)]
            outsites = [s for s in range(NSITE) if s not in S]
            for hS in range(1 << 7):
                hmask = 0
                for k, s in enumerate(S):
                    if (hS >> k) & 1:
                        hmask |= 1 << s
                wt_in = sum(W5T[(hmask >> s) & 1][nibs[s][0]]
                            + W5T[(hmask >> s) & 1][nibs[s][1]] for s in S)
                base_out = sum(W5T[0][nibs[s][0]] + W5T[0][nibs[s][1]]
                               for s in outsites)
                assert base_out == 0  # nibs vanish outside S
                for out in range(len(outsites) + 1):
                    hm = hmask | (1 << outsites[out - 1] if out else 0)
                    wt = wt_in + (10 if out else 0)
                    if wt > 14 or (e == 0 and hm == 0):
                        continue
                    n_wtpass += 1
                    recon = 0
                    for s in range(NSITE):
                        hb = (hm >> s) & 1
                        tu = sigma_fiber(hb, nibs[s][0])
                        tv = sigma_fiber(hb, nibs[s][1])
                        for a in range(5):
                            if tu[a]:
                                recon |= 1 << ucell(s, a)
                            if tv[a]:
                                recon |= 1 << vcell(s, a)
                    all_recon_masks.add(recon)
        if j % 100 == 0:
            print(f"   rep {j}: dimY={dimY} rank={rank} survs={len(survs)}")

    print(f"   dim histogram: {dim_hist}")
    print(f"   pivots total: {len(piv_pos_flat)}, gens total: {len(gen_free_flat)}, "
          f"survivors total: {len(surv_flat)}, wt-passes: {n_wtpass}, "
          f"distinct recon masks: {len(all_recon_masks)}")

    print("== TT translate table from ClassData.REPS ==")
    classdata = (QECLEAN_ROOT / "QEC/Stabilizer/Codes/BivariateBicycle/"
                 "Z5Z15F2A6/ClassData.lean").read_text()
    m = re.search(r"def REPS : Array Nat :=\s*#\[(.*?)\]", classdata, re.S)
    assert m, "REPS array not found in ClassData.lean"
    REPS = [int(x) for x in re.sub(r"\s", "", m.group(1)).split(",")]
    assert len(REPS) == 113
    tt_rows: list[tuple[int, int, int]] = []
    for i, rm in enumerate(REPS):
        for ccode in range(75):
            c1, c2 = ccode // 15, ccode % 15
            tm = 0
            for blk in range(2):
                for x in range(5):
                    for y in range(15):
                        if (rm >> cell1(x + c1, y + c2, blk)) & 1:
                            tm |= 1 << cell1(x, y, blk)
            tt_rows.append((tm, i, ccode))
    tt_rows.sort()
    TT = [r[0] for r in tt_rows]
    TT_I = [r[1] for r in tt_rows]
    TT_C = [r[2] for r in tt_rows]
    tt_set = set(TT)
    print(f"   8475 rows, {len(tt_set)} distinct masks")

    print("== SOUNDNESS: every swept recon mask is a tabulated translate ==")
    extra = all_recon_masks - tt_set
    assert not extra, f"{len(extra)} recon masks NOT in TT (sweep unsound!)"
    print(f"   {len(all_recon_masks)} distinct recon masks, all in TT")

    print("== ROUTE COMPLETENESS: the Lean assembly path re-derives every "
          "TT mask ==")
    # mirror of the Lean proof: active sites -> SHIFT_ANS -> translate ->
    # must be one of the swept recon masks.
    def active_mask_of(t: int) -> int:
        m = 0
        for s in range(NSITE):
            tu = [(t >> ucell(s, a)) & 1 for a in range(5)]
            tv = [(t >> vcell(s, a)) & 1 for a in range(5)]
            if dnib(tu) or dnib(tv):
                m |= 1 << s
        return m

    def translate_mask(t: int, c1: int, c2: int) -> int:
        out = 0
        for blk in range(2):
            for x in range(5):
                for y in range(15):
                    if (t >> cell1(x + c1, y + c2, blk)) & 1:
                        out |= 1 << cell1(x, y, blk)
        return out

    for t in tt_set:
        m = active_mask_of(t)
        assert bin(m).count("1") <= 7, "TT entry with > 7 active sites?!"
        a = int(SHIFT_ANS[m])
        tau, j = a % 15, a // 15
        tprime = translate_mask(t, tau // 3, 5 * (tau % 3))
        assert active_mask_of(tprime) & ~REP7[j] == 0
        assert tprime in all_recon_masks, \
            f"TT mask not re-derived through rep {j}"
    print(f"   all {len(tt_set)} TT masks re-derived; classification closed")

    return dict(
        COLPACK=COLPACK, DGEN=DGEN_masks, REP7=REP7,
        SHIFT_ANS=[int(x) for x in SHIFT_ANS],
        SUB_PIV_OFF=sub_piv_off, PIV_POS=piv_pos_flat, PIV_W=piv_w_flat,
        SUB_GEN_OFF=sub_gen_off, GEN_FREE=gen_free_flat,
        GEN_MASK=gen_mask_flat, GEN_PRE=gen_pre_flat,
        SUB_SURV_OFF=sub_surv_off, SURV=surv_flat,
        TT=TT, TT_I=TT_I, TT_C=TT_C,
    )


# ------------------------------------------------------------ emission

BANNER = """/-
GENERATED FILE — DO NOT HAND-EDIT.
Generator:
  qec-lab:experiments/bb_lab/scripts/a22_gen_light_classification.py
Data: computed in-script (F2-ified (eps,delta) fibering of the Z5Z15F2A6
delta-data map); every table is numpy-hard-asserted at emission, including
the global completeness cross-check recon-set == translate-table set.
Regen: cd experiments/bb_lab &&
  uv run python scripts/a22_gen_light_classification.py --force
-/
"""

CHUNK = 1024


def _wrap(body: str) -> str:
    lines = []
    cur = ""
    for tok in body.split(","):
        if cur and len(cur) + len(tok) + 1 > 93:
            lines.append(cur + ",")
            cur = tok
        else:
            cur = cur + "," + tok if cur else tok
    lines.append(cur)
    return "\n    ".join(lines)


def _wrap_doc(doc: str) -> str:
    words = doc.split(" ")
    lines: list[str] = []
    cur = "/--"
    for w in words:
        if len(cur) + len(w) + 1 > 92:
            lines.append(cur)
            cur = w
        else:
            cur = cur + " " + w
    lines.append(cur + ". -/" if not cur.endswith("-/") else cur)
    return "\n".join(lines)


def fmt_array(name: str, vals: list[int], doc: str) -> str:
    assert len(vals) <= CHUNK, f"{name}: use fmt_chunked for {len(vals)} entries"
    joined = _wrap(",".join(str(v) for v in vals))
    return f"{_wrap_doc(doc)}\ndef {name} : Array Nat :=\n  #[{joined}]\n"


def fmt_chunked(name: str, getter: str, vals: list[int], doc: str) -> str:
    """Chunked storage: 1024-entry chunk arrays + an Array-of-Arrays
    dispatcher + a total getter (large literals blow elaboration limits)."""
    out = []
    nchunks = (len(vals) + CHUNK - 1) // CHUNK
    for c in range(nchunks):
        joined = _wrap(",".join(str(v) for v in vals[c * CHUNK:(c + 1) * CHUNK]))
        out.append(f"/-- {name} chunk {c}. -/\ndef {name}_c{c} : Array Nat :=\n"
                   f"  #[{joined}]\n")
    refs = _wrap(",".join(f"{name}_c{c}" for c in range(nchunks)))
    out.append(f"/-- {name} chunk table. -/\n"
               f"def {name}_CHUNKS : Array (Array Nat) :=\n  #[{refs}]\n")
    out.append(f"{_wrap_doc(doc)}\ndef {getter} (k : ℕ) : ℕ :=\n"
               f"  ({name}_CHUNKS.getD (k / {CHUNK}) #[]).getD (k % {CHUNK}) 0\n")
    return "\n".join(out)


def emit(data: dict, force: bool) -> None:
    outdir = QECLEAN_ROOT / "QEC/Stabilizer/Codes/BivariateBicycle/Z5Z15F2A6"
    assert outdir.is_dir(), f"missing {outdir} (set QECLEAN_ROOT)"
    f_cert = outdir / "LightCertData.lean"
    f_tt = outdir / "LightTTData.lean"
    for f in (f_cert, f_tt):
        if f.exists() and not force:
            sys.exit(f"{f} exists; pass --force to overwrite")

    ns_open = ("namespace Quantum\nnamespace Stabilizer\nnamespace Homological\n"
               "namespace BB\nnamespace Z5Z15F2A6\n")
    ns_close = ("\nend Z5Z15F2A6\nend BB\nend Homological\nend Stabilizer\n"
                "end Quantum\n")

    opts = "set_option maxRecDepth 4096\n"
    cert = [BANNER,
            "import QEC.Stabilizer.Codes.BivariateBicycle.Z5Z15F2A6.Defs\n",
            "/-!\n# Z5Z15F2A6 light-classification certificate data (A22)\n\n"
            "Pivot certificates, kernel generators, survivor lists, and the\n"
            "orbit/shift tables for the analytic `lightClassification` proof.\n"
            "See the generator docstring for the exact conventions.\n-/\n",
            ns_open, opts]
    cert.append(fmt_array(
        "COLPACK", data["COLPACK"],
        "Transposed delta-data matrix: bit `g` of `COLPACK[p]` is position `p` "
        "of the delta-data of `∂₂ (Pi.single g 1)` (75-bit columns)"))
    cert.append(fmt_array(
        "DGEN", data["DGEN"],
        "Delta-data rows: `DGEN[cell2Idx g]` is the packed 120-bit delta-data "
        "of `∂₂ (Pi.single g 1)`"))
    cert.append(fmt_array(
        "REP7", data["REP7"],
        "The 429 canonical size-7 site subsets (translation-orbit reps, 15-bit)"))
    cert.append(fmt_chunked(
        "SHIFT_ANS", "shiftAnsGet", data["SHIFT_ANS"],
        "Per popcount-≤7 site mask `m`: `τ + 15·j` with `shift(m, τ) ⊆ REP7[j]`"))
    cert.append(fmt_array(
        "SUB_PIV_OFF", data["SUB_PIV_OFF"],
        "Per-rep offsets into the pivot arrays (430 entries)"))
    cert.append(fmt_chunked(
        "PIV_POS", "pivPosGet", data["PIV_POS"],
        "Pivot positions (`< 120`), grouped by rep"))
    cert.append(fmt_chunked(
        "PIV_W", "pivWGet", data["PIV_W"],
        "Pivot W-functionals (120-bit masks in the left null space of the "
        "delta-data map; RREF rows: 1 at own pivot, 0 at the rep's other "
        "pivots)"))
    cert.append(fmt_array(
        "SUB_GEN_OFF", data["SUB_GEN_OFF"],
        "Per-rep offsets into the generator arrays (430 entries)"))
    cert.append(fmt_chunked(
        "GEN_FREE", "genFreeGet", data["GEN_FREE"],
        "Generator free positions (`< 120`), grouped by rep"))
    cert.append(fmt_chunked(
        "GEN_MASK", "genMaskGet", data["GEN_MASK"],
        "Delta-normalized kernel generators (120-bit), grouped by rep"))
    cert.append(fmt_chunked(
        "GEN_PRE", "genPreGet", data["GEN_PRE"],
        "Generator `∂₂`-preimages (75-bit 2-chain masks), grouped by rep"))
    cert.append(fmt_array(
        "SUB_SURV_OFF", data["SUB_SURV_OFF"],
        "Per-rep offsets into the survivor array (430 entries)"))
    cert.append(fmt_chunked(
        "SURV", "survGet", data["SURV"],
        "Survivor coefficient vectors (mincost ≤ 14), grouped by rep"))
    cert.append(ns_close)
    f_cert.write_text("\n".join(cert))
    print(f"wrote {f_cert} ({f_cert.stat().st_size // 1024} KB)")

    tt = [BANNER,
          "import QEC.Stabilizer.Codes.BivariateBicycle.Z5Z15F2A6.ClassData\n",
          "/-!\n# Z5Z15F2A6 light-class translate table (A22)\n\n"
          "The 8,475 packed masks of `translate1 c (repChain i)` over all\n"
          "`i < 113`, `c = (ccode/15, ccode%15)`, sorted ascending for binary\n"
          "search, with aligned decode arrays.\n-/\n",
          ns_open, opts]
    tt.append(fmt_chunked(
        "TT", "ttGet", data["TT"],
        "Sorted packed masks of all 113 × 75 class-rep translates (150-bit)"))
    tt.append(fmt_chunked(
        "TT_I", "ttIGet", data["TT_I"], "Class index `i` of each `TT` row"))
    tt.append(fmt_chunked(
        "TT_C", "ttCGet", data["TT_C"],
        "Translate code `ccode = 15·c₁ + c₂` of each `TT` row"))
    tt.append(ns_close)
    f_tt.write_text("\n".join(tt))
    print(f"wrote {f_tt} ({f_tt.stat().st_size // 1024} KB)")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()
    assert QECLEAN_ROOT.is_dir(), \
        f"QECLean checkout not found at {QECLEAN_ROOT}; set QECLEAN_ROOT"
    data = build_all()
    emit(data, args.force)
    print("ALL ASSERTS PASSED; data emitted.")
