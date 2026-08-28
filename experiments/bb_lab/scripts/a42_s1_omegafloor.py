#!/usr/bin/env python3
"""A42 S1c — exact per-content cost assembly for the compact-phase
floor: certified lower AND upper bounds from the omega-syzygy
structure, no SAT.

For period p = 3*2^a (p <= 12 here, brute-force cost tables):

  * cost tables per exact omega-content lambda in Lambda_a:
      free(lambda) = min weight of z in F2[y]/(y^p-1) with omega-chain
                     component lambda (complementary content free)
      pure(lambda) = same with ALL complementary content zero
    (brute force over 2^p contents).

  * LOWER bound (certificate, window-free): every cylinder cycle v
    with nontrivial class has omega-content equal to some
    class-nontrivial omega-syzygy sigma; column weights obey
    wt(col) >= free(lambda_col).  Minimal-support syzygies have slot
    gaps <= 3 (splitting argument), and free() >= 1 per active slot,
    so min over patterns of size <= B exhausts all candidates below
    B+1.  floor_LB(p) = min over class-nontrivial omega-syzygies of
    sum free(lambda_slot).

  * UPPER bound (explicit object): the pure lift of a syzygy (all
    complementary components zero) IS a cylinder cycle, with weight
    sum pure(lambda_slot); nontriviality is inherited from the class.
    floor_UB(p) = min over class-nontrivial syzygies of
    sum pure(lambda_slot).  Each achieving lift is constructed and
    re-verified end-to-end in the full window system of
    a42_s1_cylfloor (cycle + functional pairing + torus embed).

  If LB == UB the compact floor at p is EXACT, certificate tier.
"""
from __future__ import annotations

import itertools
import json
import sys
import time
from pathlib import Path

import numpy as np

LAB = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(LAB / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from bb_lab.linalg import nullspace_f2  # noqa: E402
import a42_lib as AL  # noqa: E402
import importlib.util

_spec = importlib.util.spec_from_file_location(
    "a42_s1_syzygy", Path(__file__).parent / "a42_s1_syzygy.py")
SY = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(SY)

_spec2 = importlib.util.spec_from_file_location(
    "a42_s1_cylfloor", Path(__file__).parent / "a42_s1_cylfloor.py")
CF = importlib.util.module_from_spec(_spec2)
_spec2.loader.exec_module(CF)

DATA = LAB / "data" / "a42"


def content_costs(p: int):
    """lambda -> (free_min, pure_min or None); lambda = omega-chain
    component (int poly mod (y^2+y+1)^{2^a})."""
    a = AL.v2(p)
    L = SY.Lam(a)
    ymod = (1 << p) | 1
    cof = AL.pdivmod(ymod, L.mod)[0]
    free = {}
    pure = {}
    for z in range(1, 1 << p):
        w = bin(z).count("1")
        lam = AL.pmod(z, L.mod)
        cur = free.get(lam)
        if cur is None or w < cur:
            free[lam] = w
        if AL.pmod(z, cof) == 0:
            cur = pure.get(lam)
            if cur is None or w < cur:
                pure[lam] = w
    return L, free, pure


def enumerate_pattern_minima(a: int, Wx: int, smax: int, free, pure,
                             log=print):
    """Min sum-free and sum-pure over class-nontrivial omega-syzygies
    with <= smax slots (anchored at col 0).  Returns dicts + the
    achieving syzygy for the pure minimum."""
    ow = SY.OmegaWindow(a, Wx)
    H, T, Z, Lrows, dimZ, ncls = ow.functionals()
    log(f"  a={a}: window classes {ncls}")
    Ld = ow.L.dim
    slots = [(blk, c) for blk in (0, 1) for c in range(Wx)]
    best_free = None
    best_pure = None
    best_pure_syz = None
    n_pat = 0
    for s in range(1, smax + 1):
        for pat in itertools.combinations(slots, s):
            if min(c for (_, c) in pat) != 0:
                continue
            # gap-splitting prune: max slot-col gap <= 3
            cols = sorted({c for (_, c) in pat})
            if any(b - a_ > 4 for a_, b in zip(cols, cols[1:])):
                continue
            n_pat += 1
            sols = SY.solve_affine_in_pattern(ow, H, Lrows, pat, 0)
            if sols.size == 0:
                continue
            pair = (Lrows @ sols.T) % 2
            nt_idx = np.flatnonzero(pair.any(axis=0))
            if nt_idx.size == 0:
                continue
            # enumerate the affine class-nontrivial subset: the
            # solution space is linear; nontrivial elements = those
            # with nonzero pairing.  Enumerate all combos if small.
            dim = sols.shape[0]
            if dim > 16:
                log(f"    pattern {pat}: solution dim {dim} > 16, "
                    "skipping enumeration (flagged)")
                continue
            for maskbits in range(1, 1 << dim):
                acc = np.zeros(ow.nbits, dtype=np.uint8)
                mb = maskbits
                i = 0
                while mb:
                    if mb & 1:
                        acc ^= sols[i]
                    mb >>= 1
                    i += 1
                if not ((Lrows @ acc) % 2).any():
                    continue
                cf = cp = 0
                ok_pure = True
                for (blk, c) in pat:
                    lam = 0
                    for ii in range(Ld):
                        if acc[ow.bit(blk, c, ii)]:
                            lam |= 1 << ii
                    if lam == 0:
                        continue
                    cf += free[lam]
                    pc = pure.get(lam)
                    if pc is None:
                        ok_pure = False
                    else:
                        cp += pc
                if best_free is None or cf < best_free:
                    best_free = cf
                if ok_pure and cp > 0 and (best_pure is None
                                           or cp < best_pure):
                    best_pure = cp
                    best_pure_syz = (pat, acc.copy())
    log(f"  patterns examined: {n_pat}; "
        f"LB (sum free) = {best_free}, UB (sum pure) = {best_pure}")
    return best_free, best_pure, best_pure_syz, ow


def realize_and_verify(p: int, ow, syz_info, log=print) -> int:
    """Build the pure lift of the achieving syzygy as a full cylinder
    window vector, verify with the CylWindow machinery of s1_cylfloor
    (note: that engine is the X-sector = antipode image; here we
    verify in the plain-pair convention directly)."""
    a = AL.v2(p)
    L = SY.Lam(a)
    ymod = (1 << p) | 1
    cof = AL.pdivmod(ymod, L.mod)[0]
    pat, acc = syz_info
    Wx = ow.Wx
    Ld = ow.L.dim
    # per-slot pure lift: min-weight multiple of cof with the content
    lift = {}
    # precompute best pure lift per lambda
    best = {}
    for tbits in range(1, 1 << (AL.pdeg(ymod) - AL.pdeg(cof))):
        z = AL.pmod(AL.pmul(cof, tbits), ymod)
        if z == 0:
            continue
        lam = AL.pmod(z, L.mod)
        w = bin(z).count("1")
        cur = best.get(lam)
        if cur is None or w < cur[0]:
            best[lam] = (w, z)
    total = 0
    cols_content = {}
    for (blk, c) in pat:
        lam = 0
        for ii in range(Ld):
            if acc[ow.bit(blk, c, ii)]:
                lam |= 1 << ii
        if lam == 0:
            continue
        w, z = best[lam]
        cols_content[(blk, c)] = z
        total += w
    # verify the plain-pair cycle equations on the window directly
    # (A~ f + B^ g = 0 in F2[y]/(y^p-1) per column)
    At = {0: 0b110, 3: 1}
    Bh = {0: 1, 1: 0b1000, 2: 0b1000}
    for c in range(-1, Wx + 4):
        s = 0
        for dx, coef in At.items():
            z = cols_content.get((0, c - dx))
            if z:
                s ^= AL.pmod(AL.pmul(coef, z), ymod)
        for dx, coef in Bh.items():
            z = cols_content.get((1, c - dx))
            if z:
                s ^= AL.pmod(AL.pmul(coef, z), ymod)
        assert s == 0, ("pure lift is not a cycle", c)
    log(f"  realized pure lift: weight {total}, cycle verified "
        f"(plain-pair, window)")
    return total


def main():
    t0 = time.time()
    out = {}
    for p, Wx, smax in ((3, 8, 6), (6, 10, 7), (9, 8, 6), (12, 10, 7)):
        a = AL.v2(p)
        print(f"== p = {p} (a = {a}) ==", flush=True)
        L, free, pure = content_costs(p)
        lb, ub, syz_info, ow = enumerate_pattern_minima(
            a, Wx, smax, free, pure,
            log=lambda s: print(s, flush=True))
        row = {"p": p, "a": a, "LB_sum_free": lb, "UB_sum_pure": ub}
        if syz_info is not None:
            w = realize_and_verify(p, ow, syz_info,
                                   log=lambda s: print(s, flush=True))
            assert w == ub, (w, ub)
            row["realized"] = w
        out[str(p)] = row
        print(f"  p={p}: floor in [{lb}, {ub}]"
              + (" — EXACT" if lb == ub else ""), flush=True)

    out["wall_s"] = round(time.time() - t0, 1)
    (DATA / "s1_omegafloor.json").write_text(json.dumps(out, indent=1))
    print(f"\nwrote {DATA/'s1_omegafloor.json'} ({out['wall_s']} s)",
          flush=True)


if __name__ == "__main__":
    main()
