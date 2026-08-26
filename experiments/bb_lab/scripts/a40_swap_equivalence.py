#!/usr/bin/env python3
"""A40 P1 follow-up — are the swapped-frame literal quotients the smaller
family members?

The lattice pass found three member pairs where a group quotient exists but
no axis-aligned literal lift does:

    (2,1) -> (1,1)  deg 3   [(18,12)/<x^6>  = (6,12)-code,  swap of (12,6)]
    (3,1) -> (2,1)  deg 2   [(24,18)/<x^12> = (12,18)-code, swap of (18,12)]
    (4,1) -> (1,1)  deg 10

Each hinges on ONE question: is the family swap-symmetric as CODES — i.e.
is the literal quotient on the swapped frame (m', l') equivalent to the
member code on (l', m') under the standard presentation-move group?

Move set searched (the repo's anchorable-presentation currency):
    Aut(Z_l x Z_m)  (all bijective monomial substitutions x -> x^a y^c,
                     y -> x^d y^b, enumerated and verified bijective)
  x per-block monomial shifts (A -> g*A, B -> h*B, normalized away)
  x block swap (A,B) -> (B,A)          [L/R qubit block relabel]
  x transpose duality (A,B) -> (Abar,Bbar) via inversion (in Aut already)
  x X<->Z swap (A,B) -> (Bbar,Abar)    [CSS Hadamard swap]

A found move = a verified equivalence (checked end-to-end on the parity
matrices).  Not found = no equivalence under this move set (recorded as
such; NOT a proof of inequivalence).
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from bb_lab.tower import TowerCode, support_str  # noqa: E402

DATA = ROOT / "data" / "a40"
A_L = [(0, 0), (0, 1), (3, -1)]
B_L = [(0, 0), (1, 0), (-1, -3)]


def red(supp, lm):
    return frozenset((e[0] % lm[0], e[1] % lm[1]) for e in supp)


def swap_supp(supp):
    return frozenset((e[1], e[0]) for e in supp)


def shift_normal(supp, lm):
    """Lex-min monomial-shift normal form of a support."""
    best = None
    for t in supp:
        cand = tuple(sorted(((e[0] - t[0]) % lm[0], (e[1] - t[1]) % lm[1])
                            for e in supp))
        if best is None or cand < best:
            best = cand
    return best


def all_auts(lm):
    """All automorphisms of Z_l x Z_m as (image of x, image of y) pairs,
    verified bijective by image-set cardinality."""
    l, m = lm
    els = [(i, j) for i in range(l) for j in range(m)]
    auts = []
    for ux in els:
        # order of ux must be a multiple of ... just check hom+bijective
        for uy in els:
            # hom well-defined: l*ux = 0 automatically; m*uy = 0 needs
            # (m*uy) == 0 mod (l,m)
            if (m * uy[0]) % l or (m * uy[1]) % m:
                continue
            if (l * ux[0]) % l or (l * ux[1]) % m:
                continue
            # bijectivity: the image of the generator grid spans; test via
            # matrix determinant-free closure count on a coarse grid
            seen = np.zeros((l, m), dtype=bool)
            # image of (s,t) = s*ux + t*uy
            xs = (np.outer(np.arange(l), [ux[0]]) % l)
            # vectorized: e_(s,t) -> (s*ux0 + t*uy0 mod l, s*ux1 + t*uy1 mod m)
            S = np.arange(l).reshape(-1, 1)
            T = np.arange(m).reshape(1, -1)
            gx = (S * ux[0] + T * uy[0]) % l
            gy = (S * ux[1] + T * uy[1]) % m
            seen[gx.ravel(), gy.ravel()] = True
            if seen.all():
                auts.append((ux, uy))
    return auts


def apply_aut(supp, aut, lm):
    ux, uy = aut
    l, m = lm
    return frozenset((((e[0] * ux[0] + e[1] * uy[0]) % l),
                      ((e[0] * ux[1] + e[1] * uy[1]) % m)) for e in supp)


def bar(supp, lm):
    return frozenset(((-e[0]) % lm[0], (-e[1]) % lm[1]) for e in supp)


def find_equivalence(pairA, pairB, lm, auts):
    """Search for a move taking pairA=(As,Bs) to pairB under the move set.
    Returns (variant, aut) or None."""
    tgtA = shift_normal(pairB[0], lm)
    tgtB = shift_normal(pairB[1], lm)
    variants = {
        "identity": lambda A, B: (A, B),
        "block_swap": lambda A, B: (B, A),
        "xz_swap": lambda A, B: (bar(B, lm), bar(A, lm)),
        "xz_swap+block": lambda A, B: (bar(A, lm), bar(B, lm)),
    }
    for vname, vf in variants.items():
        A0, B0 = vf(*pairA)
        for aut in auts:
            A1 = apply_aut(A0, aut, lm)
            if shift_normal(A1, lm) != tgtA:
                continue
            B1 = apply_aut(B0, aut, lm)
            if shift_normal(B1, lm) == tgtB:
                return vname, aut
    return None


def verify_move(pairA, pairB, lm, hit):
    """End-to-end: build both codes, apply the move as an explicit qubit
    permutation, assert rowspace equality of the parity matrices."""
    vname, aut = hit
    A1, B1 = pairA
    if vname in ("block_swap",):
        A1, B1 = B1, A1
    elif vname == "xz_swap":
        A1, B1 = bar(B1, lm), bar(A1, lm)
    elif vname == "xz_swap+block":
        A1, B1 = bar(A1, lm), bar(B1, lm)
    A1 = apply_aut(A1, aut, lm)
    B1 = apply_aut(B1, aut, lm)
    # now A1,B1 should equal pairB up to per-block shifts
    sA = None
    for t in [(0, 0)] + [tuple(np.subtract(b, a) % lm) for a in A1
                         for b in pairB[0]]:
        if frozenset(((e[0] + t[0]) % lm[0], (e[1] + t[1]) % lm[1])
                     for e in A1) == pairB[0]:
            sA = t
            break
    sB = None
    for t in [(0, 0)] + [tuple(np.subtract(b, a) % lm) for a in B1
                         for b in pairB[1]]:
        if frozenset(((e[0] + t[0]) % lm[0], (e[1] + t[1]) % lm[1])
                     for e in B1) == pairB[1]:
            sB = t
            break
    assert sA is not None and sB is not None
    c1 = TowerCode("m1", lm, A1, B1)          # pre-shift
    c2 = TowerCode("m2", lm, pairB[0], pairB[1])
    # per-block monomial shift = qubit relabeling inside each block:
    # M_{A+s}[g, h] = A((g-h)-s) = M_A[g, h+s], so column h of c2's block
    # is column (h + s_blk) of c1's block.
    ng = c1.ng
    cols = np.zeros(c1.n, dtype=int)
    for i, g in enumerate(c1.G):
        gA = c1.G.index(c1.G.reduce((g[0] + sA[0], g[1] + sA[1])))
        gB = c1.G.index(c1.G.reduce((g[0] + sB[0], g[1] + sB[1])))
        cols[i] = gA
        cols[ng + i] = ng + gB
    HX1 = c1.HX[:, cols]
    # rowspace equality (rref basis is canonical as a SET, insertion-ordered)
    from bb_lab.tower import rref_ints, v2i
    b1, _ = rref_ints([v2i(r) for r in HX1])
    b2, _ = rref_ints([v2i(r) for r in c2.HX])
    assert sorted(b1) == sorted(b2), "X rowspaces differ after the move"
    HZ1 = c1.HZ[:, cols]
    z1, _ = rref_ints([v2i(r) for r in HZ1])
    z2, _ = rref_ints([v2i(r) for r in c2.HZ])
    assert sorted(z1) == sorted(z2), "Z rowspaces differ after the move"
    return {"variant": vname, "aut_x": list(aut[0]), "aut_y": list(aut[1]),
            "shift_A": list(sA), "shift_B": list(sB),
            "verified": "HX and HZ rowspaces equal after explicit qubit "
                        "permutation"}


def _pyify(o):
    if isinstance(o, dict):
        return {k: _pyify(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [_pyify(v) for v in o]
    if isinstance(o, (np.integer,)):
        return int(o)
    return o


def iso_frames(big_lm, small_lm):
    """All divisor frames (a, c) of big_lm with Z_a x Z_c iso
    Z_{small_lm} — the candidate literal-quotient frames.  Isomorphism by
    per-prime invariant comparison (has_quotient both ways = iso here
    since orders are equal)."""
    from math import prod

    def divs(n):
        return [d for d in range(1, n + 1) if n % d == 0]

    tgt = small_lm[0] * small_lm[1]
    outs = []
    for a in divs(big_lm[0]):
        for c in divs(big_lm[1]):
            if a * c != tgt:
                continue
            # iso test: same per-prime sorted exponent pairs
            import math as _m

            def pparts(u, v):
                pp = {}
                for p in {q for n in (u, v) for q in _pr(n)}:
                    pp[p] = tuple(sorted((_vp(u, p), _vp(v, p)),
                                         reverse=True))
                return pp
            if pparts(a, c) == pparts(*small_lm):
                outs.append((a, c))
    return outs


def _pr(n):
    s, d = set(), 2
    while d * d <= n:
        while n % d == 0:
            s.add(d)
            n //= d
        d += 1
    if n > 1:
        s.add(n)
    return s


def _vp(n, p):
    v = 0
    while n % p == 0:
        n //= p
        v += 1
    return v


def group_iso(src_lm, dst_lm):
    """One explicit isomorphism Z_{src} -> Z_{dst} as generator images,
    or None.  Brute-force over hom images, verified bijective."""
    l1, m1 = src_lm
    l2, m2 = dst_lm
    els = [(i, j) for i in range(l2) for j in range(m2)]
    for ux in els:
        if (l1 * ux[0]) % l2 or (l1 * ux[1]) % m2:
            continue
        for uy in els:
            if (m1 * uy[0]) % l2 or (m1 * uy[1]) % m2:
                continue
            S = np.arange(l1).reshape(-1, 1)
            T = np.arange(m1).reshape(1, -1)
            gx = (S * ux[0] + T * uy[0]) % l2
            gy = (S * ux[1] + T * uy[1]) % m2
            seen = np.zeros((l2, m2), bool)
            seen[gx.ravel(), gy.ravel()] = True
            if seen.all():
                return (ux, uy)
    return None


def transport(supp, iso, dst_lm):
    ux, uy = iso
    return frozenset((((e[0] * ux[0] + e[1] * uy[0]) % dst_lm[0]),
                      ((e[0] * ux[1] + e[1] * uy[1]) % dst_lm[1]))
                     for e in supp)


def main():
    t0 = time.time()
    out = {}
    # The three member pairs the lattice pass left open (group quotient
    # exists, no axis-aligned literal lift): resolve via ALL isomorphic
    # divisor frames + the standard monomial move set.
    CASES = [
        ("(2,1)->(1,1)", (18, 12), (12, 6)),
        ("(3,1)->(2,1)", (24, 18), (18, 12)),
        ("(4,1)->(1,1)", (30, 24), (12, 6)),
    ]
    aut_cache: dict = {}
    for tag, big_lm, small_lm in CASES:
        frames = iso_frames(big_lm, small_lm)
        rows = []
        verdict = None
        for qlm in frames:
            Aq = red(A_L, qlm)
            Bq = red(B_L, qlm)
            assert len(Aq) == 3 and len(Bq) == 3, "support collision"
            iso = group_iso(qlm, small_lm)
            assert iso is not None
            At = transport(Aq, iso, small_lm)
            Bt = transport(Bq, iso, small_lm)
            Am = red(A_L, small_lm)
            Bm = red(B_L, small_lm)
            if small_lm not in aut_cache:
                aut_cache[small_lm] = all_auts(small_lm)
            auts = aut_cache[small_lm]
            hit = find_equivalence((At, Bt), (Am, Bm), small_lm, auts)
            row = {"quotient_frame": list(qlm),
                   "iso_gen_images": [list(iso[0]), list(iso[1])],
                   "transported_pair": [support_str(At), support_str(Bt)],
                   "member_pair": [support_str(Am), support_str(Bm)],
                   "n_auts": len(auts)}
            if hit:
                row["equivalence"] = verify_move((At, Bt), (Am, Bm),
                                                 small_lm, hit)
                verdict = {"frame": list(qlm), **row["equivalence"]}
                print(f"{tag} via frame {qlm}: EQUIVALENT — "
                      f"{row['equivalence']}")
            else:
                row["equivalence"] = None
                print(f"{tag} via frame {qlm}: no equivalence "
                      f"({len(auts)} auts x 4 variants x shifts)")
            rows.append(row)
        out[tag] = {
            "frames_tried": [list(f) for f in frames],
            "rows": rows,
            "verdict": ("COVER (literal quotient equivalent to the member "
                        "under monomial moves)" if verdict else
                        "NOT a cover under literal quotients + the "
                        "standard monomial move set (twisted-descent "
                        "question left open)"),
            "witness_move": verdict,
        }
        print(f"  => {tag}: {out[tag]['verdict']}")
    out["wall_s"] = round(time.time() - t0, 1)
    (DATA / "swap_equivalence.json").write_text(
        json.dumps(_pyify(out), indent=1))
    print(f"wrote {DATA/'swap_equivalence.json'} ({out['wall_s']} s)")


if __name__ == "__main__":
    main()
