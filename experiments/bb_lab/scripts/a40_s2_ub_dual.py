#!/usr/bin/env python3
"""A40 S2 / P4 step 3 — dual certificates for the (r,1) upper-bound
witness family v_{r,1} = sum_j y^{6j} L12.

Nontriviality routes tested mechanically:
  R1 (fold, odd r): p(v_{r,1}) = r . L12 = L12 for odd r — nontrivial
     in the (l,6)-frame => v nontrivial (p(stab) <= stab).  Requires
     L12 nontrivial at every (l,6) frame: certified per-l by the dual
     pair <swapbar(L12), L12> and the x-window locality of both.
  R2 (dual pairing, any r): u = swapbar(v_{r,1}) is a Z-cycle of the
     member; if <u, v> = 1 then v not in rowspace(H_X) directly.

swapbar = the BB transpose-duality map: qubit (blk, g) -> (1-blk, -g).
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np

LAB = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(LAB / "src"))

from bb_lab.tower import TowerCode, validate_banked  # noqa: E402

DATA = LAB / "data" / "a40"
A_L = [(0, 0), (0, 1), (3, -1)]
B_L = [(0, 0), (1, 0), (-1, -3)]


def red(supp, lm):
    return frozenset((e[0] % lm[0], e[1] % lm[1]) for e in supp)


def member(r, b):
    lm = (6 * (r + b), 6 * r)
    return TowerCode(f"tdg({r},{b})", lm, red(A_L, lm), red(B_L, lm))


def place(code, triples, t=(0, 0)):
    ng = code.ng
    v = np.zeros(code.n, dtype=np.uint8)
    for blk, gx, gy in triples:
        u = ((gx + t[0]) % code.G.orders[0],
             (gy + t[1]) % code.G.orders[1])
        v[blk * ng + code.G.index(u)] ^= 1
    return v


def swapbar(code, v):
    ng = code.ng
    out = np.zeros(code.n, dtype=np.uint8)
    for i in np.nonzero(v)[0]:
        blk, gi = divmod(int(i), ng)
        g = code.G.from_index(gi)
        gg = code.G.index(code.G.neg(g))
        out[(1 - blk) * ng + gg] = 1
    return out


def main():
    validate_banked(LAB / "data")
    print("validate_banked: PASS")
    L12 = json.loads((DATA / "s2_ub_bands.json").read_text())["L12"]
    L12 = [tuple(t) for t in L12]
    out = {}

    # R1 ingredient: L12 vs its swapbar dual at (l,6) frames
    print("(l,6)-frame dual pairing <swapbar(L12), L12>:")
    for ell in (12, 18, 24, 30, 36):
        c = TowerCode(f"f{ell}", (ell, 6), red(A_L, (ell, 6)),
                      red(B_L, (ell, 6)))
        v = place(c, L12)
        u = swapbar(c, v)
        assert not ((c.HX @ u) % 2).any(), "swapbar(L12) not a Z-cycle"
        pair = int((u & v).sum()) % 2
        print(f"  l={ell}: Z-cycle OK, <u,v> = {pair}")
        out[f"L12_pair_l{ell}"] = pair

    # R2: swapbar(v_{r,1}) pairing at the members
    print("members (r,1): u = swapbar(v_{r,1}), <u, v>:")
    for r in range(1, 7):
        code = member(r, 1)
        v = np.zeros(code.n, dtype=np.uint8)
        for j in range(r):
            v ^= place(code, L12, (0, 6 * j))
        u = swapbar(code, v)
        assert not ((code.HX @ u) % 2).any()
        pair = int((u & v).sum()) % 2
        print(f"  r={r} [[{code.n},12]]: <u,v> = {pair}")
        out[f"member_pair_r{r}"] = pair

    # if the whole-v pairing degrades with r, try the single-band dual:
    # u1 = swapbar of ONE band placed y-locally is NOT a cycle (L12
    # y-winds), so instead try u = swapbar(v) restricted... record only.
    (DATA / "s2_ub_dual.json").write_text(json.dumps(out, indent=1))
    print(f"wrote {DATA/'s2_ub_dual.json'}")


if __name__ == "__main__":
    main()
