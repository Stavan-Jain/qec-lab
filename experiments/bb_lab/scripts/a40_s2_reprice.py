#!/usr/bin/env python3
"""A40 S2 / P5 — corrected pricing for (3,0) and (3,1) floors under the
v2 direct-primary architecture (session-1's [[432]] execution pattern):
pick the deepest BZ-able level (n <= 192) as the direct-census anchor,
composite H1 ranks to shrink the seam-shadow fan-out, kernel-shift for
deep caps.  COST verdicts only.

Also prices the b=1 members' PARTIAL floor questions now that their
upper bounds are witness-banked (s2_ub_bands): the gap each partial
floor would close.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

LAB = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(LAB / "src"))

from bb_lab.tower import (  # noqa: E402
    AxisDeck, TowerCode, census_nodes, colspace, gf2_rank, h1_map,
    v2i, validate_banked,
)

DATA = LAB / "data" / "a40"
A_L = [(0, 0), (0, 1), (3, -1)]
B_L = [(0, 0), (1, 0), (-1, -3)]


def red(supp, lm):
    return frozenset((e[0] % lm[0], e[1] % lm[1]) for e in supp)


def main():
    validate_banked(LAB / "data")
    print("validate_banked: PASS")
    out = {}

    # ---------- (3,0) tower (18,18) -> (9,18) -> (9,9)
    L0 = TowerCode("t30", (18, 18), red(A_L, (18, 18)),
                   red(B_L, (18, 18)))
    L1 = TowerCode("t30m", (9, 18), red(A_L, (9, 18)),
                   red(B_L, (9, 18)))
    L2 = TowerCode("t30b", (9, 9), red(A_L, (9, 9)), red(B_L, (9, 9)))
    d0 = AxisDeck(L0, L1, 0)
    d1 = AxisDeck(L1, L2, 1)
    Mp0, Mp1 = h1_map(d0), h1_map(d1)
    comp = (Mp1 @ Mp0) % 2
    rc = gf2_rank([v2i(c) for c in comp.T])
    print(f"(3,0): k chain {L0.k}/{L1.k}/{L2.k}; kappa(L2) = {L2.kappa}; "
          f"rank p0* = {gf2_rank([v2i(c) for c in Mp0.T])}, "
          f"rank p1* = {gf2_rank([v2i(c) for c in Mp1.T])}, "
          f"rank composite = {rc} -> S1' classes = {2**rc - 1}")
    rows = []
    for W in (14, 16, 18, 20, 22):
        n_s1p = 2 ** rc - 1
        nodes = (census_nodes(L2.kappa, W)              # stab
                 + census_nodes(L2.kappa, W, n_s1p)     # S1' cosets
                 + census_nodes(L2.kappa, W - 4, 2 ** L1.k // 2 ** 6)
                 )
        cap = (W - 6) // 2
        rows.append({"W": W, "floor": W + 2,
                     "log10_bottom_nodes": round(np.log10(nodes), 1),
                     "cap": cap})
        print(f"  W={W} (d >= {W+2}): bottom walks ~1e"
              f"{np.log10(nodes):.1f} nodes, cap {cap}")
    out["t30"] = {"k": [L0.k, L1.k, L2.k], "rank_composite": rc,
                  "rows": rows}

    # ---------- (3,1) tower (24,18) -> (12,18) -> (12,9) -> (6,9) -> (3,9)
    chain = [(24, 18), (12, 18), (12, 9), (6, 9), (3, 9)]
    codes = [TowerCode(f"t31_{i}", lm, red(A_L, lm), red(B_L, lm))
             for i, lm in enumerate(chain)]
    ks = [c.k for c in codes]
    print(f"(3,1): chain {chain} k {ks} "
          f"kappas {[c.kappa for c in codes]}")
    # composite rank from the top to each level
    decks = []
    for i in range(len(chain) - 1):
        ax = 0 if chain[i][0] != chain[i + 1][0] else 1
        decks.append(AxisDeck(codes[i], codes[i + 1], ax))
    M = None
    comps = []
    for i, dk in enumerate(decks):
        Mi = h1_map(dk)
        M = Mi if M is None else (Mi @ M) % 2
        comps.append(gf2_rank([v2i(c) for c in M.T]))
    print(f"  composite ranks top->level: {comps}")
    rows = []
    for W in (16, 22):
        # BZ-able levels: n <= 192: (6,9) n=108, (3,9) n=54; the direct
        # anchor = (6,9); levels above need stacked descent (2 layers).
        nodes = (census_nodes(codes[3].kappa, W)
                 + census_nodes(codes[3].kappa, W, 2 ** comps[2] - 1))
        rows.append({"W": W, "floor": W + 2,
                     "log10_anchor_nodes": round(np.log10(nodes), 1),
                     "descent_layers_above_anchor": 2,
                     "cap": (W - 6) // 2})
        print(f"  W={W} (d >= {W+2}): anchor (6,9) walks ~1e"
              f"{np.log10(nodes):.1f}, TWO stacked descent layers above "
              f"(n=216, 432), cap {(W-6)//2}")
    out["t31"] = {"k": ks, "comps": comps, "rows": rows}

    (DATA / "s2_reprice.json").write_text(json.dumps(out, indent=1))
    print(f"wrote {DATA/'s2_reprice.json'}")


if __name__ == "__main__":
    main()
