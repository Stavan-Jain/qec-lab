#!/usr/bin/env python3
"""A40 S4 gate 3 — winding characters: the b-bit exclusion mechanism.

Facts under test:
  W1: the a36 two-gross w18 witness pattern USES both wraps: its naive
      re-placement (same exponents, no stacking) at (24,12) and at
      (12,24) is NOT a cycle — the pattern x-winds and y-winds; it
      exists only at frames with 12 | l (x-side), so its x-stacks pay
      a factor l/12: x-winding cheap phases cost Theta(l) per period
      and cannot threaten the rate-2 floor at large l.
  W2: L12 is x-local: naive re-placement at (l,6) is a cycle for
      l = 12..42 step 6 (re-pin of the S2 uniformity, wider range).
  W3: the b-bit closure arithmetic: a phase with y-period 12 and
      x-period 12 (witness species) tiles a member (6(r+b), 6r) only
      if 12 | 6r and 12 | 6(r+b): for b = 1 this forces r even AND
      r odd — IMPOSSIBLE (pure arithmetic, asserted): the witness
      species never fits any b = 1 member; for b = 0 it requires
      r even, and its 2D stack cost is 18*(r/2)^2 = Theta(r^2) >> 12r.
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


def code_at(lm, name=None):
    return TowerCode(name or f"tdg{lm}", lm, red(A_L, lm), red(B_L, lm))


def place(pts, code):
    v = np.zeros(code.n, dtype=np.uint8)
    for blk, gx, gy in pts:
        v[blk * code.ng
          + code.G.index((gx % code.G.orders[0],
                          gy % code.G.orders[1]))] ^= 1
    return v


def main():
    t0 = time.time()
    validate_banked(LAB / "data")
    print("validate_banked: PASS")
    out = {}

    # witness points in the paper (12,12) frame (transport as in s3/s4)
    tg = code_at((12, 12))
    wit = json.loads((LAB / "data" / "a36" /
                      "w18_witness_banked.json").read_text())
    tg_s = TowerCode("tg/stored", (12, 12), "x^3 + y^2 + y^7",
                     "y^3 + x + x^2")
    v_s = np.zeros(tg_s.n, dtype=np.uint8)
    v_s[wit["v_support"]] = 1
    pts = []
    for i in np.nonzero(v_s)[0]:
        blk, gi = divmod(int(i), tg.ng)
        h = tg_s.G.from_index(gi)
        s = (0, 7) if blk == 0 else (1, 0)
        pts.append((blk, (h[0] + s[0]) % 12, (7 * (h[1] + s[1])) % 12))
    v0 = place(pts, tg)
    assert tg.is_cycle(v0) and not tg.is_stab(v0) and v0.sum() == 18

    # ---- W1: naive re-placements break
    for lm in [(24, 12), (12, 24)]:
        c = code_at(lm)
        v = place(pts, c)
        is_c = c.is_cycle(v)
        print(f"W1 witness naive re-placement at {lm}: cycle = {is_c}")
        assert not is_c, (lm, "witness should wind this axis")
    out["W1"] = "witness x-winds and y-winds (naive re-placements at "\
                "(24,12) and (12,24) are not cycles)"

    # ---- W2: L12 x-locality across the family's l values
    L12 = [tuple(t) for t in
           json.loads((DATA / "s2_ub_bands.json").read_text())["L12"]]
    for l in range(12, 43, 6):
        c = code_at((l, 6))
        v = place(L12, c)
        assert c.is_cycle(v) and not c.is_stab(v) and v.sum() == 12, l
    print("W2 L12 naive placement is a nontrivial w12 logical at (l,6) "
          "for l = 12,18,24,30,36,42")
    out["W2"] = "L12 x-local nontrivial at (l,6), l = 12..42 step 6"

    # ---- W3: the b-bit arithmetic (pure)
    for r in range(1, 200):
        b = 1
        fits = (6 * r) % 12 == 0 and (6 * (r + b)) % 12 == 0
        assert not fits, r
    print("W3 witness-species tiling (12 | m and 12 | l) is impossible "
          "at every b = 1 member r < 200 (and in general: r even and "
          "r odd)")
    out["W3"] = "witness species excluded from all b=1 members; b=0 "\
                "needs r even with Theta(r^2) 2D-stack cost"

    out["wall_s"] = round(time.time() - t0, 1)
    (DATA / "s4_winding_gate.json").write_text(json.dumps(out, indent=1))
    print(f"wrote {DATA/'s4_winding_gate.json'} ({out['wall_s']} s)")


if __name__ == "__main__":
    main()
