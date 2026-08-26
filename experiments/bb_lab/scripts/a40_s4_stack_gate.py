#!/usr/bin/env python3
"""A40 S4 gate 1 — quotient-stack phases: cheap nontrivial logicals
from y-deck transfer, and the frame-dependence of any transfer bound.

Facts under test:
  S1: the (18,6) weight-12 minimal logical (L12 reduced) stacks along
      the literal y-deck to (18,12) as a weight-24 NONTRIVIAL logical
      (= the certified tau0-witness of d([[432,12]]) = 24: the b=1
      minimum IS a stack of the rate-2 base phase).
  S2: the a36 two-gross w18 witness (rate 1.5 at (12,12)) stacks to
      (12,24) and (12,36): if the stacks are nontrivial, rate-1.5
      nontrivial y-spanning logicals exist at (12,12k) — so NO
      frame-free y-transfer bound at rate 2 can hold: the transfer
      theorem must be class-aware (cheap phases = quotient stacks,
      killed by base-frame floors or seam arguments).
  S3: the gross (12,6) w12 minimal stacks to (12,12) at weight 24
      > d = 18: stacks need not be minimal at b=0 (the cheap phase at
      (12,12) is the witness pattern itself, not the L12 stack).
Every vector is re-verified end-to-end (cycle, stabilizer membership,
weight, y-gap structure).
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


def max_cyclic_gap(code, v):
    m = code.G.orders[1]
    occ = sorted({code.G.from_index(int(i) % code.ng)[1]
                  for i in np.nonzero(v)[0]})
    if not occ:
        return m
    if len(occ) == m:
        return 0
    return max((occ[(i + 1) % len(occ)] - occ[i]) % m - 1
               for i in range(len(occ)))


def stack(vbase, cbase, ccover):
    """Literal y-deck transfer tau: sum over sheets of the lift."""
    m0 = cbase.G.orders[1]
    m1 = ccover.G.orders[1]
    assert m1 % m0 == 0 and cbase.G.orders[0] == ccover.G.orders[0]
    v = np.zeros(ccover.n, dtype=np.uint8)
    for i in np.nonzero(vbase)[0]:
        blk, gi = divmod(int(i), cbase.ng)
        gx, gy = cbase.G.from_index(gi)
        for j in range(m1 // m0):
            v[blk * ccover.ng + ccover.G.index((gx, gy + m0 * j))] ^= 1
    return v


def report(tag, code, v, expect_nontrivial):
    w = int(v.sum())
    cyc = code.is_cycle(v)
    stab = code.is_stab(v)
    gap = max_cyclic_gap(code, v)
    m = code.G.orders[1]
    assert cyc, f"{tag}: stack is not a cycle"
    nontriv = not stab
    print(f"{tag}: weight {w} (rate {w/m:g}/row), y-spanning "
          f"(max gap {gap}), nontrivial = {nontriv}")
    assert gap <= 3
    assert nontriv == expect_nontrivial, \
        f"{tag}: nontriviality {nontriv} != expected {expect_nontrivial}"
    return dict(tag=tag, lm=list(code.G.orders), weight=w,
                rate=w / m, max_y_gap=int(gap), nontrivial=bool(nontriv))


def main():
    t0 = time.time()
    validate_banked(LAB / "data")
    print("validate_banked: PASS")
    out = {"rows": []}

    L12 = [tuple(t) for t in
           json.loads((DATA / "s2_ub_bands.json").read_text())["L12"]]

    # ---- S1: (18,6) minimal -> (18,12) stack = the tau0-witness ------
    c186 = code_at((18, 6))
    v6 = np.zeros(c186.n, dtype=np.uint8)
    for blk, gx, gy in L12:
        v6[blk * c186.ng + c186.G.index((gx % 18, gy % 6))] ^= 1
    assert c186.is_cycle(v6) and not c186.is_stab(v6) and v6.sum() == 12
    print("base (18,6): L12 weight-12 nontrivial logical re-verified "
          "(d((18,6)) = 12 census-complete, a40 §3.2)")
    c1812 = code_at((18, 12))
    out["rows"].append(report("S1 stack (18,6)->(18,12) x2", c1812,
                              stack(v6, c186, c1812), True))
    print("  -> the certified d = 24 value at [[432,12]] is achieved "
          "by a stack of the rate-2 base phase (consistency PASS)")

    # ---- S2: a36 witness (12,12) -> (12,24), (12,36) stacks ----------
    tg = code_at((12, 12), "two-gross-paper")
    wit = json.loads((LAB / "data" / "a36" /
                      "w18_witness_banked.json").read_text())
    tg_s = TowerCode("tg/stored", (12, 12), "x^3 + y^2 + y^7",
                     "y^3 + x + x^2")
    v_s = np.zeros(tg_s.n, dtype=np.uint8)
    v_s[wit["v_support"]] = 1
    v18 = np.zeros(tg.n, dtype=np.uint8)
    for i in np.nonzero(v_s)[0]:
        blk, gi = divmod(int(i), tg.ng)
        h = tg_s.G.from_index(gi)
        s = (0, 7) if blk == 0 else (1, 0)
        u = ((h[0] + s[0]) % 12, (7 * (h[1] + s[1])) % 12)
        v18[blk * tg.ng + tg.G.index(u)] = 1
    assert tg.is_cycle(v18) and not tg.is_stab(v18) and v18.sum() == 18
    print("base (12,12): a36 w18 witness re-verified in the paper frame")
    for k in (2, 3):
        cK = code_at((12, 12 * k))
        out["rows"].append(report(
            f"S2 stack (12,12)->(12,{12*k}) x{k}", cK,
            stack(v18, tg, cK), True))
    print("  -> rate-1.5 NONTRIVIAL y-spanning logicals exist at "
          "(12,12k): no frame-free rate-2 transfer bound can hold; "
          "the theorem must be class-aware")

    # ---- S3: gross (12,6) w12 -> (12,12) stack -----------------------
    c126 = code_at((12, 6), "gross")
    vg = np.zeros(c126.n, dtype=np.uint8)
    for blk, gx, gy in L12:
        vg[blk * c126.ng + c126.G.index((gx % 12, gy % 6))] ^= 1
    assert c126.is_cycle(vg) and not c126.is_stab(vg) and vg.sum() == 12
    out["rows"].append(report("S3 stack (12,6)->(12,12) x2", tg,
                              stack(vg, c126, tg), True))
    print("  -> nontrivial but weight 24 > d = 18: at b=0 the stack is "
          "NOT minimal (the witness pattern is the cheaper phase)")

    out["wall_s"] = round(time.time() - t0, 1)
    (DATA / "s4_stack_gate.json").write_text(json.dumps(out, indent=1))
    print(f"\nwrote {DATA/'s4_stack_gate.json'} ({out['wall_s']} s)")


if __name__ == "__main__":
    main()
