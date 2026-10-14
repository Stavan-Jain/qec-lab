#!/usr/bin/env python3
"""A42 S2c — cheap fills: new-period machine checks of Theorem A,
Theorem H, and the parity lemma's per-p scope.

* Barren periods p in {10, 11, 13, 14} (3 nmid p, 127 nmid p):
  Theorem A predicts H = 0.  Machine check: dim Z_W = dim B_W on
  window systems at widths W in {6, 9, 12} — no window-supported
  nontrivial compact cycle exists (exact linear algebra, certificate
  for extent <= 12 at ANY weight; Theorem A covers all extents).

* omega-periods p = 15 (= 3*5, a = 0) and p = 18 (= 3*3*2, a = 1):
  Theorem H predicts dim H = 4 (a = 0: F4^2) and 6 (a >= 1:
  F4 + F4[pi]/pi^2), independent of the odd dilation m.  Machine
  check: dim Z_W - dim B_W == prediction, stable at two widths.
  p = 18 is the r = 3 member's own compact period (m = 6r = 18) —
  the first machine contact with the r = 3 row of the b = 1 column.

* Parity scope at every p above: every qubit column of the (slack-
  complete) window check system has weight exactly 3 (odd), which is
  the entire hypothesis of the parity lemma (sum of all checks = the
  all-ones functional => compact cycles have even weight).
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

LAB = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(LAB / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import importlib.util

_spec = importlib.util.spec_from_file_location(
    "a42_s1_cylfloor", Path(__file__).parent / "a42_s1_cylfloor.py")
CF = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(CF)

DATA = LAB / "data" / "a42"

EXPECT = {10: 0, 11: 0, 13: 0, 14: 0,   # barren (Theorem A)
          15: 4,                        # a = 0, m = 5   (Theorem H)
          18: 6}                        # a = 1, m = 3   (Theorem H)


def main():
    t0 = time.time()
    out = {}
    ok = True
    for p, expect in EXPECT.items():
        row = {"expect_dimH": expect, "widths": {}}
        for W in (6, 9, 12):
            cw = CF.CylWindow(p, W)
            H = cw.build_H()
            colw = H.sum(axis=0)
            assert (colw == 3).all(), (p, W, "column weight != 3")
            Bnd = cw.build_boundary_in_window()
            L, dimZ, dimB = cw.class_functionals(H, Bnd)
            got = dimZ - dimB
            row["widths"][str(W)] = {"dimZ": int(dimZ),
                                     "dimB": int(dimB),
                                     "dimH": int(got)}
            tag = "OK" if got == expect else "MISMATCH"
            if got != expect:
                ok = False
            print(f"p={p:2d} W={W:2d}: dimZ={dimZ:3d} dimB={dimB:3d} "
                  f"dimH={got} (expect {expect}) {tag}; "
                  f"parity columns all weight 3", flush=True)
        out[str(p)] = row
    out["all_match"] = ok
    out["wall_s"] = round(time.time() - t0, 1)
    (DATA / "s2_fills.json").write_text(json.dumps(out, indent=1))
    print(f"\nall_match={ok}; wrote {DATA/'s2_fills.json'} "
          f"({out['wall_s']} s)", flush=True)
    assert ok


if __name__ == "__main__":
    main()
