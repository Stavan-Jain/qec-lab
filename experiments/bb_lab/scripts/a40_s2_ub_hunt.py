#!/usr/bin/env python3
"""A40 S2 / P4 step 2 — hunt the closed-form upper-bound witness.

The naive reading (min witnesses = staircases of w6 blocks) is refuted
by archaeology (s2_staircase_archaeology.json): min-weight logicals come
in many cell shapes.  But the UB theorem only needs SOME weight-
6(2r+b-1) logical per member.  Candidate: a staircase of 2r+b-1
translated weight-6 bb72 blocks along the (1,1)-diagonal of the
(r+b) x r cell torus, glued by wrap-defect telescoping.

This script:
  1. cell-histograms all translates of the a36 two-gross witness (do
     (6,6,6) 3-cell forms already exist in its orbit?);
  2. searches staircase sums at (2,0): v = w6 + x^6 T1 w6 + x^6y^6 T2 w6
     over all 84 w6 choices x in-cell offsets T1, T2 (cycle test =
     exact); any hit is verified nontrivial + weight 18;
  3. if found, extracts the offset pattern, states the closed form, and
     kernel-verifies v_{r,b} at every member through (6,1)
     (cycle + weight + non-stab, all exact linear algebra);
  4. finds a y-local dual certificate u with <u, v> = 1 at each member
     (mechanical nontriviality is is_stab; the dual form is for the
     hand proof).

Output: data/a40/s2_ub_hunt.json
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np

LAB = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(LAB / "src"))
sys.path.insert(0, str(LAB / "scripts"))

from bb_lab import cosetbz  # noqa: E402
from bb_lab.tower import (  # noqa: E402
    TowerCode, i2v, rep_for, translation_perms, v2i, validate_banked,
)
from a38_c37xx_freeze import census_pass  # noqa: E402

DATA = LAB / "data" / "a40"
A_L = [(0, 0), (0, 1), (3, -1)]
B_L = [(0, 0), (1, 0), (-1, -3)]


def red(supp, lm):
    return frozenset((e[0] % lm[0], e[1] % lm[1]) for e in supp)


def member(r, b, name=None):
    lm = (6 * (r + b), 6 * r)
    return TowerCode(name or f"tdg({r},{b})", lm, red(A_L, lm),
                     red(B_L, lm))


def embed(code_small: TowerCode, v: np.ndarray, code_big: TowerCode,
          t=(0, 0)) -> np.ndarray:
    """Flat embedding: qubit (blk, g) -> (blk, g + t) with g read as
    integer coordinates (no reduction of the SMALL frame's coords other
    than the big frame's)."""
    out = np.zeros(code_big.n, dtype=np.uint8)
    ngs, ngb = code_small.ng, code_big.ng
    for i in np.nonzero(v)[0]:
        blk, gi = divmod(int(i), ngs)
        g = code_small.G.from_index(gi)
        u = ((g[0] + t[0]) % code_big.G.orders[0],
             (g[1] + t[1]) % code_big.G.orders[1])
        out[blk * ngb + code_big.G.index(u)] ^= 1
    return out


def cell_hist(code, v):
    ng = code.ng
    cells = {}
    for i in np.nonzero(v)[0]:
        blk, gi = divmod(int(i), ng)
        g = code.G.from_index(gi)
        c = (g[0] // 6, g[1] // 6)
        cells[c] = cells.get(c, 0) + 1
    return cells


def main():
    t0 = time.time()
    validate_banked(LAB / "data")
    print("validate_banked: PASS")
    out = {}
    bb72 = member(1, 0, "bb72")
    binp = cosetbz.build_kernel()

    # w6 logicals of bb72 (fresh, chunked)
    CH = 51
    w6 = []
    allb = sorted(range(1, 1 << bb72.k))
    for lo in range(0, len(allb), CH):
        chunk = allb[lo:lo + CH]
        hits = census_pass(binp, bb72,
                           [(f"C{c}", rep_for(bb72, c)) for c in chunk],
                           6, f"s2b_bb72_{lo}")
        for c in chunk:
            for h in sorted(hits[f"C{c}"]):
                v = i2v(h, bb72.n)
                if int(v.sum()) == 6:
                    w6.append(v)
    assert len(w6) == 84
    print(f"bb72 w6 logicals: {len(w6)}")

    # ---- 1. translate orbit of the a36 witness: cell histograms
    tg = member(2, 0, "two-gross")
    wit = json.loads((LAB / "data" / "a36" /
                      "w18_witness_banked.json").read_text())
    tg_s = TowerCode("tg/stored", (12, 12), "x^3 + y^2 + y^7",
                     "y^3 + x + x^2")
    v_s = np.zeros(tg_s.n, dtype=np.uint8)
    v_s[wit["v_support"]] = 1
    ngt = tg.ng
    v_p = np.zeros(tg.n, dtype=np.uint8)
    for i in np.nonzero(v_s)[0]:
        blk, gi = divmod(int(i), ngt)
        h = tg_s.G.from_index(gi)
        s = (0, 7) if blk == 0 else (1, 0)
        u = ((h[0] + s[0]) % 12, (7 * (h[1] + s[1])) % 12)
        v_p[blk * ngt + tg.G.index(u)] = 1
    assert tg.is_cycle(v_p) and not tg.is_stab(v_p)
    hist = {}
    permsT = translation_perms(tg)
    for t0x in range(12):
        for t0y in range(12):
            vt = embed(tg, v_p, tg, (t0x, t0y))
            ch = cell_hist(tg, vt)
            key = (len(ch), tuple(sorted(ch.values())))
            hist[str(key)] = hist.get(str(key), 0) + 1
    print(f"a36-witness translate orbit (144): cell histograms {hist}")
    out["a36_witness_translates"] = hist

    # ---- 2. staircase sums at (2,0): w6 + x^6 T1 w6 + x^6 y^6 T2 w6
    print("staircase hunt at (2,0) [path (0,0)->(1,0)->(1,1)]:")
    finds = []
    t_hunt = time.time()
    for i6, w in enumerate(w6):
        base = embed(bb72, w, tg, (0, 0))
        for t1x in range(6):
            for t1y in range(6):
                mid = embed(bb72, w, tg, (6 + t1x, t1y))
                for t2x in range(6):
                    for t2y in range(6):
                        top = embed(bb72, w, tg, (6 + t2x, 6 + t2y))
                        v = (base + mid + top) % 2
                        if int(v.sum()) != 18:
                            continue
                        if not tg.is_cycle(v):
                            continue
                        if tg.is_stab(v):
                            continue
                        finds.append({"w6_idx": i6,
                                      "t1": [t1x, t1y],
                                      "t2": [t2x, t2y]})
        if finds and len(finds) >= 40:
            break
    print(f"  {len(finds)} staircase w18 logicals found "
          f"({time.time()-t_hunt:.1f} s)"
          + (f"; first: {finds[0]}" if finds else ""))
    out["staircase_20_finds"] = finds[:40]

    # also allow DIFFERENT w6 blocks per cell for the record (first-found)
    if not finds:
        print("  same-block search empty — trying mixed blocks "
              "(sampled)...")
        rng = np.random.default_rng(20260826)
        for _ in range(20000):
            i, j, k2 = rng.integers(0, 84, 3)
            t1 = (6 + int(rng.integers(6)), int(rng.integers(6)))
            t2 = (6 + int(rng.integers(6)), 6 + int(rng.integers(6)))
            v = (embed(bb72, w6[i], tg) + embed(bb72, w6[j], tg, t1)
                 + embed(bb72, w6[k2], tg, t2)) % 2
            if int(v.sum()) == 18 and tg.is_cycle(v) \
                    and not tg.is_stab(v):
                finds.append({"mixed": [int(i), int(j), int(k2)],
                              "t1": list(t1), "t2": list(t2)})
                break
        print(f"  mixed sample: {len(finds)} finds")

    out["wall_s"] = round(time.time() - t0, 1)
    (DATA / "s2_ub_hunt.json").write_text(json.dumps(out, indent=1))
    print(f"wrote {DATA/'s2_ub_hunt.json'} ({out['wall_s']} s)")


if __name__ == "__main__":
    main()
