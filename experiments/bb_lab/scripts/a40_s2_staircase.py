#!/usr/bin/env python3
"""A40 S2 / P4 step 1 — witness archaeology: are the minimum-weight
logicals of the four proven family members staircases of translated
weight-6 bb72 blocks over the deck-cell grid?

Cell frame: member (r,b) on Z_{6(r+b)} x Z_{6r} (paper Laurent
presentation) covers bb72 by (x mod 6, y mod 6); deck cells =
(gx div 6, gy div 6) on an (r+b) x r grid; every qubit is
(block, cell, bb72-position).

Witness sources (all verified in-run before analysis):
  bb72   : complete <= 6 coset census (fresh, instant) — 84 expected;
  gross  : complete <= 12 nontrivial census (fresh BZ walk, cheap);
  2-gross: the banked a36 w18 witness (stored presentation), mapped to
           the paper presentation through the session-1-verified move
           (y -> y^7; shifts), then re-verified;
  [[432]]: the session-1 tau0-witness (w24 = TAU of the w12 non-SEAM
           (18,6)-logical from the W22 checkpoint), re-verified; plus
           the 12 w12 (18,6)-logicals themselves on their 3x1 grid.

Output: data/a40/s2_staircase_archaeology.json + stdout analysis.
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
    AxisDeck, TowerCode, batch_keys, i2v, in_span, rep_for,
    translation_perms, v2i, validate_banked,
)
from a38_c37xx_freeze import census_pass, whist  # noqa: E402

DATA = LAB / "data" / "a40"
A_L = [(0, 0), (0, 1), (3, -1)]
B_L = [(0, 0), (1, 0), (-1, -3)]


def red(supp, lm):
    return frozenset((e[0] % lm[0], e[1] % lm[1]) for e in supp)


def member(r, b, name=None):
    lm = (6 * (r + b), 6 * r)
    return TowerCode(name or f"tdg({r},{b})", lm, red(A_L, lm),
                     red(B_L, lm))


def cells_of(code: TowerCode, v: np.ndarray, bb72: TowerCode):
    """Decompose v over deck cells; per-cell restriction as a bb72
    vector (block structure preserved)."""
    lm = code.G.orders
    ng = code.ng
    out: dict[tuple[int, int], np.ndarray] = {}
    for i in np.nonzero(v)[0]:
        blk, gi = divmod(int(i), ng)
        g = code.G.from_index(gi)
        cell = (g[0] // 6, g[1] // 6)
        pos = (g[0] % 6, g[1] % 6)
        rho = out.setdefault(cell, np.zeros(bb72.n, dtype=np.uint8))
        rho[blk * bb72.ng + bb72.G.index(pos)] ^= 1
    return out


def analyze(code, v, bb72, perms72, tag):
    cells = cells_of(code, v, bb72)
    grid = sorted(cells)
    ws = {str(c): int(cells[c].sum()) for c in grid}
    kinds = {}
    keys = {}
    for c in grid:
        rho = cells[c]
        cyc = bb72.is_cycle(rho)
        st = bb72.is_stab(rho) if cyc else False
        kinds[str(c)] = ("stab" if cyc and st else
                         "logical" if cyc and not st else "non-cycle")
        keys[str(c)] = int(batch_keys(rho.reshape(1, -1), perms72)[0][0]) \
            if rho.any() else None
    # translate-equivalence classes among the cell restrictions
    kk = [bytes(batch_keys(cells[c].reshape(1, -1), perms72)[0])
          for c in grid]
    n_classes = len(set(kk))
    # staircase test: cells sorted by (x+y) monotone in both coords?
    xs = [c[0] for c in grid]
    ys = [c[1] for c in grid]
    stair = (len(grid) == len(set(grid)) and
             sorted(grid) == grid and
             all(b2 >= b1 for b1, b2 in zip(ys, ys[1:])) and
             all((x2 - x1, y2 - y1) in [(1, 0), (0, 1)]
                 for (x1, y1), (x2, y2) in zip(grid, grid[1:])))
    res = {"tag": tag, "weight": int(v.sum()), "n_cells": len(grid),
           "cells": [list(c) for c in grid], "cell_weights": ws,
           "cell_kinds": kinds, "translate_classes": n_classes,
           "monotone_staircase": bool(stair)}
    print(f"  {tag}: |v|={res['weight']} cells={res['cells']} "
          f"weights={list(ws.values())} kinds={list(kinds.values())} "
          f"translate-classes={n_classes} staircase={stair}")
    return res


def main():
    t0 = time.time()
    validate_banked(LAB / "data")
    print("validate_banked: PASS")
    out = {}
    bb72 = member(1, 0, "bb72")
    perms72 = translation_perms(bb72)
    binp = cosetbz.build_kernel()

    # ---- bb72 weight-6 logicals (complete <= 6), chunked over classes
    CH = 51
    w6 = []
    allb = sorted(range(1, 1 << bb72.k))
    for lo in range(0, len(allb), CH):
        chunk = allb[lo:lo + CH]
        hits = census_pass(binp, bb72,
                           [(f"C{c}", rep_for(bb72, c)) for c in chunk],
                           6, f"s2_bb72_{lo}")
        for c in chunk:
            for h in sorted(hits[f"C{c}"]):
                v = i2v(h, bb72.n)
                assert bb72.is_cycle(v) and not bb72.is_stab(v)
                if int(v.sum()) == 6:
                    w6.append((c, v))
    print(f"bb72: {len(w6)} weight-6 logicals (banked expectation 84)")
    assert len(w6) == 84
    # orbit structure
    k72 = {bytes(batch_keys(v.reshape(1, -1), perms72)[0])
           for _, v in w6}
    print(f"  {len(k72)} translation-orbit classes of w6 logicals")
    out["bb72"] = {"n_w6": len(w6), "orbits": len(k72)}

    # ---- gross weight-12 logicals (complete <= 12 nontrivial)
    gross = member(1, 1, "gross")
    permsG = translation_perms(gross)
    # CH defined above
    w12 = []
    allc = sorted(range(1, 1 << gross.k))
    for lo in range(0, len(allc), CH):
        chunk = allc[lo:lo + CH]
        hh = census_pass(binp, gross,
                         [(f"C{c}", rep_for(gross, c)) for c in chunk],
                         12, f"s2_gross{lo}")
        for c in chunk:
            for h in sorted(hh[f"C{c}"]):
                v = i2v(h, gross.n)
                assert gross.is_cycle(v) and not gross.is_stab(v)
                w12.append((c, v))
    ws = [int(v.sum()) for _, v in w12]
    assert min(ws) == 12, f"d(gross) != 12?! {min(ws)}"
    print(f"gross: {len(w12)} nontrivial <= 12 (all weight 12: "
          f"{set(ws) == {12}}) — d(gross) = 12 re-derived "
          f"census-complete")
    out["gross"] = {"n_w12": len(w12)}
    print("gross w12 cell analysis (2x1 grid):")
    ga = [analyze(gross, v, bb72, perms72, f"gross-w12[{i}]")
          for i, (c, v) in enumerate(w12[:6])]
    # full statistics over all w12
    stats: dict[str, int] = {}
    for c, v in w12:
        cells = cells_of(gross, v, bb72)
        key = (len(cells),
               tuple(sorted(int(r.sum()) for r in cells.values())))
        stats[str(key)] = stats.get(str(key), 0) + 1
    print(f"  all-{len(w12)} (n_cells, sorted cell weights): {stats}")
    out["gross_stats"] = stats

    # ---- two-gross: map the banked a36 witness into the paper frame
    tg_p = member(2, 0, "two-gross/paper")
    tg_s = TowerCode("two-gross/stored", (12, 12), "x^3 + y^2 + y^7",
                     "y^3 + x + x^2")
    wit = json.loads((LAB / "data" / "a36" /
                      "w18_witness_banked.json").read_text())
    v_s = np.zeros(tg_s.n, dtype=np.uint8)
    v_s[wit["v_support"]] = 1
    assert int(v_s.sum()) == 18
    assert tg_s.is_cycle(v_s) and not tg_s.is_stab(v_s)
    # stored qubit (blk, h) = paper qubit (blk, phi^-1(h + s_blk)),
    # phi: y -> y^7 (an involution on Z12)
    ngt = tg_p.ng
    v_p = np.zeros(tg_p.n, dtype=np.uint8)
    for i in np.nonzero(v_s)[0]:
        blk, gi = divmod(int(i), ngt)
        h = tg_s.G.from_index(gi)
        s = (0, 7) if blk == 0 else (1, 0)
        u = ((h[0] + s[0]) % 12, (7 * (h[1] + s[1])) % 12)
        v_p[blk * ngt + tg_p.G.index(u)] = 1
    if not (tg_p.is_cycle(v_p) and not tg_p.is_stab(v_p)):
        # try the other map direction
        v_p[:] = 0
        for i in np.nonzero(v_s)[0]:
            blk, gi = divmod(int(i), ngt)
            h = tg_s.G.from_index(gi)
            s = (0, 7) if blk == 0 else (1, 0)
            u = ((h[0] - s[0]) % 12, (7 * h[1] - s[1] * 7) % 12)
            u = (u[0] % 12, u[1] % 12)
            v_p[blk * ngt + tg_p.G.index(u)] = 1
    assert int(v_p.sum()) == 18
    assert tg_p.is_cycle(v_p) and not tg_p.is_stab(v_p), \
        "presentation map failed both directions"
    print("two-gross: banked a36 w18 witness mapped to the paper "
          "presentation and re-verified (cycle, non-stab, w18)")
    out["two_gross_witness"] = analyze(tg_p, v_p, bb72, perms72,
                                       "two-gross-w18(a36)")

    # ---- [[432,12]]: the session-1 tau0-witness + its w12 source
    tdg = member(2, 1, "tdg432")
    L1 = member(1, 2, "b2frame")     # (18,6)
    deck0 = AxisDeck(tdg, L1, 1)
    seam_reps = []
    u12s = []
    for line in (DATA / "tdg432" / "ckpt_W22_ntrv1.jsonl").open():
        r = json.loads(line)
        if r["w"] != 12:
            continue
        u = np.zeros(L1.n, dtype=np.uint8)
        u[r["support"]] = 1
        assert L1.is_cycle(u) and not L1.is_stab(u)
        u12s.append(u)
    assert len(u12s) == 12
    print(f"[[432]]: {len(u12s)} w12 (18,6)-logical orbit reps loaded")
    print("(18,6) w12 cell analysis (3x1 grid over bb72):")
    out["b2frame_w12"] = [analyze(L1, u, bb72, perms72, f"(18,6)-w12[{i}]")
                          for i, u in enumerate(u12s)]
    # the tau0 witness for each u12 that is non-SEAM
    from bb_lab.tower import colspace, h1_map, rref_ints, span_points
    Mp0 = h1_map(deck0)
    seam_set = span_points(rref_ints(list(colspace(Mp0)))[0]) - {0}
    wits = []
    for i, u in enumerate(u12s):
        if v2i(L1.sig(u)) in seam_set:
            continue
        v = (deck0.TAU @ u) % 2
        if tdg.is_cycle(v) and not tdg.is_stab(v):
            wits.append((i, v))
    print(f"[[432]]: {len(wits)} verified w24 tau0-witnesses")
    out["tdg432_witness"] = [analyze(tdg, v, bb72, perms72,
                                     f"[[432]]-w24[tau0 of u12[{i}]]")
                             for i, v in wits[:4]]

    out["wall_s"] = round(time.time() - t0, 1)
    (DATA / "s2_staircase_archaeology.json").write_text(
        json.dumps(out, indent=1, default=str))
    print(f"wrote {DATA/'s2_staircase_archaeology.json'} "
          f"({out['wall_s']} s)")


if __name__ == "__main__":
    main()
