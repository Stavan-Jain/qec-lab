#!/usr/bin/env python3
"""A40 S3 — falsify-first gates for the L1 floor theorem.

Candidate theorem (the gap-dichotomy + plane-Koszul argument):

  For any BB code on Z_l x Z_m with polynomials A, B whose plane pair
  is REGULAR (equivalently: the Koszul homology H1 over
  F2[x^+-, y^+-] vanishes; certified by a nonzero resultant), every
  nontrivial X-logical v satisfies
      |v| >= min( ceil(m / s_y), ceil(l / s_x) )
  where s_x, s_y are the per-axis check spans (both = 4 for the
  tour-de-gross pair).  Proof shape: if the y-support of v has all
  cyclic gaps <= s_y - 1 it occupies >= ceil(m/s_y) rows; else v is
  y-windowed; if additionally x-windowed it is a plane cycle (the
  window margins kill aliasing), hence by Koszul a plane boundary,
  which reduces mod the torus to a member boundary — contradicting
  nontriviality; so it is x-spanning and occupies >= ceil(l/s_x)
  columns.  (Z-side identical by transpose symmetry.)

Gates (this script):
  K1: the check spans (machine, from the supports)      [done above: 4]
  K2: the regularity certificate — Res_y(y*A, x*y^3*B) != 0 in F2[x],
      computed exactly (5x5 Sylvester determinant over F2[x]).
  K3: counterexample hunt — NO banked certified nontrivial logical may
      be both-windowed (x-gap >= 4 AND y-gap >= 4); populations: bb72
      w6 (84), gross <= 12 (1,884, recomputed), (18,6) <= 18 + seam
      (5,807), the a36 two-gross witness, the [[432]] w24 witnesses.
  K4: the counting inequality direct check on every such logical
      (occupied rows/columns vs the gap structure), plus the gap-
      structure and per-row-load statistics that the L2 program needs
      (G2/G3 of the session charter: x-support widths, row loads).
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
    AxisDeck, TowerCode, colspace, h1_map, i2v, rep_for, rref_ints,
    span_points, v2i, validate_banked,
)
from a38_c37xx_freeze import census_pass  # noqa: E402

DATA = LAB / "data" / "a40"
A_L = [(0, 0), (0, 1), (3, -1)]
B_L = [(0, 0), (1, 0), (-1, -3)]
SPAN = 4


def red(supp, lm):
    return frozenset((e[0] % lm[0], e[1] % lm[1]) for e in supp)


def member(r, b, name=None):
    lm = (6 * (r + b), 6 * r)
    return TowerCode(name or f"tdg({r},{b})", lm, red(A_L, lm),
                     red(B_L, lm))


# ---------------------------------------------------- F2[x] arithmetic
def pmul(a, b):
    """Multiply F2[x] polys as int-coefficient lists mod 2."""
    if not a or not b:
        return []
    out = np.zeros(len(a) + len(b) - 1, dtype=np.int64)
    for i, c in enumerate(a):
        if c:
            out[i:i + len(b)] ^= np.array(b, dtype=np.int64)
    return [int(c % 2) for c in out]


def padd(a, b):
    n = max(len(a), len(b))
    return [((a[i] if i < len(a) else 0) ^ (b[i] if i < len(b) else 0))
            for i in range(n)]


def det_f2x(M):
    """Exact determinant of a matrix with F2[x] entries (cofactor)."""
    n = len(M)
    if n == 1:
        return M[0][0]
    out = []
    for j in range(n):
        if not any(M[0][j]):
            continue
        minor = [[M[i][jj] for jj in range(n) if jj != j]
                 for i in range(1, n)]
        out = padd(out, pmul(M[0][j], det_f2x(minor)))
    return out


def gap_structure(vals, order):
    """Cyclic gap analysis of a set of coordinate values in Z_order:
    (n_occupied, max_gap)."""
    vs = sorted(set(vals))
    if not vs:
        return 0, order
    gaps = [(vs[(i + 1) % len(vs)] - vs[i]) % order - 1
            for i in range(len(vs))]
    if len(vs) == 1:
        gaps = [order - 1]
    return len(vs), max(gaps)


def audit(code, v, tag, stats):
    ng = code.ng
    xs, ys = [], []
    row_load: dict[int, int] = {}
    for i in np.nonzero(v)[0]:
        g = code.G.from_index(int(i) % ng)
        xs.append(g[0])
        ys.append(g[1])
        row_load[g[1]] = row_load.get(g[1], 0) + 1
    l, m = code.G.orders
    nx, gx = gap_structure(xs, l)
    ny, gy = gap_structure(ys, m)
    w = int(v.sum())
    x_windowed = gx >= SPAN
    y_windowed = gy >= SPAN
    assert not (x_windowed and y_windowed), \
        f"{tag}: BOTH-WINDOWED nontrivial logical — Koszul REFUTED " \
        f"(w={w}, max gaps x={gx} y={gy})"
    if not y_windowed:
        assert ny >= -(-m // SPAN), (tag, ny, m)
        assert w >= -(-m // SPAN)
    if y_windowed:
        assert not x_windowed and nx >= -(-l // SPAN) and w >= -(-l // SPAN)
    stats["n"] += 1
    stats["min_row_load"] = min(stats["min_row_load"],
                                min(row_load.values()))
    stats.setdefault("x_width_hist", {})
    xw = l - gx - 1 if gx < l else 0     # occupied x-extent
    stats["x_width_hist"][str(xw)] = \
        stats["x_width_hist"].get(str(xw), 0) + 1
    stats.setdefault("sector_hist", {})
    sec = ("y-wrap+x-wrap" if not y_windowed and not x_windowed else
           "y-wrap only" if not y_windowed else "x-wrap only")
    stats["sector_hist"][sec] = stats["sector_hist"].get(sec, 0) + 1


def main():
    t0 = time.time()
    validate_banked(LAB / "data")
    print("validate_banked: PASS")
    out = {}

    # ---- K2: the regularity certificate
    # ytA = y*A = x^3 + y + y^2 : y-coeffs [x^3, 1, 1]
    # xy3B = x*y^3*B = 1 + (x + x^2) y^3 : y-coeffs [1, 0, 0, x+x^2]
    a0, a1, a2 = [0, 0, 0, 1], [1], [1]          # F2[x] coeff lists
    b0, b1, b2, b3 = [1], [], [], [0, 1, 1]
    Z = []
    Syl = [
        [a0, a1, a2, Z, Z],
        [Z, a0, a1, a2, Z],
        [Z, Z, a0, a1, a2],
        [b0, b1, b2, b3, Z],
        [Z, b0, b1, b2, b3],
    ]
    res = det_f2x(Syl)
    res_str = " + ".join(f"x^{i}" for i, c in enumerate(res) if c) \
        or "0"
    nonzero = any(res)
    print(f"K2 resultant Res_y(y*A, x*y^3*B) = {res_str}  "
          f"(nonzero: {nonzero})")
    assert nonzero, "REGULARITY FAILS — the Koszul route is dead"
    out["K2_resultant"] = res_str

    # sanity cross-check of K2 numerically: count common zeros of A,B
    # over F_{2^12}* x F_{2^12}* via the member quotient dims already
    # banked (k = 12 = 2*dim quotient bounded) — recorded, not needed.

    # ---- K3/K4 populations
    binp = cosetbz.build_kernel()
    stats = {"n": 0, "min_row_load": 99}

    bb72 = member(1, 0, "bb72")
    CH = 51
    allb = sorted(range(1, 1 << bb72.k))
    n6 = 0
    for lo in range(0, len(allb), CH):
        chunk = allb[lo:lo + CH]
        hits = census_pass(binp, bb72,
                           [(f"C{c}", rep_for(bb72, c)) for c in chunk],
                           6, f"s3_bb72_{lo}")
        for c in chunk:
            for h in sorted(hits[f"C{c}"]):
                v = i2v(h, bb72.n)
                audit(bb72, v, "bb72-w6", stats)
                n6 += 1
    print(f"K3/K4 bb72: {n6} w6 logicals audited")

    gross = member(1, 1, "gross")
    n12 = 0
    for lo in range(0, len(allb), CH):
        chunk = allb[lo:lo + CH]
        hits = census_pass(binp, gross,
                           [(f"C{c}", rep_for(gross, c)) for c in chunk],
                           12, f"s3_gross_{lo}")
        for c in chunk:
            for h in sorted(hits[f"C{c}"]):
                audit(gross, i2v(h, gross.n), "gross-w12", stats)
                n12 += 1
    print(f"K3/K4 gross: {n12} nontrivial <= 12 audited")

    b2f = member(1, 2, "b2frame")
    n18 = 0
    for fn in ("ckpt_W22_ntrv1.jsonl", "ckpt_W22_seam1.jsonl"):
        for line in (DATA / "tdg432" / fn).open():
            r = json.loads(line)
            v = np.zeros(b2f.n, dtype=np.uint8)
            v[r["support"]] = 1
            assert b2f.is_cycle(v) and not b2f.is_stab(v)
            audit(b2f, v, "(18,6)", stats)
            n18 += 1
    print(f"K3/K4 (18,6): {n18} nontrivial <= 22 audited")

    # a36 witness (paper frame, as in S2)
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
    audit(tg, v_p, "two-gross-w18", stats)

    # [[432]] w24 witnesses (from the L12 closed form)
    L12 = [tuple(t) for t in
           json.loads((DATA / "s2_ub_bands.json").read_text())["L12"]]
    t21 = member(2, 1)
    v = np.zeros(t21.n, dtype=np.uint8)
    for j in range(2):
        for blk, gx, gy in L12:
            v[blk * t21.ng + t21.G.index((gx, (gy + 6 * j) % 12))] ^= 1
    assert t21.is_cycle(v) and not t21.is_stab(v)
    audit(t21, v, "[[432]]-w24", stats)

    print(f"\nK3 PASS: {stats['n']} certified nontrivial logicals, "
          f"ZERO both-windowed (Koszul-consistent), all counting "
          f"inequalities hold")
    print(f"K4 stats: sectors {stats['sector_hist']}; min per-row "
          f"load = {stats['min_row_load']}; x-extent histogram "
          f"(top): "
          f"{sorted(stats['x_width_hist'].items(), key=lambda kv: -kv[1])[:8]}")
    out["K3_n_audited"] = stats["n"]
    out["K4"] = {k: stats[k] for k in
                 ("sector_hist", "min_row_load", "x_width_hist")}
    out["wall_s"] = round(time.time() - t0, 1)
    (DATA / "s3_l1_gates.json").write_text(json.dumps(out, indent=1))
    print(f"wrote {DATA/'s3_l1_gates.json'} ({out['wall_s']} s)")


if __name__ == "__main__":
    main()
