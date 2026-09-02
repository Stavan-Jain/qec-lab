#!/usr/bin/env python3
"""A40 S2 / P4 steps 2-3 — the band-form upper-bound witnesses.

Archaeology (s2_staircase_archaeology + s2_ub_hunt) says the proven
minimum witnesses have y-band weights:

    (2,0) w18: bands (12, 6);   (2,1) w24: bands (12, 12);
    (1,1) w12: band  (12) [x-local];  (1,0) w6: band (6).

Candidate closed forms (bands as (l, 6)-frame vectors placed at
y-offsets 6j):

    v_{r,1} = sum_{j=0}^{r-1} y^{6j} . L12     (the pure y-transfer of
              the x-local weight-12 pattern; cycle by pullback)
    v_{r,0} = sum_{j=0}^{r-2} y^{6j} . U0  +  y^{6(r-1)} . U1
              (U0/U1 extracted from the a36 two-gross witness; the
              telescoping glue is what this script tests)

This script extracts L12/U0/U1 explicitly, tests the generalized forms
at every member with r <= 6 (kernel checks: cycle + weight + non-stab,
all exact linear algebra), and finds y-local dual certificates
(Z-side cycles u with <u, v> = 1) whose form is r-independent.

Output: data/a40/s2_ub_bands.json (+ the explicit supports for the
note's theorem statement).
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

from bb_lab.tower import (  # noqa: E402
    TowerCode, in_span, rref_ints, v2i, validate_banked,
)

DATA = LAB / "data" / "a40"
A_L = [(0, 0), (0, 1), (3, -1)]
B_L = [(0, 0), (1, 0), (-1, -3)]


def red(supp, lm):
    return frozenset((e[0] % lm[0], e[1] % lm[1]) for e in supp)


def member(r, b, name=None):
    lm = (6 * (r + b), 6 * r)
    return TowerCode(name or f"tdg({r},{b})", lm, red(A_L, lm),
                     red(B_L, lm))


def to_pairs(code, v):
    """Support as (block, gx, gy) triples."""
    ng = code.ng
    out = []
    for i in np.nonzero(v)[0]:
        blk, gi = divmod(int(i), ng)
        g = code.G.from_index(gi)
        out.append((blk, g[0], g[1]))
    return sorted(out)


def place(code_big, triples, t=(0, 0)):
    """Place (blk, gx, gy) triples into code_big shifted by t."""
    ngb = code_big.ng
    v = np.zeros(code_big.n, dtype=np.uint8)
    for blk, gx, gy in triples:
        u = ((gx + t[0]) % code_big.G.orders[0],
             (gy + t[1]) % code_big.G.orders[1])
        v[blk * ngb + code_big.G.index(u)] ^= 1
    return v


def main():
    t0 = time.time()
    validate_banked(LAB / "data")
    print("validate_banked: PASS")
    out = {}

    # ---- L12: an x-local w12 (single-cell rep from the (18,6) ckpt)
    L1 = member(1, 2, "b2frame")   # (18,6)
    u12s = []
    for line in (DATA / "tdg432" / "ckpt_W22_ntrv1.jsonl").open():
        r = json.loads(line)
        if r["w"] != 12:
            continue
        u = np.zeros(L1.n, dtype=np.uint8)
        u[r["support"]] = 1
        u12s.append(u)
    # pick a single-cell one (all support x in [6,12), y in [0,6))
    L12_triples = None
    for u in u12s:
        tr = to_pairs(L1, u)
        if all(6 <= gx < 12 for _, gx, gy in tr):
            L12_triples = [(blk, gx - 6, gy) for blk, gx, gy in tr]
            break
    assert L12_triples is not None
    print(f"L12 (x-local w12, normalized to x in [0,6)): {L12_triples}")
    out["L12"] = L12_triples

    # sanity: L12 placed at x-offset 0 is a logical of gross and of any
    # (l>=12, 6) frame flat
    for lm in [(12, 6), (18, 6), (24, 6), (30, 6)]:
        c = TowerCode(f"f{lm}", lm, red(A_L, lm), red(B_L, lm))
        v = place(c, L12_triples)
        assert c.is_cycle(v) and not c.is_stab(v), lm
    print("  L12 flat is a nontrivial logical of (l,6) for "
          "l = 12,18,24,30  [x-window check]")

    # ---- U0/U1 from the a36 two-gross witness (paper frame)
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
    tr = to_pairs(tg, v_p)
    U0 = [(blk, gx, gy) for blk, gx, gy in tr if gy < 6]
    U1 = [(blk, gx, gy - 6) for blk, gx, gy in tr if gy >= 6]
    print(f"a36 witness bands: |U0| = {len(U0)}, |U1| = {len(U1)}")
    assert {len(U0), len(U1)} == {12, 6}
    if len(U0) == 6:
        U0, U1 = U1, U0    # normalize: U0 = the 12-band
    out["U0"] = U0
    out["U1"] = U1

    # is U0 a translate of L12?  (x-local & same pattern?)
    xs = [gx for _, gx, gy in U0]
    span_ok = (max(xs) - min(xs)) < 6
    print(f"  U0 x-span within one cell: {span_ok}; "
          f"U0 x-range [{min(xs)},{max(xs)}]")

    # ---- generalized forms, kernel-checked at every member r <= 6
    print("\nkernel checks of the closed forms:")
    results = []
    for r in range(1, 7):
        for b in (0, 1):
            code = member(r, b)
            d_hat = 6 * (2 * r + b - 1)
            if b == 1:
                v = np.zeros(code.n, dtype=np.uint8)
                for j in range(r):
                    v ^= place(code, L12_triples, (0, 6 * j))
                form = f"sum_j y^6j L12 (r={r} bands)"
            else:
                if r == 1:
                    # w6: band form (6) — use U1 alone?  U1 alone at
                    # (6,6) is the seam band; test it:
                    v = place(code, U1)
                    form = "U1 alone at (6,6)"
                else:
                    v = np.zeros(code.n, dtype=np.uint8)
                    for j in range(r - 1):
                        v ^= place(code, U0, (0, 6 * j))
                    v ^= place(code, U1, (0, 6 * (r - 1)))
                    form = f"(r-1) U0 bands + U1 (r={r})"
            w = int(v.sum())
            cyc = code.is_cycle(v)
            nst = (not code.is_stab(v)) if cyc else False
            ok = (w == d_hat) and cyc and nst
            results.append({"r": r, "b": b, "n": code.n,
                            "d_conj": d_hat, "form": form,
                            "weight": w, "cycle": bool(cyc),
                            "nontrivial": bool(nst), "ok": bool(ok)})
            print(f"  ({r},{b}) [[{code.n},12]] d^={d_hat}: |v|={w} "
                  f"cycle={cyc} nontrivial={nst}  "
                  f"{'OK' if ok else '-- FAIL'}   ({form})")
    out["kernel_checks"] = results

    out["wall_s"] = round(time.time() - t0, 1)
    (DATA / "s2_ub_bands.json").write_text(json.dumps(out, indent=1))
    print(f"\nwrote {DATA/'s2_ub_bands.json'} ({out['wall_s']} s)")


if __name__ == "__main__":
    main()
