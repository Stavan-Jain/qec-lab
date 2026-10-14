#!/usr/bin/env python3
"""A42 S5 Stage 4 — the Z_9-fibre weight table (the r = 3 / 9 | p case).

p = 9 q' with 3 coprime to q'.  Z_{9q'} = Z_9 x Z_{q'} and
F2[Z_9] = F2[y]/(y^9+1) = F2 x F4 x F64  (y^9+1 = (y+1)(y^2+y+1)(y^6+y^3+1)),
so the fibre over a cell of the period-q' cylinder is a Z_9-coset, and
per fibre the 9 bits (v_0..v_8) (v_i = the cell y = j + i q') are the
DFT triple (s, mu, nu) in F2 x F4 x F64:

    s  = sum_i v_i                         (the order-1 character, barren)
    mu = sum_i v_i zeta^{i}   (zeta^3 = 1)  (the order-3 character, OMEGA)
    nu = sum_i v_i eta^{i}    (eta^9 = 1, eta^3 = zeta) (order-9, BARREN by
                                            Theorem A: 9 not in {3,127})

Homology lives in the mu-factor alone; (s, nu) are FREE barren data
(boundaries of an exact barren complex).  The per-fibre Hamming weight
w(s, mu, nu) is what the Z_9 hiding-mass inequality must control:

    wt(v) = sum_fibres w(s_f, mu_f, nu_f),   HM_9:  wt >= 18 q' = 2p.

This script tabulates w by (s, mu != 0, nu != 0) and by the finer
invariants, checks the pure-fibre weight (s = 0, nu = 0, mu != 0) = 6,
the minimal mixed fibre weight 1, and the fibre-count law, and writes
data/a42/s5_z9fibre.json.
"""
from __future__ import annotations

import itertools
import json
import sys
from pathlib import Path

LAB = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(Path(__file__).resolve().parent))
import a42_lib as AL  # noqa: E402

DATA = LAB / "data" / "a42"

# F64 = F2[eta]/(eta^6 + eta^3 + 1); eta has order 9; zeta = eta^3.
F64 = AL.F2k(0b1001001)
ETA = 0b10
ZETA = F64.pow(ETA, 3)
assert F64.pow(ETA, 9) == 1 and F64.pow(ETA, 3) != 1
assert F64.add(F64.add(F64.mul(ZETA, ZETA), ZETA), 1) == 0


def dft(v):
    s = sum(v) % 2
    mu = 0
    nu = 0
    for i, b in enumerate(v):
        if b:
            mu = F64.add(mu, F64.pow(ZETA, i))
            nu = F64.add(nu, F64.pow(ETA, i))
    return s, mu, nu


def main():
    tab = {}
    by_mu = {}
    full = {}
    for bits in itertools.product((0, 1), repeat=9):
        s, mu, nu = dft(bits)
        w = sum(bits)
        key = (s, int(mu != 0), int(nu != 0))
        tab.setdefault(key, []).append(w)
        by_mu.setdefault((int(mu != 0), w), 0)
        by_mu[(int(mu != 0), w)] += 1
        full[(s, mu, nu)] = w
    assert len(full) == 512, "DFT not bijective?!"
    print("(s, mu!=0, nu!=0) -> weight multiset")
    out = {}
    for key in sorted(tab):
        ws = sorted(tab[key])
        hist = {}
        for w in ws:
            hist[w] = hist.get(w, 0) + 1
        print(f"  s={key[0]} mu{'!=' if key[1] else '=='}0 "
              f"nu{'!=' if key[2] else '=='}0 : n={len(ws)} min={ws[0]} "
              f"max={ws[-1]} hist={hist}")
        out[str(key)] = {"n": len(ws), "min": ws[0], "max": ws[-1],
                         "hist": {str(k): v for k, v in sorted(hist.items())}}
    pure = tab[(0, 1, 0)]
    assert set(pure) == {6}, pure
    print("pure omega fibre (s=0, nu=0, mu!=0): weight 6 exactly, "
          f"{len(pure)} fibres (= 3 nonzero mu values)")
    mixed_min = min(w for key, ws in tab.items() if key[1] for w in ws)
    assert mixed_min == 1
    # the omega-support weight law: weight of a fibre given mu != 0
    print("weight histogram over fibres with mu != 0 vs mu == 0:")
    for m in (1, 0):
        hist = {w: c for (mm, w), c in sorted(by_mu.items()) if mm == m}
        print(f"  mu{'!=' if m else '=='}0: {hist}")
    # the Z_3-quotient view: fold Z_9 -> Z_3 (sum over the order-3 coset)
    # gives the (s, mu) data of the Z_3 fibre; nu is the extra barren
    # information.  Discount table: for each (s, mu != 0) what nu buys.
    out["pure_fibre_weight"] = 6
    out["min_mixed_fibre_weight"] = mixed_min
    (DATA / "s5_z9fibre.json").write_text(json.dumps(out, indent=1))
    print("->", DATA / "s5_z9fibre.json")


if __name__ == "__main__":
    main()
