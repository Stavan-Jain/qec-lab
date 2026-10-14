#!/usr/bin/env python3
"""A40 S2 / P4 step 3c — the cylinder dual certificate.

Target: u with period 6 in x (so it exists identically at every
l = 6(r+b)), supported in a bounded y-window [0, 6h) with NO y-wrap
usage (so it is a Z-cycle at every y-order >= 6h + reach), and
<u, v_{r,1}> = 1 (the pairing meets only v's bands inside the window,
so it is one finite computation).

Construction space: Z-side cycles of the CYLINDER code (x in Z_6,
y in the window, checks free in y).  Any kernel element placed
x-periodically into a member is a Z-cycle there (periodicity absorbs
the x-wrap; y-locality absorbs the y-wrap) — asserted mechanically at
every member r <= 6 below.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

LAB = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(LAB / "src"))

from bb_lab.tower import (  # noqa: E402
    TowerCode, kernel_basis, validate_banked,
)

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


def v_r1(code, L12, r):
    v = np.zeros(code.n, dtype=np.uint8)
    for j in range(r):
        v ^= place(code, L12, (0, 6 * j))
    return v


def cylinder_kernel(h):
    """Z-cycles of the x-periodic (period 6), y-window [0, 6h) strip:
    columns (blk, x in Z6, y in [0, 6h)); X-check rows (x in Z6,
    y in [-3, 6h+2)) with x mod 6, y FREE (no wrap).  Returns kernel
    basis as lists of (blk, x, y)."""
    ys = list(range(6 * h))
    cols = [(blk, x, y) for blk in (0, 1) for x in range(6) for y in ys]
    cidx = {c: i for i, c in enumerate(cols)}
    rows = []
    for gx in range(6):
        for gy in range(-3, 6 * h + 3):
            row = np.zeros(len(cols), dtype=np.uint8)
            touched = False
            for (sx, sy) in A_L:            # L-qubits: h = g - s
                q = (0, (gx - sx) % 6, gy - sy)
                if q in cidx:
                    row[cidx[q]] ^= 1
                    touched = True
            for (sx, sy) in B_L:            # R-qubits
                q = (1, (gx - sx) % 6, gy - sy)
                if q in cidx:
                    row[cidx[q]] ^= 1
                    touched = True
            if touched:
                rows.append(row)
    M = np.array(rows, dtype=np.uint8)
    ker = kernel_basis(M)
    out = []
    for kv in ker:
        out.append([cols[i] for i in np.nonzero(kv)[0]])
    return out


def main():
    validate_banked(LAB / "data")
    print("validate_banked: PASS")
    L12 = [tuple(t) for t in
           json.loads((DATA / "s2_ub_bands.json").read_text())["L12"]]
    out = {}
    found = None
    for h in (1, 2, 3):
        ker = cylinder_kernel(h)
        print(f"cylinder h={h}: kernel dim {len(ker)}")
        if not ker:
            continue
        # pairing with v's window content: bands j < h of L12 at
        # x-offsets 0 (v is x-local in [0,12); u is 6-periodic in x)
        best = None
        for kv in ker:
            pat = {(blk, x, y) for blk, x, y in kv}
            pair = 0
            for j in range(min(h, 99)):
                for blk, gx, gy in L12:
                    if (blk, gx % 6, gy + 6 * j) in pat:
                        pair ^= 1
            if pair == 1 and (best is None or len(kv) < len(best)):
                best = kv
        if best is not None:
            found = {"h": h, "cyl_weight": len(best),
                     "triples": [list(t) for t in best]}
            print(f"  FOUND pairing-odd cylinder Z-cycle: h={h}, "
                  f"cylinder weight {len(best)}")
            break
        print(f"  all {len(ker)} kernel elements pair EVEN with the "
              f"L12 bands")
    out["dual"] = found
    if not found:
        print("NO cylinder dual through h = 3 — recorded; the all-r "
              "nontriviality certificate remains open")
        (DATA / "s2_ub_dual3.json").write_text(json.dumps(out, indent=1))
        return

    # verify at every member r = 1..6: periodic placement is a Z-cycle
    # and pairs odd with v_{r,1}
    print("member verification (periodic placement):")
    rows = {}
    for r in range(1, 7):
        c = member(r, 1)
        if 6 * r < 6 * found["h"]:
            continue
        ng = c.ng
        u = np.zeros(c.n, dtype=np.uint8)
        for blk, x, y in found["triples"]:
            for rep in range(c.G.orders[0] // 6):
                u[blk * ng + c.G.index(((x + 6 * rep) % c.G.orders[0],
                                        y))] ^= 1
        zc = not ((c.HX @ u) % 2).any()
        vv = v_r1(c, L12, r)
        pr = int((u & vv).sum()) % 2
        rows[r] = {"z_cycle": bool(zc), "pairing": pr,
                   "u_weight": int(u.sum())}
        print(f"  r={r} [[{c.n},12]]: Z-cycle={zc} <u,v>={pr} "
              f"|u|={int(u.sum())}")
    out["members"] = rows
    (DATA / "s2_ub_dual3.json").write_text(json.dumps(out, indent=1))
    print(f"wrote {DATA/'s2_ub_dual3.json'}")


if __name__ == "__main__":
    main()
