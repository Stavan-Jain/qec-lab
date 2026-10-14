#!/usr/bin/env python3
"""A40 S2 / P4 step 3b — hunt a Y-LOCAL dual certificate for v_{r,1}.

Want u with: H_X u = 0 (Z-side cycle), supp(u) inside a bounded
y-window, <u, v_{r,1}> = 1.  Y-locality gives r-independence of both
the cycle condition and the pairing (v's bands repeat with y-period 6,
u touches finitely many bands: parity = a fixed finite computation),
so verifying it once + small-r checks closes nontriviality for all r
in the hand proof.

Method: restrict H_X's columns to a y-window of h bands
(y in [0, 6h)), compute the kernel of the restricted matrix (exact),
and search that kernel for an element with odd pairing against
v_{r,1}'s window content.  Window h = 1, 2, 3 tried at (2,1); any hit
is then re-verified at (r,1), r = 1..6 by flat placement.
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np

LAB = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(LAB / "src"))

from bb_lab.tower import (  # noqa: E402
    TowerCode, kernel_basis, v2i, validate_banked,
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


def main():
    validate_banked(LAB / "data")
    print("validate_banked: PASS")
    L12 = [tuple(t) for t in
           json.loads((DATA / "s2_ub_bands.json").read_text())["L12"]]
    out = {}

    code = member(2, 1)          # (18,12)
    v = v_r1(code, L12, 2)
    ng = code.ng
    found = None
    for h in (1, 2):
        # window columns: y in [0, 6h) both blocks
        cols = [blk * ng + code.G.index((gx, gy))
                for blk in (0, 1)
                for gx in range(18) for gy in range(6 * h)]
        cols = np.array(sorted(cols))
        Hw = code.HX[:, cols]
        ker = kernel_basis(Hw)
        print(f"window h={h}: {len(cols)} qubits, kernel dim {len(ker)}")
        if not ker:
            continue
        # search the kernel span greedily: need odd pairing with v|window
        vw = v[cols]
        pair = np.array([int((kv & vw).sum()) % 2 for kv in ker],
                        dtype=np.uint8)
        idx = np.nonzero(pair)[0]
        if len(idx) == 0:
            print(f"  all {len(ker)} kernel basis elements pair EVEN")
            continue
        # pick the sparsest pairing-odd basis element (greedy)
        best = min((int(ker[i].sum()), int(i)) for i in idx)
        kv = ker[best[1]]
        u = np.zeros(code.n, dtype=np.uint8)
        u[cols] = kv
        assert not ((code.HX @ u) % 2).any()
        assert int((u & v).sum()) % 2 == 1
        tr = sorted((int(i) // ng, *code.G.from_index(int(i) % ng))
                    for i in np.nonzero(u)[0])
        found = {"h": h, "weight": int(u.sum()), "triples": tr}
        print(f"  FOUND y-local dual: window h={h}, weight "
              f"{found['weight']}")
        break
    out["dual"] = found
    if not found:
        print("NO y-local dual in windows h <= 2 — record and stop")
        (DATA / "s2_ub_dual2.json").write_text(json.dumps(out, indent=1))
        return

    # re-verify the SAME pattern at r = 1..6
    print("re-verification at (r,1), r = 1..6:")
    U = [tuple(t) for t in found["triples"]]
    rows = {}
    for r in range(1, 7):
        c = member(r, 1)
        if 6 * r < 6 * found["h"]:
            continue
        u = place(c, U)
        vv = v_r1(c, L12, r)
        zc = not ((c.HX @ u) % 2).any()
        pr = int((u & vv).sum()) % 2
        rows[r] = {"z_cycle": bool(zc), "pairing": pr}
        print(f"  r={r} [[{c.n},12]]: Z-cycle={zc} <u,v>={pr}")
    out["members"] = rows
    (DATA / "s2_ub_dual2.json").write_text(json.dumps(out, indent=1))
    print(f"wrote {DATA/'s2_ub_dual2.json'}")


if __name__ == "__main__":
    main()
