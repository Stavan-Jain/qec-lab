#!/usr/bin/env python3
"""A42 S0b — anchor identifications and the 127-line frame witness.

1. Identify the tangent F4 orbit: at alpha = omega the Jacobian
   determinant of (A~, B^) is 1 + beta*omega, which vanishes iff
   beta = omega^2 = alpha^2.  Verify the engine's mult-2 orbit is
   exactly the beta = alpha^2 (antidiagonal) orbit and the mult-1
   orbit is beta = alpha (diagonal).

2. The 127-orbit's discrete log: beta = alpha^g in F_128; report g.
   A cyclic frame Z_127 = Z^2/<(127,0),(d,1)> carries the 127-line
   iff d + g == 0 mod 127 (character condition alpha^d beta = 1).
   Verify with TowerCode: k(127,1,127-g) == 2*7 = 14, and k == 0 at
   two control shears.  This certifies the 127-line as real
   code-level structure at an affordable size (n = 254), grounding
   the k = 26 prediction at member (762,762).

3. The W7 chirality anchor: at frames (l,7,l-2) the contributing
   orbit is the TANGENT one (beta = alpha^2 forced by the character
   condition); at the +2 mirror it is the transverse orbit.  Assert
   both selections symbolically.
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

LAB = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(LAB / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from bb_lab.tower import TowerCode  # noqa: E402

import a42_lib as L  # noqa: E402

DATA = LAB / "data" / "a42"


def main():
    t0 = time.time()
    out = {}
    orbits = L.variety_orbits()
    for o in orbits:
        o.plane_mult = L.plane_multiplicity(o, N=8)

    # --- 1: tangency identification ---
    print("== tangency identification ==")
    f4 = [o for o in orbits if o.ord_a == 3]
    assert len(f4) == 2
    ident = []
    for o in f4:
        F = o.field
        is_diag = F.eq(o.beta, o.alpha)
        is_anti = F.eq(o.beta, F.mul(o.alpha, o.alpha))
        assert is_diag != is_anti
        name = "diagonal (beta=alpha)" if is_diag else \
            "antidiagonal (beta=alpha^2)"
        print(f"  orbit mult={o.plane_mult}: {name}")
        ident.append({"mult": o.plane_mult, "type": name})
        if o.plane_mult == 2:
            assert is_anti, "tangent orbit must be antidiagonal"
        else:
            assert is_diag, "transverse orbit must be diagonal"
    out["f4_identification"] = ident
    print("  matches the hand Jacobian computation "
          "(det = 1 + beta*omega at alpha = omega): PASS")

    # --- 2: the 127-line at a cyclic frame ---
    print("\n== 127-line frame witness ==")
    o127 = next(o for o in orbits if o.ord_a == 127)
    F = o127.field
    # discrete log g: beta = alpha^g
    g = None
    pw = F.one
    for e in range(127):
        if F.eq(pw, o127.beta):
            g = e
            break
        pw = F.mul(pw, o127.alpha)
    assert g is not None
    print(f"  beta = alpha^{g} in F_128")
    out["dlog_g"] = g
    d_star = (-g) % 127
    ks = L.spectral_k(orbits, 127, 1, d_star)
    sa, sb, fd, coll = L.transported_supports(127, 1, d_star)
    c = TowerCode("w127", (fd.o1, fd.o2), sa, sb)
    print(f"  frame (127,1,{d_star}): spectral k={ks}, TowerCode k={c.k}")
    assert ks == c.k == 14, (ks, c.k)
    controls = []
    for d in (d_star + 1) % 127, (d_star + 63) % 127:
        ks2 = L.spectral_k(orbits, 127, 1, d)
        sa, sb, fd, _ = L.transported_supports(127, 1, d)
        c2 = TowerCode("ctrl", (fd.o1, fd.o2), sa, sb)
        assert ks2 == c2.k == 0, (d, ks2, c2.k)
        controls.append(d)
    print(f"  controls d={controls}: k=0 both, spectral==TowerCode: PASS")
    out["frame_127"] = {"d_star": d_star, "k": 14, "controls": controls}

    # --- 3: W7 chirality anchor ---
    print("\n== W7 chirality anchor ==")
    rows = []
    for l in (12, 18, 24):
        for d, expect in ((l - 2, "antidiagonal"), (2, "diagonal")):
            k, contrib, fd = L.spectral_k(orbits, l, 7, d, detail=True)
            assert len(contrib) == 1, (l, d, contrib)
            o, m = contrib[0]
            Fo = o.field
            typ = "diagonal" if Fo.eq(o.beta, o.alpha) else "antidiagonal"
            assert typ == expect, (l, d, typ, expect)
            assert k == 4
            rows.append({"l": l, "d": d, "orbit": typ, "mult": o.plane_mult})
            print(f"  (l,7,{d}): k=4 via the {typ} orbit "
                  f"(plane mult {o.plane_mult})")
    print("  chirality: the weight-8 species side (d=l-2) is the TANGENT"
          " orbit; the empty mirror (d=+2) is the transverse orbit")
    out["w7_chirality"] = rows

    out["wall_s"] = round(time.time() - t0, 1)
    (DATA / "s0b_anchors.json").write_text(json.dumps(out, indent=1))
    print(f"\nwrote {DATA/'s0b_anchors.json'} ({out['wall_s']} s)")


if __name__ == "__main__":
    main()
