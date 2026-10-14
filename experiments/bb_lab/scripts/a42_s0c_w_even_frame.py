#!/usr/bin/env python3
"""A42 S0c — second-method grounding for the k = 26 member claim.

The (762,762) member is too large for direct rank; instead verify the
W-line's contribution at an EVEN-2-part sheared frame of feasible
size: L = <(254,0),(108,3)>, |G| = 762 = 2*3*127 — the smallest kind
of frame where the W-orbit switches on together with a nontrivial
2-part.  The spectral formula's predicted k is compared against
TowerCode on the honestly-cancelled transported supports (n = 1524).
Controls: two neighboring shears.

Combined with s0b's odd cyclic frame (127,1,36) (k = 14 = pure W),
every ingredient of k(762,762) = 26 — the W switch-on conditions,
its local dim 1 at even 2-parts, and the F4 dims — is then
TowerCode-grounded at some feasible frame; the member value itself
remains formula-tier (the formula being exact by the CRT theorem).
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
    orbits = L.variety_orbits()
    for o in orbits:
        o.plane_mult = L.plane_multiplicity(o, N=8)
    out = {"frames": []}
    for (l, p, d) in ((254, 3, 108), (254, 3, 109), (254, 3, 107)):
        ks, contrib, fd = L.spectral_k(orbits, l, p, d, detail=True)
        det = [{"fx": L.poly_str(o.fx), "ord": (o.ord_a, o.ord_b),
                "D": o.D, "localdim": m} for (o, m) in contrib]
        print(f"frame ({l},{p},{d}): spectral k = {ks} via {det}; "
              f"orders ({fd.o1},{fd.o2}) 2-part ({fd.a1},{fd.a2})",
              flush=True)
        sa, sb, fdd, coll = L.transported_supports(l, p, d)
        c = TowerCode(f"wf{d}", (fdd.o1, fdd.o2), sa, sb)
        print(f"  TowerCode k = {c.k} (collision={coll}) "
              f"[{time.time()-t0:.0f} s]", flush=True)
        assert ks == c.k, (l, p, d, ks, c.k)
        out["frames"].append({"l": l, "p": p, "d": d, "k": ks,
                              "contrib": det,
                              "two_part": [fd.a1, fd.a2]})
    out["wall_s"] = round(time.time() - t0, 1)
    (DATA / "s0c_w_even_frame.json").write_text(json.dumps(out, indent=1))
    print(f"wrote {DATA/'s0c_w_even_frame.json'} ({out['wall_s']} s)")


if __name__ == "__main__":
    main()
