#!/usr/bin/env python3
"""A40 S4 — theorem B6, the finite burden: W = 11 all-class censuses
at the four remaining small family frames (l,6), l in {24,30,36,42}.

Assembly of B6 (d((l,6)) = 12 for every l = 0 mod 6, l >= 12):
  - upper: L12 places nontrivially at every (l,6) (gate W2, l-uniform
    by x-locality: no x-wrap is used, S2 §7.1);
  - lower, x-windowed logicals (any l): the compact-phase atlas at
    p = 6 is EMPTY of nontrivial cycles <= 11, and every sub-12
    compact cycle is COMPACTLY trivial (gate 4), so an x-windowed
    nontrivial (l,6)-logical of weight <= 11 cannot exist at ANY l
    (its window lift would be such a cycle: margins >= 4 kill x-wrap
    aliasing — the one-axis Lemma-K windowing move);
  - lower, x-spanning logicals: weight >= ceil(l/4) >= 12 for l >= 45
    (Theorem L1); for l in {12, 18} certified (gross Lean-grade /
    a40 §3.2 census-complete); for l in {24, 30, 36, 42}: THIS
    script — complete W = 11 coset-BZ censuses over all 2^k - 1
    nontrivial classes, all required EMPTY.
Control: at (24,6), one W = 12 run of the L12 class must be NONEMPTY
(engine-wiring positive control).

STATUS (S4, recorded in note §9.6): the direct route ABORTS — the
compiled walk kernel caps at n <= 192 qubits (cosetbz.NMAX = 3x64-bit
words) and these frames have n = 288..504.  The four certificates
were then CLOSED the same session by Z2-descent with walks at n <=
192 frames only — see a40_s4_b6_close.py.  This script guards the cap
explicitly and exits; it is kept as the assembly record.
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
    TowerCode, i2v, rep_for, v2i, validate_banked,
)
from a38_c37xx_freeze import census_pass  # noqa: E402

DATA = LAB / "data" / "a40"
A_L = [(0, 0), (0, 1), (3, -1)]
B_L = [(0, 0), (1, 0), (-1, -3)]


def red(supp, lm):
    return frozenset((e[0] % lm[0], e[1] % lm[1]) for e in supp)


def code_at(lm, name=None):
    return TowerCode(name or f"tdg{lm}", lm, red(A_L, lm), red(B_L, lm))


def main():
    t0 = time.time()
    validate_banked(LAB / "data")
    print("validate_banked: PASS", flush=True)
    if cosetbz.NMAX < 2 * 24 * 6:
        print(f"walk kernel NMAX = {cosetbz.NMAX} < n = 288..504: the "
              f"direct censuses cannot run (see docstring STATUS). "
              f"The four B6 x-spanning certificates remain S5 residue "
              f"(descent route or NMAX lift).")
        sys.exit(2)
    binp = cosetbz.build_kernel()
    out = {"frames": []}

    # ---- positive control: L12's class at (24,6), W = 12, must hit
    c24 = code_at((24, 6))
    assert c24.k == 12, c24.k
    L12 = [tuple(t) for t in
           json.loads((DATA / "s2_ub_bands.json").read_text())["L12"]]
    vL = np.zeros(c24.n, dtype=np.uint8)
    for blk, gx, gy in L12:
        vL[blk * c24.ng + c24.G.index((gx % 24, gy % 6))] ^= 1
    assert c24.is_cycle(vL) and not c24.is_stab(vL)
    sig = c24.sig(vL)
    sig_int = int("".join(map(str, sig[::-1])), 2)
    hits = census_pass(binp, c24, [("L12class", rep_for(c24, sig_int))],
                       12, "s4b6_ctrl")
    n12 = sum(1 for h in hits["L12class"]
              if int(i2v(h, c24.n).sum()) == 12)
    assert n12 > 0, "control failed: W=12 census misses L12's class"
    print(f"control PASS: (24,6) L12-class W=12 census nonempty "
          f"({n12} weight-12 elements)", flush=True)

    # ---- the four frames, W = 11, all classes empty
    CH = 51
    for l in (24, 30, 36, 42):
        t1 = time.time()
        c = code_at((l, 6))
        assert c.k == 12, (l, c.k)
        allb = sorted(range(1, 1 << c.k))
        n_found = 0
        for lo in range(0, len(allb), CH):
            chunk = allb[lo:lo + CH]
            hits = census_pass(
                binp, c,
                [(f"C{cc}", rep_for(c, cc)) for cc in chunk],
                11, f"s4b6_{l}_{lo}")
            for cc in chunk:
                for h in hits[f"C{cc}"]:
                    v = i2v(h, c.n)
                    assert c.is_cycle(v) and not c.is_stab(v)
                    n_found += 1
        wall = round(time.time() - t1, 1)
        print(f"({l},6): all {len(allb)} nontrivial classes censused "
              f"at W = 11: {n_found} elements found ({wall} s)",
              flush=True)
        assert n_found == 0, (l, "d < 12 element found!!")
        out["frames"].append(dict(l=l, W=11, classes=len(allb),
                                  nontrivial_found=0, wall_s=wall))

    print("\nVERDICT: d((l,6)) >= 12 for l in {24,30,36,42} at "
          "certificate tier (complete empty W=11 class censuses); "
          "with L12 (<= 12), L1 (l >= 45), the p=6 atlas + compact "
          "triviality (x-windowed, all l), and gross/(18,6): "
          "THEOREM B6: d((l,6)) = 12 for every l = 0 mod 6, l >= 12.")
    out["wall_s"] = round(time.time() - t0, 1)
    (DATA / "s4_b6_frames.json").write_text(json.dumps(out, indent=1))
    print(f"wrote {DATA/'s4_b6_frames.json'} ({out['wall_s']} s)")


if __name__ == "__main__":
    main()
