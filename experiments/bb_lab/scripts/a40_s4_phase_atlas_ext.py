#!/usr/bin/env python3
"""A40 S4 — compact-phase atlas extension: periods 7, 8 below rate 2,
and the p = 6 nontrivial spectrum just above the floor (Wcap 14).

Same engine and soundness as a40_s4_phase_atlas.py (imported).
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

LAB = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(LAB / "src"))
sys.path.insert(0, str(LAB / "scripts"))

from a40_s4_phase_atlas import DATA, atlas  # noqa: E402
from bb_lab.tower import validate_banked  # noqa: E402


def main():
    t0 = time.time()
    validate_banked(LAB / "data")
    print("validate_banked: PASS")
    summary, nontrivial = [], []
    jobs = [("AB", 7, 13), ("BAbar", 7, 13),
            ("AB", 8, 15), ("BAbar", 8, 15),
            ("AB", 6, 14)]
    for pair, p, W in jobs:
        t1 = time.time()
        rows, npop = atlas(pair, p, W, max_states=30_000_000,
                           max_paths=3_000_000)
        wall = round(time.time() - t1, 1)
        ntv = sorted((r["weight"] for r in rows if r["nontrivial"]))
        import collections
        spec = dict(collections.Counter(ntv))
        summary.append(dict(pair=pair, p=p, Wcap=W, n_cycles=len(rows),
                            n_nontrivial=len(ntv),
                            nontrivial_spectrum=spec,
                            states_popped=npop, wall_s=wall))
        nontrivial += [r for r in rows if r["nontrivial"]]
        print(f"{pair} p={p} W<={W}: {len(rows)} compact cycles, "
              f"nontrivial spectrum {spec}, {npop} states, {wall} s",
              flush=True)
    below = [s for s in summary
             if s["p"] in (7, 8) and s["n_nontrivial"]]
    if below:
        print("!! nontrivial phases below rate 2 at p in {7,8}:", below)
    else:
        print("p in {7,8}: no nontrivial x-compact phase below rate 2 "
              "(sweep now covers p <= 8, both pairs)")
    out = dict(summary=summary, wall_s=round(time.time() - t0, 1))
    (DATA / "s4_phase_atlas_ext.json").write_text(
        json.dumps(out, indent=1))
    print(f"wrote {DATA/'s4_phase_atlas_ext.json'} ({out['wall_s']} s)")


if __name__ == "__main__":
    main()
