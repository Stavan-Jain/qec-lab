#!/usr/bin/env python3
"""A40 S5 — straight compact-phase atlas extension to p = 9 (both
lanes, Wcap = 2p - 1 = 17): the §9.9 item-3 reach test.  Same engine
and soundness as the S4 atlas; state cap raised with a hard abort
recorded honestly if exceeded (memory guard)."""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

LAB = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(LAB / "src"))
sys.path.insert(0, str(LAB / "scripts"))

from bb_lab.tower import validate_banked  # noqa: E402

_argv = sys.argv
sys.argv = [_argv[0], "18", "8"]
from a40_s4_phase_atlas import DATA, atlas  # noqa: E402
sys.argv = _argv


def main():
    t0 = time.time()
    validate_banked(LAB / "data")
    print("validate_banked: PASS", flush=True)
    out = {"summary": []}
    for pair in ("AB", "BAbar"):
        t1 = time.time()
        try:
            rows, npop = atlas(pair, 9, 17, max_states=30_000_000,
                               max_paths=3_000_000)
        except RuntimeError as ex:
            out["summary"].append(dict(pair=pair, p=9, Wcap=17,
                                       aborted=str(ex)))
            print(f"{pair} p=9: ABORT {ex} (recorded)", flush=True)
            continue
        ntv = sorted(r["weight"] for r in rows if r["nontrivial"])
        import collections
        spec = dict(collections.Counter(ntv))
        out["summary"].append(dict(
            pair=pair, p=9, Wcap=17, n_cycles=len(rows),
            n_nontrivial=len(ntv), nontrivial_spectrum=spec,
            states_popped=npop, wall_s=round(time.time() - t1, 1)))
        print(f"{pair} p=9 W<=17: {len(rows)} compact cycles, "
              f"nontrivial spectrum {spec}, {npop} states, "
              f"{time.time()-t1:.1f}s", flush=True)
    out["wall_s"] = round(time.time() - t0, 1)
    (DATA / "s5_atlas_p9.json").write_text(json.dumps(out, indent=1))
    print(f"wrote {DATA/'s5_atlas_p9.json'} ({out['wall_s']} s)")


if __name__ == "__main__":
    main()
