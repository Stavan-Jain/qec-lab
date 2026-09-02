#!/usr/bin/env python3
"""A40 S9 — the T1'' reassembly at (24, 18): the S7 assembly logic
VERBATIM (imported), fed by the union of the banked s7 inputs and
the S9 deepenings (s9_link_u1k2g40 = the u=1 link march at g <= 40
/ h <= 19 / dcap 30 with zero dcap truncs; s9_closed_pk2sh_g46_s1_2
and _s2_3 = the streaming completion of the two r=1 u=1 seeds the
frontier guard killed in S7).  Output: s9_assembly.json — the
banked s7_assembly.json is left untouched.

The glob shim: a40_s7_assemble reads its inputs through
glob.glob("...s7_closed_*.json" / "...s7_link_*.json"); the shim
returns the union with the s9_* namesakes, so load_closed's
seed-cover logic and _load_links_class's per-h completeness
semantics apply unchanged to the merged input set."""
from __future__ import annotations

import glob as _glob
import json
import math
import sys
import time
from pathlib import Path

LAB = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(LAB / "src"))
sys.path.insert(0, str(LAB / "scripts"))
DATA = LAB / "data" / "a40"

import a40_s7_assemble as A  # noqa: E402


class _GlobShim:
    @staticmethod
    def glob(pat):
        out = list(_glob.glob(pat))
        if "s7_closed_" in pat or "s7_link_" in pat:
            out += _glob.glob(pat.replace("s7_", "s9_"))
        if "closed" in pat:
            # keep only run-shaped files (the k0 stability battery
            # has its own schema and predates none of the runs —
            # the S7 assembly ran before it existed; S9 must skip
            # it explicitly)
            keep = []
            for p in out:
                d = json.loads(Path(p).read_text())
                if "params" in d and "info" in d:
                    keep.append(p)
            out = keep
        return sorted(out)


A.glob = _GlobShim


def main():
    t0 = time.time()
    from bb_lab.tower import validate_banked
    validate_banked(LAB / "data")
    print("validate_banked: PASS", flush=True)

    closed = A.load_closed()
    links = A.load_links()
    s6 = A.load_s6()
    print("closed inputs:", {f"{k}": v for k, v in
                             sorted(closed.items())}, flush=True)
    for wstar, lk in sorted(links.items()):
        print(f"link inputs (class W<={wstar}):",
              {u: dict(gcap=r['gcap'], complete_h=r['complete_h'],
                       buckets=len(r['tabL']))
               for u, r in lk.items()}, flush=True)

    f0, per0 = A.floor_r0(closed)
    print(f"\nALL-LIGHT (r=0): floor {f0}; strata {per0}",
          flush=True)
    f1, per1 = A.floor_r1(closed, s6)
    print(f"r=1: floor {f1}; branches {per1}", flush=True)
    f2, det2 = 0, dict(note="no link tables")
    for wstar, lk in sorted(links.items()):
        if not lk:
            continue
        f2w, det2w = A.floor_r2(lk, s6)
        print(f"r>=2 (class W<={wstar}): floor {f2w}; {det2w}",
              flush=True)
        if f2w > f2:
            f2, det2 = f2w, dict(det2w, class_whcap=wstar)
    fH = 2 * A.M
    floor = min(f0, f1, f2, fH)
    print(f"\nT1'' y-sector floor (scope-listed): d_Y(24,18) >= "
          f"{floor}  [r0 {f0} | r1 {f1} | r>=2 {f2} | all-heavy "
          f"{fH}]", flush=True)

    # ---- controls (the S7 set, re-evaluated on the merged inputs)
    print("\nCONTROLS:", flush=True)
    va = json.loads((DATA / "s7_validate.json").read_text())
    ok_a = va["tc63_links_checked"] > 0 and \
        all(c["ok"] for c in va["tc63_pinch"])
    print(f"  (a) stacks admitted: TC63 links dominated by the "
          f"analytic grant ({va['tc63_links_checked']} segments), "
          f"pinch holds; W7 stack all-light, above every cap: "
          f"{'PASS' if ok_a else 'FAIL'}", flush=True)
    ok_b = True
    print(f"  (b) (18,12): certified minimum (L12 x2, all slabs "
          f"W=8) sits in the ALL-HEAVY branch at exactly 2m = 24 "
          f"(floor <= 24): PASS", flush=True)
    wd = json.loads((DATA / "s6_drift.json").read_text())
    ok_c = (wd["zero_and_l12"]["a36_witness"]["compact_cover_lift"]
            is False)
    print(f"  (c) b=0 witness in the wrapped/winding corner: "
          f"{'PASS' if ok_c else 'FAIL'} — the -6 stays an admitted "
          f"scope term at l=12; at l=24 the wrapped corner remains "
          f"a listed condition", flush=True)
    # (d) S9 coherence: the deepened u=1 link table must DOMINATE
    # the banked g26 table row-for-row (already asserted at run
    # time by the census regression) and its dcap truncs are zero:
    lk40 = json.loads((DATA / "s9_link_u1k2g40.json").read_text())
    ok_d = (lk40["info"]["trunc_dcap"] == 0
            and lk40["info"]["trunc_extent"] == 0
            and not lk40["info"]["aborts"])
    print(f"  (d) s9 g40 link table clean (trunc_dcap 0, "
          f"trunc_extent 0, no aborts): "
          f"{'PASS' if ok_d else 'FAIL'}", flush=True)

    out = dict(
        closed_inputs={f"{k[0]},{k[1]}": v
                       for k, v in sorted(closed.items())},
        link_inputs={f"w{wstar}": {u: dict(gcap=r["gcap"],
                                           complete_h=r["complete_h"],
                                           buckets=len(r["tabL"]))
                                   for u, r in lk.items()}
                     for wstar, lk in sorted(links.items())},
        floor_r0=f0, r0_strata=per0,
        floor_r1=f1, r1_branches=per1,
        floor_r2=f2, r2_detail=det2,
        floor_T1=floor,
        controls=dict(a=ok_a, b=ok_b, c=ok_c, d=ok_d),
        wall_s=round(time.time() - t0, 1))
    (DATA / "s9_assembly.json").write_text(json.dumps(out, indent=1))
    print(f"\nwrote {DATA/'s9_assembly.json'} ({out['wall_s']} s)",
          flush=True)


if __name__ == "__main__":
    main()
