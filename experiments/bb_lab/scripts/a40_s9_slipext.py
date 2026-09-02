#!/usr/bin/env python3
"""A40 S9 — the rightward-relocation inequality: census extension
(g 32-34) + parent-logged specimen replay.

Census mode reuses the S8 SlipLinkMarch UNCHANGED (so the g26
regression stays byte-comparable); this driver only adds an output
name in the s9_* namespace and records the full cap set in params.

Replay mode: a parent-logged rerun of the same march below a small
gcap, harvesting every completed crossing of a requested
(L, gH, slip) class and re-verifying the winner end-to-end through
the independent CoverFragment (every E, slab classes, window rule,
weight, drift) — the S7 specimen discipline (§12.6)."""
from __future__ import annotations

import argparse
import json
import sys
import time
from itertools import combinations
from pathlib import Path

LAB = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(LAB / "src"))
sys.path.insert(0, str(LAB / "scripts"))

DATA = LAB / "data" / "a40"

from a40_s6_frontier import (  # noqa: E402
    wt, lsb, norm, dilate, tooth_ok, seeds_full,
)
from a40_s8_slip import SlipLinkMarch, derive_caps, _rss_mb  # noqa: E402

PRE, HEAVY, POST = 0, 1, 2


def run_census(args):
    import gc
    gc.disable()
    t0 = time.time()
    from bb_lab.tower import validate_banked
    validate_banked(LAB / "data")
    print("validate_banked: PASS", flush=True)
    out = {"caps": derive_caps()}

    print(f"s9 slip census: u={args.u} kmax={args.kmax} "
          f"whcap={args.whcap} gcap={args.gcap} hcap={args.hcap} "
          f"dcap={args.dcap} (byseed)", flush=True)
    m = SlipLinkMarch(args.u, kmax=args.kmax, whcap=args.whcap,
                      gcap=args.gcap, hcap=args.hcap,
                      dcap=args.dcap)
    sds = seeds_full(args.u)
    info = m.run_byseed(sds)
    print(f"info: {info}", flush=True)
    reg = {}
    if args.ref:
        ref = json.loads((DATA / args.ref).read_text())
        ref_pre = {tuple(map(int, k.split(","))): g
                   for k, g in ref["tab_pre"].items()}
        ref_link = {tuple(map(int, k.split(","))): g
                    for k, g in ref["tab_link"].items()}
        my_pre_ok = all(m.tab_pre.get(k) == g
                        for k, g in ref_pre.items()
                        if g <= args.gcap and k[0] <= args.hcap)
        my_link_ok = all(m.tab_link.get(k) == g
                         for k, g in ref_link.items()
                         if g <= args.gcap and k[0] <= args.hcap)
        reg = dict(ref=args.ref, pre_equal=bool(my_pre_ok),
                   link_equal=bool(my_link_ok),
                   pre_rows=len(ref_pre), link_rows=len(ref_link))
        ref_lL = {tuple(map(int, k.split(","))): g
                  for k, g in ref.get("tab_link_L", {}).items()}
        if ref_lL:
            my_lL_ok = all(m.tab_link_L.get(k) == g
                           for k, g in ref_lL.items()
                           if g <= args.gcap and k[0] <= args.hcap)
            reg["link_L_equal"] = bool(my_lL_ok)
            assert my_lL_ok, "LINK_L REGRESSION FAILED"
        # internal consistency: tab_link == min over L of tab_link_L
        coll = {}
        for (h, L, dl), g in m.tab_link_L.items():
            k = (h, dl)
            if k not in coll or coll[k] > g:
                coll[k] = g
        assert coll == m.tab_link, "tab_link_L collapse mismatch"
        reg["linkL_collapse_ok"] = True
        print(f"regression vs {args.ref}: pre {my_pre_ok}, "
              f"link {my_link_ok}, link_L "
              f"{reg.get('link_L_equal', 'n/a')}, collapse ok",
              flush=True)
        assert my_pre_ok and my_link_ok, "REGRESSION FAILED"
    # slip-table regression vs the banked g30 census: every banked
    # (L, gH, slip) class must reappear with min_g <= banked value
    # (deeper caps can only match or improve completed-crossing
    # minima; h/d caps only widen).
    g30 = json.loads((DATA / "s8_slip_u1_g30.json").read_text())
    slipreg_ok = True
    for r in g30["slip_rows"]:
        k = (r["L"], r["gH"], r["slip"])
        v = m.tab_slip.get(k)
        if v is None or v > r["min_g"]:
            slipreg_ok = False
            print(f"  SLIPREG MISS {k}: banked {r['min_g']} "
                  f"got {v}", flush=True)
    reg["slip_vs_g30"] = bool(slipreg_ok)
    print(f"slip-table regression vs g30: {slipreg_ok}", flush=True)
    assert slipreg_ok, "SLIP REGRESSION FAILED"
    out["regression"] = reg

    rows = []
    for (L, gH, slip), g in sorted(m.tab_slip.items()):
        rows.append(dict(L=L, gH=gH, slip=slip, min_g=g))
    pos = [r for r in rows if r["slip"] > 0]
    if rows:
        c_half_c0 = max(abs(r["slip"]) - 0.5 * (r["gH"] - 8)
                        for r in rows)
        mx = max(rows, key=lambda r: abs(r["slip"]))
        out["census"] = dict(
            n=len(rows), max_abs_slip=abs(mx["slip"]),
            argmax=mx, n_positive=len(pos),
            fit_c_half=dict(c=0.5, c0=round(c_half_c0, 3)))
        print(f"census: {len(rows)} (L, gH, slip) classes; max "
              f"|slip| {abs(mx['slip'])} at {mx}; POSITIVE slips: "
              f"{len(pos)} {pos if pos else '(none)'}", flush=True)
    out["slip_rows"] = rows
    # the rightward-transport wall, PRE/LINK phase (measured):
    # every (h, dlt >= +1) row of the pre and link tables
    out["pre_pos_rows"] = sorted(
        [dict(h=k[0], dlt=k[1], min_g=g)
         for k, g in m.tab_pre.items() if k[1] >= 1],
        key=lambda r: (r["dlt"], r["h"]))
    out["link_pos_rows"] = sorted(
        [dict(h=k[0], dlt=k[1], min_g=g)
         for k, g in m.tab_link.items() if k[1] >= 0],
        key=lambda r: (r["dlt"], r["h"]))
    out["pre_max_dlt"] = max((k[1] for k in m.tab_pre), default=None)
    out["link_max_dlt"] = max((k[1] for k in m.tab_link),
                              default=None)
    print(f"pre_max_dlt {out['pre_max_dlt']}; link_max_dlt "
          f"{out['link_max_dlt']}; pre rows dlt>=+1: "
          f"{len(out['pre_pos_rows'])}", flush=True)
    out["info"] = info
    out["params"] = dict(u=args.u, kmax=args.kmax,
                         whcap=args.whcap, gcap=args.gcap,
                         hcap=args.hcap, dcap=args.dcap,
                         extent_cap=m.extent_cap, smax=m.smax,
                         dil=m.dil)
    p = DATA / f"s9_slip_u{args.u}_g{args.gcap}.json"
    p.write_text(json.dumps(out, indent=1))
    print(f"wrote {p} ({round(time.time() - t0, 1)} s)", flush=True)

    if args.save_linkl:
        # assembly-consumable link-table file (s7_link schema):
        # sound ONLY as a complete-below-caps table, so demand a
        # clean run (no aborts, no extent truncs; dcap truncs are
        # reported and must be judged at consume time).
        assert not info["aborts"], "aborted run: no link file"
        assert info["trunc_extent"] == 0
        lout = dict(
            params=dict(u=args.u, kmax=args.kmax, whcap=args.whcap,
                        gcap=args.gcap, hcap=args.hcap, smax=3,
                        dil=4, dcap=args.dcap),
            info=dict(trunc_g=info["trunc_g"], trunc_extent=0,
                      trunc_dcap=info["trunc_dcap"],
                      complete_h=args.hcap, nodes=info["nodes"],
                      wall_s=info["wall_s"], aborts=[]),
            tab_pre={f"{h},{d}": g
                     for (h, d), g in m.tab_pre.items()},
            tab_link={f"{h},{d}": g
                      for (h, d), g in m.tab_link.items()},
            tab_link_L={f"{h},{L},{d}": g
                        for (h, L, d), g in m.tab_link_L.items()})
        lp = DATA / f"s9_link_u{args.u}k{args.kmax}g{args.gcap}.json"
        lp.write_text(json.dumps(lout, indent=1))
        print(f"wrote {lp}", flush=True)


# ---------------------------------------------------------------------
# specimen: parented march -> fragment-verified (L, gH, slip) witness
# ---------------------------------------------------------------------

def run_specimen(args):
    """Harvest crossings of a requested slip class and verify the
    cheapest end-to-end through CoverFragment (S7 §12.6 machinery
    reused verbatim: ParentedLinkMarch + replay).  The (L, gH,
    slip) of each candidate is recomputed FROM THE FRAGMENT (block
    location by slab weights, anchors by column minima), so the
    verdict is independent of the march's slip bookkeeping.
    Narrowing (dcap/hcap) is sound for the specimen hunt: any found
    object is fully verified (S7 precedent)."""
    sys.argv = [sys.argv[0]]
    from bb_lab.tower import validate_banked
    validate_banked(LAB / "data")
    print("validate_banked: PASS", flush=True)
    from a40_s7_tax import ParentedLinkMarch, replay
    from a40_s6_drift import CoverFragment
    want_L, want_gH, want_slip = map(int, args.want.split(","))
    sds = [w8 for w8 in seeds_full(args.u)
           if sum(wt(r) for r in w8) == args.u]
    lo, hi = 0, len(sds)
    if args.seeds:
        lo, hi = map(int, args.seeds.split(":"))
    found = []
    t0 = time.time()
    for si in range(lo, hi):
        w8 = sds[si]
        m = ParentedLinkMarch(args.u, kmax=args.kmax,
                              whcap=args.whcap, gcap=args.gcap,
                              hcap=args.hcap, dcap=args.dcap)
        info = m.run([w8], log=False)
        rss = _rss_mb()
        print(f"  seed {si}: {info['popped']} nodes, rss {rss}MB",
              flush=True)
        for key in list(m.parents):
            dyn, anch, phase, L, h, dlt = key
            if phase != POST or L != want_L:
                continue
            g = m.tab_link_L.get((h, L, dlt))
            if g is None or g > args.gcap:
                continue
            rows, chain = replay(m, key, [w8])
            lo2, hi2 = min(rows), max(rows)
            fr = CoverFragment(
                [rows[j] for j in range(lo2, hi2 + 1)], lo2)
            if not fr.admissible():
                continue
            sl = fr.slabs()
            heavy_idx = [i for i, w in enumerate(sl) if w >= 8]
            if not heavy_idx:
                continue
            a, b = heavy_idx[0], heavy_idx[-1]
            if b - a + 1 != want_L or b != a + want_L - 1:
                continue
            gH = sum(sl[a:b + 1])
            anchors = fr.anchors()
            if not (0 < a and b + 1 < len(anchors)):
                continue
            slip = anchors[b + 1] - anchors[a - 1]
            if gH == want_gH and slip == want_slip:
                found.append(dict(
                    sid=si, g=int(sum(sl) + 0),
                    slabs=[int(x) for x in sl],
                    anchors=[int(x) for x in anchors],
                    slip=int(slip), gH=int(gH), L=int(b - a + 1),
                    weight=int(fr.weight()),
                    drift=int(fr.drift()),
                    rows=[[sorted(rows[j][0]), sorted(rows[j][1])]
                          for j in range(lo2, hi2 + 1)]))
        del m
        if found and args.first:
            break
    best = min(found, key=lambda r: r["g"]) if found else None
    if best:
        print(f"SPECIMEN VERIFIED: seed {best['sid']}, g "
              f"{best['g']}, slabs {best['slabs']}, anchors "
              f"{best['anchors']}, slip {best['slip']} across "
              f"gH={best['gH']} L={best['L']}; fragment "
              f"E-admissible, weight {best['weight']}, drift "
              f"{best['drift']}", flush=True)
    else:
        print(f"NO specimen of class (L={want_L}, gH={want_gH}, "
              f"slip={want_slip}) at g <= {args.gcap} within "
              f"hcap {args.hcap} dcap {args.dcap}", flush=True)
    out = dict(params=dict(u=args.u, kmax=args.kmax,
                           whcap=args.whcap, gcap=args.gcap,
                           hcap=args.hcap, dcap=args.dcap,
                           want=args.want, seeds=[lo, hi]),
               n_found=len(found), best=best,
               wall_s=round(time.time() - t0, 1))
    p = DATA / f"s9_specimen_{args.tag}.json"
    p.write_text(json.dumps(out, indent=1))
    print(f"wrote {p}", flush=True)


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    s = sub.add_parser("census")
    s.add_argument("--u", type=int, default=1)
    s.add_argument("--kmax", type=int, default=2)
    s.add_argument("--whcap", type=int, default=14)
    s.add_argument("--gcap", type=int, required=True)
    s.add_argument("--hcap", type=int, required=True)
    s.add_argument("--dcap", type=int, required=True)
    s.add_argument("--ref", type=str, default="")
    s.add_argument("--save-linkl", action="store_true")
    r = sub.add_parser("specimen")
    r.add_argument("--u", type=int, default=1)
    r.add_argument("--kmax", type=int, default=2)
    r.add_argument("--whcap", type=int, default=14)
    r.add_argument("--gcap", type=int, required=True)
    r.add_argument("--hcap", type=int, default=10)
    r.add_argument("--dcap", type=int, default=12)
    r.add_argument("--want", type=str, required=True,
                   help="L,gH,slip")
    r.add_argument("--seeds", type=str, default="")
    r.add_argument("--first", action="store_true")
    r.add_argument("--tag", type=str, required=True)
    args = ap.parse_args()
    if args.cmd == "census":
        run_census(args)
    else:
        run_specimen(args)


if __name__ == "__main__":
    main()
