#!/usr/bin/env python3
"""A40 S10 — Stage 2: the ladder rate J' (the per-slab growth of the
certified u = 1 link deficit) — census deepening to the h = 11/12
rungs, the h-stratified readout, and the refutation-witness lanes.

UNITS.  A link of h slabs and slab-weight sum g has deficit
J = 8h - g quarters (2 per slab is the all-heavy rate).  The S9
ladder observation: min_g(h) = 25, 30, 35, 40 at h = 7..10 in the
g <= 40 census, i.e. +5 g per slab = +3 quarters of deficit per slab
("J' = 3q/slab"), and the conjecture (§14.6 item 3) extrapolates it
to every h, which would give the drift-blind member floor
2m - 3m/4 = 5m/4 (~7.5 r).

`census`: the S8 SlipLinkMarch UNCHANGED (S9 driver logic verbatim,
s10_* output names), run to a deeper g-cap so the h = 11 (and 12)
rungs become visible: +5/slab predicts min_g(11) = 45, min_g(12) =
50; any h = 11 link at g <= 44 REFUTES the rung.
`strat`: the h-stratified ladder readout of a banked link table
(min g, J, the marginal rates, tight (L, dlt) configurations).
`w7link`: the asymptotic refutation witness — a u = 1 seed, one
heavy block and a light transient INTO the W7 species coast (rate
32/7 g per slab < 5), found by a backward cost-ordered march from
the banked W7 cover lift and verified end to end (CoverFragment)."""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

LAB = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(LAB / "src"))
sys.path.insert(0, str(LAB / "scripts"))
DATA = LAB / "data" / "a40"

from a40_s6_frontier import wt, lsb, seeds_full  # noqa: E402
from a40_s8_slip import SlipLinkMarch, derive_caps, _rss_mb  # noqa: E402

PRE, HEAVY, POST = 0, 1, 2


def _load_tab(d, name):
    return {tuple(map(int, k.split(","))): g for k, g in d[name].items()}


def ladder_rows(tab_link, tab_link_L, gcap):
    """h-stratified readout: min g, J = 8h - g, marginal rates."""
    byh = {}
    for (h, dl), g in tab_link.items():
        if h not in byh or byh[h][0] > g:
            byh[h] = (g, dl)
    rows = []
    prev = None
    for h in sorted(byh):
        g, dl = byh[h]
        tight = sorted([(L, d) for (hh, L, d), gg in tab_link_L.items()
                        if hh == h and gg == g])
        J = 8 * h - g
        rows.append(dict(h=h, min_g=g, J_q=J, at_cap=(g == gcap),
                         dg=(g - prev[1]) if prev and prev[0] == h - 1
                         else None,
                         dJ_q=(J - prev[2]) if prev and prev[0] == h - 1
                         else None,
                         tight=tight))
        prev = (h, g, J)
    return rows


def run_census(args):
    import gc
    gc.disable()
    t0 = time.time()
    from bb_lab.tower import validate_banked
    validate_banked(LAB / "data")
    print("validate_banked: PASS", flush=True)
    out = {"caps": derive_caps()}
    print(f"s10 ladder census: u={args.u} kmax={args.kmax} "
          f"whcap={args.whcap} gcap={args.gcap} hcap={args.hcap} "
          f"dcap={args.dcap} (byseed, S8 engine unchanged)",
          flush=True)
    m = SlipLinkMarch(args.u, kmax=args.kmax, whcap=args.whcap,
                      gcap=args.gcap, hcap=args.hcap, dcap=args.dcap)
    sds = seeds_full(args.u)
    info = m.run_byseed(sds, rss_cap=args.rss_cap)
    print(f"info: {info}", flush=True)
    reg = {}
    if args.ref:
        ref = json.loads((DATA / args.ref).read_text())
        ref_pre = _load_tab(ref, "tab_pre")
        ref_link = _load_tab(ref, "tab_link")
        ref_lL = _load_tab(ref, "tab_link_L")
        ok_pre = all(m.tab_pre.get(k) == g for k, g in ref_pre.items()
                     if g <= args.gcap and k[0] <= args.hcap)
        ok_link = all(m.tab_link.get(k) == g
                      for k, g in ref_link.items()
                      if g <= args.gcap and k[0] <= args.hcap)
        ok_lL = all(m.tab_link_L.get(k) == g for k, g in ref_lL.items()
                    if g <= args.gcap and k[0] <= args.hcap)
        coll = {}
        for (h, L, dl), g in m.tab_link_L.items():
            k = (h, dl)
            if k not in coll or coll[k] > g:
                coll[k] = g
        reg = dict(ref=args.ref, pre_equal=bool(ok_pre),
                   link_equal=bool(ok_link), link_L_equal=bool(ok_lL),
                   linkL_collapse_ok=(coll == m.tab_link))
        print(f"regression vs {args.ref}: {reg}", flush=True)
        assert ok_pre and ok_link and ok_lL and coll == m.tab_link, \
            "REGRESSION FAILED"
    out["regression"] = reg
    rows = [dict(L=L, gH=gH, slip=s, min_g=g)
            for (L, gH, s), g in sorted(m.tab_slip.items())]
    pos = [r for r in rows if r["slip"] > 0]
    out["slip_rows"] = rows
    out["n_positive_slips"] = len(pos)
    print(f"slip classes {len(rows)}; POSITIVE slips: "
          f"{len(pos)} {pos if pos else '(none)'}", flush=True)
    out["pre_max_dlt"] = max((k[1] for k in m.tab_pre), default=None)
    out["link_max_dlt"] = max((k[1] for k in m.tab_link), default=None)
    lad = ladder_rows(m.tab_link, m.tab_link_L, args.gcap)
    out["ladder"] = lad
    print("LADDER (h, min_g, J_q, dg, dJ_q, at_cap):", flush=True)
    for r in lad:
        print(f"  h={r['h']:2d} min_g={r['min_g']:2d} J={r['J_q']:3d}q "
              f"dg={r['dg']} dJ={r['dJ_q']} cap={r['at_cap']} "
              f"tight={r['tight']}", flush=True)
    out["info"] = info
    out["params"] = dict(u=args.u, kmax=args.kmax, whcap=args.whcap,
                         gcap=args.gcap, hcap=args.hcap, dcap=args.dcap,
                         extent_cap=m.extent_cap, smax=m.smax, dil=m.dil)
    p = DATA / f"s10_slip_u{args.u}_g{args.gcap}.json"
    p.write_text(json.dumps(out, indent=1))
    print(f"wrote {p} ({round(time.time() - t0, 1)} s)", flush=True)
    if args.save_linkl:
        assert not info["aborts"], "aborted run: no link file"
        assert info["trunc_extent"] == 0
        lout = dict(
            params=dict(u=args.u, kmax=args.kmax, whcap=args.whcap,
                        gcap=args.gcap, hcap=args.hcap, smax=3, dil=4,
                        dcap=args.dcap),
            info=dict(trunc_g=info["trunc_g"], trunc_extent=0,
                      trunc_dcap=info["trunc_dcap"],
                      complete_h=args.hcap, nodes=info["nodes"],
                      wall_s=info["wall_s"], aborts=[]),
            tab_pre={f"{h},{d}": g for (h, d), g in m.tab_pre.items()},
            tab_link={f"{h},{d}": g
                      for (h, d), g in m.tab_link.items()},
            tab_link_L={f"{h},{L},{d}": g
                        for (h, L, d), g in m.tab_link_L.items()})
        lp = DATA / f"s10_link_u{args.u}k{args.kmax}g{args.gcap}.json"
        lp.write_text(json.dumps(lout, indent=1))
        print(f"wrote {lp}", flush=True)


def run_strat(args):
    from bb_lab.tower import validate_banked
    validate_banked(LAB / "data")
    print("validate_banked: PASS", flush=True)
    d = json.loads((DATA / args.table).read_text())
    gcap = d["params"]["gcap"]
    lad = ladder_rows(_load_tab(d, "tab_link"), _load_tab(d, "tab_link_L"),
                      gcap)
    pre = _load_tab(d, "tab_pre")
    prem = {}
    for (h, dl), g in pre.items():
        if h not in prem or prem[h] > g:
            prem[h] = g
    print(f"table {args.table} (gcap {gcap})", flush=True)
    print("PRE min_g by h:", dict(sorted(prem.items())), flush=True)
    for r in lad:
        print(f"  h={r['h']:2d} min_g={r['min_g']:2d} J={r['J_q']:3d}q "
              f"dg={r['dg']} dJ={r['dJ_q']} cap={r['at_cap']} "
              f"tight={r['tight']}", flush=True)
    # the affine fits: J <= 3h + C over the certified rows
    C3 = max(r["J_q"] - 3 * r["h"] for r in lad)
    print(f"J(h) <= 3h + {C3} quarters over the certified rows "
          f"(g(h) >= 5h - {C3 + 0})", flush=True)
    out = dict(table=args.table, gcap=gcap, ladder=lad,
               pre_min_by_h=prem, fit_J_le_3h_plus=C3)
    p = DATA / f"s10_ladder_strat_{Path(args.table).stem}.json"
    p.write_text(json.dumps(out, indent=1))
    print(f"wrote {p}", flush=True)


def run_bwdprobe(args):
    """Does the u = 1 BACKWARD light tree (the post piece of a
    u1 -> u1 link seen from the next seed) exhaust like the forward
    one (§12.4)?  The S6 March, bwd direction, u = 1, run to a
    deeper cap; per-h minima and the truncation flags are the
    readout.  Output in the s10 namespace (the S6 tables untouched)."""
    from a40_s6_frontier import March
    from bb_lab.tower import validate_banked
    validate_banked(LAB / "data")
    print("validate_banked: PASS", flush=True)
    t0 = time.time()
    sds = seeds_full(1)
    res = {}
    for direction in args.dirs.split(","):
        m = March(1, direction, smax=3, dil=4, hcap=args.hcap,
                  gcap=args.gcap)
        info = m.run(sds, log=True, max_nodes=args.maxnodes)
        byh = {}
        for (h, d), g in m.table.items():
            if h not in byh or byh[h][0] > g:
                byh[h] = (g, d)
        hmax = max(byh) if byh else 0
        print(f"{direction} u=1 gcap {args.gcap} hcap {args.hcap}: {info};"
              f" min g by h: {dict(sorted(byh.items()))}; max h reached "
              f"{hmax}; EXHAUSTED: {not info['trunc_g'] and not info['trunc_nodes'] and hmax < args.hcap}",
              flush=True)
        res[direction] = dict(info=info,
                              min_g_by_h={h: dict(g=g, dlt=d)
                                          for h, (g, d) in sorted(byh.items())},
                              hmax=hmax,
                              exhausted=(not info["trunc_g"]
                                         and not info["trunc_nodes"]
                                         and hmax < args.hcap))
    res["params"] = dict(gcap=args.gcap, hcap=args.hcap,
                         maxnodes=args.maxnodes, dirs=args.dirs)
    res["wall_s"] = round(time.time() - t0, 1)
    p = DATA / f"s10_bwdprobe_g{args.gcap}.json"
    p.write_text(json.dumps(res, indent=1))
    print(f"wrote {p}", flush=True)


def run_w7fwd(args):
    """The asymptotic-ladder witness, FORWARD form: the S8
    SlipLinkMarch (layered, streaming, unchanged) in stratum u = 2
    seeded at the W7 species' own minimum window (weight 2).  Any
    completed crossing (POST state at (h0, dlt) with cost g0) is a
    u = 2 link [W7 seed -> light -> block -> post]; prepending k
    further W7 periods keeps it a link (the seed stays the run's
    first weight-2 window) with h = h0 + 7k, g = g0 + 32k, so
    J(h) grows by 24 quarters per 7 slabs = 3.43 q/slab > 3: the
    per-slab ladder extrapolation fails in the u = 2 table for every
    h beyond the first crossing.  The found crossing is replayed and
    verified (ParentedLinkMarch + CoverFragment) before any claim."""
    import gc
    gc.disable()
    t0 = time.time()
    from bb_lab.tower import validate_banked
    validate_banked(LAB / "data")
    print("validate_banked: PASS", flush=True)
    from a40_s10_w7link import w7_lift, rows_to_masks
    from a40_s6_frontier import norm
    fr = w7_lift(args.l, n_periods=6)
    sl = fr.slabs()
    OFF = 64
    mrows = rows_to_masks(fr.rows, OFF)
    # the weight-2 windows: slab index j (rows j..j+3 of the lift)
    seeds = []
    for j, w in enumerate(sl):
        if w == 2:
            w8 = (mrows[j][0], mrows[j + 1][0], mrows[j + 2][0],
                  mrows[j + 3][0], mrows[j][1], mrows[j + 1][1],
                  mrows[j + 2][1], mrows[j + 3][1])
            seeds.append(norm(w8)[0])
    seeds = sorted(set(seeds))
    print(f"W7 lift l={args.l}: slabs {sl[:14]}...; {len(seeds)} distinct "
          f"weight-2 windows (normalized): {seeds}", flush=True)
    m = SlipLinkMarch(2, kmax=args.kmax, whcap=args.whcap, gcap=args.gcap,
                      hcap=args.hcap, dcap=args.dcap)
    info = m.run_byseed(seeds, rss_cap=args.rss_cap)
    print(f"info: {info}", flush=True)
    lad = ladder_rows(m.tab_link, m.tab_link_L, args.gcap)
    prem = {}
    for (h, dl), g in m.tab_pre.items():
        if h not in prem or prem[h] > g:
            prem[h] = g
    print(f"PRE min_g by h (from the W7 window): {dict(sorted(prem.items()))}",
          flush=True)
    print("LINK ladder from the W7 seed (h, min_g, J_q):", flush=True)
    for r in lad:
        print(f"  h={r['h']:2d} min_g={r['min_g']:2d} J={r['J_q']:3d}q "
              f"dg={r['dg']} cap={r['at_cap']} tight={r['tight']}",
              flush=True)
    out = dict(params=dict(u=2, kmax=args.kmax, whcap=args.whcap,
                           gcap=args.gcap, hcap=args.hcap, dcap=args.dcap,
                           l=args.l, seeds=[list(s) for s in seeds]),
               info=info, pre_min_by_h=prem, ladder=lad,
               tab_link={f"{h},{d}": g for (h, d), g in m.tab_link.items()},
               tab_link_L={f"{h},{L},{d}": g
                           for (h, L, d), g in m.tab_link_L.items()},
               w7_period=dict(dh=7, dg=32, dJ_q=24))
    best = min(lad, key=lambda r: r["min_g"]) if lad else None
    if best:
        print(f"\ncheapest crossing from the W7 window: h={best['h']} "
              f"g={best['min_g']} J={best['J_q']}q; family h+7k, g+32k, "
              f"J+24k q => asymptotic J' >= 24/7 = 3.43 q/slab in the "
              f"u=2 table", flush=True)
        # replay + verify through the S7 parented march (heap-based,
        # small at this cap) — the crossing must be E-admissible with
        # the slab profile the table claims
        from a40_s7_tax import ParentedLinkMarch, replay
        from a40_s6_drift import CoverFragment
        pm_ = ParentedLinkMarch(2, kmax=args.kmax, whcap=args.whcap,
                                gcap=best["min_g"], hcap=best["h"] + 1,
                                dcap=args.dcap)
        pm_.run(seeds, log=False)
        ver = None
        for key in list(pm_.parents):
            dyn, anch, phase, L, h, dlt = key
            if phase != POST or h != best["h"]:
                continue
            g = pm_.tab_link_L.get((h, L, dlt))
            if g != best["min_g"]:
                continue
            rows_d, chain = replay(pm_, key, seeds)
            lo2, hi2 = min(rows_d), max(rows_d)
            rws = [rows_d[j] for j in range(lo2, hi2 + 1)]
            frag = CoverFragment(rws, lo2)
            if not frag.admissible():
                continue
            sl2 = frag.slabs()
            ver = dict(slabs=sl2, g=sum(sl2), h=len(sl2),
                       anchors=frag.anchors(), weight=frag.weight(),
                       rows=[[sorted(a), sorted(b)] for a, b in rws],
                       heavy=[i for i, w in enumerate(sl2) if w >= 8])
            break
        out["verified_crossing"] = ver
        print(f"verified crossing: {ver}", flush=True)
    out["wall_s"] = round(time.time() - t0, 1)
    p = DATA / f"s10_w7fwd_g{args.gcap}.json"
    p.write_text(json.dumps(out, indent=1))
    print(f"wrote {p}", flush=True)


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    f = sub.add_parser("w7fwd")
    f.add_argument("--gcap", type=int, default=40)
    f.add_argument("--hcap", type=int, default=14)
    f.add_argument("--dcap", type=int, default=30)
    f.add_argument("--kmax", type=int, default=2)
    f.add_argument("--whcap", type=int, default=14)
    f.add_argument("--l", type=int, default=24)
    f.add_argument("--rss-cap", type=int, default=2500)
    b = sub.add_parser("bwdprobe")
    b.add_argument("--gcap", type=int, default=28)
    b.add_argument("--hcap", type=int, default=16)
    b.add_argument("--maxnodes", type=int, default=12_000_000)
    b.add_argument("--dirs", type=str, default="bwd")
    s = sub.add_parser("census")
    s.add_argument("--u", type=int, default=1)
    s.add_argument("--kmax", type=int, default=2)
    s.add_argument("--whcap", type=int, default=14)
    s.add_argument("--gcap", type=int, required=True)
    s.add_argument("--hcap", type=int, default=19)
    s.add_argument("--dcap", type=int, default=30)
    s.add_argument("--rss-cap", type=int, default=2500)
    s.add_argument("--ref", type=str, default="")
    s.add_argument("--save-linkl", action="store_true")
    t = sub.add_parser("strat")
    t.add_argument("--table", type=str, default="s9_link_u1k2g40.json")
    w = sub.add_parser("w7link")
    w.add_argument("--gcap", type=int, default=60)
    w.add_argument("--l", type=int, default=24)
    w.add_argument("--rss-cap", type=int, default=2500)
    w.add_argument("--node-cap", type=int, default=30_000_000)
    ap.add_argument("--log", type=str, default="",
                    help="tee stdout/stderr into this file (repo data "
                         "dir; shell redirection is unavailable here)")
    args = ap.parse_args()
    if args.log:
        fh = open(args.log, "a", buffering=1)
        sys.stdout = sys.stderr = fh
    if args.cmd == "census":
        run_census(args)
    elif args.cmd == "strat":
        run_strat(args)
    elif args.cmd == "bwdprobe":
        run_bwdprobe(args)
    elif args.cmd == "w7fwd":
        run_w7fwd(args)
    else:
        from a40_s10_w7link import run_w7link
        run_w7link(args)


if __name__ == "__main__":
    main()
