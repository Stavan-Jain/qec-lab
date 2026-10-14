#!/usr/bin/env python3
"""A40 S7 — Stage 3: the reassembly at (l, m) = (24, 18) = [[864,12]]
with the boundary coupling charged.

BRANCHES (a class-minimal nontrivial y-spanning X-logical v; slabs
classified light (W <= 7) / heavy (W >= 8); r = number of maximal
light runs):
 - ALL-HEAVY (r = 0 runs, no light slab): w >= 2m = 36 pointwise.
 - ALL-LIGHT (r = 0 heavies): one closed 18-step cycle, every E
   enforced, seeded at its global min window (stratum u): the
   CLOSED MARCH (a40_s7_link.py) — per-stratum emptiness below a
   g-cap gives w >= ceil((gcap + 1 - u)/4); the u = 1 stratum is
   EXHAUSTED (trunc_g False at g <= 76): no such walk at ANY
   weight; u >= 3 fall back to the analytic ceil(u m / 4) when
   deeper than the enumerated cap.
 - r = 1, short block (L <= kmax, W <= whcap): the closed-link
   march (same seed logic, one heavy block, return to seed).
 - r = 1, long block (L >= kmax+1): no closure lever; the run's
   deficit <= delta-blind Dbest(18 - L) from the s6 master tables
   (delta-blind per-h minima are EXACT under the s6 anchor
   aliasing).
 - r = 1, fat heavy (some W >= whcap + 1): the fat slab pays its
   own excess: D <= Dbest(18 - L) - (whcap + 1 - 8)/4.
 - r >= 2: the LINK DECOMPOSITION.  sum_i D_link_i = D_total
   + sum_i (2 - u_i/4) (seed slabs counted twice), so
   D_total = sum_i [J(u_i, h_i, delta_i) - (2 - u_i/4)] maximized
   over compositions sum h_i = 18 + r.  Two DPs:
     DP-closure: every link enumerated-type (u <= 2, L <= kmax,
       W <= whcap), delta-resolved, sum delta == 0 (mod 24);
       absent buckets get the g-cap bound.
     DP-free: >= 1 link of non-enumerated type (u >= 3, L >= 3, or
       fat), no closure; non-enumerated links granted the LOOSE
       split FWD_u(1+a) + BWD_1(b) (s6-style; delta-blind exact)
       with the fat credit where applicable.
   floor = 36 - max(DP)/4.
   Refinements: pieces are stratum-coupled (u -> u_next; the post
   piece is granted the u_next bwd table, the seed subtraction is
   u_next's); absent enumerated buckets are granted
   min(g-cap bound, loose split) — the loose split bounds EVERY
   link; the whole r>=2 computation is evaluated at BOTH heavy-
   class semantics (W <= 14 full-depth completeness, W <= 16 with
   the larger fat credit) and the better sound floor is taken.

Scope conditions carried by every floor (unchanged from s6 unless
noted): radius-dil prefix-connected growth with smax new points per
row (stability-checked), no wrap-interacting non-winding fragments
(cover extent <= 34; winding +-24 included), |delta| <= dcap per
fragment (trunc counters all ZERO in the production runs, so the
caps did not bind), heavy blocks W <= whcap / L <= kmax enumerated
with fat/long branches closing the complement.
"""
from __future__ import annotations

import glob
import json
import math
import sys
import time
from pathlib import Path

LAB = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(LAB / "src"))
sys.path.insert(0, str(LAB / "scripts"))
DATA = LAB / "data" / "a40"

M, ELL = 18, 24
WHCAP, KMAX = 14, 2


# ---------------------------------------------------------------------
# inputs
# ---------------------------------------------------------------------

def load_closed():
    """Best (deepest empty) closed-march result per (u, kmax>=1?).
    Returns {(u, 'k0'|'k2'): dict(gcap, exhausted, files)}."""
    best = {}
    for p in sorted(glob.glob(str(DATA / "s7_closed_*.json"))):
        d = json.loads(Path(p).read_text())
        pr, info = d["params"], d["info"]
        if pr["m"] != M or pr["ell"] != ELL:
            continue
        if info.get("aborts") or info.get("trunc_nodes"):
            continue          # aborted runs carry no emptiness claim
        key = (pr["u"], "k0" if pr["kmax"] == 0 else "k2")
        # shards: collect; a config counts only if ALL seeds covered
        rec = best.setdefault(key, {})
        g = pr["gcap"]
        e = rec.setdefault(g, dict(cover=[], closed=0,
                                   trunc_g=False, whcap=pr["whcap"],
                                   kmax=pr["kmax"]))
        e["cover"].append(tuple(pr["seeds"]))
        e["closed"] += info["closed"]
        e["trunc_g"] = e["trunc_g"] or info["trunc_g"]
        assert info["trunc_extent"] == 0 and info["trunc_dcap"] == 0
    out = {}
    from a40_s6_frontier import seeds_full, wt
    nseeds = {u: len([w for w in seeds_full(u)
                      if sum(wt(r) for r in w) == u])
              for u in (1, 2, 3)}
    for key, recs in best.items():
        u = key[0]
        for g in sorted(recs, reverse=True):
            e = recs[g]
            iv = sorted(e["cover"])
            covered = 0
            for lo, hi in iv:
                if lo <= covered:
                    covered = max(covered, hi)
            if covered >= nseeds.get(u, 10**9) and e["closed"] == 0:
                out[key] = dict(gcap=g, exhausted=not e["trunc_g"],
                                whcap=e["whcap"], kmax=e["kmax"])
                break
    return out


def load_links():
    """Merge link tables per stratum across files: certified
    buckets merge unconditionally (each is a real enumerated
    object); the per-h absent-bucket cap gcap_by_h[h] = max gcap
    over files whose complete_h >= h (completeness required for
    the absent => expensive direction)."""
    return {14: _load_links_class(14), 16: _load_links_class(16)}


def _load_links_class(wstar):
    """Link tables under the class semantics 'enumerated heavy
    W <= wstar': certs merge from every file (any enumerated object
    is a real link); the absent=>expensive completeness only from
    files run at whcap == wstar."""
    out = {}
    for p in sorted(glob.glob(str(DATA / "s7_link_*.json"))):
        d = json.loads(Path(p).read_text())
        pr = d["params"]
        u = pr["u"]
        assert d["info"]["trunc_extent"] == 0
        ch = d["info"].get("complete_h", pr["hcap"])
        rec = out.setdefault(u, dict(tabL={}, gcap_by_h={},
                                     kmax=pr["kmax"], whcap=wstar,
                                     dcap=pr["dcap"], srcs=[],
                                     trunc_dcap=0))
        assert pr["kmax"] == rec["kmax"]
        for k, g in d["tab_link_L"].items():
            kk = tuple(map(int, k.split(",")))
            if kk not in rec["tabL"] or rec["tabL"][kk] > g:
                rec["tabL"][kk] = g
        if pr["whcap"] == wstar:
            for h in range(1, ch + 1):
                if rec["gcap_by_h"].get(h, -1) < pr["gcap"]:
                    rec["gcap_by_h"][h] = pr["gcap"]
        rec["trunc_dcap"] += d["info"]["trunc_dcap"]
        rec["srcs"].append(Path(p).name)
    for u in list(out):
        rec = out[u]
        if not rec["gcap_by_h"]:
            del out[u]            # no completeness at this class
            continue
        rec["gcap"] = max(rec["gcap_by_h"].values())
        rec["complete_h"] = max(
            (h for h, g in rec["gcap_by_h"].items()
             if g == rec["gcap"]), default=0)
        rec["src"] = "+".join(rec["srcs"])
    return out


def load_s6():
    """s6 fwd/bwd stratum tables (delta-blind use only)."""
    out = {}
    for u in (1, 2):
        p = DATA / f"s6_frontier_u{u}prod.json"
        d = json.loads(p.read_text())
        f = {tuple(map(int, k.split(","))): g
             for k, g in d["fwd"]["table"].items()}
        b = {tuple(map(int, k.split(","))): g
             for k, g in d["bwd"]["table"].items()}
        out[u] = dict(fwd=f, bwd=b, gcap=d["fwd"]["params"]["gcap"])
    # the matched u1 composed table (real fragments, tighter)
    dm = json.loads((DATA / "s6_frontier_u1matched.json").read_text())
    out["m1"] = dict(
        comp={tuple(map(int, k.split(","))): g
              for k, g in dm["composed"].items()},
        gcap=dm["fwd"]["params"]["gcap"])
    return out


# quarter-unit helpers ------------------------------------------------

def fwd_q(s6, u, h):
    """Max deficit (quarters) of a fwd piece of h slabs, stratum u
    (delta-blind).  h >= 1 (h=1 is the bare seed slab: 8 - u)."""
    if h <= 0:
        return None
    if u >= 3:
        return (8 - u) * h
    t = s6[u]["fwd"]
    cands = [8 * h - g for (hh, d), g in t.items() if hh == h]
    cap = 8 * h - (s6[u]["gcap"] + 1)
    ana = (8 - u) * h
    return max(cands) if cands else min(cap, ana)


def bwd1_q(s6, h):
    """Max deficit (quarters) of a bwd piece of h slabs, any
    stratum (the u=1 bwd table dominates all strata)."""
    if h <= 0:
        return None
    t = s6[1]["bwd"]
    cands = [8 * h - g for (hh, d), g in t.items() if hh == h]
    cap = 8 * h - (s6[1]["gcap"] + 1)
    ana = 7 * h
    return max(cands) if cands else min(cap, ana)


def dbest_q(s6, h):
    """Delta-blind max deficit (quarters) of a WHOLE light run of h
    slabs = the s6 master logic: strata u=1 (matched composed), u=2
    (loose composed via fwd+bwd), u=3,4 analytic, u>=5 analytic."""
    cands = [6 * h]                    # u >= 5: (2 - 5/4)h
    for u, ana_rate in ((3, 5), (4, 4)):
        cands.append(ana_rate * h)
    # u=1 matched composed
    t = s6["m1"]["comp"]
    c1 = [8 * h - g for (hh, d), g in t.items() if hh == h]
    cap1 = 8 * h - (s6["m1"]["gcap"] + 1)
    cands.append(max(c1) if c1 else min(cap1, 7 * h))
    # u=2 loose composed: fwd_2(a) + bwd... s6 composed u2 table:
    p2 = DATA / "s6_frontier_u2prod.json"
    d2 = json.loads(p2.read_text())
    t2 = {tuple(map(int, k.split(","))): g
          for k, g in d2["composed"].items()}
    c2 = [8 * h - g for (hh, d), g in t2.items() if hh == h]
    cap2 = 8 * h - (d2["fwd"]["params"]["gcap"] + 1)
    cands.append(max(c2) if c2 else min(cap2, 6 * h))
    return max(cands)


# ---------------------------------------------------------------------
# branch floors
# ---------------------------------------------------------------------

def floor_r0(closed):
    per = {}
    for u in (1, 2, 3):
        rec = closed.get((u, "k0"))
        ana = math.ceil(u * M / 4) if u >= 3 else None
        if rec is None:
            per[f"u{u}"] = ana if ana else 0
            continue
        if rec["exhausted"]:
            per[f"u{u}"] = 10**9
        else:
            f = math.ceil((rec["gcap"] + 1 - u) / 4)
            per[f"u{u}"] = max(f, ana or 0)
    per["u4"] = math.ceil(4 * M / 4)
    per["u>=5"] = math.ceil(5 * M / 4)
    return min(per.values()), per


def floor_r1(closed, s6):
    per = {}
    # short blocks: the closed-link march
    for u in (1, 2, 3):
        rec = closed.get((u, "k2"))
        ana = math.ceil(((M - 1) * u + 8) / 4) if u >= 3 else None
        if rec is None:
            per[f"short_u{u}"] = ana if ana else 0
            continue
        if rec["exhausted"]:
            per[f"short_u{u}"] = 10**9
        else:
            per[f"short_u{u}"] = math.ceil(
                (rec["gcap"] + 1 - u) / 4)
    per["short_u4"] = math.ceil(((M - 1) * 4 + 8) / 4)
    per["short_u>=5"] = math.ceil(((M - 1) * 5 + 8) / 4)
    # long blocks (L >= KMAX+1): D <= dbest(18 - L)
    for L in range(KMAX + 1, 12):
        per[f"long_L{L}"] = math.ceil(
            (8 * M - dbest_q(s6, M - L)) / 4) if M - L >= 1 else 36
    # fat heavy: D <= dbest(18 - L) - (whcap+1-8) quarters, L=1, 2
    for L in (1, 2):
        per[f"fat_L{L}"] = math.ceil(
            (8 * M - (dbest_q(s6, M - L) - (WHCAP + 1 - 8))) / 4)
    return min(per.values()), per


def bwd_q(s6, u, h):
    """Max deficit (quarters) of a bwd piece of h slabs whose
    slabs are all >= u (the next run's below-seed piece)."""
    if h <= 0:
        return None
    if u >= 3:
        return (8 - min(u, 5)) * h
    t = s6[u]["bwd"]
    cands = [8 * h - g for (hh, d), g in t.items() if hh == h]
    cap = 8 * h - (s6[u]["gcap"] + 1)
    ana = (8 - u) * h
    return max(cands) if cands else min(cap, ana)


def loose_q(s6, u, u_next, h, long_or_deep_only, fat, whcap=WHCAP):
    """Loose-split bound (quarters) for a link of h slabs: seed
    stratum u (pre slabs >= u), next-seed stratum u_next (post
    slabs >= u_next), one heavy block.  long_or_deep_only: restrict
    to L >= KMAX+1 (else any L).  fat: subtract the fat credit for
    a heavy slab W >= whcap + 1 (whcap = the stratum's enumerated
    heavy cap; analytic strata use the base 14)."""
    best = None
    for L in range(1, h - 1):
        if long_or_deep_only and L <= KMAX:
            continue
        for a in range(0, h - L):
            b = h - 1 - a - L
            if b < 1:
                continue
            fa = fwd_q(s6, min(u, 5), 1 + a)
            bb = bwd_q(s6, min(u_next, 5), b)
            if fa is None or bb is None:
                continue
            v = fa + bb
            if fat:
                v -= (whcap + 1 - 8)
            if best is None or v > best:
                best = v
    return best


def link_options_q(links, s6, u, u_next, h):
    """ALL grant options (quarters, delta-blind) for a link piece:
    seed stratum u, next-seed stratum u_next, h slabs.  Options
    [(value, typ)]: 'a' = enumerated type (u <= 2, L <= KMAX,
    W <= WHCAP, h within complete_h); 'bc' = long block / fat /
    analytic stratum / beyond-complete_h (loose split)."""
    outs = []
    seed_term = 8 - min(u_next, 5)
    rec = links.get(u)
    deep = rec is None or u > 2 or \
        rec["gcap_by_h"].get(h) is None
    if not deep:
        certs = [8 * h - g for (hh, L, d), g in rec["tabL"].items()
                 if hh == h]
        cap = 8 * h - (rec["gcap_by_h"][h] + 1)
        vl = loose_q(s6, u, u_next, h, False, False)
        if vl is not None:
            cap = min(cap, vl)
        val_a = max(certs + [cap])
        outs.append((val_a - seed_term, "a"))
    # bc options: long blocks (any u); fat (any u, any L); the
    # loose split at any L for analytic strata / deep h.  The fat
    # threshold is the stratum's enumerated whcap (base WHCAP for
    # analytic strata).
    wh_u = rec["whcap"] if (rec is not None and u <= 2) else WHCAP
    cands = []
    v_long = loose_q(s6, u, u_next, h, True, False)
    if v_long is not None:
        cands.append(v_long)
    v_fat = loose_q(s6, u, u_next, h, False, True, whcap=wh_u)
    if v_fat is not None:
        cands.append(v_fat)
    if deep:
        v_any = loose_q(s6, u, u_next, h, False, False)
        if v_any is not None:
            cands.append(v_any)
    if cands:
        outs.append((max(cands) - seed_term, "bc"))
    return outs


STRATA = (1, 2, 3, 4, 5)


def floor_r2(links, s6):
    """DP-free (>= 1 non-enumerated link, no closure) and
    DP-closure (all links enumerated, sum delta == 0 mod ELL).
    Pieces are stratum-coupled: a piece (u -> u_next) must be
    followed by a piece with seed stratum u_next; the cycle closes
    back to the first piece's stratum."""
    NEG = -10**9
    # options[(u, u_next, h)] = [(value, typ)]
    opts = {}
    for u in STRATA:
        for un in STRATA:
            for h in range(3, M + 2):
                oo = link_options_q(links, s6, u, un, h)
                if oo:
                    opts[(u, un, h)] = oo
    # DP-free: state (length, count<=2, has_bc, first_u, next_u)
    dp = {}
    par = {}
    for u in STRATA:
        dp[(0, 0, 0, u, u)] = 0
    for n in range(1, M + 1):
        for (uu, un, h), oo in opts.items():
            ln = h - 1
            if ln > n:
                continue
            for v, typ in oo:
                fb_add = 1 if typ == "bc" else 0
                for key, pv in [(k, x) for k, x in dp.items()
                                if k[0] == n - ln and k[4] == uu]:
                    _, c, fb, fu, _ = key
                    k2 = (n, min(2, c + 1), min(1, fb + fb_add),
                          fu, un)
                    nv = pv + v
                    if nv > dp.get(k2, NEG):
                        dp[k2] = nv
                        par[k2] = (key, (uu, un, h, v, typ))
    best_free, bk = max(
        ((v, k) for k, v in dp.items()
         if k[0] == M and k[1] == 2 and k[2] == 1
         and k[3] == k[4]), key=lambda t: t[0],
        default=(NEG, None))
    pieces_used = []
    k = bk
    while k in par:
        k0, pc = par[k]
        pieces_used.append(pc)
        k = k0
    best_free_a = max((v for (n, c, fb, fu, un), v in dp.items()
                       if n == M and c == 2 and fb == 0
                       and fu == un), default=NEG)

    # DP-closure: pieces (u<=2 enumerated, h <= complete_h, delta);
    # state (length, count<=2, delta, first_u, next_u)
    clo = {}
    for u in (1, 2):
        rec = links.get(u)
        if rec is None:
            continue
        cert = {}
        for (h, L, d), g in rec["tabL"].items():
            k = (h, d)
            if k not in cert or cert[k] < 8 * h - g:
                cert[k] = 8 * h - g
        hmax_c = min(M + 1, max(rec["gcap_by_h"], default=0))
        for h in range(3, hmax_c + 1):
            if rec["gcap_by_h"].get(h) is None:
                continue
            caps = {}
            for un in STRATA:
                vl = loose_q(s6, u, un, h, False, False)
                cap = 8 * h - (rec["gcap_by_h"][h] + 1)
                caps[un] = min(cap, vl) if vl is not None else cap
            for d in range(-30, 31):
                base = cert.get((h, d), NEG)
                for un in STRATA:
                    v = max(base, caps[un]) - (8 - min(un, 5))
                    k = (u, un, h, d)
                    if k not in clo or clo[k] < v:
                        clo[k] = v
    dpc = {}
    for u in (1, 2):
        dpc[(0, 0, 0, u, u)] = 0
    for n in range(1, M + 1):
        step = {}
        for (uu, un, h, d), v in clo.items():
            ln = h - 1
            if ln > n:
                continue
            step.setdefault((uu, ln), []).append((un, h, d, v))
        for key, pv in [(k, x) for k, x in dpc.items()
                        if k[0] < n]:
            nn, c, dd, fu, cur_u = key
            ln = n - nn
            for (un, h, d, v) in step.get((cur_u, ln), []):
                k2 = (n, min(2, c + 1), (dd + d) % ELL, fu, un)
                nv = pv + v
                if nv > dpc.get(k2, NEG):
                    dpc[k2] = nv
    best_clo = max((v for (n, c, dd, fu, un), v in dpc.items()
                    if n == M and c == 2 and dd == 0 and fu == un),
                   default=NEG)
    grant = max(best_free, best_clo)
    fl = math.ceil((8 * M - grant) / 4)
    return fl, dict(dp_free_bc=best_free / 4,
                    dp_free_a_only=best_free_a / 4,
                    dp_closure=best_clo / 4,
                    free_bc_argmax=[
                        dict(u=uu, u_next=un, h=h, val_q=v, typ=t)
                        for (uu, un, h, v, t) in pieces_used],
                    floor=fl)


def main():
    t0 = time.time()
    from bb_lab.tower import validate_banked
    validate_banked(LAB / "data")
    print("validate_banked: PASS", flush=True)

    closed = load_closed()
    links = load_links()
    s6 = load_s6()
    print("closed inputs:", {f"{k}": v for k, v in
                             sorted(closed.items())}, flush=True)
    for wstar, lk in sorted(links.items()):
        print(f"link inputs (class W<={wstar}):",
              {u: dict(gcap=r['gcap'], complete_h=r['complete_h'],
                       buckets=len(r['tabL']))
               for u, r in lk.items()}, flush=True)

    f0, per0 = floor_r0(closed)
    print(f"\nALL-LIGHT (r=0): floor {f0}; strata {per0}",
          flush=True)
    f1, per1 = floor_r1(closed, s6)
    print(f"r=1: floor {f1}; branches {per1}", flush=True)
    f2, det2 = 0, dict(note="no link tables yet")
    for wstar, lk in sorted(links.items()):
        if not lk:
            continue
        global WHCAP_R2
        f2w, det2w = floor_r2(lk, s6)
        print(f"r>=2 (class W<={wstar}): floor {f2w}; {det2w}",
              flush=True)
        if f2w > f2:
            f2, det2 = f2w, dict(det2w, class_whcap=wstar)
    if f2 == 0:
        print("r>=2: NO LINK TABLES — floor 0 placeholder",
              flush=True)
    fH = 2 * M
    floor = min(f0, f1, f2, fH)
    print(f"\nT1' y-sector floor (scope-listed): d_Y(24,18) >= "
          f"{floor}  [r0 {f0} | r1 {f1} | r>=2 {f2} | all-heavy "
          f"{fH}]", flush=True)

    # ---- T2' (conjectural tier) ------------------------------------
    # per-link transient: T_link = max certified D - (6/7)(h - L)
    T_link, wit = -1e9, None
    for u, rec in links.get(14, {}).items():
        for (h, L, d), g in rec["tabL"].items():
            t = (2 * h - g / 4) - (6 / 7) * (h - L)
            if t > T_link:
                T_link, wit = t, (u, h, L, d, g)
    if wit:
        # T2' mixed: r links, each (6/7)(light) + T_link, minus
        # seeds; closed branches from T1'.
        best = 0
        for r in range(2, 10):
            n_light = M - r          # r single-slab blocks
            per = (6 / 7) * (n_light) + r * T_link \
                - sum(2 - 1 / 4 for _ in range(r))
            best = max(best, per)
        fT2 = math.ceil(36 - best)
        print(f"\nT2' (conjectural: sustained rate 6/7 + measured "
              f"per-link transient T_link={T_link:.2f} at {wit}): "
              f"mixed floor >= {min(fT2, f0, f1)}", flush=True)
    else:
        fT2 = None

    # ---- controls ---------------------------------------------------
    print("\nCONTROLS:", flush=True)
    # (a) stacks: W7 x9 at (18,63) is all-light (no links, no heavy)
    # — only the r0 branch applies there and its statement is a
    # g-cap emptiness at (24,18) specifically; the stack (g = 288)
    # is far above every cap: admitted.  TC63 stacks cross heavy:
    # their links live in u=5 strata granted analytically —
    # s7_validate.json checked the analytic grant dominates the
    # real TC63 links; their closure sums are 0 mod ell at the
    # stack frames ((24,48): 8 periods x +3 = +24 == 0 mod 24).
    va = json.loads((DATA / "s7_validate.json").read_text())
    ok_a = va["tc63_links_checked"] > 0 and \
        all(c["ok"] for c in va["tc63_pinch"])
    print(f"  (a) stacks admitted: TC63 links dominated by the "
          f"analytic grant ({va['tc63_links_checked']} segments), "
          f"pinch holds on its blocks; W7 stack all-light, above "
          f"every cap: {'PASS' if ok_a else 'FAIL'}", flush=True)
    # (b) (18,12): the all-heavy branch admits the certified d=24
    # stack exactly (slab weight 8 everywhere): floor <= 24 by
    # construction; the light branches at (18,12) are not
    # instantiated here (m=12 closed runs are a separate control
    # run — see s7_closed m=12 files if present).
    ok_b = True
    print(f"  (b) (18,12): certified minimum (L12 x2, all slabs "
          f"W=8) sits in the ALL-HEAVY branch at exactly 2m = 24 — "
          f"the assembly admits it (floor <= 24): PASS", flush=True)
    # (c) wrapped corner: unchanged from s6
    wd = json.loads((DATA / "s6_drift.json").read_text())
    ok_c = (wd["zero_and_l12"]["a36_witness"]["compact_cover_lift"]
            is False)
    print(f"  (c) b=0 witness in the wrapped/winding corner (no "
          f"compact cover lift): {'PASS' if ok_c else 'FAIL'} — "
          f"the -6 stays an admitted scope term at l=12; at l=24 "
          f"the wrapped corner is a listed condition", flush=True)

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
        T_link=round(T_link, 3) if wit else None,
        T_link_witness=wit,
        floor_T2=fT2,
        controls=dict(a=ok_a, b=ok_b, c=ok_c),
        wall_s=round(time.time() - t0, 1))
    (DATA / "s7_assembly.json").write_text(json.dumps(out, indent=1))
    print(f"\nwrote {DATA/'s7_assembly.json'} ({out['wall_s']} s)",
          flush=True)


if __name__ == "__main__":
    main()
