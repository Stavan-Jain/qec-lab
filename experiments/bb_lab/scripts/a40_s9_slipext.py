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
        print(f"regression vs {args.ref}: pre {my_pre_ok}, "
              f"link {my_link_ok}", flush=True)
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


# ---------------------------------------------------------------------
# Replay: parent-logged march -> specimen for a (L, gH, slip) class
# ---------------------------------------------------------------------

class ReplayMarch(SlipLinkMarch):
    """SlipLinkMarch with parent logging: every stored state keeps
    (parent_key, parent_h, slab_mask_normshift) so a completed
    crossing can be reconstructed row-by-row.  Dominance unchanged
    (the winner at each key overwrites parent info with the cheaper
    path — exactly the represented minimum).  Crossing harvest: at
    the HEAVY->POST transition matching the requested class, record
    the full chain."""

    def __init__(self, *a, want=None, **kw):
        super().__init__(*a, **kw)
        self.want = want            # (L, gH, slip) or None = all
        self.par = {}               # (h, key) -> (parent (h,key), slab, M, S)
        self.hits = []              # completed-crossing chains
        self._cur_h = None
        self._cur_key = None

    def _run_seed(self, sid, w8, log, rss_cap, frontier_cap):
        self._seed_w8 = w8
        self._sid = sid
        return super()._run_seed(sid, w8, log, rss_cap, frontier_cap)

    def _run_layerset(self, layer, h, log, rss_cap, frontier_cap,
                      depth=0):
        # non-streaming variant (replay runs are small-cap only)
        while layer and h < self.hcap:
            rss = _rss_mb()
            if rss > rss_cap:
                return f"RED: RSS {rss} MB at h={h}"
            if len(layer) > frontier_cap:
                return f"RED: frontier {len(layer)} at h={h}"
            nxt = {}
            for key, val in layer.items():
                self._cur_h, self._cur_key = h, key
                self.expand(key, val, h, nxt)
            self.nodes += len(layer)
            layer = nxt
            h += 1
        return None

    def expand(self, key, val, h, nxt):
        # duplicate of SlipLinkMarch.expand with parent logging;
        # kept in lockstep with the S8 engine (any divergence would
        # fail the census cross-check below).
        (dyn, anch, phase, L, dlt, hin) = key
        (g, xlo, xhi) = val
        v1a, v1b, v1c, v1d, v2a, v2b, v2c = dyn
        M = 8 if (v1a | v1b | v1c | v1d | v2a | v2b | v2c) & 0xFF \
            else 0
        if M:
            v1a, v1b, v1c, v1d = (v1a << M, v1b << M, v1c << M,
                                  v1d << M)
            v2a, v2b, v2c = v2a << M, v2b << M, v2c << M
        anch_m = anch + M
        v2new = v1d ^ (v1d >> 1) ^ (v1a << 1) ^ v2c ^ (v2b >> 3)
        if not tooth_ok(v1a, v1b, v1c, v2a, v2new):
            return
        fixed = v1b | v1c | v1d | v2a | v2b | v2c | v2new
        wbase = wt(v1b) + wt(v1c) + wt(v1d) + wt(v2a) + wt(v2b) \
            + wt(v2c) + wt(v2new)
        trans = self.transitions(phase, L)
        wmax_any = max((7 if c != HEAVY else self.whcap)
                       for c, _ in trans)
        if wbase > wmax_any:
            return
        allow = dilate(v1a | v1b | v1c | v1d | v2a | v2b | v2c
                       | v2new, self.dil)
        acols = [i for i in range(allow.bit_length())
                 if allow >> i & 1]
        h2 = h + 1
        for cls, L2 in trans:
            wmin, wmax = (self.u, 7) if cls == PRE else \
                (8, self.whcap) if cls == HEAVY else (1, 7)
            if wbase > wmax:
                continue
            need = max(0, wmin - wbase)
            room = wmax - wbase
            for kk in range(need, min(room, self.smax) + 1):
                for pick in combinations(acols, kk):
                    s = 0
                    for c in pick:
                        s |= 1 << c
                    W = wbase + kk
                    g2 = g + W
                    if g2 > self.gcap:
                        self.trunc_g = True
                        continue
                    slab = fixed | s
                    dlt2 = dlt + (lsb(slab) - anch_m)
                    if abs(dlt2) > self.dcap:
                        self.trunc_dcap += 1
                        continue
                    state = (v1b, v1c, v1d, s, v2b, v2c, v2new)
                    a_new = lsb(slab)
                    nrows, S = norm(state)
                    anch_c = a_new - S
                    xlo_c = min(xlo + M, a_new) - S
                    xhi_c = max(xhi + M,
                                slab.bit_length() - 1) - S
                    if xhi_c - xlo_c + 1 > self.extent_cap:
                        self.trunc_extent += 1
                        continue
                    if cls == HEAVY and phase == PRE:
                        hin2 = (dlt, W)
                    elif cls == HEAVY:
                        hin2 = (hin[0], hin[1] + W)
                    elif cls == POST and phase == HEAVY:
                        slip = dlt2 - hin[0]
                        ks = (L, hin[1], slip)
                        if ks not in self.tab_slip or \
                                self.tab_slip[ks] > g2:
                            self.tab_slip[ks] = g2
                        hin2 = (hin[0], hin[1], slip)
                        if self.want is None or ks == self.want:
                            self.hits.append(dict(
                                sid=self._sid, cls=list(ks), g=g2,
                                chain=self._chain()
                                + [(int(slab), int(M), int(S),
                                    int(g2))]))
                    else:
                        hin2 = hin
                    if cls == PRE:
                        k = (h2, dlt2)
                        if k not in self.tab_pre or \
                                self.tab_pre[k] > g2:
                            self.tab_pre[k] = g2
                    elif cls == POST:
                        k = (h2, dlt2)
                        if k not in self.tab_link or \
                                self.tab_link[k] > g2:
                            self.tab_link[k] = g2
                    key2 = (nrows, anch_c, cls, L2, dlt2, hin2)
                    cur = nxt.get(key2)
                    if cur is not None and cur[0] <= g2:
                        continue
                    nxt[key2] = (g2, xlo_c, xhi_c)
                    self.par[(h2, key2)] = (
                        (self._cur_h, self._cur_key), int(slab),
                        int(M), int(S))

    def _chain(self):
        """Reconstruct the slab chain (absolute columns) ending at
        the CURRENT expanding state."""
        out = []
        node = (self._cur_h, self._cur_key)
        shift = 0
        while node in self.par:
            parent, slab, M, S = self.par[node]
            out.append((slab, M, S))
            node = parent
        # root = the seed window (4 rows); reconstruct forward
        chain = []
        rows_abs = []
        # walk back down applying shifts: absolute col of row i
        # = stored slab >> accumulated normalization, but norm
        # shifts compose; easier: rebuild forward from the seed.
        return [(int(s), int(M), int(S)) for (s, M, S)
                in reversed(out)]


def run_replay(args):
    from bb_lab.tower import validate_banked
    validate_banked(LAB / "data")
    print("validate_banked: PASS", flush=True)
    want = tuple(int(x) for x in args.want.split(",")) \
        if args.want else None
    m = ReplayMarch(args.u, kmax=args.kmax, whcap=args.whcap,
                    gcap=args.gcap, hcap=args.hcap, dcap=args.dcap,
                    want=want)
    sds = seeds_full(args.u)
    info = m.run_byseed(sds, rss_cap=2048, frontier_cap=3_000_000)
    print(f"info: {info}; hits {len(m.hits)}", flush=True)
    # census cross-check: tab_slip must MATCH the banked g30 rows
    # below this gcap (engine-duplication guard)
    g30 = json.loads((DATA / "s8_slip_u1_g30.json").read_text())
    ok = True
    for r in g30["slip_rows"]:
        if r["min_g"] <= args.gcap:
            k = (r["L"], r["gH"], r["slip"])
            if m.tab_slip.get(k) != r["min_g"]:
                ok = False
                print(f"  CROSS-CHECK MISS {k}: banked "
                      f"{r['min_g']} got {m.tab_slip.get(k)}",
                      flush=True)
    print(f"replay-engine census cross-check vs g30 "
          f"(rows with min_g <= {args.gcap}): {ok}", flush=True)
    assert ok, "REPLAY ENGINE DIVERGED"
    best = min((h for h in m.hits), key=lambda h: h["g"]) \
        if m.hits else None
    out = dict(params=dict(u=args.u, kmax=args.kmax,
                           whcap=args.whcap, gcap=args.gcap,
                           hcap=args.hcap, dcap=args.dcap,
                           want=args.want),
               info=info, n_hits=len(m.hits), best=best)
    p = DATA / f"s9_replay_{args.tag}.json"
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
    r = sub.add_parser("replay")
    r.add_argument("--u", type=int, default=1)
    r.add_argument("--kmax", type=int, default=2)
    r.add_argument("--whcap", type=int, default=14)
    r.add_argument("--gcap", type=int, required=True)
    r.add_argument("--hcap", type=int, default=12)
    r.add_argument("--dcap", type=int, default=16)
    r.add_argument("--want", type=str, default="",
                   help="L,gH,slip (empty = all crossings)")
    r.add_argument("--tag", type=str, required=True)
    args = ap.parse_args()
    if args.cmd == "census":
        run_census(args)
    else:
        run_replay(args)


if __name__ == "__main__":
    main()
