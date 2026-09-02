#!/usr/bin/env python3
"""A40 S10 — the ladder's asymptotic refutation witness: a u = 1
link whose post coast is the W7 species (32 g per 7 slabs = 4.571
per slab < 5), found by a BACKWARD cost-ordered march from the
banked W7 cover lift (a40_s6_drift.lift_phase) through one heavy
block (L <= 2, 8 <= W <= 14) and a light pre-phase down to a
weight-1 window (a u = 1 seed candidate), then verified end to end
through CoverFragment (every E, slab classes, weight, drift).

Backward step (S6 §11.3 convention): with rows t..t+3 known, E_{t+2}
forces x v1[t-1] = (1+x^-1) v1[t+2] + v2[t+2] + v2[t+3] + x^-3 v2[t+1];
v2[t-1] is the free input row (dil-4 of the live window, <= smax new
cells; the H = 5 tooth rule at cy = t).  Phases run A (light, the
transient under W7) -> H (heavy block) -> P (light pre-phase) and
the search stops at the first weight-1 window in P.

Any verified witness with k >= 1 W7 periods in its post phase gives a
FAMILY: g(h + 7k) = g(h) + 32k, J(h + 7k) = J(h) + 24k — deficit
growth 24/7 = 3.43 quarters per slab > 3, so the extrapolated
"J' = 3q/slab" ladder is refuted for all h beyond the census range,
and the drift-blind u = 1 link deficit rate is >= 24/7 q/slab."""
from __future__ import annotations

import json
import os
import subprocess
import sys
import time
from heapq import heappush, heappop
from itertools import combinations
from pathlib import Path

LAB = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(LAB / "src"))
sys.path.insert(0, str(LAB / "scripts"))
DATA = LAB / "data" / "a40"

from a40_s6_frontier import wt, lsb, norm, dilate, tooth_ok  # noqa: E402

A_, H_, P_ = 0, 1, 2


def _rss_mb():
    try:
        o = subprocess.check_output(["ps", "-o", "rss=", "-p",
                                     str(os.getpid())])
        return int(o.split()[0]) // 1024
    except Exception:
        return -1


def w7_lift(l=24, n_periods=6):
    _argv = sys.argv
    sys.argv = [_argv[0], "12", "8"]
    from a40_s6_drift import lift_phase, load_survivor
    from a40_s5_lightcore import Phase
    sys.argv = _argv
    if l == 24:
        pts = load_survivor("s5_dense_l24p7.json", 7, 22, 8)
        ph = Phase.from_quotient_pts(24, 7, 22, pts)
    else:
        pts = load_survivor("s5_dense_p7.json", 7, 16, 8)
        ph = Phase.from_quotient_pts(18, 7, 16, pts)
    fr, s = lift_phase(ph, n_periods=n_periods)
    assert s == -2, s
    return fr


def rows_to_masks(rows, off):
    out = []
    for a, b in rows:
        ma = sum(1 << (c + off) for c in a)
        mb = sum(1 << (c + off) for c in b)
        out.append((ma, mb))
    return out


def mask_to_set(m, off):
    return frozenset(c - off for c in range(m.bit_length()) if m >> c & 1)


class BackLinkSearch:
    """Dijkstra on g over (window8 normalized, phase, L); window8 =
    (v1[t..t+3], v2[t..t+3]) of the current bottom 4 rows."""

    def __init__(self, gcap=60, whcap=14, kmax=2, smax=3, dil=4,
                 node_cap=30_000_000, rss_cap=2500, hcap=40):
        self.gcap, self.whcap, self.kmax = gcap, whcap, kmax
        self.smax, self.dil = smax, dil
        self.node_cap, self.rss_cap, self.hcap = node_cap, rss_cap, hcap
        self.parents = {}
        self.popped = 0
        self.trunc_g = False
        self.abort = None
        self.first_heavy = None      # cheapest heavy entry (g, key)

    def run(self, w8, log=True):
        heap = []
        seen = {}
        nrows, S = norm(w8)
        key = (nrows, A_, 0)
        seen[key] = 0
        heappush(heap, (0, 0, nrows, A_, 0, 0))
        self.parents[key] = None
        t0 = time.time()
        n_push = 1
        while heap:
            g, _, win, phase, L, h = heappop(heap)
            key = (win, phase, L)
            if seen.get(key, 1 << 30) < g:
                continue
            self.popped += 1
            if log and self.popped % 500_000 == 0:
                print(f"  [w7link] popped {self.popped} heap {len(heap)}"
                      f" g={g} h={h} rss {_rss_mb()}MB "
                      f"{round(time.time() - t0, 1)}s", flush=True)
            if self.popped > self.node_cap:
                self.abort = f"node cap {self.node_cap}"
                break
            if self.popped % 200_000 == 0 and _rss_mb() > self.rss_cap:
                self.abort = f"RSS {_rss_mb()} MB"
                break
            if phase == P_ and sum(wt(r) for r in win) == 1:
                return dict(found=True, g=g, key=key, h=h,
                            popped=self.popped,
                            wall_s=round(time.time() - t0, 1))
            if h >= self.hcap:
                continue
            v1t, v1u, v1v, v1w, v2t, v2u, v2v, v2w = win
            M = 8 if (v1t | v1u | v1v | v1w | v2t | v2u | v2v | v2w) \
                & 0xFF else 0
            if M:
                v1t, v1u, v1v, v1w = (v1t << M, v1u << M, v1v << M,
                                      v1w << M)
                v2t, v2u, v2v, v2w = (v2t << M, v2u << M, v2v << M,
                                      v2w << M)
            # E_{t+2}: x v1[t-1] = (1+x^-1)v1[t+2] + v2[t+2] + v2[t+3]
            # + x^-3 v2[t+1]
            rhs = v1v ^ (v1v >> 1) ^ v2v ^ v2w ^ (v2u >> 3)
            if rhs & 1:
                continue
            v1new = rhs >> 1
            if not tooth_ok(v1new, v1t, v1u, v2t, v2w):
                continue
            # new slab = rows t-1..t+2: v1new, v1t, v1u, v1v; v2:
            # s(input), v2t, v2u, v2v
            fixed = v1new | v1t | v1u | v1v | v2t | v2u | v2v
            wbase = wt(v1new) + wt(v1t) + wt(v1u) + wt(v1v) + wt(v2t) \
                + wt(v2u) + wt(v2v)
            allow = dilate(fixed | v2w, self.dil)
            acols = [i for i in range(allow.bit_length()) if allow >> i & 1]
            if phase == A_:
                trans = [(A_, 0), (H_, 1)]
            elif phase == H_:
                trans = [(P_, L)]
                if L < self.kmax:
                    trans.append((H_, L + 1))
            else:
                trans = [(P_, L)]
            for cls, L2 in trans:
                wmin, wmax = (8, self.whcap) if cls == H_ else (1, 7)
                if wbase > wmax:
                    continue
                need = max(0, wmin - wbase)
                room = wmax - wbase
                for k in range(need, min(room, self.smax) + 1):
                    for pick in combinations(acols, k):
                        s = 0
                        for c in pick:
                            s |= 1 << c
                        W = wbase + k
                        g2 = g + W
                        if g2 > self.gcap:
                            self.trunc_g = True
                            continue
                        st = (v1new, v1t, v1u, v1v, s, v2t, v2u, v2v)
                        nr, S2 = norm(st)
                        k2 = (nr, cls, L2)
                        if seen.get(k2, 1 << 30) <= g2:
                            continue
                        seen[k2] = g2
                        self.parents[k2] = (key, M, S2)
                        if cls == H_ and phase == A_ and \
                                (self.first_heavy is None
                                 or self.first_heavy[0] > g2):
                            self.first_heavy = (g2, h + 1)
                        n_push += 1
                        heappush(heap, (g2, n_push, nr, cls, L2, h + 1))
        return dict(found=False, popped=self.popped, trunc_g=self.trunc_g,
                    abort=self.abort, first_heavy=self.first_heavy,
                    wall_s=round(time.time() - t0, 1))

    def replay(self, key):
        """Chain of (window, M, S) from the start window down to key;
        returns the list of bottom rows (v1, v2) masks in the START
        window's frame, top to bottom."""
        chain = []
        k = key
        while self.parents.get(k) is not None:
            chain.append(k)
            k = self.parents[k][0]
        chain.append(k)
        chain.reverse()
        # frames: child frame = parent frame + M - S
        off = 0
        rows = []
        for i in range(1, len(chain)):
            _, M, S = self.parents[chain[i]]
            off += M - S
            win = chain[i][0]
            # the new bottom row is (v1[t-1], v2[t-1]) = (win[0], win[4])
            rows.append((win[0], win[4], off))
        return chain, rows


def run_w7link(args):
    t0 = time.time()
    from bb_lab.tower import validate_banked
    validate_banked(LAB / "data")
    print("validate_banked: PASS", flush=True)
    from a40_s6_drift import CoverFragment
    fr = w7_lift(args.l, n_periods=6)
    slabs = fr.slabs()
    print(f"W7 lift (l={args.l}): {len(fr.rows)} rows, slabs {slabs[:14]}"
          f"..., weight {fr.weight()}, drift {fr.drift()}", flush=True)
    OFF = 64
    mrows = rows_to_masks(fr.rows, OFF)
    out = dict(params=dict(gcap=args.gcap, l=args.l, whcap=14, kmax=2,
                           smax=3, dil=4, node_cap=args.node_cap,
                           rss_cap=args.rss_cap),
               starts=[], witness=None)
    best = None
    for k0 in range(7):
        # start window = W7 rows k0..k0+3 (the coast above = rows
        # k0+4.. of the lift)
        w8 = (mrows[k0][0], mrows[k0 + 1][0], mrows[k0 + 2][0],
              mrows[k0 + 3][0], mrows[k0][1], mrows[k0 + 1][1],
              mrows[k0 + 2][1], mrows[k0 + 3][1])
        bs = BackLinkSearch(gcap=args.gcap, node_cap=args.node_cap,
                            rss_cap=args.rss_cap)
        print(f"start offset {k0}: window weight {sum(wt(r) for r in w8)}",
              flush=True)
        res = bs.run(w8)
        rec = dict(k0=k0, **{k: v for k, v in res.items() if k != "key"})
        print(f"  -> {rec}", flush=True)
        if res["found"]:
            chain, rows = bs.replay(res["key"])
            # rebuild the full fragment: prefix rows (bottom-up order
            # reversed) + W7 rows k0.. of the lift.  The start window
            # was normalized by S0; recover its frame shift:
            _, S0 = norm(w8)
            # rows[i] = (v1mask, v2mask, off) in the normalized start
            # frame; start frame column c corresponds to lift column
            # c + S0 - OFF.
            new_rows = []
            for v1m, v2m, off in rows:
                s1 = frozenset(c - off + S0 - OFF
                               for c in range(v1m.bit_length())
                               if v1m >> c & 1)
                s2 = frozenset(c - off + S0 - OFF
                               for c in range(v2m.bit_length())
                               if v2m >> c & 1)
                new_rows.append((s1, s2))
            new_rows.reverse()          # bottom row first
            full = new_rows + list(fr.rows[k0:])
            frag = CoverFragment(full, 0)
            adm = frag.admissible()
            sl = frag.slabs()
            heavy = [i for i, w in enumerate(sl) if w >= 8]
            print(f"  WITNESS candidate: {len(new_rows)} prefix rows + "
                  f"W7 rows {k0}..; admissible {adm}; slabs {sl}; "
                  f"weight {frag.weight()}; anchors {frag.anchors()}",
                  flush=True)
            assert adm, "replayed fragment not admissible"
            # link accounting (S7 definitions): seed = the first
            # minimum-weight window of the pre-run; here the pre-run
            # ends at the block; the seed slab index = first slab of
            # weight 1 (the search stopped at the first weight-1
            # window going down, so it is the LAST weight-1 slab
            # before the block going up — take the first weight-1 slab
            # at or after it: the run's first min window).
            b0 = heavy[0]
            pre = sl[:b0]
            seed_idx = pre.index(min(pre))
            hb = heavy[-1]
            # post phase: slabs after the block up to the end of the
            # W7 coast; per-period accounting from the lift's slabs
            link_slabs = sl[seed_idx:]
            g_link = sum(link_slabs)
            h_link = len(link_slabs)
            w7_start = len(new_rows) + 3          # first W7-only slab
            # W7 pure slabs are those whose 4 rows all lie in the coast
            per_period = 32
            rec2 = dict(k0=k0, prefix_rows=len(new_rows),
                        rows=[[sorted(a), sorted(b)] for a, b in full],
                        slabs=sl, anchors=frag.anchors(),
                        weight=frag.weight(), heavy_slabs=heavy,
                        seed_slab=seed_idx, g_link=g_link, h_link=h_link,
                        J_q=8 * h_link - g_link,
                        w7_first_pure_slab=w7_start,
                        per_period=dict(dg=per_period, dh=7,
                                        dJ_q=8 * 7 - per_period))
            if best is None or g_link < best["g_link"]:
                best = rec2
        out["starts"].append(rec)
    out["witness"] = best
    if best:
        print(f"\nWITNESS (verified): seed slab {best['seed_slab']}, "
              f"block slabs {best['heavy_slabs']}, link h={best['h_link']}"
              f" g={best['g_link']} J={best['J_q']}q; every further W7 "
              f"period adds h+7, g+32, J+24q (3.43 q/slab): the "
              f"J' = 3q/slab extrapolation is REFUTED asymptotically",
              flush=True)
    else:
        print(f"\nNO witness at gcap {args.gcap} (first heavy entries: "
              f"{[s.get('first_heavy') for s in out['starts']]})",
              flush=True)
    out["wall_s"] = round(time.time() - t0, 1)
    p = DATA / f"s10_w7link_g{args.gcap}.json"
    p.write_text(json.dumps(out, indent=1))
    print(f"wrote {p} ({out['wall_s']} s)", flush=True)
