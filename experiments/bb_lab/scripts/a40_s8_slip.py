#!/usr/bin/env python3
"""A40 S8 — Stage 2: the heavy-slip weight bound — machine-derived
per-step anchor caps, an EXHAUSTIVE slip census of enumerated links,
and the verdict on the strong local lemma.

PER-STEP ANCHOR CAPS (the provable piece; constants derived from
the Laurent supports, both lanes, verified here on random states):

  Forced-row reach.  The forward-forced row is an F2 combination of
  window rows with column shifts; its support therefore lies within
  [min(window) + smin, max(window) + smax_r] where smin/smax_r are
  the extreme shifts of the solved recurrence:
    y-lane  v2[t+1] = (1+x^-1) v1[t] + x v1[t-3] + v2[t]
            + x^-3 v2[t-1]:            shifts {-3..+1}
    mirror  u1[t+1] = x^-3 ((1+x^-1) u2[t] + x u2[t-3] + u1[t]
            + u1[t-1]):                shifts {-4..-2}
  Input-row reach (SCOPE, dil-D growth): within [min - D, max + D]
  of the live window.

  LEMMA S8.2a (per-step, scope dil-D, both lanes).  With D = 4:
    A_{j+1} >= A_j - 4   (min-anchor gauge: left slip <= 4/step)
    Amax_{j+1} <= Amax_j + 4  (max-anchor gauge: right <= 4/step).
  Proof: every new cell (forced or input) sits >= min(window) - 4
  and <= max(window) + 4; the entering slab's min/max are over old
  rows (within the window envelope) and the new rows.  QED.
  NOTE the caps are ONE-SIDED PER GAUGE: the min-anchor can jump
  RIGHT by up to the slab span (a strand death), and the max-anchor
  can jump LEFT likewise — those are not growth events and no
  footprint constant bounds them by local weight.

EXHAUSTIVE SLIP CENSUS (links, y-lane).  The S7 link march measured
slips only on table-improving representatives.  Here the heavy-entry
drift d_in is part of the dominance key from the HEAVY phase onward,
so every (block shape, slip) pair reachable below the caps appears:
tab_slip[(L, gH, slip)] = min g over links realizing it.  The PRE
phase is unchanged (regression-checked against the banked s7 pre
table).  slip := dlt(entry into POST) - dlt(at heavy entry), i.e.
the anchor movement across the heavy block including its exit row.

VERDICT TARGET: is |slip| <= c * gH + c0 with small constants on
ALL enumerated links, and does the bound survive as a lemma (the
two-strand adversary: a cheap left excursion dying INTO the block
makes slip ~ excursion length at fixed block weight — the census
range and the per-(g, slip) frontier measure how much slip the caps
actually admit)."""
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

PRE, HEAVY, POST = 0, 1, 2


def _rss_mb():
    """CURRENT process RSS in MB (ps-based; ru_maxrss is a
    lifetime peak — see a40_s8_xlane and §13.6)."""
    import os
    import subprocess
    try:
        out = subprocess.check_output(
            ["ps", "-o", "rss=", "-p", str(os.getpid())])
        return int(out.split()[0]) // 1024
    except Exception:
        import resource
        return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss \
            // (1024 * 1024)


# ---------------------------------------------------------------------
# Stage 2a: the per-step caps — derive + verify
# ---------------------------------------------------------------------

def derive_caps():
    """Machine-derive the forced-row shift extremes for both lanes
    from the recurrence supports, and verify per-step anchor caps on
    random march states of both engines."""
    import random
    out = {}
    # y-lane solved form: v2[t+1] = sum shifts of window rows
    y_shifts = [0, -1, +1, 0, -3]
    m_shifts = [a - 3 for a in (0, -1, +1, 0, 0)]   # mirror: then -3
    out["y_forced_reach"] = [min(y_shifts), max(y_shifts)]
    out["m_forced_reach"] = [min(m_shifts), max(m_shifts)]
    assert out["y_forced_reach"] == [-3, 1]
    assert out["m_forced_reach"] == [-4, -2]
    # verify on random states: y-lane forced row (s6 kernel)
    rng = random.Random(2)
    for _ in range(3000):
        rows = [sum(1 << c for c in rng.sample(range(0, 12),
                                               rng.randint(1, 4)))
                for _ in range(7)]
        v1a, v1b, v1c, v1d, v2a, v2b, v2c = [r << 8 for r in rows]
        v2new = v1d ^ (v1d >> 1) ^ (v1a << 1) ^ v2c ^ (v2b >> 3)
        allm = v1a | v1b | v1c | v1d | v2a | v2b | v2c
        if v2new:
            assert lsb(v2new) >= lsb(allm) - 3
            assert v2new.bit_length() <= allm.bit_length() + 1
    from a40_s8_xlane import m_forced
    for _ in range(3000):
        rows = [sum(1 << c for c in rng.sample(range(0, 12),
                                               rng.randint(1, 4)))
                for _ in range(8)]
        p = [r << 8 for r in rows]
        u1new = m_forced(p[2], p[3], p[4], p[7])
        allm = 0
        for r in p:
            allm |= r
        if u1new:
            assert lsb(u1new) >= lsb(allm) - 4
            assert u1new.bit_length() <= allm.bit_length() - 2 + 1
    out["per_step_caps"] = dict(
        left_min_gauge=4, right_max_gauge=4, scope="dil-4 growth",
        note="one-sided per gauge; strand deaths unbounded by "
             "local weight")
    print("caps: y forced reach [-3,+1], mirror [-4,-2]; with "
          "dil-4 inputs => per-step left slip <= 4 (min gauge), "
          "right <= 4 (max gauge): verified on 6000 random states",
          flush=True)
    return out


# ---------------------------------------------------------------------
# Stage 2b: the slip-resolved link march (y-lane)
# ---------------------------------------------------------------------

class SlipLinkMarch:
    """The S7 LinkMarch (layered, byseed) with the heavy-entry drift
    d_in carried in the dominance key from HEAVY onward, and the
    completed-crossing slip carried through POST.  Tables:
      tab_pre[(h, dlt)]           (regression target vs banked s7)
      tab_slip[(L, gH, slip)] -> min g   (EXHAUSTIVE below caps)
      tab_link[(h, dlt)]          (regression: min over slips)
    slip = dlt at POST entry - dlt at HEAVY entry."""

    def __init__(self, u, kmax=2, whcap=14, smax=3, dil=4, hcap=19,
                 gcap=24, extent_cap=34, dcap=30):
        self.u, self.kmax, self.whcap = u, kmax, whcap
        self.smax, self.dil = smax, dil
        self.hcap, self.gcap = hcap, gcap
        self.extent_cap, self.dcap = extent_cap, dcap
        self.tab_pre = {}
        self.tab_link = {}
        self.tab_link_L = {}   # (h, L, dlt) -> min g (S9 additive:
        #                        assembly-consumable resolution; no
        #                        dominance-key change, so the S8
        #                        regressions are unaffected)
        self.tab_slip = {}
        self.nodes = 0
        self.trunc_g = False
        self.trunc_extent = 0
        self.trunc_dcap = 0
        self.aborts = []

    def run_byseed(self, seeds, log=True, rss_cap=2048,
                   frontier_cap=3_000_000):
        t0 = time.time()
        nseed = 0
        for sid, w8 in enumerate(seeds):
            if sum(wt(r) for r in w8) != self.u:
                continue
            nseed += 1
            ab = self._run_seed(sid, w8, log, rss_cap, frontier_cap)
            if ab:
                self.aborts.append(dict(sid=sid, why=ab))
                print(f"  seed {sid}: ABORT {ab}", flush=True)
        return dict(n_seeds=nseed, nodes=self.nodes,
                    trunc_g=self.trunc_g,
                    trunc_extent=self.trunc_extent,
                    trunc_dcap=self.trunc_dcap,
                    aborts=self.aborts,
                    wall_s=round(time.time() - t0, 1))

    def _run_seed(self, sid, w8, log, rss_cap, frontier_cap):
        u = self.u
        dyn = (w8[0], w8[1], w8[2], w8[3], w8[5], w8[6], w8[7])
        allm = 0
        for r in w8:
            allm |= r
        anchor = lsb(allm)
        # key: (dyn, anch, phase, L, dlt, hin) with hin = None in
        # PRE, (d_in, gH) in HEAVY, (d_in, gH, slip) in POST
        layer = {(dyn, anchor, PRE, 0, 0, None):
                 (u, anchor, allm.bit_length() - 1)}
        k = (1, 0)
        if k not in self.tab_pre or self.tab_pre[k] > u:
            self.tab_pre[k] = u
        return self._run_layerset(layer, 1, log, rss_cap,
                                  frontier_cap)

    # recursive layer chunking (S8): a layer above split_thresh is
    # partitioned and each chunk marched independently to
    # completion.  Sound: cross-chunk dominance is only an
    # optimization, so splitting can only OVER-explore; the tables
    # (updated during expansion) merge automatically.  Guards run
    # per chunk, with a mid-expansion RSS check.
    SPLIT_THRESH = 180_000

    def _run_layerset(self, layer, h, log, rss_cap, frontier_cap,
                      depth=0):
        import os
        import pickle
        import tempfile
        t0 = time.time()
        while layer and h < self.hcap:
            rss = _rss_mb()
            if rss > rss_cap:
                return f"RED: RSS {rss} MB at h={h}"
            nxt = {}
            parts = []

            def spill(d):
                fd, pth = tempfile.mkstemp(
                    suffix=f"_slip_d{depth}_h{h}.pkl")
                with os.fdopen(fd, "wb") as f:
                    pickle.dump(list(d.items()), f, protocol=4)
                parts.append(pth)

            nexp = 0
            for key, val in layer.items():
                self.expand(key, val, h, nxt)
                nexp += 1
                if len(nxt) > self.SPLIT_THRESH:
                    spill(nxt)
                    nxt = {}
                if nexp % 65536 == 0 and _rss_mb() > rss_cap:
                    for pth in parts:
                        os.remove(pth)
                    return f"RED: RSS {_rss_mb()} MB mid-layer " \
                        f"h={h}"
            self.nodes += len(layer)
            if parts:
                if nxt:
                    spill(nxt)
                nxt = None
                layer = None
                if log:
                    print(f"[slip d{depth}] h={h + 1} streamed "
                          f"into {len(parts)} parts rss "
                          f"{_rss_mb()}MB "
                          f"{round(time.time() - t0, 1)}s",
                          flush=True)
                try:
                    for pi, pth in enumerate(parts):
                        with open(pth, "rb") as f:
                            chunk = dict(pickle.load(f))
                        os.remove(pth)
                        ab = self._run_layerset(
                            chunk, h + 1, log, rss_cap,
                            frontier_cap, depth + 1)
                        chunk = None
                        if ab:
                            return ab
                finally:
                    for pth in parts:
                        if os.path.exists(pth):
                            os.remove(pth)
                return None
            if log and (h % 4 == 0 or len(nxt) > 100000):
                print(f"[slip d{depth}] h={h + 1} states "
                      f"{len(nxt)} rss {_rss_mb()}MB "
                      f"{round(time.time() - t0, 1)}s", flush=True)
            layer = nxt
            h += 1
        return None

    def transitions(self, phase, L):
        if phase == PRE:
            out = [(PRE, 0)]
            if self.kmax >= 1:
                out.append((HEAVY, 1))
            return out
        if phase == HEAVY:
            out = [(POST, L)]
            if L < self.kmax:
                out.append((HEAVY, L + 1))
            return out
        return [(POST, L)]

    def expand(self, key, val, h, nxt):
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
                    # hin bookkeeping
                    if cls == HEAVY and phase == PRE:
                        hin2 = (dlt, W)          # (d_in, gH so far)
                    elif cls == HEAVY:
                        hin2 = (hin[0], hin[1] + W)
                    elif cls == POST and phase == HEAVY:
                        slip = dlt2 - hin[0]
                        ks = (L, hin[1], slip)
                        if ks not in self.tab_slip or \
                                self.tab_slip[ks] > g2:
                            self.tab_slip[ks] = g2
                        hin2 = (hin[0], hin[1], slip)
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
                        kL = (h2, L2, dlt2)
                        if kL not in self.tab_link_L or \
                                self.tab_link_L[kL] > g2:
                            self.tab_link_L[kL] = g2
                    key2 = (nrows, anch_c, cls, L2, dlt2, hin2)
                    cur = nxt.get(key2)
                    if cur is not None and cur[0] <= g2:
                        continue
                    nxt[key2] = (g2, xlo_c, xhi_c)


def run_census(args):
    import gc
    gc.disable()
    t0 = time.time()
    from bb_lab.tower import validate_banked
    validate_banked(LAB / "data")
    print("validate_banked: PASS", flush=True)
    out = {"caps": derive_caps()}

    print(f"slip census: u={args.u} kmax={args.kmax} "
          f"whcap={args.whcap} gcap={args.gcap} (byseed)",
          flush=True)
    m = SlipLinkMarch(args.u, kmax=args.kmax, whcap=args.whcap,
                      gcap=args.gcap, hcap=args.hcap,
                      dcap=args.dcap)
    sds = seeds_full(args.u)
    info = m.run_byseed(sds)
    print(f"info: {info}", flush=True)
    # regression vs the banked s7 byseed link tables (same caps)
    reg = {}
    ref_name = args.ref
    if ref_name:
        ref = json.loads((DATA / ref_name).read_text())
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
        # the slip-resolved keys can only ADD states, never change
        # the (h, dlt) minima; demand exact equality on in-cap rows
        reg = dict(ref=ref_name, pre_equal=bool(my_pre_ok),
                   link_equal=bool(my_link_ok),
                   pre_rows=len(ref_pre), link_rows=len(ref_link))
        print(f"regression vs {ref_name}: pre {my_pre_ok}, "
              f"link {my_link_ok}", flush=True)
        assert my_pre_ok and my_link_ok, "REGRESSION FAILED"
    out["regression"] = reg

    # the census verdict: slip range per (L, gH), and the affine
    # test |slip| <= c*(gH - 8) + c0
    rows = []
    viol = []
    best_c, best_c0 = None, None
    for (L, gH, slip), g in sorted(m.tab_slip.items()):
        rows.append(dict(L=L, gH=gH, slip=slip, min_g=g))
    if rows:
        # fit the minimal c0 at c = 1/2 (charter's suggested shape)
        # and the minimal c at c0 = 1: report both
        import math
        c_half_c0 = max(abs(r["slip"]) - 0.5 * (r["gH"] - 8)
                        for r in rows)
        c0_one_c = max((abs(r["slip"]) - 1) / max(r["gH"] - 8, 1e-9)
                       for r in rows if r["gH"] > 8) \
            if any(r["gH"] > 8 for r in rows) else 0.0
        mx = max(rows, key=lambda r: abs(r["slip"]))
        out["census"] = dict(
            n=len(rows), max_abs_slip=abs(mx["slip"]),
            argmax=mx,
            fit_c_half=dict(c=0.5, c0=round(c_half_c0, 3)),
            fit_c0_one=dict(c=round(c0_one_c, 3), c0=1))
        print(f"census: {len(rows)} (L, gH, slip) classes; max "
              f"|slip| {abs(mx['slip'])} at {mx}; "
              f"|slip| <= (gH-8)/2 + {c_half_c0}", flush=True)
    out["slip_rows"] = rows
    out["info"] = info
    out["params"] = dict(u=args.u, kmax=args.kmax,
                         whcap=args.whcap, gcap=args.gcap,
                         hcap=args.hcap, dcap=args.dcap)
    p = DATA / (f"s8_slip_u{args.u}_g{args.gcap}.json")
    p.write_text(json.dumps(out, indent=1))
    print(f"wrote {p} ({round(time.time() - t0, 1)} s)", flush=True)


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("caps")
    sub.add_parser("species")
    s = sub.add_parser("census")
    s.add_argument("--u", type=int, default=1)
    s.add_argument("--kmax", type=int, default=2)
    s.add_argument("--whcap", type=int, default=14)
    s.add_argument("--gcap", type=int, required=True)
    s.add_argument("--hcap", type=int, default=19)
    s.add_argument("--dcap", type=int, default=30)
    s.add_argument("--ref", type=str, default="")
    args = ap.parse_args()
    if args.cmd == "caps":
        from bb_lab.tower import validate_banked
        validate_banked(LAB / "data")
        print("validate_banked: PASS", flush=True)
        derive_caps()
    elif args.cmd == "species":
        species_check()
    else:
        run_census(args)




# ---------------------------------------------------------------------
# Stage 2c: per-step cap verification on the banked species lifts
# (both lanes) — appended after the census run; invoked as `species`.
# ---------------------------------------------------------------------

def species_check():
    from bb_lab.tower import validate_banked
    validate_banked(LAB / "data")
    print("validate_banked: PASS", flush=True)
    out = {}
    _argv = sys.argv
    sys.argv = [_argv[0], "12", "8"]
    from a40_s6_drift import lift_phase, load_survivor
    from a40_s5_lightcore import Phase
    sys.argv = _argv
    rows = []
    for name, l, p, d, w in [("W7_l18", 18, 7, 16, 8),
                             ("W7_l24", 24, 7, 22, 8),
                             ("TC63_l18", 18, 6, 3, 10),
                             ("TC63_l24", 24, 6, 3, 10)]:
        fname = ("s5_dense_p7.json" if (l, p) == (18, 7) else
                 "s5_dense_p6.json" if (l, p) == (18, 6) else
                 f"s5_dense_l24p{p}.json")
        ph = Phase.from_quotient_pts(
            l, p, d, load_survivor(fname, p, d, w))
        fr, s = lift_phase(ph, n_periods=3)
        anch = fr.anchors()
        incs = [anch[i + 1] - anch[i] for i in range(len(anch) - 1)]
        assert min(incs) >= -4, (name, incs)
        rows.append(dict(name=name, min_step=min(incs),
                         max_step=max(incs)))
        print(f"  {name}: per-step anchor increments in "
              f"[{min(incs)}, {max(incs)}] — left cap -4 respected",
              flush=True)
    # mirror species TC63' via the S8 control lift
    from a40_s8_xlane import MirrorFragment, P_M, Q_M
    from a40_s5_twisted_atlas import (
        transform_class, tr_supp, atlas_run)
    import numpy as np
    w1, w2, gg = transform_class(6, 3)
    trP = tr_supp(P_M, w1, w2, gg)
    trQ = tr_supp(Q_M, w1, w2, gg)
    rws, _ = atlas_run(trP, trQ, gg, 11)
    pick = next(r for r in rws if r["nontrivial"]
                and r["weight"] == 10)
    Wm = np.array([[w1[0], w1[1]], [w2[0], w2[1]]], dtype=np.int64)
    Wi = np.array([[Wm[1, 1], -Wm[0, 1]], [-Wm[1, 0], Wm[0, 0]]],
                  dtype=np.int64)
    base_pts = [(int(Wi[0, 0] * c + Wi[0, 1] * y),
                 int(Wi[1, 0] * c + Wi[1, 1] * y), blk)
                for (c, y, blk) in pick["pts"]]
    rws2 = []
    for t in range(0, 6 * 3 + 5):
        s1, s2 = set(), set()
        for (e0, e1, blk) in base_pts:
            if (t - e1) % 6 == 0:
                k = (t - e1) // 6
                (s1 if blk == 0 else s2).add(e0 + 3 * k)
        rws2.append((frozenset(s1), frozenset(s2)))
    frT = MirrorFragment(rws2, 0)
    assert frT.admissible()
    anch = frT.anchors()
    incs = [anch[i + 1] - anch[i] for i in range(len(anch) - 1)]
    assert min(incs) >= -4, incs
    rows.append(dict(name="TC63p_mirror", min_step=min(incs),
                     max_step=max(incs)))
    print(f"  TC63' (mirror): per-step increments in "
          f"[{min(incs)}, {max(incs)}] — left cap -4 respected",
          flush=True)
    out["species_rows"] = rows
    (DATA / "s8_slip_species.json").write_text(
        json.dumps(out, indent=1))
    print(f"wrote {DATA/'s8_slip_species.json'}", flush=True)

if __name__ == "__main__":
    main()
