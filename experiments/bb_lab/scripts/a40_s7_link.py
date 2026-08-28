#!/usr/bin/env python3
"""A40 S7 — Stage 1: the boundary-coupling engine.

LINK MARCH: forward-only min-window-seeded march with a slab-class
schedule — light pre-phase (W in [u,7]), heavy block (W in
[8, whcap], at most kmax slabs), light post-phase (W in [1,7]) —
enforcing every interface E THROUGH the heavy block.  CLOSED MARCH:
per-seed forward march whose readout demands return to the seed
window (normalized 4-row-window state equality) with total drift
== 0 (mod ell) after exactly m steps — the r in {0,1} branches
(all-light r=0 at phase 0; one-block r=1 at phase 2).

LINK DECOMPOSITION (the s7 bookkeeping; note §12).  Split every
maximal light run at its first minimum-weight window (the seed).
The cyclic walk = seeds + LINKS, where link_i runs from seed_i
(inclusive) upward through run i's fwd piece, heavy block i, run
(i+1)'s lower piece, ENDING AT seed_{i+1} (inclusive, so its final
slab is the next seed's window).  Then:
 - every slab lies in exactly one link except seed slabs (exactly
   two), so sum_i g_link_i = 4|v| + sum_i u_i and
   sum_i D_link_i = D_total + sum_i (2 - u_i/4);
 - drift telescopes with NO uncounted anchor steps:
   sum_i delta_link_i == 0 (mod ell) (winding included) — the s6
   assembly's free heavy slips are now MEASURED inside links;
 - every E_j is enforced by exactly one link (a fwd march reaching
   final row s enforces E_{s-1}; the next link, seeded on rows
   [s-3, s], enforces from E_s).
 - the forward recurrence is monic in v2[t+1], so the fwd march
   from a full-content seed window enumerates EVERY admissible
   continuation: the long backward ghost coasts of run i+1 appear
   as post-phase content driven from the block below — with the
   cross-boundary E-consistency the s6 loose IP ignored now
   enforced.
Tables are relaxations of true links (free far end, scope caps), so
deficits are overstated and floors sound.

THE PINCH LEMMA (combinatorial, proven; mechanical check in
selftest).  For a heavy block of L <= 3 slabs between light slabs b
(exit) and b+L+1 (entry): rows(slab b+1) = {b-2..b+1} is contained
in rows(slab b) u rows(slab b+L+1) = {b-3..b} u {b+L-2..b+L+1}
exactly when b+1 >= b+L-2, i.e. L <= 3.  Hence
    W_b + W_{b+L+1} >= W_{b+1} >= 8   (L <= 3):
double-vacuum interfaces are impossible at short heavy blocks.  At
L >= 4 the covering fails (an interior fat row can feed every heavy
slab), so the lemma is scoped L <= 3.
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from heapq import heappush, heappop
from pathlib import Path

LAB = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(LAB / "src"))
sys.path.insert(0, str(LAB / "scripts"))

DATA = LAB / "data" / "a40"

from a40_s6_frontier import (  # noqa: E402
    wt, lsb, norm, dilate, tooth_ok, seeds_full,
)

# phases
PRE, HEAVY, POST = 0, 1, 2


def _rss_mb():
    import resource
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss \
        // (1024 * 1024)


class LinkMarch:
    """Forward march with slab-class schedule.  Node state =
    (dyn, phase, L): dyn = (v1[t-3..t], v2[t-2..t]) normalized;
    phase in {PRE, HEAVY, POST}; L = heavy slabs used so far.

    Tables (min g):
      tab_pre[(h, dlt)]        — phase PRE (the plain fwd table;
                                 regression target vs s6)
      tab_link[(h, dlt)]       — phase POST (a completed crossing:
                                 >= 1 heavy slab, >= 1 post slab —
                                 every true link ends POST because
                                 its final slab is the next seed)
      tab_link_L[(h, L, dlt)]  — POST, heavy-length-resolved
      slip[(L, gH)]            — max |anchor movement across the
                                 heavy block| observed per (L, block
                                 weight): Stage-2 instrumentation
                                 (dlt at heavy exit minus dlt at
                                 heavy entry, carried in the node).
    """

    def __init__(self, u, kmax=2, whcap=14, smax=3, dil=4, hcap=19,
                 gcap=32, extent_cap=34, dcap=30, slips=False):
        self.do_slips = slips
        self.u, self.kmax, self.whcap = u, kmax, whcap
        self.smax, self.dil = smax, dil
        self.hcap, self.gcap = hcap, gcap
        self.extent_cap, self.dcap = extent_cap, dcap
        self.tab_pre = {}
        self.tab_link = {}
        self.tab_link_L = {}
        self.slip = {}
        self.nodes_popped = 0
        self.nodes_pushed = 0
        self.trunc_g = False
        self.trunc_extent = 0
        self.trunc_dcap = 0
        self.trunc_nodes = False

    # weight window of a slab in class `cls`
    def wwin(self, cls):
        if cls == PRE:
            return self.u, 7
        if cls == HEAVY:
            return 8, self.whcap
        return 1, 7

    def transitions(self, phase, L):
        """Allowed class of the NEXT slab given the current phase."""
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

    def record(self, phase, L, h, dlt, g, hin):
        if phase == PRE:
            k = (h, dlt)
            if k not in self.tab_pre or self.tab_pre[k] > g:
                self.tab_pre[k] = g
        elif phase == POST:
            k = (h, dlt)
            if k not in self.tab_link or self.tab_link[k] > g:
                self.tab_link[k] = g
            kL = (h, L, dlt)
            if kL not in self.tab_link_L or self.tab_link_L[kL] > g:
                self.tab_link_L[kL] = g
            # slip instrumentation: hin = (dlt at heavy entry,
            # dlt at heavy exit, block weight)
            if hin is not None:  # only when slips instrumentation on
                d_in, d_out, gH = hin
                s = abs(d_out - d_in)
                key = (L, gH)
                if s > self.slip.get(key, -1):
                    self.slip[key] = s

    def run(self, seed_states, log=True, log_every=2_000_000,
            max_nodes=120_000_000):
        u = self.u
        heap = []
        seen = {}
        n_seed = 0
        for w8 in seed_states:
            if sum(wt(r) for r in w8) != u:
                continue
            n_seed += 1
            dyn = (w8[0], w8[1], w8[2], w8[3], w8[5], w8[6], w8[7])
            allm = 0
            for r in w8:
                allm |= r
            anchor = lsb(allm)
            key = (dyn, anchor, PRE, 0, 1, 0)
            if key in seen and seen[key] <= u:
                continue
            seen[key] = u
            heappush(heap, (u, self.nodes_pushed, dyn, PRE, 0, 1, 0,
                            anchor, anchor, allm.bit_length() - 1,
                            None))
            self.nodes_pushed += 1
            self.record(PRE, 0, 1, 0, u, None)
        t0 = time.time()
        while heap:
            (g, _, dyn, phase, L, h, dlt, anch, xlo, xhi,
             hin) = heappop(heap)
            key = (dyn, anch, phase, L, h, dlt)
            if seen.get(key, 1 << 30) < g:
                continue
            self.nodes_popped += 1
            if log and self.nodes_popped % log_every == 0:
                print(f"[link u={u}] popped {self.nodes_popped} "
                      f"pushed {self.nodes_pushed} heap {len(heap)} "
                      f"g={g} {round(time.time()-t0, 1)}s",
                      flush=True)
            if self.nodes_popped > max_nodes:
                self.trunc_nodes = True
                break
            if h >= self.hcap:
                continue
            self.expand(g, dyn, phase, L, h, dlt, anch, xlo, xhi,
                        hin, heap, seen)
        return dict(n_seeds=n_seed, popped=self.nodes_popped,
                    pushed=self.nodes_pushed, trunc_g=self.trunc_g,
                    trunc_extent=self.trunc_extent,
                    trunc_dcap=self.trunc_dcap,
                    trunc_nodes=self.trunc_nodes,
                    wall_s=round(time.time() - t0, 1))

    def expand(self, g, dyn, phase, L, h, dlt, anch, xlo, xhi, hin,
               heap, seen):
        v1a, v1b, v1c, v1d, v2a, v2b, v2c = dyn
        M = 8 if (v1a | v1b | v1c | v1d | v2a | v2b | v2c) & 0xFF \
            else 0
        if M:
            v1a, v1b, v1c, v1d = (v1a << M, v1b << M, v1c << M,
                                  v1d << M)
            v2a, v2b, v2c = v2a << M, v2b << M, v2c << M
        anch_m = anch + M
        # forced v2[t+1] from E_t (identical to the s6 engine)
        v2new = v1d ^ (v1d >> 1) ^ (v1a << 1) ^ v2c ^ (v2b >> 3)
        if not tooth_ok(v1a, v1b, v1c, v2a, v2new):
            return
        fixed = v1b | v1c | v1d | v2a | v2b | v2c | v2new
        wbase = wt(v1b) + wt(v1c) + wt(v1d) + wt(v2a) + wt(v2b) \
            + wt(v2c) + wt(v2new)
        trans = self.transitions(phase, L)
        wmax_any = max(self.wwin(c)[1] for c, _ in trans)
        if wbase > wmax_any:
            return
        allow = dilate(v1a | v1b | v1c | v1d | v2a | v2b | v2c
                       | v2new, self.dil)
        acols = [i for i in range(allow.bit_length())
                 if allow >> i & 1]
        from itertools import combinations
        for cls, L2 in trans:
            wmin, wmax = self.wwin(cls)
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
                    slab = fixed | s
                    dlt2 = dlt + (lsb(slab) - anch_m)
                    if abs(dlt2) > self.dcap:
                        self.trunc_dcap += 1
                        continue
                    # slip instrumentation bookkeeping (opt-in)
                    if not self.do_slips:
                        hin2 = None
                    elif cls == HEAVY and phase == PRE:
                        hin2 = (dlt, dlt2, W)        # entry
                    elif cls == HEAVY:
                        hin2 = (hin[0], dlt2, hin[2] + W)
                    elif cls == POST and phase == HEAVY:
                        hin2 = (hin[0], hin[1], hin[2])  # frozen
                    else:
                        hin2 = hin
                    self.push(g2, (v1b, v1c, v1d, s, v2b, v2c,
                                   v2new), slab,
                              slab.bit_length() - 1, cls if cls !=
                              PRE else PRE, L2, h, dlt2, M, anch,
                              xlo, xhi, hin2, heap, seen)

    # -- layered engine (production): identical tables, no heap ----
    def run_layers(self, seed_states, log=True):
        """Level-synchronous BFS over h-layers.  Dominance = min g
        per (dyn, anch, phase, L, dlt) within each layer — exactly
        the heap engine's key (h was part of it), so the resulting
        tables are identical; memory holds only two layers."""
        u = self.u
        layer = {}
        n_seed = 0
        for w8 in seed_states:
            if sum(wt(r) for r in w8) != u:
                continue
            n_seed += 1
            dyn = (w8[0], w8[1], w8[2], w8[3], w8[5], w8[6], w8[7])
            allm = 0
            for r in w8:
                allm |= r
            anchor = lsb(allm)
            key = (dyn, anchor, PRE, 0, 0)
            cur = layer.get(key)
            if cur is None or cur[0] > u:
                layer[key] = (u, anchor, allm.bit_length() - 1)
            self.record(PRE, 0, 1, 0, u, None)
        t0 = time.time()
        h = 1
        self.complete_h = 1
        aborted = None
        while layer and h < self.hcap:
            # ops guards (checked BETWEEN layers, so layers are
            # atomic and the tables stay complete through
            # complete_h): RSS <= 2 GB, frontier <= 3M states.
            rss = _rss_mb()
            if rss > 2048:
                aborted = f"RED: RSS {rss} MB > 2048"
                break
            if len(layer) > 3_000_000:
                aborted = f"RED: frontier {len(layer)} > 3M"
                break
            nxt = {}
            for (dyn, anch, phase, L, dlt), (g, xlo, xhi) \
                    in layer.items():
                self.expand_layer(g, dyn, phase, L, h, dlt, anch,
                                  xlo, xhi, nxt)
            self.nodes_popped += len(layer)
            self.nodes_pushed += len(nxt)
            if log:
                print(f"[link u={u} layers] h={h + 1} states "
                      f"{len(nxt)} (cum {self.nodes_pushed}) "
                      f"rss {_rss_mb()}MB "
                      f"{round(time.time() - t0, 1)}s", flush=True)
            layer = nxt
            h += 1
            self.complete_h = h
        if aborted:
            self.trunc_nodes = True
            print(f"[link u={u}] ABORT {aborted} — tables complete "
                  f"through h={self.complete_h}", flush=True)
        else:
            # loop ended by hcap or by a dead frontier: the tables
            # are complete at every h (an empty layer has no
            # descendants).
            self.complete_h = self.hcap
        return dict(n_seeds=n_seed, popped=self.nodes_popped,
                    pushed=self.nodes_pushed, trunc_g=self.trunc_g,
                    trunc_extent=self.trunc_extent,
                    trunc_dcap=self.trunc_dcap,
                    trunc_nodes=self.trunc_nodes,
                    complete_h=self.complete_h,
                    aborted=aborted,
                    wall_s=round(time.time() - t0, 1))

    def expand_layer(self, g, dyn, phase, L, h, dlt, anch, xlo,
                     xhi, nxt):
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
        wmax_any = max(self.wwin(c)[1] for c, _ in trans)
        if wbase > wmax_any:
            return
        allow = dilate(v1a | v1b | v1c | v1d | v2a | v2b | v2c
                       | v2new, self.dil)
        acols = [i for i in range(allow.bit_length())
                 if allow >> i & 1]
        from itertools import combinations
        h2 = h + 1
        for cls, L2 in trans:
            wmin, wmax = self.wwin(cls)
            if wbase > wmax:
                continue
            need = max(0, wmin - wbase)
            room = wmax - wbase
            for k in range(need, min(room, self.smax) + 1):
                for pick in combinations(acols, k):
                    s = 0
                    for c in pick:
                        s |= 1 << c
                    g2 = g + wbase + k
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
                    xhi_c = max(xhi + M, slab.bit_length() - 1) - S
                    if xhi_c - xlo_c + 1 > self.extent_cap:
                        self.trunc_extent += 1
                        continue
                    key = (nrows, anch_c, cls, L2, dlt2)
                    cur = nxt.get(key)
                    if cur is not None and cur[0] <= g2:
                        continue
                    nxt[key] = (g2, xlo_c, xhi_c)
                    self.record(cls, L2, h2, dlt2, g2, None)

    def push(self, g2, state_rows, slab_mask, slab_hi, phase2, L2,
             h, dlt_new, M, anch, xlo, xhi, hin2, heap, seen):
        a_new = lsb(slab_mask)
        nrows, S = norm(state_rows)
        anch_c = a_new - S
        xlo_c = min(xlo + M, a_new) - S
        xhi_c = max(xhi + M, slab_hi) - S
        if xhi_c - xlo_c + 1 > self.extent_cap:
            self.trunc_extent += 1
            return
        h2 = h + 1
        # anch_c is part of the dominance key: the last slab's
        # anchor is NOT a function of the 7-row state (v2[t-3] was
        # dropped), and future drift increments depend on it — the
        # s6 March merged such states (anchor aliasing), which is
        # exact for delta-blind readouts but incomplete for
        # delta-resolved buckets.  Here the tables are consumed
        # delta-resolved, so the key carries the anchor.
        key = (nrows, anch_c, phase2, L2, h2, dlt_new)
        if seen.get(key, 1 << 30) <= g2:
            return
        seen[key] = g2
        self.record(phase2, L2, h2, dlt_new, g2, hin2)
        heappush(heap, (g2, self.nodes_pushed, nrows, phase2, L2,
                        h2, dlt_new, anch_c, xlo_c, xhi_c, hin2))
        self.nodes_pushed += 1


class ClosedMarch:
    """Per-seed forward march for the r in {0, 1} branches.  State =
    the full last-4-rows window (8 masks: v1[t-3..t], v2[t-3..t]),
    so the closure readout is exact: after exactly m steps, the
    normalized state equals the normalized seed window and the
    accumulated drift is == 0 (mod ell) (|delta| = ell covered:
    winding once).  Phase/L as in LinkMarch with kmax = 0 (r=0,
    all-light) or kmax >= 1 (r=1, one heavy block).  Every closed
    y-spanning walk whose minimum window has weight u is counted
    from the seed equal to that window, so per-stratum emptiness
    below the g-cap is a sound floor: g_cycle = g_march - u (the
    seed slab is counted twice)."""

    def __init__(self, u, m=18, ell=24, kmax=0, whcap=14, smax=3,
                 dil=4, gcap=40, extent_cap=34, dcap=30,
                 record_open=False, complete_prune=True):
        self.u, self.m, self.ell = u, m, ell
        self.kmax, self.whcap = kmax, whcap
        self.smax, self.dil = smax, dil
        self.gcap = gcap
        self.complete_prune = complete_prune
        self.extent_cap, self.dcap = extent_cap, dcap
        self.hcap = m + 1
        self.closed = []          # (g, phase, L, dlt, sid)
        self.record_open = record_open
        self.tab_open = {}        # (phase, h, dlt) -> min g
        self.nodes_popped = 0
        self.nodes_pushed = 0
        self.trunc_g = False
        self.trunc_extent = 0
        self.trunc_dcap = 0
        self.trunc_nodes = False

    def wwin(self, cls):
        if cls == HEAVY:
            return 8, self.whcap
        return self.u, 7          # global min: post slabs >= u too

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

    def run(self, seeds, sid0=0, log=True, log_every=2_000_000,
            max_nodes=120_000_000):
        t0 = time.time()
        info_seeds = 0
        for sid, w8 in enumerate(seeds):
            if sum(wt(r) for r in w8) != self.u:
                continue
            info_seeds += 1
            self._run_seed(sid0 + sid, w8, log, log_every, max_nodes)
            if self.trunc_nodes:
                break
        return dict(n_seeds=info_seeds, popped=self.nodes_popped,
                    pushed=self.nodes_pushed, trunc_g=self.trunc_g,
                    trunc_extent=self.trunc_extent,
                    trunc_dcap=self.trunc_dcap,
                    trunc_nodes=self.trunc_nodes,
                    closed=len(self.closed),
                    wall_s=round(time.time() - t0, 1))

    def run_layers(self, seeds, sid0=0, log=True):
        """Sequential per-seed layered engine with the ops guards
        (RSS <= 2 GB, frontier <= 3M).  Aborts are per-seed and
        recorded; a seed aborted at layer h leaves its emptiness
        claim INVALID (the run reports it)."""
        t0 = time.time()
        info_seeds = 0
        self.aborts = []
        for sid, w8 in enumerate(seeds):
            if sum(wt(r) for r in w8) != self.u:
                continue
            info_seeds += 1
            ab = self._run_seed_layers(sid0 + sid, w8, log)
            if ab:
                self.aborts.append(dict(sid=sid0 + sid, why=ab))
        return dict(n_seeds=info_seeds, popped=self.nodes_popped,
                    pushed=self.nodes_pushed, trunc_g=self.trunc_g,
                    trunc_extent=self.trunc_extent,
                    trunc_dcap=self.trunc_dcap,
                    trunc_nodes=bool(self.aborts),
                    aborts=self.aborts,
                    closed=len(self.closed),
                    wall_s=round(time.time() - t0, 1))

    def _run_seed_layers(self, sid, w8, log):
        u = self.u
        seed_norm, _ = norm(tuple(w8))
        dyn = tuple(w8)
        allm = 0
        for r in w8:
            allm |= r
        layer = {(dyn, PRE, 0, 0): (u, lsb(allm),
                                    allm.bit_length() - 1)}
        t0 = time.time()
        h = 1
        while layer and h < self.m + 1:
            rss = _rss_mb()
            if rss > 2048:
                return f"RED: RSS {rss} MB > 2048 at h={h}"
            if len(layer) > 3_000_000:
                return f"RED: frontier {len(layer)} > 3M at h={h}"
            nxt = {}
            for (dyn, phase, L, dlt), (g, xlo, xhi) in layer.items():
                self.expand_layer(g, dyn, phase, L, h, dlt, xlo,
                                  xhi, nxt)
            self.nodes_popped += len(layer)
            self.nodes_pushed += len(nxt)
            if log and (h % 6 == 0 or len(nxt) > 300000):
                print(f"[closed u={u} k{self.kmax} sid{sid} lyr] "
                      f"h={h + 1} states {len(nxt)} rss "
                      f"{_rss_mb()}MB {round(time.time()-t0, 1)}s",
                      flush=True)
            layer = nxt
            h += 1
        # readout at h = m + 1
        for (dyn, phase, L, dlt), (g, xlo, xhi) in layer.items():
            nd, _ = norm(dyn)
            if nd == seed_norm and dlt % self.ell == 0 and \
                    (phase == PRE or phase == POST):
                self.closed.append((g, phase, L, dlt, sid))
        return None

    def expand_layer(self, g, dyn, phase, L, h, dlt, xlo, xhi, nxt):
        v1a, v1b, v1c, v1d, v2z, v2a, v2b, v2c = dyn
        M = 8 if (v1a | v1b | v1c | v1d | v2z | v2a | v2b
                  | v2c) & 0xFF else 0
        if M:
            v1a, v1b, v1c, v1d = (v1a << M, v1b << M, v1c << M,
                                  v1d << M)
            v2z, v2a, v2b, v2c = (v2z << M, v2a << M, v2b << M,
                                  v2c << M)
        # the anchor of the current slab IS determined by dyn (the
        # 8 rows are the whole slab): no aliasing.
        allm = v1a | v1b | v1c | v1d | v2z | v2a | v2b | v2c
        anch_m = lsb(allm)
        v2new = v1d ^ (v1d >> 1) ^ (v1a << 1) ^ v2c ^ (v2b >> 3)
        if not tooth_ok(v1a, v1b, v1c, v2a, v2new):
            return
        fixed = v1b | v1c | v1d | v2a | v2b | v2c | v2new
        wbase = wt(v1b) + wt(v1c) + wt(v1d) + wt(v2a) + wt(v2b) \
            + wt(v2c) + wt(v2new)
        trans = self.transitions(phase, L)
        wmax_any = max(self.wwin(c)[1] for c, _ in trans)
        if wbase > wmax_any:
            return
        allow = dilate(v1a | v1b | v1c | v1d | v2a | v2b | v2c
                       | v2new, self.dil)
        acols = [i for i in range(allow.bit_length())
                 if allow >> i & 1]
        from itertools import combinations
        for cls, L2 in trans:
            wmin, wmax = self.wwin(cls)
            if wbase > wmax:
                continue
            need = max(0, wmin - wbase)
            room = wmax - wbase
            for k in range(need, min(room, self.smax) + 1):
                for pick in combinations(acols, k):
                    s = 0
                    for c in pick:
                        s |= 1 << c
                    g2 = g + wbase + k
                    rem = (self.m - h) * self.u \
                        if self.complete_prune else 0
                    if g2 + rem > self.gcap:
                        self.trunc_g = True
                        continue
                    slab = fixed | s
                    dlt2 = dlt + (lsb(slab) - anch_m)
                    if abs(dlt2) > self.dcap:
                        self.trunc_dcap += 1
                        continue
                    state = (v1b, v1c, v1d, s, v2a, v2b, v2c, v2new)
                    a_new = lsb(slab)
                    nrows, S = norm(state)
                    xlo_c = min(xlo + M, a_new) - S
                    xhi_c = max(xhi + M, slab.bit_length() - 1) - S
                    if xhi_c - xlo_c + 1 > self.extent_cap:
                        self.trunc_extent += 1
                        continue
                    key = (nrows, cls, L2, dlt2)
                    cur = nxt.get(key)
                    if cur is not None and cur[0] <= g2:
                        continue
                    nxt[key] = (g2, xlo_c, xhi_c)

    def _run_seed(self, sid, w8, log, log_every, max_nodes):
        u = self.u
        seed_norm, _ = norm(tuple(w8))
        heap = []
        seen = {}
        dyn = tuple(w8)           # (v1[0..3], v2[0..3]) — 8 masks
        allm = 0
        for r in w8:
            allm |= r
        anchor = lsb(allm)
        key = (dyn, PRE, 0, 1, 0)
        seen[key] = u
        heappush(heap, (u, self.nodes_pushed, dyn, PRE, 0, 1, 0,
                        anchor, anchor, allm.bit_length() - 1))
        self.nodes_pushed += 1
        t0 = time.time()
        while heap:
            g, _, dyn, phase, L, h, dlt, anch, xlo, xhi = \
                heappop(heap)
            key = (dyn, phase, L, h, dlt)
            if seen.get(key, 1 << 30) < g:
                continue
            self.nodes_popped += 1
            if log and self.nodes_popped % log_every == 0:
                print(f"[closed u={u} k{self.kmax} sid{sid}] popped "
                      f"{self.nodes_popped} heap {len(heap)} g={g} "
                      f"{round(time.time()-t0, 1)}s", flush=True)
            if self.nodes_popped > max_nodes:
                self.trunc_nodes = True
                return
            if h == self.m + 1:
                # readout: exactly m steps done; closure?
                nd, _ = norm(dyn)
                if nd == seed_norm and dlt % self.ell == 0 and \
                        (phase == PRE or phase == POST):
                    self.closed.append((g, phase, L, dlt, sid))
                continue
            self.expand(g, dyn, phase, L, h, dlt, anch, xlo, xhi,
                        heap, seen)

    def expand(self, g, dyn, phase, L, h, dlt, anch, xlo, xhi, heap,
               seen):
        v1a, v1b, v1c, v1d, v2z, v2a, v2b, v2c = dyn
        M = 8 if (v1a | v1b | v1c | v1d | v2z | v2a | v2b
                  | v2c) & 0xFF else 0
        if M:
            v1a, v1b, v1c, v1d = (v1a << M, v1b << M, v1c << M,
                                  v1d << M)
            v2z, v2a, v2b, v2c = (v2z << M, v2a << M, v2b << M,
                                  v2c << M)
        anch_m = anch + M
        v2new = v1d ^ (v1d >> 1) ^ (v1a << 1) ^ v2c ^ (v2b >> 3)
        if not tooth_ok(v1a, v1b, v1c, v2a, v2new):
            return
        fixed = v1b | v1c | v1d | v2a | v2b | v2c | v2new
        wbase = wt(v1b) + wt(v1c) + wt(v1d) + wt(v2a) + wt(v2b) \
            + wt(v2c) + wt(v2new)
        trans = self.transitions(phase, L)
        wmax_any = max(self.wwin(c)[1] for c, _ in trans)
        if wbase > wmax_any:
            return
        allow = dilate(v1a | v1b | v1c | v1d | v2a | v2b | v2c
                       | v2new, self.dil)
        acols = [i for i in range(allow.bit_length())
                 if allow >> i & 1]
        from itertools import combinations
        for cls, L2 in trans:
            wmin, wmax = self.wwin(cls)
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
                    # completion prune (exact, not a scope change):
                    # each of the remaining m+1-(h+1) slabs costs
                    # >= u, so a closed walk through this state
                    # already exceeds the cap.
                    rem = (self.m - h) * self.u \
                        if self.complete_prune else 0
                    if g2 + rem > self.gcap:
                        self.trunc_g = True
                        continue
                    slab = fixed | s
                    dlt2 = dlt + (lsb(slab) - anch_m)
                    if abs(dlt2) > self.dcap:
                        self.trunc_dcap += 1
                        continue
                    state = (v1b, v1c, v1d, s, v2a, v2b, v2c, v2new)
                    a_new = lsb(slab)
                    nrows, S = norm(state)
                    anch_c = a_new - S
                    xlo_c = min(xlo + M, a_new) - S
                    xhi_c = max(xhi + M, slab.bit_length() - 1) - S
                    if xhi_c - xlo_c + 1 > self.extent_cap:
                        self.trunc_extent += 1
                        continue
                    h2 = h + 1
                    key = (nrows, cls, L2, h2, dlt2)
                    if seen.get(key, 1 << 30) <= g2:
                        continue
                    seen[key] = g2
                    if self.record_open:
                        ko = (cls, h2, dlt2)
                        if self.tab_open.get(ko, 1 << 30) > g2:
                            self.tab_open[ko] = g2
                    heappush(heap, (g2, self.nodes_pushed, nrows,
                                    cls, L2, h2, dlt2, anch_c,
                                    xlo_c, xhi_c))
                    self.nodes_pushed += 1


# ---------------------------------------------------------------------
# selftest: regression, pinch, small link run
# ---------------------------------------------------------------------

def selftest():
    t0 = time.time()
    print("SELFTEST 1: pinch-lemma row covering (mechanical) ...",
          flush=True)
    # slab rows: slab j = rows [j-3, j].  For heavy block b+1..b+L
    # between light slabs b and b+L+1, every row of slab b+1 must
    # lie in rows(slab b) u rows(slab b+L+1) iff L <= 3.
    for L in range(1, 7):
        rows_b1 = set(range(-2, 2))                 # slab b+1, b=0
        rows_exit = set(range(-3, 1))               # slab b
        rows_entry = set(range(L - 2, L + 2))       # slab b+L+1
        covered = rows_b1 <= (rows_exit | rows_entry)
        assert covered == (L <= 3), (L, covered)
    print("  row covering holds exactly for L <= 3: PASS",
          flush=True)

    print("SELFTEST 2: kmax=0 regression vs s6 fwd table (u=1, "
          "gcap=24) ...", flush=True)
    m = LinkMarch(1, kmax=0, gcap=24, hcap=19)
    info = m.run(seeds_full(1), log=False)
    ref = json.loads((DATA / "s6_frontier_u1prod.json").read_text())
    ref_tab = {tuple(map(int, k.split(","))): g
               for k, g in ref["fwd"]["table"].items()}
    # the s7 anchor-keyed dominance is FINER than s6's (which
    # aliased anchors): s7 must dominate every s6 bucket and agree
    # on the delta-blind per-h minima (where s6 is exact).
    for k, g in ref_tab.items():
        assert k in m.tab_pre and m.tab_pre[k] <= g, (k, g)
    for h in {h for (h, d) in ref_tab}:
        m6 = min(g for (hh, d), g in ref_tab.items() if hh == h)
        m7 = min(g for (hh, d), g in m.tab_pre.items() if hh == h)
        assert m6 == m7, (h, m6, m7)
    extra = len(m.tab_pre) - len(ref_tab)
    assert not m.tab_link
    print(f"  dominates s6 ({len(ref_tab)} buckets, +{extra} "
          f"de-aliased, per-h minima equal): PASS", flush=True)

    print("SELFTEST 3: small link run (u=1, kmax=1, gcap=24) + "
          "sanity ...", flush=True)
    m = LinkMarch(1, kmax=1, whcap=14, gcap=24, hcap=19)
    info = m.run(seeds_full(1), log=False)
    print(f"  {info}", flush=True)
    print(f"  link buckets: {len(m.tab_link)}", flush=True)
    # sanity: every link bucket has g >= (h-2)*1 + 8 + 1 (h-2 light
    # slabs at >= 1, one heavy >= 8, one post >= 1  — h counts seed)
    bad = [(h, d, g) for (h, d), g in m.tab_link.items()
           if g < (h - 2) * 1 + 8 + 1]
    assert not bad, bad[:5]
    # sanity: no link at h < 3 (seed + heavy + >= 1 post)
    assert all(h >= 3 for (h, d) in m.tab_link), \
        sorted(m.tab_link)[:5]
    print("  weight/shape sanity: PASS", flush=True)

    print("SELFTEST 4: closed march smoke (u=1, all-light, "
          "gcap=26) ...", flush=True)
    c = ClosedMarch(1, m=18, ell=24, kmax=0, gcap=26)
    info = c.run(seeds_full(1), log=False)
    print(f"  {info}", flush=True)
    print(f"  closed walks found: {len(c.closed)} (expect 0 at "
          f"this cap)", flush=True)

    print(f"selftest done ({round(time.time()-t0, 1)} s)",
          flush=True)


# ---------------------------------------------------------------------
# production runs
# ---------------------------------------------------------------------

def run_link(args):
    import gc
    gc.disable()
    t0 = time.time()
    sds = seeds_full(args.u)
    print(f"link u={args.u} kmax={args.kmax} whcap={args.whcap} "
          f"gcap={args.gcap} smax={args.smax} dil={args.dil}: "
          f"{len(sds)} raw seeds", flush=True)
    if args.byseed:
        # sequential per-seed layered runs (memory-bounded);
        # dominance across seeds is only an optimization, so the
        # merged tables are identical to the joint run's.
        m = LinkMarch(args.u, kmax=args.kmax, whcap=args.whcap,
                      smax=args.smax, dil=args.dil, hcap=args.hcap,
                      gcap=args.gcap, dcap=args.dcap)
        infos = []
        complete_h = args.hcap
        from a40_s6_frontier import wt as _wt
        good = [w8 for w8 in sds
                if sum(_wt(r) for r in w8) == args.u]
        for i, w8 in enumerate(good):
            mi = LinkMarch(args.u, kmax=args.kmax,
                           whcap=args.whcap, smax=args.smax,
                           dil=args.dil, hcap=args.hcap,
                           gcap=args.gcap, dcap=args.dcap)
            inf = mi.run_layers([w8], log=False)
            infos.append(inf)
            complete_h = min(complete_h, inf["complete_h"])
            for tab, mtab in ((mi.tab_pre, m.tab_pre),
                              (mi.tab_link, m.tab_link),
                              (mi.tab_link_L, m.tab_link_L)):
                for k, g in tab.items():
                    if k not in mtab or mtab[k] > g:
                        mtab[k] = g
            m.trunc_g = m.trunc_g or mi.trunc_g
            m.trunc_extent += mi.trunc_extent
            m.trunc_dcap += mi.trunc_dcap
            m.trunc_nodes = m.trunc_nodes or mi.trunc_nodes
            m.nodes_pushed += inf["pushed"]
            print(f"  seed {i + 1}/{len(good)}: pushed "
                  f"{inf['pushed']} complete_h "
                  f"{inf['complete_h']}"
                  + (f" ABORT {inf['aborted']}"
                     if inf['aborted'] else ""), flush=True)
        info = dict(n_seeds=len(good), popped=m.nodes_pushed,
                    pushed=m.nodes_pushed, trunc_g=m.trunc_g,
                    trunc_extent=m.trunc_extent,
                    trunc_dcap=m.trunc_dcap,
                    trunc_nodes=m.trunc_nodes,
                    complete_h=complete_h,
                    aborted=[i_["aborted"] for i_ in infos
                             if i_["aborted"]] or None,
                    wall_s=round(sum(i_["wall_s"]
                                     for i_ in infos), 1))
    else:
        m = LinkMarch(args.u, kmax=args.kmax, whcap=args.whcap,
                      smax=args.smax, dil=args.dil, hcap=args.hcap,
                      gcap=args.gcap, dcap=args.dcap)
        info = m.run_layers(sds)
    print(f"info: {info}", flush=True)
    print(f"pre buckets {len(m.tab_pre)}, link buckets "
          f"{len(m.tab_link)}", flush=True)
    best = {}
    for (h, d), g in m.tab_link.items():
        D = 2 * h - g / 4
        if D > best.get(h, (-99, None, None))[0]:
            best[h] = (D, d, g)
    for h, v in sorted(best.items()):
        print(f"  link h={h}: D={v[0]:.2f} delta={v[1]} g={v[2]}",
              flush=True)
    out = dict(
        params=dict(u=args.u, kmax=args.kmax, whcap=args.whcap,
                    gcap=args.gcap, hcap=args.hcap, smax=args.smax,
                    dil=args.dil, dcap=args.dcap),
        info=info,
        tab_pre={f"{h},{d}": g for (h, d), g
                 in sorted(m.tab_pre.items())},
        tab_link={f"{h},{d}": g for (h, d), g
                  in sorted(m.tab_link.items())},
        tab_link_L={f"{h},{L},{d}": g for (h, L, d), g
                    in sorted(m.tab_link_L.items())},
        slip={f"{L},{gH}": s for (L, gH), s
              in sorted(m.slip.items())},
        best_link_D={h: dict(D=round(v[0], 3), delta=v[1], g=v[2])
                     for h, v in sorted(best.items())})
    p = DATA / f"s7_link_{args.tag or f'u{args.u}_g{args.gcap}'}.json"
    p.write_text(json.dumps(out, indent=1))
    print(f"wrote {p} ({round(time.time()-t0, 1)} s)", flush=True)


def run_closed(args):
    import gc
    gc.disable()
    t0 = time.time()
    sds = seeds_full(args.u)
    sds = [w8 for w8 in sds if sum(wt(r) for r in w8) == args.u]
    lo, hi = 0, len(sds)
    if args.seeds:
        lo, hi = map(int, args.seeds.split(":"))
    print(f"closed u={args.u} m={args.m} ell={args.ell} "
          f"kmax={args.kmax} gcap={args.gcap}: seeds [{lo}, {hi}) "
          f"of {len(sds)}", flush=True)
    c = ClosedMarch(args.u, m=args.m, ell=args.ell, kmax=args.kmax,
                    whcap=args.whcap, smax=args.smax, dil=args.dil,
                    gcap=args.gcap, dcap=args.dcap)
    info = c.run_layers(sds[lo:hi], sid0=lo)
    print(f"info: {info}", flush=True)
    for (g, phase, L, dlt, sid) in sorted(c.closed)[:20]:
        print(f"  CLOSED: g={g} phase={phase} L={L} dlt={dlt} "
              f"sid={sid}", flush=True)
    out = dict(
        params=dict(u=args.u, m=args.m, ell=args.ell,
                    kmax=args.kmax, whcap=args.whcap,
                    gcap=args.gcap, smax=args.smax, dil=args.dil,
                    dcap=args.dcap, seeds=[lo, hi]),
        info=info,
        closed=[dict(g=g, phase=p_, L=L, dlt=d, sid=s)
                for (g, p_, L, d, s) in sorted(c.closed)])
    p = DATA / (f"s7_closed_{args.tag or f'u{args.u}_k{args.kmax}'}"
                f"_g{args.gcap}_s{lo}_{hi}.json")
    p.write_text(json.dumps(out, indent=1))
    print(f"wrote {p} ({round(time.time()-t0, 1)} s)", flush=True)


def control():
    """End-to-end validation of the CLOSED lane on a real closed
    walk: the L12 rate-2 species at (18,6) (drift 0, crosses heavy,
    d((18,6)) = 12 certified).  The march, seeded at the walk's own
    minimum window, must detect its closure at exactly
    g = sum(slabs) + u with dlt == 0 (mod 18).  Also cross-checks
    ClosedMarch's expansion against the regression-validated
    LinkMarch on open tables at equal caps."""
    t0 = time.time()
    print("CONTROL 1: ClosedMarch-vs-LinkMarch open-table "
          "cross-check (u=1, kmax=1, gcap=22) ...", flush=True)
    lm = LinkMarch(1, kmax=1, whcap=14, gcap=22, hcap=19)
    lm.run(seeds_full(1), log=False)
    cm = ClosedMarch(1, m=18, ell=10**9, kmax=1, whcap=14, gcap=22,
                     record_open=True, complete_prune=False)
    cm.run(seeds_full(1), log=False)
    pre_c = {(h, d): g for (ph, h, d), g in cm.tab_open.items()
             if ph == PRE}
    post_c = {(h, d): g for (ph, h, d), g in cm.tab_open.items()
              if ph == POST}
    pre_l = {k: v for k, v in lm.tab_pre.items() if k != (1, 0)}
    assert pre_c == pre_l, (
        len(pre_c), len(pre_l), "PRE tables differ")
    assert post_c == lm.tab_link, (
        len(post_c), len(lm.tab_link), "POST tables differ")
    print(f"  identical open tables (PRE {len(pre_c)}, POST "
          f"{len(post_c)} buckets): PASS", flush=True)

    print("CONTROL 2: L12 closure detection at (18,6) ...",
          flush=True)
    sys.argv = [sys.argv[0]]
    import numpy as np  # noqa: F401
    from a40_s6_drift import lift_phase
    from a40_s5_lightcore import Phase
    from a40_s4_phase_atlas import atlas
    rows_at, _ = atlas("AB", 6, 12, keep_pts=True)
    cands = [r for r in rows_at
             if r["nontrivial"] and r["weight"] == 12
             and r["extent"] <= 13]
    done = False
    for cand in cands:
        ph = Phase(18, 6, 0,
                   [(blk, c, y) for (c, y, blk) in cand["pts"]])
        fr, s = lift_phase(ph, n_periods=3)
        if s != 0:
            continue
        sl_all = fr.slabs()
        prof = sl_all[3:9]        # one cyclic period, mid-lift
        # cyclic heavy blocks
        hv = [i for i, w in enumerate(prof) if w >= 8]
        if not hv:
            continue
        blocks = []
        cur = [hv[0]]
        for i in hv[1:]:
            if i == cur[-1] + 1:
                cur.append(i)
            else:
                blocks.append(cur)
                cur = [i]
        blocks.append(cur)
        if prof[0] >= 8 and prof[-1] >= 8 and len(blocks) > 1:
            blocks[0] = blocks.pop() + blocks[0]
        if len(blocks) != 1 or len(blocks[0]) > 2:
            continue
        Lblk = len(blocks[0])
        # seed at the min-weight slab of the period
        lights = [i for i, w in enumerate(prof) if w <= 7]
        j0 = min(lights, key=lambda i: prof[i])
        u0 = prof[j0]
        # absolute slab index in the lift: slabs()[k] is slab of
        # rows [k, k+3] + t0... slab index k covers rows k..k+3
        # (t0=0, slabs start at j=t0+3: slabs()[k] = slab j=k+3,
        # rows k..k+3).  Profile index i -> lift slab index 3+i.
        kk = 3 + j0
        w8 = []
        for blkk in (0, 1):
            for rr in range(kk, kk + 4):
                w8.append(sorted(fr.s(blkk, rr)))
        cols = [c for row in w8 for c in row]
        c0 = min(cols)
        w8m = tuple(
            sum(1 << (c - c0) for c in row) for row in w8)
        gsum = sum(prof)
        print(f"  L12 member: profile {prof}, block L={Lblk}, "
              f"seed slab {j0} (u={u0}), g_target={gsum + u0}",
              flush=True)
        c2 = ClosedMarch(u0, m=6, ell=18, kmax=Lblk,
                         whcap=max(prof), smax=8, dil=8,
                         gcap=gsum + u0, dcap=30)
        info = c2.run([w8m], log=False)
        hits = [(g, phse, L, d) for (g, phse, L, d, _) in c2.closed]
        print(f"  march info: {info}; closed hits: {hits}",
              flush=True)
        ok = any(g <= gsum + u0 and d % 18 == 0
                 for (g, phse, L, d) in hits)
        assert ok, "L12 closure NOT detected"
        exact = any(g == gsum + u0 and d == 0
                    for (g, phse, L, d) in hits)
        print(f"  closure detected (exact-g hit: {exact}): PASS",
              flush=True)
        done = True
        break
    assert done, "no single-block L12 member found for the control"
    print(f"control done ({round(time.time()-t0, 1)} s)",
          flush=True)


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("selftest")
    sub.add_parser("control")
    for name in ("link", "closed"):
        s = sub.add_parser(name)
        s.add_argument("--u", type=int, required=True)
        s.add_argument("--gcap", type=int, required=True)
        s.add_argument("--kmax", type=int, default=2)
        s.add_argument("--whcap", type=int, default=14)
        s.add_argument("--smax", type=int, default=3)
        s.add_argument("--dil", type=int, default=4)
        s.add_argument("--hcap", type=int, default=19)
        s.add_argument("--dcap", type=int, default=30)
        s.add_argument("--maxnodes", type=int, default=120_000_000)
        s.add_argument("--byseed", action="store_true")
        s.add_argument("--tag", type=str, default="")
        if name == "closed":
            s.add_argument("--m", type=int, default=18)
            s.add_argument("--ell", type=int, default=24)
            s.add_argument("--seeds", type=str, default="")
    args = ap.parse_args()
    if args.cmd == "selftest":
        selftest()
    elif args.cmd == "control":
        control()
    elif args.cmd == "link":
        run_link(args)
    else:
        run_closed(args)


if __name__ == "__main__":
    main()
