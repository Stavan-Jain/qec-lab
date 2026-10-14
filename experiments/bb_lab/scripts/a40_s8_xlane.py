#!/usr/bin/env python3
"""A40 S8 — Stage 1: the x-sector mirror (the theta'-image lane).

THE MIRROR REDUCTION (exact, weight-preserving).  theta': (x, y) |->
(y^-1, x) is a ring automorphism of F2[x^pm, y^pm] carrying the pair
(A, B) to (A o theta', B o theta') = (B, Abar): the AB code on
Z_l x Z_m is isomorphic to the BAbar code on Z_m x Z_l (support map
(a, b) |-> (b, -a), blocks preserved, weights/cycles/stabilizers/
nontriviality preserved).  Under this map the X-SECTOR of AB-(l, m)
— nontrivial X-logicals with a cyclic y-gap >= 4, which by Lemma K
have NO x-gap >= 4 (K transports: theta' is a ring automorphism, so
(B, Abar) is plane-regular because (A, B) is) — corresponds exactly
to the y-walk sector of BAbar-(m, l) whose CONTENT has a cyclic
x-gap >= 4:

  - walk rows = the l y-rows of BAbar-(m, l) (= x-columns of the
    original); every 4-row slab NONEMPTY (no original x-gap >= 4);
  - content = x-columns on an INTERVAL of extent <= m - 4 (the
    original y-gap >= 4 cuts the content circle): the interval
    embedding is the canonical cover lift, there is NO wrapped
    corner and NO winding branch, and closure demands total drift
    == 0 EXACTLY (|dlt| <= extent - 1 <= m - 5 < m anyway);
  - slab telescope: 4|v| = sum of the l slab weights (cyclic in the
    walk direction).

MIRROR RECURRENCE (derived generically from the pair supports in
selftest; hand form): X-cycles of (P', Q') = (B, Abar) satisfy
Q'bar u1 + P'bar u2 = A u1 + Bbar u2 = 0, i.e. per walk row t
  E'_t: u1[t] + u1[t-1] + x^3 u1[t+1] + (1+x^-1) u2[t] + x u2[t-3] = 0
— monic in u1[t+1] up to the unit x^3: BLOCK 1 is forced forward
(u1[t+1] = x^-3 (...)), block-2 rows are the free inputs (S4 §9.1's
mirrored convolutional structure, offset 3).  Footprint spans (4,4);
the X-stabilizer generator ("tooth") is (Bbar, A):
  blk1 (cx,cy), (cx-1,cy), (cx+1,cy+3);
  blk2 (cx,cy), (cx,cy+1), (cx+3,cy-1).

COUNTING FLOORS (exact, from the telescope; n_H = # heavy slabs,
walk length l): 4|v| >= l + 7 n_H, and min-slab u >= 2 gives
4|v| >= 2 l.  At (l, m) = (24, 18), target floor 11 (resp. 12), only
u = 1 with n_H <= 2 (<= 3) and per-heavy W <= 19 (23) survives the
counting — enumerated by per-seed closed marches (full-window state,
alias-free) at gcap 42 (46) with the completion prune.  March
emptiness at gcap G gives 4|v| + u > G (seed slab counted twice).

SOUNDNESS NOTES.  (i) The extent cap is the SECTOR DEFINITION here
(not a scope condition), so the running envelope (xlo, xhi) is part
of the dominance key — merging paths with different envelopes under
a binding extent cap would under-enumerate.  (ii) Closure = exact
dlt == 0 (the interval lift of any x-sector object has zero total
anchor drift; winding is unreachable below extent m - 4).  (iii)
The window(tooth) prune is sound for class-minimal logicals (local
reduction, §9.2, pair-generic).  (iv) Scope conditions that DO
remain: smax new points per input row, dilation-radius growth
(prefix-connected placement) — production smax 3 / dil 4, stability
at 4 / 6, same discipline as S7.

Runs: selftest | control | coherence | prod (see main()).
"""
from __future__ import annotations

import argparse
import json
import random
import sys
import time
from itertools import combinations
from pathlib import Path

LAB = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(LAB / "src"))
sys.path.insert(0, str(LAB / "scripts"))

DATA = LAB / "data" / "a40"

from a40_s6_frontier import (  # noqa: E402
    wt, lsb, norm, dilate, seeds_full,
)
from a40_s4_phase_atlas import (  # noqa: E402
    A_L, B_L, bar, atlas,
)

# the mirror pair (P', Q') = (B, Abar); cycle eq Q'bar u1 + P'bar u2
P_M = list(B_L)
Q_M = bar(A_L)
TERMS_U1 = bar(Q_M)         # = A_L: [(0,0), (0,1), (3,-1)]
TERMS_U2 = bar(P_M)         # = [(0,0), (-1,0), (1,3)]
# generator (tooth) of the mirror system: (P'bar, Q'bar) = (Bbar, A)
TOOTH_B1 = bar(P_M)         # blk1 offsets: (0,0), (-1,0), (1,3)
TOOTH_B2 = bar(Q_M)         # = A_L: (0,0), (0,1), (3,-1)


def bit(mask, i):
    return (mask >> i) & 1 if i >= 0 else 0


def _rss_mb():
    """CURRENT process RSS in MB (ps-based; ru_maxrss is a
    lifetime peak and would trip per-seed guards forever after one
    spike — the R1 incident, §13.6)."""
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
# the mirror cover fragment (generic from the supports — the
# independent verifier for everything the marches produce)
# ---------------------------------------------------------------------

class MirrorFragment:
    """Walk rows on the y-cover of the mirror (BAbar) system:
    rows[j] = (frozenset u1, frozenset u2) of INTEGER content
    columns, j in [t0, t1]."""

    def __init__(self, rows, t0):
        self.t0 = t0
        self.rows = rows
        self.t1 = t0 + len(rows) - 1

    def s(self, blk, j):
        if j < self.t0 or j > self.t1:
            return frozenset()
        return self.rows[j - self.t0][blk]

    def weight(self):
        return sum(len(a) + len(b) for a, b in self.rows)

    def e_residual(self, j):
        """Support of E'_j, built GENERICALLY from the term lists:
        term (a, b) of block i reads u_i[j - b] shifted by +a."""
        acc = set()

        def add(cols):
            for c in cols:
                if c in acc:
                    acc.remove(c)
                else:
                    acc.add(c)
        for (a, b) in TERMS_U1:
            add(c + a for c in self.s(0, j - b))
        for (a, b) in TERMS_U2:
            add(c + a for c in self.s(1, j - b))
        return acc

    def admissible(self):
        """E'_j = 0 for every fully-supported j: E'_j reads rows
        j-3 .. j+1, so j in [t0+3, t1-1]."""
        return all(not self.e_residual(j)
                   for j in range(self.t0 + 3, self.t1))

    def slab(self, j):
        return sum(len(self.s(b, t)) for b in (0, 1)
                   for t in range(j - 3, j + 1))

    def slabs(self):
        return [self.slab(j) for j in range(self.t0 + 3, self.t1 + 1)]

    def anchor(self, j):
        cols = [c for b in (0, 1) for t in range(j - 3, j + 1)
                for c in self.s(b, t)]
        return min(cols) if cols else None

    def anchors(self):
        return [self.anchor(j) for j in range(self.t0 + 3, self.t1 + 1)]

    def drift(self):
        a = self.anchors()
        assert a and a[0] is not None and a[-1] is not None
        return a[-1] - a[0]

    def extent(self):
        cols = [c for a, b in self.rows for c in a | b]
        return (max(cols) - min(cols) + 1) if cols else 0

    def window_prune_events(self):
        """Fully-visible mirror-tooth alignments held > half (6
        cells; overlap > 3).  Tooth rows: blk1 cy, cy+3; blk2
        cy-1, cy, cy+1 — visible iff rows cy-1..cy+3 inside."""
        ev = []
        cols = sorted({c for a, b in self.rows for c in a | b})
        if not cols:
            return ev
        for cy in range(self.t0 + 1, self.t1 - 2):
            for cx in range(cols[0] - 4, cols[-1] + 5):
                cells = [(0, cx + a, cy + b) for (a, b) in TOOTH_B1] \
                    + [(1, cx + a, cy + b) for (a, b) in TOOTH_B2]
                ov = sum(1 for blk, x, y in cells
                         if x in self.s(blk, y))
                if 2 * ov > 6:
                    ev.append((cy, cx, ov))
        return ev

    def subfragment(self, a, b):
        return MirrorFragment(
            self.rows[a - self.t0:b - self.t0 + 1], a)


def forced_row_sets(u1_t, u1_tm1, u2_t, u2_tm3):
    """Generic set-arithmetic solve of E'_t for u1[t+1] =
    x^-3 (u1[t] + u1[t-1] + (1+x^-1) u2[t] + x u2[t-3])."""
    acc = set()

    def add(cols):
        for c in cols:
            if c in acc:
                acc.remove(c)
            else:
                acc.add(c)
    add(u1_t)
    add(u1_tm1)
    add(u2_t)
    add(c - 1 for c in u2_t)
    add(c + 1 for c in u2_tm3)
    return frozenset(c - 3 for c in acc)


# ---------------------------------------------------------------------
# the march kernel (bitmask form; checked against the generic solve
# in selftest)
# ---------------------------------------------------------------------

def m_forced(p1b, p1c, p2z, p2c):
    """u1new in the MARGIN frame (caller pre-shifts so all content
    sits at bits >= 8 or has none below bit 8; min result bit >= 4,
    no underflow)."""
    return (p1c ^ p1b ^ p2c ^ (p2c >> 1) ^ (p2z << 1)) >> 3


def m_tooth_ok(p1a, u1new, p2z, p2a, p2b):
    """Mirror tooth at cy = t-2 (decided by the forced row alone):
    cells blk1 (cx, cx-1) on row t-2 [p1a], (cx+1) on row t+1
    [u1new]; blk2 (cx+3) on row t-3 [p2z], (cx) on rows t-2, t-1
    [p2a, p2b].  Require overlap <= 3 at every alignment."""
    m = p1a | u1new | p2z | p2a | p2b
    if m == 0:
        return True
    lo, hi = lsb(m) - 4, m.bit_length() + 4
    for cx in range(lo, hi + 1):
        ov = bit(p1a, cx) + bit(p1a, cx - 1) + bit(u1new, cx + 1) \
            + bit(p2z, cx + 3) + bit(p2a, cx) + bit(p2b, cx)
        if ov > 3:
            return False
    return True


LIGHT, HEAVY = 0, 1


class MClosedMarch:
    """Per-seed forward closed march of the mirror system.

    State = the full last-4-rows window (8 masks: u1[t-3..t],
    u2[t-3..t]) — anchors alias-free; dominance key =
    (state, nH, dlt, xlo, xhi): the envelope is IN the key because
    the extent cap is the sector definition (see module docstring).
    Slab classes: light W in [u, 7] (u = the walk minimum: every
    slab >= u), heavy W in [8, whcap], at most nHmax heavy slabs in
    ANY arrangement (adjacent blocks, separated blocks).  Readout
    after exactly m steps: normalized state == normalized seed and
    dlt == 0 (exact_dlt) or dlt == 0 mod ell (control mode).
    Emptiness below gcap gives 4|v| + u > gcap for the covered
    branch set."""

    def __init__(self, u, m, ell, nHmax=0, whcap=19, smax=3, dil=4,
                 gcap=42, extent_cap=14, dcap=20, exact_dlt=True):
        self.u, self.m, self.ell = u, m, ell
        self.nHmax, self.whcap = nHmax, whcap
        self.smax, self.dil = smax, dil
        self.gcap = gcap
        self.extent_cap, self.dcap = extent_cap, dcap
        self.exact_dlt = exact_dlt
        self.closed = []          # dicts
        self.nodes_pushed = 0
        self.nodes_popped = 0
        self.trunc_g = False
        self.trunc_extent = 0
        self.trunc_dcap = 0
        self.aborts = []
        self.layer_peak = 0

    def run_layers(self, seeds, sid0=0, log=True, rss_cap=2048,
                   frontier_cap=3_000_000):
        t0 = time.time()
        nseed = 0
        for sid, w8 in enumerate(seeds):
            if sum(wt(r) for r in w8) != self.u:
                continue
            nseed += 1
            ab = self._run_seed(sid0 + sid, w8, log, rss_cap,
                                frontier_cap)
            if ab:
                self.aborts.append(dict(sid=sid0 + sid, why=ab))
        return dict(n_seeds=nseed, pushed=self.nodes_pushed,
                    popped=self.nodes_popped, trunc_g=self.trunc_g,
                    trunc_extent=self.trunc_extent,
                    trunc_dcap=self.trunc_dcap,
                    aborts=self.aborts, closed=len(self.closed),
                    layer_peak=self.layer_peak,
                    wall_s=round(time.time() - t0, 1))

    SPLIT_THRESH = 180_000

    def _run_seed(self, sid, w8, log, rss_cap, frontier_cap):
        u = self.u
        seed_norm, _ = norm(tuple(w8))
        allm = 0
        for r in w8:
            allm |= r
        layer = {(tuple(w8), 0, 0, lsb(allm),
                  allm.bit_length() - 1): u}
        return self._run_layerset(layer, 1, sid, seed_norm, log,
                                  rss_cap, frontier_cap)

    # recursive layer chunking (see a40_s8_slip): a layer above
    # SPLIT_THRESH is partitioned and each chunk marched to its own
    # readout.  Sound: cross-chunk dominance is only an
    # optimization (over-exploration), closures merge in
    # self.closed; guards run per chunk with a mid-expansion RSS
    # check.
    def _run_layerset(self, layer, h, sid, seed_norm, log, rss_cap,
                      frontier_cap, depth=0):
        import os
        import pickle
        import tempfile
        u = self.u
        t0 = time.time()
        while layer and h < self.m + 1:
            rss = _rss_mb()
            if rss > rss_cap:
                return f"RED: RSS {rss} MB > {rss_cap} at h={h}"
            nxt = {}
            parts = []

            def spill(d):
                fd, pth = tempfile.mkstemp(
                    suffix=f"_s8_sid{sid}_d{depth}_h{h}.pkl")
                with os.fdopen(fd, "wb") as f:
                    pickle.dump(list(d.items()), f, protocol=4)
                parts.append(pth)

            nexp = 0
            for (dyn, nH, dlt, xlo, xhi), g in layer.items():
                self.expand(g, dyn, nH, h, dlt, xlo, xhi, nxt)
                nexp += 1
                if len(nxt) > self.SPLIT_THRESH:
                    spill(nxt)
                    nxt = {}
                if nexp % 65536 == 0 and _rss_mb() > rss_cap:
                    for pth in parts:
                        os.remove(pth)
                    return f"RED: RSS {_rss_mb()} MB mid-layer " \
                        f"h={h}"
            self.nodes_popped += len(layer)
            if parts:
                if nxt:
                    spill(nxt)
                nxt = None
                layer = None
                if log:
                    print(f"[s8 sid{sid} d{depth}] h={h + 1} "
                          f"streamed into {len(parts)} parts rss "
                          f"{_rss_mb()}MB "
                          f"{round(time.time() - t0, 1)}s",
                          flush=True)
                try:
                    for pth in parts:
                        with open(pth, "rb") as f:
                            chunk = dict(pickle.load(f))
                        os.remove(pth)
                        self.nodes_pushed += len(chunk)
                        self.layer_peak = max(self.layer_peak,
                                              len(chunk))
                        ab = self._run_layerset(
                            chunk, h + 1, sid, seed_norm, log,
                            rss_cap, frontier_cap, depth + 1)
                        chunk = None
                        if ab:
                            return ab
                finally:
                    for pth in parts:
                        if os.path.exists(pth):
                            os.remove(pth)
                return None
            self.nodes_pushed += len(nxt)
            self.layer_peak = max(self.layer_peak, len(nxt))
            if log and (h % 6 == 0 or len(nxt) > 100000):
                print(f"[s8 u={u} nH<={self.nHmax} sid{sid} "
                      f"d{depth}] h={h + 1} states {len(nxt)} rss "
                      f"{_rss_mb()}MB {round(time.time() - t0, 1)}s",
                      flush=True)
            layer = nxt
            h += 1
        for (dyn, nH, dlt, xlo, xhi), g in layer.items():
            nd, _ = norm(dyn)
            ok = (dlt == 0) if self.exact_dlt \
                else (dlt % self.ell == 0)
            if nd == seed_norm and ok:
                self.closed.append(dict(g=g, nH=nH, dlt=dlt,
                                        extent=xhi - xlo + 1,
                                        sid=sid))
        return None

    def expand(self, g, dyn, nH, h, dlt, xlo, xhi, nxt):
        p1z, p1a, p1b, p1c, p2z, p2a, p2b, p2c = dyn
        M = 8 if (p1z | p1a | p1b | p1c | p2z | p2a | p2b
                  | p2c) & 0xFF else 0
        if M:
            p1z, p1a, p1b, p1c = (p1z << M, p1a << M, p1b << M,
                                  p1c << M)
            p2z, p2a, p2b, p2c = (p2z << M, p2a << M, p2b << M,
                                  p2c << M)
        allm = p1z | p1a | p1b | p1c | p2z | p2a | p2b | p2c
        anch_m = lsb(allm)
        u1new = m_forced(p1b, p1c, p2z, p2c)
        if not m_tooth_ok(p1a, u1new, p2z, p2a, p2b):
            return
        fixed = p1a | p1b | p1c | u1new | p2a | p2b | p2c
        wbase = wt(p1a) + wt(p1b) + wt(p1c) + wt(u1new) + wt(p2a) \
            + wt(p2b) + wt(p2c)
        wmax_any = self.whcap if nH < self.nHmax else 7
        if wbase > wmax_any:
            return
        allow = dilate(allm | u1new, self.dil)
        acols = [i for i in range(allow.bit_length())
                 if allow >> i & 1]
        classes = [(LIGHT, self.u, 7)]
        if nH < self.nHmax:
            classes.append((HEAVY, 8, self.whcap))
        for cls, wmin, wmax in classes:
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
                    if cls == LIGHT and W > 7:
                        continue
                    if cls == HEAVY and W < 8:
                        continue
                    g2 = g + W
                    rem = (self.m - h) * self.u
                    if g2 + rem > self.gcap:
                        self.trunc_g = True
                        continue
                    slab = fixed | s
                    dlt2 = dlt + (lsb(slab) - anch_m)
                    if abs(dlt2) > self.dcap:
                        self.trunc_dcap += 1
                        continue
                    state = (p1a, p1b, p1c, u1new, p2a, p2b, p2c, s)
                    a_new = lsb(slab)
                    nrows, S = norm(state)
                    xlo_c = min(xlo + M, a_new) - S
                    xhi_c = max(xhi + M,
                                slab.bit_length() - 1) - S
                    if xhi_c - xlo_c + 1 > self.extent_cap:
                        self.trunc_extent += 1
                        continue
                    nH2 = nH + (1 if cls == HEAVY else 0)
                    key = (nrows, nH2, dlt2, xlo_c, xhi_c)
                    cur = nxt.get(key)
                    if cur is not None and cur <= g2:
                        continue
                    nxt[key] = g2

    # -- parent-tracked replay (used only to reconstruct found
    # closures for independent verification) -------------------------
    def replay_seed(self, sid, w8):
        """Re-run one seed keeping all layers + parents; returns the
        reconstructed closures as absolute row lists."""
        u = self.u
        seed_norm, _ = norm(tuple(w8))
        allm = 0
        for r in w8:
            allm |= r
        key0 = (tuple(w8), 0, 0, lsb(allm), allm.bit_length() - 1)
        layers = [{key0: (u, None, None, 0)}]  # key -> (g, pkey,
        #                (u1new_c, s_c), base)
        h = 1
        while h < self.m + 1:
            nxt = {}
            for key, (g, _, _, base) in layers[-1].items():
                dyn, nH, dlt, xlo, xhi = key
                self._expand_rec(g, dyn, nH, h, dlt, xlo, xhi,
                                 base, key, nxt)
            layers.append(nxt)
            h += 1
        out = []
        for key, (g, pk, newrows, base) in layers[-1].items():
            dyn, nH, dlt, xlo, xhi = key
            nd, _ = norm(dyn)
            ok = (dlt == 0) if self.exact_dlt \
                else (dlt % self.ell == 0)
            if not (nd == seed_norm and ok):
                continue
            # walk parents back
            steps = []
            kk, hh = key, len(layers) - 1
            while hh > 0:
                g_, pk_, nr_, base_ = layers[hh][kk]
                steps.append((nr_, base_))
                kk = pk_
                hh -= 1
            steps.reverse()
            rows = [(set(), set()) for _ in range(4 + len(steps))]
            for r in range(4):
                for c in range(w8[r].bit_length()):
                    if w8[r] >> c & 1:
                        rows[r][0].add(c)
                for c in range(w8[4 + r].bit_length()):
                    if w8[4 + r] >> c & 1:
                        rows[r][1].add(c)
            for i, ((u1n, sn), base_) in enumerate(steps):
                t = 4 + i
                for c in range(u1n.bit_length()):
                    if u1n >> c & 1:
                        rows[t][0].add(c + base_)
                for c in range(sn.bit_length()):
                    if sn >> c & 1:
                        rows[t][1].add(c + base_)
            out.append(dict(g=g, nH=nH, dlt=dlt, sid=sid,
                            rows=[(frozenset(a), frozenset(b))
                                  for a, b in rows]))
        return out

    def _expand_rec(self, g, dyn, nH, h, dlt, xlo, xhi, base, pkey,
                    nxt):
        p1z, p1a, p1b, p1c, p2z, p2a, p2b, p2c = dyn
        M = 8 if (p1z | p1a | p1b | p1c | p2z | p2a | p2b
                  | p2c) & 0xFF else 0
        if M:
            p1z, p1a, p1b, p1c = (p1z << M, p1a << M, p1b << M,
                                  p1c << M)
            p2z, p2a, p2b, p2c = (p2z << M, p2a << M, p2b << M,
                                  p2c << M)
        allm = p1z | p1a | p1b | p1c | p2z | p2a | p2b | p2c
        anch_m = lsb(allm)
        u1new = m_forced(p1b, p1c, p2z, p2c)
        if not m_tooth_ok(p1a, u1new, p2z, p2a, p2b):
            return
        fixed = p1a | p1b | p1c | u1new | p2a | p2b | p2c
        wbase = wt(p1a) + wt(p1b) + wt(p1c) + wt(u1new) + wt(p2a) \
            + wt(p2b) + wt(p2c)
        wmax_any = self.whcap if nH < self.nHmax else 7
        if wbase > wmax_any:
            return
        allow = dilate(allm | u1new, self.dil)
        acols = [i for i in range(allow.bit_length())
                 if allow >> i & 1]
        classes = [(LIGHT, self.u, 7)]
        if nH < self.nHmax:
            classes.append((HEAVY, 8, self.whcap))
        for cls, wmin, wmax in classes:
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
                    if (cls == LIGHT and W > 7) or \
                            (cls == HEAVY and W < 8):
                        continue
                    g2 = g + W
                    rem = (self.m - h) * self.u
                    if g2 + rem > self.gcap:
                        continue
                    slab = fixed | s
                    dlt2 = dlt + (lsb(slab) - anch_m)
                    if abs(dlt2) > self.dcap:
                        continue
                    state = (p1a, p1b, p1c, u1new, p2a, p2b, p2c, s)
                    a_new = lsb(slab)
                    nrows, S = norm(state)
                    xlo_c = min(xlo + M, a_new) - S
                    xhi_c = max(xhi + M,
                                slab.bit_length() - 1) - S
                    if xhi_c - xlo_c + 1 > self.extent_cap:
                        continue
                    nH2 = nH + (1 if cls == HEAVY else 0)
                    key = (nrows, nH2, dlt2, xlo_c, xhi_c)
                    base2 = base - M + S
                    cur = nxt.get(key)
                    if cur is not None and cur[0] <= g2:
                        continue
                    nxt[key] = (g2, pkey, (u1new >> S, s >> S),
                                base2)


# ---------------------------------------------------------------------
# verification helpers: torus embedding + classification of the
# mirror system, and the theta' point map back to AB coordinates
# ---------------------------------------------------------------------

def mirror_code(ellc, mw, cache={}):
    """The BAbar TowerCode on (x order ellc, y order mw)."""
    key = (ellc, mw)
    if key not in cache:
        from bb_lab.tower import TowerCode
        cache[key] = TowerCode(
            f"BAbar({ellc},{mw})", (ellc, mw),
            frozenset((a % ellc, b % mw) for a, b in P_M),
            frozenset((a % ellc, b % mw) for a, b in Q_M))
    return cache[key]


def classify_rows(rows_abs, ellc, mw):
    """Embed absolute mirror rows (list over walk rows 0..len-1 of
    (u1 set, u2 set)) on the BAbar (ellc, mw) torus; the walk index
    is y mod mw, content is x mod ellc.  Returns weight/cycle/stab
    plus the theta'-mapped AB check."""
    import numpy as np
    code = mirror_code(ellc, mw)
    v = np.zeros(code.n, dtype=np.uint8)
    for t, (s1, s2) in enumerate(rows_abs):
        for c in s1:
            v[0 * code.ng + code.G.index((c % ellc, t % mw))] ^= 1
        for c in s2:
            v[1 * code.ng + code.G.index((c % ellc, t % mw))] ^= 1
    w = int(v.sum())
    isc = bool(code.is_cycle(v))
    nt = bool(isc and not code.is_stab(v))
    # theta' map back to AB-(mw, ellc): u-point (a, b) -> v-point
    # (-b mod mw, a mod ellc)
    from bb_lab.tower import TowerCode
    ab = TowerCode(f"AB({mw},{ellc})", (mw, ellc),
                   frozenset((a % mw, b % ellc) for a, b in A_L),
                   frozenset((a % mw, b % ellc) for a, b in B_L))
    va = np.zeros(ab.n, dtype=np.uint8)
    for t, (s1, s2) in enumerate(rows_abs):
        for c in s1:
            va[0 * ab.ng + ab.G.index(((-t) % mw, c % ellc))] ^= 1
        for c in s2:
            va[1 * ab.ng + ab.G.index(((-t) % mw, c % ellc))] ^= 1
    isc_ab = bool(ab.is_cycle(va))
    nt_ab = bool(isc_ab and not ab.is_stab(va))
    return dict(weight=w, is_cycle=isc, nontrivial=nt,
                ab_is_cycle=isc_ab, ab_nontrivial=nt_ab)


# ---------------------------------------------------------------------
# selftest
# ---------------------------------------------------------------------

def selftest():
    t0 = time.time()
    from bb_lab.tower import validate_banked
    validate_banked(LAB / "data")
    print("validate_banked: PASS", flush=True)
    out = {}

    print("SELFTEST 1: derivation — terms, unique forced pivot, "
          "footprint ...", flush=True)
    # forced term: unique max-row (min b) term of TERMS_U1 u TERMS_U2
    rows_read1 = sorted(-b for _, b in TERMS_U1)
    rows_read2 = sorted(-b for _, b in TERMS_U2)
    assert max(rows_read1) == 1 and max(rows_read2) == 0
    top = [(a, b) for (a, b) in TERMS_U1 if b == -1]
    assert top == [(3, -1)], top
    # footprint: cells read by E'_t at position q
    cells = [(0, -a, -b) for (a, b) in TERMS_U1] \
        + [(1, -a, -b) for (a, b) in TERMS_U2]
    xs = [c for _, c, _ in cells]
    ys = [y for _, _, y in cells]
    assert max(xs) - min(xs) == 4 and max(ys) - min(ys) == 4
    out["footprint"] = dict(x_span=4, y_span=4)
    print("  pivot x^3 u1[t+1] unique; footprint spans (4,4): PASS",
          flush=True)

    print("SELFTEST 2: bit kernel vs generic set solve (2000 "
          "random states) ...", flush=True)
    rng = random.Random(8)
    for _ in range(2000):
        rows = [frozenset(rng.sample(range(0, 14),
                                     rng.randint(0, 4)))
                for _ in range(8)]
        p = [sum(1 << c for c in r) for r in rows]
        M = 8
        got = m_forced(p[2] << M, p[3] << M, p[4] << M, p[7] << M)
        want = forced_row_sets(rows[3], rows[2], rows[7], rows[4])
        wantm = sum(1 << (c + M) for c in want)
        assert got == wantm, (rows, got, wantm)
    print("  forced-row bit kernel == generic solve: PASS",
          flush=True)

    print("SELFTEST 3: random march walks re-verified through "
          "MirrorFragment ...", flush=True)
    rng = random.Random(17)
    nwalk = 0
    for trial in range(400):
        seed = rng.choice(seeds_full(1))
        rows = [(set(), set()) for _ in range(4)]
        for r in range(4):
            for c in range(20):
                if seed[r] >> c & 1:
                    rows[r][0].add(c)
                if seed[4 + r] >> c & 1:
                    rows[r][1].add(c)
        ok = True
        for step in range(rng.randint(3, 9)):
            t = len(rows) - 1
            u1n = forced_row_sets(rows[t][0], rows[t - 1][0],
                                  rows[t][1], rows[t - 3][1])
            supp = {c for rr in rows[-4:] for c in rr[0] | rr[1]} \
                | set(u1n)
            if not supp:
                break
            loc = sorted({c + d for c in supp for d in range(-4, 5)})
            k = rng.randint(0, 3)
            new = set(rng.sample(loc, min(k, len(loc))))
            rows.append((set(u1n), new))
        fr = MirrorFragment(
            [(frozenset(a), frozenset(b)) for a, b in rows], 0)
        assert fr.admissible(), f"walk {trial} INADMISSIBLE"
        nwalk += 1
    print(f"  {nwalk} random forced walks admissible under the "
          f"generic E': PASS", flush=True)

    print("SELFTEST 4: drift additivity + gauge (on a march walk) "
          "...", flush=True)
    # build one long deterministic walk and check telescope
    rng = random.Random(5)
    seed = seeds_full(1)[3]
    rows = [(set(), set()) for _ in range(4)]
    for r in range(4):
        for c in range(20):
            if seed[r] >> c & 1:
                rows[r][0].add(c)
            if seed[4 + r] >> c & 1:
                rows[r][1].add(c)
    for step in range(12):
        t = len(rows) - 1
        u1n = forced_row_sets(rows[t][0], rows[t - 1][0],
                              rows[t][1], rows[t - 3][1])
        supp = {c for rr in rows[-4:] for c in rr[0] | rr[1]} \
            | set(u1n)
        loc = sorted({c + d for c in supp for d in range(-2, 3)})
        new = set(rng.sample(loc, 1))
        rows.append((set(u1n), new))
    fr = MirrorFragment([(frozenset(a), frozenset(b))
                         for a, b in rows], 0)
    assert fr.admissible()
    T = fr.t1
    Amid = 4 + (T - 4) // 2
    d1 = fr.subfragment(0, Amid + 3).drift()
    d2 = fr.subfragment(Amid, T).drift()
    dall = fr.drift()
    assert d1 + d2 == dall, (d1, d2, dall)
    out["additivity"] = dict(d1=d1, d2=d2, dall=dall)
    print(f"  glue telescope {d1} + {d2} = {dall}: PASS",
          flush=True)

    print("SELFTEST 5: pinch row-covering (walk-direction "
          "structure, pair-free) ...", flush=True)
    for L in range(1, 7):
        rows_b1 = set(range(-2, 2))
        rows_exit = set(range(-3, 1))
        rows_entry = set(range(L - 2, L + 2))
        assert (rows_b1 <= (rows_exit | rows_entry)) == (L <= 3)
    print("  covering holds exactly for L <= 3: PASS", flush=True)

    out["wall_s"] = round(time.time() - t0, 1)
    (DATA / "s8_xlane_selftest.json").write_text(
        json.dumps(out, indent=1))
    print(f"wrote {DATA/'s8_xlane_selftest.json'} "
          f"({out['wall_s']} s)", flush=True)


# ---------------------------------------------------------------------
# controls
# ---------------------------------------------------------------------

def control():
    t0 = time.time()
    from bb_lab.tower import validate_banked
    validate_banked(LAB / "data")
    print("validate_banked: PASS", flush=True)
    out = {}

    # ---- B: the L12' mirror species — lift + closed-march
    # detection at exact cost --------------------------------------
    print("CONTROL B: L12' (BAbar p=6 compact, w12) lift + closure "
          "detection ...", flush=True)
    rows_at, _ = atlas("BAbar", 6, 12, keep_pts=True)
    cands = [r for r in rows_at
             if r["nontrivial"] and r["weight"] == 12]
    assert cands, "no BAbar w12 nontrivial phase?!"
    done = False
    for cand in cands:
        # pts are (c, y, blk) on the BAbar cylinder (y-period 6)
        pts = [tuple(t) for t in cand["pts"]]
        ys = [y for _, y, _ in pts]
        assert min(ys) >= 0
        rowsD = {}
        for (c, y, blk) in pts:
            for k in range(0, 4):
                rowsD.setdefault(y + 6 * k, [set(), set()])[
                    blk].add(c)
        T1 = 6 * 3 + 3
        rows = [(frozenset(rowsD.get(t, [set(), set()])[0]),
                 frozenset(rowsD.get(t, [set(), set()])[1]))
                for t in range(0, T1 + 1)]
        fr = MirrorFragment(rows, 0)
        if not fr.admissible():
            continue
        anch = fr.anchors()
        per = [anch[i + 6] - anch[i] for i in range(len(anch) - 6)]
        assert all(x == 0 for x in per), per
        prof = fr.subfragment(4, 4 + 6 + 3).slabs()[:6]
        hv = [i for i, w in enumerate(prof) if w >= 8]
        # need a single-block profile for the kmax analog; count
        # heavy slabs cyclically
        nH = len(hv)
        lights = [i for i, w in enumerate(prof) if w <= 7]
        if not lights:
            continue
        j0 = min(lights, key=lambda i: prof[i])
        u0 = prof[j0]
        # slab index kk covers rows kk-3 .. kk (fr.s takes absolute
        # row indices); prof[i] is the absolute slab 7 + i of the
        # lift (subfragment(4, .) slabs start at j = 7)
        kk = 7 + j0
        w8rows = []
        for blkk in (0, 1):
            for rr in range(kk - 3, kk + 1):
                w8rows.append(sorted(fr.s(blkk, rr)))
        cols = [c for row in w8rows for c in row]
        c0 = min(cols)
        w8m = tuple(sum(1 << (c - c0) for c in row)
                    for row in w8rows)
        gsum = sum(prof)
        print(f"  L12' member: profile {prof}, heavy slabs {nH}, "
              f"seed u={u0}, g_target={gsum + u0}, extent "
              f"{fr.extent()}", flush=True)
        cm = MClosedMarch(u0, m=6, ell=cand["lstar"], nHmax=nH,
                          whcap=max(prof), smax=8, dil=8,
                          gcap=gsum + u0,
                          extent_cap=fr.extent() + 4, dcap=30,
                          exact_dlt=True)
        info = cm.run_layers([w8m], log=False)
        hits = sorted((c_["g"], c_["nH"], c_["dlt"])
                      for c_ in cm.closed)
        print(f"  march info: {info}", flush=True)
        print(f"  closed hits: {hits[:6]}", flush=True)
        ok = any(g <= gsum + u0 and d == 0 for (g, nh, d) in hits)
        assert ok, "L12' closure NOT detected"
        exact = any(g == gsum + u0 and d == 0
                    for (g, nh, d) in hits)
        # replay + independent verification of the detected closure
        reps = cm.replay_seed(0, w8m)
        ver = None
        for rp in reps:
            if rp["dlt"] != 0:
                continue
            fr2 = MirrorFragment(rp["rows"], 0)
            assert fr2.admissible()
            cl = classify_rows(rp["rows"][:6], cand["lstar"], 6)
            ver = dict(g=rp["g"], weight=cl["weight"],
                       nontrivial=cl["nontrivial"],
                       ab_nontrivial=cl["ab_nontrivial"])
            if cl["nontrivial"]:
                break
        assert ver and ver["nontrivial"], ver
        assert ver["ab_nontrivial"], \
            "theta' image must be AB-nontrivial too"
        out["L12p"] = dict(profile=prof, nH=nH, u=u0,
                           g_target=gsum + u0, exact_hit=exact,
                           verified=ver)
        print(f"  closure detected (exact-g {exact}); replay "
              f"re-verified: nontrivial on BAbar torus AND its "
              f"theta' image nontrivial on AB: PASS", flush=True)
        done = True
        break
    assert done, "no admissible single-lift L12' candidate"

    # ---- C: TC63' drift verification (the mirrored species) ------
    print("CONTROL C: TC63' (BAbar <(3,6)>, w10) — mirror cover "
          "drift +3/6 rows ...", flush=True)
    from a40_s5_twisted_atlas import (
        transform_class, tr_supp, atlas_run, verify_on_torus)
    w1, w2, gg = transform_class(6, 3)
    assert gg == 3
    trP = tr_supp(P_M, w1, w2, gg)
    trQ = tr_supp(Q_M, w1, w2, gg)
    rows_tw, _ = atlas_run(trP, trQ, gg, 11)
    ntv = [r for r in rows_tw if r["nontrivial"]]
    wmin = min(r["weight"] for r in ntv)
    assert wmin == 10, wmin
    pick = next(r for r in ntv if r["weight"] == 10)
    tv = verify_on_torus(6, 3, w1, w2, gg, pick["pts"], P_M, Q_M)
    assert tv["is_cycle"] and tv["nontrivial"] and \
        tv["weight"] == 10, tv
    # back-map to original BAbar coords, build the periodic lift
    import numpy as np
    Wm = np.array([[w1[0], w1[1]], [w2[0], w2[1]]], dtype=np.int64)
    Wi = np.array([[Wm[1, 1], -Wm[0, 1]], [-Wm[1, 0], Wm[0, 0]]],
                  dtype=np.int64)
    base_pts = []
    for (c, y, blk) in pick["pts"]:
        e0 = int(Wi[0, 0] * c + Wi[0, 1] * y)
        e1 = int(Wi[1, 0] * c + Wi[1, 1] * y)
        base_pts.append((e0, e1, blk))
    T0, T1 = 0, 6 * 3 + 4
    rows = []
    for t in range(T0, T1 + 1):
        s1, s2 = set(), set()
        for (e0, e1, blk) in base_pts:
            if (t - e1) % 6 == 0:
                k = (t - e1) // 6
                (s1 if blk == 0 else s2).add(e0 + 3 * k)
        rows.append((frozenset(s1), frozenset(s2)))
    frT = MirrorFragment(rows, T0)
    assert frT.admissible(), "TC63' lift INADMISSIBLE"
    anch = frT.anchors()
    per = [anch[i + 6] - anch[i] for i in range(len(anch) - 6)]
    assert all(x == 3 for x in per), per
    sl = frT.subfragment(T0 + 4, T0 + 4 + 6 + 3).slabs()[:6]
    out["TC63p"] = dict(drift_per_period=3, one_period_slabs=sl,
                        weight_per_period=10,
                        deficit_per_period=2 * 6 - 10,
                        window_pruned=bool(
                            frT.window_prune_events()),
                        torus_verify=tv)
    print(f"  TC63' drift +3/6 rows (gauge-free per-period "
          f"anchors), slabs {sl}, deficit 2, torus re-verified "
          f"w10 nontrivial: PASS", flush=True)

    # ---- E: the b=0 witness is NOT in the mirror sector ----------
    print("CONTROL E: a36 witness (12,12) y-gaps ...", flush=True)
    _argv = sys.argv
    sys.argv = [_argv[0], "12", "8"]
    from a40_s5_lightcore import load_witness
    sys.argv = _argv
    wit = load_witness()
    ys = sorted({p[2] % 12 for p in wit})
    occ = set(ys)
    mx = 0
    for y in range(12):
        if y in occ:
            continue
        run = 0
        yy = y
        while yy not in occ and run <= 12:
            run += 1
            yy = (yy + 1) % 12
        mx = max(mx, run)
    assert mx <= 3, (ys, mx)
    out["b0_witness"] = dict(occupied_rows=len(occ),
                             max_y_gap=mx)
    print(f"  witness occupies {len(occ)}/12 rows, max cyclic "
          f"y-gap {mx} <= 3: it is a Y-SECTOR object — the mirror "
          f"claims nothing about it (and needs no wrapped-corner "
          f"term): PASS", flush=True)

    # ---- F: mirrored stack arithmetic ----------------------------
    print("CONTROL F: TC63' closure arithmetic at the member ...",
          flush=True)
    # at AB-(24,18): mirror walk 24 rows, content circle 18.
    # A TC63' tiling needs 6 | 24 (yes, 4 periods) AND total drift
    # 4*3 = 12 == 0 mod 18 — FALSE; and any x-sector object has
    # zero drift exactly.  The species cannot close the member's
    # mirror walk.
    assert (24 % 6 == 0) and ((24 // 6) * 3) % 18 != 0
    out["member_protection"] = dict(
        periods=4, total_drift=12, content_circle=18,
        closes=False)
    print("  4 periods x drift 3 = 12 != 0 (mod 18): TC63' cannot "
          "wrap the (24,18) mirror walk (and winding is excluded "
          "in-sector): PASS", flush=True)

    out["wall_s"] = round(time.time() - t0, 1)
    (DATA / "s8_xlane_control.json").write_text(
        json.dumps(out, indent=1, default=str))
    print(f"wrote {DATA/'s8_xlane_control.json'} "
          f"({out['wall_s']} s)", flush=True)


# ---------------------------------------------------------------------
# production / coherence runs
# ---------------------------------------------------------------------

def run_prod(args):
    import gc
    gc.disable()
    t0 = time.time()
    from bb_lab.tower import validate_banked
    validate_banked(LAB / "data")
    print("validate_banked: PASS", flush=True)
    sds = [w8 for w8 in seeds_full(args.u)
           if sum(wt(r) for r in w8) == args.u]
    lo, hi = 0, len(sds)
    if args.seeds:
        lo, hi = map(int, args.seeds.split(":"))
    print(f"s8 mirror closed: u={args.u} m={args.m} ell={args.ell} "
          f"nHmax={args.nhmax} whcap={args.whcap} gcap={args.gcap} "
          f"smax={args.smax} dil={args.dil} extent<={args.extent}: "
          f"seeds [{lo}, {hi}) of {len(sds)}", flush=True)
    cm = MClosedMarch(args.u, m=args.m, ell=args.ell,
                      nHmax=args.nhmax, whcap=args.whcap,
                      smax=args.smax, dil=args.dil, gcap=args.gcap,
                      extent_cap=args.extent, dcap=args.dcap,
                      exact_dlt=not args.modell)
    info = cm.run_layers(sds[lo:hi], sid0=lo)
    print(f"info: {info}", flush=True)
    verified = []
    if cm.closed:
        print(f"  !! {len(cm.closed)} closures found — replaying "
              f"for classification", flush=True)
        for c_ in sorted(cm.closed, key=lambda d: d["g"])[:12]:
            sid = c_["sid"]
            reps = cm.replay_seed(sid, sds[sid])
            for rp in reps:
                fr = MirrorFragment(rp["rows"], 0)
                assert fr.admissible()
                cl = classify_rows(rp["rows"][:args.m], args.ell,
                                   args.m)
                verified.append(dict(g=rp["g"], nH=rp["nH"],
                                     dlt=rp["dlt"], sid=sid, **cl))
                print(f"  closure g={rp['g']} nH={rp['nH']} "
                      f"dlt={rp['dlt']}: weight={cl['weight']} "
                      f"cycle={cl['is_cycle']} "
                      f"nontrivial={cl['nontrivial']} "
                      f"(AB: {cl['ab_nontrivial']})", flush=True)
    out = dict(
        params=dict(u=args.u, m=args.m, ell=args.ell,
                    nHmax=args.nhmax, whcap=args.whcap,
                    gcap=args.gcap, smax=args.smax, dil=args.dil,
                    extent_cap=args.extent, dcap=args.dcap,
                    exact_dlt=not args.modell, seeds=[lo, hi]),
        info=info,
        closed=sorted(cm.closed, key=lambda d: d["g"]),
        verified=verified)
    p = DATA / (f"s8_{args.tag or 'prod'}_u{args.u}"
                f"_nh{args.nhmax}_g{args.gcap}_s{lo}_{hi}.json")
    p.write_text(json.dumps(out, indent=1))
    print(f"wrote {p} ({round(time.time() - t0, 1)} s)", flush=True)


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("selftest")
    sub.add_parser("control")
    s = sub.add_parser("prod")
    s.add_argument("--u", type=int, default=1)
    s.add_argument("--m", type=int, default=24)
    s.add_argument("--ell", type=int, default=18)
    s.add_argument("--nhmax", type=int, default=2)
    s.add_argument("--whcap", type=int, default=19)
    s.add_argument("--gcap", type=int, required=True)
    s.add_argument("--smax", type=int, default=3)
    s.add_argument("--dil", type=int, default=4)
    s.add_argument("--extent", type=int, default=14)
    s.add_argument("--dcap", type=int, default=20)
    s.add_argument("--seeds", type=str, default="")
    s.add_argument("--modell", action="store_true",
                   help="closure mod ell (control mode)")
    s.add_argument("--tag", type=str, default="")
    args = ap.parse_args()
    if args.cmd == "selftest":
        selftest()
    elif args.cmd == "control":
        control()
    else:
        run_prod(args)


if __name__ == "__main__":
    main()
