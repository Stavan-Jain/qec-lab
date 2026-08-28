#!/usr/bin/env python3
"""A40 S6 — Stage 1b: the deficit-vs-drift frontier of light pruned
walk fragments, by min-window-stratified cost-ordered BnB (cover
mode; NEVER graph materialization).

Object: admissible fragments of the y-walk system on the universal
x-cover — rows v1[t], v2[t] (integer columns), E_j enforced whenever
fully supported, every 4-row slab weight in [1, 7] (light,
y-spanning), the H = 5 window-prune rule (no fully visible tooth held
> half), (4,4)-connected new growth.  Per fragment: h = #slabs,
g = sum of slab weights (the exact telescope cost: 4|v| on closed
walks), drift delta = A_last - A_first where A_j = min occupied
column of rows [j-3, j].  The frontier: minG(h, delta), equivalently
D(h, delta) = 2h - minG/4 (deficit; per light slab the deficit gains
2 - W/4 >= 1/4, so the frontier is a min-cost problem).

STRATIFICATION (the pricing move that beats the seed wall).  Split
every fragment at its FIRST minimum-weight window u:
- u >= 5: every slab >= 5 gives D <= (2 - u/4) h <= 0.75 h < (6/7) h
  analytically — below the W7 rate; no enumeration needed.
- u <= 4: the split window has weight u <= 4; full-content connected
  u <= 4 windows number ~4e5 (vs ~2e8 for u <= 7 — the 4-row
  materialization wall, avoided).  The piece above the split is a
  BACKWARD fragment (v1-forced march down), the piece below a FORWARD
  fragment (v2-forced march up), both in stratum u (all slabs >= u),
  glued on the shared window: g = gB + gF - u, h = hB + hF - 1,
  delta = deltaB + deltaF (anchor telescope, exact).
The tables here use the LOOSE join (min over seeds independently on
each side) — it can only OVERSTATE deficits, hence is sound for the
assembly; the species runs measure its tightness.

Search: Dijkstra on g with dominance memo[(state, h, delta)] -> min
g; states = the dynamic window (forward: v1[t-3..t], v2[t-2..t];
backward: v1[t..t+2], v2[t..t+4]) normalized by x-translation.
Free-input branching capped at SMAX new points per row within the
+-DILATE dilation of the live support (scope caps, stability-checked
by re-runs at larger caps).  Forced rows alone decide the tooth
check, so pruning happens before input branching.
"""
from __future__ import annotations

import json
import sys
import time
from heapq import heappush, heappop
from itertools import combinations
from pathlib import Path

LAB = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(LAB / "src"))
sys.path.insert(0, str(LAB / "scripts"))

DATA = LAB / "data" / "a40"

# rows are bitmasks; bit c = column c (columns shifted >= 0 locally).
# support shift by +s = mask << s.  Forced rows:
#   forward:  v2[t+1] = v1[t] ^ (v1[t]>>1) ^ (v1[t-3]<<1) ^ v2[t]
#             ^ (v2[t-1]>>3)
#   backward: v1[t-1] = ( (v1[t+2] ^ (v1[t+2]>>1) ^ v2[t+2] ^ v2[t+3]
#             ^ (v2[t+1]>>3) ) >> 1     [E_{t+2} solved for x*v1[t-1]]


def wt(x):
    return bin(x).count("1")


def lsb(x):
    return (x & -x).bit_length() - 1


def span(x):
    return 0 if x == 0 else x.bit_length() - lsb(x)


def norm(rows):
    """Shift a tuple of masks so the least set bit overall is 0;
    returns (shifted rows, shift applied)."""
    m = 0
    for r in rows:
        m |= r
    if m == 0:
        return rows, 0
    s = lsb(m)
    if s == 0:
        return rows, 0
    return tuple(r >> s for r in rows), s


def dilate(mask, r):
    out = mask
    for _ in range(r):
        out |= (out << 1) | (out >> 1)
    return out


def tooth_ok(v1a, v1b, v1c, v2a, v2new):
    """Tooth check at cy: blk1 rows cy-1,cy,cy+1 = (v1a, v1b, v1c);
    blk2 row cy = v2a (cells cx, cx-1), row cy+3 = v2new (cell cx+1).
    Overlap at alignment cx =
      v1a[cx]? wait — cells: (0,cx,cy)=v1b, (0,cx,cy-1)=v1a,
      (0,cx-3,cy+1)=v1c, (1,cx,cy)=v2a, (1,cx-1,cy)=v2a,
      (1,cx+1,cy+3)=v2new.
    Require overlap <= 3 for every cx."""
    m = v1a | v1b | v1c | v2a | v2new
    if m == 0:
        return True
    lo, hi = lsb(m) - 4, m.bit_length() + 4
    for cx in range(lo, hi + 1):
        ov = ((v1b >> cx) & 1) + ((v1a >> cx) & 1) \
            + ((v1c >> (cx - 3)) & 1 if cx >= 3 else 0) \
            + ((v2a >> cx) & 1) \
            + ((v2a >> (cx - 1)) & 1 if cx >= 1 else 0) \
            + ((v2new >> (cx + 1)) & 1)
        if ov > 3:
            return False
    return True


def seeds_full(u, max_span=25):
    """All connected full-content 4-row windows of weight exactly u:
    columns with consecutive gaps <= 4 (the (4,4) adjacency), 8 cells
    per column (v1 rows 0..3 = bits 0..3, v2 rows 0..3 = bits 4..7),
    min column = 0.  Yields (v1a,v1b,v1c,v1d, v2a,v2b,v2c,v2d)."""
    out = []

    def build(cols_masks, rem, last_col):
        if rem == 0:
            rows = [0] * 8
            for c, cm in cols_masks:
                for b in range(8):
                    if cm >> b & 1:
                        rows[b] |= 1 << c
            out.append(tuple(rows))
            return
        # extend current column set with a new column
        for c in range(last_col + 1, last_col + 5):
            if c > max_span:
                continue
            for w in range(1, rem + 1):
                for cells in combinations(range(8), w):
                    cm = 0
                    for b in cells:
                        cm |= 1 << b
                    build(cols_masks + [(c, cm)], rem - w, c)

    # first column at 0
    for w in range(1, u + 1):
        for cells in combinations(range(8), w):
            cm = 0
            for b in cells:
                cm |= 1 << b
            build([(0, cm)], u - w, 0)
    return out


class March:
    """Cost-ordered stratified march (forward or backward).

    Node = (g, tie, dyn, h, dlt, anch, xlo, xhi): dyn = normalized
    state rows (masks in the node's LOCAL frame, min set bit = 0
    across the state); anch = anchor (min col) of the node's last
    completed slab, xlo/xhi = fragment extent so far — all three in
    the node's local frame (they may be negative: dropped rows can
    hold the old minimum).  Frame passing: expansion may pre-shift
    everything by a margin M (so >> terms cannot underflow), then
    the child normalization shift S re-bases; anch/xlo/xhi are
    translated with the same M and S, so anchor DIFFERENCES (drift
    increments) are frame-exact."""

    def __init__(self, u, direction="fwd", smax=3, dil=4, hcap=18,
                 gcap=60, extent_cap=34, dcap=30):
        self.u, self.dir = u, direction
        self.smax, self.dil = smax, dil
        self.hcap, self.gcap = hcap, gcap
        self.extent_cap, self.dcap = extent_cap, dcap
        self.table = {}          # (h, delta) -> min g  (pure)
        self.table_pen = {}      # (h, delta) -> min g + boundary pen
        self.nodes_popped = 0
        self.nodes_pushed = 0
        self.trunc_g = False       # benign: defines the g-cap scope
        self.trunc_extent = 0      # scope caps hit (report!)
        self.trunc_dcap = 0
        self.trunc_nodes = False

    def debt(self, rows):
        """Weight of the E-forced extension row at this march's OUTER
        end (fwd: the forced v2 of the next row up; bwd: the forced
        v1 of the next row down).  A maximal run's flanking slab is
        heavy AND must hold this row: W_flank >= max(8, debt).  The
        penalty max(0, debt-8)//2 (halved: the flanking slab may be
        shared with the neighbouring run's debt; avg <= max keeps the
        halved sum sound) charges ghost boundary relief against heavy
        excess.  Pure function of the endpoint state, so pure-g
        dominance remains sound for the penalized table."""
        if self.dir == "fwd":
            v1a, v1b, v1c, v1d, v2a, v2b, v2c = rows
            # E_t forces v2[t+1] = (1+x^-1)v1[t] + x v1[t-3] + v2[t]
            # + x^-3 v2[t-1]; margin-shift to protect >>:
            M = 8 if (v1a | v1d | v2b | v2c) & 0xFF else 0
            a, d_, b, c = v1a << M, v1d << M, v2b << M, v2c << M
            nxt = d_ ^ (d_ >> 1) ^ (a << 1) ^ c ^ (b >> 3)
            return wt(nxt)
        v1t, v1u_, v1v, v2t, v2u_, v2v, v2w, v2x = rows
        # E_{t+2} forces x v1[t-1] = (1+x^-1)v1[t+2] + v2[t+2]
        # + v2[t+3] + x^-3 v2[t+1]
        M = 8 if (v1v | v2u_ | v2v | v2w) & 0xFF else 0
        vv, u2, v2_, w2 = v1v << M, v2u_ << M, v2v << M, v2w << M
        rhs = vv ^ (vv >> 1) ^ v2_ ^ w2 ^ (u2 >> 3)
        return wt(rhs)

    def record(self, h, dlt, g, pen):
        k = (h, dlt)
        if k not in self.table or self.table[k] > g:
            self.table[k] = g
        gp = g + pen
        if k not in self.table_pen or self.table_pen[k] > gp:
            self.table_pen[k] = gp

    def run(self, seed_states, log=None, log_every=200000,
            max_nodes=30_000_000):
        u = self.u
        heap = []
        seen = {}
        n_seed = 0
        for w8 in seed_states:
            if sum(wt(r) for r in w8) != u:
                continue
            n_seed += 1
            if self.dir == "fwd":
                # state after row t=3: (v1[0..3], v2[1..3])
                dyn = (w8[0], w8[1], w8[2], w8[3], w8[5], w8[6],
                       w8[7])
            else:
                # state-down at t = w-3: (v1[w-3..w-1],
                # v2[w-3..w+1]); v2[w+1] = 0 (outside the piece)
                dyn = (w8[0], w8[1], w8[2], w8[4], w8[5], w8[6],
                       w8[7], 0)
            allm = 0
            for r in w8:
                allm |= r
            anchor = lsb(allm)     # 0 for normalized seeds
            key = (dyn, 1, 0)
            if key in seen and seen[key] <= u:
                continue
            seen[key] = u
            heappush(heap, (u, self.nodes_pushed, dyn, 1, 0, anchor,
                            anchor, allm.bit_length() - 1))
            self.nodes_pushed += 1
            self.record(1, 0, u, max(0, self.debt(dyn) - 8) // 2)
        t0 = time.time()
        while heap:
            g, _, dyn, h, dlt, anch, xlo, xhi = heappop(heap)
            key = (dyn, h, dlt)
            if seen.get(key, 1 << 30) < g:
                continue
            self.nodes_popped += 1
            if log and self.nodes_popped % log_every == 0:
                print(f"[{self.dir} u={u}] popped {self.nodes_popped}"
                      f" pushed {self.nodes_pushed} heap {len(heap)}"
                      f" g={g} {round(time.time()-t0, 1)}s",
                      flush=True)
            if self.nodes_popped > max_nodes:
                self.trunc_nodes = True
                break
            if h >= self.hcap:
                continue
            if self.dir == "fwd":
                self.expand_fwd(g, dyn, h, dlt, anch, xlo, xhi, heap,
                                seen)
            else:
                self.expand_bwd(g, dyn, h, dlt, anch, xlo, xhi, heap,
                                seen)
        return dict(n_seeds=n_seed, popped=self.nodes_popped,
                    pushed=self.nodes_pushed,
                    trunc_g=self.trunc_g,
                    trunc_extent=self.trunc_extent,
                    trunc_dcap=self.trunc_dcap,
                    trunc_nodes=self.trunc_nodes,
                    wall_s=round(time.time() - t0, 1))

    def push(self, g2, state_rows, slab_mask, slab_hi, h, dlt_new,
             M, anch, xlo, xhi, heap, seen):
        """Common push: state_rows and slab data in the margin frame
        (parent local + M); anch/xlo/xhi still in parent local."""
        a_new = lsb(slab_mask)                 # margin frame
        nrows, S = norm(state_rows)
        # child frame = margin frame - S
        anch_c = a_new - S
        xlo_c = min(xlo + M, a_new) - S
        xhi_c = max(xhi + M, slab_hi) - S
        if xhi_c - xlo_c + 1 > self.extent_cap:
            self.trunc_extent += 1
            return
        h2 = h + 1
        key = (nrows, h2, dlt_new)
        if seen.get(key, 1 << 30) <= g2:
            return
        seen[key] = g2
        self.record(h2, dlt_new, g2,
                    max(0, self.debt(nrows) - 8) // 2)
        heappush(heap, (g2, self.nodes_pushed, nrows, h2, dlt_new,
                        anch_c, xlo_c, xhi_c))
        self.nodes_pushed += 1

    # -- forward ------------------------------------------------------
    def expand_fwd(self, g, dyn, h, dlt, anch, xlo, xhi, heap, seen):
        v1a, v1b, v1c, v1d, v2a, v2b, v2c = dyn
        M = 8 if (v1a | v1b | v1c | v1d | v2a | v2b | v2c) & 0xFF \
            else 0
        if M:
            v1a, v1b, v1c, v1d = (v1a << M, v1b << M, v1c << M,
                                  v1d << M)
            v2a, v2b, v2c = v2a << M, v2b << M, v2c << M
        anch_m = anch + M                     # margin frame
        # forced v2[t+1] from E_t
        v2new = v1d ^ (v1d >> 1) ^ (v1a << 1) ^ v2c ^ (v2b >> 3)
        # tooth at cy = t-2: (v1[t-3], v1[t-2], v1[t-1], v2[t-2],
        # v2new)
        if not tooth_ok(v1a, v1b, v1c, v2a, v2new):
            return
        # new slab j = t+1 covers rows t-2..t+1 BOTH blocks:
        # v1: v1b, v1c, v1d, s(input); v2: v2a, v2b, v2c, v2new.
        fixed = v1b | v1c | v1d | v2a | v2b | v2c | v2new
        wbase = wt(v1b) + wt(v1c) + wt(v1d) + wt(v2a) + wt(v2b) \
            + wt(v2c) + wt(v2new)
        if wbase > 7:
            return
        allow = dilate(v1a | v1b | v1c | v1d | v2a | v2b | v2c
                       | v2new, self.dil)
        acols = [i for i in range(allow.bit_length())
                 if allow >> i & 1]
        need = max(0, self.u - wbase)
        room = 7 - wbase
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
                self.push(g2, (v1b, v1c, v1d, s, v2b, v2c, v2new),
                          slab, slab.bit_length() - 1, h, dlt2, M,
                          anch, xlo, xhi, heap, seen)

    # -- backward -----------------------------------------------------
    def expand_bwd(self, g, dyn, h, dlt, anch, xlo, xhi, heap, seen):
        v1t, v1u_, v1v, v2t, v2u_, v2v, v2w, v2x = dyn
        M = 8 if (v1t | v1u_ | v1v | v2t | v2u_ | v2v | v2w
                  | v2x) & 0xFF else 0
        if M:
            v1t, v1u_, v1v = v1t << M, v1u_ << M, v1v << M
            v2t, v2u_, v2v, v2w, v2x = (v2t << M, v2u_ << M,
                                        v2v << M, v2w << M,
                                        v2x << M)
        anch_m = anch + M
        # forced v1[t-1] via E_{t+2}: x v1[t-1] = (1+x^-1)v1[t+2]
        # + v2[t+2] + v2[t+3] + x^-3 v2[t+1]
        rhs = v1v ^ (v1v >> 1) ^ v2v ^ v2w ^ (v2u_ >> 3)
        if rhs & 1:
            return
        v1new = rhs >> 1
        # tooth at cy = t: blk1 rows t-1..t+1 = (v1new, v1t, v1u_);
        # blk2 rows t (v2t) and t+3 (v2w)
        if not tooth_ok(v1new, v1t, v1u_, v2t, v2w):
            return
        # new slab j = t+2 covers rows t-1..t+2:
        # v1: v1new, v1t, v1u_, v1v; v2: s(input), v2t, v2u_, v2v.
        fixed = v1new | v1t | v1u_ | v1v | v2t | v2u_ | v2v
        wbase = wt(v1new) + wt(v1t) + wt(v1u_) + wt(v1v) + wt(v2t) \
            + wt(v2u_) + wt(v2v)
        if wbase > 7:
            return
        allow = dilate(v1new | v1t | v1u_ | v1v | v2t | v2u_ | v2v
                       | v2w, self.dil)
        acols = [i for i in range(allow.bit_length())
                 if allow >> i & 1]
        need = max(0, self.u - wbase)
        room = 7 - wbase
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
                # marching down: dlt tracks A_seed - A_bottom
                dlt2 = dlt - (lsb(slab) - anch_m)
                if abs(dlt2) > self.dcap:
                    self.trunc_dcap += 1
                    continue
                self.push(g2, (v1new, v1t, v1u_, s, v2t, v2u_, v2v,
                               v2w), slab, slab.bit_length() - 1, h,
                          dlt2, M, anch, xlo, xhi, heap, seen)


def compose(tabF, tabB, u, hmax=18):
    """Loose min-window join: minG(h, delta) <= min over splits of
    gB + gF - u.  Returns {(h, delta): minG}."""
    out = {}
    for (hf, df), gf in tabF.items():
        for (hb, db), gb in tabB.items():
            h = hf + hb - 1
            if h > hmax:
                continue
            d = df + db
            g = gf + gb - u
            k = (h, d)
            if k not in out or out[k] > g:
                out[k] = g
    return out


def species_subfragments(hmax=18):
    """(h, delta, g, min_slab) for every contiguous slab-window of
    the two light species' long runs (W7 at l=24 lifted long; the
    l=12 W7 twin) — the validation set the tables must dominate."""
    import numpy as np  # noqa: F401
    sys.argv = [sys.argv[0]]
    from a40_s6_drift import lift_phase, load_survivor
    from a40_s5_lightcore import Phase
    out = []
    pts = load_survivor("s5_dense_l24p7.json", 7, 22, 8)
    ph = Phase.from_quotient_pts(24, 7, 22, pts)
    fr, s = lift_phase(ph, n_periods=4)
    assert s == -2
    slabs = fr.slabs()
    anchors = fr.anchors()
    n = len(slabs)
    for i in range(n):
        for j in range(i, min(n, i + hmax)):
            h = j - i + 1
            g = sum(slabs[i:j + 1])
            d = anchors[j] - anchors[i]
            out.append((h, d, g, min(slabs[i:j + 1])))
    return out


def main():
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--u", type=int, required=True)
    ap.add_argument("--gcap", type=int, default=40)
    ap.add_argument("--hcap", type=int, default=18)
    ap.add_argument("--smax", type=int, default=3)
    ap.add_argument("--dil", type=int, default=4)
    ap.add_argument("--maxnodes", type=int, default=30_000_000)
    ap.add_argument("--tag", type=str, default="")
    args = ap.parse_args()

    t0 = time.time()
    u = args.u
    print(f"stratum u={u}: generating full-content seeds...",
          flush=True)
    sds = seeds_full(u)
    print(f"  {len(sds)} seeds ({round(time.time()-t0, 1)} s)",
          flush=True)
    res = {}
    tabs, tabs_pen = {}, {}
    for direction in ("fwd", "bwd"):
        m = March(u, direction, smax=args.smax, dil=args.dil,
                  hcap=args.hcap, gcap=args.gcap)
        info = m.run(sds, log=True, max_nodes=args.maxnodes)
        tabs[direction] = dict(m.table)
        tabs_pen[direction] = dict(m.table_pen)
        res[direction] = dict(info=info,
                              table={f"{h},{d}": g for (h, d), g
                                     in sorted(m.table.items())},
                              table_pen={f"{h},{d}": g for (h, d), g
                                         in sorted(
                                             m.table_pen.items())})
        print(f"{direction} u={u}: {info} — {len(m.table)} (h,delta) "
              f"buckets", flush=True)
        res[direction]["params"] = vars(args)
    comp = compose(tabs["fwd"], tabs["bwd"], u, args.hcap)
    comp_pen = compose(tabs_pen["fwd"], tabs_pen["bwd"], u, args.hcap)
    res["composed"] = {f"{h},{d}": g for (h, d), g
                       in sorted(comp.items())}
    res["composed_pen"] = {f"{h},{d}": g for (h, d), g
                           in sorted(comp_pen.items())}
    # deficit view
    best = {}
    for (h, d), g in comp.items():
        D = 2 * h - g / 4
        if D > best.get(h, (-99, None))[0]:
            best[h] = (D, d)
    best_pen = {}
    for (h, d), g in comp_pen.items():
        D = 2 * h - g / 4
        if D > best_pen.get(h, (-99, None))[0]:
            best_pen[h] = (D, d)
    res["best_pen_deficit_per_h"] = {
        h: dict(D=round(v[0], 3), delta=v[1])
        for h, v in sorted(best_pen.items())}
    res["best_deficit_per_h"] = {h: dict(D=round(v[0], 3), delta=v[1])
                                 for h, v in sorted(best.items())}
    print("best deficit per h (composed pure | penalized):",
          flush=True)
    for h in sorted(best):
        v, vp = best[h], best_pen.get(h, (None, None))
        print(f"  h={h}: D={v[0]:.2f}@delta={v[1]}  |  "
              f"Dpen={vp[0]:.2f}@delta={vp[1]}" if vp[0] is not None
              else f"  h={h}: D={v[0]:.2f}@delta={v[1]}", flush=True)
    # validation: the composed table must dominate every species
    # sub-fragment whose min slab lands in this stratum and whose
    # cost is within the g-cap.
    n_checked = n_viol = 0
    for (h, d, g, mins) in species_subfragments(args.hcap):
        if mins != u or g > args.gcap or h > args.hcap:
            continue
        n_checked += 1
        got = comp.get((h, d))
        if got is None or got > g:
            n_viol += 1
            print(f"  VIOLATION: species sub-fragment (h={h}, "
                  f"delta={d}, g={g}) not dominated (table: {got})",
                  flush=True)
    res["species_validation"] = dict(checked=n_checked,
                                    violations=n_viol)
    print(f"species validation (stratum {u}): {n_checked} "
          f"sub-fragments checked, {n_viol} violations", flush=True)
    tag = args.tag or f"u{u}_g{args.gcap}"
    out_p = DATA / f"s6_frontier_{tag}.json"
    out_p.write_text(json.dumps(res, indent=1))
    print(f"wrote {out_p} ({round(time.time()-t0, 1)} s)", flush=True)


if __name__ == "__main__":
    main()
