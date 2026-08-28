#!/usr/bin/env python3
"""A40 S6 — ghost-passage verification: (1) MATCHED (per-seed)
min-window join at u = 1 — are the deep cheap composed buckets real
fragments or loose-join chimeras?  (2) end-to-end reconstruction of
one deep backward ghost coast and its independent verification
through the CoverFragment machinery (E_j, slabs, window rule,
anchors, drift).
"""
from __future__ import annotations

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
    March, seeds_full, wt, lsb, norm, tooth_ok, dilate,
)


class MarchSeeded(March):
    """March with per-seed attribution: dominance and tables keyed
    by seed id.  Also keeps parent pointers for path readout."""

    def __init__(self, *a, keep_parents=False, **kw):
        super().__init__(*a, **kw)
        self.tables = {}          # sid -> {(h, delta): g}
        self.parents = {} if keep_parents else None
        self.keep_parents = keep_parents

    def record_sid(self, sid, h, dlt, g):
        t = self.tables.setdefault(sid, {})
        k = (h, dlt)
        if k not in t or t[k] > g:
            t[k] = g

    def run_seeded(self, seed_states, log_every=2_000_000,
                   max_nodes=60_000_000):
        u = self.u
        heap = []
        seen = {}
        sids = []
        for sid, w8 in enumerate(seed_states):
            if sum(wt(r) for r in w8) != u:
                continue
            sids.append(sid)
            if self.dir == "fwd":
                dyn = (w8[0], w8[1], w8[2], w8[3], w8[5], w8[6],
                       w8[7])
            else:
                dyn = (w8[0], w8[1], w8[2], w8[4], w8[5], w8[6],
                       w8[7], 0)
            allm = 0
            for r in w8:
                allm |= r
            anchor = lsb(allm)
            key = (sid, dyn, 1, 0)
            seen[key] = u
            heappush(heap, (u, self.nodes_pushed,
                            (sid, dyn, 1, 0, anchor, anchor,
                             allm.bit_length() - 1)))
            self.nodes_pushed += 1
            self.record_sid(sid, 1, 0, u)
            if self.keep_parents:
                self.parents[key] = None
        t0 = time.time()
        while heap:
            g, _, (sid, dyn, h, dlt, anch, xlo, xhi) = heappop(heap)
            key = (sid, dyn, h, dlt)
            if seen.get(key, 1 << 30) < g:
                continue
            self.nodes_popped += 1
            if self.nodes_popped % log_every == 0:
                print(f"[{self.dir} u={self.u} seeded] popped "
                      f"{self.nodes_popped} heap {len(heap)} g={g} "
                      f"{round(time.time()-t0, 1)}s", flush=True)
            if self.nodes_popped > max_nodes:
                self.trunc_nodes = True
                break
            if h >= self.hcap:
                continue
            self._sid = sid
            self._pkey = key
            if self.dir == "fwd":
                self.expand_fwd(g, dyn, h, dlt, anch, xlo, xhi, heap,
                                seen)
            else:
                self.expand_bwd(g, dyn, h, dlt, anch, xlo, xhi, heap,
                                seen)
        return dict(popped=self.nodes_popped,
                    pushed=self.nodes_pushed,
                    trunc_nodes=self.trunc_nodes,
                    wall_s=round(time.time() - t0, 1))

    # override push to carry sid + parents
    def push(self, g2, state_rows, slab_mask, slab_hi, h, dlt_new,
             M, anch, xlo, xhi, heap, seen):
        a_new = lsb(slab_mask)
        nrows, S = norm(state_rows)
        anch_c = a_new - S
        xlo_c = min(xlo + M, a_new) - S
        xhi_c = max(xhi + M, slab_hi) - S
        if xhi_c - xlo_c + 1 > self.extent_cap:
            self.trunc_extent += 1
            return
        h2 = h + 1
        sid = self._sid
        key = (sid, nrows, h2, dlt_new)
        if seen.get(key, 1 << 30) <= g2:
            return
        seen[key] = g2
        self.record_sid(sid, h2, dlt_new, g2)
        if self.keep_parents:
            # store the un-normalized new state + norm shift so the
            # walk can be replayed: parent key + the new last row
            # pair in the CHILD frame
            self.parents[key] = (self._pkey, S, M)
        heappush(heap, (g2, self.nodes_pushed,
                        (sid, nrows, h2, dlt_new, anch_c, xlo_c,
                         xhi_c)))
        self.nodes_pushed += 1


def matched_join(tabF, tabB, u, hmax=19):
    out = {}
    for sid in tabF:
        if sid not in tabB:
            continue
        for (hf, df), gf in tabF[sid].items():
            for (hb, db), gb in tabB[sid].items():
                h = hf + hb - 1
                if h > hmax:
                    continue
                k = (h, df + db)
                g = gf + gb - u
                if k not in out or out[k] > g:
                    out[k] = g
    return out


def main():
    t0 = time.time()
    gcap = int(sys.argv[1]) if len(sys.argv) > 1 else 24
    sds = seeds_full(1)
    print(f"u=1 matched join at gcap {gcap}: {len(sds)} seeds",
          flush=True)
    tabs = {}
    for direction in ("fwd", "bwd"):
        m = MarchSeeded(1, direction, hcap=19, gcap=gcap)
        info = m.run_seeded(sds)
        tabs[direction] = m.tables
        print(f"{direction}: {info}; {sum(len(t) for t in m.tables.values())} "
              f"seeded entries", flush=True)
    comp = matched_join(tabs["fwd"], tabs["bwd"], 1)
    best = {}
    for (h, d), g in comp.items():
        D = 2 * h - g / 4
        if D > best.get(h, (-99, None, None))[0]:
            best[h] = (D, d, g)
    print("MATCHED best-D per h (u=1):", flush=True)
    for h, v in sorted(best.items()):
        print(f"  h={h}: D={v[0]:.2f} delta={v[1]} g={v[2]}",
              flush=True)
    # compare against the loose join from the prod file
    loose = json.loads((DATA / "s6_frontier_u1prod.json").read_text()
                       )["composed"]
    loose = {tuple(map(int, k.split(","))): g for k, g in loose.items()}
    diffs = []
    for (h, d), g in comp.items():
        lg = loose.get((h, d))
        if lg is not None and lg < g:
            diffs.append((h, d, lg, g))
    print(f"loose < matched on {len(diffs)} buckets (loose join "
          f"looseness); largest gaps:",
          sorted(diffs, key=lambda t: t[3] - t[2])[-5:], flush=True)

    # ---- explicit deep-ghost reconstruction + verification --------
    m = MarchSeeded(1, "bwd", hcap=12, gcap=20, keep_parents=True)
    info = m.run_seeded(sds)
    print(f"parented bwd run: {info}", flush=True)
    # pick max-h min-g among parent keys
    cand = {}
    for key in m.parents:
        sid, dyn, h, dlt = key
        g = m.tables[sid].get((h, dlt), 10**9)
        if h not in cand or g < cand[h][1]:
            cand[h] = (key, g)
    hstar = max(h for h in cand if cand[h][1] < 10**9)
    key, gstar = cand[hstar]
    print(f"replaying bwd ghost: h={hstar}, g={gstar}, "
          f"key-delta={key[3]}", flush=True)
    # replay: walk parents up to the seed, collecting norm shifts;
    # then re-simulate the backward march applying the recorded
    # free inputs is complex — instead reconstruct rows directly
    # from the chain of states: each backward step's new bottom row
    # is (v1new = dyn[0], v2new = dyn[3]) in that node's frame.
    chain = []
    k = key
    while k is not None:
        chain.append(k)
        pk = m.parents[k]
        if pk is None:
            k = None
        else:
            k = pk[0]
    chain.reverse()          # seed ... deepest
    # frames: child frame = margin(parent) - S; accumulate offset of
    # each node's frame relative to the seed frame:
    offs = [0]
    for i in range(1, len(chain)):
        pk = m.parents[chain[i]]
        S, M = pk[1], pk[2]
        offs.append(offs[-1] + M - S)
    # wait: child cols = parent cols + M - S  => to express child
    # masks in seed frame: add offs (col_seed = col_child + off?)
    # col_child = col_parent + M - S  => col_parent = col_child -
    # (M - S): seed frame = parent-most: col_seed = col_child -
    # sum(M - S) along the chain: offs[i] holds the sum; so
    # col_seed = col_child - offs[i]... sign check done empirically
    # below via E-verification.
    sys.argv = [sys.argv[0]]
    from a40_s6_drift import CoverFragment
    # rows: the seed contributes rows [w-3, w] (frame 0); each step i
    # adds one row BELOW: row w-3-i with v1 = dyn[0], v2 = dyn[3]:
    sid0, dyn0, _, _ = chain[0]
    w8 = sds[sid0]
    rows = {}                # global row index -> (set1, set2)
    for r in range(4):
        rows[r] = (frozenset(i for i in range(64) if w8[r] >> i & 1),
                   frozenset(i for i in range(64)
                             if w8[4 + r] >> i & 1))
    for i in range(1, len(chain)):
        sid, dyn, h, dlt = chain[i]
        off = offs[i]
        v1new, v2new = dyn[0], dyn[3]
        rows[-i] = (
            frozenset((c - off) for c in range(v1new.bit_length())
                      if v1new >> c & 1),
            frozenset((c - off) for c in range(v2new.bit_length())
                      if v2new >> c & 1))
    lo = min(rows)
    fr = CoverFragment([rows[j] for j in range(lo, 4)], lo)
    ok = fr.admissible()
    sl = fr.slabs()
    print(f"reconstructed ghost fragment rows [{lo}, 3]: "
          f"admissible={ok}, weight={fr.weight()}, slabs={sl}, "
          f"drift={fr.drift() if sl else None}, window_pruned="
          f"{bool(fr.window_prune_events())}", flush=True)
    assert ok, "ghost fragment fails independent E verification!"
    assert max(sl) <= 7 and min(sl) >= 1
    out = dict(
        matched_best={f"{h}": dict(D=round(v[0], 3), delta=v[1],
                                   g=v[2])
                      for h, v in sorted(best.items())},
        loose_vs_matched_diffs=len(diffs),
        ghost=dict(h=hstar, g=gstar, admissible=bool(ok),
                   weight=fr.weight(), slabs=sl,
                   drift=fr.drift(),
                   window_pruned=bool(fr.window_prune_events())),
        wall_s=round(time.time() - t0, 1))
    (DATA / "s6_ghost_verify.json").write_text(
        json.dumps(out, indent=1))
    print(f"wrote {DATA/'s6_ghost_verify.json'} "
          f"({out['wall_s']} s)", flush=True)
    # also dump the matched composed table in the frontier schema so
    # the assembly can consume it as the (tighter, still sound)
    # u = 1 table: matched join = real fragments only.
    front = dict(
        fwd=dict(info=dict(note="seeded"), params=dict(gcap=gcap),
                 table={}),
        bwd=dict(info=dict(note="seeded"), params=dict(gcap=gcap),
                 table={}),
        composed={f"{h},{d}": g for (h, d), g in sorted(comp.items())},
        species_validation=dict(checked=0, violations=0,
                                note="u=1: no species sub-fragment "
                                     "has min slab 1"))
    (DATA / "s6_frontier_u1matched.json").write_text(
        json.dumps(front, indent=1))
    print(f"wrote {DATA/'s6_frontier_u1matched.json'}", flush=True)


if __name__ == "__main__":
    main()
