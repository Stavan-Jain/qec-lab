#!/usr/bin/env python3
"""A40 S9 — streaming port of the S7 ClosedMarch (the r=1 u=1
g46 completion: seeds 1-2, the two the frontier guard killed).

The S8 streaming rework (§13.6) applied to the closed lane: flush
the child layer to disk whenever it crosses SPLIT_THRESH during
expansion, then recurse per spilled part.  Cross-part dominance is
lost = over-exploration only (sound for emptiness: exploring MORE
paths can only find MORE closures, and the horizon readout runs on
every part that reaches h = m + 1).  The ru_maxrss lifetime-peak
guard flaw of the S7 module (§13.6) is patched out: the inherited
guard calls the module-level _rss_mb, replaced here with the
ps-based CURRENT RSS.

Port controls (subcommand `portcontrol`):
  C1 exact-identity: seed 3 at (24,18) u=1 k2 g46 with the split
     threshold disabled — popped/closed must EQUAL the banked
     s7_closed_pk2sh_g46_s3_4.json (same tree, same dominance).
  C2 spill-stability: same seed with SPLIT_THRESH = 20k — spills
     engaged, closed must stay 0, popped >= C1 (over-exploration).
  C3 closure-through-spills: the L12 (18,6) positive control
     (§12.4 control (i)) run through run_layers with SPLIT_THRESH
     = 64 — the certified closure must be DETECTED across spilled
     layers (readout path validated).
"""
from __future__ import annotations

import argparse
import json
import os
import pickle
import sys
import tempfile
import time
from pathlib import Path

LAB = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(LAB / "src"))
sys.path.insert(0, str(LAB / "scripts"))
DATA = LAB / "data" / "a40"

from a40_s6_frontier import wt, lsb, norm, seeds_full  # noqa: E402
import a40_s7_link as S7  # noqa: E402
from a40_s8_slip import _rss_mb  # noqa: E402  (ps-based)

# patch the flawed lifetime-peak guard out of the S7 module
S7._rss_mb = _rss_mb

PRE, HEAVY, POST = 0, 1, 2


class StreamClosedMarch(S7.ClosedMarch):
    SPLIT_THRESH = 180_000

    def _run_seed_layers(self, sid, w8, log):
        self._seed_norm, _ = norm(tuple(w8))
        self._sid = sid
        dyn = tuple(w8)
        allm = 0
        for r in w8:
            allm |= r
        layer = {(dyn, PRE, 0, 0): (self.u, lsb(allm),
                                    allm.bit_length() - 1)}
        return self._stream(layer, 1, sid, log, 0)

    def _stream(self, layer, h, sid, log, depth):
        t0 = time.time()
        while layer and h < self.m + 1:
            rss = _rss_mb()
            if rss > 2048:
                return f"RED: RSS {rss} MB at h={h}"
            nxt = {}
            parts = []

            def spill(d):
                fd, pth = tempfile.mkstemp(
                    suffix=f"_s9cl_d{depth}_h{h}.pkl")
                with os.fdopen(fd, "wb") as f:
                    pickle.dump(list(d.items()), f, protocol=4)
                parts.append(pth)

            nexp = 0
            for (dyn, phase, L, dlt), (g, xlo, xhi) in layer.items():
                self.expand_layer(g, dyn, phase, L, h, dlt, xlo,
                                  xhi, nxt)
                nexp += 1
                if len(nxt) > self.SPLIT_THRESH:
                    spill(nxt)
                    nxt = {}
                if nexp % 65536 == 0 and _rss_mb() > 2048:
                    for pth in parts:
                        os.remove(pth)
                    return f"RED: RSS {_rss_mb()} MB mid-layer h={h}"
            self.nodes_popped += len(layer)
            if parts:
                if nxt:
                    spill(nxt)
                nxt = None
                layer = None
                if log:
                    print(f"[s9cl sid{sid} d{depth}] h={h + 1} "
                          f"streamed into {len(parts)} parts rss "
                          f"{_rss_mb()}MB "
                          f"{round(time.time() - t0, 1)}s",
                          flush=True)
                try:
                    for pth in parts:
                        with open(pth, "rb") as f:
                            chunk = dict(pickle.load(f))
                        os.remove(pth)
                        ab = self._stream(chunk, h + 1, sid, log,
                                          depth + 1)
                        chunk = None
                        if ab:
                            return ab
                finally:
                    for pth in parts:
                        if os.path.exists(pth):
                            os.remove(pth)
                return None
            self.nodes_pushed += len(nxt)
            if log and (h % 6 == 0 or len(nxt) > 300000):
                print(f"[s9cl sid{sid} d{depth}] h={h + 1} states "
                      f"{len(nxt)} rss {_rss_mb()}MB "
                      f"{round(time.time() - t0, 1)}s", flush=True)
            layer = nxt
            h += 1
        if layer and h == self.m + 1:
            for (dyn, phase, L, dlt), (g, xlo, xhi) in layer.items():
                nd, _ = norm(dyn)
                if nd == self._seed_norm and \
                        dlt % self.ell == 0 and \
                        (phase == PRE or phase == POST):
                    self.closed.append((g, phase, L, dlt, sid))
        return None


def portcontrol():
    from bb_lab.tower import validate_banked
    validate_banked(LAB / "data")
    print("validate_banked: PASS", flush=True)
    sds = seeds_full(1)
    sds = [w8 for w8 in sds if sum(wt(r) for r in w8) == 1]

    ref = json.loads(
        (DATA / "s7_closed_pk2sh_g46_s3_4.json").read_text())
    print("C1 exact-identity: seed 3, g46, threshold off ...",
          flush=True)
    m1 = StreamClosedMarch(1, m=18, ell=24, kmax=2, whcap=14,
                           gcap=46, dcap=30)
    m1.SPLIT_THRESH = 10**9
    i1 = m1.run_layers(sds[3:4], sid0=3, log=False)
    ok1 = (i1["popped"] == ref["info"]["popped"]
           and i1["closed"] == ref["info"]["closed"] == 0
           and not i1["aborts"])
    print(f"  popped {i1['popped']} vs banked "
          f"{ref['info']['popped']}; closed {i1['closed']}: "
          f"{'PASS' if ok1 else 'FAIL'}", flush=True)
    assert ok1

    print("C2 spill-stability: seed 3, g46, threshold 20k ...",
          flush=True)
    m2 = StreamClosedMarch(1, m=18, ell=24, kmax=2, whcap=14,
                           gcap=46, dcap=30)
    m2.SPLIT_THRESH = 20_000
    i2 = m2.run_layers(sds[3:4], sid0=3, log=False)
    ok2 = (i2["closed"] == 0 and not i2["aborts"]
           and i2["popped"] >= i1["popped"])
    print(f"  popped {i2['popped']} (>= C1 {i1['popped']}: "
          f"over-exploration), closed {i2['closed']}: "
          f"{'PASS' if ok2 else 'FAIL'}", flush=True)
    assert ok2

    print("C3 closure-through-spills: L12 at (18,6), threshold 64 "
          "...", flush=True)
    _argv = sys.argv
    sys.argv = [sys.argv[0]]
    import numpy as np  # noqa: F401
    from a40_s6_drift import lift_phase
    from a40_s5_lightcore import Phase
    from a40_s4_phase_atlas import atlas
    sys.argv = _argv
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
        prof = sl_all[3:9]
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
        lights = [i for i, w in enumerate(prof) if w <= 7]
        j0 = min(lights, key=lambda i: prof[i])
        u0 = prof[j0]
        kk = 3 + j0
        w8 = []
        for blkk in (0, 1):
            for rr in range(kk, kk + 4):
                w8.append(sorted(fr.s(blkk, rr)))
        cols = [c for row in w8 for c in row]
        c0 = min(cols)
        w8m = tuple(sum(1 << (c - c0) for c in row) for row in w8)
        gsum = sum(prof)
        print(f"  L12 member: profile {prof}, block L={Lblk}, "
              f"seed slab {j0} (u={u0}), g_target={gsum + u0}",
              flush=True)
        c2 = StreamClosedMarch(u0, m=6, ell=18, kmax=Lblk,
                               whcap=max(prof), smax=8, dil=8,
                               gcap=gsum + u0, dcap=30)
        c2.SPLIT_THRESH = 64
        info = c2.run_layers([w8m], log=False)
        hits = [(g, phse, L, d) for (g, phse, L, d, _) in c2.closed]
        print(f"  march info: {info}; closed hits: {hits}",
              flush=True)
        ok = any(g <= gsum + u0 and d % 18 == 0
                 for (g, phse, L, d) in hits)
        assert ok, "L12 closure NOT detected through spills"
        exact = any(g == gsum + u0 and d == 0
                    for (g, phse, L, d) in hits)
        print(f"  closure detected across spilled layers (exact-g "
              f"hit: {exact}): PASS", flush=True)
        done = True
        break
    assert done, "no L12 control member found"
    print("portcontrol: ALL PASS", flush=True)


def prod(args):
    import gc
    gc.disable()
    t0 = time.time()
    from bb_lab.tower import validate_banked
    validate_banked(LAB / "data")
    print("validate_banked: PASS", flush=True)
    sds = seeds_full(args.u)
    sds = [w8 for w8 in sds if sum(wt(r) for r in w8) == args.u]
    lo, hi = map(int, args.seeds.split(":"))
    print(f"s9 closed (streaming) u={args.u} m=18 ell=24 kmax=2 "
          f"gcap={args.gcap}: seeds [{lo}, {hi}) of {len(sds)}",
          flush=True)
    c = StreamClosedMarch(args.u, m=18, ell=24, kmax=2, whcap=14,
                          gcap=args.gcap, dcap=30)
    info = c.run_layers(sds[lo:hi], sid0=lo)
    print(f"info: {info}", flush=True)
    for (g, phase, L, dlt, sid) in sorted(c.closed)[:20]:
        print(f"  CLOSED: g={g} phase={phase} L={L} dlt={dlt} "
              f"sid={sid}", flush=True)
    out = dict(
        params=dict(u=args.u, m=18, ell=24, kmax=2, whcap=14,
                    gcap=args.gcap, smax=3, dil=4, dcap=30,
                    seeds=[lo, hi]),
        info=info,
        closed=[dict(g=g, phase=p_, L=L, dlt=d, sid=s)
                for (g, p_, L, d, s) in sorted(c.closed)])
    p = DATA / f"s9_closed_pk2sh_g{args.gcap}_s{lo}_{hi}.json"
    p.write_text(json.dumps(out, indent=1))
    print(f"wrote {p} ({round(time.time() - t0, 1)} s)", flush=True)


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("portcontrol")
    s = sub.add_parser("prod")
    s.add_argument("--u", type=int, default=1)
    s.add_argument("--gcap", type=int, required=True)
    s.add_argument("--seeds", type=str, required=True)
    args = ap.parse_args()
    if args.cmd == "portcontrol":
        portcontrol()
    else:
        prod(args)


if __name__ == "__main__":
    main()
