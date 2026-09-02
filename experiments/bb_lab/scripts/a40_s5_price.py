#!/usr/bin/env python3
"""A40 S5 — pricing gate for the l = 18 periodic-leg completion
(p = 6, 7, 8 sheared frames, dense branch) and the p = 9, 10 reach.

Per S4 §9.7.1 the periodic leg at (l, p) splits: orbits with an x-gap
>= 4 lift to x-compact cylinder phases (the atlas closed p <= 8, both
lanes, forall l), and DENSE orbits (all x-gaps <= 3) need per-frame
censuses.  At l = 18, p in {6,7,8} the frames have n = 36p in
{216, 252, 288} > 192 (walk-kernel cap) — the descent lane is the only
node-exact route: fold an even SNF axis by 2, census the base at
W = 2p-1, then enumerate ALL cover cycles <= W as lifts (tau-lane +
deep fibers).  This script prices that plan BEFORE any run:

  - recompute k of every shear (must match the S4 triage — asserted);
  - detect duplicate normalized codes among shears (same orders +
    transported supports => same census, run once);
  - per k > 0 frame: fold options (even axes), base k, base kappa,
    census_nodes(kappa, W) x class count = the node bill;
  - calibrate nodes/s with one real stab census on the smallest base;
  - extend the k-map to p = 9, 10 (new triage rows) and price those.

Verdicts: GREEN < 30 min est., AMBER < 3 h, RED beyond (per frame and
total).  No claims here — cost verdicts only.
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np

LAB = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(LAB / "src"))
sys.path.insert(0, str(LAB / "scripts"))

from bb_lab import cosetbz  # noqa: E402
from bb_lab.tower import (  # noqa: E402
    TowerCode, census_nodes, validate_banked,
)
from a38_c37xx_freeze import census_pass  # noqa: E402
from a40_s4_phase_triage import quotient_code  # noqa: E402

DATA = LAB / "data" / "a40"
L = 18


def fold_options(code, orders):
    """Even axes and the folded base code for each."""
    outs = []
    for ax in (0, 1):
        if orders[ax] % 2:
            continue
        newo = list(orders)
        newo[ax] //= 2
        A2 = frozenset((e[0] % newo[0], e[1] % newo[1])
                       for e in code.A.support)
        B2 = frozenset((e[0] % newo[0], e[1] % newo[1])
                       for e in code.B.support)
        base = TowerCode(f"{code.name}_f{ax}", tuple(newo), A2, B2)
        try:
            cosetbz.disjoint_info_sets(base.HX)
            info = True
        except RuntimeError:
            info = False
        outs.append(dict(axis=ax, orders=newo, n=base.n, k=base.k,
                         kappa=base.kappa, info_sets=info))
    return outs


def main():
    t0 = time.time()
    validate_banked(LAB / "data")
    print("validate_banked: PASS", flush=True)

    triage = json.load(open(DATA / "s4_phase_triage.json"))
    tri_k = {}
    for r in triage["rows"]:
        if r["l"] == L:
            for d in r["shears"]:
                tri_k[(r["p"], d)] = r["k"]

    out = {"frames": [], "p910_kmap": []}
    seen = {}
    frames = []
    for p in (6, 7, 8):
        for d in range(L):
            code, o = quotient_code(L, p, d)
            assert tri_k[(p, d)] == code.k, (p, d, code.k, tri_k[(p, d)])
            if code.k == 0:
                continue
            key = (tuple(o), tuple(sorted(code.A.support)),
                   tuple(sorted(code.B.support)))
            if key in seen:
                out["frames"].append(dict(p=p, d=d, k=code.k,
                                          duplicate_of=seen[key]))
                print(f"p={p} d={d:2d}: k={code.k} DUPLICATE of "
                      f"{seen[key]} (same normalized code)", flush=True)
                continue
            seen[key] = [p, d]
            frames.append((p, d, code, o))
    print(f"\n{len(frames)} distinct k>0 frames at p=6..8 (of 30 "
          f"triage rows)", flush=True)

    # calibrate: one stab census on the smallest foldable base
    cal = None
    binp = cosetbz.build_kernel()
    for p, d, code, o in frames:
        for fo in fold_options(code, o):
            if fo["info_sets"] and fo["n"] <= 192:
                newo = fo["orders"]
                A2 = frozenset((e[0] % newo[0], e[1] % newo[1])
                               for e in code.A.support)
                B2 = frozenset((e[0] % newo[0], e[1] % newo[1])
                               for e in code.B.support)
                base = TowerCode("cal", tuple(newo), A2, B2)
                W = 2 * p - 1
                t1 = time.time()
                census_pass(binp, base, [("S", np.zeros(base.n, np.uint8))],
                            W, "s5cal")
                wall = time.time() - t1
                nodes = census_nodes(base.kappa, W)
                cal = dict(p=p, d=d, base_n=base.n, kappa=base.kappa,
                           W=W, nodes=nodes, wall_s=round(wall, 2),
                           nodes_per_s=round(nodes / max(wall, 1e-9)))
                print(f"\ncalibration: base of (18,{p},{d}) n={base.n} "
                      f"kappa={base.kappa} W={W}: {nodes:.3g} nodes in "
                      f"{wall:.2f}s -> {cal['nodes_per_s']:.3g} nodes/s",
                      flush=True)
                break
        if cal:
            break
    nps = cal["nodes_per_s"]

    tot_s = 0.0
    for p, d, code, o in frames:
        W = 2 * p - 1
        fos = fold_options(code, o)
        best = None
        for fo in fos:
            if not fo["info_sets"] or fo["n"] > 192:
                continue
            ncls = (1 << fo["k"]) - 1 if fo["k"] else 0
            nodes = census_nodes(fo["kappa"], W) * (ncls + 1)
            est = nodes / nps
            if best is None or est < best["est_s"]:
                best = dict(fo, classes=ncls, nodes=nodes,
                            est_s=round(est, 1))
        row = dict(p=p, d=d, k=code.k, orders=list(o), n=code.n,
                   folds=fos, best=best)
        out["frames"].append(row)
        if best is None:
            print(f"p={p} d={d:2d}: k={code.k} orders={o} n={code.n} "
                  f"-> NO usable fold (info sets / cap) — needs depth-2 "
                  f"or other route", flush=True)
        else:
            tot_s += best["est_s"]
            print(f"p={p} d={d:2d}: k={code.k} orders={o} n={code.n} "
                  f"-> fold ax{best['axis']} base n={best['n']} "
                  f"k={best['k']} kappa={best['kappa']}; "
                  f"{best['classes']+1} census bases x W={W}: "
                  f"{best['nodes']:.3g} nodes ~ {best['est_s']}s",
                  flush=True)

    # p = 9, 10 k-map extension + pricing sketch
    print("\n--- p = 9, 10 k-map (new triage rows) ---", flush=True)
    for p in (9, 10):
        ks = {}
        for d in range(L):
            code, o = quotient_code(L, p, d)
            ks.setdefault(code.k, []).append(d)
        for k, ds in sorted(ks.items()):
            out["p910_kmap"].append(dict(l=L, p=p, k=k, shears=ds))
            print(f"l={L} p={p}: k={k} at d={ds}", flush=True)
        for k, ds in sorted(ks.items()):
            if k == 0:
                continue
            d = ds[0]
            code, o = quotient_code(L, p, d)
            fos = fold_options(code, o)
            for fo in fos:
                if fo["info_sets"] and fo["n"] <= 192:
                    W = 2 * p - 1
                    ncls = (1 << fo["k"]) - 1
                    nodes = census_nodes(fo["kappa"], W) * (ncls + 1)
                    print(f"   p={p} d={d} k={k}: base n={fo['n']} "
                          f"kappa={fo['kappa']} W={W}: {nodes:.3g} "
                          f"nodes ~ {nodes/nps:.0f}s", flush=True)
                    break

    v = ("GREEN" if tot_s < 1800 else
         "AMBER" if tot_s < 10800 else "RED")
    print(f"\nTOTAL census bill p=6..8 (distinct frames, best folds): "
          f"~{tot_s:.0f}s of kernel walks (+ fiber lanes, typically "
          f"comparable) -> {v}", flush=True)
    out["calibration"] = cal
    out["total_est_s"] = round(tot_s, 1)
    out["verdict"] = v
    out["wall_s"] = round(time.time() - t0, 1)
    (DATA / "s5_price.json").write_text(json.dumps(out, indent=1))
    print(f"wrote {DATA/'s5_price.json'} ({out['wall_s']} s)")


if __name__ == "__main__":
    main()
