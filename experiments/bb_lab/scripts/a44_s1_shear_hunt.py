#!/usr/bin/env python3
"""A44 S1 — which family reaches l + m - 6 at the next b = 0 frames?

Witness lane (SAT existence hunts on SMALL quotient frames, every witness
pulled back to the torus and kernel-verified there; witness weights are
UPPER bounds, UNSAT / timeouts are observations).

For the b = 0 frame (L, L), every sublattice Lambda' of Z^2 containing
Lambda = <(L,0),(0,L)> gives a shear frame Z^2/Lambda' = Z^2/<(L,0),(d,p)>
(p | L, 0 <= d < L) covered by the torus with index L/p; a nontrivial
weight-w cycle there pulls back to a Lambda'-periodic weight (L/p) w torus
cycle (odd index => nontriviality is inherited through the pullback of a
pairing partner; even index is checked in the kernel).  The (12,12)
D-minimizers are the (12;4,4) pullbacks + one band (A40 §16.2, A43 §5);
this scan asks the same question at (18,18) [target 30] and (24,24)
[target 42], and reproduces (12,12) as the control.

usage: a44_s1_shear_hunt.py L [time_limit_s]   (L in {12, 18, 24})
Output data/a44/s1_shear_L{L}.{json,log}.
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from a44_common import (  # noqa: E402
    LAB, DATA, Logger, member_code, Shear, hunt, pair_of, poly_str, picture,
    rss_guard, translation_canon, C,
)
from bb_lab.tower import validate_banked, v2i  # noqa: E402


def main():
    L = int(sys.argv[1])
    tl = float(sys.argv[2]) if len(sys.argv) > 2 else 120.0
    log = Logger(DATA / f"s1_shear_L{L}.log")
    t0 = time.time()
    validate_banked(LAB / "data")
    log("validate_banked: PASS")
    target = 2 * L - 6
    torus = member_code(L, L)
    log(f"=== frame ({L},{L}) n={torus.n} k={torus.k}; target l+m-6 = {target}")
    out = {"L": L, "target": target, "frames": [], "torus_orbits": {}}
    seen_orbits: dict = {}
    for p in range(1, L + 1):
        if L % p:
            continue
        idx = L // p
        if idx == 1 or target % idx:
            continue
        W = target // idx
        for d in range(L):
            sh = Shear(L, p, d)
            if not sh.covers(torus):
                continue
            tag = f"({L};{d},{p}) n={sh.code.n} idx={idx} W<={W}"
            if sh.code.k == 0:
                log(f"  {tag}: k=0, skip")
                out["frames"].append(dict(p=p, d=d, n=sh.code.n, k=0, idx=idx, W=W,
                                          status="k=0"))
                continue
            sols, status, dt = hunt(sh.code, W, nsol=6, time_limit=tl, tag=tag)
            rec = dict(p=p, d=d, n=sh.code.n, k=sh.code.k, orders=list(sh.orders),
                       idx=idx, W=W, status=status, t=dt, witnesses=[])
            for s in sols:
                v = C.vec_of(sh.code, [tuple(c) for c in s["cells"]])
                assert sh.code.is_cycle(v) and not sh.code.is_stab(v)
                vt = sh.pullback(v, torus)
                nontriv = not torus.is_stab(vt)
                cells = C.cells_of(torus, vt)
                nx, gx = C.gap_structure([c[1] for c in cells], L)
                ny, gy = C.gap_structure([c[2] for c in cells], L)
                f, g = pair_of(torus, vt)
                wrec = dict(w_shear=s["w"], w_torus=int(vt.sum()), nontrivial=bool(nontriv),
                            n_x=nx, n_y=ny, gap_x=gx, gap_y=gy,
                            blocks=[len(f), len(g)], f=poly_str(f), g=poly_str(g),
                            cells_torus=sorted(cells))
                if nontriv:
                    key = translation_canon(torus, vt)
                    wrec["orbit_key"] = str(key)
                    if key not in seen_orbits:
                        seen_orbits[key] = dict(p=p, d=d, w=int(vt.sum()), gap_x=gx, gap_y=gy,
                                                blocks=[len(f), len(g)], f=wrec["f"], g=wrec["g"])
                        log(f"    NEW torus orbit: w={int(vt.sum())} gaps=({gx},{gy}) "
                            f"blocks={[len(f), len(g)]}\n      f = {wrec['f']}\n      g = {wrec['g']}")
                        log(picture(torus, cells))
                rec["witnesses"].append(wrec)
            log(f"  {tag}: k={sh.code.k} orders={sh.orders} status={status} "
                f"witnesses={len(sols)} nontrivial-on-torus="
                f"{sum(1 for w in rec['witnesses'] if w['nontrivial'])} ({dt} s)")
            out["frames"].append(rec)
            rss_guard(tag)
            (DATA / f"s1_shear_L{L}.json").write_text(json.dumps(out, indent=1, default=str))
    out["torus_orbits"] = {str(k): v for k, v in seen_orbits.items()}
    out["wall_s"] = round(time.time() - t0, 1)
    log(f"=== ({L},{L}): {len(seen_orbits)} distinct nontrivial torus translation-orbits "
        f"at weight {target} from shear pullbacks; {out['wall_s']} s, RSS {rss_guard('end'):.2f} GB")
    (DATA / f"s1_shear_L{L}.json").write_text(json.dumps(out, indent=1, default=str))


if __name__ == "__main__":
    main()
