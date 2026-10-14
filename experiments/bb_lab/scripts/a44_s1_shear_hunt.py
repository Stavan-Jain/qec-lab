#!/usr/bin/env python3
"""A44 S1 — which family reaches the target weight at the next frames?

Witness lane (SAT existence hunts on SMALL quotient frames, every witness
pulled back to the torus and kernel-verified there; witness weights are
UPPER bounds, UNSAT / timeouts are observations).

For the frame (L, M), every sublattice Lambda' of Z^2 containing
Lambda = <(L,0),(0,M)> has a Hermite basis <(a,0),(d,p)> with a | L,
p | M, 0 <= d < a, and (0,M) in Lambda' iff a | (M/p) d.  The shear frame
Z^2/Lambda' is covered by the torus with index LM/(ap); a nontrivial
weight-w cycle there pulls back to a Lambda'-periodic torus cycle of
weight (LM/ap) w (odd index => nontriviality is inherited; every pullback
is checked in the kernel anyway).  The (12,12) D-minimizers are the
(12;4,4) pullbacks + one band (A40 §16.2, A43 §5); this scan asks the same
question at (18,18) [target 30], (24,24) [42], (18,12) [26, 28], (24,18).

usage: a44_s1_shear_hunt.py L M target [time_limit_s] [amin]
   amin: scan only a >= amin (default 1); the a = L, p | M slice is the
   original (L;d,p) scan.
Output data/a44/s1_shear_{L}x{M}_W{target}.{json,log}.
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from a44_common import (  # noqa: E402
    LAB, DATA, Logger, member_code, hunt, pair_of, poly_str, picture,
    rss_guard, translation_canon, C,
)
from a40_s4_phase_triage import snf2  # noqa: E402
from bb_lab.tower import TowerCode, validate_banked  # noqa: E402

A_L = [(0, 0), (0, 1), (3, -1)]
B_L = [(0, 0), (1, 0), (-1, -3)]


class HShear:
    """Z^2 / <(a,0),(d,p)> with the fixed pair, SNF coordinates."""

    def __init__(self, a, p, d):
        self.a, self.p, self.d = a, p, d
        M = [[a, 0], [d, p]]
        D, U, V = snf2(M)
        self.V = V
        self.orders = (abs(D[0][0]), abs(D[1][1]))
        assert self.orders[0] * self.orders[1] == a * p

        def tr(supp):
            cnt: dict = {}
            for e in supp:
                k = self.phi(e)
                cnt[k] = cnt.get(k, 0) ^ 1
            return frozenset(k for k, v in cnt.items() if v)
        self.code = TowerCode(f"shear(a={a},p={p},d={d})", self.orders, tr(A_L), tr(B_L))

    def phi(self, e):
        V, (o1, o2) = self.V, self.orders
        return ((e[0] * V[0][0] + e[1] * V[1][0]) % o1,
                (e[0] * V[0][1] + e[1] * V[1][1]) % o2)

    def covers(self, target):
        L, M = target.G.orders
        return self.phi((L, 0)) == (0, 0) and self.phi((0, M)) == (0, 0)

    def pullback(self, v, target):
        L, M = target.G.orders
        import numpy as np
        vt = np.zeros(target.n, dtype=np.uint8)
        code = self.code
        for blk in range(2):
            for x in range(L):
                for y in range(M):
                    if v[blk * code.ng + code.G.index(self.phi((x, y)))]:
                        vt[blk * target.ng + target.G.index((x, y))] = 1
        assert target.is_cycle(vt)
        return vt


def main():
    L, M, target = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
    tl = float(sys.argv[4]) if len(sys.argv) > 4 else 120.0
    amin = int(sys.argv[5]) if len(sys.argv) > 5 else 1
    amax = int(sys.argv[6]) if len(sys.argv) > 6 else L
    log = Logger(DATA / f"s1_shear_{L}x{M}_W{target}_a{amin}-{amax}.log")
    t0 = time.time()
    validate_banked(LAB / "data")
    log("validate_banked: PASS")
    torus = member_code(L, M)
    log(f"=== frame ({L},{M}) n={torus.n} k={torus.k}; target {target}; a >= {amin}")
    out = {"L": L, "M": M, "target": target, "frames": [], "torus_orbits": {}}
    seen_orbits: dict = {}
    jpath = DATA / f"s1_shear_{L}x{M}_W{target}_a{amin}-{amax}.json"
    for a in range(amin, amax + 1):
        if L % a:
            continue
        for p in range(1, M + 1):
            if M % p:
                continue
            idx = (L * M) // (a * p)
            if idx == 1 or target % idx:
                continue
            W = target // idx
            for d in range(a):
                if ((M // p) * d) % a:
                    continue
                sh = HShear(a, p, d)
                assert sh.covers(torus)
                tag = f"<({a},0),({d},{p})> n={sh.code.n} idx={idx} W<={W}"
                if sh.code.k == 0:
                    log(f"  {tag}: k=0, skip")
                    out["frames"].append(dict(a=a, p=p, d=d, n=sh.code.n, k=0, idx=idx,
                                              W=W, status="k=0"))
                    continue
                sols, status, dt = hunt(sh.code, W, nsol=6, time_limit=tl, tag=tag)
                rec = dict(a=a, p=p, d=d, n=sh.code.n, k=sh.code.k, orders=list(sh.orders),
                           idx=idx, W=W, status=status, t=dt, witnesses=[])
                for s in sols:
                    v = C.vec_of(sh.code, [tuple(c) for c in s["cells"]])
                    assert sh.code.is_cycle(v) and not sh.code.is_stab(v)
                    vt = sh.pullback(v, torus)
                    nontriv = not torus.is_stab(vt)
                    cells = C.cells_of(torus, vt)
                    nx, gx = C.gap_structure([c[1] for c in cells], L)
                    ny, gy = C.gap_structure([c[2] for c in cells], M)
                    f, g = pair_of(torus, vt)
                    wrec = dict(w_shear=s["w"], w_torus=int(vt.sum()), nontrivial=bool(nontriv),
                                n_x=nx, n_y=ny, gap_x=gx, gap_y=gy, blocks=[len(f), len(g)],
                                f=poly_str(f), g=poly_str(g), cells_torus=sorted(cells))
                    if nontriv:
                        key = translation_canon(torus, vt)
                        wrec["orbit_key"] = str(key)
                        if key not in seen_orbits:
                            seen_orbits[key] = dict(a=a, p=p, d=d, w=int(vt.sum()), gap_x=gx,
                                                    gap_y=gy, blocks=[len(f), len(g)],
                                                    f=wrec["f"], g=wrec["g"],
                                                    cells=sorted(cells))
                            log(f"    NEW torus orbit: w={int(vt.sum())} gaps=({gx},{gy}) "
                                f"blocks={[len(f), len(g)]}\n      f = {wrec['f']}\n"
                                f"      g = {wrec['g']}")
                            log(picture(torus, cells))
                    rec["witnesses"].append(wrec)
                log(f"  {tag}: k={sh.code.k} orders={sh.orders} status={status} "
                    f"witnesses={len(sols)} nontrivial-on-torus="
                    f"{sum(1 for w in rec['witnesses'] if w['nontrivial'])} ({dt} s)")
                out["frames"].append(rec)
                rss_guard(tag)
                out["torus_orbits"] = {str(k): v for k, v in seen_orbits.items()}
                jpath.write_text(json.dumps(out, indent=1, default=str))
    out["wall_s"] = round(time.time() - t0, 1)
    log(f"=== ({L},{M}) target {target}: {len(seen_orbits)} distinct nontrivial torus "
        f"translation-orbits from shear pullbacks; {out['wall_s']} s, "
        f"RSS {rss_guard('end'):.2f} GB")
    jpath.write_text(json.dumps(out, indent=1, default=str))


if __name__ == "__main__":
    main()
