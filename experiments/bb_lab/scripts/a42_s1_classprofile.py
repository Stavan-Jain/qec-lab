#!/usr/bin/env python3
"""A42 S1g — the class-weight profile of the compact cylinder:
min realization weight PER nonzero H-class (63 classes, dim H = 6),
at p = 6 and p = 12, by the sigma-conditioned h-DP.

The profile refines the floor (its min is the floor) and is the
object Stage 2's member assembly consumes: member logicals pair
x- and y-cylinder classes, so the per-class weights (not just the
global min) drive the torus bounds.  Also identifies WHICH class
ray is extremal: the tangency tower's socle (pi-multiples), the
unit part, or the transverse F4 summand (Theorem H inventory).

Scope: omega-support <= smax slots per sigma (as in s1_jointdp).
"""
from __future__ import annotations

import itertools
import json
import sys
import time
from pathlib import Path

import numpy as np

LAB = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(LAB / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import a42_lib as AL  # noqa: E402
import importlib.util

_spec = importlib.util.spec_from_file_location(
    "a42_s1_syzygy", Path(__file__).parent / "a42_s1_syzygy.py")
SY = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(SY)

_spec2 = importlib.util.spec_from_file_location(
    "a42_s1_jointdp", Path(__file__).parent / "a42_s1_jointdp.py")
JD = importlib.util.module_from_spec(_spec2)
_spec2.loader.exec_module(JD)

DATA = LAB / "data" / "a42"


def profile(p: int, Wx: int, smax: int, log=print):
    a = AL.v2(p)
    t0 = time.time()
    tab, mprime, L = JD.joint_cost_table(p)
    ow = SY.OmegaWindow(a, Wx)
    H, T, Z, Lrows, dimZ, ncls = ow.functionals()
    log(f"p={p}: classes {ncls}")
    assert ncls == (6 if a >= 1 else 4)
    Ld = ow.L.dim
    slots = [(blk, c) for blk in (0, 1) for c in range(Wx)]
    best = {}
    n_sigma = 0
    for s in range(1, smax + 1):
        for pat in itertools.combinations(slots, s):
            if min(c for (_, c) in pat) != 0:
                continue
            cols = sorted({c for (_, c) in pat})
            if any(b_ - a_ > 4 for a_, b_ in zip(cols, cols[1:])):
                continue
            sols = SY.solve_affine_in_pattern(ow, H, Lrows, pat, 0)
            if sols.size == 0:
                continue
            pair = (Lrows @ sols.T) % 2
            if not pair.any():
                continue
            dimS = sols.shape[0]
            if dimS > 14:
                continue
            for mb in range(1, 1 << dimS):
                acc = np.zeros(ow.nbits, dtype=np.uint8)
                t = mb
                i = 0
                while t:
                    if t & 1:
                        acc ^= sols[i]
                    t >>= 1
                    i += 1
                cv = (Lrows @ acc) % 2
                key = int("".join(map(str, cv)), 2)
                if key == 0:
                    continue
                n_sigma += 1
                sigma_cols = {}
                for (blk, c) in pat:
                    lam = 0
                    for ii in range(Ld):
                        if acc[ow.bit(blk, c, ii)]:
                            lam |= 1 << ii
                    if lam:
                        sigma_cols[(blk, c)] = lam
                if not sigma_cols:
                    continue
                cur = best.get(key)
                # pre-filter vs this class's current best
                nh_ = 1 << (1 << a)
                fsum = sum(min(tab[(zp, lam)] for zp in range(nh_))
                           for lam in sigma_cols.values())
                if cur is not None and fsum >= cur:
                    continue
                cmin = min(c for (_, c) in sigma_cols)
                cmaxx = max(c for (_, c) in sigma_cols)
                val = JD.dp_min_for_sigma(p, sigma_cols, tab, a,
                                          cmin - 4, cmaxx + 4)
                if cur is None or val < cur:
                    best[key] = val
    # histogram
    from collections import Counter
    hist = Counter(best.values())
    dt = round(time.time() - t0, 1)
    log(f"p={p}: {n_sigma} sigma priced; classes reached "
        f"{len(best)}/{(1 << ncls) - 1}; weight histogram {dict(hist)} "
        f"({dt} s)")
    return {"p": p, "classes": ncls, "n_sigma": n_sigma,
            "profile": {str(k): v for k, v in sorted(best.items())},
            "histogram": {str(k): v for k, v in sorted(hist.items())},
            "wall_s": dt}


def main():
    out = {}
    for p, Wx, smax in ((6, 10, 7), (12, 10, 7)):
        out[str(p)] = profile(p, Wx, smax,
                              log=lambda s: print(s, flush=True))
    (DATA / "s1_classprofile.json").write_text(json.dumps(out, indent=1))
    print(f"wrote {DATA/'s1_classprofile.json'}", flush=True)


if __name__ == "__main__":
    main()
