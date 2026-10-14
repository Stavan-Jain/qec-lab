#!/usr/bin/env python3
"""A40 S4 gate 4 — compact triviality: every sub-rate-2 x-compact
cycle has a COMPACT stabilizer generator, so it is trivial at EVERY
torus (l, p) it fits — closing the forall-l hole left by the atlas's
single-embedding classification.

Method: for pair (P,Q), stabilizers are v = (Pbar s, Qbar s).  Given a
compact cycle v, a compact generator (if one exists) is UNIQUE and is
produced by the deterministic x-march on the block-1 equation
conv(Pbar, s) = v1 (solve for the most-advanced column; unique top
term asserted), then verified against conv(Qbar, s) = v2 exactly.  If
the march flushes to zero and the block-2 equation checks, v = d(s)
with compact s: reducing mod x^l - 1 shows v is a stabilizer at every
(l >= extent(s) + span, p).  Any cycle whose march fails is reported
and tested per-l directly over the family-relevant list (FALLBACK —
none expected).
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

from a40_s4_phase_atlas import (  # noqa: E402
    DATA, PAIRS, atlas, bar, rot,
)
from bb_lab.tower import TowerCode, validate_banked  # noqa: E402


def conv_cols(T, cols, p, lo, hi):
    """conv of a support T (list of (i,j)) with columns dict c->int,
    returned as dict over [lo, hi]."""
    out = {}
    for c in range(lo, hi + 1):
        acc = 0
        for (i, j) in T:
            acc ^= rot(cols.get(c - i, 0), j, p)
        if acc:
            out[c] = acc
    return out


def march_generator(Tsolve, Tcheck, v1cols, v2cols, p, margin=12):
    """Solve conv(Tsolve, s) = v1 by the forced x-march; verify
    conv(Tcheck, s) = v2.  Returns s columns dict or None."""
    imin = min(i for i, _ in Tsolve)          # most-advanced column c-imin
    tops = [(i, j) for (i, j) in Tsolve if i == imin]
    assert len(tops) == 1, ("non-unique top term", Tsolve)
    j0 = tops[0][1]
    rest = [(i, j) for (i, j) in Tsolve if i != imin]
    cs = ([c for c, x in v1cols.items() if x]
          + [c for c, x in v2cols.items() if x])
    lo, hi = min(cs) - margin, max(cs) + margin
    s = {}
    for c in range(lo, hi + 1):
        acc = v1cols.get(c, 0)
        for (i, j) in rest:
            acc ^= rot(s.get(c - i, 0), j, p)
        val = rot(acc, -j0, p)
        if val:
            s[c - imin] = val
    # flush check: s must vanish beyond hi - margin/2
    if any(c > hi - margin // 2 for c in s):
        return None
    # verify both equations exactly on a window covering everything
    ss = sorted(s) or [0]
    LO = min(lo, ss[0]) - 6
    HI = max(hi, ss[-1]) + 6
    if conv_cols(Tsolve, s, p, LO, HI) != {c: x for c, x in
                                           v1cols.items() if x}:
        return None
    if conv_cols(Tcheck, s, p, LO, HI) != {c: x for c, x in
                                           v2cols.items() if x}:
        return None
    return s


def verify_on_torus(pair_name, p, pts, s, cache={}):
    """Independent check: on (lstar, p) with lstar >= all extents + 8,
    the placed generator's stabilizer row-combination equals v."""
    Psupp, Qsupp = PAIRS[pair_name]
    cs = [c for c, _, _ in pts] + sorted(s)
    lstar = max(24, max(cs) - min(cs) + 10)
    key = (pair_name, p, lstar)
    if key not in cache:
        cache[key] = TowerCode(
            f"vt{key}", (lstar, p),
            frozenset((i % lstar, j % p) for i, j in Psupp),
            frozenset((i % lstar, j % p) for i, j in Qsupp))
    code = cache[key]
    c0 = min(cs)
    v = np.zeros(code.n, dtype=np.uint8)
    for c, y, blk in pts:
        v[blk * code.ng + code.G.index(((c - c0) % lstar, y % p))] ^= 1
    w = np.zeros(code.n, dtype=np.uint8)
    for c, col in s.items():
        for y in range(p):
            if col >> y & 1:
                w ^= code.HX[code.G.index(((c - c0) % lstar, y % p))]
    return (v == w).all()


def main():
    t0 = time.time()
    validate_banked(LAB / "data")
    print("validate_banked: PASS")
    out = {"rows": []}
    total, solved = 0, 0
    for pair_name in ("AB", "BAbar"):
        Psupp, Qsupp = PAIRS[pair_name]
        Pb, Qb = bar(Psupp), bar(Qsupp)
        for p in (2, 3, 4, 5, 6, 7, 8):
            rows, _ = atlas(pair_name, p, 2 * p - 1,
                            max_states=30_000_000,
                            max_paths=3_000_000, keep_pts=True)
            n_ok, smax = 0, 0
            for r in rows:
                assert not r["nontrivial"]
                v1cols, v2cols = {}, {}
                for c, y, blk in r["pts"]:
                    d = v1cols if blk == 0 else v2cols
                    d[c] = d.get(c, 0) ^ (1 << y)
                s = march_generator(Pb, Qb, v1cols, v2cols, p)
                total += 1
                if s is not None:
                    assert verify_on_torus(pair_name, p, r["pts"], s)
                    n_ok += 1
                    solved += 1
                    if s:
                        smax = max(smax, max(s) - min(s) + 1)
                else:
                    print(f"  FALLBACK needed: {pair_name} p={p} "
                          f"weight {r['weight']} extent {r['extent']}")
            out["rows"].append(dict(pair=pair_name, p=p,
                                    n_cycles=len(rows),
                                    n_compactly_trivial=n_ok,
                                    max_generator_extent=smax))
            print(f"{pair_name} p={p}: {n_ok}/{len(rows)} sub-rate-2 "
                  f"cycles have COMPACT generators "
                  f"(max s-extent {smax})", flush=True)
    print(f"\nTOTAL: {solved}/{total} compactly trivial")
    assert solved == total, "some cycles need per-l fallback — " \
        "the forall-l claim is NOT closed for those"
    print("VERDICT: every x-compact cycle below rate 2 (p <= 8, both "
          "pairs) is COMPACTLY trivial => trivial at every (l, p) "
          "that fits it: the atlas's below-rate-2 verdict holds at "
          "EVERY l, not just the test embedding.")
    out["wall_s"] = round(time.time() - t0, 1)
    (DATA / "s4_compact_triviality.json").write_text(
        json.dumps(out, indent=1))
    print(f"wrote {DATA/'s4_compact_triviality.json'} "
          f"({out['wall_s']} s)")


if __name__ == "__main__":
    main()
