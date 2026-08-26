#!/usr/bin/env python3
"""A40 S4 — the local-reduction pruning pilot (L-W's first mechanical
leg): the COMPLETE periodic-orbit landscape at l = 12, periods p <= 8,
all shears, with the class-minimality filter applied.

Universe: a period-p closed walk of the y-transfer system at l = 12
with x-drift d is an X-cycle of the sheared frame Z^2/<(12,0),(d,p)>
(the §9.4 correspondence; x-winding orbits INCLUDED, unlike the
compact atlas).  For every such frame we census ALL cycles of weight
<= 2p (node-exact coset-BZ over every class, stabilizers included)
and keep the y-spanning ones (phase row-support has all cyclic gaps
<= 3 in the ORIGINAL y coordinate, recovered through the SNF
transform — tiling makes walk-spanning equal phase-spanning).

Filter (sound, one-sided): a cycle v cannot be class-minimal if some
LOCAL stabilizer z (catalog: all single H_X rows and all row-pair
sums of weight <= 10 of the frame) has 2|v cap z| > |z| (then
|v + z| < |v| in the same class).  Trivial cycles are never
class-minimal (class 0's minimum is 0), so any trivial cycle the
catalog MISSES is a measured gap of the catalog, not a soundness
error.

Question: does any sub-rate-2 (weight < 2p) y-spanning orbit SURVIVE
the filter?  Expected: none — combs and their relatives all contain
more than half of a tooth.  Surviving rate-2 orbits (weight = 2p) are
the tight population (the L12/stack species and frames' minima).
Output: the survivor spectrum per frame + the pruning statistics.
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
    i2v, rep_for, v2i, validate_banked,
)
from a38_c37xx_freeze import census_pass  # noqa: E402
from a40_s4_phase_triage import quotient_code, snf2  # noqa: E402

DATA = LAB / "data" / "a40"
L = int(sys.argv[1]) if len(sys.argv) > 1 else 12
P_MAX = int(sys.argv[2]) if len(sys.argv) > 2 else 8


def inv2(V):
    """Inverse of a unimodular 2x2 integer matrix."""
    det = V[0][0] * V[1][1] - V[0][1] * V[1][0]
    assert det in (1, -1)
    return [[det * V[1][1], -det * V[0][1]],
            [-det * V[1][0], det * V[0][0]]]


def y_coords(code, p, d, v):
    """Original-cylinder y coordinates (mod p) of a vector's support.

    The quotient's normalized coordinates f relate to exponent coords
    e by f = e V (a40_s4_phase_triage.quotient_code); so e = f V^-1,
    and the original y is e[1] mod p (the phase period)."""
    M = [[L, 0], [d, p]]
    D, U, V = snf2(M)
    Vi = inv2(V)
    ys = set()
    for i in np.nonzero(v)[0]:
        f = code.G.from_index(int(i) % code.ng)
        y = (f[0] * Vi[0][1] + f[1] * Vi[1][1]) % p
        ys.add(y)
    return ys


def y_spanning(code, p, d, v):
    ys = sorted(y_coords(code, p, d, v))
    if not ys:
        return False
    if len(ys) == p:
        return True
    if len(ys) == 1:
        return p - 1 <= 3
    gaps = [(ys[(i + 1) % len(ys)] - ys[i]) % p - 1
            for i in range(len(ys))]
    return max(gaps) <= 3


def build_catalog(code, wmax=10):
    """All single H_X rows + all row-pair sums of weight <= wmax, as
    packed ints (the local-generator catalog)."""
    rows = [v2i(r) for r in code.HX]
    cat = set()
    for r in rows:
        cat.add(r)
    for i in range(len(rows)):
        for j in range(i + 1, len(rows)):
            z = rows[i] ^ rows[j]
            if z and bin(z).count("1") <= wmax:
                cat.add(z)
    return [(z, bin(z).count("1")) for z in cat]


def prunable(vi, catalog):
    for z, wz in catalog:
        if 2 * bin(vi & z).count("1") > wz:
            return True
    return False


def b1_closure_exists(p, d, rmax=2000):
    """Does the (p,d)-phase tile ANY b=1 member (6r+6, 6r)?  Needs
    p | 6r and (6r/p)*d = 0 mod 6r+6.  (The x-order constraint of an
    x-winding phase is conservatively ignored — one-sided: closure
    chances are reported a fortiori.)"""
    for r in range(1, rmax + 1):
        m = 6 * r
        if m % p:
            continue
        if ((m // p) * d) % (6 * r + 6) == 0:
            return r
    return None


def x_winds(p, d, pts):
    """Does the phase use the x-wrap?  Re-place the point set on the
    doubled-l shear frame; not-a-cycle => x-winding."""
    code2, _ = quotient_code(2 * L, p, d)
    M = [[L, 0], [d, p]]
    D, U, V = snf2(M)
    Vi = inv2(V)
    M2 = [[2 * L, 0], [d, p]]
    D2, U2, V2 = snf2(M2)
    v = np.zeros(code2.n, dtype=np.uint8)
    for blk, f1, f2 in pts:
        e = (f1 * Vi[0][0] + f2 * Vi[1][0],
             f1 * Vi[0][1] + f2 * Vi[1][1])
        g = ((e[0] * V2[0][0] + e[1] * V2[1][0]) % code2.G.orders[0],
             (e[0] * V2[0][1] + e[1] * V2[1][1]) % code2.G.orders[1])
        v[blk * code2.ng + code2.G.index(g)] ^= 1
    return not code2.is_cycle(v)


def main():
    t0 = time.time()
    validate_banked(LAB / "data")
    print("validate_banked: PASS", flush=True)
    binp = cosetbz.build_kernel()
    out = {"frames": [], "L": L}
    n_orbits = n_below = n_below_pruned = 0
    survivors_below = []
    tight = {}
    for p in range(2, P_MAX + 1):
        for d in range(L):
            code, o = quotient_code(L, p, d)
            W = 2 * p
            if code.n > cosetbz.NMAX:
                print(f"p={p} d={d}: SKIPPED (n = {code.n} > kernel "
                      f"cap {cosetbz.NMAX}) — needs the descent lane",
                      flush=True)
                out["frames"].append(dict(p=p, d=d, k=code.k,
                                          skipped="n>NMAX"))
                continue
            cat = build_catalog(code)
            try:
                cosetbz.disjoint_info_sets(code.HX)
            except RuntimeError:
                print(f"p={p} d={d}: SKIPPED (no disjoint info-set "
                      f"pair) — frame recorded as uncovered", flush=True)
                out["frames"].append(dict(p=p, d=d, k=code.k,
                                          skipped="no-info-sets"))
                continue
            # stabilizers (class 0) + all nontrivial classes
            offsets = [("S", np.zeros(code.n, np.uint8))]
            offsets += [(f"C{c}", rep_for(code, c))
                        for c in range(1, 1 << code.k)]
            found = []
            CH = 51
            for lo in range(0, len(offsets), CH):
                hits = census_pass(binp, code, offsets[lo:lo + CH], W,
                                   f"pp_{p}_{d}_{lo}")
                for lab, _ in offsets[lo:lo + CH]:
                    for h in sorted(hits[lab]):
                        found.append((lab, h))
            rows = []
            for lab, h in found:
                v = i2v(h, code.n)
                w = int(v.sum())
                if w == 0 or not code.is_cycle(v):
                    continue
                if not y_spanning(code, p, d, v):
                    continue
                triv = code.is_stab(v)
                pr = prunable(h, cat)
                n_orbits += 1
                if w < 2 * p:
                    n_below += 1
                    if pr:
                        n_below_pruned += 1
                    else:
                        pts = []
                        for i in np.nonzero(v)[0]:
                            blk, gi = divmod(int(i), code.ng)
                            pts.append((blk,)
                                       + tuple(code.G.from_index(gi)))
                        survivors_below.append(
                            dict(p=p, d=d, weight=w,
                                 trivial=bool(triv),
                                 x_winds=bool(x_winds(p, d, pts)),
                                 b1_closure_r=b1_closure_exists(p, d),
                                 pts=sorted(pts)))
                elif w == 2 * p and not pr:
                    tight[(p, d)] = tight.get((p, d), 0) + 1
                rows.append((w, triv, pr))
            below = [r for r in rows if r[0] < 2 * p]
            if below or any(w == 2 * p for w, _, _ in rows):
                surv_t = tight.get((p, d), 0)
                print(f"p={p} d={d:2d} (k={code.k}): {len(rows)} "
                      f"y-spanning orbits <= {W}; below-rate-2: "
                      f"{len(below)} (pruned {sum(1 for r in below if r[2])}"
                      f"/{len(below)}); tight rate-2 survivors: {surv_t}",
                      flush=True)
            out["frames"].append(dict(
                p=p, d=d, k=code.k, n_spanning=len(rows),
                n_below=len(below),
                n_below_pruned=sum(1 for r in below if r[2]),
                tight_survivors=tight.get((p, d), 0)))
    print(f"\nTOTALS: {n_orbits} y-spanning orbits (all p <= 8, all "
          f"shears, l = {L}); below rate 2: {n_below}, of which "
          f"PRUNED by the local catalog: {n_below_pruned}")
    if survivors_below:
        import collections
        c = collections.Counter(
            (s["p"], s["d"], s["weight"], s["trivial"], s["x_winds"],
             s["b1_closure_r"]) for s in survivors_below)
        print("SUB-RATE-2 SURVIVORS by "
              "(p, d, w, trivial, x_winds, b1_closure_r):")
        for kk, njm in sorted(c.items()):
            print(f"   {kk}: {njm}")
        all_wind = all(s["x_winds"] for s in survivors_below)
        none_close = all(s["b1_closure_r"] is None
                         for s in survivors_below)
        print(f"ALL survivors x-wind: {all_wind}; NONE closes around "
              f"any b=1 member (r <= 2000): {none_close}")
        if all_wind and none_close:
            print("=> consistent with the architecture: cheap phases "
                  "exist but are x-winding non-closers — the wall "
                  "accounting must consume the closure arithmetic; no "
                  "b=1 member can ride any of them as a pure phase.")
    else:
        print(f"VERDICT: every y-spanning periodic orbit below rate 2 "
              f"(period <= {P_MAX}, any shear, l = {L}) is PRUNED by the "
              "local-reduction catalog — no class-minimal logical can "
              "ride any of them. The L-W certificate's periodic-orbit "
              "leg holds on the swept range.")
    out["totals"] = dict(orbits=n_orbits, below=n_below,
                         below_pruned=n_below_pruned,
                         survivors_below=survivors_below)
    out["wall_s"] = round(time.time() - t0, 1)
    (DATA / f"s4_prune_pilot_l{L}.json").write_text(json.dumps(out, indent=1))
    print(f"wrote {DATA}/s4_prune_pilot_l{L}.json ({out['wall_s']} s)")


if __name__ == "__main__":
    main()
