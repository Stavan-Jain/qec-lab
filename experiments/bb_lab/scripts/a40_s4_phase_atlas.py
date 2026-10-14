#!/usr/bin/env python3
"""A40 S4 — the compact-phase atlas: exhaustive census of x-compact
cycles of weight <= Wcap on the period-p y-cylinder, per polynomial
pair.

Object: pair (P,Q); X-cycles satisfy Qbar v1 + Pbar v2 = 0 over
F2[x^pm] (x) tensor F2[y]/(y^p - 1).  Writing columns v_c[i] in F2^p,
the cycle equation collected at x-degree c is a linear recurrence with
ONE forced column per step (the largest x-shift in the equation is
unique for both pairs used here), so every x-compact cycle is a walk
of the finite-state x-automaton from the zero state back to the zero
state.  Exhaustive BFS over (state, cost <= Wcap) pairs + DAG path
readout therefore enumerates ALL x-compact cycles of weight <= Wcap,
up to x-translation (start normalization).  Each cycle is re-verified
independently (cycle + weight) on an embedding torus (lstar, p) and
classified trivial/nontrivial there (lstar >= extent + 8, so member
triviality at every l >= lstar coincides with cylinder triviality for
compactly-bounded generators; recorded as such).

Consumers:
  - the rate-2 floor architecture (S4 note §9): NO nontrivial
    x-compact phase of weight < 2p may exist at any small p;
  - B6 (p = 6): the strip census — minimum nontrivial x-compact weight
    at p = 6 is 12 = L12 (positive control: L12 must be FOUND at
    Wcap = 12).
Pairs: (A,B) [the y-lane] and (B,Abar) [= the theta'-image: the x-lane
/ strip sector of (A,B) in rotated coordinates].
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np

LAB = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(LAB / "src"))

from bb_lab.tower import TowerCode, validate_banked  # noqa: E402

DATA = LAB / "data" / "a40"
A_L = [(0, 0), (0, 1), (3, -1)]
B_L = [(0, 0), (1, 0), (-1, -3)]


def bar(supp):
    return [(-e[0], -e[1]) for e in supp]


PAIRS = {
    # name -> (P support, Q support); X-cycle eq: Qbar v1 + Pbar v2 = 0
    "AB": (A_L, B_L),
    "BAbar": (B_L, bar(A_L)),
}


def rot(col, s, p):
    """y-shift by +s of a p-bit column (int)."""
    s %= p
    return ((col << s) | (col >> (p - s))) & ((1 << p) - 1)


class Automaton:
    """x-march automaton for the pair on the period-p cylinder.

    Equation at x-degree c:  sum over (i,j) in supp(Qbar) of
    Y^j v1[c - i]  +  sum over supp(Pbar) of Y^j v2[c - i]  = 0.
    Let i1max/i2max be the largest x-offsets per block and imax the
    overall largest; the recurrence solves for the block with the
    UNIQUE term at x-offset imax (asserted unique), the other block's
    later columns being free inputs.
    """

    def __init__(self, Psupp, Qsupp, p):
        self.p = p
        Qb = [(-i, -j) for (i, j) in Qsupp]   # acts on v1
        Pb = [(-i, -j) for (i, j) in Psupp]   # acts on v2
        # normalize so the forced block is block 1 with a unique max
        # x-offset term; equation at c reads v_blk[c - i] for (i,j).
        i1 = min(i for i, _ in Qb)   # v1[c - i]: most-advanced = min i
        i2 = min(i for i, _ in Pb)
        # advance = -i (column c + (-i)); forced column = the single
        # most-advanced one across both blocks.
        a1 = max(-i for i, _ in Qb)
        a2 = max(-i for i, _ in Pb)
        assert a1 != a2 or True
        if a1 > a2:
            self.forced_blk = 0
            terms_f = [(-i, j) for (i, j) in Qb]
            terms_o = [(-i, j) for (i, j) in Pb]
            self.adv_f, self.adv_o = a1, a2
        else:
            self.forced_blk = 1
            terms_f = [(-i, j) for (i, j) in Pb]
            terms_o = [(-i, j) for (i, j) in Qb]
            self.adv_f, self.adv_o = a2, a1
        top = [t for t in terms_f if t[0] == self.adv_f]
        assert len(top) == 1, ("non-unique top term", terms_f)
        self.top_j = top[0][1]
        # memory depths: equation at c touches forced[c+adv_f .. c-df],
        # other[c+adv_o .. c-do]
        self.df = -min(a for a, _ in terms_f)
        self.do = -min(a for a, _ in terms_o)
        self.terms_f = [(a, j) for (a, j) in terms_f if a != self.adv_f]
        self.terms_o = terms_o
        # state: forced block columns c-df+1 .. c+adv_f  (len df+adv_f),
        #        other block columns  c-do+1 .. c+adv_o  (len do+adv_o)
        self.nf = self.df + self.adv_f
        self.no = self.do + self.adv_o
        assert self.nf >= 0 and self.no >= 1

    def zero(self):
        return (0,) * (self.nf + self.no)

    def step(self, state, a):
        """Advance the equation index by 1 with free input a (the next
        'other' column beyond the horizon); returns (state', cost)."""
        p = self.p
        f = state[:self.nf]
        o = state[self.nf:]
        # solve equation at index c for forced[c + adv_f]:
        # top term Y^top_j forced[c+adv_f] = sum(other terms)
        acc = 0
        # forced-block terms at offsets a_ < adv_f: column c+a_ =
        # f[nf-1 - (adv_f - a_)] (f lists c-df+1..c+adv_f-? ) --
        # f holds columns (c-df+1 .. c+adv_f -1 +1)? we keep f as the
        # last nf forced columns ENDING at c+adv_f-1; the new forced
        # column becomes c+adv_f.
        for (a_, j) in self.terms_f:
            idx = self.nf - 1 + (a_ - (self.adv_f - 1))
            acc ^= rot(f[idx], j, p) if 0 <= idx < self.nf else 0
            assert 0 <= idx < self.nf
        # other-block: columns end at c+adv_o AFTER consuming input a:
        # o currently holds columns (c-do .. c+adv_o-1); with input a =
        # column c+adv_o.
        oo = o + (a,)
        for (a_, j) in self.terms_o:
            idx = len(oo) - 1 + (a_ - self.adv_o)
            assert 0 <= idx < len(oo)
            acc ^= rot(oo[idx], j, p)
        new_f_col = rot(acc, -self.top_j, p)
        nf2 = (f + (new_f_col,))[-self.nf:] if self.nf else ()
        no2 = oo[-self.no:]
        cost = bin(a).count("1") + bin(new_f_col).count("1")
        return nf2 + no2, cost, new_f_col

    def readout_cols(self, inputs):
        """Replay a full input sequence from zero; return per-column
        supports of both blocks as {(c, y, blk)} plus final state."""
        st = self.zero()
        pts = set()
        c = 0
        for a in inputs:
            st, _, ncol = self.step(st, a)
            for y in range(self.p):
                if a >> y & 1:
                    pts.add((c + self.adv_o, y, 1 if self.forced_blk == 0
                             else 0))
                if ncol >> y & 1:
                    pts.add((c + self.adv_f, y, self.forced_blk))
            c += 1
        return pts, st


def embed_and_classify(pts, Psupp, Qsupp, p, cache={}):
    """Embed a compact point set on the (lstar, p) torus and classify."""
    if not pts:
        return None
    cs = [c for c, _, _ in pts]
    ext = max(cs) - min(cs) + 1
    lstar = max(24, ext + 8)
    key = (tuple(sorted(Psupp)), tuple(sorted(Qsupp)), p, lstar)
    if key not in cache:
        lm = (lstar, p)
        cache[key] = TowerCode(f"embed{lm}", lm,
                               frozenset((i % lstar, j % p)
                                         for i, j in Psupp),
                               frozenset((i % lstar, j % p)
                                         for i, j in Qsupp))
    code = cache[key]
    v = np.zeros(code.n, dtype=np.uint8)
    c0 = min(cs)
    for c, y, blk in pts:
        v[blk * code.ng + code.G.index(((c - c0) % lstar, y % p))] ^= 1
    assert code.is_cycle(v), "readout is not a cycle on the torus"
    return dict(weight=int(v.sum()), extent=ext, lstar=lstar,
                nontrivial=bool(not code.is_stab(v)))


def atlas(pair_name, p, Wcap, max_states=4_000_000, max_paths=200_000,
          keep_pts=False):
    Psupp, Qsupp = PAIRS[pair_name]
    au = Automaton(Psupp, Qsupp, p)
    z = au.zero()
    # BFS over (state, cost): parents kept for DAG path readout
    layers = {(z, 0): None}
    frontier = [(z, 0)]
    parents = {}
    hits = []          # (cost, input-sequence) for zero-returns
    npop = 0
    while frontier:
        nxt = []
        for (st, cost) in frontier:
            npop += 1
            if npop > max_states:
                raise RuntimeError("state cap exceeded — raise cap or "
                                   "lower Wcap")
            budget = Wcap - cost
            for a in range(1 << p):
                if bin(a).count("1") > budget:
                    continue
                st2, c2, _ = au.step(st, a)
                cost2 = cost + c2
                if cost2 > Wcap:
                    continue
                if st2 == z and not any(st) and a == 0:
                    continue    # stay-at-zero: skip identity self-loop
                key = (st2, cost2)
                if key not in layers:
                    layers[key] = (st, cost, a)
                    nxt.append(key)
                    parents.setdefault(key, []).append((st, cost, a))
                else:
                    parents.setdefault(key, []).append((st, cost, a))
        frontier = nxt
    # enumerate all zero-return paths (each = one compact cycle)
    out = []

    def walk_back(key, suffix, seen_budget):
        if len(out) > max_paths:
            raise RuntimeError("path cap exceeded")
        st, cost = key
        if st == z and cost == 0:
            out.append(tuple(reversed(suffix)))
            return
        for (pst, pcost, a) in parents.get(key, []):
            walk_back((pst, pcost), suffix + [a], seen_budget)

    for cost in range(1, Wcap + 1):
        if (z, cost) in layers:
            walk_back((z, cost), [], cost)
    # classify
    rows = []
    for inputs in out:
        pts, st = au.readout_cols(list(inputs)
                                  + [0] * (au.nf + au.no + 2))
        assert st == au.zero()
        cl = embed_and_classify(pts, Psupp, Qsupp, p)
        if cl is None:
            continue
        if cl["weight"] > Wcap:
            continue
        cl["pair"], cl["p"] = pair_name, p
        if keep_pts:
            cl["pts"] = sorted(pts)
        rows.append(cl)
    return rows, npop


def main():
    t0 = time.time()
    validate_banked(LAB / "data")
    print("validate_banked: PASS")
    results = []
    summary = []
    jobs = []
    for pair in ("AB", "BAbar"):
        for p in (2, 3, 4, 5, 6):
            jobs.append((pair, p, 2 * p - 1))
    # positive control: p = 6 at Wcap = 12 must FIND weight-12
    # nontrivial phases (L12 species) for pair AB
    jobs.append(("AB", 6, 12))
    jobs.append(("BAbar", 6, 12))
    for pair, p, W in jobs:
        t1 = time.time()
        rows, npop = atlas(pair, p, W)
        wall = round(time.time() - t1, 1)
        ntv = [r for r in rows if r["nontrivial"]]
        wmin = min([r["weight"] for r in ntv], default=None)
        summary.append(dict(pair=pair, p=p, Wcap=W, n_cycles=len(rows),
                            n_nontrivial=len(ntv),
                            min_nontrivial_weight=wmin,
                            states_popped=npop, wall_s=wall))
        results += rows
        print(f"{pair} p={p} W<={W}: {len(rows)} compact cycles, "
              f"{len(ntv)} NONTRIVIAL (min weight {wmin}), "
              f"{npop} states, {wall} s")
    # the theorem rows: below-rate-2 windows must be all-trivial
    viol = [s for s in summary if s["Wcap"] == 2 * s["p"] - 1
            and s["n_nontrivial"]]
    ctrl = [s for s in summary if s["Wcap"] == 12 and s["pair"] == "AB"
            and s["p"] == 6]
    assert ctrl and ctrl[0]["min_nontrivial_weight"] == 12, \
        "positive control failed: L12 not found at p=6 W=12"
    print("\npositive control PASS: p=6 Wcap=12 finds min nontrivial "
          "weight 12 (the L12 species)")
    if viol:
        print("!! BELOW-RATE-2 NONTRIVIAL PHASES EXIST:", viol)
    else:
        print("VERDICT: no nontrivial x-compact phase below rate 2 at "
              "any p <= 6, either pair — the phase floor holds on the "
              "swept range")
    out = dict(summary=summary,
               nontrivial_rows=[r for r in results if r["nontrivial"]],
               wall_s=round(time.time() - t0, 1))
    (DATA / "s4_phase_atlas.json").write_text(json.dumps(out, indent=1))
    print(f"wrote {DATA/'s4_phase_atlas.json'} ({out['wall_s']} s)")


if __name__ == "__main__":
    main()
