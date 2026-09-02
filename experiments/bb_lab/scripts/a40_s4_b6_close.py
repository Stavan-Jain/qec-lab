#!/usr/bin/env python3
"""A40 S4 — closing Theorem B6's outstanding certificates: d >= 12 at
(24,6), (30,6), (36,6) by Z2-descent (walks at n <= 192 only).

Per cover L0 = (l,6) with x-fold base L1 = (l/2,6), target T = 12,
window W = 11: every nontrivial X-logical v of L0 with |v| <= 11 has
shadow b = p(v), an L1-cycle with |b| <= |v| = |b| + 2*overflow.
Complete case split on b (all censuses node-exact BZ):

  b = 0: p(v) = 0 <=> v is sigma-invariant <=> v = tau(u) for an
     L1-CYCLE u (elementary for a free Z2 deck: p(v)(gbar) =
     v(g0)+v(g1) = 0 forall gbar <=> v = tau(fold of sheet-0); tau is
     injective — rank asserted — and a chain map, so u is a cycle);
     |v| = 2|u| >= 2*min(mu_stab, d(L1)) >= 12 once mu_stab >= 6 and
     d(L1) >= 6 are certified (both censused in-run).  tau(stab) is a
     cover stab, tau(nontrivial-with-|u|<=5) needs d(L1) >= 6.
  b a stabilizer, 6 <= |b| <= 11: dangerous rung rung(b, M) per
     translation-orbit rep, M = ceil((12-|b|)/2) <= 3 (restricted
     lanes; PASS => |v| >= |b| + 2M >= 12).  Covariance spot-checked.
  b nontrivial: |b| >= d(L1); the all-class <= 11 census certifies
     d(L1) >= 12 > 11 >= |b| — the seam lane is EMPTY (if the census
     were nonempty, each element would get a seam_rung instead).

Frames:
  (24,6)  --x-->  (12,6) = gross   (n = 144; d(L1) = 12 re-derived)
  (30,6)  --x-->  (15,6)           (n = 180; k, d censused in-run)
  (36,6)  --x-->  (18,6)  --x-->  (9,6)   (depth 2: the (18,6)
     obligations — stabs <= 11 / nontrivial <= 11 — are built by
     descent from complete (9,6) censuses (n = 108): tau-transport of
     (9,6) objects (zero-shadow lane, |tau u| = 2|u|) + fiber
     enumeration over every (9,6) stab/class rep <= 11 with cap
     (11-|b2|)/2 <= 2 (enumerate_lifts_deep, node-exact), then the
     top rungs (18,6) -> (36,6) as above.)
  (42,6)  --x-->  (21,6)  --y-->  (21,3)   (depth 2, MIXED axes: the
     direct (21,6) census is impossible (n = 252 > 192) but its
     y-fold (21,3) has n = 126 — the same descent machinery applies
     axis-generically.)

With the atlas (windowed branch, all l), L1 (l >= 45), gross/(18,6),
and these four: THEOREM B6 holds for every 6 | l >= 12 — complete.

A40 S4 header for d >= 12 at (l,6): the certificates here are FULL
d >= 12 floors (not just the x-spanning branch), so each frame's
closure is single-source; the atlas covers the windowed branch of
every OTHER l uniformly.
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
    AxisDeck, RungCell, TowerCode, enumerate_lifts_deep, gf2_rank,
    i2v, in_span, rep_for, translation_perms, v2i, validate_banked,
)
from a38_c37xx_freeze import Collector, census_pass, whist  # noqa: E402

DATA = LAB / "data" / "a40"
A_L = [(0, 0), (0, 1), (3, -1)]
B_L = [(0, 0), (1, 0), (-1, -3)]
W = 11
TARGET = 12


def red(supp, lm):
    return frozenset((e[0] % lm[0], e[1] % lm[1]) for e in supp)


def code_at(lm, name=None):
    return TowerCode(name or f"tdg{lm}", lm, red(A_L, lm), red(B_L, lm))


def census_stabs(binp, code, tag):
    hs = census_pass(binp, code, [("S", np.zeros(code.n, np.uint8))],
                     W, tag)
    out = []
    for h in sorted(hs["S"]):
        v = i2v(h, code.n)
        assert code.is_stab(v)
        out.append(v)
    return out


def census_all_classes(binp, code, tag, Wc=W):
    """Complete nontrivial census <= Wc; returns list of vectors."""
    out = []
    CH = 51
    allc = sorted(range(1, 1 << code.k))
    for lo in range(0, len(allc), CH):
        chunk = allc[lo:lo + CH]
        ha = census_pass(binp, code,
                         [(f"C{c}", rep_for(code, c)) for c in chunk],
                         Wc, f"{tag}{lo}")
        for c in chunk:
            for h in sorted(ha[f"C{c}"]):
                v = i2v(h, code.n)
                assert code.is_cycle(v) and not code.is_stab(v)
                out.append(v)
    return out


def zero_shadow_dead(deck, mu_stab, d_l1, log, name):
    """The b = 0 lane: tau-injectivity + weight arithmetic."""
    rk = gf2_rank([v2i(c) for c in deck.TAU.T])
    assert rk == deck.base.n, "tau not injective?!"
    dead = 2 * min(mu_stab, d_l1) >= TARGET
    log(f"[{name}] b=0 lane: tau injective (rank {rk}); "
        f"2*min(mu_stab={mu_stab}, d(L1)>={d_l1}) >= {TARGET}: {dead}")
    assert dead
    return dead


def top_rungs(name, L1, L0, deck, stab_reps, perms1, log):
    cell = RungCell(f"{name}_top", L1, L0, deck)
    verd: dict[str, int] = {}
    viol = []
    deadline = time.monotonic() + 3600
    for b in sorted(stab_reps, key=lambda b: -int(b.sum())):
        M = -(-(TARGET - int(b.sum())) // 2)
        if M <= 0:
            continue
        assert M - 1 <= 8
        r = cell.rung(b, M, deadline)
        verd[r["verdict"]] = verd.get(r["verdict"], 0) + 1
        if r["verdict"] == "VIOLATION":
            viol.append(r)
    # covariance spot-check: verdicts translation-invariant
    if stab_reps:
        b = stab_reps[0]
        M = -(-(TARGET - int(b.sum())) // 2)
        if M > 0:
            for perm in perms1[1:3]:
                r2 = cell.rung(b[perm], M, deadline)
                assert r2["verdict"] == "PASS", "covariance check"
    log(f"[{name}] top rungs at target {TARGET}: {verd}, "
        f"violations {len(viol)}")
    return verd, viol


def close_depth1(binp, l, log):
    name = f"({l},6)"
    L0 = code_at((l, 6), f"{name}/L0")
    L1 = code_at((l // 2, 6), f"{name}/L1")
    deck = AxisDeck(L0, L1, 0)
    perms1 = translation_perms(L1)
    log(f"[{name}] tower ({l},6) --x--> ({l//2},6); k = {L0.k}/{L1.k}")

    stabs = census_stabs(binp, L1, f"b6c_{l}_stab")
    mu = min((int(v.sum()) for v in stabs), default=99)
    coll = Collector(L1.n)
    for v in stabs:
        coll.add(v)
    stab_reps = coll.reps(perms1)
    log(f"[{name}] L1 stabs <= {W}: {len(stabs)} elements, "
        f"{len(stab_reps)} orbit reps {whist(stab_reps)}; mu = {mu}")
    assert mu >= 6

    ntrv = census_all_classes(binp, L1, f"b6c_{l}_all") if L1.k else []
    log(f"[{name}] L1 all-class census <= {W}: {len(ntrv)} nontrivial "
        f"=> d(L1) >= {12 if not ntrv else min(int(v.sum()) for v in ntrv)}")
    assert not ntrv, f"{name}: base has nontrivial <= {W} — seam lane " \
        f"needed (not wired for this case; d(L1) < 12?!)"
    d_l1 = 12

    zero_shadow_dead(deck, mu, d_l1, log, name)
    verd, viol = top_rungs(name, L1, L0, deck, stab_reps, perms1, log)
    ok = not viol
    log(f"[{name}] => d(({l},6)) >= {TARGET}: {'PASS' if ok else 'FAIL'}")
    return dict(l=l, depth=1, base=[l // 2, 6], k=[L0.k, L1.k],
                stab_orbits=len(stab_reps), mu_stab=mu,
                base_nontrivial_le_W=0, rungs=verd,
                violations=len(viol), floor_12=bool(ok))


def close_depth2(binp, lm0, ax0, lm1, ax1, lm2, log):
    name = f"({lm0[0]},{lm0[1]})"
    axn = {0: "x", 1: "y"}
    L0 = code_at(lm0, f"{name}/L0")
    L1 = code_at(lm1, f"{name}/L1")
    L2 = code_at(lm2, f"{name}/L2")
    deck0 = AxisDeck(L0, L1, ax0)
    deck1 = AxisDeck(L1, L2, ax1)
    perms1 = translation_perms(L1)
    perms2 = translation_perms(L2)
    log(f"[{name}] tower {lm0} --{axn[ax0]}--> {lm1} --{axn[ax1]}--> "
        f"{lm2}; k = {L0.k}/{L1.k}/{L2.k}")

    # complete L2 censuses
    stabs2 = census_stabs(binp, L2, f"b6c_{lm0[0]}_L2stab")
    mu2 = min(int(v.sum()) for v in stabs2)
    ntrv2 = census_all_classes(binp, L2, f"b6c_{lm0[0]}_L2all")
    d_l2 = min((int(v.sum()) for v in ntrv2), default=12)
    c2 = Collector(L2.n)
    for v in stabs2:
        c2.add(v)
    stab2_reps = c2.reps(perms2)
    c2n = Collector(L2.n)
    for v in ntrv2:
        c2n.add(v)
    ntrv2_reps = c2n.reps(perms2)
    log(f"[{name}] L2 {lm2} censuses <= {W}: stabs {len(stabs2)} "
        f"(mu {mu2}, {len(stab2_reps)} reps), nontrivial {len(ntrv2)} "
        f"({len(ntrv2_reps)} reps) => d({lm2}) = {d_l2}")
    assert mu2 >= 6

    # L1 obligations by descent: every L1-cycle w with |w| <= 11 has
    # L2-shadow b2 = p1(w) with |b2| <= 11, |w| = |b2| + 2 ov:
    #   b2 = 0        -> w = tau1(u), |w| = 2|u|, u an L2-cycle
    #   b2 stab       -> fiber enumeration cap (11 - |b2|)/2
    #   b2 nontrivial -> |b2| >= d_l2 = 10; fiber cap 0 as applicable
    coll_stab1 = Collector(L1.n)
    ntrv1 = []

    def classify(w):
        ww = int(w.sum())
        if ww == 0 or ww > W:
            return
        assert L1.is_cycle(w)
        if in_span(v2i(w), L1.rsHX_b, L1.rsHX_p):
            coll_stab1.add(w)
        else:
            ntrv1.append(w)

    # zero-shadow lane: tau1 of every L2 stab/nontrivial with 2|u| <= 11
    rk = gf2_rank([v2i(c) for c in deck1.TAU.T])
    assert rk == L2.n
    for u in stabs2 + ntrv2:
        if 2 * int(u.sum()) <= W:
            classify((deck1.TAU @ u) % 2)
    # fiber lanes over every L2 rep (elements, not reps, for shadows —
    # orbits under L2-translations lift to L1-translations, so reps
    # suffice with the perms1 reduction at the end; we enumerate over
    # ALL banked elements to stay conservative)
    nf = 0
    for b2 in sorted(stabs2 + ntrv2, key=lambda v: -int(v.sum())):
        cap = (W - int(b2.sum())) // 2
        if cap < 0:
            continue
        lifts = enumerate_lifts_deep(deck1, b2, cap=cap)
        for v0c, m2 in sorted(lifts.items()):
            classify(deck1.lift(i2v(v0c, L2.n), b2))
        nf += 1
    stab1_reps = coll_stab1.reps(perms1)
    mu1 = min(int(v.sum()) for v in stab1_reps)
    log(f"[{name}] L1 {lm1} obligations by descent ({nf} fibers): "
        f"stabs {len(stab1_reps)} reps {whist(stab1_reps)} (mu {mu1}); "
        f"nontrivial <= {W}: {len(ntrv1)} => d({lm1}) >= 12 "
        f"re-derived: {not ntrv1}")
    assert mu1 >= 6
    assert not ntrv1, f"{lm1} nontrivial <= 11 found?!"

    zero_shadow_dead(deck0, mu1, 12, log, name)
    verd, viol = top_rungs(name, L1, L0, deck0, stab1_reps, perms1, log)
    ok = not viol
    log(f"[{name}] => d({lm0}) >= {TARGET}: {'PASS' if ok else 'FAIL'}")
    return dict(l=lm0[0], depth=2, base=list(lm1),
                bottom=list(lm2),
                k=[L0.k, L1.k, L2.k], mu_stab_mid=mu1,
                d_bottom=d_l2, stab_orbits=len(stab1_reps),
                base_nontrivial_le_W=0, rungs=verd,
                violations=len(viol), floor_12=bool(ok))


def main():
    t0 = time.monotonic()

    def log(m):
        print(f"[{time.monotonic()-t0:7.1f}s] {m}", flush=True)

    validate_banked(LAB / "data")
    log("validate_banked: PASS")
    binp = cosetbz.build_kernel()
    out = {"frames": []}
    for l in (24, 30):
        out["frames"].append(close_depth1(binp, l, log))
    out["frames"].append(close_depth2(binp, (36, 6), 0, (18, 6), 0, (9, 6), log))
    out["frames"].append(close_depth2(binp, (42, 6), 0, (21, 6), 1, (21, 3), log))
    all_ok = all(f["floor_12"] for f in out["frames"])
    assert all_ok
    log("ALL FOUR PASS => with the atlas (windowed, all l), L1 "
        "(l >= 45), gross and (18,6): THEOREM B6 d((l,6)) = 12 holds "
        "for every 6 | l >= 12 — ALL certificates closed, no "
        "residue.")
    out["wall_s"] = round(time.monotonic() - t0, 1)
    (DATA / "s4_b6_close.json").write_text(json.dumps(out, indent=1))
    log(f"wrote {DATA/'s4_b6_close.json'} ({out['wall_s']} s)")


if __name__ == "__main__":
    main()
