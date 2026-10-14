#!/usr/bin/env python3
"""A40 P2: the tour-de-gross (2,1) member [[432,12,?]] — the W-sweep.

Tower (the yxy ordering: top rung (R), the A36 shape; P1 lattice +
pricing GREEN at W = 22, `data/a40/pricing.json`):

  L0 (18,12) [[432,12,?]]   d_conj = 24 (arXiv:2506.03094, conjecture)
     |  y-fold (deck sigma = y^6), twisted, (R): k 12->12, exact both,
     |  sigma* = id, rank p0* = 6  => SEAM = im p0* (63 classes)
  L1 (18,6) [[216,12,?]]    d(L1) = 12 EXACT (measured here, W16 run)
     |  x-fold (x^9), twisted, non-(R): k 12->8, exact_cover only,
     |  rank p1* = 6; rank(p1* o p0*) = 2 => S1' = p1*(SEAM) has just
     |  3 nontrivial classes — the seam shadows
  L2 (9,6)  [[108,8,?]]     d(L2) = 10 EXACT (measured here)
     |  y-fold (y^3), twisted, (R)
  L3 (9,3)  [[54,8,?]]      d(L3) = 6 EXACT (direct BZ)

No banked distances existed for L1/L2/L3 — every input is certified
in-run by COMPLETE weight-window census generation (node-exact BZ walks
at n <= 108; composed fibers cross-gated).

ARCHITECTURE (v2, direct-L2-primary — the v1 descent-primary ran the
banked W = 16 shakedown; at W = 22 the L3 fiber layer explodes, and
rank(p1* o p0*) = 2 makes the direct route GREEN):

  L2 censuses DIRECT (coset-BZ, exact node asserts):
     stab <= W (1 base);  S1' 3 classes <= W;  im-p1* 63 classes
     <= WNT = W-4 (feeds the nontrivial-L1 window censuses);
     ALL 255 classes <= WALL = 16 (kernel-shift windows + tau2 sources
     <= W/2 + the <= 14 gate slice + d(L2) exact).
  L2 GATE: the <= 14 slice re-derived by DESCENT from L3 (tau3-family +
     composed fibers over the L3 <= 14 censuses) — orbit key sets must
     equal the direct slices exactly (stab AND nontrivial).
  L1 obligations by descent from the L2 censuses (completeness by the
     shadow-class law [P v] = p1*[v]):
     stab1 <= W    from stab2 fibers + tau2   (stab shadows are stabs);
     seam1 <= W    from stab2 + S1' fibers + tau2  (seam shadows lie in
                   S1' u {0} u stab);
     ntrv1 <= WNT  from stab2 + im-p1* fibers + tau2 (all classes) —
                   the top kernel-shift windows + d(L1) exact.
  L1 GATE: the <= 12 L1-stab census re-derived through the INDEPENDENT
     y-quotient (18,3) deck.
  Rungs at the top only (target W+2): dangerous per stab1 orbit rep,
     seam per seam1 orbit rep; kernel-shift lane windows <= WNT
     (all-cycle censuses complete there); every candidate re-verified
     in-line; covariance; X<->Z duality spot-check.

Questions (TARGET = W + 2; all cycles even — parity scope):
  W = 16 (banked, v1): d >= 18. + tau0-witness 24 => 18 <= d <= 24.
  W = 22: all rungs PASS => d >= 24; with the weight-24 tau0-witness
          (= tau0 of a weight-12 non-SEAM L1-logical, re-verified) =>
          **d([[432,12]]) = 24 EXACT, certificate tier** — the
          conjectured value at (2,1). Any violation is an explicit
          verified logical <= 22 => d < 24 (conjecture refuted at
          (2,1)); the complete sweep makes the minimum exact.

Usage: python a40_tdg432_close.py 22 --census-only   (phase 1)
       python a40_tdg432_close.py 22 --rungs-only    (phase 2)
       python a40_tdg432_close.py 16                 (one shot)

Output: data/a40/tdg432/sweep_W{W}*.json + ckpt_W{W}_*.jsonl
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
    AxisDeck, RungCell, TowerCode, batch_keys, colspace,
    enumerate_lifts_deep, fold_support, gf2_rank, h1_map, i2v, in_span,
    kernel_ints, perm_for, rep_for, rref_ints, span_eq, span_points,
    translation_mat, translation_perms, v2i, validate_banked,
)
from a38_c37xx_freeze import (  # noqa: E402
    Collector, KernelShift, census_pass, row_lift_v0, whist,
)

DATA = LAB / "data" / "a40" / "tdg432"
DATA.mkdir(parents=True, exist_ok=True)

A_L = [(0, 0), (0, 1), (3, -1)]
B_L = [(0, 0), (1, 0), (-1, -3)]
CH = 51                      # coset-class offsets per BZ pass


def red(supp, lm):
    return frozenset((e[0] % lm[0], e[1] % lm[1]) for e in supp)


def key_set(vecs, perms):
    if not vecs:
        return set()
    return {bytes(k) for k in
            batch_keys(np.array(vecs, dtype=np.uint8), perms)}


def main() -> None:
    W = int(sys.argv[1]) if len(sys.argv) > 1 else 22
    assert W in (16, 18, 22)
    census_only = "--census-only" in sys.argv
    rungs_only = "--rungs-only" in sys.argv
    assert not (census_only and rungs_only)
    TARGET = W + 2
    WNT = W - 4      # nontrivial-window bound (top kernel-shift windows)
    WALL = min(16, W)  # all-class census bound (L1-fiber ks windows)
    t0 = time.monotonic()
    out: dict = {"W": W, "target": TARGET, "WNT": WNT, "WALL": WALL,
                 "phase": ("census" if census_only else
                           "rungs" if rungs_only else "full"),
                 "architecture": "v2 direct-L2-primary"}

    def log(msg: str) -> None:
        print(f"[{time.monotonic()-t0:7.1f}s] {msg}", flush=True)

    # ------------------------------------------------- -1. the A38 gate
    if not rungs_only:
        validate_banked(LAB / "data")
        log("validate_banked: PASS (charter s6.0 gate)")

    # ------------------------------------------------ 0. tower + structure
    LMS = [(18, 12), (18, 6), (9, 6), (9, 3)]
    L = [TowerCode(f"L{i}", lm, red(A_L, lm), red(B_L, lm))
         for i, lm in enumerate(LMS)]
    assert [c.k for c in L] == [12, 12, 8, 8]
    assert [c.n for c in L] == [432, 216, 108, 54]
    deck0 = AxisDeck(L[0], L[1], 1)    # y: 12 -> 6
    deck1 = AxisDeck(L[1], L[2], 0)    # x: 18 -> 9
    deck2 = AxisDeck(L[2], L[3], 1)    # y: 6 -> 3
    for c in L:
        assert not any(int(kv.sum()) % 2 for kv in c.kerHZ), "parity!"
    Mp0, Mp1, Mp2 = h1_map(deck0), h1_map(deck1), h1_map(deck2)
    Mt0 = h1_map(deck0, tau=True)
    SEAMb, _ = rref_ints(list(colspace(Mp0)))
    P2b, _ = rref_ints(list(colspace(Mp2)))
    comp = (Mp1 @ Mp0) % 2
    S1pb, _ = rref_ints(list(colspace(comp)))
    assert len(SEAMb) == 6 and len(P2b) == 4 and len(S1pb) == 2
    assert gf2_rank([v2i(c) for c in Mp1.T]) == 6
    seam_set = span_points(SEAMb) - {0}
    imp1_set = span_points(rref_ints(list(colspace(Mp1)))[0]) - {0}
    S1pset = span_points(S1pb) - {0}
    P2set = span_points(P2b) - {0}
    assert len(seam_set) == 63 and len(imp1_set) == 63 \
        and len(S1pset) == 3 and len(P2set) == 15
    assert S1pset <= imp1_set
    assert span_eq(kernel_ints(Mt0), list(colspace(Mp0))), \
        "top rung exact_base fails?!"
    assert (translation_mat(L[0], deck0.sigma)
            == np.eye(L[0].k, dtype=np.uint8)).all(), "sigma* != id"
    assert deck0.twisted() and deck1.twisted() and deck2.twisted()
    for gi in range(0, L[1].ng, 27):
        tst = (deck0.TAU @ L[1].HX[gi]) % 2
        assert in_span(v2i(tst), L[0].rsHX_b, L[0].rsHX_p), \
            "tau0(stab) not stab: b=0 branch argument breaks"
    log("tower + structure: k 12/12/8/8; SEAM dim 6; rank(p1*p0*) = 2 "
        "=> S1' = 3 classes; top rung (R)+exact+sigma*=id; "
        "tau0(stab)=stab sampled")
    out["structure"] = {"k": [12, 12, 8, 8], "seam_dim": 6,
                        "s1prime_classes": 3}
    perms1 = translation_perms(L[1])
    perms2 = translation_perms(L[2])
    perms3 = translation_perms(L[3])
    binp = cosetbz.build_kernel()

    ckpt_stab = DATA / f"ckpt_W{W}_stab1.jsonl"
    ckpt_seam = DATA / f"ckpt_W{W}_seam1.jsonl"
    ckpt_ntrv = DATA / f"ckpt_W{W}_ntrv1.jsonl"
    ckpt_meta = DATA / f"ckpt_W{W}_meta.json"

    if not rungs_only:
        # ------------------------- 1. L2 primary censuses (direct BZ)
        WG = 14                      # gate window
        hits_s = census_pass(binp, L[2],
                             [("S", np.zeros(L[2].n, np.uint8))], W,
                             "tdg432v2_L2stab")
        s2v = [i2v(h, L[2].n) for h in sorted(hits_s["S"])]
        for v in s2v[:: max(1, len(s2v) // 50)]:
            assert L[2].is_stab(v)
        mu2 = min(int(v.sum()) for v in s2v)
        assert mu2 >= 6
        c2s = Collector(L[2].n)
        for v in s2v:
            c2s.add(v)
        stab2_reps = c2s.reps(perms2)
        log(f"L2 stab <= {W} (direct): {len(s2v)} vectors, "
            f"{len(stab2_reps)} orbit reps {whist(stab2_reps)}, "
            f"mu2 = {mu2}")

        hits_p = census_pass(binp, L[2],
                             [(f"C{c}", rep_for(L[2], c))
                              for c in sorted(S1pset)], W,
                             "tdg432v2_L2s1p")
        c2p = Collector(L[2].n)
        n_s1p = 0
        for c in sorted(S1pset):
            for h in sorted(hits_p[f"C{c}"]):
                v = i2v(h, L[2].n)
                assert L[2].is_cycle(v) and not L[2].is_stab(v)
                assert v2i(L[2].sig(v)) == c
                c2p.add(v)
                n_s1p += 1
        s1p_reps = c2p.reps(perms2)
        log(f"L2 S1'-cosets <= {W} (direct, 3 classes): {n_s1p} "
            f"elements, {len(s1p_reps)} orbit reps {whist(s1p_reps)}")

        c2i = Collector(L[2].n)
        n_imp1 = 0
        imp1_sorted = sorted(imp1_set)
        for lo in range(0, len(imp1_sorted), CH):
            chunk = imp1_sorted[lo:lo + CH]
            hh = census_pass(binp, L[2],
                             [(f"C{c}", rep_for(L[2], c)) for c in chunk],
                             WNT, f"tdg432v2_L2imp1_{lo}")
            for c in chunk:
                for h in sorted(hh[f"C{c}"]):
                    v = i2v(h, L[2].n)
                    assert L[2].is_cycle(v) and not L[2].is_stab(v)
                    c2i.add(v)
                    n_imp1 += 1
        imp1_reps = c2i.reps(perms2)
        log(f"L2 im-p1*-cosets <= {WNT} (direct, 63 classes): "
            f"{n_imp1} elements, {len(imp1_reps)} orbit reps")

        c2a = Collector(L[2].n)       # ALL cycles <= WALL (stab + ntrv)
        d_l2 = None
        alll2 = sorted(set(range(1, 1 << L[2].k)))
        for v in s2v:
            if int(v.sum()) <= WALL:
                c2a.add(v)
        for lo in range(0, len(alll2), CH):
            chunk = alll2[lo:lo + CH]
            hh = census_pass(binp, L[2],
                             [(f"C{c}", rep_for(L[2], c)) for c in chunk],
                             WALL, f"tdg432v2_L2all_{lo}")
            for c in chunk:
                for h in sorted(hh[f"C{c}"]):
                    v = i2v(h, L[2].n)
                    w_ = int(v.sum())
                    d_l2 = w_ if d_l2 is None else min(d_l2, w_)
                    c2a.add(v)
        all2_reps = c2a.reps(perms2)
        assert d_l2 == 10, f"d(L2) = {d_l2} != the banked W16-run 10"
        log(f"L2 all-class <= {WALL} (direct, 255 classes): "
            f"{len(all2_reps)} orbit reps; d(L2) = {d_l2} EXACT")
        out["L2"] = {"stab_orbits": len(stab2_reps),
                     "stab_whist": whist(stab2_reps),
                     "s1p_orbits": len(s1p_reps),
                     "imp1_orbits": len(imp1_reps),
                     "all_orbits_WALL": len(all2_reps),
                     "d_L2": d_l2, "mu2": mu2}

        # ------------------- 2. L3 censuses + the L2 descent GATE <= 14
        hits3 = census_pass(binp, L[3],
                            [("S", np.zeros(L[3].n, np.uint8))], WG,
                            "tdg432v2_L3stab")
        s3v = [i2v(h, L[3].n) for h in sorted(hits3["S"])]
        for v in s3v[:: max(1, len(s3v) // 40)]:
            assert L[3].is_stab(v)
        mu3 = min(int(v.sum()) for v in s3v)
        assert mu3 >= 6
        c3 = Collector(L[3].n)
        for v in s3v:
            c3.add(v)
        stab3_reps = c3.reps(perms3)
        c3p = Collector(L[3].n)
        hits3p = census_pass(
            binp, L[3],
            [(f"C{c}", rep_for(L[3], c)) for c in sorted(P2set)], WG,
            "tdg432v2_L3p2cos")
        for c in sorted(P2set):
            for h in sorted(hits3p[f"C{c}"]):
                v = i2v(h, L[3].n)
                assert L[3].is_cycle(v) and not L[3].is_stab(v)
                c3p.add(v)
        cos3_reps = c3p.reps(perms3)
        # all-class <= 7 for tau3 sources (2w <= 14) + d(L3)
        tau3_src = [v for v in s3v if 2 * int(v.sum()) <= WG]
        d_l3 = None
        allcls3 = sorted(set(range(1, 1 << L[3].k)))
        for lo in range(0, len(allcls3), CH):
            chunk = allcls3[lo:lo + CH]
            hh = census_pass(binp, L[3],
                             [(f"C{c}", rep_for(L[3], c)) for c in chunk],
                             WG, f"tdg432v2_L3all_{lo}")
            for c in chunk:
                for h in sorted(hh[f"C{c}"]):
                    v = i2v(h, L[3].n)
                    w_ = int(v.sum())
                    d_l3 = w_ if d_l3 is None else min(d_l3, w_)
                    if 2 * w_ <= WG:
                        tau3_src.append(v)
        assert d_l3 == 6, f"d(L3) = {d_l3} != the banked W16-run 6"
        log(f"L3 censuses <= {WG}: {len(stab3_reps)} stab + "
            f"{len(cos3_reps)} im-p2* coset orbit reps; d(L3) = {d_l3} "
            f"EXACT; tau3 sources {len(tau3_src)}")
        out["L3"] = {"stab_orbits": len(stab3_reps),
                     "p2cos_orbits": len(cos3_reps), "d_L3": d_l3}

        gate_stab = Collector(L[2].n)
        gate_ntrv = Collector(L[2].n)

        def classify_gate(b2: np.ndarray) -> None:
            w_ = int(b2.sum())
            if w_ == 0 or w_ > WG:
                return
            if in_span(v2i(b2), L[2].rsHX_b, L[2].rsHX_p):
                gate_stab.add(b2)
            else:
                gate_ntrv.add(b2)

        for src in tau3_src:
            b2 = (deck2.TAU @ src) % 2
            assert L[2].is_cycle(b2)
            classify_gate(b2)
        nfg = 0
        for beta in stab3_reps + cos3_reps:
            cap = (WG - int(beta.sum())) // 2
            if cap < 0:
                continue
            lifts = enumerate_lifts_deep(deck2, beta, cap=min(cap, 8))
            nfg += 1
            for v0c, m2 in sorted(lifts.items()):
                classify_gate(deck2.lift(i2v(v0c, L[3].n), beta))
        kdesc_s = key_set(gate_stab.reps(perms2), perms2)
        kdir_s = key_set([b for b in stab2_reps if int(b.sum()) <= WG],
                         perms2)
        assert kdesc_s == kdir_s, \
            f"L2 stab gate FAIL {len(kdesc_s)} != {len(kdir_s)}"
        kdesc_n = key_set(gate_ntrv.reps(perms2), perms2)
        kdir_n = key_set(
            [b for b in all2_reps
             if int(b.sum()) <= WG
             and not in_span(v2i(b), L[2].rsHX_b, L[2].rsHX_p)],
            perms2)
        assert kdesc_n == kdir_n, \
            f"L2 ntrv gate FAIL {len(kdesc_n)} != {len(kdir_n)}"
        log(f"L2 GATE: descent-from-L3 <= {WG} == direct slices EXACTLY "
            f"(stab {len(kdir_s)} orbits, nontrivial {len(kdir_n)} "
            f"orbits; {nfg} fibers) — composed fibers census-complete")
        out["L2_gate"] = {"stab_orbits": len(kdir_s),
                          "ntrv_orbits": len(kdir_n), "equal": True}

        # --------------------------------- 3. L1 obligations (descent)
        ks_fib1 = KernelShift(deck1, all2_reps, complete_to=WALL)
        coll_stab1 = Collector(L[1].n)
        coll_seam1 = Collector(L[1].n)
        coll_ntrv1 = Collector(L[1].n)
        d_l1 = [None]

        def classify_l1(b1: np.ndarray) -> None:
            w_ = int(b1.sum())
            if w_ == 0 or w_ > W:
                return
            if in_span(v2i(b1), L[1].rsHX_b, L[1].rsHX_p):
                coll_stab1.add(b1)
                return
            d_l1[0] = w_ if d_l1[0] is None else min(d_l1[0], w_)
            if v2i(L[1].sig(b1)) in seam_set:
                coll_seam1.add(b1)
            if w_ <= WNT:
                coll_ntrv1.add(b1)

        n_tau2 = 0
        for u in ([v for v in s2v if 2 * int(v.sum()) <= W] +
                  [v for v in c2a.vecs
                   if not in_span(v2i(v), L[2].rsHX_b, L[2].rsHX_p)
                   and 2 * int(v.sum()) <= W]):
            b1 = (deck1.TAU @ u) % 2
            assert L[1].is_cycle(b1)
            classify_l1(b1)
            n_tau2 += 1
        log(f"L1 tau2-family: {n_tau2} sources (all classes <= {W//2})")
        nfd = nfk = 0
        fiber_plan = ([(b, W) for b in stab2_reps] +
                      [(b, W) for b in s1p_reps] +
                      [(b, WNT) for b in imp1_reps])
        for b2, Wb in sorted(fiber_plan, key=lambda t: -int(t[0].sum())):
            wb2 = int(b2.sum())
            cap = (Wb - wb2) // 2
            if cap < 0:
                continue
            if cap <= 6:
                lifts = enumerate_lifts_deep(deck1, b2, cap=cap)
                for v0c, m2 in sorted(lifts.items()):
                    classify_l1(deck1.lift(i2v(v0c, L[2].n), b2))
                nfd += 1
            else:
                v0p, ovp = row_lift_v0(deck1, b2)
                B = wb2 + cap + ovp
                if B <= WALL:
                    rhs = (deck1.RHS @ b2) % 2
                    seen_v0: set[int] = set()
                    bmask = v2i(b2)
                    for v0i in ks_fib1.candidates(b2, v0p, cap):
                        canon = min(v0i, v0i ^ bmask)
                        if canon in seen_v0:
                            continue
                        seen_v0.add(canon)
                        v0 = i2v(v0i, L[2].n)
                        assert not (((deck1.E @ v0) + rhs) % 2).any()
                        classify_l1(deck1.lift(v0, b2))
                    nfk += 1
                else:
                    assert cap <= 8
                    lifts = enumerate_lifts_deep(deck1, b2, cap=cap)
                    for v0c, m2 in sorted(lifts.items()):
                        classify_l1(deck1.lift(i2v(v0c, L[2].n), b2))
                    nfd += 1
            if (nfd + nfk) % 2000 == 0:
                log(f"  ... L2 fibers {nfd+nfk}/{len(fiber_plan)}")
        stab1_reps = coll_stab1.reps(perms1)
        seam1_reps = coll_seam1.reps(perms1)
        ntrv1_reps = coll_ntrv1.reps(perms1)
        mu1 = min(int(b.sum()) for b in stab1_reps)
        assert mu1 >= 6
        assert d_l1[0] == 12, \
            f"d(L1) = {d_l1[0]} != the banked W16-run 12"
        log(f"L1 descent censuses: stab reps <= {W}: {len(stab1_reps)} "
            f"{whist(stab1_reps)} (mu1 {mu1}); SEAM-coset reps <= {W}: "
            f"{len(seam1_reps)} {whist(seam1_reps)}; nontrivial reps "
            f"<= {WNT}: {len(ntrv1_reps)} {whist(ntrv1_reps)}; "
            f"d(L1) = {d_l1[0]} EXACT; fibers {nfd} direct + {nfk} "
            f"kernel-shift")
        out["L1"] = {"stab_orbits": len(stab1_reps),
                     "stab_whist": whist(stab1_reps),
                     "seam_orbits": len(seam1_reps),
                     "seam_whist": whist(seam1_reps),
                     "ntrv_orbits_WNT": len(ntrv1_reps),
                     "d_L1": d_l1[0], "mu1": mu1}

        # 3G. GATE: <= 12 L1-stab census via the INDEPENDENT y-quotient
        WG1 = 12
        LY = TowerCode("LY", (18, 3),
                       fold_support(red(A_L, (18, 6)), 1, 3),
                       fold_support(red(B_L, (18, 6)), 1, 3))
        assert (1 << LY.k) - 1 <= 256, f"k(LY) = {LY.k}"
        decky = AxisDeck(L[1], LY, 1)
        permsy = translation_perms(LY)
        hitsy = census_pass(binp, LY,
                            [("S", np.zeros(LY.n, np.uint8))], WG1,
                            "tdg432v2_LYgS")
        syv = [i2v(h, LY.n) for h in sorted(hitsy["S"])]
        cy = Collector(LY.n)
        for v in syv:
            assert LY.is_stab(v)
            cy.add(v)
        ysrc = [v for v in syv if 2 * int(v.sum()) <= WG1]
        ally = sorted(set(range(1, 1 << LY.k)))
        for lo in range(0, len(ally), CH):
            chunk = ally[lo:lo + CH]
            hh = census_pass(binp, LY,
                             [(f"C{c}", rep_for(LY, c)) for c in chunk],
                             WG1 // 2, f"tdg432v2_LYgC{lo}")
            for c in chunk:
                for h in sorted(hh[f"C{c}"]):
                    v = i2v(h, LY.n)
                    if 2 * int(v.sum()) <= WG1:
                        ysrc.append(v)
        coll_y = Collector(L[1].n)
        for u in ysrc:
            b1 = (decky.TAU @ u) % 2
            assert L[1].is_cycle(b1)
            if 0 < int(b1.sum()) <= WG1 and \
                    in_span(v2i(b1), L[1].rsHX_b, L[1].rsHX_p):
                coll_y.add(b1)
        nyf = 0
        for beta in cy.reps(permsy):
            cap = (WG1 - int(beta.sum())) // 2
            if cap < 0:
                continue
            lifts = enumerate_lifts_deep(decky, beta, cap=min(cap, 8))
            nyf += 1
            for v0c, m2 in sorted(lifts.items()):
                b1 = decky.lift(i2v(v0c, LY.n), beta)
                w_ = int(b1.sum())
                if 0 < w_ <= WG1 and \
                        in_span(v2i(b1), L[1].rsHX_b, L[1].rsHX_p):
                    coll_y.add(b1)
        kdy = key_set(coll_y.reps(perms1), perms1)
        kmx = key_set([b for b in stab1_reps if int(b.sum()) <= WG1],
                      perms1)
        assert kdy == kmx, f"L1 gate FAIL {len(kdy)} != {len(kmx)}"
        log(f"L1 GATE: <= {WG1} L1-stab census re-derived through the "
            f"independent y-quotient (18,3) [[{LY.n},{LY.k}]] "
            f"({nyf} fibers) == x-route slice EXACTLY ({len(kdy)} "
            f"orbits)")
        out["L1_gate"] = {"orbits": len(kdy), "equal": True,
                          "k_LY": LY.k}

        # ------------------------------------------------ checkpoint
        for path, reps in [(ckpt_stab, stab1_reps),
                           (ckpt_seam, seam1_reps),
                           (ckpt_ntrv, ntrv1_reps)]:
            with path.open("w") as f:
                for b in reps:
                    f.write(json.dumps({
                        "w": int(b.sum()),
                        "support": sorted(int(j)
                                          for j in np.nonzero(b)[0]),
                    }) + "\n")
        ckpt_meta.write_text(json.dumps(
            {"d_L1": d_l1[0], "d_L2": d_l2, "d_L3": d_l3,
             "mu1": mu1, "WNT": WNT}, indent=1))
        log(f"checkpoint: {len(stab1_reps)} stab + {len(seam1_reps)} "
            f"seam + {len(ntrv1_reps)} nontrivial orbit reps + meta")
        if census_only:
            out["wall_s"] = round(time.monotonic() - t0, 1)
            (DATA / f"sweep_W{W}_census.json").write_text(
                json.dumps(out, indent=1))
            log(f"census phase complete ({out['wall_s']}s)")
            return
    else:
        meta = json.loads(ckpt_meta.read_text())
        d_l1 = [meta["d_L1"]]
        d_l2 = meta["d_L2"]
        d_l3 = meta["d_L3"]
        mu1 = meta["mu1"]
        assert meta["WNT"] == WNT

        def load(path, checker):
            reps = []
            for line in path.open():
                r = json.loads(line)
                v = np.zeros(L[1].n, dtype=np.uint8)
                v[r["support"]] = 1
                assert int(v.sum()) == r["w"]
                checker(v)
                reps.append(v)
            return reps

        def chk_stab(v):
            assert L[1].is_cycle(v)
            assert in_span(v2i(v), L[1].rsHX_b, L[1].rsHX_p)

        def chk_seam(v):
            assert L[1].is_cycle(v)
            assert not in_span(v2i(v), L[1].rsHX_b, L[1].rsHX_p)
            assert v2i(L[1].sig(v)) in seam_set

        def chk_ntrv(v):
            assert L[1].is_cycle(v)
            assert not in_span(v2i(v), L[1].rsHX_b, L[1].rsHX_p)

        stab1_reps = load(ckpt_stab, chk_stab)
        seam1_reps = load(ckpt_seam, chk_seam)
        ntrv1_reps = load(ckpt_ntrv, chk_ntrv)
        log(f"checkpoint reloaded + re-verified: {len(stab1_reps)} stab "
            f"+ {len(seam1_reps)} seam + {len(ntrv1_reps)} nontrivial; "
            f"d(L1) = {d_l1[0]}, d(L2) = {d_l2}, d(L3) = {d_l3}")

    # ------------------------------------------- 4. the witness ladder
    wit = None
    cand = sorted(seam1_reps + ntrv1_reps, key=lambda v: int(v.sum()))
    for u in cand:
        if v2i(L[1].sig(u)) in seam_set:
            continue
        v = (deck0.TAU @ u) % 2
        assert L[0].is_cycle(v)
        if in_span(v2i(v), L[0].rsHX_b, L[0].rsHX_p):
            continue
        wv = int(v.sum())
        assert wv == 2 * int(u.sum())
        wit = {"w_L0": wv, "w_L1": int(u.sum()),
               "class_L1": int(v2i(L[1].sig(u)))}
        break
    if wit:
        log(f"tau0-witness: nontrivial L0-logical of weight "
            f"{wit['w_L0']} (= 2 x L1 weight {wit['w_L1']}, non-SEAM "
            f"class) VERIFIED end-to-end => d([[432,12]]) <= "
            f"{wit['w_L0']}")
    else:
        log("tau0-witness: none in the window")
    out["tau_witness"] = wit

    # --------------------------------------------------- 5. the rungs
    d_L1_eff = d_l1[0] if d_l1[0] is not None else W + 1
    cell = RungCell("tdg432_top", L[1], L[0], deck0)
    ks_top = KernelShift(deck0, stab1_reps + ntrv1_reps,
                         complete_to=WNT)

    def kernel_shift_rung(b: np.ndarray, M: int) -> dict:
        wb = int(b.sum())
        v0p, ovp = row_lift_v0(deck0, b)
        B = wb + (M - 1) + ovp
        assert B <= WNT, \
            f"top kernel-shift window {B} > WNT (cell |b|={wb} M={M})"
        rhs = (deck0.RHS @ b) % 2
        viols = []
        seen_v0: set[int] = set()
        bmask = v2i(b)
        for v0i in ks_top.candidates(b, v0p, M - 1):
            canon = min(v0i, v0i ^ bmask)
            if canon in seen_v0:
                continue
            seen_v0.add(canon)
            v0 = i2v(v0i, L[1].n)
            assert not (((deck0.E @ v0) + rhs) % 2).any()
            ch = (deck0.EMB[0] @ v0 + deck0.EMB[1] @ ((v0 + b) % 2)) % 2
            if in_span(v2i(ch), L[0].rsHX_b, L[0].rsHX_p):
                continue
            ov = bin(v0i & ~bmask).count("1")
            wt = int(ch.sum())
            assert wt == wb + 2 * ov, "slice identity violated"
            viols.append({"overflow": ov, "weight": wt,
                          "v0_hex": f"{v0i:x}"})
        if viols:
            return {"verdict": "VIOLATION", "lane": "kernel-shift",
                    "w_b": wb, "M": M, "n_viol": len(viols),
                    "violations": viols[:5],
                    "min_weight": min(x["weight"] for x in viols)}
        return {"verdict": "PASS", "lane": "kernel-shift", "w_b": wb,
                "M": M}

    verd: dict[str, int] = {}
    lanes: dict[str, int] = {}
    viol_finds: list[dict] = []
    first3: list[tuple[np.ndarray, int, str]] = []
    n_xval = 0
    n_done = 0
    tR = time.monotonic()
    rung_rows = []
    for b in sorted(stab1_reps, key=lambda b: -int(b.sum())):
        wb = int(b.sum())
        M = (TARGET - wb) // 2
        if M <= 0:
            continue
        if (M - 1) <= 4:
            r = cell.rung(b, M, time.monotonic() + 3600)
            if n_xval < 3:
                v0p, ovp = row_lift_v0(deck0, b)
                if wb + (M - 1) + ovp <= WNT:
                    rks = kernel_shift_rung(b, M)
                    assert rks["verdict"] == r["verdict"], "LANE MISMATCH"
                    n_xval += 1
        else:
            v0p, ovp = row_lift_v0(deck0, b)
            B = wb + (M - 1) + ovp
            if B <= WNT:
                r = kernel_shift_rung(b, M)
                if n_xval < 6 and (M - 1) <= 6:
                    rd = cell.rung(b, M, time.monotonic() + 3600)
                    assert rd["verdict"] == r["verdict"], "LANE MISMATCH"
                    if r["verdict"] == "VIOLATION":
                        assert rd["n_viol"] == r["n_viol"]
                    n_xval += 1
            else:
                assert (M - 1) <= 6, \
                    f"cell |b|={wb} M={M}: no sound lane (B={B})"
                r = cell.rung(b, M, time.monotonic() + 3600)
        verd[r["verdict"]] = verd.get(r["verdict"], 0) + 1
        lanes[r["lane"]] = lanes.get(r["lane"], 0) + 1
        rung_rows.append({"species": "dangerous", "w": wb, "M": M,
                          "verdict": r["verdict"], "lane": r["lane"]})
        if r["verdict"] == "VIOLATION":
            viol_finds.append(r)
        if len(first3) < 3:
            first3.append((b, M, r["verdict"]))
        n_done += 1
        if n_done % 5000 == 0:
            log(f"  ... dangerous rungs {n_done}/{len(stab1_reps)} "
                f"({time.monotonic()-tR:.0f}s)")
    dtR = time.monotonic() - tR
    log(f"dangerous rungs: {sum(verd.values())} at target {TARGET} "
        f"({dtR:.1f}s): {verd}, lanes {lanes}, xval {n_xval}")
    out["dangerous_rungs"] = {"rungs": sum(verd.values()),
                              "verdicts": verd, "lanes": lanes,
                              "wall_s": round(dtR, 1)}

    verd2: dict[str, int] = {}
    for w_el in sorted(seam1_reps, key=lambda v: -int(v.sum())):
        ww = int(w_el.sum())
        M = (TARGET - ww) // 2
        if M <= 0:
            continue
        r = cell.seam_rung(w_el, M)
        verd2[r["verdict"]] = verd2.get(r["verdict"], 0) + 1
        rung_rows.append({"species": "seam", "w": ww, "M": M,
                          "verdict": r["verdict"], "lane": r["lane"]})
        if r["verdict"] == "VIOLATION":
            viol_finds.append(r)
    log(f"seam rungs: {sum(verd2.values())} at target {TARGET}: {verd2}")
    out["seam_rungs"] = {"rungs": sum(verd2.values()), "verdicts": verd2}

    if len(rung_rows) <= 50000:
        with (DATA / f"rungs_W{W}.jsonl").open("w") as f:
            for rr in rung_rows:
                f.write(json.dumps(rr) + "\n")
    else:
        agg: dict[tuple, int] = {}
        for rr in rung_rows:
            key = (rr["species"], rr["w"], rr["M"], rr["verdict"],
                   rr["lane"])
            agg[key] = agg.get(key, 0) + 1
        with (DATA / f"rungs_W{W}.jsonl").open("w") as f:
            for (sp, w_, M_, vd, ln), cnt in sorted(agg.items()):
                f.write(json.dumps({"species": sp, "w": w_, "M": M_,
                                    "verdict": vd, "lane": ln,
                                    "count": cnt}) + "\n")

    # covariance
    g = (5, 7)
    perm_g = perm_for(L[1], g)
    for b, M, v0_verd in first3:
        bt = b[perm_g]
        v0p, ovp = row_lift_v0(deck0, bt)
        B = int(bt.sum()) + (M - 1) + ovp
        rt = (kernel_shift_rung(bt, M) if B <= WNT
              else cell.rung(bt, M, time.monotonic() + 600))
        assert rt["verdict"] == v0_verd, "covariance broken"
    log("covariance: 3 translated reps re-rung, verdicts transport")

    # ------------------------------------------------- 6. the assembly
    b0_dead = 2 * d_L1_eff > W
    if not viol_finds and b0_dead:
        log(f"ASSEMBLY (W = {W}): b = 0 branch dead (|v| = 2|u| >= "
            f"2 d(L1) = {2*d_L1_eff} > {W}); dangerous "
            f"({sum(verd.values())} PASS) + seam ({sum(verd2.values())} "
            f"PASS) + G-transport => NO nontrivial X-logical of weight "
            f"<= {W}: d([[432,12]]) >= {TARGET} at certificate tier "
            f"(consuming d(L1) = {d_l1[0]}, d(L2) = "
            f"{d_l2 if not rungs_only else d_l2}, both measured "
            f"census-complete in-run). Z side by transpose duality.")
        out["verdict"] = {"floor": TARGET, "all_pass": True}
        if wit and wit["w_L0"] == TARGET:
            log(f"WITH the tau0-witness at weight {wit['w_L0']}: "
                f"d([[432,12]]) = {TARGET} EXACT, certificate tier.")
            out["verdict"]["exact"] = TARGET
    elif not viol_finds:
        ub = 2 * d_L1_eff
        log(f"ASSEMBLY (W = {W}): rungs all PASS but the b = 0 branch "
            f"is LIVE: d(L1) = {d_L1_eff} <= {W//2} — the tau0-witness "
            f"(if non-SEAM) gives d <= {ub}. See tau_witness.")
        out["verdict"] = {"rungs_pass": True, "b0_live": True,
                          "d_L1": d_L1_eff}
    else:
        wts = [x["weight"] for r in viol_finds
               for x in r.get("violations", [])]
        wts += [r["min_weight"] for r in viol_finds if "min_weight" in r]
        wmin = min(wts)
        if wit:
            wmin = min(wmin, wit["w_L0"])
        log(f"ASSEMBLY (W = {W}): VIOLATIONS — explicit nontrivial "
            f"L0-logicals found and re-verified, min weight {wmin} "
            f"({len(viol_finds)} cells). With the sweep complete to "
            f"{W}: d([[432,12]]) = {wmin} EXACT if {wmin} <= {W}.")
        out["verdict"] = {"upper": wmin, "all_pass": False,
                          "viol_cells": len(viol_finds)}

    # X<->Z duality spot-check
    ng0 = L[0].ng
    iota = np.zeros(L[0].n, dtype=np.int64)
    for i, e in enumerate(L[0].G):
        j = L[0].G.index(L[0].G.neg(e))
        iota[i] = j
        iota[ng0 + i] = ng0 + j
    swap = np.concatenate([np.arange(ng0, 2 * ng0), np.arange(0, ng0)])
    hzb, hzp = rref_ints([v2i(r) for r in L[0].HZ])
    for kv in L[0].kerHZ[:10]:
        d_ = kv[iota][swap]
        assert not ((L[0].HX @ d_) % 2).any()
    for row in L[0].HX[::100]:
        d_ = row[iota][swap]
        assert in_span(v2i(d_), hzb, hzp)
    log("X<->Z duality spot-check OK => Z side follows")

    out["wall_s"] = round(time.monotonic() - t0, 1)
    (DATA / f"sweep_W{W}.json").write_text(json.dumps(out, indent=1))
    log(f"total {out['wall_s']}s -> {DATA / f'sweep_W{W}.json'}")


if __name__ == "__main__":
    main()
