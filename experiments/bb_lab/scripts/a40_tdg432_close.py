#!/usr/bin/env python3
"""A40 P2: the tour-de-gross (2,1) member [[432,12,?]] — the W-sweep.

Tower (the yxy ordering: top rung (R), the A36 shape; P1 lattice +
pricing GREEN at W = 22, `data/a40/pricing.json`):

  L0 (18,12) [[432,12,?]]   d_conj = 24 (arXiv:2506.03094, conjecture)
     |  y-fold (deck sigma = y^6), twisted, (R): k 12->12, exact both,
     |  sigma* = id, rank p0* = 6  => SEAM = im p0* (63 classes)
  L1 (18,6) [[216,12,?]]    d(L1) DISCOVERED here (complete <= W census)
     |  x-fold (x^9), twisted, non-(R): k 12->8, exact_cover only,
     |  rank p1* = 6  => S1 = im p1* (63 of 255 classes)
  L2 (9,6)  [[108,8,?]]     d(L2) DISCOVERED here
     |  y-fold (y^3), twisted, (R): k 8->8, rank p2* = 4 (15 classes)
  L3 (9,3)  [[54,8,?]]      d(L3) DISCOVERED here (direct BZ censuses)

No banked distances exist for L1/L2/L3 — everything is certified in-run,
bottom-up, by COMPLETE weight-<= W census generation:

  every Li-cycle v with |v| <= W has shadow P v in the (complete)
  L(i+1) censuses with class in im p(i+1)* (+ stab), and overflow
  <= (W - |shadow|)/2; the b = 0 branch is the tau-family over the
  complete <= W/2 all-class census one level down.  Fibers enumerate
  the lift sets exhaustively (deep MITM cap <= 8 at n <= 108; the
  kernel-shift lane with ALL-CYCLE windows <= WC at n = 216, sound
  unconditionally because ker E = Z(base) exactly and the window census
  is complete to WC).

Questions (TARGET = W + 2; all cycles even — parity scope):
  W = 16: all rungs PASS  =>  d([[432,12]]) >= 18 (matches two-gross).
  W = 22: all rungs PASS  =>  d >= 24 = the conjectured value; any
          violation/tau-find is an explicit verified logical => d <= w
          (and the complete sweep makes the minimum found EXACT).

The witness side (either W): m_ns := min weight of a nontrivial
L1-cycle <= W with class OUTSIDE SEAM gives the explicit L0-logical
tau0(u) of weight 2 m_ns (exact_base at the top rung: ker tau0* =
im p0* = SEAM, so non-SEAM classes lift to NONTRIVIAL tau-images) =>
d <= 2 m_ns.  d(L1) = 12 with a non-SEAM minimizer + a PASS sweep at
W = 22 would pin d = 24 EXACT.

Gates (falsify-first, charter A38 s6.0): validate_banked() first; the
L2 descent censuses == direct BZ censuses <= 14 (stab + all 255
classes, orbit key sets); the <= 12 L1-stab census re-derived through
the INDEPENDENT y-quotient (18,3); every census pass carries the exact
node-count assert (cosetbz); every rung candidate re-verified in-line
(E-system, membership, slice identity); covariance spot-checks;
X<->Z duality spot-check.

Usage: python a40_tdg432_close.py 16                 (shakedown, 1 shot)
       python a40_tdg432_close.py 22 --census-only   (phase 1)
       python a40_tdg432_close.py 22 --rungs-only    (phase 2)

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
    enumerate_lifts_deep, fold_support, h1_map, i2v, in_span,
    kernel_ints, perm_for, rep_for, rref_ints, span_points,
    translation_perms, v2i, validate_banked,
)
from a38_c37xx_freeze import (  # noqa: E402
    Collector, KernelShift, census_pass, row_lift_v0, whist,
)

DATA = LAB / "data" / "a40" / "tdg432"
DATA.mkdir(parents=True, exist_ok=True)

A_L = [(0, 0), (0, 1), (3, -1)]
B_L = [(0, 0), (1, 0), (-1, -3)]


def red(supp, lm):
    return frozenset((e[0] % lm[0], e[1] % lm[1]) for e in supp)


def main() -> None:
    W = int(sys.argv[1]) if len(sys.argv) > 1 else 16
    assert W in (16, 18, 22)
    census_only = "--census-only" in sys.argv
    rungs_only = "--rungs-only" in sys.argv
    assert not (census_only and rungs_only)
    TARGET = W + 2
    WC = W - 4               # all-class window bound (>= 2*ceil(W/4))
    t0 = time.monotonic()
    out: dict = {"W": W, "target": TARGET, "WC": WC,
                 "phase": ("census" if census_only else
                           "rungs" if rungs_only else "full")}

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
    # screen-measured structure, hard-asserted (data/a40/pricing.json)
    Mp0, Mp1, Mp2 = h1_map(deck0), h1_map(deck1), h1_map(deck2)
    Mt0 = h1_map(deck0, tau=True)
    SEAMb, _ = rref_ints(list(colspace(Mp0)))
    S1b, _ = rref_ints(list(colspace(Mp1)))
    P2b, _ = rref_ints(list(colspace(Mp2)))
    assert len(SEAMb) == 6 and len(S1b) == 6 and len(P2b) == 4
    seam_set = span_points(SEAMb) - {0}
    S1set = span_points(S1b) - {0}
    P2set = span_points(P2b) - {0}
    assert len(seam_set) == 63 and len(S1set) == 63 and len(P2set) == 15
    # top rung (R) + exactness (ker tau0* = im p0*): the witness machine
    kt0b, _ = rref_ints(list(kernel_ints(Mt0)))
    from bb_lab.tower import span_eq
    assert span_eq(kernel_ints(Mt0), list(colspace(Mp0))), \
        "top rung exact_base fails?!"
    St = np.eye(L[0].k, dtype=np.uint8)
    from bb_lab.tower import translation_mat
    assert (translation_mat(L[0], deck0.sigma) == St).all(), "sigma*!=id"
    assert deck0.twisted() and deck1.twisted() and deck2.twisted()
    # tau0(L1-stab) is an L0-STAB (b=0 branch transport), sampled
    for gi in range(0, L[1].ng, 27):
        tst = (deck0.TAU @ L[1].HX[gi]) % 2
        assert in_span(v2i(tst), L[0].rsHX_b, L[0].rsHX_p), \
            "tau0(stab) not stab: b=0 branch argument breaks"
    log("tower + structure: k 12/12/8/8; SEAM dim 6 / S1 dim 6 / "
        "im p2* dim 4; top rung (R)+exact+sigma*=id; tau0(stab)=stab "
        "sampled — screen values reproduced, b=0 transport verified")
    out["structure"] = {"k": [12, 12, 8, 8], "seam_dim": 6, "s1_dim": 6}
    perms1 = translation_perms(L[1])
    perms2 = translation_perms(L[2])
    perms3 = translation_perms(L[3])
    binp = cosetbz.build_kernel()

    ckpt_stab = DATA / f"ckpt_W{W}_stab1.jsonl"
    ckpt_seam = DATA / f"ckpt_W{W}_seam1.jsonl"
    ckpt_ntrv = DATA / f"ckpt_W{W}_ntrv1.jsonl"
    ckpt_meta = DATA / f"ckpt_W{W}_meta.json"

    if not rungs_only:
        # --------------------------------------------- 1. L3 (direct BZ)
        hits3 = census_pass(binp, L[3],
                            [("S", np.zeros(L[3].n, np.uint8))], W,
                            "tdg432_L3stab")
        s3v = [i2v(h, L[3].n) for h in sorted(hits3["S"])]
        for v in s3v[:: max(1, len(s3v) // 40)]:
            assert L[3].is_stab(v)
        mu3 = min(int(v.sum()) for v in s3v)
        assert mu3 >= 6, f"mu3 = {mu3} < 6"
        c3 = Collector(L[3].n)
        for v in s3v:
            c3.add(v)
        stab3_reps = c3.reps(perms3)
        log(f"L3 stab <= {W}: {len(s3v)} vectors, "
            f"{len(stab3_reps)} orbit reps, mu3 = {mu3}")
        # im p2* cosets <= W (shadow classes of L2-cycles)
        hits3p = census_pass(
            binp, L[3],
            [(f"C{c}", rep_for(L[3], c)) for c in sorted(P2set)], W,
            "tdg432_L3p2cos")
        cos3: dict[int, list[np.ndarray]] = {}
        for c in sorted(P2set):
            vs = [i2v(h, L[3].n) for h in sorted(hits3p[f"C{c}"])]
            for v in vs[:: max(1, len(vs) // 10)]:
                assert L[3].is_cycle(v) and not L[3].is_stab(v)
                assert v2i(L[3].sig(v)) == c
            cos3[c] = vs
        n_cos3 = sum(len(v) for v in cos3.values())
        # all-class <= W/2 (tau3 sources) + d(L3)
        WT3 = W // 2
        allcls = sorted(set(range(1, 1 << L[3].k)))
        d_l3 = None
        tau3_src: list[np.ndarray] = [v for v in s3v
                                      if 2 * int(v.sum()) <= W]
        CH = 51                     # offsets per pass (kernel cap 256)
        for lo in range(0, len(allcls), CH):
            chunk = allcls[lo:lo + CH]
            hh = census_pass(
                binp, L[3],
                [(f"C{c}", rep_for(L[3], c)) for c in chunk],
                max(WT3, 14), f"tdg432_L3all{lo}")
            for c in chunk:
                for h in sorted(hh[f"C{c}"]):
                    v = i2v(h, L[3].n)
                    w = int(v.sum())
                    d_l3 = w if d_l3 is None else min(d_l3, w)
                    if 2 * w <= W:
                        tau3_src.append(v)
        log(f"L3: im-p2* coset elements <= {W}: {n_cos3} over 15 "
            f"classes; d(L3) = {d_l3} EXACT (all 255 classes censused "
            f"<= {max(WT3, 14)}); tau3 sources {len(tau3_src)}")
        out["L3"] = {"stab_vecs": len(s3v), "stab_orbits": len(stab3_reps),
                     "p2cos_elems": n_cos3, "d_L3": d_l3}

        # ------------------------------------- 2. L2 (descent from L3)
        coll_stab2 = Collector(L[2].n)
        coll_s1cos2 = Collector(L[2].n)
        coll_all2 = Collector(L[2].n)     # all cycles <= WC (windows)
        d_l2 = [None]

        def classify_l2(b2: np.ndarray) -> None:
            w = int(b2.sum())
            if w == 0 or w > W:
                return
            if in_span(v2i(b2), L[2].rsHX_b, L[2].rsHX_p):
                if w <= W:
                    coll_stab2.add(b2)
                if w <= WC:
                    coll_all2.add(b2)
                return
            d_l2[0] = w if d_l2[0] is None else min(d_l2[0], w)
            if v2i(L[2].sig(b2)) in S1set:
                coll_s1cos2.add(b2)
            if w <= WC:
                coll_all2.add(b2)

        n_tau3 = 0
        for src in tau3_src:
            b2 = (deck2.TAU @ src) % 2
            assert L[2].is_cycle(b2)
            classify_l2(b2)
            n_tau3 += 1
        nf = nl = 0
        fib_srcs = ([(v, "stab") for v in stab3_reps] +
                    [(Collector(L[3].n), c) for c in []])
        # orbit reps for the 15 coset classes
        cos3_reps = []
        for c in sorted(P2set):
            cc = Collector(L[3].n)
            for v in cos3[c]:
                cc.add(v)
            cos3_reps.extend(cc.reps(perms3))
        for beta in stab3_reps + cos3_reps:
            cap = (W - int(beta.sum())) // 2
            if cap < 0:
                continue
            lifts = enumerate_lifts_deep(deck2, beta, cap=min(cap, 8))
            nf += 1
            for v0c, m2 in sorted(lifts.items()):
                classify_l2(deck2.lift(i2v(v0c, L[3].n), beta))
                nl += 1
            if nf % 500 == 0:
                log(f"  ... L3 fibers {nf}/{len(stab3_reps)+len(cos3_reps)}"
                    f" ({nl} lifts)")
        stab2_reps = coll_stab2.reps(perms2)
        s1cos2_reps = coll_s1cos2.reps(perms2)
        all2_reps = coll_all2.reps(perms2)
        mu2 = min(int(b.sum()) for b in stab2_reps)
        assert mu2 >= 6
        log(f"L2 descent: stab reps <= {W}: {len(stab2_reps)} "
            f"{whist(stab2_reps)} (mu2 {mu2}); S1-coset reps <= {W}: "
            f"{len(s1cos2_reps)} {whist(s1cos2_reps)}; all-cycle reps "
            f"<= {WC}: {len(all2_reps)}; d(L2) = {d_l2[0]} EXACT "
            f"(complete <= {W}); fibers {nf}, lifts {nl}, tau3 {n_tau3}")
        out["L2"] = {"stab_orbits": len(stab2_reps),
                     "stab_whist": whist(stab2_reps),
                     "s1_orbits": len(s1cos2_reps),
                     "all_orbits": len(all2_reps), "d_L2": d_l2[0],
                     "mu2": mu2}

        # 2G. GATE: direct BZ censuses at L2 <= 14 == the descent slice
        WG = 14
        dg_stab = census_pass(binp, L[2],
                              [("S", np.zeros(L[2].n, np.uint8))], WG,
                              "tdg432_L2gS")
        gs = [i2v(h, L[2].n) for h in sorted(dg_stab["S"])]
        cg = Collector(L[2].n)
        for v in gs:
            assert L[2].is_stab(v)
            cg.add(v)
        kd = {bytes(k) for k in batch_keys(
            np.array(cg.vecs, np.uint8), perms2)} if cg.vecs else set()
        m14 = [b for b in stab2_reps if int(b.sum()) <= WG]
        km = {bytes(k) for k in batch_keys(
            np.array(m14, np.uint8), perms2)} if m14 else set()
        assert kd == km, f"L2 stab gate FAIL {len(kd)} != {len(km)}"
        # all 255 classes <= 14, chunked
        cg2 = Collector(L[2].n)
        n_dir_c = 0
        alll2 = sorted(set(range(1, 1 << L[2].k)))
        for lo in range(0, len(alll2), CH):
            chunk = alll2[lo:lo + CH]
            hh = census_pass(binp, L[2],
                             [(f"C{c}", rep_for(L[2], c)) for c in chunk],
                             WG, f"tdg432_L2gC{lo}")
            for c in chunk:
                for h in sorted(hh[f"C{c}"]):
                    v = i2v(h, L[2].n)
                    assert L[2].is_cycle(v) and not L[2].is_stab(v)
                    cg2.add(v)
                    n_dir_c += 1
        # nontrivial comparison at the honest common window: the descent
        # all-cycle collection stops at WC, so compare <= min(WG, WC)
        WGC = min(WG, WC)
        cgc = [v for v in cg2.vecs if int(v.sum()) <= WGC]
        kd2 = {bytes(k) for k in batch_keys(
            np.array(cgc, np.uint8), perms2)} if cgc else set()
        m14c = [b for b in all2_reps
                if int(b.sum()) <= WGC
                and not in_span(v2i(b), L[2].rsHX_b, L[2].rsHX_p)]
        km2 = {bytes(k) for k in batch_keys(
            np.array(m14c, np.uint8), perms2)} if m14c else set()
        assert kd2 == km2, \
            f"L2 all-class gate FAIL at <= {WGC}: " \
            f"{len(kd2)} != {len(km2)}"
        # d(L2) cross-check from the direct pass (complete to WG >= WC)
        if n_dir_c:
            dmin_dir = min(int(v.sum()) for v in cg2.vecs)
            assert dmin_dir == d_l2[0], (dmin_dir, d_l2[0])
        log(f"L2 GATE: direct BZ <= {WG} == descent slice EXACTLY "
            f"(stab {len(kd)} orbits; nontrivial {len(kd2)} orbits, "
            f"d(L2) cross-checked) — composed fibers census-complete")
        out["L2_gate"] = {"stab_orbits": len(kd),
                          "ntrv_orbits": len(kd2), "equal": True}

        # --------------------------------- 3. L1 obligations (descent)
        ks_fib1 = KernelShift(deck1, all2_reps, complete_to=WC)
        coll_stab1 = Collector(L[1].n)
        coll_seam1 = Collector(L[1].n)
        coll_ntrv1 = Collector(L[1].n)    # nontrivial <= WC (top windows)
        d_l1 = [None]

        def classify_l1(b1: np.ndarray) -> None:
            w = int(b1.sum())
            if w == 0 or w > W:
                return
            if in_span(v2i(b1), L[1].rsHX_b, L[1].rsHX_p):
                coll_stab1.add(b1)
                return
            d_l1[0] = w if d_l1[0] is None else min(d_l1[0], w)
            if v2i(L[1].sig(b1)) in seam_set:
                coll_seam1.add(b1)
            if w <= WC:
                coll_ntrv1.add(b1)

        n_tau2 = 0
        for w_, bucket in ks_fib1.by_w.items():
            if 2 * w_ <= W:
                for zi in bucket:
                    b1 = (deck1.TAU @ i2v(zi, L[2].n)) % 2
                    assert L[1].is_cycle(b1)
                    classify_l1(b1)
                    n_tau2 += 1
        log(f"L1 tau2-family: {n_tau2} sources")
        nfd = nfk = 0
        for b2 in sorted(stab2_reps + s1cos2_reps,
                         key=lambda b: -int(b.sum())):
            wb2 = int(b2.sum())
            cap = (W - wb2) // 2
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
                assert B <= WC, \
                    f"L1-fiber kernel-shift window {B} > WC = {WC}"
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
            if (nfd + nfk) % 1000 == 0:
                log(f"  ... L2 fibers {nfd+nfk}/"
                    f"{len(stab2_reps)+len(s1cos2_reps)}")
        stab1_reps = coll_stab1.reps(perms1)
        seam1_reps = coll_seam1.reps(perms1)
        ntrv1_reps = coll_ntrv1.reps(perms1)
        mu1 = min(int(b.sum()) for b in stab1_reps)
        assert mu1 >= 6
        log(f"L1 descent censuses <= {W}: stab reps {len(stab1_reps)} "
            f"{whist(stab1_reps)} (mu1 {mu1}); SEAM-coset reps "
            f"{len(seam1_reps)} {whist(seam1_reps)}; nontrivial reps "
            f"<= {WC}: {len(ntrv1_reps)}; d(L1) = "
            f"{d_l1[0] if d_l1[0] is not None else f'> {W}'} "
            f"(complete <= {W}); fibers {nfd} direct + {nfk} "
            f"kernel-shift")
        out["L1"] = {"stab_orbits": len(stab1_reps),
                     "stab_whist": whist(stab1_reps),
                     "seam_orbits": len(seam1_reps),
                     "seam_whist": whist(seam1_reps),
                     "ntrv_orbits_WC": len(ntrv1_reps),
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
                            "tdg432_LYgS")
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
                             WG1 // 2, f"tdg432_LYgC{lo}")
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
        yreps = coll_y.reps(perms1)
        kdy = {bytes(k) for k in batch_keys(
            np.array(yreps, np.uint8), perms1)} if yreps else set()
        m12 = [b for b in stab1_reps if int(b.sum()) <= WG1]
        kmx = {bytes(k) for k in batch_keys(
            np.array(m12, np.uint8), perms1)} if m12 else set()
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
            {"d_L1": d_l1[0], "d_L2": d_l2[0], "d_L3": d_l3,
             "mu1": mu1}, indent=1))
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
        d_l2 = [meta["d_L2"]]
        d_l3 = meta["d_L3"]
        mu1 = meta["mu1"]

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
            f"d(L1) = {d_l1[0]}, d(L2) = {d_l2[0]}, d(L3) = {d_l3}")

    # ------------------------------------------- 4. the witness ladder
    # m_ns = min weight of a nontrivial L1-cycle <= W with class outside
    # SEAM: tau0(u) is then a verified NONTRIVIAL L0-logical, wt 2|u|.
    wit = None
    cand = sorted(seam1_reps + ntrv1_reps, key=lambda v: int(v.sum()))
    seen_w = set()
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
        log(f"tau0-witness: no non-SEAM nontrivial L1-cycle <= "
            f"{WC if not rungs_only else WC} found -> no tau upper "
            f"bound from this window")
    out["tau_witness"] = wit

    # --------------------------------------------------- 5. the rungs
    d_L1_eff = d_l1[0] if d_l1[0] is not None else W + 1
    cell = RungCell("tdg432_top", L[1], L[0], deck0)
    ks_top = KernelShift(deck0, stab1_reps + ntrv1_reps,
                         complete_to=WC)

    def kernel_shift_rung(b: np.ndarray, M: int) -> dict:
        wb = int(b.sum())
        v0p, ovp = row_lift_v0(deck0, b)
        B = wb + (M - 1) + ovp
        assert B <= WC, \
            f"top kernel-shift window {B} > WC (cell |b|={wb} M={M})"
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
                if wb + (M - 1) + ovp <= WC:
                    rks = kernel_shift_rung(b, M)
                    assert rks["verdict"] == r["verdict"], "LANE MISMATCH"
                    n_xval += 1
        else:
            v0p, ovp = row_lift_v0(deck0, b)
            B = wb + (M - 1) + ovp
            if B <= WC:
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
        rt = (kernel_shift_rung(bt, M) if B <= WC
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
            f"(consuming d(L1) = {d_l1[0]}, d(L2) = {d_l2[0]}, both "
            f"measured census-complete in-run). Z side by transpose "
            f"duality.")
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
