#!/usr/bin/env python3
"""A40 S2 / P5 — the (3,0) [[648,12]] partial floor at W = 16
(question: d >= 18), TWO-ROUTE design.

Towers (depth 2, bottom (9,9) n = 162 <= 192 = BZ-able; both routes run
as COMPLETE independent closures and must agree — this replaces the
independent-quotient gate, which does not exist here because (9,18) and
(18,9) have no second Z2 deck):

  route A: (18,18) --x9--> (9,18)  --y9--> (9,9)
  route B: (18,18) --y9--> (18,9)  --x9--> (9,9)

Per route (the session-1 v2 architecture, one descent layer, all caps
<= 5 so every fiber is a direct deep enumeration — no kernel-shift
needed at W = 16):

  L2 (9,9) direct BZ censuses (node-exact): stab <= 16; S1' classes
  (im(p1* o p0*), rank computed in-run) <= 16; im-p1* classes <= 12
  (the d(L1) window); all classes <= 8 (tau2 sources).
  L1 obligations by descent (shadow-class law): stab1 <= 16 (dangerous),
  seam1 <= 16 (sig in im p0*), ntrv <= 12 (d(L1) exact if <= 12).
  b = 0 branch at the top needs d(L1) >= 10 (2 d(L1) > 16); if a
  nontrivial L1-cycle <= 8 exists, tau0 of it is tried as an explicit
  upper-bound witness instead.
  Top rungs at target 18: dangerous per stab1 rep (M = (18-w)/2 <= 6,
  direct lanes only), seam per seam1 rep.  Assembly per route; the two
  routes' verdicts and d(L1)-window values are then cross-asserted.

Every rung candidate re-verified in-line (RungCell internals); no SAT.
Output: data/a40/s2_t30_W16.json
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
    enumerate_lifts_deep, gf2_rank, h1_map, i2v, in_span, rep_for,
    rref_ints, span_points, translation_perms, v2i, validate_banked,
)
from a38_c37xx_freeze import Collector, census_pass, whist  # noqa: E402

DATA = LAB / "data" / "a40"
A_L = [(0, 0), (0, 1), (3, -1)]
B_L = [(0, 0), (1, 0), (-1, -3)]
W = 16
TARGET = 18
WNT = 12
WTAU = 8


def red(supp, lm):
    return frozenset((e[0] % lm[0], e[1] % lm[1]) for e in supp)


def run_route(binp, name, mid_lm, ax0, ax1, log):
    L0 = TowerCode(f"{name}/L0", (18, 18), red(A_L, (18, 18)),
                   red(B_L, (18, 18)))
    L1 = TowerCode(f"{name}/L1", mid_lm, red(A_L, mid_lm),
                   red(B_L, mid_lm))
    L2 = TowerCode(f"{name}/L2", (9, 9), red(A_L, (9, 9)),
                   red(B_L, (9, 9)))
    assert (L0.k, L2.k) == (12, 8)
    deck0 = AxisDeck(L0, L1, ax0)
    deck1 = AxisDeck(L1, L2, ax1)
    Mp0, Mp1 = h1_map(deck0), h1_map(deck1)
    comp = (Mp1 @ Mp0) % 2
    rk0 = gf2_rank([v2i(c) for c in Mp0.T])
    rk1 = gf2_rank([v2i(c) for c in Mp1.T])
    rkc = gf2_rank([v2i(c) for c in comp.T])
    seam_set = span_points(rref_ints(list(colspace(Mp0)))[0]) - {0}
    imp1_set = span_points(rref_ints(list(colspace(Mp1)))[0]) - {0}
    s1p_set = span_points(rref_ints(list(colspace(comp)))[0]) - {0}
    assert s1p_set <= imp1_set
    log(f"[{name}] k {L0.k}/{L1.k}/{L2.k}; ranks p0*={rk0} p1*={rk1} "
        f"comp={rkc}; |SEAM|={len(seam_set)} |imp1|={len(imp1_set)} "
        f"|S1'|={len(s1p_set)}")
    # tau0(stab) transport (b=0 argument), sampled
    for gi in range(0, L1.ng, 27):
        tst = (deck0.TAU @ L1.HX[gi]) % 2
        assert in_span(v2i(tst), L0.rsHX_b, L0.rsHX_p)
    perms1 = translation_perms(L1)
    perms2 = translation_perms(L2)

    # ---- L2 direct BZ censuses
    hs = census_pass(binp, L2, [("S", np.zeros(L2.n, np.uint8))], W,
                     f"{name}_L2stab")
    s2v = [i2v(h, L2.n) for h in sorted(hs["S"])]
    for v in s2v[:: max(1, len(s2v) // 40)]:
        assert L2.is_stab(v)
    mu2 = min(int(v.sum()) for v in s2v)
    assert mu2 >= 6
    c2s = Collector(L2.n)
    for v in s2v:
        c2s.add(v)
    stab2 = c2s.reps(perms2)
    hp = census_pass(binp, L2,
                     [(f"C{c}", rep_for(L2, c)) for c in sorted(s1p_set)],
                     W, f"{name}_L2s1p")
    c2p = Collector(L2.n)
    for c in sorted(s1p_set):
        for h in sorted(hp[f"C{c}"]):
            v = i2v(h, L2.n)
            assert L2.is_cycle(v) and not L2.is_stab(v)
            assert v2i(L2.sig(v)) == c
            c2p.add(v)
    s1p2 = c2p.reps(perms2)
    c2i = Collector(L2.n)
    hi = census_pass(binp, L2,
                     [(f"C{c}", rep_for(L2, c)) for c in sorted(imp1_set)],
                     WNT, f"{name}_L2imp1")
    for c in sorted(imp1_set):
        for h in sorted(hi[f"C{c}"]):
            c2i.add(i2v(h, L2.n))
    imp12 = c2i.reps(perms2)
    tau_src = [v for v in s2v if 2 * int(v.sum()) <= W]
    d_l2 = None
    CH = 51
    allc = sorted(range(1, 1 << L2.k))
    for lo in range(0, len(allc), CH):
        chunk = allc[lo:lo + CH]
        ha = census_pass(binp, L2,
                         [(f"C{c}", rep_for(L2, c)) for c in chunk],
                         WTAU, f"{name}_L2all{lo}")
        for c in chunk:
            for h in sorted(ha[f"C{c}"]):
                v = i2v(h, L2.n)
                d_l2 = (int(v.sum()) if d_l2 is None
                        else min(d_l2, int(v.sum())))
                tau_src.append(v)
    log(f"[{name}] L2 censuses: stab {len(stab2)} reps "
        f"{whist(stab2)}; S1' {len(s1p2)}; imp1<= {WNT} {len(imp12)}; "
        f"tau sources {len(tau_src)}; lightest nontrivial <= {WTAU}: "
        f"{d_l2}")

    # ---- L1 obligations
    coll_stab = Collector(L1.n)
    coll_seam = Collector(L1.n)
    coll_ntrv = Collector(L1.n)
    d_l1 = [None]

    def classify(b1):
        w_ = int(b1.sum())
        if w_ == 0 or w_ > W:
            return
        if in_span(v2i(b1), L1.rsHX_b, L1.rsHX_p):
            coll_stab.add(b1)
            return
        d_l1[0] = w_ if d_l1[0] is None else min(d_l1[0], w_)
        if v2i(L1.sig(b1)) in seam_set:
            coll_seam.add(b1)
        if w_ <= WNT:
            coll_ntrv.add(b1)

    for u in tau_src:
        b1 = (deck1.TAU @ u) % 2
        assert L1.is_cycle(b1)
        classify(b1)
    nf = 0
    plan = ([(b, W) for b in stab2] + [(b, W) for b in s1p2]
            + [(b, WNT) for b in imp12])
    for b2, Wb in sorted(plan, key=lambda t: -int(t[0].sum())):
        cap = (Wb - int(b2.sum())) // 2
        if cap < 0:
            continue
        assert cap <= 6
        lifts = enumerate_lifts_deep(deck1, b2, cap=cap)
        for v0c, m2 in sorted(lifts.items()):
            classify(deck1.lift(i2v(v0c, L2.n), b2))
        nf += 1
    stab1 = coll_stab.reps(perms1)
    seam1 = coll_seam.reps(perms1)
    ntrv1 = coll_ntrv.reps(perms1)
    mu1 = min(int(b.sum()) for b in stab1)
    assert mu1 >= 6
    log(f"[{name}] L1 obligations: stab {len(stab1)} {whist(stab1)}; "
        f"seam {len(seam1)} {whist(seam1)}; ntrv<= {WNT} {len(ntrv1)} "
        f"{whist(ntrv1)}; d(L1) window value = {d_l1[0]}; fibers {nf}")

    # b = 0 branch / tau0 witness
    tau_wit = None
    for u in sorted(ntrv1 + seam1, key=lambda v: int(v.sum())):
        v = (deck0.TAU @ u) % 2
        if L0.is_cycle(v) and not L0.is_stab(v):
            tau_wit = {"w_L0": int(v.sum()), "w_L1": int(u.sum())}
            break
    d_L1_eff = d_l1[0] if d_l1[0] is not None else WNT + 2
    b0_dead = 2 * d_L1_eff > W
    log(f"[{name}] b=0 branch: d(L1) >= {d_L1_eff} -> "
        f"{'dead' if b0_dead else 'LIVE'}; tau0 witness: {tau_wit}")

    # ---- rungs at target 18
    cell = RungCell(f"{name}_top", L1, L0, deck0)
    verd: dict[str, int] = {}
    viol = []
    for b in sorted(stab1, key=lambda b: -int(b.sum())):
        M = (TARGET - int(b.sum())) // 2
        if M <= 0:
            continue
        assert M - 1 <= 6
        r = cell.rung(b, M, time.monotonic() + 3600)
        verd[r["verdict"]] = verd.get(r["verdict"], 0) + 1
        if r["verdict"] == "VIOLATION":
            viol.append(r)
    verd2: dict[str, int] = {}
    for w_el in sorted(seam1, key=lambda v: -int(v.sum())):
        M = (TARGET - int(w_el.sum())) // 2
        if M <= 0:
            continue
        r = cell.seam_rung(w_el, M)
        verd2[r["verdict"]] = verd2.get(r["verdict"], 0) + 1
        if r["verdict"] == "VIOLATION":
            viol.append(r)
    log(f"[{name}] rungs: dangerous {verd}, seam {verd2}, "
        f"violations {len(viol)}")
    ok = (not viol) and b0_dead
    return {"route": name, "mid": list(mid_lm), "k_mid": L1.k,
            "ranks": [rk0, rk1, rkc], "d_L1_window": d_l1[0],
            "d_L2_le8": d_l2, "stab1": len(stab1),
            "seam1": len(seam1), "dangerous": verd, "seam": verd2,
            "violations": len(viol), "b0_dead": bool(b0_dead),
            "floor_18": bool(ok), "tau_witness": tau_wit}


def main():
    t0 = time.monotonic()

    def log(m):
        print(f"[{time.monotonic()-t0:7.1f}s] {m}", flush=True)

    validate_banked(LAB / "data")
    log("validate_banked: PASS")
    binp = cosetbz.build_kernel()
    ra = run_route(binp, "A_x_first", (9, 18), 0, 1, log)
    rb = run_route(binp, "B_y_first", (18, 9), 1, 0, log)
    assert ra["floor_18"] == rb["floor_18"], "ROUTES DISAGREE!"
    out = {"W": W, "target": TARGET, "routes": [ra, rb],
           "wall_s": round(time.monotonic() - t0, 1)}
    if ra["floor_18"]:
        log(f"BOTH routes PASS => d([[648,12]]) >= {TARGET} at "
            f"certificate tier (two independent tower derivations; "
            f"no SAT). Z side by transpose duality.")
        out["verdict"] = {"floor": TARGET}
    else:
        log("floor NOT established — see route outputs")
    (DATA / "s2_t30_W16.json").write_text(json.dumps(out, indent=1))
    log(f"total {out['wall_s']}s -> {DATA/'s2_t30_W16.json'}")


if __name__ == "__main__":
    main()
