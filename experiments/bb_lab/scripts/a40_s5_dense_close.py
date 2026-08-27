#!/usr/bin/env python3
"""A40 S5 — the l = 18 periodic-leg completion: complete censuses of
every k > 0 sheared frame (18, p, d), p in {6, 7, 8}, at W = 2p - 1,
by Z2-descent (node-exact, no SAT).

Statement produced per frame: the COMPLETE list of nontrivial X-cycles
of weight <= 2p - 1 of the quotient code on Z^2/<(18,0),(d,p)> — the
sub-rate-2 nontrivial periodic phases of the l = 18 y-walk system with
period p and x-drift d.  Expected (extending the S4 pilot's p <= 5
verdict): EMPTY at every frame; any survivor is classified exactly as
in the pilot (y-spanning / prunable / x-winding / b=1 closure).

Method (the b6_close fiber-complete pattern, frame-generic):
  cover (n = 36p > 192) --fold even SNF axis--> base (n = 18p <= 192,
  kernel-censusable), and for p = 8 a second fold to a bottom at
  n = 72 (depth 2 — the depth-1 kernel bill at kappa = 70, W = 15 is
  ~35 min/frame; the bottom censuses at kappa ~ 35 are seconds).
  Every cover cycle v with |v| <= W has base shadow b = P(v), a base
  cycle with |b| <= |v| and |v| = |b| + 2*overflow, so
    - b = 0: v = tau(u), u a base cycle, |v| = 2|u| (free-Z2-deck
      elementary; tau injectivity asserted per deck) — the tau lane
      over all base cycles <= W//2;
    - b != 0: v in the deep fiber over b with cap (W - |b|)//2
      (enumerate_lifts_deep, node-exact, cap <= 8 asserted).
  Fibers run over base translation-ORBIT REPS: the fold is
  translation-equivariant and cover translations surject onto base
  translations, so any cheap nontrivial cover cycle translates to one
  whose shadow is a rep — emptiness over reps is emptiness outright,
  and survivor lists are complete up to cover translation.  A
  covariance spot-check (lift-count equality on a translated shadow)
  runs per frame.  Every produced vector is re-verified end-to-end
  (is_cycle; weight arithmetic |v| = |b| + 2*off implicitly by cap).

Controls (falsify-first, run in phase "ctrl" BEFORE the new frames):
  A. (18,5,1) k=4 (n = 180, direct kernel census possible): the
     descent census, expanded under cover translations, must equal the
     direct census EXACTLY (as element sets, stabs and nontrivial).
  B. (18,6,0) = the rectangular (18,6) frame up to SNF: its <= 11
     nontrivial census must be EMPTY (banked: d((18,6)) = 12,
     census-complete, S4 note (S)§3.2/§9.6).

Phases (argv[1]): ctrl | p6 | p7 | p8probe | p8rest.  Per-frame
checkpoints under data/a40/s5_dense/; summary s5_dense_<phase>.json.
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
    AxisDeck, TowerCode, census_nodes, enumerate_lifts_deep,
    gf2_rank, i2v, in_span, rep_for, translation_perms, v2i,
    validate_banked,
)

_argv = sys.argv
sys.argv = [_argv[0], "18", "8"]  # a40_s4_prune_pilot reads argv at import
from a38_c37xx_freeze import Collector, census_pass, whist  # noqa: E402
from a40_s4_prune_pilot import (  # noqa: E402
    b1_closure_exists, build_catalog, prunable, x_winds, y_spanning,
)
from a40_s4_phase_triage import quotient_code  # noqa: E402
sys.argv = _argv

DATA = LAB / "data" / "a40"
CKPT = DATA / "s5_dense"
L = 18


def fold_code(code, o, ax, name):
    newo = list(o)
    assert newo[ax] % 2 == 0
    newo[ax] //= 2
    A2 = frozenset((e[0] % newo[0], e[1] % newo[1])
                   for e in code.A.support)
    B2 = frozenset((e[0] % newo[0], e[1] % newo[1])
                   for e in code.B.support)
    return TowerCode(name, tuple(newo), A2, B2), tuple(newo)


def has_info(code):
    try:
        cosetbz.disjoint_info_sets(code.HX)
        return True
    except RuntimeError:
        return False


def kernel_censuses(binp, code, W, tag, log):
    """Complete element censuses <= W: (stab_elements, ntrv_elements)."""
    stabs = []
    hs = census_pass(binp, code, [("S", np.zeros(code.n, np.uint8))],
                     W, f"{tag}S")
    for h in sorted(hs["S"]):
        v = i2v(h, code.n)
        assert code.is_stab(v)
        stabs.append(v)
    ntrv = []
    if code.k:
        allc = sorted(range(1, 1 << code.k))
        CH = 51
        for lo in range(0, len(allc), CH):
            chunk = allc[lo:lo + CH]
            offs = [(f"C{c}", rep_for(code, c)) for c in chunk]
            ha = census_pass(binp, code, offs, W, f"{tag}C{lo}")
            for c in chunk:
                for h in sorted(ha[f"C{c}"]):
                    v = i2v(h, code.n)
                    assert code.is_cycle(v) and not code.is_stab(v)
                    ntrv.append(v)
    log(f"    [{tag}] base censuses <= {W}: {len(stabs)} stabs "
        f"{whist(stabs)}, {len(ntrv)} nontrivial {whist(ntrv)}")
    return stabs, ntrv


def lift_stage(deck, b_stabs, b_ntrv, W, tag, log):
    """All cover cycles <= W, complete up to cover translation.

    Returns dict packed-int -> vector.  tau lane over all given base
    cycle elements with 2|u| <= W; fiber lane over base orbit reps."""
    L0, L1 = deck.cover, deck.base
    rk = gf2_rank([v2i(c) for c in deck.TAU.T])
    assert rk == L1.n, "tau not injective"
    found: dict[int, np.ndarray] = {}

    def take(v):
        w = int(v.sum())
        if 0 < w <= W:
            assert L0.is_cycle(v)
            found[v2i(v)] = v

    n_tau = 0
    for u in b_stabs + b_ntrv:
        if 2 * int(u.sum()) <= W:
            take((deck.TAU @ u) % 2)
            n_tau += 1
    perms1 = translation_perms(L1)
    coll = Collector(L1.n)
    for u in b_stabs + b_ntrv:
        coll.add(u)
    reps = coll.reps(perms1)
    t1 = time.monotonic()
    for i, b in enumerate(sorted(reps, key=lambda v: -int(v.sum()))):
        cap = (W - int(b.sum())) // 2
        if cap < 0:
            continue
        assert cap <= 8
        lifts = enumerate_lifts_deep(deck, b, cap=cap)
        for v0c, m2 in sorted(lifts.items()):
            take(deck.lift(i2v(v0c, L1.n), b))
        if i == 0 and len(perms1) > 1:  # covariance spot-check
            tb = b[perms1[1]]
            assert len(enumerate_lifts_deep(deck, tb, cap=cap)) \
                == len(lifts), "fiber covariance fails"
        if (i + 1) % 50 == 0:
            log(f"    [{tag}] fibers {i+1}/{len(reps)} "
                f"({time.monotonic()-t1:.0f}s)")
    log(f"    [{tag}] lift stage: tau lane {n_tau} in-range base "
        f"cycles, {len(reps)} fiber reps -> {len(found)} cover cycles "
        f"<= {W} (up to cover translation)")
    return found


def classify_frame(cover, p, d, found, W, log):
    """Split found cover cycles; classify nontrivial survivors AND
    verify the local-reduction catalog prunes every TRIVIAL y-spanning
    sub-rate-2 orbit (the pilot's claim, extended to p = 6, 7, 8)."""
    stabs, ntrv = [], []
    for v in found.values():
        (stabs if cover.is_stab(v) else ntrv).append(v)
    log(f"    census by descent <= {W}: {len(stabs)} trivial "
        f"{whist(stabs)}, {len(ntrv)} NONTRIVIAL {whist(ntrv)}")
    cat = build_catalog(cover) if (stabs or ntrv) else []
    triv_yspan = [v for v in stabs if y_spanning(cover, p, d, v)]
    triv_unpruned = [v for v in triv_yspan
                     if not prunable(v2i(v), cat)]
    log(f"    trivial y-spanning sub-rate-2 orbits: {len(triv_yspan)} "
        f"(reps), catalog-UNPRUNED: {len(triv_unpruned)}"
        + (f"  !! catalog gap {whist(triv_unpruned)}"
           if triv_unpruned else " — all pruned"))
    survivors = []
    if ntrv:
        for v in ntrv:
            pts = []
            for i in np.nonzero(v)[0]:
                blk, gi = divmod(int(i), cover.ng)
                pts.append((blk,) + tuple(cover.G.from_index(gi)))
            survivors.append(dict(
                weight=int(v.sum()),
                y_spanning=bool(y_spanning(cover, p, d, v)),
                prunable=bool(prunable(v2i(v), cat)),
                x_winds=bool(x_winds(p, d, pts)),
                b1_closure_r=b1_closure_exists(p, d),
                pts=sorted(pts)))
        for s in survivors:
            log(f"      SURVIVOR w={s['weight']} yspan={s['y_spanning']} "
                f"pruned={s['prunable']} xwinds={s['x_winds']} "
                f"b1_r={s['b1_closure_r']}")
    return stabs, ntrv, survivors, len(triv_yspan), len(triv_unpruned)


def close_frame(binp, p, d, log, depth2=False):
    W = 2 * p - 1
    cover, o = quotient_code(L, p, d)
    name = f"({L},{p},{d})"
    log(f"[{name}] k={cover.k} orders={o} n={cover.n} W={W} "
        f"depth={2 if depth2 else 1}")
    t0 = time.monotonic()

    # choose the first fold: even axis minimizing base census cost
    fopts = []
    for ax in (0, 1):
        if o[ax] % 2 == 0:
            base, bo = fold_code(cover, o, ax, f"{name}b{ax}")
            fopts.append((ax, base, bo))
    assert fopts, f"{name}: no even axis?!"

    if not depth2:
        cands = [(ax, base, bo) for ax, base, bo in fopts
                 if base.n <= 192 and has_info(base)]
        assert cands, f"{name}: no kernel-censusable fold"
        ax, base, bo = min(
            cands, key=lambda t: census_nodes(t[1].kappa, W))
        deck = AxisDeck(cover, base, ax)
        log(f"    fold ax{ax} -> {bo} n={base.n} k={base.k} "
            f"kappa={base.kappa}")
        b_stabs, b_ntrv = kernel_censuses(binp, base, W,
                                          f"d{p}_{d}", log)
        mu = min((int(v.sum()) for v in b_stabs), default=None)
        found = lift_stage(deck, b_stabs, b_ntrv, W, name, log)
        route = dict(depth=1, fold_axis=ax, base_orders=list(bo),
                     base_k=base.k, mu_stab_base=mu,
                     base_ntrv=len(b_ntrv))
    else:
        # depth 2: pick (mid, bottom) minimizing bottom census cost
        best = None
        for ax0, mid, mo in fopts:
            for ax1 in (0, 1):
                if mo[ax1] % 2 == 0:
                    bot, boo = fold_code(mid, mo, ax1, f"{name}bb{ax1}")
                    if bot.n <= 192 and has_info(bot):
                        c = census_nodes(bot.kappa, W)
                        if best is None or c < best[0]:
                            best = (c, ax0, mid, mo, ax1, bot, boo)
        assert best, f"{name}: no depth-2 tower"
        _, ax0, mid, mo, ax1, bot, boo = best
        deck0 = AxisDeck(cover, mid, ax0)
        deck1 = AxisDeck(mid, bot, ax1)
        log(f"    tower ax{ax0} -> {mo} (n={mid.n}, k={mid.k}) "
            f"--ax{ax1}--> {boo} (n={bot.n}, k={bot.k}, "
            f"kappa={bot.kappa})")
        bt_stabs, bt_ntrv = kernel_censuses(binp, bot, W,
                                            f"d{p}_{d}L2", log)
        mid_found = lift_stage(deck1, bt_stabs, bt_ntrv, W,
                               f"{name}/mid", log)
        m_stabs, m_ntrv = [], []
        for v in mid_found.values():
            if in_span(v2i(v), mid.rsHX_b, mid.rsHX_p):
                m_stabs.append(v)
            else:
                m_ntrv.append(v)
        mu = min((int(v.sum()) for v in m_stabs), default=None)
        log(f"    mid censuses by descent <= {W}: {len(m_stabs)} "
            f"stabs {whist(m_stabs)} (mu {mu}), {len(m_ntrv)} "
            f"nontrivial {whist(m_ntrv)} (up to mid translation)")
        found = lift_stage(deck0, m_stabs, m_ntrv, W, name, log)
        route = dict(depth=2, mid_orders=list(mo), bottom_orders=list(boo),
                     mid_k=mid.k, bottom_k=bot.k, mu_stab_mid=mu,
                     mid_ntrv=len(m_ntrv), bottom_ntrv=len(bt_ntrv))

    stabs, ntrv, survivors, n_ty, n_tu = classify_frame(
        cover, p, d, found, W, log)
    wall = round(time.monotonic() - t0, 1)
    log(f"[{name}] {'EMPTY — no nontrivial <= ' + str(W) if not ntrv else str(len(ntrv)) + ' NONTRIVIAL'} ({wall}s)")
    row = dict(p=p, d=d, k=cover.k, orders=list(o), W=W, route=route,
               n_cover_cycles=len(found), n_trivial=len(stabs),
               trivial_whist=whist(stabs), n_nontrivial=len(ntrv),
               nontrivial_whist=whist(ntrv), survivors=survivors,
               trivial_yspan_reps=n_ty, trivial_yspan_unpruned=n_tu,
               wall_s=wall)
    CKPT.mkdir(parents=True, exist_ok=True)
    (CKPT / f"p{p}_d{d}.json").write_text(json.dumps(row, indent=1))
    return row


def expand_translations(found, code):
    perms = translation_perms(code)
    out = set()
    for v in found.values():
        for P in perms:
            out.add(v2i(v[P]))
    return out


def control_A(binp, log):
    """(18,5,d) direct vs descent census equality."""
    p = 5
    for d in (1, 2, 4):
        cover, o = quotient_code(L, p, d)
        if cover.k and has_info(cover) and cover.n <= 192:
            break
    W = 2 * p - 1
    name = f"ctrlA({L},{p},{d})"
    log(f"[{name}] k={cover.k} orders={o} n={cover.n} W={W}")
    dstabs, dntrv = kernel_censuses(binp, cover, W, "cA_dir", log)
    direct = {v2i(v) for v in dstabs} | {v2i(v) for v in dntrv}
    # descent side
    fopts = [(ax,) + fold_code(cover, o, ax, f"cA_b{ax}")
             for ax in (0, 1) if o[ax] % 2 == 0]
    cands = [(ax, b, bo) for ax, b, bo in fopts
             if b.n <= 192 and has_info(b)]
    ax, base, bo = cands[0]
    deck = AxisDeck(cover, base, ax)
    b_stabs, b_ntrv = kernel_censuses(binp, base, W, "cA_desc", log)
    found = lift_stage(deck, b_stabs, b_ntrv, W, name, log)
    desc = expand_translations(found, cover)
    assert desc == direct, \
        f"descent != direct: {len(desc)} vs {len(direct)}"
    log(f"[{name}] PASS: descent census (translation-expanded, "
        f"{len(desc)} elements) == direct census exactly")
    return dict(frame=[L, p, d], n_elements=len(direct), passed=True)


def main():
    t0 = time.monotonic()
    phase = sys.argv[1] if len(sys.argv) > 1 else "ctrl"

    logf = open(DATA / f"s5_dense_{phase}.log", "w")

    def log(m):
        line = f"[{time.monotonic()-t0:8.1f}s] {m}"
        print(line, flush=True)
        logf.write(line + "\n")
        logf.flush()

    validate_banked(LAB / "data")
    log("validate_banked: PASS")
    binp = cosetbz.build_kernel()
    out = {"phase": phase, "frames": []}

    if phase == "ctrl":
        out["control_A"] = control_A(binp, log)
        row = close_frame(binp, 6, 0, log)   # control B
        assert row["n_nontrivial"] == 0, \
            "(18,6,0) must be empty (banked d((18,6)) = 12)"
        log("[ctrlB (18,6,0)] PASS: nontrivial census <= 11 EMPTY — "
            "matches banked d((18,6)) = 12")
        out["frames"].append(row)
    elif phase == "p6":
        for d in (3, 6, 9, 12, 15):
            out["frames"].append(close_frame(binp, 6, d, log))
    elif phase == "p7":
        for d in (1, 2, 4, 5, 7, 8, 10, 11, 13, 14, 16, 17):
            out["frames"].append(close_frame(binp, 7, d, log))
    elif phase == "p8probe":
        # depth-1 repriced GREEN once census_pass's shared-walk chunking
        # is priced right (ctrl phase evidence); depth-2 kept as escape
        out["frames"].append(close_frame(binp, 8, 1, log))
    elif phase == "p8rest":
        for d in (2, 4, 5, 7, 8, 10, 11, 13, 14, 16, 17):
            out["frames"].append(close_frame(binp, 8, d, log))
    else:
        raise SystemExit(f"unknown phase {phase}")

    n_surv = sum(f["n_nontrivial"] for f in out["frames"])
    log(f"phase {phase} done: {len(out['frames'])} frames, "
        f"{n_surv} nontrivial survivors total")
    out["wall_s"] = round(time.monotonic() - t0, 1)
    (DATA / f"s5_dense_{phase}.json").write_text(json.dumps(out, indent=1))
    log(f"wrote {DATA}/s5_dense_{phase}.json ({out['wall_s']} s)")


if __name__ == "__main__":
    main()
