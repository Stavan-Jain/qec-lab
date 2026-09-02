#!/usr/bin/env python3
"""A40 S10 — Stage 1: the FINANCED HALF of the mu-echo lemma, as a
theorem.  R0 (§14.2) proved that zero-input dynamics never raises
the pair-min mu_t := min supp(v2[t] u v2[t+1]).  This module extends
R0's bookkeeping to count the K1 input interventions.

THE STEP ALPHABET (y-lane; solved recurrence
    v2[t+1] = (1+x^-1) v1[t] + x v1[t-3] + v2[t] + x^-3 v2[t-1]).
Read E at (t+2, c-3) with c := m_t := min supp v2[t]:
    v2[t+2](c-3) = v1[t+1](c-3) + v1[t+1](c-2) + v1[t-2](c-4)
                   + v2[t+1](c-3) + v2[t](c),
and v2[t](c) = 1.  Hence

LEMMA F (echo-charge; the financed half).  For every row t with
v2[t] nonempty, c = m_t:  EITHER  mu_{t+1} <= c - 3  (the echo fires
at (t+2, c-3), or v2[t+1] already holds a cell there),  OR an ODD
number of the three K1 cells  (t+1, c-3), (t+1, c-2), (t-2, c-4)
is occupied.  Proof: if the kill parity is even then
v2[t+2](c-3) = v2[t+1](c-3) + 1, so one of v2[t+1], v2[t+2] holds
c-3, and mu_{t+1} = min(v2[t+1] u v2[t+2]) <= c - 3.  QED.

The three kill cells lie at columns <= c-2, strictly LEFT of the
current v2-minimum: a step at which mu fails to drop by 3 below m_t
is CHARGED to an input cell left of the frontier.  Consequences:
  (F1) zero-input: mu_{t+1} <= m_t - 3 whenever v2[t] != 0, and
       mu is non-increasing — R0, with the drop RATE made explicit
       (>= 3 columns per two rows: 1.5 columns/row of free
       leftward drift).
  (F2, holding) over T consecutive steps in which mu never drops
       below its starting value, no two consecutive steps are free,
       so >= floor(T/2) steps are charged; a K1 cell (tau, c') can
       charge at most two steps (t = tau-1 via m_{tau-1} in
       {c'+2, c'+3}; t = tau+2 via m_{tau+2} = c'+4), so the hold
       costs >= ceil(T/4) input cells, all placed LEFT of the v2
       frontier.
  (F3, the window anchor) a charging cell at (t+1, <= c-2) sits in
       the 4-row windows of slabs t+1..t+4, so the WINDOW anchor
       (v1 u v2) cannot exceed c-2 before slab t+5: a K1-financed
       rise of the pair-min is invisible to the slab anchor for four
       slabs, and visible only after the killer's own descendants
       ((1+x^-1) debris at c-3/c-2, the x-echo at c-1 four rows up)
       are dealt with in turn — the recursion S9 called the K1
       diagonal.

WHAT IS NOT TRUE (the naive statement, re-scoped).  "Net anchor rise
<= 3 + #K1 cells in the stretch" is FALSE as a local lemma: once a
right strand exists, ONE echo-kill lets the pair-min jump onto it,
and the jump size is the strand's distance, not the input count.
The `naive` lane exhibits a verified zero-input fragment (rise 8 with
NO input cell in the stretch) and the smallest such windows.  The
census wall "+3" (§14.1) is therefore a COST statement (the right
strand and the kills must be paid for from the seed) — consistent
with F, not implied by a per-cell count.  From a weight-1 seed the
max-gauge cap gives the only local per-cell bound: A_t - A_seed
<= 4 * #K1 + 1 (dil-4 growth scope).

LANES.  `lemma`: the identity behind F verified EXHAUSTIVELY over
the four-row step alphabet (v1[t-2], v1[t+1], v2[t], v2[t+1]) at
width 6 (2^24 states) and on 10^6 random width-14 states; the F2
two-step drop verified on zero-input states; the mirror lane's
analogue.  `naive`: witness search.  `verify`: F checked row by row
on the S9 slip -8 witness, on the banked species lifts (W7, TC63 at
both ell), and on every completed crossing replayed from a parented
link march below a small g-cap (zero violations required)."""
from __future__ import annotations

import argparse
import json
import random
import sys
import time
from pathlib import Path

LAB = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(LAB / "src"))
sys.path.insert(0, str(LAB / "scripts"))
DATA = LAB / "data" / "a40"

from a40_s6_frontier import wt, lsb, tooth_ok, seeds_full  # noqa: E402

PRE, HEAVY, POST = 0, 1, 2


def v2_forced(v1_tm3, v1_t, v2_tm1, v2_t):
    """v2[t+1] from the y-lane recurrence (bitmask rows, columns
    shifted so no >> underflows)."""
    return v1_t ^ (v1_t >> 1) ^ (v1_tm3 << 1) ^ v2_t ^ (v2_tm1 >> 3)


def min_col(mask):
    return lsb(mask) if mask else None


# ---------------------------------------------------------------------
# F on a cover fragment (rows as (set v1, set v2)), the reference
# implementation used by every verification lane
# ---------------------------------------------------------------------

def f_steps(rows, t0):
    """rows[j - t0] = (set v1[j], set v2[j]).  For every t with
    v2[t] nonempty and rows t-2..t+2 inside the fragment, classify
    the step: drop (mu_{t+1} <= m_t - 3), charged (odd kill parity),
    and flag a VIOLATION if neither.  Returns per-step records."""
    n = len(rows)

    def s(blk, j):
        if j < t0 or j >= t0 + n:
            return frozenset()
        return rows[j - t0][blk]
    recs = []
    for t in range(t0 + 2, t0 + n - 2):
        v2t = s(1, t)
        if not v2t:
            continue
        c = min(v2t)
        mu_next_cells = s(1, t + 1) | s(1, t + 2)
        mu_next = min(mu_next_cells) if mu_next_cells else None
        drop = mu_next is not None and mu_next <= c - 3
        kill = ((c - 3) in s(0, t + 1)) + ((c - 2) in s(0, t + 1)) \
            + ((c - 4) in s(0, t - 2))
        charged = (kill % 2 == 1)
        # the identity itself: v2[t+2](c-3) == v2[t+1](c-3) + 1 + kill
        lhs = (c - 3) in s(1, t + 2)
        rhs = (((c - 3) in s(1, t + 1)) + 1 + kill) % 2 == 1
        recs.append(dict(t=t, m=c, mu_next=mu_next, drop=drop,
                         kill=kill, charged=charged,
                         identity_ok=(lhs == rhs),
                         violation=(not drop and not charged)))
    return recs


# ---------------------------------------------------------------------
# lane: lemma
# ---------------------------------------------------------------------

def lane_lemma(out):
    import numpy as np
    t0 = time.time()
    W = 6
    OFF = 8                       # frame offset: columns live in [8, 14)
    n_states = 0
    n_case = 0
    viol_identity = 0
    viol_F = 0
    n_two_step = 0
    viol_two_step = 0
    # exhaustive over (v1[t-2], v1[t+1], v2[t+1]) as numpy vectors
    # per value of v2[t]
    A = np.arange(1 << W, dtype=np.int64)
    v1m2 = (A[:, None, None] << OFF)
    v1p1 = (A[None, :, None] << OFF)
    v2p1 = (A[None, None, :] << OFF)
    for r0 in range(1, 1 << W):
        v2t = np.int64(r0 << OFF)
        c = lsb(int(v2t))
        v2p2 = v1p1 ^ (v1p1 >> 1) ^ (v1m2 << 1) ^ v2p1 ^ (v2t >> 3)
        kill = ((v1p1 >> (c - 3)) & 1) ^ ((v1p1 >> (c - 2)) & 1) \
            ^ ((v1m2 >> (c - 4)) & 1)
        lhs = (v2p2 >> (c - 3)) & 1
        rhs = ((v2p1 >> (c - 3)) & 1) ^ 1 ^ kill
        bad = np.count_nonzero(lhs != rhs)
        viol_identity += int(bad)
        # F: kill even  =>  min(v2[t+1] u v2[t+2]) <= c-3
        union = v2p1 | v2p2
        low = union & ((1 << (c - 2)) - 1)     # cells at columns <= c-3
        f_bad = np.count_nonzero((kill == 0) & (low == 0))
        viol_F += int(f_bad)
        n_states += union.size
        n_case += int(np.count_nonzero(kill == 0))
        # F2 two-step drop in ZERO-input dynamics: mu_{t+2} <= mu_t - 3
        # (v1 rows all zero: slice [0, 0, :])
        u0 = union[0, 0, :]
        mu_t = np.minimum.reduce([np.full_like(u0, c),
                                  np.where(v2p1[0, 0, :] > 0,
                                           (v2p1[0, 0, :]
                                            & -v2p1[0, 0, :]), 1 << 62)])
        # (min column of a mask m is log2(m & -m); compare masks by
        # low bit instead: compute mu via explicit loop — small)
        for r1 in range(1 << W):
            b = r1 << OFF
            a = int(v2t)
            mu0 = lsb(a | b)
            cc = b ^ (a >> 3)
            dd = cc ^ (b >> 3)
            if not (cc | dd):
                continue
            mu2 = lsb(cc | dd)
            n_two_step += 1
            if mu2 > mu0 - 3:
                viol_two_step += 1
    assert viol_identity == 0, viol_identity
    assert viol_F == 0, viol_F
    assert viol_two_step == 0, viol_two_step
    out["F_exhaustive"] = dict(
        width=W, states=n_states, even_kill_states=n_case,
        identity_violations=viol_identity, F_violations=viol_F,
        zero_input_two_step_pairs=n_two_step,
        two_step_violations=viol_two_step)
    print(f"LEMMA F exhaustive (width {W}): {n_states} four-row states,"
          f" identity v2[t+2](c-3) = v2[t+1](c-3) + 1 + kill holds in "
          f"ALL; F (even kill => mu_(t+1) <= c-3) holds in all "
          f"{n_case} even-kill states; zero-input two-step drop "
          f"mu_(t+2) <= mu_t - 3 holds on all {n_two_step} pairs: PASS "
          f"({round(time.time() - t0, 1)} s)", flush=True)

    # random wide states (width 14, 10^6): identity + F + full-window
    # version through the reference f_steps
    rng = random.Random(10)
    n_rand = 0
    n_ch = 0
    for _ in range(1_000_000):
        rows = [sum(1 << c for c in rng.sample(range(0, 14),
                                               rng.randint(0, 5)))
                for _ in range(4)]
        v1m2, v1p1, v2t, v2p1 = [r << 8 for r in rows]
        if not v2t:
            continue
        c = lsb(v2t)
        v2p2 = v1p1 ^ (v1p1 >> 1) ^ (v1m2 << 1) ^ v2p1 ^ (v2t >> 3)
        kill = (v1p1 >> (c - 3) & 1) + (v1p1 >> (c - 2) & 1) \
            + (v1m2 >> (c - 4) & 1)
        assert (v2p2 >> (c - 3) & 1) == ((v2p1 >> (c - 3) & 1) + 1
                                         + kill) % 2
        if kill % 2 == 0:
            assert lsb(v2p1 | v2p2) <= c - 3
        else:
            n_ch += 1
        n_rand += 1
    out["F_random"] = dict(width=14, states=n_rand, charged=n_ch)
    print(f"LEMMA F random (width 14): {n_rand} states, {n_ch} charged"
          f" (odd kill), zero violations: PASS", flush=True)

    # F2 per-cell charge count: a cell (tau, c') can charge at most
    # the steps t = tau-1 (m_t in {c'+2, c'+3}) and t = tau+2
    # (m_t = c'+4) — enumerate the kill-set membership relation
    memb = set()
    for (dt, dc) in ((1, -3), (1, -2), (-2, -4)):
        memb.add((dt, dc))          # cell at (t+dt, m_t+dc)
    steps_per_cell = {}
    for (dt, dc) in memb:
        steps_per_cell.setdefault(-dt, []).append(-dc)
    out["F2_cell_charges"] = {str(k): v for k, v in steps_per_cell.items()}
    assert set(steps_per_cell) == {-1, 2}
    print(f"F2: a K1 cell charges at most 2 steps (t = tau-1 with m_t "
          f"in {sorted(steps_per_cell[-1])} above it, t = tau+2 with "
          f"m_t = {steps_per_cell[2]} above it): holding T rows costs "
          f">= ceil(T/4) input cells left of the frontier", flush=True)

    # mirror lane (x-sector): u1[t+1] = x^-3((1+x^-1)u2[t] + x u2[t-3]
    # + u1[t] + u1[t-1]).  Zero-input: u1[t+1] = x^-3(u1[t]+u1[t-1])
    # — every shift is -3, so the mirror needs no echo argument at
    # all: the WHOLE new row lies >= 3 columns left of the pair's
    # max (a gauge-free statement; the pair-MIN can be held for a
    # few rows by cancellation, e.g. ({m},{m,m+3}) -> ({m,m+3},{m})
    # -> ({m},{m}) -> ({m},0) -> ({m-3}), so the y-lane's two-row
    # pair-min drop is NOT the right mirror invariant — recorded).
    # Exhaustive: max(u1[t+1]) <= max(u1[t] u u1[t-1]) - 3, and
    # min(u1[t+1]) >= min(u1[t] u u1[t-1]) - 3, zero-input.
    n_m = 0
    viol_m = 0
    held = 0
    for a in range(1 << W):
        for b in range(1 << W):
            if not (a | b):
                continue
            A_ = a << OFF
            B_ = b << OFF
            C_ = (A_ ^ B_) >> 3
            if not C_:
                continue
            n_m += 1
            if C_.bit_length() - 1 > (A_ | B_).bit_length() - 1 - 3:
                viol_m += 1
            if lsb(C_) < lsb(A_ | B_) - 3:
                viol_m += 1
            D_ = (B_ ^ C_) >> 3
            if (C_ | D_) and lsb(C_ | D_) > lsb(A_ | B_) - 3:
                held += 1
    assert viol_m == 0
    out["mirror_zero_input"] = dict(
        width=W, pairs=n_m, max_gauge_violations=0,
        pair_min_held_two_rows=held,
        note="mirror zero-input rows shift by exactly -3 as a whole; "
             "the pair-min two-row drop is FALSE there (cancellation "
             "holds), the max gauge carries the mirror R0")
    print(f"mirror lane: zero-input new row lies entirely in "
          f"[min-3, max-3] of the pair on all {n_m} width-{W} pairs "
          f"(max-gauge R0, gauge-free); NOTE the y-lane two-row "
          f"pair-min drop FAILS in the mirror on {held} pairs "
          f"(cancellation holds the min) — the mirror financed half "
          f"lives in the max gauge, listed residue", flush=True)


# ---------------------------------------------------------------------
# lane: naive — the local per-cell count is FALSE
# ---------------------------------------------------------------------

def march_zero_input(w8, T):
    """Forward march T rows with NO input (v1 rows above the window
    empty); returns rows list [(set v1, set v2)] from the window's
    bottom row, in the window's own columns."""
    v1 = [w8[0], w8[1], w8[2], w8[3]]
    v2 = [w8[4], w8[5], w8[6], w8[7]]
    for _ in range(T):
        nxt = v2_forced(v1[-4], v1[-1], v2[-2], v2[-1])
        v1.append(0)
        v2.append(nxt)
    rows = []
    for a, b in zip(v1, v2):
        rows.append((frozenset(c for c in range(a.bit_length())
                               if a >> c & 1),
                     frozenset(c for c in range(b.bit_length())
                               if b >> c & 1)))
    return rows


def lane_naive(out):
    from a40_s6_drift import CoverFragment
    t0 = time.time()
    OFF = 16
    # (i) the hand construction (§15): window rows t-3..t
    c = OFF
    w8 = (0, 1 << (c - 1), 0, 0,                        # v1[t-3..t]
          sum(1 << x for x in (c + 1, c + 2, c + 12, c + 16, c + 19)),
          sum(1 << x for x in (c + 6, c + 9, c + 16)),
          sum(1 << x for x in (c + 6, c + 13)),
          1 << (c + 3))                                  # v2[t-3..t]
    rows = march_zero_input(w8, 5)
    fr = CoverFragment(rows, 0)
    assert fr.admissible()
    anch = fr.anchors()
    rise = anch[-1] - anch[0]
    k1 = sum(len(r[0]) for r in rows[4:])
    print(f"naive lane (i) hand fragment: slabs {fr.slabs()}, anchors "
          f"{[a - OFF for a in anch]}, net rise {rise} over "
          f"{len(anch) - 1} slabs with {k1} input cells in the stretch "
          f"— rise > 3 + #K1 = {3 + k1}: the local per-cell count is "
          f"FALSE (E-admissible: {fr.admissible()})", flush=True)
    assert rise > 3 + k1
    out["hand_witness"] = dict(
        rows=[[sorted(x - OFF for x in a), sorted(x - OFF for x in b)]
              for a, b in rows], slabs=fr.slabs(),
        anchors=[a - OFF for a in anch], rise=rise, k1_in_stretch=k1,
        window_weight=sum(wt(r) for r in w8), admissible=True)
    # (ii) random light/heavy windows, zero-input continuation: the
    # distribution of max rise; the lightest window with rise >= 4
    rng = random.Random(11)
    best = {}
    n = 0
    for _ in range(300_000):
        wgt = rng.randint(2, 9)
        cells = set()
        while len(cells) < wgt:
            cells.add((rng.randint(0, 7), rng.randint(0, 11)))
        w8 = [0] * 8
        for r, col in cells:
            w8[r] |= 1 << (col + OFF)
        w8 = tuple(w8)
        allm = 0
        for r in w8:
            allm |= r
        rows = march_zero_input(w8, 6)
        fr = CoverFragment(rows, 0)
        anch = fr.anchors()
        if any(a is None for a in anch):
            continue
        n += 1
        r_max = max(anch[j] - anch[0] for j in range(len(anch)))
        if r_max >= 4:
            key = wgt
            tooth = tooth_ok(w8[0], w8[1], w8[2], w8[4],
                             v2_forced(w8[0], w8[3], w8[6], w8[7]))
            rec = dict(window_weight=wgt, rise=r_max, tooth_ok=tooth,
                       rows=[[sorted(x - OFF for x in a),
                              sorted(x - OFF for x in b)]
                             for a, b in rows],
                       slabs=fr.slabs(),
                       anchors=[a - OFF for a in anch])
            if key not in best or best[key]["rise"] < r_max:
                best[key] = rec
    lightest = min(best) if best else None
    out["random_zero_input"] = dict(
        windows=n, by_weight={str(k): dict(rise=v["rise"],
                                           tooth_ok=v["tooth_ok"])
                              for k, v in sorted(best.items())},
        lightest_weight_with_rise_ge4=lightest,
        lightest_witness=best.get(lightest))
    print(f"naive lane (ii): {n} random windows marched zero-input 6 "
          f"rows; max net rise >= 4 realized at window weights "
          f"{sorted(best)} (max rises "
          f"{[best[k]['rise'] for k in sorted(best)]}); lightest "
          f"{lightest} — every one is a zero-input (#K1 = 0) rise "
          f"onto a pre-existing right strand ({round(time.time() - t0, 1)}"
          f" s)", flush=True)
    if lightest is not None:
        wrec = best[lightest]
        fr = CoverFragment([(frozenset(x + OFF for x in a),
                             frozenset(x + OFF for x in b))
                            for a, b in wrec["rows"]], 0)
        assert fr.admissible()
        print(f"  lightest witness (weight {lightest}, rise "
              f"{wrec['rise']}, tooth_ok {wrec['tooth_ok']}): rows "
              f"{wrec['rows']}, slabs {wrec['slabs']}, anchors "
              f"{wrec['anchors']}", flush=True)


# ---------------------------------------------------------------------
# lane: verify — F on real objects
# ---------------------------------------------------------------------

def _summarize(recs):
    return dict(steps=len(recs), drops=sum(r["drop"] for r in recs),
                charged=sum(r["charged"] for r in recs),
                both=sum(r["drop"] and r["charged"] for r in recs),
                identity_fail=sum(not r["identity_ok"] for r in recs),
                violations=sum(r["violation"] for r in recs))


def lane_verify(out, gcap, seeds_lo, seeds_hi, hcap=10, dcap=16):
    from a40_s6_drift import CoverFragment
    t0 = time.time()
    res = {}
    # (a) the S9 slip -8 witness
    sp = json.loads((DATA / "s9_specimen_slipm8.json").read_text())
    rows = [(frozenset(a), frozenset(b)) for a, b in sp["best"]["rows"]]
    fr = CoverFragment(rows, 0)
    assert fr.admissible()
    recs = f_steps(rows, 0)
    s = _summarize(recs)
    assert s["violations"] == 0 and s["identity_fail"] == 0
    res["s9_slipm8_witness"] = s
    print(f"verify (a) S9 slip -8 witness: {s} — F holds on every "
          f"step", flush=True)
    # (b) the species lifts
    _argv = sys.argv
    sys.argv = [_argv[0], "12", "8"]
    from a40_s6_drift import lift_phase, load_survivor
    from a40_s5_lightcore import Phase
    sys.argv = _argv
    for name, l, p, d, w in [("W7_l18", 18, 7, 16, 8),
                             ("W7_l24", 24, 7, 22, 8),
                             ("TC63_l18", 18, 6, 3, 10),
                             ("TC63_l24", 24, 6, 3, 10)]:
        fname = ("s5_dense_p7.json" if (l, p) == (18, 7) else
                 "s5_dense_p6.json" if (l, p) == (18, 6) else
                 f"s5_dense_l24p{p}.json")
        ph = Phase.from_quotient_pts(l, p, d, load_survivor(fname, p, d, w))
        frs, s_ = lift_phase(ph, n_periods=3)
        recs = f_steps(frs.rows, frs.t0)
        s = _summarize(recs)
        assert s["violations"] == 0 and s["identity_fail"] == 0, name
        res[name] = dict(s, drift_per_period=s_)
        print(f"verify (b) {name}: {s} — F holds on every step",
              flush=True)
    # (c) parented replay of EVERY completed crossing below gcap
    from a40_s7_tax import ParentedLinkMarch, replay
    sds = [w8 for w8 in seeds_full(1) if sum(wt(r) for r in w8) == 1]
    tot = dict(fragments=0, steps=0, drops=0, charged=0, both=0,
               violations=0, identity_fail=0, naive_from_seed_max=None,
               naive_from_seed_violations=0)
    hold_max = 0
    for si in range(seeds_lo, min(seeds_hi, len(sds))):
        w8 = sds[si]
        m = ParentedLinkMarch(1, kmax=2, whcap=14, gcap=gcap, hcap=hcap,
                              dcap=dcap)
        info = m.run([w8], log=False)
        nkeys = 0
        for key in list(m.parents):
            dyn, anch, phase, L, h, dlt = key
            if phase != POST:
                continue
            rows_d, chain = replay(m, key, [w8])
            lo2, hi2 = min(rows_d), max(rows_d)
            rows = [rows_d[j] for j in range(lo2, hi2 + 1)]
            fr = CoverFragment(rows, lo2)
            if not fr.admissible():
                continue
            nkeys += 1
            recs = f_steps(rows, lo2)
            s = _summarize(recs)
            tot["fragments"] += 1
            for k in ("steps", "drops", "charged", "both",
                      "violations", "identity_fail"):
                tot[k] += s[k]
            # the naive count FROM THE SEED: A_h - A_seed vs 3 + K1
            anch_ = fr.anchors()
            k1 = 0
            worst = None
            for j in range(len(anch_)):
                if j > 0:
                    k1 += len(rows[j + 3][0])       # v1 cells of row
                rise = anch_[j] - anch_[0]
                slack = 3 + k1 - rise
                if worst is None or slack < worst:
                    worst = slack
                if rise > 3 + k1:
                    tot["naive_from_seed_violations"] += 1
            if tot["naive_from_seed_max"] is None or \
                    worst < tot["naive_from_seed_max"]:
                tot["naive_from_seed_max"] = worst
            # longest hold (no drop) run
            run = 0
            for r in recs:
                run = run + 1 if not r["drop"] else 0
                hold_max = max(hold_max, run)
        print(f"  seed {si}: {info['popped']} nodes, {nkeys} POST "
              f"fragments replayed+verified, rss {_rss()}MB; running "
              f"totals {tot}", flush=True)
        del m
    tot["longest_no_drop_run"] = hold_max
    assert tot["violations"] == 0 and tot["identity_fail"] == 0
    res["census_replay"] = dict(gcap=gcap, hcap=hcap, dcap=dcap,
                                seeds=[seeds_lo, seeds_hi], **tot)
    print(f"verify (c) census replay g<={gcap}: {tot} — F holds on "
          f"every step of every completed crossing; the naive "
          f"from-seed count holds with min slack "
          f"{tot['naive_from_seed_max']} (the +3 wall)", flush=True)
    res["wall_s"] = round(time.time() - t0, 1)
    out["verify"] = res


def _rss():
    import os
    import subprocess
    try:
        o = subprocess.check_output(["ps", "-o", "rss=", "-p",
                                     str(os.getpid())])
        return int(o.split()[0]) // 1024
    except Exception:
        return -1


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("lanes", nargs="*", default=["lemma", "naive",
                                                 "verify"])
    ap.add_argument("--gcap", type=int, default=26)
    ap.add_argument("--seeds", type=str, default="0:8")
    ap.add_argument("--hcap", type=int, default=10)
    ap.add_argument("--dcap", type=int, default=16)
    ap.add_argument("--log", type=str, default="")
    args = ap.parse_args()
    if args.log:
        fh = open(args.log, "a", buffering=1)
        sys.stdout = sys.stderr = fh
    t0 = time.time()
    from bb_lab.tower import validate_banked
    validate_banked(LAB / "data")
    print("validate_banked: PASS", flush=True)
    out = {}
    if "lemma" in args.lanes:
        lane_lemma(out)
    if "naive" in args.lanes:
        lane_naive(out)
    if "verify" in args.lanes:
        lo, hi = map(int, args.seeds.split(":"))
        lane_verify(out, args.gcap, lo, hi, args.hcap, args.dcap)
    out["wall_s"] = round(time.time() - t0, 1)
    tag = "_".join(args.lanes)
    p = DATA / f"s10_financed_{tag}.json"
    p.write_text(json.dumps(out, indent=1))
    print(f"wrote {p} ({out['wall_s']} s)", flush=True)


if __name__ == "__main__":
    main()
