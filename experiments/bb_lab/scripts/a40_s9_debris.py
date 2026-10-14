#!/usr/bin/env python3
"""A40 S9 — the debris accounting behind the rightward-relocation
inequality: machine-verified local identities (V1-V4) + the derived
constants.  These are the checkable steps of the hand argument
(§14); the census (s9_slip_u1_g{32,34,40}.json) is the measured
side.

THE MECHANISM (y-lane; solved recurrence
  v2[t+1] = (1+x^-1) v1[t] + x v1[t-3] + v2[t] + x^-3 v2[t-1]):

* BORROW RADIUS c0 = 3 (structural): reading E at column c, the
  only sources STRICTLY RIGHT of c are v2[t-1](c+3) (the x^-3
  term) and v1[t](c+1); the farthest is +3.  A dying left column
  can be financed by existing content at distance <= 3 — and no
  farther.  (Mirror lane: radius 4, from x^-4.)
* V1 (persistence): a v2 cell at the local left frontier persists
  upward unless one of exactly four positions is occupied; at a
  true frontier the left option is empty, so the killer is at
  {(c,t) blk1, (c+1,t) blk1, (c+3,t-1) blk2}.
* V2 (leftward spawn): a v1 cell at a true left frontier forces
  v2[t+1](c-1) = 1 (frontier advances LEFT) unless financed by
  v2[t-1](c+2) — leftward motion is the free direction.
* V3 (anchor rise = exiting-row event): the window anchor can rise
  only when every window cell left of the new anchor sits in the
  single exiting row.
* V4 (the K1 retreat step, in vivo): killing a frontier v2 cell at
  c by the input v1[t] = {c+1} moves the frontier EXACTLY one
  column right at the cost of exactly ONE input cell, and re-seeds
  v2 at c+1 (the retreat diagonal: <= 1 column/row, >= 1 cell/
  column, cells pairwise distinct by column).

CONSEQUENCE (hand, §14): a net rightward anchor transport of
delta >= 4 columns cannot be financed by borrows alone (a borrow
cell for column c sits at c+3 < A_end when delta >= 4, so it must
itself die — the recursion terminates only in K1 inputs); each
column beyond the borrow radius is charged >= 1 input cell, at
pairwise-distinct positions, concentrated on the transition rows
(which therefore go heavy — the measured ~8/row holding wall)."""
from __future__ import annotations

import json
import random
import sys
import time
from pathlib import Path

LAB = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(LAB / "src"))
sys.path.insert(0, str(LAB / "scripts"))
DATA = LAB / "data" / "a40"

from a40_s6_frontier import wt, lsb  # noqa: E402


def v2new_kernel(v1a, v1b, v1c, v1d, v2a, v2b, v2c):
    """The S6/S7/S8 forced-row bit kernel (v1a = v1[t-3], v1d =
    v1[t], v2a..v2c = v2[t-2..t])."""
    return v1d ^ (v1d >> 1) ^ (v1a << 1) ^ v2c ^ (v2b >> 3)


def main():
    t0 = time.time()
    from bb_lab.tower import validate_banked
    validate_banked(LAB / "data")
    print("validate_banked: PASS", flush=True)
    out = {}
    rng = random.Random(9)

    # ---- V0: the borrow radius, derived from the term supports ----
    # E at column c reads: v1[t]{c, c+1}, v1[t-3]{c-1}, v2[t]{c},
    # v2[t+1]{c}, v2[t-1]{c+3}.  Source offsets relative to c:
    src_offsets = dict(v1_t=(0, 1), v1_tm3=(-1,), v2_t=(0,),
                       v2_tp1=(0,), v2_tm1=(3,))
    radius = max(o for offs in src_offsets.values() for o in offs)
    assert radius == 3
    # mirror: u1[t+1] = x^-3((1+x^-1)u2[t] + x u2[t-3] + u1[t]
    #         + u1[t-1]) — source offsets at column c:
    m_offsets = (3, 4, 2, 3, 3)
    assert max(m_offsets) == 4
    out["borrow_radius"] = dict(y=3, mirror=4)
    print("V0 borrow radius: y-lane 3 (the x^-3 term of B), "
          "mirror 4 — the ONLY rightward reaches in E", flush=True)

    # ---- V1: v2 frontier persistence (200k random states) --------
    n_checked = 0
    for _ in range(200_000):
        rows = [sum(1 << c for c in rng.sample(range(2, 14),
                                               rng.randint(0, 4)))
                for _ in range(7)]
        v1a, v1b, v1c, v1d, v2a, v2b, v2c = rows
        v2n = v2new_kernel(*rows)
        # every column c occupied in v2c and vacated in v2n must
        # have odd parity among the four E-sources at c:
        for c in range(0, 16):
            if v2c >> c & 1 and not (v2n >> c & 1):
                k = ((v1d >> c & 1) + (v1d >> (c + 1) & 1)
                     + ((v1a >> (c - 1) & 1) if c >= 1 else 0)
                     + (v2b >> (c + 3) & 1))
                assert k % 2 == 1, (c, rows)
                n_checked += 1
    out["V1_persistence_checks"] = n_checked
    print(f"V1 persistence: {n_checked} vacation events, every one "
          f"killed by an odd source set among v1[t](c), v1[t](c+1),"
          f" v1[t-3](c-1), v2[t-1](c+3): PASS", flush=True)

    # ---- V2: v1 frontier leftward spawn --------------------------
    n2 = 0
    for _ in range(200_000):
        rows = [sum(1 << c for c in rng.sample(range(3, 14),
                                               rng.randint(0, 4)))
                for _ in range(7)]
        v1a, v1b, v1c, v1d, v2a, v2b, v2c = rows
        allm = v1a | v1b | v1c | v1d | v2a | v2b | v2c
        if not v1d:
            continue
        c = lsb(v1d)
        if c != lsb(allm):
            continue                       # v1[t] must hold the min
        v2n = v2new_kernel(*rows)
        spawn = (v2n >> (c - 1)) & 1
        borrow = (v2b >> (c + 2)) & 1
        # at a global window min in v1[t]: v1[t](c-1) = 0,
        # v1[t-3](c-2) = 0, v2[t](c-1) = 0 hold automatically, so
        # spawn XOR borrow = 1 exactly:
        assert spawn ^ borrow == 1, (c, rows)
        n2 += 1
    out["V2_spawn_checks"] = n2
    print(f"V2 leftward spawn: {n2} frontier-v1 states, "
          f"v2[t+1](c-1) = 1 XOR borrowed at v2[t-1](c+2) in every "
          f"one (left motion is the free direction): PASS",
          flush=True)

    # ---- V3: anchor rise = exiting-row event ---------------------
    n3 = 0
    for _ in range(300_000):
        rows = [sum(1 << c for c in rng.sample(range(2, 12),
                                               rng.randint(0, 3)))
                for _ in range(7)]
        v1a, v1b, v1c, v1d, v2a, v2b, v2c = rows
        allm = v1a | v1b | v1c | v1d | v2a | v2b | v2c
        if not allm:
            continue
        A_old = lsb(allm)
        v2n = v2new_kernel(*rows)
        s_in = sum(1 << c for c in rng.sample(range(2, 12),
                                              rng.randint(0, 3)))
        new_win = v1b | v1c | v1d | s_in | v2a | v2b | v2c | v2n
        if not new_win:
            continue
        A_new = lsb(new_win)
        if A_new > A_old:
            passed = allm & ((1 << A_new) - 1)
            # every passed cell lies ONLY in the exiting row v1a
            # (v2[t-3] is already outside the 7-row state), and in
            # particular the passed columns are a v1a-only feature:
            assert (passed & (v1b | v1c | v1d | v2a | v2b | v2c)) \
                == 0, (rows, s_in)
            assert passed & v1a == passed, (rows, s_in)
            n3 += 1
    out["V3_anchor_rise_checks"] = n3
    print(f"V3 anchor rise: {n3} rising steps, in every one the "
          f"passed columns lived ONLY in the exiting row: PASS",
          flush=True)

    # ---- V4: the K1 retreat step, in vivo ------------------------
    # A pure v2 strand at column c, no borrow content.  (i) With no
    # input: it PERSISTS (v2new(c) = 1).  (ii) With the single
    # input v1[t] = {c+1}: it DIES at c and re-seeds at c+1 —
    # exactly one input cell per column of retreat.
    c = 6
    base = (0, 0, 0, 0, 0, 0, 1 << c)        # single v2[t] cell
    v2n_hold = v2new_kernel(*base)
    assert v2n_hold == 1 << c, "strand failed to persist cleanly"
    withk1 = (0, 0, 0, 1 << (c + 1), 0, 0, 1 << c)
    v2n_kill = v2new_kernel(*withk1)
    assert not (v2n_kill >> c & 1), "K1 kill failed"
    assert v2n_kill == 1 << (c + 1), "K1 re-seed at c+1 missing"
    # the borrow variant: v2[t-1](c+3) kills WITHOUT re-seed
    withk2 = (0, 0, 0, 0, 0, 1 << (c + 3), 1 << c)
    v2n_b = v2new_kernel(*withk2)
    assert v2n_b == 0, "K2 borrow kill failed"
    # and the B-term's leftward drive: a 2-row vertical v2 strand
    # sheds x^-3 debris (leftward motion is E-POWERED):
    tall = (0, 0, 0, 0, 0, 1 << c, 1 << c)
    v2n_t = v2new_kernel(*tall)
    assert v2n_t == (1 << c) | (1 << (c - 3)), "x^-3 shed missing"
    out["V4_k1_k2"] = dict(
        persist=True, k1_kill_and_reseed=True, k2_borrow_kill=True,
        xm3_leftward_shed=True)
    print("V4 in vivo: single v2 cell PERSISTS verbatim (1 weight/"
          "row); the 1-cell K1 input at c+1 kills it and re-seeds "
          "EXACTLY at c+1 (retreat = 1 column per row per input "
          "cell); the K2 borrow at (c+3, t-1) kills with NO "
          "re-seed; a 2-row v2 strand sheds x^-3 debris leftward "
          "(left motion E-powered, right motion priced): PASS",
          flush=True)

    # ---- V5: the K2-only channel, EXHAUSTED (pure-v2 dynamics) ---
    # With v1 == 0 the evolution v2[t+1] = v2[t] + x^-3 v2[t-1] is
    # DETERMINISTIC, so the zero-input retreat question is finite:
    # enumerate ALL initial pairs (v2[0], v2[1]) with columns in
    # [0, W) and min column 0 (translation-reduced), march T rows,
    # and record the max sustained rise of m_t = min(v2[t] u
    # v2[t+1]) while both rows stay inside a width guard (paths
    # that die to zero are strand extinctions, not relocations —
    # excluded by requiring nonzero v2 through the readout row).
    W5, T5 = 11, 12
    OFF = 3 * (T5 + 2)              # frame offset so x^-3 shifts
    max_rise = 0                    # never go negative
    arg = None
    n5 = 0
    for r0 in range(1 << W5):
        for r1 in range(1 << W5):
            lo = r0 | r1
            if not lo or lo & 1 == 0:
                continue            # translation-reduce: min col 0
            n5 += 1
            a = r0 << OFF
            b = r1 << OFF
            m_start = lsb(a | b)
            best = 0
            for t in range(T5):
                c = b ^ (a >> 3)
                if not c or not (b | c):
                    break           # extinction, not relocation
                best = max(best, lsb(b | c) - m_start)
                a, b = b, c
            if best > max_rise:
                max_rise = best
                arg = (r0, r1, best)
    out["V5_pure_v2"] = dict(W=W5, T=T5, n_pairs=n5,
                             max_net_rise=max_rise, argmax=arg)
    print(f"V5 pure-v2 exhaustion: {n5} translation-reduced "
          f"initial pairs (width {W5}), {T5} rows — max net "
          f"anchor rise of the zero-input dynamics = {max_rise} "
          f"(argmax {arg})", flush=True)
    assert max_rise == 0, f"pure-v2 rise {max_rise} != 0"

    # ---- LEMMA R0 (mu-echo; the PROOF the exhaustion suggested) --
    # mu_t := min(v2[t] u v2[t+1]).  In zero-input dynamics
    # (v2[t+1] = v2[t] + x^-3 v2[t-1]) mu NEVER rises:
    #   case 1: the min cell c = mu_t lies in v2[t+1] — it is in
    #     the next pair too: mu_{t+1} <= c.
    #   case 2: c in v2[t] only.  Then, reading E at (t+2, c-3):
    #     v2[t+2](c-3) = v2[t+1](c-3) + v2[t](c) = 0 + 1 = 1
    #     (v2[t+1] has nothing < mu_t), UNCONDITIONALLY — so
    #     mu_{t+1} <= c - 3: mu strictly DROPS by >= 3.  QED.
    # (Case 2 is the x^-3 echo of V4's "leftward shed", promoted
    # to the min cell.)  Machine-check of the case-2 identity:
    n_echo = 0
    for _ in range(100_000):
        a = sum(1 << c for c in rng.sample(range(0, 14),
                                           rng.randint(1, 5)))
        b = sum(1 << c for c in rng.sample(range(0, 14),
                                           rng.randint(0, 5)))
        a <<= OFF
        b <<= OFF
        c_ = b ^ (a >> 3)
        cmin = lsb(a | b)
        if not (b >> cmin & 1):         # min cell in v2[t] only
            assert c_ >> (cmin - 3) & 1, "echo missing"
            n_echo += 1
    out["R0_echo_checks"] = n_echo
    print(f"LEMMA R0 (mu-echo, PROVEN 3 lines): zero-input "
          f"evolution never raises the pair-min mu (case 2: an "
          f"a-only min cell echoes at c-3 two rows up, "
          f"unconditionally — {n_echo} random case-2 events "
          f"checked); rightward anchor relocation is "
          f"INPUT-CHARGED, always", flush=True)

    out["wall_s"] = round(time.time() - t0, 1)
    (DATA / "s9_debris.json").write_text(json.dumps(out, indent=1))
    print(f"wrote {DATA/'s9_debris.json'} ({out['wall_s']} s)",
          flush=True)


if __name__ == "__main__":
    main()
