#!/usr/bin/env python3
"""A42 S2a — the meet-in-the-middle sandwich, analyzed mechanically.

Session 1 named "forward level-11 + backward level-11, join on states
(11 + 11 = 22)" as the p = 12 unrestricted closer.  This script grounds
the join geometry end-to-end and measures the two failure modes:

  (1) BOTH-STANDARD join: fwd cost P_F(g) = |v1[<=g]| + |v2[<=g+2]|
      and bwd cost Q_B(g) = |v1[>=g-1]| + |v2[>=g]| DOUBLE-PAY the
      5-column state block s(g); P_F + Q_B = w + s(g) at every cut, so
      caps (11, 11) reach only w <= 22 - s(cut): real cycles with no
      light cut escape.  Measured here on every atlas cycle.
  (2) EXACT-PARTITION join (bwd charged on exit, "exclusive"):
      P_F(g) + Bx(g) = w with Bx = Q_B - s, and the max-cut argument
      covers all w <= Wcap with caps (CF, Wcap - CF - 1) — but the
      exclusive frontier contains the free-suffix shell (the last
      three v2-inputs and induced forced columns ride unpaid):
      |shell| = 2^{3p} states at exclusive cost 0.  Counted here at
      p = 3, 4, 5 (exact) and extrapolated.

Also banked: the reversed automaton (x -> x^-1 substitution into the
same generic Automaton class) is verified mechanically — every atlas
cycle at p in {3, 6} replays through BOTH directions (fwd consumes
v1-columns left-to-right, bwd consumes v2-columns right-to-left, each
emitting the other block's columns exactly, both returning to the
zero state with accumulated cost = weight), and the per-cut state
content + cost formulas are asserted against the replay accumulators
at every step.
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

LAB = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(LAB / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import importlib.util

_spec = importlib.util.spec_from_file_location(
    "a40_s4_phase_atlas", Path(__file__).parent / "a40_s4_phase_atlas.py")
AT = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(AT)

DATA = LAB / "data" / "a42"


def xinv(supp):
    return [(-i, j) for (i, j) in supp]


def pc(x):
    return bin(x).count("1")


def cols_from_pts(pts):
    """pts: {(c, y, blk)} -> (v1cols, v2cols) dicts c -> p-bit ints.
    blk 0 = v1 (input block of the fwd automaton), blk 1 = v2."""
    v1, v2 = {}, {}
    for (c, y, blk) in pts:
        d = v1 if blk == 0 else v2
        d[c] = d.get(c, 0) | (1 << y)
    return v1, v2


def check_cycle_conv(v1, v2, p):
    """Independent re-verification: Bbar v1 + Abar v2 = 0 over
    F2[x^pm] tensor F2[y]/(y^p-1), directly from the supports."""
    acc = {}
    for (i, j) in [(-i, -j) for (i, j) in AT.B_L]:      # Bbar
        for c, col in v1.items():
            acc[c + i] = acc.get(c + i, 0) ^ AT.rot(col, j, p)
    for (i, j) in [(-i, -j) for (i, j) in AT.A_L]:      # Abar
        for c, col in v2.items():
            acc[c + i] = acc.get(c + i, 0) ^ AT.rot(col, j, p)
    assert all(v == 0 for v in acc.values()), "not a cycle (conv)"


def fwd_replay_with_cuts(au, v1, v2, p):
    """Replay the fwd automaton over the cycle's own input columns;
    assert state content + P_F formula at every cut; return the
    per-cut table [(gamma, P_F)] and the total cost."""
    los = min(list(v1) + list(v2)) - 6
    his = max(list(v1) + list(v2)) + 6
    st = au.zero()
    cost = 0
    cuts = []
    for g in range(los, his + 1):
        # state before step at equation-index g:
        # (v2[g], v2[g+1], v2[g+2], v1[g-1], v1[g])
        want = (v2.get(g, 0), v2.get(g + 1, 0), v2.get(g + 2, 0),
                v1.get(g - 1, 0), v1.get(g, 0))
        assert st == want, (g, st, want)
        pf = sum(w for c, col in v1.items() if c <= g
                 for w in [pc(col)]) + \
            sum(w for c, col in v2.items() if c <= g + 2
                for w in [pc(col)])
        assert cost == pf, (g, cost, pf)
        cuts.append((g, cost))
        st, c2, newf = au.step(st, v1.get(g + 1, 0))
        assert newf == v2.get(g + 3, 0), (g, "forced v2 mismatch")
        cost += c2
    assert st == au.zero() and all(x == 0 for x in st)
    return cuts, cost


def bwd_replay_with_cuts(bau, v1, v2, p):
    """Replay the REVERSED automaton (x -> x^-1 pair) over the
    mirrored data w_i[d] = v_i[-d]: inputs are the v2-columns in
    decreasing original c; forced outputs must equal the v1-columns.
    Assert state content + Q_B formula at every cut."""
    w1 = {-c: col for c, col in v1.items()}
    w2 = {-c: col for c, col in v2.items()}
    lod = min(list(w1) + list(w2)) - 6
    hid = max(list(w1) + list(w2)) + 6
    st = bau.zero()
    cost = 0
    qb = {}
    for d in range(lod, hid + 1):
        # state before step at index d:
        # f = (w1[d-1], w1[d]); o = (w2[d-3], w2[d-2], w2[d-1])
        want = (w1.get(d - 1, 0), w1.get(d, 0),
                w2.get(d - 3, 0), w2.get(d - 2, 0), w2.get(d - 1, 0))
        assert st == want, (d, st, want)
        # paid so far: w1[<= d-1] + w2[<= d-1]  (inputs w2 consumed to
        # d-1, forced w1 to d... f ends at d: w1[d] emitted by the
        # PREVIOUS step) -> original coords: v1[>= -(d-1)... ] etc.
        paid = sum(pc(c2) for dd, c2 in w1.items() if dd <= d) + \
            sum(pc(c2) for dd, c2 in w2.items() if dd <= d - 1)
        assert cost == paid, (d, cost, paid)
        # original-coordinate cut: block(e) with e = -d matches the
        # fwd block at gamma = e + 1 = -d + 1:
        # Q_B(gamma) = |v1[>= gamma-1]| + |v2[>= gamma]|
        gamma = -d + 1
        qbv = sum(pc(c2) for c, c2 in v1.items() if c >= gamma - 1) + \
            sum(pc(c2) for c, c2 in v2.items() if c >= gamma)
        assert qbv == paid, (d, gamma, qbv, paid)
        qb[gamma] = qbv
        st, c2, newf = bau.step(st, w2.get(d, 0))
        assert newf == w1.get(d + 1, 0), (d, "forced w1 mismatch")
        cost += c2
    assert st == bau.zero()
    return qb, cost


def cut_diagnostics(v1, v2, p, w):
    """Per-cut P_F, Q_B, s, Bx; verify identities; return balance
    stats for the two join schemes."""
    los = min(list(v1) + list(v2)) - 4
    his = max(list(v1) + list(v2)) + 4
    bal_std = None
    bal_part = None
    smin_int = None
    for g in range(los, his + 1):
        pf = sum(pc(c2) for c, c2 in v1.items() if c <= g) + \
            sum(pc(c2) for c, c2 in v2.items() if c <= g + 2)
        qbv = sum(pc(c2) for c, c2 in v1.items() if c >= g - 1) + \
            sum(pc(c2) for c, c2 in v2.items() if c >= g)
        s = pc(v1.get(g - 1, 0)) + pc(v1.get(g, 0)) + \
            pc(v2.get(g, 0)) + pc(v2.get(g + 1, 0)) + pc(v2.get(g + 2, 0))
        assert pf + qbv == w + s, (g, pf, qbv, w, s)
        bx = qbv - s
        assert pf + bx == w
        m1 = max(pf, qbv)
        m2 = max(pf, bx)
        bal_std = m1 if bal_std is None else min(bal_std, m1)
        bal_part = m2 if bal_part is None else min(bal_part, m2)
        if 0 < pf < w:                      # interior cut
            smin_int = s if smin_int is None else min(smin_int, s)
    return bal_std, bal_part, smin_int


def exclusive_shell(p: int, cap: int = 2_000_000):
    """Count the zero-exclusive-cost reachable shell of the REVERSED
    automaton (exit-charged): all states reachable from zero through
    edges whose source has f[0] = o[0] = 0 (exit cost 0)."""
    bau = AT.Automaton(xinv(AT.A_L), xinv(AT.B_L), p)
    z = bau.zero()
    seen = {z}
    frontier = [z]
    while frontier:
        nxt = []
        for st in frontier:
            if pc(st[0]) + pc(st[2]) != 0:      # exits of source
                continue
            for a in range(1 << p):
                st2, _, _ = bau.step(st, a)
                if st2 not in seen:
                    seen.add(st2)
                    if len(seen) > cap:
                        return len(seen), False
                    nxt.append(st2)
        frontier = nxt
    return len(seen), True


def main():
    t0 = time.time()
    out = {"identities_checked": 0, "cycles": []}

    # -- the reversed automaton: derived shape --------------------
    bau6 = AT.Automaton(xinv(AT.A_L), xinv(AT.B_L), 6)
    shape = dict(nf=bau6.nf, no=bau6.no, forced_blk=bau6.forced_blk,
                 adv_f=bau6.adv_f, adv_o=bau6.adv_o, top_j=bau6.top_j,
                 terms_f=sorted(bau6.terms_f),
                 terms_o=sorted(bau6.terms_o))
    print(f"reversed automaton shape: {shape}", flush=True)
    assert (bau6.nf, bau6.no, bau6.forced_blk) == (2, 3, 0)
    out["reversed_shape"] = {k: v for k, v in shape.items()
                             if k not in ("terms_f", "terms_o")}

    # -- harvest real cycles and run every check ------------------
    for p, W in ((3, 8), (6, 12)):
        au = AT.Automaton(AT.A_L, AT.B_L, p)
        bau = AT.Automaton(xinv(AT.A_L), xinv(AT.B_L), p)
        rows, npop = AT.atlas("AB", p, W, keep_pts=True)
        print(f"p={p} W<={W}: {len(rows)} cycles harvested "
              f"({npop} states)", flush=True)
        stats = {}
        for r in rows:
            v1, v2 = cols_from_pts(r["pts"])
            w = r["weight"]
            assert sum(pc(c) for c in v1.values()) + \
                sum(pc(c) for c in v2.values()) == w
            check_cycle_conv(v1, v2, p)
            cuts, cf = fwd_replay_with_cuts(au, v1, v2, p)
            assert cf == w, (cf, w)
            qb, cb = bwd_replay_with_cuts(bau, v1, v2, p)
            assert cb == w, (cb, w)
            b_std, b_part, smin = cut_diagnostics(v1, v2, p, w)
            out["identities_checked"] += 1
            key = (w, r["nontrivial"])
            ex = b_std - (w + 1) // 2
            st = stats.setdefault(key, dict(n=0, bal_std_min=99,
                                            bal_std_max=0, excess_max=0,
                                            bal_part_max=0, smin_min=99))
            st["n"] += 1
            st["bal_std_min"] = min(st["bal_std_min"], b_std)
            st["bal_std_max"] = max(st["bal_std_max"], b_std)
            st["excess_max"] = max(st["excess_max"], ex)
            st["bal_part_max"] = max(st["bal_part_max"], b_part)
            st["smin_min"] = min(st["smin_min"], 99 if smin is None
                                 else smin)
        for (w, nt), st in sorted(stats.items()):
            print(f"  w={w} {'NONTRIVIAL' if nt else 'trivial':10s} "
                  f"n={st['n']:5d}  bal_std in "
                  f"[{st['bal_std_min']},{st['bal_std_max']}] "
                  f"(excess over w/2 up to +{st['excess_max']})  "
                  f"bal_exactpart <= {st['bal_part_max']}  "
                  f"min interior block s >= {st['smin_min']}",
                  flush=True)
            out["cycles"].append(dict(
                p=p, w=w, nontrivial=nt, **st))

    # -- the exclusive-frontier shell (the monster) ---------------
    out["shell"] = {}
    for p in (3, 4, 5):
        n, complete = exclusive_shell(p)
        print(f"exclusive zero-cost shell p={p}: {n} states "
              f"(2^(3p) = {1 << (3 * p)}; complete={complete})",
              flush=True)
        out["shell"][str(p)] = dict(count=n, pow23p=1 << (3 * p),
                                    complete=complete)

    out["wall_s"] = round(time.time() - t0, 1)
    (DATA / "s2_sandwich_analysis.json").write_text(
        json.dumps(out, indent=1))
    print(f"wrote {DATA/'s2_sandwich_analysis.json'} "
          f"({out['wall_s']} s)", flush=True)


if __name__ == "__main__":
    main()
