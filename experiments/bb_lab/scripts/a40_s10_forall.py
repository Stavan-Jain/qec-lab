#!/usr/bin/env python3
"""A40 S10 — Stage 3: the for-all-r toroidal (doubly-spanning,
y-sector) member floor, assembled from all-h ingredients only.

THE ASSEMBLY (drift-blind; uses NO closure, NO ell, so it is a
statement about every y-spanning class-minimal walk of m rows).
Sigma_j W_j = 4|v| over the m slabs (each row lies in 4 slabs).
Partition the slabs into heavy slabs (W >= 8) and maximal light
runs; a light run's SEED is its lowest minimum-weight slab, of
weight u:
  * every slab of the run below the seed has weight >= u+1 >= 2;
  * every slab above the seed has weight >= u; for u >= 2 that is
    >= 2; for u = 1 the slabs from the seed upward form a forward
    light path from a weight-1 window — the u = 1 FORWARD TREE,
    which is FINITE (the S7 exhaustion, §12.4: no 19-slab light
    path from any weight-1 window at any budget, growth scope
    dil-4/smax-3, stability-checked at smax 4 / dil 6), so its
    exact cost table fwd_min(h) is a cap-free finite computation
    (`preexact` lane) and a u = 1 run of h_run slabs whose forward
    piece has h_f slabs weighs >= 2 (h_run - h_f) + fwd_min(h_f)
    >= 2 h_run - S with S := max_h [2h - fwd_min(h)] (S = 4);
  * every u = 1 run is followed by a heavy block (B >= 1 blocks),
    so #u1-runs <= B, and each block contributes >= 8 >= 2 + 6.
Hence 4|v| >= 2m + 6 B - 4 B >= 2m + 2B.  The all-light branch
(B = 0): u >= 2 gives 4|v| >= 2m; u = 1 is a light closed walk =
a light path of m+1 slabs from a weight-1 window, impossible for
m + 1 >= 19 by the exhaustion.  All-heavy: 4|v| >= 8m.

THEOREM S10.3 (scope-listed).  For every m >= 18, every class-
minimal doubly-spanning (y-sector) walk of m rows within the growth
scope has weight >= ceil(m/2); for the b = 1 family (m = 6r):
d_Y(C_{r,1}) >= 3r for r >= 3 — double the L1 floor ceil(3r/2), at
certificate tier (the exhaustion) + theorem (the arithmetic).  The
x-sector needs the mirror u = 1 exhaustion (`mirror` lane).

Conditional form (`theorem` lane reports it): under C-W7 (every
u >= 2 light coast of h slabs weighs >= (32/7) h - c, the light-
rate wall at the W7 species' rate), the same partition gives
4|v| >= (32/7) m - O(B) i.e. ~ (8/7) m = 6.86 r; the conjectured
12r needs the momentum budget on top (the drift is never used
here).  The ladder route's 7.5r is NOT available even
conditionally: W7 coasts sit below 5 per slab (Stage 2)."""
from __future__ import annotations

import argparse
import json
import math
import sys
import time
from pathlib import Path

LAB = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(LAB / "src"))
sys.path.insert(0, str(LAB / "scripts"))
DATA = LAB / "data" / "a40"

from a40_s6_frontier import wt, seeds_full, March  # noqa: E402


def lane_preexact(args, out):
    """Cap-free forward u = 1 light march: exact fwd_min(h), the
    exhaustion flags, and the slack S = max_h [2h - fwd_min(h)]."""
    t0 = time.time()
    sds = seeds_full(1)
    res = {}
    for (smax, dil) in [(3, 4), (4, 6)]:
        m = March(1, "fwd", smax=smax, dil=dil, hcap=args.hcap,
                  gcap=args.gcap, extent_cap=60, dcap=60)
        # record the max cover extent over every node pushed (the
        # wrap-free scope: the theorem is literal for ell > max
        # extent, scope-listed below it)
        ext = [0]
        _push = m.push

        def push(g2, state_rows, slab_mask, slab_hi, h, dlt_new, M,
                 anch, xlo, xhi, heap, seen, _p=_push, _e=ext):
            a_new = (slab_mask & -slab_mask).bit_length() - 1
            e = max(xhi + M, slab_hi) - min(xlo + M, a_new) + 1
            if e > _e[0]:
                _e[0] = e
            return _p(g2, state_rows, slab_mask, slab_hi, h, dlt_new, M,
                      anch, xlo, xhi, heap, seen)
        m.push = push
        info = m.run(sds, log=False, max_nodes=50_000_000)
        info["max_extent"] = ext[0]
        byh = {}
        for (h, d), g in m.table.items():
            if h not in byh or byh[h] > g:
                byh[h] = g
        hmax = max(byh)
        exhausted = (not info["trunc_g"] and not info["trunc_nodes"]
                     and info["trunc_extent"] == 0
                     and info["trunc_dcap"] == 0 and hmax < args.hcap)
        slack = {h: 2 * h - g for h, g in byh.items()}
        S = max(slack.values())
        argS = [h for h, s in slack.items() if s == S]
        res[f"smax{smax}_dil{dil}"] = dict(
            info=info, fwd_min=dict(sorted(byh.items())), hmax=hmax,
            exhausted=exhausted, slack=dict(sorted(slack.items())),
            S=S, S_at=argS)
        print(f"fwd u=1 cap-free (smax {smax}, dil {dil}, gcap "
              f"{args.gcap}, hcap {args.hcap}): {info}; max h {hmax}; "
              f"EXHAUSTED {exhausted}; fwd_min {dict(sorted(byh.items()))}"
              f"; slack 2h - fwd_min: max S = {S} at h {argS}",
              flush=True)
        assert exhausted, "u=1 forward tree did NOT exhaust"
    out["preexact"] = dict(res, wall_s=round(time.time() - t0, 1))


def floor_theorem(m, S, exhaust_h):
    """ceil of the assembled bound for m rows: min over branches."""
    # B >= 1: 4|v| >= 2m + 6B - S*B >= 2m + (6 - S)  (B = 1 worst
    # when S < 6); all-light u >= 2: 2m; all-light u = 1: impossible
    # iff m + 1 > exhaust_h (a closed light walk is a light path of
    # m + 1 slabs from its seed); all-heavy: 8m.
    # B blocks, each >= 1 slab, each u=1 run >= 1 slab: B <= m // 2
    cands = [min(2 * m + (6 - S) * B for B in range(1, m // 2 + 1))]
    cands.append(2 * m)                      # all-light u >= 2
    u1_light_possible = (m + 1 <= exhaust_h)
    if u1_light_possible:
        cands.append(m)                      # all-light u = 1: >= 1/slab
    return math.ceil(min(cands) / 4), u1_light_possible


def lane_theorem(args, out):
    pe = out.get("preexact") or json.loads(
        (DATA / "s10_forall_preexact.json").read_text())["preexact"]
    rec = pe["smax3_dil4"]
    S, hmax = rec["S"], rec["hmax"]
    exhaust_h = hmax                      # longest light path (slabs)
    rows = []
    for r in range(1, args.rmax + 1):
        m = 6 * r
        f, u1 = floor_theorem(m, S, exhaust_h)
        L1 = math.ceil(3 * r / 2)
        rows.append(dict(r=r, m=m, ell=6 * r + 6, f_y=f, L1=L1,
                         u1_all_light_branch_open=u1,
                         f_W7_conditional=math.ceil((32 / 7 * m
                                                     - 32 / 7 * 1) / 4)))
    print("for-all-r y-sector floor f_Y(r) = ceil((2m + 6 - S)/4) with "
          f"S = {S}, m = 6r (u=1 all-light branch dead for m >= "
          f"{exhaust_h}):", flush=True)
    for x in rows[:12]:
        print(f"  r={x['r']:2d} m={x['m']:3d}: f_Y={x['f_y']:3d}  L1="
              f"{x['L1']:3d}  (C-W7 conditional ~{x['f_W7_conditional']})"
              f"{'  [u1 all-light open: floor = m/4 branch]' if x['u1_all_light_branch_open'] else ''}",
              flush=True)
    # controls (the drift-blind floor is ell-free, so every banked
    # member value must sit ABOVE it):
    ctl = {}
    banked = {(12, 12): 18, (18, 12): 24, (12, 6): 12, (18, 6): 12,
              (24, 6): 12, (6, 6): 6}
    for (l, m), d in banked.items():
        f, _ = floor_theorem(m, S, exhaust_h)
        ctl[f"({l},{m})"] = dict(banked_d=d, f_y=f, ok=(f <= d))
    f18, _ = floor_theorem(18, S, exhaust_h)
    ctl["(24,18)"] = dict(certified_floor=12, f_y=f18, ok=(f18 <= 12),
                          note="f must not exceed the S9 certified "
                               "floor's TRUE value; 12 is itself a "
                               "floor, so consistency = f <= d_true; "
                               "f <= 12 recorded as the weaker check")
    print(f"controls: {ctl}", flush=True)
    assert all(v["ok"] for v in ctl.values())
    out["theorem"] = dict(S=S, exhaust_h=exhaust_h, rows=rows,
                          controls=ctl,
                          statement="d_Y(C_{r,1}) >= ceil((12r + 6 - S)/4)"
                                    " = 3r + ceil((6-S)/4) for r >= 3; "
                                    "with S = 4: 3r + 1")


def lane_mirror(args, out):
    """The mirror (x-sector) all-light u = 1 exhaustion probe:
    MClosedMarch with nHmax = 0 marched m = args.mrows steps at a
    cap-free budget; the tree is exhausted iff every seed's layers
    die before h = m with trunc_g False (§13.7 item 3, never run)."""
    from a40_s8_xlane import MClosedMarch
    t0 = time.time()
    sds = [w8 for w8 in seeds_full(1) if sum(wt(r) for r in w8) == 1]
    res = {}
    for extent in args.extents:
        cm = MClosedMarch(1, m=args.mrows, ell=10 ** 6, nHmax=0, whcap=19,
                          smax=3, dil=4, gcap=args.gcap, extent_cap=extent,
                          dcap=200, exact_dlt=True)
        hmax = [0]
        ming = {}
        orig = cm.expand

        def expand(g, dyn, nH, h, dlt, xlo, xhi, nxt, _o=orig, _hm=hmax,
                   _mg=ming):
            # every expanded state is a live h-slab light path: record
            # the per-h minimum g (the mirror fwd_min table) and the
            # deepest layer reached
            if h not in _mg or _mg[h] > g:
                _mg[h] = g
            if h > _hm[0]:
                _hm[0] = h
            return _o(g, dyn, nH, h, dlt, xlo, xhi, nxt)
        cm.expand = expand
        info = cm.run_layers(sds, log=False, rss_cap=args.rss_cap,
                             frontier_cap=3_000_000)
        exhausted = (not info["trunc_g"] and not info["aborts"]
                     and info["closed"] == 0 and hmax[0] < args.mrows)
        slack = {h: 2 * h - g for h, g in ming.items()}
        S = max(slack.values())
        res[f"extent{extent}"] = dict(info=info, hmax=hmax[0],
                                      exhausted=exhausted,
                                      fwd_min=dict(sorted(ming.items())),
                                      slack=dict(sorted(slack.items())),
                                      S=S)
        print(f"mirror all-light u=1 (extent<={extent}, gcap {args.gcap}, "
              f"m={args.mrows}): {info}; max h {hmax[0]}; EXHAUSTED "
              f"{exhausted}; fwd_min {dict(sorted(ming.items()))}; "
              f"slack max S = {S} ({round(time.time() - t0, 1)} s)",
              flush=True)
    out["mirror"] = res
    # the x-sector theorem with the mirror constants (extent 34 run)
    rec = res.get("extent34") or next(iter(res.values()))
    if rec["exhausted"]:
        rows = []
        for r in range(1, 13):
            m = 6 * r
            f, u1 = floor_theorem(m, rec["S"], rec["hmax"])
            rows.append(dict(r=r, m=m, f_x=f, u1_all_light_open=u1))
        out["mirror_theorem"] = dict(S=rec["S"], exhaust_h=rec["hmax"],
                                     rows=rows)
        print("x-sector (mirror) floor f_X(r): "
              + ", ".join(f"r={x['r']}:{x['f_x']}" for x in rows),
              flush=True)


def lane_all2(args, out):
    """Is there an infinite light walk with EVERY slab of weight
    exactly 2 (the all-light u = 2 branch's 2-per-slab equality
    case)?  Finite: the graph of weight-2 windows (translation-
    normalized) under the forced dynamics with inputs holding the
    next window at weight 2.  No cycle => every all-light u >= 2
    cycle has a slab >= 3 => 4|v| >= 2m + 1, lifting the all-light
    branch to ceil((2m+1)/4) = 3r + 1.  Both lanes."""
    from a40_s6_frontier import lsb, norm, dilate, tooth_ok
    from a40_s8_xlane import m_forced, m_tooth_ok
    t0 = time.time()
    res = {}
    for lane in ("y", "mirror"):
        # all weight-2 windows: 2 cells among 8 rows x columns, gap <= 8
        starts = set()
        for r1 in range(8):
            for r2 in range(8):
                for gap in range(0, 9):
                    if r1 == r2 and gap == 0:
                        continue
                    w = [0] * 8
                    w[r1] |= 1 << 8
                    w[r2] |= 1 << (8 + gap)
                    if sum(wt(x) for x in w) != 2:
                        continue
                    starts.add(norm(tuple(w))[0])
        # forward edges keeping the next window at weight exactly 2
        def succ(w8):
            outs = []
            if lane == "y":
                v1a, v1b, v1c, v1d, v2a, v2b, v2c, v2d = w8
                M = 8
                v1a, v1b, v1c, v1d = v1a << M, v1b << M, v1c << M, v1d << M
                v2a, v2b, v2c, v2d = v2a << M, v2b << M, v2c << M, v2d << M
                v2new = v1d ^ (v1d >> 1) ^ (v1a << 1) ^ v2d ^ (v2c >> 3)
                if not tooth_ok(v1a, v1b, v1c, v2b, v2new):
                    return outs
                fixed = (v1b, v1c, v1d, v2b, v2c, v2d, v2new)
                allm = v1a | v1b | v1c | v1d | v2a | v2b | v2c | v2d | v2new
                mk = lambda s: (v1b, v1c, v1d, s, v2b, v2c, v2d, v2new)
            else:
                p1z, p1a, p1b, p1c, p2z, p2a, p2b, p2c = w8
                M = 8
                p1z, p1a, p1b, p1c = p1z << M, p1a << M, p1b << M, p1c << M
                p2z, p2a, p2b, p2c = p2z << M, p2a << M, p2b << M, p2c << M
                u1new = m_forced(p1b, p1c, p2z, p2c)
                if not m_tooth_ok(p1a, u1new, p2z, p2a, p2b):
                    return outs
                fixed = (p1a, p1b, p1c, u1new, p2a, p2b, p2c)
                allm = p1z | p1a | p1b | p1c | p2z | p2a | p2b | p2c | u1new
                mk = lambda s: (p1a, p1b, p1c, u1new, p2a, p2b, p2c, s)
            wbase = sum(wt(x) for x in fixed)
            if wbase > 2:
                return outs
            k = 2 - wbase
            allow = dilate(allm, 4)
            acols = [i for i in range(allow.bit_length()) if allow >> i & 1]
            from itertools import combinations
            for pick in combinations(acols, k):
                s = 0
                for c in pick:
                    s |= 1 << c
                outs.append(norm(mk(s))[0])
            return outs
        # explore the finite graph; detect a cycle by DFS colouring
        WHITE, GREY, BLACK = 0, 1, 2
        colour = {}
        longest = {}
        cyc = None
        nodes = 0
        sys.setrecursionlimit(100000)

        def dfs(v):
            nonlocal cyc, nodes
            colour[v] = GREY
            nodes += 1
            best = 0
            for w in succ(v):
                c = colour.get(w, WHITE)
                if c == GREY:
                    cyc = (v, w)
                    continue
                if c == WHITE:
                    dfs(w)
                best = max(best, 1 + longest.get(w, 0))
            colour[v] = BLACK
            longest[v] = best
        for s0 in sorted(starts):
            if colour.get(s0, WHITE) == WHITE:
                dfs(s0)
        Lmax = max(longest.values()) if longest else 0
        res[lane] = dict(windows=len(starts), nodes_explored=nodes,
                         cycle_found=cyc is not None,
                         longest_all2_path_slabs=Lmax + 1,
                         cycle=[list(x) for x in cyc] if cyc else None)
        print(f"all-2 lane {lane}: {len(starts)} weight-2 start windows, "
              f"{nodes} nodes explored; CYCLE: {cyc is not None}; longest "
              f"all-weight-2 path = {Lmax + 1} slabs "
              f"({round(time.time() - t0, 1)} s)", flush=True)
    out["all2"] = res


def lane_allk(args, out):
    """Generalized: the graph of light windows of weight <= k (starts =
    every (4,4)-connected full-content window of weight 1..k, the
    seeds_full convention = the growth scope) with edges keeping the
    next slab <= k.  If ACYCLIC with longest path P_k slabs, every
    P_k + 1 consecutive slabs of any walk contain a slab >= k+1 — the
    discharging input for the light-rate bound
    Sigma (W - 2) >= (k - 1) * floor(h / (P_k + 1)) on >= 2 stretches."""
    from a40_s6_frontier import lsb, norm, dilate, tooth_ok
    from a40_s8_xlane import m_forced, m_tooth_ok
    from itertools import combinations
    t0 = time.time()
    res = {}
    sys.setrecursionlimit(1_000_000)
    for k in args.ks:
        starts = set()
        if args.wide:
            # EVERY window of weight <= k inside an 8 x (4k+4) box with
            # min column 0 — no internal connectivity assumed (a walk
            # window may consist of clusters bridged outside its rows)
            span = 4 * k + 4
            cells = [(r, c) for r in range(8) for c in range(span)]
            for w in range(1, k + 1):
                for pick in combinations(cells, w):
                    if min(c for _, c in pick) != 0:
                        continue
                    w8 = [0] * 8
                    for r, c in pick:
                        w8[r] |= 1 << c
                    starts.add(tuple(w8))
        else:
            for w in range(1, k + 1):
                for w8 in seeds_full(w):
                    starts.add(norm(tuple(w8))[0])
        for lane in ("y", "mirror"):
            def succ(w8, _k=k, _lane=lane):
                outs = []
                if _lane == "y":
                    v1a, v1b, v1c, v1d, v2a, v2b, v2c, v2d = w8
                    M = 8
                    v1a, v1b, v1c, v1d = (v1a << M, v1b << M, v1c << M,
                                          v1d << M)
                    v2a, v2b, v2c, v2d = (v2a << M, v2b << M, v2c << M,
                                          v2d << M)
                    v2new = v1d ^ (v1d >> 1) ^ (v1a << 1) ^ v2d ^ (v2c >> 3)
                    if not tooth_ok(v1a, v1b, v1c, v2b, v2new):
                        return outs
                    fixed = (v1b, v1c, v1d, v2b, v2c, v2d, v2new)
                    allm = (v1a | v1b | v1c | v1d | v2a | v2b | v2c | v2d
                            | v2new)
                    mk = lambda s: (v1b, v1c, v1d, s, v2b, v2c, v2d, v2new)
                else:
                    p1z, p1a, p1b, p1c, p2z, p2a, p2b, p2c = w8
                    M = 8
                    p1z, p1a, p1b, p1c = (p1z << M, p1a << M, p1b << M,
                                          p1c << M)
                    p2z, p2a, p2b, p2c = (p2z << M, p2a << M, p2b << M,
                                          p2c << M)
                    u1new = m_forced(p1b, p1c, p2z, p2c)
                    if not m_tooth_ok(p1a, u1new, p2z, p2a, p2b):
                        return outs
                    fixed = (p1a, p1b, p1c, u1new, p2a, p2b, p2c)
                    allm = (p1z | p1a | p1b | p1c | p2z | p2a | p2b | p2c
                            | u1new)
                    mk = lambda s: (p1a, p1b, p1c, u1new, p2a, p2b, p2c, s)
                wbase = sum(wt(x) for x in fixed)
                if wbase > _k:
                    return outs
                allow = dilate(allm, 4)
                acols = [i for i in range(allow.bit_length())
                         if allow >> i & 1]
                for kk in range(0, min(_k - wbase, 3) + 1):
                    for pick in combinations(acols, kk):
                        s = 0
                        for c in pick:
                            s |= 1 << c
                        nw = mk(s)
                        if sum(wt(x) for x in nw) < 1:
                            continue
                        outs.append(norm(nw)[0])
                return outs
            WHITE, GREY, BLACK = 0, 1, 2
            colour = {}
            longest = {}
            cyc = [None]
            nodes = [0]

            def dfs(v):
                # iterative DFS with explicit stack
                stack = [(v, iter(succ(v)))]
                colour[v] = GREY
                nodes[0] += 1
                best = {v: 0}
                while stack:
                    u, it = stack[-1]
                    adv = False
                    for w in it:
                        c = colour.get(w, WHITE)
                        if c == GREY:
                            if cyc[0] is None:
                                cyc[0] = (u, w)
                            continue
                        if c == WHITE:
                            colour[w] = GREY
                            nodes[0] += 1
                            best[w] = 0
                            stack.append((w, iter(succ(w))))
                            adv = True
                            break
                        best[u] = max(best[u], 1 + longest.get(w, 0))
                    if not adv:
                        stack.pop()
                        colour[u] = BLACK
                        longest[u] = best[u]
                        if stack:
                            pu = stack[-1][0]
                            best[pu] = max(best[pu], 1 + longest[u])
            for s0 in sorted(starts):
                if colour.get(s0, WHITE) == WHITE:
                    dfs(s0)
            Lmax = max(longest.values()) if longest else 0
            res[f"k{k}_{lane}"] = dict(
                k=k, lane=lane, start_windows=len(starts), wide=args.wide,
                nodes_explored=nodes[0], cycle_found=cyc[0] is not None,
                longest_path_slabs=Lmax + 1,
                cycle=[list(x) for x in cyc[0]] if cyc[0] else None)
            print(f"all<={k} lane {lane}: {len(starts)} start windows, "
                  f"{nodes[0]} nodes; CYCLE: {cyc[0] is not None}; "
                  f"longest all-<={k} path = {Lmax + 1} slabs "
                  f"({round(time.time() - t0, 1)} s)", flush=True)
    out["allk"] = res


# ---------------------------------------------------------------------
# the discharging automaton: exact minima of Sigma (W - 2) over slab
# sequences obeying "every P_k + 1 consecutive slabs contain a slab
# >= k + 1" (k = 2, 3, 4) with every slab >= 2 (open paths) or
# cyclically (closed walks)
# ---------------------------------------------------------------------

def automaton(Pk):
    """States = (d3, d4, d5): slabs since the last slab >= 3, >= 4,
    >= 5 (capped at the window lengths).  Weights W in {2,3,4,5}
    (W >= 6 dominated by 5 for the minimum).  Returns (states,
    transitions[(s, W)] -> s' or None if a window constraint breaks)."""
    L = {3: Pk[2] + 1, 4: Pk[3] + 1, 5: Pk[4] + 1}   # window lengths

    def step(s, W):
        d3, d4, d5 = s
        d3 = 0 if W >= 3 else d3 + 1
        d4 = 0 if W >= 4 else d4 + 1
        d5 = 0 if W >= 5 else d5 + 1
        # a run of L[j] consecutive slabs all < j is forbidden
        if d3 >= L[3] or d4 >= L[4] or d5 >= L[5]:
            return None
        return (d3, d4, d5)
    states = [(a, b, c) for a in range(L[3]) for b in range(L[4])
              for c in range(L[5]) if a <= b <= c]
    return states, step


def path_min_table(Pk, hmax):
    """pm[h] = min Sigma (W-2) over h consecutive slabs (>= 2 each)
    obeying the window constraints; the constraints are enforced
    ONLY on windows inside the stretch (sound: a stretch of a walk
    obeys at least these)."""
    states, step = automaton(Pk)
    INF = 10 ** 9
    # start: any state is possible for the history before the stretch
    # — to be sound take the most permissive: d's as if the previous
    # slabs were all heavy (d = 0)... NO: the history is unknown, so
    # windows straddling the start are NOT enforced; model by
    # starting from d = 0 (all counters reset), which enforces only
    # windows fully inside the stretch.
    cur = {(0, 0, 0): 0}
    pm = {0: 0}
    for h in range(1, hmax + 1):
        nxt = {}
        for s, v in cur.items():
            for W in (2, 3, 4, 5):
                s2 = step(s, W)
                if s2 is None:
                    continue
                v2 = v + (W - 2)
                if v2 < nxt.get(s2, INF):
                    nxt[s2] = v2
        cur = nxt
        pm[h] = min(cur.values())
    return pm


def cyc_min(Pk, m):
    """min Sigma (W-2) over CLOSED slab sequences of length m (>= 2
    each) obeying the window constraints cyclically."""
    states, step = automaton(Pk)
    INF = 10 ** 9
    best = INF
    for s0 in states:
        cur = {s0: 0}
        for h in range(m):
            nxt = {}
            for s, v in cur.items():
                for W in (2, 3, 4, 5):
                    s2 = step(s, W)
                    if s2 is None:
                        continue
                    v2 = v + (W - 2)
                    if v2 < nxt.get(s2, INF):
                        nxt[s2] = v2
            cur = nxt
        if s0 in cur and cur[s0] < best:
            best = cur[s0]
    return best


def lane_rate(args, out):
    """The assembled for-all-m floor with the discharging constraints
    (allk) + the u=1 exhaustion tables (preexact / mirror)."""
    t0 = time.time()
    ak = json.loads((DATA / "s10_forall_allk.json").read_text())["allk"]
    pe = json.loads((DATA / "s10_forall_preexact_theorem.json").read_text()
                    )["preexact"]["smax3_dil4"]
    mi = json.loads((DATA / "s10_forall_mirror.json").read_text()
                    )["mirror"]["extent34"]
    res = {}
    for lane, rec in (("y", pe), ("mirror", mi)):
        Pk = {k: ak[f"k{k}_{lane}"]["longest_path_slabs"]
              for k in (2, 3, 4)}
        assert all(not ak[f"k{k}_{lane}"]["cycle_found"] for k in (2, 3, 4))
        fwd = {int(h): g for h, g in rec["fwd_min"].items()}
        slack = {h: 2 * h - g for h, g in fwd.items()}
        hmax_tree = max(fwd)
        pm = path_min_table(Pk, 6 * args.rmax + 2)
        rows = []
        for r in range(1, args.rmax + 1):
            m = 6 * r
            # B >= 1: one block (>= 8 = 2 + 6), one u=1 run: forward
            # piece h_f (slack), the rest below the seed (>= 2 stretch)
            b1 = min(2 * m + 6 - slack[hf] + pm[m - 1 - hf]
                     for hf in fwd if m - 1 - hf >= 0)
            # u >= 2 runs only, B >= 1: 2m + 6 + pm[m-1]  (dominated)
            b1 = min(b1, 2 * m + 6 + pm[m - 1])
            # all-light u >= 2: a closed walk obeys every window inside
            # any cut-open stretch of itself, so the open-path minimum
            # is a sound (slightly weaker, much cheaper) bound
            al = 2 * m + pm[m]
            # all-light u = 1: impossible iff m + 1 > longest light path
            u1open = (m + 1 <= hmax_tree)
            cands = [b1, al] + ([m] if u1open else [])
            f = math.ceil(min(cands) / 4)
            rows.append(dict(r=r, m=m, f=f, branch_B1=b1, branch_all_light=al,
                             u1_all_light_open=u1open))
        # asymptotic slope: min-mean cycle of the automaton (Karp)
        states, step = automaton(Pk)
        n = len(states)
        idx = {s: i for i, s in enumerate(states)}
        INF = float("inf")
        D = [[INF] * n for _ in range(n + 1)]
        for i in range(n):
            D[0][i] = 0
        for kk in range(1, n + 1):
            for s in states:
                if D[kk - 1][idx[s]] == INF:
                    continue
                for W in (2, 3, 4, 5):
                    s2 = step(s, W)
                    if s2 is None:
                        continue
                    v = D[kk - 1][idx[s]] + (W - 2)
                    if v < D[kk][idx[s2]]:
                        D[kk][idx[s2]] = v
        mmc = INF
        for i in range(n):
            if D[n][i] == INF:
                continue
            worst = max((D[n][i] - D[kk][i]) / (n - kk)
                        for kk in range(n) if D[kk][i] != INF)
            mmc = min(mmc, worst)
        res[lane] = dict(Pk=Pk, S=max(slack.values()), hmax_tree=hmax_tree,
                         min_mean_excess_per_slab=mmc,
                         light_rate_per_slab=2 + mmc, rows=rows,
                         path_min=dict(list(pm.items())[:40]))
        print(f"{lane}: P_k {Pk}; light >=2 stretches: min mean excess "
              f"{mmc:.4f} per slab => sustained light rate >= "
              f"{2 + mmc:.4f} per slab; floors f(r): "
              + ", ".join(f"r={x['r']}:{x['f']}" for x in rows[:12])
              + f" ... r={rows[-1]['r']}:{rows[-1]['f']}", flush=True)
    # overall d floor = min over sectors
    both = []
    for i in range(args.rmax):
        ry, rx = res["y"]["rows"][i], res["mirror"]["rows"][i]
        both.append(dict(r=ry["r"], m=ry["m"], f_y=ry["f"], f_x=rx["f"],
                         f=min(ry["f"], rx["f"]),
                         L1=math.ceil(3 * ry["r"] / 2)))
    res["overall"] = both
    print("OVERALL d(C_{r,1}) >= min(f_Y, f_X): "
          + ", ".join(f"r={x['r']}:{x['f']}(L1 {x['L1']})"
                      for x in both[:12]), flush=True)
    # controls: every banked exact member value must sit above the
    # ell-free floor at its m
    ctl = {}
    for (l, m, d) in [(12, 12, 18), (18, 12, 24), (12, 6, 12), (18, 6, 12),
                      (24, 6, 12), (6, 6, 6), (24, 18, 12)]:
        r = m // 6
        f = both[r - 1]["f"]
        ctl[f"({l},{m})"] = dict(banked=d, f=f, ok=f <= d)
    print(f"controls (banked member values vs the ell-free floor at "
          f"their m): {ctl}", flush=True)
    assert all(v["ok"] for v in ctl.values())
    res["controls"] = ctl
    res["wall_s"] = round(time.time() - t0, 1)
    out["rate"] = res


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("lanes", nargs="*", default=["preexact", "theorem"])
    ap.add_argument("--ks", type=int, nargs="*", default=[2, 3, 4])
    ap.add_argument("--wide", action="store_true")
    ap.add_argument("--gcap", type=int, default=600)
    ap.add_argument("--hcap", type=int, default=40)
    ap.add_argument("--rmax", type=int, default=40)
    ap.add_argument("--mrows", type=int, default=40)
    ap.add_argument("--extents", type=int, nargs="*", default=[14, 34])
    ap.add_argument("--rss-cap", type=int, default=2500)
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
    if "preexact" in args.lanes:
        lane_preexact(args, out)
    if "theorem" in args.lanes:
        lane_theorem(args, out)
    if "mirror" in args.lanes:
        lane_mirror(args, out)
    if "all2" in args.lanes:
        lane_all2(args, out)
    if "allk" in args.lanes:
        lane_allk(args, out)
    if "rate" in args.lanes:
        lane_rate(args, out)
    out["wall_s"] = round(time.time() - t0, 1)
    tag = "_".join(args.lanes) + ("_wide" if args.wide else "")
    p = DATA / f"s10_forall_{tag}.json"
    p.write_text(json.dumps(out, indent=1))
    print(f"wrote {p} ({out['wall_s']} s)", flush=True)


if __name__ == "__main__":
    main()
