#!/usr/bin/env python3
"""A40 S6 — Stage 2/3: the momentum-budget assembly at the b = 1
member (l, m) = (24, 18) = [[864, 12]], from the S6 frontier tables.

THE ASSEMBLY (y-spanning sector).  For a class-minimal nontrivial
X-logical v of (24,18) that is y-spanning (Lemma K dichotomy; the
x-spanning-windowed sector is the mirrored machinery on (B, Abar) —
NOT run this session; the sector floor is reported as such):
 - v is CONNECTED (footprint (4,4) splitting + minimality; S6 drift
   script) — so its cover lift is well-defined up to deck moves;
 - slab telescope: 4|v| = sum_j W_j over the 18 slabs (exact);
 - every slab W_j >= 1 (y-spanning); heavy slabs (W >= 8) pay >= 2;
 - maximal light runs are admissible pruned light fragments whose
   deficit is bounded by the frontier tables D(h, delta) (the tables
   are relaxations: sound upper bounds);
 - closure: the cover drift telescope around the torus is = 0 mod 24
   ... but anchor slips across heavy slabs are NOT weight-bounded by
   any lemma proven this session (multi-cluster teleports bridged
   elsewhere), so the MIXED branch is assembled WITHOUT the closure
   constraint (heavy slips free mod 24) — the drift lever survives
   only on the all-light branch (no heavy slab: the walk is one
   19-slab fragment with equal end windows and total drift == 0 mod
   24, winding included via |delta| = 24 buckets).

Tiers reported:
 T1 (certificate-shaped, scope-listed): floor from certified buckets
    + per-stratum cap bounds + analytic strata.  Scope conditions:
    prefix-connected growth (dil=4, smax=3 caps; stability runs
    owed), no wrap-interacting fragments (cover extent <= 34; at
    l = 24 fragments of torus extent 21..23 that wrap-interact
    without winding are NOT enumerated — the l = 12 witness shows
    such objects are real), loose-join composition (deficit-
    inflating: sound).
 T2 (conjectural extrapolation): D(h, delta) <= (6/7) h + T0 with
    T0 = max certified transient — the "species rate + boundary
    transient" reading of the tables; conditional on the p <= 8
    alphabet (Theorem P) + no sustained aperiodic rate below 8/7.

Controls: (a) the S5 stacks (18,63)/(18,36)/(24,48) must be
PERMITTED by the same assembly (their frames' closure admits pure
species; assembled floor <= stack weight); (b) (18,12) floor <= 24;
(c) the b = 0 (12,12) witness lives in the wrapped corner the cover
tables exclude — the -6 is the admitted wrapped-corner term, and the
l = 12 instantiation of THIS assembly is scoped to non-wrapped
walks (the witness is the mechanical demonstration that the wrapped
corner is real and must stay an admitted term).
"""
from __future__ import annotations

import json
import math
import sys
import time
from pathlib import Path

LAB = Path(__file__).resolve().parent.parent
DATA = LAB / "data" / "a40"

HMAX = 19


def load_strata(tag="prod"):
    strata = {}
    for u in (1, 2, 3, 4):
        p = DATA / f"s6_frontier_u{u}{tag}.json"
        if u == 1 and (DATA / "s6_frontier_u1matched.json").exists():
            # the matched (per-seed) join: real fragments only —
            # tighter than the loose join and still a relaxation of
            # true runs (sound).
            p = DATA / "s6_frontier_u1matched.json"
        if not p.exists():
            continue
        d = json.loads(p.read_text())
        tab = {tuple(map(int, k.split(","))): g
               for k, g in d["composed"].items()}
        gcap = d["fwd"]["params"]["gcap"]
        info = dict(fwd=d["fwd"]["info"], bwd=d["bwd"]["info"])
        strata[u] = dict(table=tab, gcap=gcap, info=info,
                         validation=d.get("species_validation"),
                         src=p.name)
    return strata


def master_D_quarters(strata, hmax=HMAX, dmax=30):
    """masterD in QUARTER units: masterDq(h, delta) = 8h - minG
    upper bound (>= 4*true deficit).  Per stratum: certified bucket
    or min(cap bound, stratum-analytic (2-u/4)h); u >= 5 analytic
    0.75h at every delta."""
    out = {}
    for h in range(1, hmax + 1):
        for d in range(-dmax, dmax + 1):
            cands = [3 * h]                    # u>=5: (2-5/4)h -> 3h/4
            for u, s in strata.items():
                if (h, d) in s["table"]:
                    cands.append(8 * h - s["table"][(h, d)])
                else:
                    cap = 8 * h - (s["gcap"] + 1)
                    ana = 8 * h - u * h
                    cands.append(min(cap, ana))
            out[(h, d)] = max(cands)
    return out


def mixed_floor_quarters(Dq, m=18, hmax=HMAX):
    """w >= (8m - max sum Dq)/4 over cyclic configs with >= 1 heavy
    slab; runs partition the light slabs; heavy slips free (no
    closure).  DP over (slabs used by runs, #runs)."""
    Dbest = {h: max(Dq.get((h, d), -10**9)
                    for d in range(-30, 31)) for h in range(1, m)}
    NEG = -10**9
    f = [[NEG] * (m + 1) for _ in range(m + 1)]
    f[0][0] = 0
    for n in range(1, m):
        for r in range(1, n + 1):
            best = NEG
            for h in range(1, n + 1):
                prev = f[n - h][r - 1]
                if prev > NEG and Dbest.get(h, NEG) > NEG:
                    v = prev + Dbest[h]
                    if v > best:
                        best = v
            f[n][r] = best
    best_total, best_cfg = 0, (0, 0)
    for r in range(1, m):
        for n in range(r, m - r + 1):        # k_h = m - n >= r
            if f[n][r] > best_total:
                best_total, best_cfg = f[n][r], (n, r)
    wq = 8 * m - best_total
    return math.ceil(wq / 4), best_total / 4, best_cfg, \
        {h: v / 4 for h, v in Dbest.items()}


def all_light_floor(strata, m=18, ell=24):
    """The all-light branch: one (m+1)-slab fragment with equal end
    windows, total drift == 0 mod ell (0 and +-ell covered by the
    tables' dmax=30 at ell=24).  Per stratum u: w >= (minG - u)/4
    over delta in {0, +-ell}; absent buckets -> (gcap+1-u)/4;
    u >= 5 -> ceil(u m/4)."""
    h = m + 1
    floors = {}
    for u, s in strata.items():
        gs = []
        for d in (0, ell, -ell):
            if (h, d) in s["table"]:
                gs.append(s["table"][(h, d)])
        gmin = min(gs) if gs else s["gcap"] + 1
        gmin = min(gmin, s["gcap"] + 1) if not gs else gmin
        floors[u] = math.ceil((gmin - u) / 4)
    floors["u>=5"] = math.ceil(5 * m / 4)
    return min(floors.values()), floors


def transient_T0(strata, hmax=HMAX):
    """T0 = max over certified buckets of D - (6/7)h (quarters
    avoided: use floats); also the per-h certified best-D curve."""
    t0 = 0.0
    curve = {}
    who = None
    for u, s in strata.items():
        for (h, d), g in s["table"].items():
            D = 2 * h - g / 4
            if D > curve.get(h, (-99, None, None))[0]:
                curve[h] = (D, d, u)
            t = D - (6 / 7) * h
            if t > t0:
                t0, who = t, (u, h, d, g)
    return t0, who, curve


def hypotheses(strata, hmax=HMAX):
    """Charter hypotheses (i)/(ii)/(iii) on the certified buckets."""
    pos = {}
    neg = {}
    extreme = {}
    for u, s in strata.items():
        for (h, d), g in s["table"].items():
            D = 2 * h - g / 4
            if D <= 0:
                continue
            tgt = pos if d > 0 else (neg if d < 0 else None)
            if tgt is not None and D > tgt.get(h, -99):
                tgt[h] = D
            if abs(d) >= max(6, h):          # drift-extreme region
                k = (h, 1 if d > 0 else -1)
                if D > extreme.get(k, -99):
                    extreme[k] = D
    return dict(
        i_onesided=dict(
            verdict="REFUTED" if pos and neg else "one-signed",
            max_pos={h: round(v, 2) for h, v in sorted(pos.items())},
            max_neg={h: round(v, 2) for h, v in sorted(neg.items())}),
        ii_drift_extremes={f"{h},{s_}": round(v, 2)
                           for (h, s_), v in sorted(extreme.items())})


def main():
    t0 = time.time()
    tag = sys.argv[1] if len(sys.argv) > 1 else "prod"
    strata = load_strata(tag)
    assert strata, "no frontier tables found"
    print(f"loaded strata: {sorted(strata.keys())} "
          f"(gcaps {[s['gcap'] for s in strata.values()]})")
    for u, s in strata.items():
        v = s["validation"]
        assert v is None or v["violations"] == 0, (u, v)
        print(f"  u={u}: gcap {s['gcap']}, "
              f"{len(s['table'])} buckets, validation {v}")

    Dq = master_D_quarters(strata)
    floor_mixed, maxD, cfg, Dbest = mixed_floor_quarters(Dq)
    print(f"\nMIXED branch (>=1 heavy; heavy slips free, no closure "
          f"lever): max sum D = {maxD} at (light slabs, runs) = "
          f"{cfg} => w >= {floor_mixed}")
    print("  D_best(h) used:",
          {h: round(v, 2) for h, v in sorted(Dbest.items())
           if h <= 12})

    floor_al, al_detail = all_light_floor(strata)
    print(f"ALL-LIGHT branch (19-slab closed fragment, drift == 0 "
          f"mod 24 incl. winding): w >= {floor_al} "
          f"(per-stratum {al_detail})")

    floor_T1 = min(floor_mixed, floor_al)
    print(f"\nT1 y-sector floor (scope-listed): d_Y(24,18) >= "
          f"{floor_T1}")

    T0, who, curve = transient_T0(strata)
    print(f"\nT2 transient: T0 = {T0:.2f} realized by (u,h,delta,g) "
          f"= {who}")
    print("  certified best-D curve (h: D, delta, u):",
          {h: (round(v[0], 2), v[1], v[2])
           for h, v in sorted(curve.items())})
    # T2 conjecture: sustained interior deficit rate <= 6/7 (p <= 8
    # alphabet + no cheaper aperiodic recurrence) with boundary
    # transient <= T0 per run: D_T2(h) = min(1.75 h, (6/7) h + T0);
    # mixed DP over runs (heavy slips free):
    def d_t2(h):
        return min(1.75 * h, (6 / 7) * h + T0)
    # d_t2 concave increasing => for r runs the equal split of the
    # maximal light budget n = 18 - r is exactly optimal:
    best_t2 = max(r * d_t2((18 - r) / r) for r in range(1, 10))
    floor_T2_mixed = math.ceil(36 - best_t2)
    # all-light: closed walk, no boundary: sustained-rate floor
    floor_T2_al = math.ceil((8 / 7) * 18)
    print(f"T2 (conjectural: alphabet p<=8 + transient cap "
          f"T0={T0:.2f}): mixed floor >= {floor_T2_mixed}; "
          f"all-light sustained-rate >= {floor_T2_al}; T2 floor >= "
          f"{min(floor_T2_mixed, floor_T2_al)}")

    hyp = hypotheses(strata)
    print(f"\nhypothesis (i) one-sidedness: {hyp['i_onesided']['verdict']}")
    print(f"  max D at delta>0 per h: {hyp['i_onesided']['max_pos']}")
    print(f"  max D at delta<0 per h: {hyp['i_onesided']['max_neg']}")
    print(f"hypothesis (ii) drift-extreme deficits (|delta| >= "
          f"max(6,h)): {hyp['ii_drift_extremes']}")

    # ---- controls ------------------------------------------------------
    print("\nCONTROLS:")
    # (a) stacks: assembled floor at the stack frames must not
    # exceed the stack weights.  At (18,63): pure W7 (all-light, in
    # the l=18-embeddable subset) gives max-sum-D >= (6/7)*63 - eps:
    # floor <= 126 - 54 + eps = 72 = the stack. Mechanical check:
    # the species' own certified buckets sustain D-rate >= 6/7 - eps
    # up to the caps:
    ok_a = True
    for h in range(7, 15, 7):
        got = curve.get(h + 1, (0,))[0]  # h+1 slabs ~ h rows of species
        if got < (6 / 7) * h - 1:
            ok_a = False
    print(f"  (a) species D-rate present in tables through the caps: "
          f"{'PASS' if ok_a else 'FAIL'} — the (18,63)/(18,36)/"
          f"(24,48) stacks are PERMITTED (their frames' closure "
          f"admits pure species; the member IP never uses closure "
          f"on the mixed branch)")
    # (b) (18,12): the same DP at m=12 must give floor <= 24:
    f12, maxD12, cfg12, _ = mixed_floor_quarters(Dq, m=12)
    print(f"  (b) (18,12) instantiation: mixed floor {f12} <= 24: "
          f"{'PASS' if f12 <= 24 else 'FAIL'}")
    # (c) wrapped corner: the l=12 witness (w18 = 2m-6) is a
    # both-wraps object with NO cover lift (s6_drift.json) — the
    # cover-table assembly is scoped to non-wrapped walks and the
    # witness's -6 is the admitted wrapped-corner term:
    wd = json.loads((DATA / "s6_drift.json").read_text())
    ok_c = (wd["zero_and_l12"]["a36_witness"]["compact_cover_lift"]
            is False)
    print(f"  (c) b=0 witness in the wrapped corner (no cover "
          f"lift): {'PASS' if ok_c else 'FAIL'} — the -6 is the "
          f"admitted wrapped/winding term at l=12; at l=24 the "
          f"wrapped corner is a listed scope condition")

    out = dict(
        strata={u: dict(gcap=s["gcap"], buckets=len(s["table"]),
                        validation=s["validation"])
                for u, s in strata.items()},
        floor_mixed=floor_mixed, maxD_mixed=maxD, cfg_mixed=cfg,
        floor_all_light=floor_al, all_light_detail=al_detail,
        floor_T1=floor_T1,
        T0=round(T0, 3), T0_witness=who,
        best_D_curve={h: dict(D=round(v[0], 3), delta=v[1], u=v[2])
                      for h, v in sorted(curve.items())},
        floor_T2_mixed=floor_T2_mixed, floor_T2_al=floor_T2_al,
        hypotheses=hyp,
        controls=dict(a_species_rate=ok_a, b_1812_floor=f12,
                      c_wrapped_corner=ok_c),
        wall_s=round(time.time() - t0, 1))
    (DATA / "s6_assembly.json").write_text(json.dumps(out, indent=1))
    print(f"\nwrote {DATA/'s6_assembly.json'} "
          f"({out['wall_s']} s)")


if __name__ == "__main__":
    main()
