#!/usr/bin/env python3
"""A17 E20 — the P-33 zone/free endgame: finite decisive checks.

ZONES (both-side torsion-confined terminals, exponent d̃): every
class-frame realization lies in T_d̃(G) = Z_{gcd(d̃,ℓ)} × Z_{gcd(d̃,m)}
which embeds in Z_d̃² (injective hom; all violation checks are
equality patterns, preserved both ways). Z_d̃² = Z_d̃ ⊕ Z_d̃
componentwise, so the solution set of the terminal's lattice in
Z_d̃² is S × S with S = {x ∈ Z_d̃⁸ : Lx ≡ 0}, enumerated via the
Smith diagonalization (x = Vz, diag·z ≡ 0). A pair (x, x′) is a
VALID realization iff no battery form vanishes on it — i.e. iff
W(x) ∩ W(x′) = ∅ where W = the set of battery forms vanishing mod
d̃ (d̃ odd ⟹ the d = 2 torsion caveat collapses to vanishing).
Exhausting S × S with empty intersection everywhere KILLS the zone
family on every class frame.

FREE (rank-7 terminals): if the Smith diagonal is (1,…,1,0), the
solution set in ANY group is the line x = t·v (v = the kernel
column of V, primitive). Then A, B ⊆ ⟨t⟩ and everything reduces to
integers: a battery form f violates iff N | f·v (N = ord t; eq
forms) or N | 2(f·v) (difference-class forms). (iii) enters through
the projections: with n_x = ord(t_x), n_y = ord(t_y) (n_x | ℓ,
n_y | m, 4∤ both), A's x-pattern = multiset of {0,v₁..v₄} mod n_x
must be MONO ({5},{1,4},{2,3},{1,2,2}), A's y-pattern mod n_y must
be ANTI ((1,1,3),(1,1,1,2),(1⁵)), and mirrored for B. Concentration
of 5 values into ≤ 3 classes forces n | (a pairwise difference) or
n ≤ 3, so the candidate (n_x, n_y) set is FINITE and explicit. An
empty search proves (iii)-necessity for the family; a survivor is
an explicit counterexample candidate (report it).

Usage: uv run python scripts/a17_e20_p33_endgame.py \
    --table data/a17/e20_p33_table.json \
    --families data/a17/e20_p33_families.json \
    --out data/a17/e20_p33_endgame.json
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from collections import Counter
from itertools import combinations
from math import gcd, lcm
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from a17_e20_overlay_engine import (ZLattice, avec, bvec,
                                    p33_grid_forms, static_battery)
from a17_e20_p33_classify import reduce_torsion


def smith_diag(M: list, m: int = 8) -> tuple:
    """Diagonalize M by unimodular row/col ops; return (diag, V)
    with x = V·z and (row ops)·M·V = diag. diag entries ≥ 0."""
    A = [list(r) for r in M]
    if not A:
        A = [[0] * m]
    n = len(A)
    V = [[1 if i == j else 0 for j in range(m)] for i in range(m)]
    for t in range(m):
        while True:
            piv = None
            for i in range(t, n):
                for j in range(t, m):
                    if A[i][j] and (piv is None or
                                    abs(A[i][j]) < abs(A[piv[0]][piv[1]])):
                        piv = (i, j)
            if piv is None:
                break
            i0, j0 = piv
            A[t], A[i0] = A[i0], A[t]
            if j0 != t:
                for r_ in A:
                    r_[t], r_[j0] = r_[j0], r_[t]
                for r_ in V:
                    r_[t], r_[j0] = r_[j0], r_[t]
            clean = True
            for i in range(t + 1, n):
                if A[i][t]:
                    q = A[i][t] // A[t][t]
                    A[i] = [x - q * y for x, y in zip(A[i], A[t])]
                    if A[i][t]:
                        clean = False
            for j in range(t + 1, m):
                if A[t][j]:
                    q = A[t][j] // A[t][t]
                    for r_ in A:
                        r_[j] -= q * r_[t]
                    for r_ in V:
                        r_[j] -= q * r_[t]
                    if A[t][j]:
                        clean = False
            if clean:
                break
        # next t
    diag = []
    for i in range(m):
        diag.append(abs(A[i][i]) if i < n else 0)
    return diag, V


def solutions_mod(diag, V, D: int):
    """All x = V z mod D with diag_i z_i ≡ 0 mod D."""
    gs = [gcd(d, D) for d in diag]
    steps = [D // g for g in gs]
    sols = []

    def rec(i, z):
        if i == 8:
            x = tuple(sum(V[r][c] * z[c] for c in range(8)) % D
                      for r in range(8))
            sols.append(x)
            return
        for k in range(gs[i]):
            z[i] = k * steps[i]
            rec(i + 1, z)
        z[i] = 0

    rec(0, [0] * 8)
    return sols




def collision_forms(t: dict) -> list:
    """The 3+3 pairwise-distinctness forms of the collision cells,
    from the terminal's witness record (steps order:
    D01, D02, D12 (B-side), G10, G20, G21 (A-side)). Vanishing =
    coinciding collision cells = a non-grid sigma — violation."""
    from a17_e20_overlay_engine import p33_geometry, vadd, vsub
    r, g, _ = p33_geometry(t["alpha"], t["beta"])
    w = t["wits"]
    # A-side overlap cells: G10: g1-g0 = a_p - a_q -> cell a_p + g0;
    # G20: cell a_p' + g0; G21: g2-g1 = a_p'' - a_q'' -> a_p'' + g1
    vA = [avec(w[3][1]), avec(w[4][1]),
          vadd(avec(w[5][1]), g[1])]
    # B-side: D01: r0-r1 = b_s-b_t -> cell b_t + r0; D02: b_t' + r0;
    # D12: r1-r2 = b_s''-b_t'' -> cell b_t'' + r1
    vB = [bvec(w[0][2]), bvec(w[1][2]),
          vadd(bvec(w[2][2]), r[1])]
    forms = []
    for u, v_ in combinations(vA, 2):
        forms.append(vsub(u, v_))
    for u, v_ in combinations(vB, 2):
        forms.append(vsub(u, v_))
    return forms


def battery_for(t: dict) -> tuple:
    stat = static_battery()
    gforms, gkinds, _, _ = p33_grid_forms(t["alpha"], t["beta"])
    cf = collision_forms(t)
    return (stat.forms + gforms + cf,
            stat.kinds + gkinds + ["eq"] * len(cf))


MONO = {(5,), (1, 4), (2, 3), (1, 2, 2)}
ANTI = {(1, 1, 3), (1, 1, 1, 2), (1, 1, 1, 1, 1)}


def pattern(vals, n: int) -> tuple:
    return tuple(sorted(Counter(v % n for v in vals).values()))


def divisors(x: int) -> set:
    x = abs(x)
    out = set()
    d = 1
    while d * d <= x:
        if x % d == 0:
            out.add(d)
            out.add(x // d)
        d += 1
    return out


def zone_endgame(zones: list, cap: int) -> dict:
    import numpy as np
    out = {"families": 0, "skipped_big": 0, "survivor_pairs": 0,
           "examples": [], "sol_hist": Counter()}
    for t, dt in zones:
        lat = ZLattice()
        lat.rows = [r[:] for r in t["lattice"]]
        diag, V = smith_diag(lat.rows)
        S = solutions_mod(diag, V, dt)
        # sanity: every solution satisfies the lattice mod dt
        for x in S[:5]:
            for row in lat.rows:
                assert sum(c * v for c, v in zip(row, x)) % dt == 0
        out["sol_hist"][len(S)] += 1
        if len(S) ** 2 > cap:
            out["skipped_big"] += 1
            continue
        forms, kinds = battery_for(t)
        X = np.array(S, dtype=np.int64)
        FM = np.array([list(f) for f in forms], dtype=np.int64)
        Z = (X @ FM.T) % dt == 0
        W = [frozenset(np.nonzero(Z[i])[0].tolist())
             for i in range(len(S))]
        n_valid = 0
        n_iii = 0
        for i in range(len(S)):
            for j in range(len(S)):
                if not (W[i] & W[j]):
                    n_valid += 1
                    VA = [0, *S[i][:4]]
                    VAp = [0, *S[j][:4]]
                    VB = [0, *S[i][4:]]
                    VBp = [0, *S[j][4:]]
                    for ax, ay, bx, by in ((VA, VAp, VB, VBp),
                                           (VAp, VA, VBp, VB)):
                        if (pattern(ax, dt) in MONO
                                and pattern(ay, dt) in ANTI
                                and pattern(by, dt) in MONO
                                and pattern(bx, dt) in ANTI):
                            n_iii += 1
                            if len(out["examples"]) < 5:
                                out["examples"].append({
                                    "dt": dt, "x": list(S[i]),
                                    "xp": list(S[j]),
                                    "alpha": t["alpha"],
                                    "beta": t["beta"],
                                    "iii": True})
                            break
        out["survivor_pairs"] += n_valid
        out["iii_pairs"] = out.get("iii_pairs", 0) + n_iii
        out["families"] += 1
    out["sol_hist"] = dict(sorted(out["sol_hist"].items()))
    return out


def free_endgame(frees: list) -> dict:
    out = {"families": 0, "pure_line": 0, "mixed_torsion": 0,
           "iii_survivors": [], "v_profiles": Counter()}
    for t in frees:
        lat = ZLattice()
        lat.rows = [r[:] for r in t["lattice"]]
        diag, V = smith_diag(lat.rows)
        nz = [d for d in diag if d != 0]
        zero_idx = [i for i, d in enumerate(diag) if d == 0]
        out["families"] += 1
        mixed = len(zero_idx) != 1 or any(d != 1 for d in nz)
        if mixed:
            out["mixed_torsion"] += 1
            prof = tuple(sorted(d for d in nz if d != 1))
            k = "torsion_profiles"
            out.setdefault(k, {})
            out[k][str((len(zero_idx), prof))] = \
                out[k].get(str((len(zero_idx), prof)), 0) + 1
            if len(zero_idx) != 1:
                out["multi_free"] = out.get("multi_free", 0) + 1
                continue
            # fall through: analyze the PURE-LINE (w0 = 0) sub-case;
            # torsion dressings (w0 != 0, e in {2,3}) remain pending
        else:
            out["pure_line"] += 1
        v = [V[r][zero_idx[0]] for r in range(8)]
        g = 0
        for x in v:
            g = gcd(g, x)
        v = [x // g for x in v]
        # sanity: v generates the kernel
        for row in lat.rows:
            assert sum(c * x for c, x in zip(row, v)) == 0
        VA = [0] + v[:4]
        VB = [0] + v[4:]
        out["v_profiles"][tuple(v)] += 1
        forms, kinds = battery_for(t)
        vals = [(sum(c * x for c, x in zip(f, v)), k)
                for f, k in zip(forms, kinds)]
        if any(c == 0 for c, _ in vals):
            # the pure line is dead (an in-span eq-form fires).
            # Dressed layer 1: x = t·v + w0·u, w0 in G[e] \ 0
            # (e = the torsion divisor, u = its V-column). Every
            # line-vanishing form f takes the value (f·u)·w0:
            #   e = 2: diff-kind f -> value is 0 or 2-torsion, BOTH
            #     violations -> family dead for every dressing;
            #     eq-kind with f·u even -> value 0 -> dead;
            #   e = 3: f·u ≡ 0 (mod 3) -> value 0 -> dead;
            #     otherwise the form survives all w0 ≠ 0.
            out["line_dressed"] = out.get("line_dressed", 0) + 1
            e = max(nz)
            i_e = diag.index(e)
            u = [V[r][i_e] for r in range(8)]
            vanish = [(f, k) for (f, k), (c, _) in
                      zip(zip(forms, kinds), vals) if c == 0]
            dead = False
            for f, k in vanish:
                fu = sum(c * x for c, x in zip(f, u)) % e
                if e == 2 and (k == "diff" or fu == 0):
                    dead = True
                    break
                if e == 3 and fu == 0:
                    dead = True
                    break
            key = (f"dressed_dead_e{e}" if dead
                   else f"dressed_open_e{e}")
            out[key] = out.get(key, 0) + 1
            if not dead:
                out.setdefault("dressed_open", []).append({
                    "e": e, "v": v, "u": u,
                    "n_vanish": len(vanish),
                    "alpha": t["alpha"], "beta": t["beta"],
                    "wits": t["wits"]})
            continue
        out["line_analyzed"] = out.get("line_analyzed", 0) + 1
        cand_x = {1, 2, 3}
        for p, q in combinations(VA, 2):
            cand_x |= {n for n in divisors(p - q) if n % 4 != 0}
        cand_y = {1, 2, 3}
        for p, q in combinations(VB, 2):
            cand_y |= {n for n in divisors(p - q) if n % 4 != 0}
        for nx in sorted(cand_x):
            if pattern(VA, nx) not in MONO:
                continue
            for ny in sorted(cand_y):
                if pattern(VB, ny) not in MONO:
                    continue
                if pattern(VA, ny) not in ANTI:
                    continue
                if pattern(VB, nx) not in ANTI:
                    continue
                N = lcm(nx, ny)
                bad = False
                for c, k in vals:
                    m = 2 * c if k == "diff" else c
                    if m % N == 0:
                        bad = True
                        break
                if not bad:
                    out["iii_survivors"].append({
                        "v": v, "nx": nx, "ny": ny, "N": N,
                        "alpha": t["alpha"], "beta": t["beta"]})
    out["v_profiles"] = {str(k): c for k, c in
                         sorted(out["v_profiles"].items(),
                                key=lambda x: -x[1])[:10]}
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--table", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--cap", type=int, default=5_000_000)
    args = ap.parse_args()
    t0 = time.time()
    table = json.loads(Path(args.table).read_text())
    zones, frees, r1 = [], [], 0
    for t in table["terminal_patterns"]:
        lat = ZLattice()
        lat.rows = [r[:] for r in t["lattice"]]
        forms, kinds = battery_for(t)
        killed = False
        for f, k in zip(forms, kinds):
            if k != "diff":
                continue
            d = lat.forced_denom(tuple(f))
            if d is not None and (d & (d - 1)) == 0:
                killed = True
                break
        if killed:
            r1 += 1
            continue
        da = [lat.forced_denom(avec(i)) for i in range(1, 5)]
        db = [lat.forced_denom(bvec(i)) for i in range(1, 5)]
        if all(x is not None for x in da) and \
                all(x is not None for x in db):
            dt = reduce_torsion(lcm(*[x for x in da + db]))
            zones.append((t, dt))
        else:
            frees.append(t)
    print(f"inline classify: R1={r1}, zones={len(zones)}, "
          f"frees={len(frees)}", flush=True)
    zres = zone_endgame(zones, args.cap)
    fres = free_endgame(frees)
    out = {"zone": zres, "free": fres,
           "secs": round(time.time() - t0, 1)}
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    Path(args.out).write_text(json.dumps(out, indent=1))
    print(json.dumps(out, indent=1))
    zdead = zres["survivor_pairs"] == 0 and zres["skipped_big"] == 0
    fdead = (len(fres["iii_survivors"]) == 0
             and fres["mixed_torsion"] == 0)
    print(f"\nZONES: {'ALL DEAD' if zdead else 'RESIDUALS'} "
          f"({zres['families']} checked); "
          f"FREE: {'ALL (iii)-DEAD' if fdead else 'RESIDUALS'} "
          f"({fres['pure_line']} pure lines, "
          f"{fres['mixed_torsion']} mixed-torsion pending)")


if __name__ == "__main__":
    main()
