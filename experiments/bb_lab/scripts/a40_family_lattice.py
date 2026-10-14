#!/usr/bin/env python3
"""A40 P1 — the tour-de-gross family's internal cover lattice, mechanical.

Primary source (arXiv:2506.03094v1, "Future directions" item 3): members are
(r, b) with r >= 1 integer, b in {0,1} a BIT; group Z_{6(r+b)} x Z_{6r};
FIXED Laurent polynomials shared by every member,

    A = 1 + y + x^3 y^-1      support {(0,0),(0,1),(3,-1)}
    B = 1 + x + x^-1 y^-3     support {(0,0),(1,0),(-1,-3)}

conjectured d = 6(2r+b-1); n = 72 r (r+b); k = 12 ("we find", numerical).

This script, after the falsify-first gate (charter A38 s6.0):

  P1a  builds every member with r <= RMAX, asserts k = 12, and asserts the
       paper presentation is unit-shift / Aut-equivalent to the repo's
       stored presentations (bb72, gross, BCGMRY two-gross);
  P1b  decides every ordered member pair mechanically:
         - covering degree n'/n integer?  (else NO cover of any kind)
         - does ANY subgroup K <= G' give G'/K iso G?  (p-part criterion +
           brute-force cross-check on small groups)
         - is the AXIS-ALIGNED projection a LITERAL polynomial lift?
           (fixed Laurent supports reduce with no collisions)
         - deck structure Z_{l'/l} x Z_{m'/m}, 2-part / odd part split;
  P1c  screens the free-Z2 sub-steps of the two known edges against the
       banked facts (k rows, (R), regime data) — validation edges;
  P1d  each member's own tower_inventory (maximal iterated-Z2 chain to an
       odd bottom, which intermediates are family members vs not);
  P1e  the A13 x-ladder correction: (24,6) literal lift has k = 12 but is
       NOT a family member, and A14 s13's SAT@12 witness makes it a
       provably different code from the (12,12) two-gross.

Output: data/a40/family_lattice.json (+ stdout tables).
"""
from __future__ import annotations

import json
import math
import sys
import time
from itertools import product
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from bb_lab.tower import (  # noqa: E402
    TowerCode, AxisDeck, screen_rung, tower_inventory, fold_support,
    support_str, validate_banked,
)

DATA = ROOT / "data" / "a40"
DATA.mkdir(parents=True, exist_ok=True)

RMAX = 4                      # members through (4,1) = n 1440
EXTRA = [(5, 0)]              # + (5,0) n=1800 for the diagonal chain

# The paper's fixed Laurent supports (exponent tuples of (x, y)).
A_LAUR = [(0, 0), (0, 1), (3, -1)]
B_LAUR = [(0, 0), (1, 0), (-1, -3)]


def member_group(r: int, b: int) -> tuple[int, int]:
    return (6 * (r + b), 6 * r)


def reduce_supp(supp, lm):
    out = frozenset((e[0] % lm[0], e[1] % lm[1]) for e in supp)
    return out


def shift_supp(supp, t, lm):
    return frozenset(((e[0] + t[0]) % lm[0], (e[1] + t[1]) % lm[1])
                     for e in supp)


def aut_y(supp, u, lm):
    """Group automorphism y -> y^u (u invertible mod m)."""
    assert math.gcd(u, lm[1]) == 1
    return frozenset((e[0], (e[1] * u) % lm[1]) for e in supp)


# ---------------------------------------------------------- group theory
def p_parts(a: int, b: int) -> dict[int, tuple[int, int]]:
    """Invariant p-exponent pairs (hi >= lo) of Z_a x Z_b, per prime."""
    out: dict[int, tuple[int, int]] = {}
    for p in {q for n in (a, b) for q in _primes(n)}:
        ea = _vp(a, p)
        eb = _vp(b, p)
        out[p] = (max(ea, eb), min(ea, eb))
    return out


def _primes(n: int) -> set[int]:
    s, d = set(), 2
    while d * d <= n:
        while n % d == 0:
            s.add(d)
            n //= d
        d += 1
    if n > 1:
        s.add(n)
    return s


def _vp(n: int, p: int) -> int:
    v = 0
    while n % p == 0:
        n //= p
        v += 1
    return v


def has_quotient(big: tuple[int, int], small: tuple[int, int]) -> bool:
    """Does Z_{big} admit ANY quotient iso Z_{small}?  For finite abelian
    groups, quotients of G' = subgroups of G' (duality); for 2-generated
    p-groups Z_{p^a} x Z_{p^b} (a>=b) contains Z_{p^c} x Z_{p^d} (c>=d)
    iff c <= a and d <= b; per-prime."""
    if (big[0] * big[1]) % (small[0] * small[1]) != 0:
        return False
    P, Q = p_parts(*big), p_parts(*small)
    for p, (c, d) in Q.items():
        a, b = P.get(p, (0, 0))
        if c > a or d > b:
            return False
    return True


def brute_quotient_check(big: tuple[int, int], small: tuple[int, int]) -> bool:
    """Brute force (|big| <= 300): enumerate 2-generated subgroups K and
    test G'/K iso small via element-order census."""
    a, b = big
    n = a * b
    deg = n // (small[0] * small[1])
    if n % (small[0] * small[1]):
        return False
    els = [(i, j) for i in range(a) for j in range(b)]

    def closure(gens):
        K = {(0, 0)}
        frontier = list(gens)
        while frontier:
            g = frontier.pop()
            if g in K:
                continue
            K.add(g)
            for h in list(K):
                s = ((g[0] + h[0]) % a, (g[1] + h[1]) % b)
                if s not in K:
                    frontier.append(s)
        return frozenset(K)

    def quot_orders(K):
        # order census of G'/K
        reps: dict[frozenset, int] = {}
        seen = set()
        census: dict[int, int] = {}
        for g in els:
            cs = frozenset(((g[0] + k[0]) % a, (g[1] + k[1]) % b) for k in K)
            if cs in seen:
                continue
            seen.add(cs)
            o = 1
            x = g
            while x not in K:
                x = ((x[0] + g[0]) % a, (x[1] + g[1]) % b)
                o += 1
                if o > n:
                    raise RuntimeError
            census[o] = census.get(o, 0) + 1
        return census

    def group_orders(lm):
        census: dict[int, int] = {}
        for i in range(lm[0]):
            for j in range(lm[1]):
                # order of (i,j) = lcm(ord_i, ord_j)
                oi = lm[0] // math.gcd(i, lm[0])
                oj = lm[1] // math.gcd(j, lm[1])
                o = oi * oj // math.gcd(oi, oj)
                census[o] = census.get(o, 0) + 1
        return census

    target = group_orders(small)
    tried = set()
    for g1 in els:
        for g2 in els:
            K = closure([g1, g2])
            if len(K) != deg or K in tried:
                continue
            tried.add(K)
            if quot_orders(K) == target:
                return True
    return False


# ------------------------------------------------------------------ main
def main() -> None:
    t0 = time.time()
    out: dict = {"when": time.strftime("%Y-%m-%d %H:%M:%S")}

    # ---------- gate (falsify-first, charter s6.0)
    print("== gate: bb_lab.tower.validate_banked() ==")
    gate = validate_banked(ROOT / "data", log=lambda s: print(s))
    out["gate"] = "PASS"
    print(f"   gate PASS ({time.time()-t0:.1f} s)")

    # ---------- P1a: members
    print("\n== P1a: family members (fixed Laurent polynomials) ==")
    members: dict[tuple[int, int], dict] = {}
    codes: dict[tuple[int, int], TowerCode] = {}
    pairs_rb = [(r, b) for r in range(1, RMAX + 1) for b in (0, 1)] + EXTRA
    for (r, b) in pairs_rb:
        lm = member_group(r, b)
        As = reduce_supp(A_LAUR, lm)
        Bs = reduce_supp(B_LAUR, lm)
        assert len(As) == 3 and len(Bs) == 3, f"support collision at {lm}"
        code = TowerCode(f"tdg({r},{b})", lm, As, Bs)
        codes[(r, b)] = code
        d_conj = 6 * (2 * r + b - 1)
        assert code.n == 72 * r * (r + b)
        assert code.k == 12, f"k({lm}) = {code.k} != 12"
        members[(r, b)] = {
            "lm": list(lm), "n": code.n, "k": code.k, "d_conj": d_conj,
            "A": support_str(As), "B": support_str(Bs),
            "inventory": tower_inventory(lm),
        }
        inv = members[(r, b)]["inventory"]
        print(f"  ({r},{b}) Z{lm[0]}xZ{lm[1]:<3} [[{code.n},{code.k}]] "
              f"d_conj={d_conj:<3} v2={inv['v2_per_axis']} "
              f"odd={inv['odd_part_per_axis']} depth={inv['depth']}")

    # presentation checks vs the repo's stored forms
    print("\n  presentation checks (paper Laurent -> stored forms):")
    # bb72 (6,6) and gross (12,6): y*A, x*B unit shifts
    for (r, b), stored_A, stored_B in [
        ((1, 0), {(3, 0), (0, 1), (0, 2)}, {(1, 0), (2, 0), (0, 3)}),
        ((1, 1), {(3, 0), (0, 1), (0, 2)}, {(1, 0), (2, 0), (0, 3)}),
    ]:
        lm = member_group(r, b)
        As = reduce_supp(A_LAUR, lm)
        Bs = reduce_supp(B_LAUR, lm)
        assert shift_supp(As, (0, 1), lm) == frozenset(stored_A), (r, b)
        assert shift_supp(Bs, (1, 0), lm) == frozenset(stored_B), (r, b)
        print(f"    ({r},{b}): y*A == {support_str(stored_A)}, "
              f"x*B == {support_str(stored_B)}  [unit shifts]  OK")
    # two-gross (12,12): Aut y->y^7 then unit shifts == BCGMRY stored form
    lm = member_group(2, 0)
    As = reduce_supp(A_LAUR, lm)
    Bs = reduce_supp(B_LAUR, lm)
    A_bc = frozenset({(3, 0), (0, 2), (0, 7)})   # x^3 + y^2 + y^7
    B_bc = frozenset({(0, 3), (1, 0), (2, 0)})   # y^3 + x + x^2
    assert shift_supp(aut_y(As, 7, lm), (0, 7), lm) == A_bc
    assert shift_supp(aut_y(Bs, 7, lm), (1, 0), lm) == B_bc
    print("    (2,0): (y->y^7; y^7*A, x*B) == BCGMRY (x^3+y^2+y^7, "
          "y^3+x+x^2)  [Aut + unit shifts]  OK")
    # and the stored form's k as a cross-check
    tg_stored = TowerCode("two-gross/stored", lm, "x^3 + y^2 + y^7",
                          "y^3 + x + x^2")
    assert tg_stored.k == 12
    out["members"] = {f"({r},{b})": v for (r, b), v in members.items()}

    # ---------- P1b: the pairwise cover lattice
    print("\n== P1b: pairwise cover relations ==")
    edges = []
    keys = sorted(members, key=lambda t: (members[t]["n"], t))
    small_brute = 0
    for big in keys:
        for small in keys:
            if big == small:
                continue
            nb, ns = members[big]["n"], members[small]["n"]
            if nb <= ns:
                continue
            LMb = tuple(members[big]["lm"])
            LMs = tuple(members[small]["lm"])
            row: dict = {"cover": f"({big[0]},{big[1]})",
                         "base": f"({small[0]},{small[1]})",
                         "lm_cover": list(LMb), "lm_base": list(LMs)}
            if nb % ns:
                row["verdict"] = "IMPOSSIBLE (non-integer degree)"
                row["degree"] = round(nb / ns, 3)
                edges.append(row)
                continue
            deg = nb // ns
            row["degree"] = deg
            anyq = has_quotient(LMb, LMs)
            if (LMb[0] * LMb[1]) <= 300:
                bq = brute_quotient_check(LMb, LMs)
                assert bq == anyq, (LMb, LMs, bq, anyq)
                small_brute += 1
            row["group_quotient_exists"] = anyq
            axis_ok = LMb[0] % LMs[0] == 0 and LMb[1] % LMs[1] == 0
            row["axis_aligned"] = axis_ok
            if not anyq:
                row["verdict"] = "IMPOSSIBLE (no group quotient)"
                edges.append(row)
                continue
            if axis_ok:
                # literal-lift test: fixed Laurent supports reduce cleanly
                Ab = reduce_supp(A_LAUR, LMb)
                Bb = reduce_supp(B_LAUR, LMb)
                As2 = reduce_supp(A_LAUR, LMs)
                Bs2 = reduce_supp(B_LAUR, LMs)
                foldA = reduce_supp(Ab, LMs)
                foldB = reduce_supp(Bb, LMs)
                lit = (foldA == As2 and foldB == Bs2 and
                       len(foldA) == 3 and len(foldB) == 3)
                row["literal_lift"] = bool(lit)
                dx, dy = LMb[0] // LMs[0], LMb[1] // LMs[1]
                row["deck"] = f"Z{dx} x Z{dy} (x,y)"
                row["deck_two_part"] = (dx // _odd(dx)) * (dy // _odd(dy))
                row["deck_odd_part"] = _odd(dx) * _odd(dy)
                row["verdict"] = ("COVER (literal lift)" if lit else
                                  "axis-aligned but NOT literal")
            else:
                row["verdict"] = ("group quotient exists but NOT "
                                  "axis-aligned (no literal lift of the "
                                  "fixed polynomials)")
            edges.append(row)
    out["edges"] = edges
    print(f"  {len(edges)} ordered pairs; brute-force quotient "
          f"cross-checks on {small_brute} small pairs all agree")
    for row in edges:
        if "COVER" in row["verdict"]:
            print(f"  COVER  {row['cover']} -> {row['base']}  deg "
                  f"{row['degree']}  deck {row['deck']}  "
                  f"2-part {row['deck_two_part']}  odd {row['deck_odd_part']}")
    n_cover = sum(1 for r in edges if "COVER" in r["verdict"])
    n_imposs = sum(1 for r in edges if "IMPOSSIBLE" in r["verdict"])
    print(f"  totals: {n_cover} cover edges, {n_imposs} impossible, "
          f"{len(edges)-n_cover-n_imposs} other")

    # consecutive members (the zigzag): which steps are covers?
    print("\n  consecutive zigzag steps:")
    zig = sorted(members, key=lambda t: members[t]["n"])
    zsteps = []
    for i in range(len(zig) - 1):
        a, c = zig[i], zig[i + 1]
        row = next(r for r in edges
                   if r["cover"] == f"({c[0]},{c[1]})"
                   and r["base"] == f"({a[0]},{a[1]})")
        zsteps.append({"step": f"{a} -> {c}", "verdict": row["verdict"],
                       "degree": row["degree"]})
        print(f"    {a} -> {c}: {row['verdict']} (deg {row['degree']})")
    out["zigzag"] = zsteps

    # ---------- P1c: screen the two known edges (validation)
    print("\n== P1c: known-edge validation (structure screens) ==")
    val = {}
    # bb72 -> gross is the x-axis Z2 fold of (12,6) onto (6,6)
    rs = screen_rung(codes[(1, 1)], codes[(1, 0)], 0, 16, do_fibers=False)
    assert rs["k_cover"] == 12 and rs["k_base"] == 12 and rs["R_holds"]
    val["bb72->gross"] = {k: v for k, v in rs.items()
                          if not k.startswith("_")}
    print(f"  (1,0)->(1,1) x-fold: k 12->12, (R) {rs['R_holds']}, "
          f"twisted {rs['twisted']}, exact(cov/base) "
          f"{rs['exact_cover']}/{rs['exact_base']}, sigma*=id "
          f"{rs['sigma_id']}  [banked: perfect doubling, d 6->12, Lean]")
    # gross -> two-gross is the y-axis Z2 fold of (12,12) onto (12,6)
    rs2 = screen_rung(codes[(2, 0)], codes[(1, 1)], 1, 16, do_fibers=False)
    assert rs2["k_cover"] == 12 and rs2["k_base"] == 12 and rs2["R_holds"]
    val["gross->two-gross"] = {k: v for k, v in rs2.items()
                               if not k.startswith("_")}
    print(f"  (1,1)->(2,0) y-fold: k 12->12, (R) {rs2['R_holds']}, "
          f"twisted {rs2['twisted']}, exact(cov/base) "
          f"{rs2['exact_cover']}/{rs2['exact_base']}, sigma*=id "
          f"{rs2['sigma_id']}  [banked: deficit 6, d 12->18 A36 cert]")
    out["known_edges"] = val

    # ---------- P1e: the A13 x-ladder correction
    print("\n== P1e: the x-ladder (24,6) is not the two-gross ==")
    x2 = TowerCode("xladder/L2", (24, 6),
                   reduce_supp([(3, 0), (0, 1), (0, 2)], (24, 6)),
                   reduce_supp([(1, 0), (2, 0), (0, 3)], (24, 6)))
    assert x2.k == 12   # A13 T3: k = 12 along the whole x-ladder
    sol = [(r, b) for r in range(1, 9) for b in (0, 1)
           if member_group(r, b) == (24, 6)]
    assert sol == [], "(24,6) solves the member equations?!"
    out["xladder_check"] = {
        "lm": [24, 6], "n": x2.n, "k": x2.k,
        "family_member": False,
        "note": "A14 s13 battery: SAT witness at weight 12 -> d_X <= 12 "
                "(witness grade); two-gross has d = 18 (A36 certificate) "
                "-> provably different codes",
    }
    print(f"  (24,6) literal x-ladder L2: [[{x2.n},{x2.k}]], NOT a family "
          f"member (b=3 required); d_X <= 12 (A14 s13 witness) != 18")

    out["wall_s"] = round(time.time() - t0, 1)
    (DATA / "family_lattice.json").write_text(json.dumps(out, indent=1))
    print(f"\nwrote {DATA/'family_lattice.json'}  ({out['wall_s']} s)")


def _odd(n: int) -> int:
    while n % 2 == 0:
        n //= 2
    return n


if __name__ == "__main__":
    main()
