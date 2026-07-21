#!/usr/bin/env python3
"""A17 E20 — the overlay/lattice enumeration engine; P-33 driver.

Generalizes the (C)-table machinery (a17_c_matching_table.py) to the
mixed dA/dB overlay systems left open by E18/E19. Abstract variables
x = (a1..a4, b1..b4) over an arbitrary abelian group (a0 = b0 = 0 by
translation gauge); every branching choice adds an INTEGER relation
Σ c·x = 0; a branch dies only on SOUND kills:

  * a form f forced with denominator d = 1 (f ∈ Z-span of the
    relations ⟹ f = 0 in EVERY abelian realization):
      - vertex collision (a_p = a_q, b_p = b_q, translate collision),
      - Sidon repeat (two distinct ±edge classes equal, either side),
      - D2 violation (an a-edge equal to a b-edge),
      - a forbidden second membership (grid cells have m = 1),
      - grid-cell collision;
  * a DIFFERENCE-CLASS form (dA/dB element) forced with d = 2
    (2f ∈ Z-span ⟹ f is 2-torsion — dead by the D1 no-2-torsion
    lemma). d ≥ 3 is NOT a kill (the class frames have odd torsion);
    full Q-rank alone is NOT a kill either (rank-8 integer systems
    can retain d-torsion solutions) — death only via d ∈ {1,2} forms.

Soundness architecture: every kill is confirmed by exact rational
back-substitution against the Z-echelon lattice (unique combo ⟹
exact denominator); the mod-p machinery (p = 2^31 − 1) is only a
screen; every TERMINAL is re-verified against the full battery with
exact arithmetic. So the incremental bookkeeping can only affect
speed, never soundness.

P-33 geometry ((3,3) closure; Theorem E19.1 + doc §11). Grid cells
y(i,c) = a_{α(i,c)} + g_i = b_{β(i,c)} + r_c; gauge g_0 = 0 and
row-0 witnesses relabeled to α(0,c) = β(0,c) = c (S5×S5 freedom;
row-0 injectivity is proven). Eliminations:
    r_c = a_c − b_c,        g_i = b_{β(i,0)} − a_{α(i,0)},
cell relations for i, c ∈ {1,2}:
    a_{α(i,c)} + g_i − b_{β(i,c)} − r_c = 0,
plus 6 membership witnesses: Δ_{cc'} = r_c − r_{c'} ∈ dB and
Γ_{ii'} = g_i − g_{i'} ∈ dA, each an (s,t) witness branching. α, β
are row- AND column-injective (all four proven necessary). P-33
holds iff NO pattern survives.

Validation: --mode p33-inject rebuilds the shift-passing triangle
grids from live members (the E19 O4 population), maps them into the
abstract gauge, ASSERTS every relation concretely in the group, and
requires the engine NOT to kill them (they are realized, collided
configurations — checked against the static + translate battery,
since the full grid battery presumes an actual match). --mode
p33-relaxed drops injectivity (collided patterns must survive).

Usage:
  uv run python scripts/a17_e20_overlay_engine.py --mode selftest
  uv run python scripts/a17_e20_overlay_engine.py --mode p33-inject \
      --members data/a17/members_7x9.jsonl,data/a17/members_6x9_6x10.jsonl
  uv run python scripts/a17_e20_overlay_engine.py --mode p33 \
      --out data/a17/e20_p33_table.json
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from fractions import Fraction
from itertools import combinations
from math import gcd, lcm
from pathlib import Path

import numpy as np

NV = 8
P = (1 << 31) - 1


def avec(i: int) -> tuple:
    v = [0] * NV
    if i:
        v[i - 1] = 1
    return tuple(v)


def bvec(i: int) -> tuple:
    v = [0] * NV
    if i:
        v[4 + i - 1] = 1
    return tuple(v)


def vsub(a, b):
    return tuple(x - y for x, y in zip(a, b))


def vadd(a, b):
    return tuple(x + y for x, y in zip(a, b))


ZERO = tuple([0] * NV)


class ZLattice:
    """Integer echelon row lattice; exact membership denominators."""

    def __init__(self):
        self.rows: list[list[int]] = []   # echelon by pivot column

    def clone(self) -> "ZLattice":
        l2 = ZLattice()
        l2.rows = [r[:] for r in self.rows]
        return l2

    @staticmethod
    def _piv(r) -> int:
        return next(i for i, x in enumerate(r) if x)

    def add(self, vec) -> bool:
        """Insert integer vector; True iff the Z-span grew."""
        before = [r[:] for r in self.rows]
        v = list(vec)
        i = 0
        while i < len(self.rows):
            r = self.rows[i]
            c = self._piv(r)
            lead = next((j for j, x in enumerate(v) if x), None)
            if lead is None:
                break
            if lead < c:
                break
            if lead > c:
                i += 1
                continue
            while v[c]:
                q = v[c] // r[c]
                v = [x - q * y for x, y in zip(v, r)]
                if v[c]:
                    self.rows[i], v = v, r[:]
                    r = self.rows[i]
        lead = next((j for j, x in enumerate(v) if x), None)
        if lead is not None:
            if v[lead] < 0:
                v = [-x for x in v]
            self.rows.insert(
                next((k for k, r in enumerate(self.rows)
                      if self._piv(r) > lead), len(self.rows)), v)
        return self.rows != before

    def rank(self) -> int:
        return len(self.rows)

    def forced_denom(self, f) -> int | None:
        """Smallest d with d·f in the Z-span (unique echelon combo),
        or None if f is outside the Q-span. Fraction-free: the state
        (v, den) represents the true vector v/den; each elimination
        multiplies through by the pivot, and the combo coefficient
        for a row is v[c]/(den·r[c]), whose reduced denominator
        feeds the lcm."""
        v = list(f)
        den = 1
        dmax = 1
        for r in self.rows:
            c = self._piv(r)
            if v[c]:
                rc = r[c]
                dk = den * rc
                dk //= gcd(v[c], dk)
                dmax = lcm(dmax, abs(dk))
                vc = v[c]
                v = [rc * x - vc * y for x, y in zip(v, r)]
                den *= rc
                g = gcd(den, *(abs(x) for x in v)) if any(v) else den
                if g > 1:
                    den //= g
                    v = [x // g for x in v]
        if any(v):
            return None
        return dmax

    def basis_modq(self, q: int) -> list:
        """Echelon basis of the row space over F_q, as
        (normalized row, pivot) pairs. Used as a NECESSITY screen:
        d·f ∈ Z-span with gcd(d, q) = 1 forces f into the mod-q
        row space."""
        rows = []
        for r in self.rows:
            v = [x % q for x in r]
            for br, bp in rows:
                if v[bp]:
                    c = v[bp]
                    v = [(x - c * y) % q for x, y in zip(v, br)]
            lead = next((j for j, x in enumerate(v) if x), None)
            if lead is None:
                continue
            inv = pow(v[lead], -1, q)
            v = [(x * inv) % q for x in v]
            rows.append((v, lead))
        return rows

    def basis_modp(self) -> tuple:
        rows, pivs = [], []
        for r in self.rows:
            v = np.array(r, dtype=np.int64) % P
            for br, bp in zip(rows, pivs):
                v = (v - v[bp] * br) % P
            nz = np.nonzero(v)[0]
            if len(nz) == 0:
                continue
            p0 = int(nz[0])
            v = (v * pow(int(v[p0]), P - 2, P)) % P
            rows.append(v)
            pivs.append(p0)
        return rows, pivs


def killable(d: int | None, kind: str) -> bool:
    return d == 1 or (d == 2 and kind == "diff")


class Battery:
    """Kill-form battery with mod-p prefilter + exact confirm."""

    def __init__(self, forms: list, kinds: list):
        self.forms = forms
        self.kinds = kinds
        self.Fz = np.array(forms, dtype=np.int64)
        self.F = self.Fz % P
        self.Fq = {q: self.Fz % q for q in (2, 3, 5)}

    def reduced(self, lat: ZLattice) -> np.ndarray:
        rows, pivs = lat.basis_modp()
        R = self.F.copy()
        for br, bp in zip(rows, pivs):
            R = (R - np.outer(R[:, bp], br)) % P
        return R

    def first_kill(self, lat: ZLattice) -> int | None:
        R = self.reduced(lat)
        for idx in np.nonzero(~R.any(axis=1))[0]:
            if killable(lat.forced_denom(self.forms[idx]),
                        self.kinds[idx]):
                return int(idx)
        return None

    def extended(self, nf: list, nk: list) -> "Battery":
        """Append forms without rebuilding the numpy block."""
        b = Battery.__new__(Battery)
        b.forms = self.forms + nf
        b.kinds = self.kinds + nk
        if nf:
            add = np.array(nf, dtype=np.int64)
            b.Fz = np.vstack([self.Fz, add])
            b.F = np.vstack([self.F, add % P])
            b.Fq = {q: np.vstack([self.Fq[q], add % q])
                    for q in (2, 3, 5)}
        else:
            b.Fz = self.Fz
            b.F = self.F
            b.Fq = self.Fq
        return b

    def kill_candidates(self, lat: ZLattice) -> tuple:
        """(indices possibly killable, degenerate flag). A form can
        be killable only if its reduced row vanishes mod p AND mod
        q for every q ∈ {2,3,5} coprime to the allowed denominator
        ("eq": d = 1 needs all three; "diff": d ∈ {1,2} needs 3, 5)
        — necessary conditions, so skipping the rest is sound.
        degenerate=True (mod-p basis lost exact rank): caller must
        fall back to the exhaustive exact scan."""
        rows, pivs = lat.basis_modp()
        if len(rows) < lat.rank():
            return [], True
        R = self.F.copy()
        for br, bp in zip(rows, pivs):
            R = (R - np.outer(R[:, bp], br)) % P
        z = ~R.any(axis=1)
        nz = int(z.sum())
        if nz == 0:
            return [], False
        zidx = np.nonzero(z)[0]
        if nz <= 12:      # small: exact checks beat the screen
            return [int(i) for i in zidx], False
        oks = {}
        for q in (2, 3, 5):
            Rq = self.Fq[q][zidx]
            for br, bp in lat.basis_modq(q):
                brv = np.array(br, dtype=np.int64)
                Rq = (Rq - np.outer(Rq[:, bp], brv)) % q
            oks[q] = ~Rq.any(axis=1)
        out = []
        for k, idx in enumerate(zidx):
            idx = int(idx)
            if self.kinds[idx] == "diff":
                if oks[3][k] and oks[5][k]:
                    out.append(idx)
            elif oks[2][k] and oks[3][k] and oks[5][k]:
                out.append(idx)
        return out, False

    def exact_all(self, lat: ZLattice) -> list:
        """All killable form indices (exact; screened fast path
        with exhaustive fallback on mod-p degeneracy)."""
        cands, degen = self.kill_candidates(lat)
        if degen:
            return [i for i, (f, k) in enumerate(
                        zip(self.forms, self.kinds))
                    if killable(lat.forced_denom(f), k)]
        return [idx for idx in cands
                if killable(lat.forced_denom(self.forms[idx]),
                            self.kinds[idx])]


def static_battery() -> Battery:
    forms, kinds = [], []
    for p, q in combinations(range(5), 2):
        forms.append(vsub(avec(p), avec(q))); kinds.append("diff")
        forms.append(vsub(bvec(p), bvec(q))); kinds.append("diff")
    aedges = [vsub(avec(p), avec(q))
              for p, q in combinations(range(5), 2)]
    bedges = [vsub(bvec(p), bvec(q))
              for p, q in combinations(range(5), 2)]
    for edges in (aedges, bedges):
        for u, v in combinations(edges, 2):
            forms.append(vsub(u, v)); kinds.append("eq")
            forms.append(vadd(u, v)); kinds.append("eq")
    for u in aedges:
        for v in bedges:
            forms.append(vsub(u, v)); kinds.append("eq")
            forms.append(vadd(u, v)); kinds.append("eq")
    return Battery(forms, kinds)


def p33_geometry(alpha, beta) -> tuple:
    """r-, g-, y-vectors for a pattern (correct signs:
    r_c = a_c − b_c, g_i = b_{β(i,0)} − a_{α(i,0)})."""
    r = [ZERO, vsub(avec(1), bvec(1)), vsub(avec(2), bvec(2))]
    g = [ZERO,
         vsub(bvec(beta[1][0]), avec(alpha[1][0])),
         vsub(bvec(beta[2][0]), avec(alpha[2][0]))]
    y = [[vadd(avec(alpha[i][c]), g[i]) for c in range(3)]
         for i in range(3)]
    return r, g, y


def p33_translate_forms(r, g) -> tuple:
    forms, kinds = [], []
    for i, j in ((1, 0), (2, 0), (2, 1)):
        forms.append(vsub(g[i], g[j])); kinds.append("diff")
        forms.append(vsub(r[i], r[j])); kinds.append("diff")
    return forms, kinds


def p33_grid_forms(alpha, beta) -> tuple:
    """Translate + cell-distinctness + m = 1 membership forms."""
    r, g, y = p33_geometry(alpha, beta)
    forms, kinds = p33_translate_forms(r, g)
    cells = [y[i][c] for i in range(3) for c in range(3)]
    for u, v in combinations(cells, 2):
        forms.append(vsub(u, v)); kinds.append("eq")
    for i in range(3):
        for c in range(3):
            for i2 in range(3):
                if i2 != i:
                    for p in range(5):
                        forms.append(
                            vsub(vsub(y[i][c], g[i2]), avec(p)))
                        kinds.append("eq")
            for c2 in range(3):
                if c2 != c:
                    for p in range(5):
                        forms.append(
                            vsub(vsub(y[i][c], r[c2]), bvec(p)))
                        kinds.append("eq")
    return forms, kinds, r, g


AEDGE = [(p, q) for p in range(5) for q in range(5) if p != q]


def p33_search(relax: bool, max_nodes: int, stop_first: bool,
               progress: bool = True) -> dict:
    t0 = time.time()
    stat = static_battery()
    stats = {"nodes": 0, "dead": 0, "terminals": 0, "capped": False}
    terminals = []

    def bump() -> bool:
        stats["nodes"] += 1
        if progress and stats["nodes"] % 500000 == 0:
            print(f"  ..{stats['nodes']} nodes, dead={stats['dead']},"
                  f" term={stats['terminals']},"
                  f" {round(time.time() - t0)}s", flush=True)
        if stats["nodes"] >= max_nodes:
            stats["capped"] = True
        return not stats["capped"]

    def stop() -> bool:
        return stats["capped"] or (stop_first
                                   and stats["terminals"] > 0)

    def memb(k, lat, battery, RF, dset, steps, alpha, beta, wits):
        """Incremental witness layers. RF = battery reduced mod p by
        lat's Q-basis; dset = form indices in the Q-span, confirmed
        unkilled (their denominators can still drop)."""
        if k == len(steps):
            if battery.exact_all(lat):    # exhaustive exact confirm
                stats["dead"] += 1
                return
            stats["terminals"] += 1
            terminals.append({
                "alpha": [row[:] for row in alpha],
                "beta": [row[:] for row in beta],
                "wits": wits[:], "rank": lat.rank(),
                "lattice": [r[:] for r in lat.rows]})
            return
        target, side = steps[k]
        vecf = avec if side == "A" else bvec
        base_rows, base_pivs = lat.basis_modp()
        for (s, t) in AEDGE:
            if not bump():
                return
            rel = vsub(target, vsub(vecf(s), vecf(t)))
            l2 = lat.clone()
            grew = l2.add(rel)
            if not grew:                   # relation already known
                wits.append((side, s, t))
                memb(k + 1, l2, battery, RF, dset, steps,
                     alpha, beta, wits)
                wits.pop()
                if stop():
                    return
                continue
            dead = False
            for idx in dset:               # denominators can drop
                if killable(l2.forced_denom(battery.forms[idx]),
                            battery.kinds[idx]):
                    dead = True
                    break
            if dead:
                stats["dead"] += 1
                continue
            rr = np.array(rel, dtype=np.int64) % P
            for br, bp in zip(base_rows, base_pivs):
                rr = (rr - rr[bp] * br) % P
            if not rr.any():               # Q-span unchanged
                nstate = (RF, dset)
            else:
                j0 = int(np.nonzero(rr)[0][0])
                rr = (rr * pow(int(rr[j0]), P - 2, P)) % P
                RF2 = (RF - np.outer(RF[:, j0], rr)) % P
                newly = np.nonzero(RF.any(axis=1)
                                   & ~RF2.any(axis=1))[0]
                d2 = list(dset)
                for idx in newly:
                    d = l2.forced_denom(battery.forms[idx])
                    if killable(d, battery.kinds[idx]):
                        dead = True
                        break
                    if d is not None:
                        d2.append(int(idx))
                if dead:
                    stats["dead"] += 1
                    continue
                nstate = (RF2, d2)
            wits.append((side, s, t))
            memb(k + 1, l2, battery, nstate[0], nstate[1],
                 steps, alpha, beta, wits)
            wits.pop()
            if stop():
                return

    def cells_done(lat, alpha, beta):
        gforms, gkinds, r, g = p33_grid_forms(alpha, beta)
        battery = Battery(stat.forms + gforms, stat.kinds + gkinds)
        if battery.first_kill(lat) is not None:
            stats["dead"] += 1
            return
        RF = battery.reduced(lat)
        dset = [int(i) for i in np.nonzero(~RF.any(axis=1))[0]]
        steps = [(vsub(r[0], r[1]), "B"), (vsub(r[0], r[2]), "B"),
                 (vsub(r[1], r[2]), "B"), (g[1], "A"), (g[2], "A"),
                 (vsub(g[1], g[2]), "A")]
        memb(0, lat, battery, RF, dset, steps, alpha, beta, [])

    CELLS = [(1, 1), (1, 2), (2, 1), (2, 2)]

    def rec_cell(k, lat, alpha, beta, used_a, used_b):
        if k == len(CELLS):
            cells_done(lat, alpha, beta)
            return
        i, c = CELLS[k]
        acol = {alpha[i2][c] for i2 in range(3) if i2 != i
                and alpha[i2][c] is not None}
        arow = {alpha[i][c2] for c2 in range(3) if c2 != c
                and alpha[i][c2] is not None}
        bcol = {beta[i2][c] for i2 in range(3) if i2 != i
                and beta[i2][c] is not None}
        brow = {beta[i][c2] for c2 in range(3) if c2 != c
                and beta[i][c2] is not None}
        g_i = vsub(bvec(beta[i][0]), avec(alpha[i][0]))
        r_c = vsub(avec(c), bvec(c))
        for av in range(5):
            if not relax and (av in acol or av in arow):
                continue
            if av == 4 and 3 not in used_a:
                continue
            for bv in range(5):
                if not relax and (bv in bcol or bv in brow):
                    continue
                if bv == 4 and 3 not in used_b:
                    continue
                if not bump():
                    return
                rel = vsub(vsub(vadd(avec(av), g_i), bvec(bv)), r_c)
                l2 = lat.clone()
                l2.add(rel)
                if stat.first_kill(l2) is not None:
                    stats["dead"] += 1
                    continue
                alpha[i][c], beta[i][c] = av, bv
                rec_cell(k + 1, l2, alpha, beta,
                         used_a | {av}, used_b | {bv})
                alpha[i][c] = beta[i][c] = None
                if stop():
                    return

    for a10 in range(1, 4):           # label 4 cannot precede 3
        for b10 in range(1, 4):
            for a20 in range(1, 5):
                if not relax and a20 == a10:
                    continue
                if a20 == 4 and a10 != 3:
                    continue
                for b20 in range(1, 5):
                    if not relax and b20 == b10:
                        continue
                    if b20 == 4 and b10 != 3:
                        continue
                    if (a10, b10) > (a20, b20):
                        continue      # row-swap canonicalization
                    alpha = [[0, 1, 2], [a10, None, None],
                             [a20, None, None]]
                    beta = [[0, 1, 2], [b10, None, None],
                            [b20, None, None]]
                    rec_cell(0, ZLattice(), alpha, beta,
                             {0, 1, 2, a10, a20},
                             {0, 1, 2, b10, b20})
                    if stop():
                        stats["secs"] = round(time.time() - t0, 1)
                        stats["terminal_patterns"] = terminals
                        return stats
    stats["secs"] = round(time.time() - t0, 1)
    stats["terminal_patterns"] = terminals
    return stats


def p33_inject(member_paths: list[str]) -> dict:
    """Validation: realized shift-passing grids must NOT be killed."""
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent
                           / "src"))
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    from a17_e18_k4_census import diffs
    from a17_e19_odd_census import triangles
    from bb_lab.group import AbelianGroup
    from bb_lab.poly import Poly

    def gmul(G, k, x):
        acc = tuple(0 for _ in range(G.rank))
        step = x if k > 0 else G.neg(x)
        for _ in range(abs(k)):
            acc = G.add(acc, step)
        return acc

    stat = static_battery()
    out = {"grids": 0, "gauged": 0, "validated": 0, "killed": [],
           "collided": 0}
    rows = []
    for path in member_paths:
        with open(path) as f:
            rows += [json.loads(l) for l in f if '"A"' in l]
    for rmem in rows:
        G = AbelianGroup(tuple(rmem["frame"]))
        A = sorted(Poly.from_string(rmem["A"], G).support)
        B = sorted(Poly.from_string(rmem["B"], G).support)
        dAs, dBs = diffs(frozenset(A), G), diffs(frozenset(B), G)
        S_map = {}
        for a in A:
            for b in B:
                S_map[G.sub(a, b)] = (a, b)  # A−B keys (r−g = a−b)
        elems = list(G)
        for Tg in triangles(dAs, G):
            for Tr0 in triangles(dBs, G):
              for tau in elems:
                Tr = tuple(G.add(t, tau) for t in Tr0)
                wit = {}
                ok = True
                for i in range(3):
                    for c in range(3):
                        w = S_map.get(G.sub(Tr[c], Tg[i]))
                        if w is None:
                            ok = False
                            break
                        wit[(i, c)] = w
                    if not ok:
                        break
                if not ok:
                    continue
                out["grids"] += 1
                if out["gauged"] >= 2000:
                    continue          # validation cap
                row0 = next(
                    (i for i in range(3)
                     if len({wit[(i, c)][0] for c in range(3)}) == 3
                     and len({wit[(i, c)][1]
                              for c in range(3)}) == 3), None)
                if row0 is None:
                    continue       # no gauge-able row; skip
                out["gauged"] += 1
                order = [row0] + [i for i in range(3) if i != row0]
                aL = {wit[(row0, c)][0]: c for c in range(3)}
                for a in A:
                    if a not in aL:
                        aL[a] = len(aL)
                bL = {wit[(row0, c)][1]: c for c in range(3)}
                for b in B:
                    if b not in bL:
                        bL[b] = len(bL)
                alpha = [[aL[wit[(order[i], c)][0]]
                          for c in range(3)] for i in range(3)]
                beta = [[bL[wit[(order[i], c)][1]]
                         for c in range(3)] for i in range(3)]
                if any(len(set(alpha[i])) < 3 for i in range(3)) \
                        or any(len({beta[i][c] for i in range(3)}) < 3
                               for c in range(3)):
                    out["collided"] += 1
                rels = []
                r, g, _ = p33_geometry(alpha, beta)
                for i in (1, 2):
                    for c in (1, 2):
                        rels.append(vsub(vsub(
                            vadd(avec(alpha[i][c]), g[i]),
                            bvec(beta[i][c])), r[c]))
                for c1, c2 in ((0, 1), (0, 2), (1, 2)):
                    delta = G.sub(Tr[c1], Tr[c2])
                    s, t = next((bL[b1], bL[b2]) for b1 in B
                                for b2 in B
                                if G.sub(b1, b2) == delta)
                    rels.append(vsub(vsub(r[c1], r[c2]),
                                     vsub(bvec(s), bvec(t))))
                for i1, i2 in ((1, 0), (2, 0), (2, 1)):
                    gam = G.sub(Tg[order[i1]], Tg[order[i2]])
                    p_, q_ = next((aL[a1], aL[a2]) for a1 in A
                                  for a2 in A
                                  if G.sub(a1, a2) == gam)
                    rels.append(vsub(vsub(g[i1], g[i2]),
                                     vsub(avec(p_), avec(q_))))
                # concrete verification of every abstract relation
                aInv = {v: k for k, v in aL.items()}
                bInv = {v: k for k, v in bL.items()}
                xs = [G.sub(aInv[l], aInv[0]) for l in (1, 2, 3, 4)]
                xs += [G.sub(bInv[l], bInv[0]) for l in (1, 2, 3, 4)]
                for rel in rels:
                    acc = tuple(0 for _ in range(G.rank))
                    for cf, xv in zip(rel, xs):
                        if cf:
                            acc = G.add(acc, gmul(G, cf, xv))
                    assert all(v == 0 for v in acc), \
                        f"relation fails concretely: {rel}"
                lat = ZLattice()
                for rel in rels:
                    lat.add(rel)
                tforms, tkinds = p33_translate_forms(r, g)
                bat = Battery(stat.forms + tforms,
                              stat.kinds + tkinds)
                hit = bat.first_kill(lat)
                if hit is None:
                    out["validated"] += 1
                else:
                    out["killed"].append({
                        "frame": rmem["frame"], "form": int(hit)})
    return out


def selftest() -> None:
    lat = ZLattice()
    lat.add(tuple([2] + [0] * 7))
    assert lat.forced_denom(avec(1)) == 2
    assert lat.forced_denom(avec(2)) is None
    lat2 = ZLattice()
    lat2.add(vsub(avec(1), avec(2)))
    assert lat2.forced_denom(vsub(avec(1), avec(2))) == 1
    assert lat2.forced_denom(avec(1)) is None
    lat3 = ZLattice()
    lat3.add(vsub(avec(1), bvec(1)))
    lat3.add(vsub(avec(1), bvec(2)))
    assert lat3.forced_denom(vsub(bvec(1), bvec(2))) == 1
    b = static_battery()
    assert b.first_kill(lat3) is not None, "b1=b2 must kill"
    lat4 = ZLattice()
    lat4.add(vsub(avec(1), bvec(1)))
    assert b.first_kill(lat4) is not None, "a1=b1 is a D2 kill"
    lat5 = ZLattice()
    lat5.add(vsub(vsub(avec(1), avec(2)), vsub(bvec(1), bvec(2))))
    assert b.first_kill(lat5) is not None, "D2 edge equality kill"
    lat6 = ZLattice()
    lat6.add(vadd(vsub(avec(1), bvec(1)), vsub(avec(2), bvec(2))))
    assert b.first_kill(lat6) is None, "allowed mixed relation lives"
    lat7 = ZLattice()
    lat7.add(tuple([2] + [0] * 7))
    assert lat7.add(avec(1)) is True, "Z-refinement must count"
    assert lat7.forced_denom(avec(1)) == 1
    lat8 = ZLattice()
    lat8.add(avec(1))
    assert lat8.add(avec(1)) is False, "re-adding must not grow"
    print("selftest OK")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--mode", required=True,
                    choices=["selftest", "p33", "p33-relaxed",
                             "p33-inject"])
    ap.add_argument("--max-nodes", type=int, default=200_000_000)
    ap.add_argument("--members", default=None)
    ap.add_argument("--out", default=None)
    args = ap.parse_args()
    if args.mode == "selftest":
        selftest()
        return
    if args.mode == "p33-inject":
        res = p33_inject(args.members.split(","))
        print(json.dumps(res, indent=1))
        ok = (res["gauged"] > 0 and res["killed"] == []
              and res["validated"] == res["gauged"])
        print(f"\nINJECTION {'OK' if ok else 'FAILED'}: "
              f"{res['validated']}/{res['gauged']} gauged grids "
              f"survive the engine ({res['collided']} collided, "
              f"{res['grids']} total)")
        sys.exit(0 if ok else 1)
    relax = args.mode == "p33-relaxed"
    res = p33_search(relax=relax, max_nodes=args.max_nodes,
                     stop_first=relax)
    summary = {k: v for k, v in res.items()
               if k != "terminal_patterns"}
    summary["n_terminals"] = len(res["terminal_patterns"])
    print(json.dumps(summary, indent=1))
    if args.out:
        Path(args.out).parent.mkdir(parents=True, exist_ok=True)
        Path(args.out).write_text(json.dumps(res, indent=1))
    if relax:
        ok = res["terminals"] > 0
        print(f"\nVALIDATION {'OK' if ok else 'FAILED'}: relaxed run"
              f" {'found' if ok else 'found NO'} terminals")
        sys.exit(0 if ok else 1)
    ok = res["terminals"] == 0 and not res["capped"]
    verdict = ("PROVEN — no injective pattern survives" if ok
               else "HAS RESIDUALS/CAPPED")
    print(f"\nP-33 {verdict}: terminals={res['terminals']}, "
          f"nodes={res['nodes']}, capped={res['capped']}")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
