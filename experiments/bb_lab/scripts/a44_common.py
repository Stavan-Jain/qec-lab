#!/usr/bin/env python3
"""A44 shared helpers — Laurent-pair view of member logicals, shear
frames and their pullbacks, and a generic SAT witness hunter.

Conventions (A40 §1, `a40_s11_compare.py`): fixed pair
A = 1 + y + x^3 y^-1, B = 1 + x + x^-1 y^-3 on Z_l x Z_m, H_X = [M_A|M_B],
H_Z = [M_B^T|M_A^T].  With alpha := A(x^-1,y^-1) = 1 + y^-1 + x^-3 y and
beta := B(x^-1,y^-1) = 1 + x^-1 + x y^3, a vector v = (f, g) (block 0 = f,
block 1 = g as group-ring elements) is an X-cycle iff

        beta * f = alpha * g        (in F2[x,y]/(x^l-1, y^m-1)),

and an X-stabilizer iff (f, g) = h * (alpha, beta) for some h.  (Both are
verified against the code's own matrices in `assert_laurent_conventions`.)

Only this session's scripts import this module.  No SAT on any
certificate path: `hunt` returns kernel-verified WITNESSES (upper bounds)
and its UNSAT / timeout statuses are observations.
"""
from __future__ import annotations

import os
import resource
import sys
import time
from pathlib import Path

import numpy as np

LAB = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(LAB / "src"))
sys.path.insert(0, str(LAB / "scripts"))

from bb_lab.tower import TowerCode, v2i, i2v  # noqa: E402
import a40_s11_compare as C  # noqa: E402
from a40_s4_phase_triage import snf2, quotient_code  # noqa: E402

DATA = LAB / "data" / "a44"
DATA.mkdir(parents=True, exist_ok=True)
RSS_CAP_GB = 2.5

member_code = C.member_code


def rss_gb() -> float:
    ru = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return ru / (1 << 30)  # darwin reports bytes


def rss_guard(tag=""):
    r = rss_gb()
    if r > RSS_CAP_GB:
        raise RuntimeError(f"RSS {r:.2f} GB > cap {RSS_CAP_GB} at {tag}")
    return r


class Logger:
    def __init__(self, path: Path, append=False):
        self.path = path
        if not append:
            path.write_text(f"=== {time.strftime('%Y-%m-%d %H:%M:%S')} {path.name}\n")

    def __call__(self, msg: str):
        print(msg, flush=True)
        with self.path.open("a") as f:
            f.write(msg + "\n")


# ------------------------------------------------------------ Laurent view
Poly = frozenset  # of (ex, ey) exponent pairs, i.e. an F2 group-ring element


def pmul(p: Poly, q: Poly, lm) -> Poly:
    l, m = lm
    cnt: dict = {}
    for a in p:
        for b in q:
            k = ((a[0] + b[0]) % l, (a[1] + b[1]) % m)
            cnt[k] = cnt.get(k, 0) ^ 1
    return frozenset(k for k, v in cnt.items() if v)


def padd(p: Poly, q: Poly) -> Poly:
    return frozenset(p ^ q)


def pshift(p: Poly, s, lm) -> Poly:
    l, m = lm
    return frozenset(((a[0] + s[0]) % l, (a[1] + s[1]) % m) for a in p)


def alpha_beta(lm):
    """(alpha, beta) = (A-bar, B-bar) reduced mod the frame."""
    l, m = lm
    alpha = frozenset({(0, 0), (0, (-1) % m), ((-3) % l, 1)})
    beta = frozenset({(0, 0), ((-1) % l, 0), (1, 3 % m)})
    return alpha, beta


def pair_of(code: TowerCode, v: np.ndarray):
    """(f, g) as Polys from a vector."""
    f, g = set(), set()
    for blk, x, y in C.cells_of(code, v):
        (f if blk == 0 else g).add((x, y))
    return frozenset(f), frozenset(g)


def vec_of_pair(code: TowerCode, f: Poly, g: Poly) -> np.ndarray:
    return C.vec_of(code, [(0,) + e for e in f] + [(1,) + e for e in g])


def is_cycle_pair(f: Poly, g: Poly, lm) -> bool:
    alpha, beta = alpha_beta(lm)
    return pmul(beta, f, lm) == pmul(alpha, g, lm)


def assert_laurent_conventions(code: TowerCode):
    """beta f = alpha g <=> H_Z v = 0 and stabilizer generator = (alpha,
    beta) shifted, checked on random vectors and on the H_X rows."""
    C.assert_conventions(code)
    lm = tuple(code.G.orders)
    alpha, beta = alpha_beta(lm)
    rng = np.random.default_rng(44)
    for _ in range(20):
        v = (rng.integers(2, size=code.n)).astype(np.uint8)
        f, g = pair_of(code, v)
        assert is_cycle_pair(f, g, lm) == code.is_cycle(v)
        c = code.random_cycle(rng)
        f, g = pair_of(code, c)
        assert is_cycle_pair(f, g, lm)
    # stabilizer at g0 = (alpha, beta) * x^g0
    for i in [0, 1, code.ng // 2, code.ng - 1]:
        f, g = pair_of(code, code.HX[i])
        g0 = code.G.from_index(i)
        assert f == pshift(alpha, g0, lm) and g == pshift(beta, g0, lm)


def poly_str(p: Poly, lm=None, centre=None) -> str:
    """Monomial list x^a y^b, sorted, optionally re-centred by
    subtracting `centre` with symmetric representatives."""
    def rep(e):
        if lm is None or centre is None:
            return e
        l, m = lm
        a = (e[0] - centre[0]) % l
        b = (e[1] - centre[1]) % m
        if a > l // 2:
            a -= l
        if b > m // 2:
            b -= m
        return (a, b)
    terms = sorted(rep(e) for e in p)

    def mono(e):
        a, b = e
        s = []
        if a:
            s.append(f"x^{a}" if a != 1 else "x")
        if b:
            s.append(f"y^{b}" if b != 1 else "y")
        return "*".join(s) if s else "1"
    return " + ".join(mono(e) for e in terms) if terms else "0"


# ------------------------------------------------------------ shear frames
class Shear:
    """Z^2 / <(l,0),(d,p)> as a TowerCode (SNF coordinates) with the
    projection phi from Z^2 and the pullback to a covering torus."""

    def __init__(self, l: int, p: int, d: int):
        self.l, self.p, self.d = l, p, d
        self.code, self.orders = quotient_code(l, p, d)
        M = [[l, 0], [d, p]]
        D, U, V = snf2(M)
        self.V = V
        assert (abs(D[0][0]), abs(D[1][1])) == tuple(self.orders)

    def phi(self, e):
        V, (o1, o2) = self.V, self.orders
        return ((e[0] * V[0][0] + e[1] * V[1][0]) % o1,
                (e[0] * V[0][1] + e[1] * V[1][1]) % o2)

    def covers(self, target: TowerCode) -> bool:
        L, M = target.G.orders
        return self.phi((L, 0)) == (0, 0) and self.phi((0, M)) == (0, 0)

    def index_in(self, target: TowerCode) -> int:
        L, M = target.G.orders
        return (L * M) // (self.l * self.p)

    def pullback(self, v: np.ndarray, target: TowerCode) -> np.ndarray:
        assert self.covers(target)
        L, M = target.G.orders
        vt = np.zeros(target.n, dtype=np.uint8)
        code = self.code
        for blk in range(2):
            for x in range(L):
                for y in range(M):
                    if v[blk * code.ng + code.G.index(self.phi((x, y)))]:
                        vt[blk * target.ng + target.G.index((x, y))] = 1
        assert target.is_cycle(vt)
        return vt

    def pullback_pair(self, v: np.ndarray, target: TowerCode):
        return pair_of(target, self.pullback(v, target))


# ------------------------------------------------------------ SAT hunter
def hunt(code: TowerCode, Wmax: int, nsol=10, time_limit=300.0, tag="",
         dense=False, forbid_stab_only=True, extra_clauses=None,
         min_weight=None):
    """Existence lane: nontrivial X-logicals of weight <= Wmax (optionally
    >= min_weight, optionally gap-dense on both axes).  Witnesses are
    re-verified in the kernel; whole translation orbits are blocked after
    each find.  Status 'unsat'/'timeout' are OBSERVATIONS only."""
    from pycryptosat import Solver
    from pysat.card import CardEnc, EncType
    l, m = code.G.orders
    n = code.n
    top = n
    clauses = list(extra_clauses or [])
    xors = []
    for row in code.HZ:
        xors.append(([int(i) + 1 for i in np.nonzero(row)[0]], False))
    pv = []
    for zr in code.zreps:
        top += 1
        pv.append(top)
        xors.append(([int(i) + 1 for i in np.nonzero(zr)[0]] + [top], False))
    clauses.append(pv)
    if dense:
        for axis in range(2):
            order = code.G.orders[axis]
            other = code.G.orders[1 - axis]
            for c0 in range(order):
                lits = []
                for blk in range(2):
                    for dc in range(C.SPAN):
                        for q in range(other):
                            e = ((c0 + dc) % order, q) if axis == 0 else \
                                (q, (c0 + dc) % order)
                            lits.append(blk * code.ng + code.G.index(e) + 1)
                clauses.append(lits)
    card = CardEnc.atmost(lits=list(range(1, n + 1)), bound=Wmax,
                          top_id=top, encoding=EncType.seqcounter)
    clauses.extend(card.clauses)
    top = max(top, card.nv)
    if min_weight is not None and min_weight > 0:
        card2 = CardEnc.atleast(lits=list(range(1, n + 1)), bound=min_weight,
                                top_id=top, encoding=EncType.seqcounter)
        clauses.extend(card2.clauses)
        top = max(top, card2.nv)
    sols, status = [], "unknown"
    t0 = time.time()
    while len(sols) < nsol:
        left = time_limit - (time.time() - t0)
        if left <= 0:
            status = "timeout"
            break
        s = Solver(time_limit=max(1.0, left))
        for cl in clauses:
            s.add_clause(cl)
        for lits, rhs in xors:
            s.add_xor_clause(lits, rhs)
        sat, sol = s.solve()
        if sat is None:
            status = "timeout"
            break
        if not sat:
            status = "unsat" if not sols else "exhausted"
            break
        v = np.array([1 if sol[i + 1] else 0 for i in range(n)], dtype=np.uint8)
        assert code.is_cycle(v) and not code.is_stab(v) and v.sum() <= Wmax
        cells = C.cells_of(code, v)
        sols.append(dict(w=int(v.sum()), cells=sorted(cells),
                         sig=int(v2i(code.sig(v))), t=round(time.time() - t0, 1)))
        print(f"  [{tag}] witness {len(sols)}: w={int(v.sum())} "
              f"({time.time() - t0:.1f} s)", flush=True)
        for dx in range(l):
            for dy in range(m):
                sh = [blk * code.ng + code.G.index(((x + dx) % l, (y + dy) % m)) + 1
                      for blk, x, y in cells]
                clauses.append([-lit for lit in sh])
        status = "found"
    return sols, status, round(time.time() - t0, 1)


def picture(code: TowerCode, cells) -> str:
    l, m = code.G.orders
    grid = [["." for _ in range(l)] for _ in range(m)]
    for blk, x, y in cells:
        ch = "a" if blk == 0 else "b"
        cur = grid[m - 1 - y][x]
        grid[m - 1 - y][x] = ch if cur == "." else "*"
    rows = ["".join(r) for r in grid]
    return "\n".join(f"y={m-1-i:2d} {r}" for i, r in enumerate(rows)) + \
        "\n     " + "".join(str(x % 10) for x in range(l))


def translation_canon(code: TowerCode, v: np.ndarray) -> int:
    """Canonical int over the translation orbit (min over shifts)."""
    l, m = code.G.orders
    cells = C.cells_of(code, v)
    best = None
    for dx in range(l):
        for dy in range(m):
            k = v2i(C.vec_of(code, [(b, (x + dx) % l, (y + dy) % m)
                                    for b, x, y in cells]))
            if best is None or k < best:
                best = k
    return best
