"""The tower slice calculus as a library (A38 front F5, session S1).

The A32/A33/A36 distance-descent machinery and the A35 generality screen,
promoted from frozen session scripts into one reusable, **rank-generic**
module: groups are ``AbelianGroup`` order tuples of any rank (the group
layer was already dimension-generic), decks are "one axis, index-2 fold"
objects, and every screen/costing entry point takes order tuples — never
scalar ``(ell, m)``.

Ancestry (the session scripts stay frozen as artifacts; this module is the
library-grade copy, following the ``a30_coset_bz.py -> bb_lab.cosetbz``
precedent):

* ``scripts/a32_tower_slice.py``   — BBCode/Deck frames, slice + carry
  machinery, H1 maps, translation action (A32 Lemmas 1-2, Thm 3).
* ``scripts/a32_subclosures.py``   — the cap<=4 restricted-MITM lift-fiber
  enumerator (complete by the exact-off-support subset-sum argument).
* ``scripts/a33_tower_cells.py``   — h1_map / translation_mat / rep_for.
* ``scripts/a32_sectorAC_full.py`` — batched translation-canonical keys.
* ``scripts/a35_generality_screen.py`` (branch
  ``claude/tower-slice-calculus-generalize-410ed1``) — the precondition
  screen: per-rung structure, per-pair regime lattice, cost gates G1-G5.

Falsify-first contract (charter A38 §6.0): every engine change re-passes
:func:`validate_banked` — the banked A32/A33 structural asserts, the
census node-count anchors (x1.00 / x1.01, and the A33 x3.00 shared-walk
accounting datum), the A19 deck-survey k-verdicts, and the bit-level
sector-C fiber-layer reproduction (397 fibers / 4,132 lifts / exact
m2-histogram) — before anything new runs.

Claim tiers: nothing in this module makes a distance claim.  Screen
verdicts GREEN/AMBER/RED are **cost** verdicts on quantified gates;
censuses/fibers/rungs built on top of it inherit the deterministic
certificate tier of A30/A32/A33/A36 (exact node-count invariants,
complete-by-construction enumerations, in-line re-verification), never
more.  SAT appears nowhere.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from itertools import combinations
from typing import Callable, Iterable, Optional, Sequence

import numpy as np

from .checks import circulant
from .group import AbelianGroup
from .poly import Poly

__all__ = [
    "v2i", "i2v", "rref_ints", "reduce_int", "in_span", "gf2_rank",
    "kernel_basis", "rref_np", "gf2_inv", "colspace", "kernel_ints",
    "apply_map", "span_eq", "span_leq", "span_points", "preimage_basis",
    "TowerCode", "AxisDeck", "fold_support", "support_twisted",
    "support_str", "h1_map", "translation_mat", "rep_for",
    "translation_action", "translation_perms", "perm_for", "canon_key",
    "batch_keys", "orbits", "census_nodes", "enumerate_lifts",
    "liftable_codim", "fiber_sample", "screen_rung", "screen_tower",
    "gate_verdict", "tower_inventory", "A35_DOCKET", "A35_NO_DECK",
    "validate_banked",
]


# --------------------------------------------------------------- GF(2) core
# Bit-int rows (LSB = coordinate 0) for RREF-style algebra; numpy uint8
# vectors for chain-level arithmetic.  Verbatim algorithms from
# a30_rung_pass.py / a32_tower_slice.py so that every downstream number is
# reproduced exactly.

def v2i(v: np.ndarray) -> int:
    """Pack a 0/1 vector into a bit-int (LSB = index 0)."""
    x = 0
    for i in np.nonzero(v)[0]:
        x |= 1 << int(i)
    return x


def i2v(x: int, n: int) -> np.ndarray:
    """Unpack a bit-int into a length-n uint8 vector."""
    return np.array([(x >> i) & 1 for i in range(n)], dtype=np.uint8)


def rref_ints(rows: list[int]) -> tuple[list[int], list[int]]:
    """RREF over GF(2) on bit-int rows; returns (basis rows, pivot bits)."""
    basis: list[int] = []
    piv: list[int] = []
    for r in rows:
        cur = r
        for b, p in zip(basis, piv):
            if (cur >> p) & 1:
                cur ^= b
        if cur:
            p = (cur & -cur).bit_length() - 1
            for i in range(len(basis)):
                if (basis[i] >> p) & 1:
                    basis[i] ^= cur
            basis.append(cur)
            piv.append(p)
    return basis, piv


def reduce_int(x: int, basis: list[int], piv: list[int]) -> int:
    for b, p in zip(basis, piv):
        if (x >> p) & 1:
            x ^= b
    return x


def in_span(x: int, basis: list[int], piv: list[int]) -> bool:
    return reduce_int(x, basis, piv) == 0


def gf2_rank(vectors: Iterable[int]) -> int:
    b, _ = rref_ints([v for v in vectors])
    return len(b)


def kernel_basis(M: np.ndarray) -> list[np.ndarray]:
    """Basis of ker M over GF(2), M uint8 (rows x cols)."""
    M = M.copy() % 2
    rows, n = M.shape
    piv, r = [], 0
    for c in range(n):
        if r >= rows:
            break
        nz = np.nonzero(M[r:, c])[0]
        if len(nz) == 0:
            continue
        M[[r, r + nz[0]]] = M[[r + nz[0], r]]
        for i in np.nonzero(M[:, c])[0]:
            if i != r:
                M[i] ^= M[r]
        piv.append(c)
        r += 1
    Mr = M[:r]
    pivset = set(piv)
    out = []
    for c in [c for c in range(n) if c not in pivset]:
        v = np.zeros(n, dtype=np.uint8)
        v[c] = 1
        for i, pc in enumerate(piv):
            if Mr[i, c]:
                v[pc] ^= 1
        out.append(v)
    return out


def rref_np(M: np.ndarray) -> tuple[np.ndarray, list[int]]:
    """RREF of a uint8 matrix; returns (reduced rows, pivot columns)."""
    M = M.copy() % 2
    rows, n = M.shape
    piv, r = [], 0
    for c in range(n):
        if r >= rows:
            break
        nz = np.nonzero(M[r:, c])[0]
        if len(nz) == 0:
            continue
        M[[r, r + nz[0]]] = M[[r + nz[0], r]]
        for i in np.nonzero(M[:, c])[0]:
            if i != r:
                M[i] ^= M[r]
        piv.append(c)
        r += 1
    return M[:r], piv


def gf2_inv(M: np.ndarray) -> np.ndarray:
    k = M.shape[0]
    aug = np.concatenate([M % 2, np.eye(k, dtype=np.uint8)], axis=1)
    R, piv = rref_np(aug)
    assert piv == list(range(k)), "matrix not invertible"
    return R[:, k:]


def colspace(M: np.ndarray) -> list[int]:
    b, _ = rref_ints([v2i(c) for c in (M % 2).T])
    return b


def kernel_ints(M: np.ndarray) -> list[int]:
    return [v2i(v) for v in kernel_basis(M % 2)]


def apply_map(M: np.ndarray, sig_int: int) -> int:
    v = i2v(sig_int, M.shape[1])
    return v2i((M @ v) % 2)


def span_eq(a: list[int], b: list[int]) -> bool:
    ab, _ = rref_ints(list(a))
    bb, _ = rref_ints(list(b))
    return sorted(ab) == sorted(bb)


def span_leq(a: list[int], b: list[int]) -> bool:
    bb, bp = rref_ints(list(b))
    return all(in_span(x, bb, bp) for x in a)


def span_points(basis: Iterable[int]) -> set[int]:
    pts = {0}
    for b in basis:
        pts |= {p ^ b for p in pts}
    return pts


def preimage_basis(M: np.ndarray, Wb: list[int], Wp: list[int]) -> list[int]:
    """Basis of {s : M s in span(Wb)} = ker(F o M), F = the functionals
    vanishing on span(Wb).  Full-codomain W means empty annihilator and
    the preimage is the whole domain (the a158028 fix, carried over)."""
    n = M.shape[0]
    Wmat = np.array([i2v(x, n) for x in Wb], dtype=np.uint8) \
        if Wb else np.zeros((0, n), dtype=np.uint8)
    F = np.array(kernel_basis(Wmat), dtype=np.uint8)
    if F.size == 0:
        return [v2i(v) for v in np.eye(M.shape[1], dtype=np.uint8)]
    QM = (F @ M) % 2
    return [v2i(v) for v in kernel_basis(QM)]


# ------------------------------------------------------- polynomial descent
def _as_support(poly, G: AbelianGroup) -> frozenset[tuple[int, ...]]:
    """Accept a poly string, a Poly, or an iterable of exponent tuples."""
    if isinstance(poly, Poly):
        return frozenset(G.reduce(g) for g in poly.support)
    if isinstance(poly, str):
        return frozenset(Poly.from_string(poly, G).support)
    return frozenset(G.reduce(tuple(int(e) for e in g)) for g in poly)


def fold_support(support: Iterable[tuple[int, ...]], axis: int,
                 newmod: int) -> frozenset[tuple[int, ...]]:
    """Descend a polynomial support along an index-2 axis fold
    (mod-2 term merging over F2)."""
    counts: dict[tuple[int, ...], int] = {}
    for t in support:
        u = list(t)
        u[axis] %= newmod
        key = tuple(u)
        counts[key] = counts.get(key, 0) + 1
    return frozenset(t for t, c in counts.items() if c % 2 == 1)


def support_twisted(support: Iterable[tuple[int, ...]], axis: int,
                    newmod: int) -> bool:
    """Does any cover term carry the deck element (exponent >= newmod)?"""
    return any(t[axis] >= newmod for t in support)


_VAR_NAMES = "xyzwvu"


def support_str(support: Iterable[tuple[int, ...]]) -> str:
    """Canonical human-readable form, any rank ('0' for the zero poly)."""
    def mono(t: tuple[int, ...]) -> str:
        fs = []
        for a, e in enumerate(t):
            if e:
                fs.append(f"{_VAR_NAMES[a]}^{e}" if e > 1 else _VAR_NAMES[a])
        return "*".join(fs) if fs else "1"
    sup = sorted(support)
    return " + ".join(mono(t) for t in sup) if sup else "0"


# ------------------------------------------------------------------- codes
class TowerCode:
    """One two-block group-algebra CSS code in the lab convention
    (H_X = [M_A|M_B], H_Z = [M_B^T|M_A^T]), over an AbelianGroup of ANY
    rank, with the homology helpers the calculus consumes.

    Rank-generic port of ``a32_tower_slice.BBCode`` (which pinned
    ``lm: tuple[int, int]``); for rank 2 every field is numerically
    identical (same row-major group enumeration, same circulant, same
    RREF order)."""

    def __init__(self, name: str, orders: Sequence[int], A, B):
        self.name = name
        self.G = AbelianGroup(tuple(int(o) for o in orders))
        self.A = Poly.from_support(_as_support(A, self.G), self.G)
        self.B = Poly.from_support(_as_support(B, self.G), self.G)
        ng = self.G.cardinality
        self.ng, self.n = ng, 2 * ng
        MA = circulant(self.A).astype(np.uint8) % 2
        MB = circulant(self.B).astype(np.uint8) % 2
        self.HX = np.concatenate([MA, MB], axis=1) % 2          # ng x n
        self.HZ = np.concatenate([MB.T, MA.T], axis=1) % 2      # ng x n
        assert not ((self.HX @ self.HZ.T) % 2).any(), "CSS fails"
        self.rsHX_b, self.rsHX_p = rref_ints([v2i(r) for r in self.HX])
        self.rsHZ_b, self.rsHZ_p = rref_ints([v2i(r) for r in self.HZ])
        self.kerHZ = kernel_basis(self.HZ)   # X-side cycles
        self.kerHX = kernel_basis(self.HX)   # Z-side cycles
        self.k = len(self.kerHZ) - len(self.rsHX_b)
        # Z-logical reps: ker HX mod rowspace HZ
        bb, pp = list(self.rsHZ_b), list(self.rsHZ_p)
        zreps = []
        for kv in self.kerHX:
            x = reduce_int(v2i(kv), bb, pp)
            if x:
                bb.append(x)
                pp.append((x & -x).bit_length() - 1)
                zreps.append(kv)
        assert len(zreps) == self.k
        self.zreps = np.array(zreps, dtype=np.uint8) if self.k else \
            np.zeros((0, self.n), dtype=np.uint8)
        # X-logical reps: ker HZ mod rowspace HX
        bb, pp = list(self.rsHX_b), list(self.rsHX_p)
        xreps = []
        for kv in self.kerHZ:
            x = reduce_int(v2i(kv), bb, pp)
            if x:
                bb.append(x)
                pp.append((x & -x).bit_length() - 1)
                xreps.append(kv)
        assert len(xreps) == self.k
        self.xreps = np.array(xreps, dtype=np.uint8) if self.k else \
            np.zeros((0, self.n), dtype=np.uint8)
        if self.k:
            P = (self.xreps @ self.zreps.T) % 2
            assert gf2_rank([v2i(r) for r in P]) == self.k, \
                "degenerate pairing"

    # -- membership / signatures -------------------------------------
    def sig(self, v: np.ndarray) -> np.ndarray:
        return (self.zreps @ v) % 2

    def is_cycle(self, v: np.ndarray) -> bool:
        return not ((self.HZ @ v) % 2).any()

    def is_stab(self, v: np.ndarray) -> bool:
        return in_span(v2i(v), self.rsHX_b, self.rsHX_p)

    def random_cycle(self, rng: np.random.Generator) -> np.ndarray:
        c = np.zeros(self.n, dtype=np.uint8)
        for kv in self.kerHZ:
            if rng.integers(2):
                c ^= kv
        return c

    @property
    def kappa(self) -> int:
        """dim rowspace H_X = ng - k/2 (the census generator count)."""
        return len(self.rsHX_b)

    def __repr__(self) -> str:  # pragma: no cover
        return (f"TowerCode({self.name}: {self.G.label()} "
                f"[[{self.n},{self.k}]])")


# ------------------------------------------------------------------- decks
def _fold_maps(cover_G: AbelianGroup, axis: int, newmod: int
               ) -> tuple[Callable, Callable]:
    """(fold, emb) closures for the index-2 fold along one axis."""
    def fold(e: tuple[int, ...]) -> tuple[int, ...]:
        u = list(e)
        u[axis] %= newmod
        return tuple(u)

    def emb(e: tuple[int, ...], s: int) -> tuple[int, ...]:
        u = list(e)
        u[axis] += newmod * s
        return tuple(u)

    return fold, emb


class AxisDeck:
    """A free Z2 deck along one axis: cover.G.orders[axis] = 2 * newmod,
    all other axes equal.  Carries the A32 Lemma-1 machinery — fold /
    sheet embeddings / carry system E v0 = RHS b / tau — with the
    constructor asserting every transport (chain map, stabilizer
    transport, twist-invariance im S <= ker E, sections).

    Rank-generic port of ``a32_tower_slice.Deck`` + the screen's
    ``make_deck``; the deck element is sigma = newmod * e_axis (central —
    abelian G), and twisted lifts cost nothing (Lemma 1 is twist-generic).
    """

    def __init__(self, cover: TowerCode, base: TowerCode, axis: int):
        assert cover.G.rank == base.G.rank, "rank mismatch"
        newmod = base.G.orders[axis]
        assert cover.G.orders[axis] == 2 * newmod, \
            f"not an index-2 fold on axis {axis}: " \
            f"{cover.G.orders} -> {base.G.orders}"
        for a in range(cover.G.rank):
            if a != axis:
                assert cover.G.orders[a] == base.G.orders[a], \
                    f"axis {a} differs: {cover.G.orders} vs {base.G.orders}"
        self.cover, self.base = cover, base
        self.axis, self.newmod = axis, newmod
        fold, emb = _fold_maps(cover.G, axis, newmod)
        self.fold, self.emb = fold, emb

        Gc, Gb = cover.G, base.G
        ngc, ngb = Gc.cardinality, Gb.cardinality
        P1 = np.zeros((ngb, ngc), dtype=np.uint8)
        for i, e in enumerate(Gc):
            P1[Gb.index(fold(e)), i] = 1
        Z = np.zeros_like(P1)
        self.P = np.block([[P1, Z], [Z, P1]]) % 2       # 1-chains fold
        self.P0 = P1.copy()                             # 0-cells fold
        self.EMB = []
        for s in (0, 1):
            E1 = np.zeros((ngc, ngb), dtype=np.uint8)
            for i, e in enumerate(Gb):
                E1[Gc.index(Gc.reduce(emb(e, s))), i] = 1
            Zb = np.zeros_like(E1)
            self.EMB.append(np.block([[E1, Zb], [Zb, E1]]) % 2)
        # sections: P o EMB_s = id
        for s in (0, 1):
            assert (self.P @ self.EMB[s] % 2 == np.eye(base.n,
                    dtype=np.uint8)).all(), "fold o emb != id"
        # chain-map transport: HZ_base o P = P0 o HZ_cover
        lhs = (base.HZ @ self.P) % 2
        rhs = (self.P0 @ cover.HZ) % 2
        assert (lhs == rhs).all(), "fold is not a chain map"
        # stabilizer transport: fold of each cover HX row = a base HX row
        for gidx in range(0, cover.ng, max(1, cover.ng // 24)):
            fb = (self.P @ cover.HX[gidx]) % 2
            gbar = base.G.index(fold(cover.G.from_index(gidx)))
            assert (fb == base.HX[gbar]).all(), "stab transport fails"
        # E-system for lifts of a base chain: E v0 = RHS beta
        self.E = (cover.HZ @ ((self.EMB[0] + self.EMB[1]) % 2)) % 2
        self.RHS = (cover.HZ @ self.EMB[1]) % 2
        # twist-invariance: tau(base stabilizer) is a cover cycle
        TAU = (self.EMB[0] + self.EMB[1]) % 2
        chk = (cover.HZ @ ((TAU @ base.HX.T) % 2)) % 2
        assert not chk.any(), "im S_base not <= ker E (twist argument fails)"
        self.TAU = TAU

    @property
    def sigma(self) -> tuple[int, ...]:
        """The deck element of the cover group (newmod on the fold axis)."""
        t = [0] * self.cover.G.rank
        t[self.axis] = self.newmod
        return tuple(t)

    def twisted(self) -> bool:
        """Does any lifted polynomial term carry the deck element?"""
        return bool(
            support_twisted(self.cover.A.support, self.axis, self.newmod)
            or support_twisted(self.cover.B.support, self.axis, self.newmod))

    def sheets(self, v: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        """(v0, v1) with v = EMB0 v0 + EMB1 v1."""
        v0 = (self.EMB[0].T @ v) % 2
        v1 = (self.EMB[1].T @ v) % 2
        rec = ((self.EMB[0] @ v0) + (self.EMB[1] @ v1)) % 2
        assert (rec == v).all(), "sheet decomposition failed"
        return v0, v1

    def slice_data(self, v: np.ndarray):
        """(shadow, overflow m, sheet0) with |v| = |shadow| + 2m."""
        v0, v1 = self.sheets(v)
        shadow = (v0 + v1) % 2
        m = int((v0 & v1).sum())
        assert int(v.sum()) == int(shadow.sum()) + 2 * m, "slice identity"
        assert (shadow == (self.P @ v) % 2).all(), "shadow != fold"
        return shadow, m, v0

    def lift(self, v0: np.ndarray, shadow: np.ndarray) -> np.ndarray:
        return ((self.EMB[0] @ v0)
                + (self.EMB[1] @ ((v0 + shadow) % 2))) % 2


# ------------------------------------------------ H1 maps and group action
def h1_map(deck: AxisDeck, tau: bool = False) -> np.ndarray:
    """Matrix of p_* (or tau_*): H1(cover) -> H1(base) (or reverse),
    in sig coordinates."""
    src = deck.cover if not tau else deck.base
    dst = deck.base if not tau else deck.cover
    S = np.array([src.sig(r) for r in src.xreps], dtype=np.uint8)
    op = (deck.P if not tau else deck.TAU)
    D = np.array([dst.sig((op @ r) % 2) for r in src.xreps], dtype=np.uint8)
    return (D.T @ gf2_inv(S.T)) % 2


def perm_for(code: TowerCode, t: tuple[int, ...]) -> np.ndarray:
    """Index array: (t . v)[i] = v[perm[i]] for the translation by t."""
    ng = code.ng
    perm = np.zeros(2 * ng, dtype=np.int64)
    for i, e in enumerate(code.G):
        j = code.G.index(code.G.sub(e, t))
        perm[i] = j
        perm[ng + i] = ng + j
    return perm


def translation_mat(code: TowerCode, t: tuple[int, ...]) -> np.ndarray:
    """Sig-space matrix of the translation by t acting on H1."""
    S = np.array([code.sig(r) for r in code.xreps], dtype=np.uint8)
    Sinv = gf2_inv(S.T)
    perm = perm_for(code, t)
    D = np.array([code.sig(r[perm]) for r in code.xreps], dtype=np.uint8)
    return (D.T @ Sinv) % 2


def translation_action(code: TowerCode) -> list[np.ndarray]:
    """Sig-space matrices of all |G| translations acting on H1."""
    S = np.array([code.sig(r) for r in code.xreps], dtype=np.uint8)
    Sinv = gf2_inv(S.T)
    mats = []
    for t in code.G:
        perm = perm_for(code, t)
        D = np.array([code.sig(r[perm]) for r in code.xreps],
                     dtype=np.uint8)
        mats.append((D.T @ Sinv) % 2)
    return mats


def translation_perms(code: TowerCode) -> list[np.ndarray]:
    return [perm_for(code, t) for t in code.G]


def rep_for(code: TowerCode, sig_int: int) -> np.ndarray:
    """A cycle with the prescribed H1 signature."""
    S = np.array([code.sig(r) for r in code.xreps], dtype=np.uint8)
    SinvT = gf2_inv(S.T)
    tvec = i2v(sig_int, code.k)
    coeff = (SinvT @ tvec) % 2
    v = np.zeros(code.n, dtype=np.uint8)
    for i in range(code.k):
        if coeff[i]:
            v ^= code.xreps[i]
    assert v2i(code.sig(v)) == sig_int
    return v


def canon_key(v: np.ndarray, perms: list[np.ndarray]) -> int:
    """Lexicographic-min bit-int over the translation orbit of v."""
    return min(v2i(v[p]) for p in perms)


def batch_keys(vecs: np.ndarray, perms: list[np.ndarray]) -> np.ndarray:
    """Lexicographic-min packed key over translation perms, per row
    (verbatim from a32_sectorAC_full.py)."""
    N, n = vecs.shape
    nwords = -(-(-(-n // 8)) // 8)  # ceil(ceil(n/8)/8)
    cur = None
    for p in perms:
        t = np.packbits(vecs[:, p], axis=1)
        pad = np.zeros((N, nwords * 8 - t.shape[1]), dtype=np.uint8)
        t64 = np.ascontiguousarray(
            np.concatenate([t, pad], axis=1)).view(">u8").reshape(N, nwords)
        if cur is None:
            cur = t64.copy()
        else:
            better = np.zeros(N, dtype=bool)
            tied = np.ones(N, dtype=bool)
            for w in range(nwords):
                better |= tied & (t64[:, w] < cur[:, w])
                tied &= t64[:, w] == cur[:, w]
            cur[better] = t64[better]
    return cur


def orbits(points: set[int], mats: list[np.ndarray]) -> list[set[int]]:
    """Orbit partition of a set of sig-ints under the given H1 matrices."""
    left = set(points)
    orbs = []
    while left:
        x = left.pop()
        orb = {x}
        stack = [x]
        while stack:
            y = stack.pop()
            for M in mats:
                z = apply_map(M, y)
                if z not in orb:
                    orb.add(z)
                    stack.append(z)
        left -= orb
        orbs.append(orb)
    return orbs


# --------------------------------------------------------- cost gate G1
def census_nodes(kappa: int, W: int, n_bases: int = 1) -> int:
    """Exact node count of the two-window coset-BZ census, balanced
    r-pair (r1 + r2 + 1 = W); per window sum_{w=1..r} C(kappa, w), per
    coset base.  Calibrated against the A30/A32/A33 anchors in
    :func:`validate_banked` (the A33 anchor reads x3.00 because that
    engine shares one walk across its 3 coset offsets — shared-walk
    kernels take multiple offsets in one walk)."""
    if W < 1:
        return 0
    r1 = (W - 1 + 1) // 2
    r2 = (W - 1) - r1
    tot = 0
    for r in (r1, r2):
        tot += sum(math.comb(kappa, w) for w in range(1, r + 1))
    return tot * n_bases


# ------------------------------------------------------------ lift fibers
def enumerate_lifts(deck: AxisDeck, beta: np.ndarray, cap: int,
                    kernel_cap: int = 16) -> dict[int, int]:
    """All v0 with E v0 = RHS beta and |v0 off supp(beta)| <= cap.

    Returns {canonical v0 int -> m2}, canonical under v0 -> v0 + beta
    (the deck translate).  Complete by the exact-off-support subset-sum
    argument (verbatim from a32_subclosures.py; the deep ordered-split
    lane for cap > 4 stays in the session scripts)."""
    assert cap <= 4, "MITM lane implemented to size 4"
    n = deck.base.n
    E_cols = [v2i(deck.E[:, j]) for j in range(n)]
    rhs = (deck.RHS @ beta) % 2
    rhs_i = v2i(rhs)
    bsupp = [int(j) for j in np.nonzero(beta)[0]]
    bmask = v2i(beta)
    bcols = [E_cols[j] for j in bsupp]
    bb, bp = rref_ints(bcols)
    rhs_res = reduce_int(rhs_i, bb, bp)
    offb = [j for j in range(n) if not (bmask >> j) & 1]
    red = {j: reduce_int(E_cols[j], bb, bp) for j in offb}
    by_val: dict[int, list[int]] = {}
    for j in offb:
        by_val.setdefault(red[j], []).append(j)
    hits_X: set[tuple[int, ...]] = set()
    if rhs_res == 0:
        hits_X.add(())
    if cap >= 1:
        for j in by_val.get(rhs_res, []):
            hits_X.add((j,))
    if cap >= 2:
        for j1 in offb:
            for j2 in by_val.get(rhs_res ^ red[j1], []):
                if j2 > j1:
                    hits_X.add((j1, j2))
    if cap >= 3:
        for j1, j2 in combinations(offb, 2):
            for j3 in by_val.get(rhs_res ^ red[j1] ^ red[j2], []):
                if j3 > j2:
                    hits_X.add((j1, j2, j3))
    if cap >= 4:
        pair_sum: dict[int, list[tuple[int, int]]] = {}
        for j1, j2 in combinations(offb, 2):
            pair_sum.setdefault(red[j1] ^ red[j2], []).append((j1, j2))
        for val, prs in pair_sum.items():
            for j3, j4 in pair_sum.get(rhs_res ^ val, []):
                for j1, j2 in prs:
                    if j2 < j3:
                        hits_X.add((j1, j2, j3, j4))
    out: dict[int, int] = {}
    for X in sorted(hits_X):
        cols = bsupp + list(X)
        b3: list[int] = []
        p3: list[int] = []
        h3: list[int] = []
        deps: list[int] = []
        for ci, j in enumerate(cols):
            cur, h = E_cols[j], 1 << ci
            for bb3, pp3, hh in zip(b3, p3, h3):
                if (cur >> pp3) & 1:
                    cur ^= bb3
                    h ^= hh
            if cur:
                b3.append(cur)
                p3.append((cur & -cur).bit_length() - 1)
                h3.append(h)
            else:
                deps.append(h)
        cur, hsel = rhs_i, 0
        for bb3, pp3, hh in zip(b3, p3, h3):
            if (cur >> pp3) & 1:
                cur ^= bb3
                hsel ^= hh
        if cur:
            continue
        assert len(deps) <= kernel_cap, f"kernel 2^{len(deps)} at X={X}"
        for kt in range(1 << len(deps)):
            sel = hsel
            for jj in range(len(deps)):
                if (kt >> jj) & 1:
                    sel ^= deps[jj]
            v0_int = 0
            for ci, j in enumerate(cols):
                if (sel >> ci) & 1:
                    v0_int |= 1 << j
            m2 = bin(v0_int & ~bmask).count("1")
            if m2 > cap:
                continue
            v0 = i2v(v0_int, n)
            assert not (((deck.E @ v0) + rhs) % 2).any(), "not a solution"
            canon = min(v0_int, v0_int ^ bmask)
            prev = out.get(canon)
            if prev is None or m2 < prev:
                out[canon] = m2
    return out


# ------------------------------------------------------------- rung screen
def liftable_codim(deck: AxisDeck) -> int:
    """Rank of the carry obstruction on H1(base): k_base - dim of the
    classes [b] whose carry system E v0 = RHS b is consistent.  Equals
    k_base - k_cover/2 on every screened rung (the delta / connecting-map
    law, A35 §3)."""
    F = np.array(kernel_basis(deck.E.T), dtype=np.uint8)
    if F.size == 0:
        return 0
    O = (F @ deck.RHS) % 2
    for ridx in (0, deck.base.ng // 2):
        assert not ((O @ deck.base.HX[ridx]) % 2).any(), \
            "obstruction does not vanish on a stabilizer row"
    M = np.array([(O @ r) % 2 for r in deck.base.xreps], dtype=np.uint8)
    return gf2_rank([v2i(r) for r in M])


def fiber_sample(deck: AxisDeck, W: int, rng: np.random.Generator,
                 want: int = 30, tries: int = 600) -> dict:
    """Restricted-fiber probe at the production regime: HEAVY stabilizer
    shadows (|beta| in [W-6, W], caps 0-3), sampled as localized row sums.
    Known limitation (A35 §7, falsified claim ledger): this shallow probe
    CANNOT estimate the census-population emptiness rates — it exists to
    exercise the enumerator and sample mu (lightest nonzero stabilizer
    seen).  RNG call order matches a35_generality_screen.fiber_sample
    exactly at rank 2 (same stream -> same numbers)."""
    base = deck.base
    ords = base.G.orders
    tried = empty = 0
    by_cap: dict[int, list[int]] = {}
    mu = int(base.HX[0].sum())
    for _ in range(tries):
        if tried >= want:
            break
        j = int(rng.integers(2, 6))
        g0 = tuple(int(rng.integers(o)) for o in ords)
        beta = np.zeros(base.n, dtype=np.uint8)
        for _ in range(j):
            deltas = [int(rng.integers(-2, 3)) for _ in ords]
            gi = base.G.index(tuple((g0[a] + deltas[a]) % ords[a]
                                    for a in range(len(ords))))
            beta = (beta + base.HX[gi]) % 2
        wb = int(beta.sum())
        if wb:
            mu = min(mu, wb)
        if not (W - 6 <= wb <= W):
            continue
        cap = (W - wb) // 2
        lifts = enumerate_lifts(deck, beta, cap)
        tried += 1
        if not lifts:
            empty += 1
        by_cap.setdefault(cap, [0, 0])
        by_cap[cap][0] += 0 if lifts else 1
        by_cap[cap][1] += 1
    return {"tried": tried, "empty": empty,
            "empty_rate": round(empty / tried, 3) if tried else None,
            "by_cap": {str(c): f"{e}/{t}" for c, (e, t)
                       in sorted(by_cap.items())},
            "mu_sampled": mu}


def screen_rung(cover: TowerCode, base: TowerCode, axis: int, W_eff: int,
                rng: Optional[np.random.Generator] = None,
                do_fibers: bool = True) -> dict:
    """Measure every structural hypothesis of one rung (A35 §1 layers
    L0-L4): Lemma-1 transports (constructor asserts), twist, (R) via
    k-preservation (A12), the H1 rank lattice + exactness + sigma*, the
    liftable-cycle codimension, and (optionally) the sampled fiber probe.
    """
    deck = AxisDeck(cover, base, axis)      # Lemma 1 asserts inside
    r: dict = {"axis": _VAR_NAMES[axis], "fold_to": base.G.orders[axis],
               "lemma1": "PASS (constructor asserts)"}
    r["twisted"] = deck.twisted()
    r["k_cover"], r["k_base"] = cover.k, base.k
    r["R_holds"] = bool(cover.k == base.k)   # A12: (R) <=> k preserved
    if cover.k > 0 and base.k > 0:
        Mp = h1_map(deck)
        Mt = h1_map(deck, tau=True)
        r["rank_p"] = gf2_rank([v2i(c) for c in Mp.T])
        r["rank_tau"] = gf2_rank([v2i(c) for c in Mt.T])
        imP = colspace(Mp)
        imT = colspace(Mt)
        kerP = kernel_ints(Mp)
        kerT = kernel_ints(Mt)
        r["exact_cover"] = bool(span_eq(imT, kerP))   # im tau* = ker p*
        r["exact_base"] = bool(span_eq(kerT, imP))    # ker tau* = im p*
        St = translation_mat(cover, deck.sigma)
        r["sigma_id"] = bool((St == np.eye(cover.k, dtype=np.uint8)).all())
        r["_imP"] = imP
        r["_kerP"] = kerP
    else:
        r["rank_p"] = r["rank_tau"] = 0
        r["exact_cover"] = r["exact_base"] = None
        r["sigma_id"] = None
        r["_imP"] = []
        r["_kerP"] = []
    r["codim_lift"] = liftable_codim(deck)
    if do_fibers:
        assert rng is not None, "fiber sampling needs an rng"
        r["fibers"] = fiber_sample(deck, W_eff, rng)
    return r


def gate_verdict(nodes_bottom: int, cap_max: int) -> str:
    """The A35 §5 cost verdict.  GREEN/AMBER/RED are COST verdicts on the
    quantified gates — never distance claims."""
    return ("GREEN" if nodes_bottom <= 2e11 and cap_max <= 8 else
            "AMBER" if nodes_bottom <= 1e14 and cap_max <= 12 else "RED")


def screen_tower(spec: dict, rng: Optional[np.random.Generator] = None,
                 log: Callable[[str], None] = lambda s: None) -> dict:
    """The A35 tower screen, rank-generic: build the level chain by
    literal polynomial descent, screen every rung and adjacent pair,
    and price the cost gates — WITHOUT running any closure.

    spec keys: name; top = (orders, A, B); folds = [(axis, newmod), ...];
    d_top (None = unknown); optional d_mid, W_eff, W_list, fibers, tag,
    notes.  Orders are tuples of any rank; folds name an axis index."""
    name = spec["name"]
    orders, As, Bs = spec["top"]
    orders = tuple(int(o) for o in orders)
    d_top = spec.get("d_top")
    W_eff = spec.get("W_eff") or ((d_top - 2) if d_top else None)
    log(f"--- {name}: top {orders} A={As!r} B={Bs!r}")

    G_top = AbelianGroup(orders)
    levels = [(orders, _as_support(As, G_top), _as_support(Bs, G_top))]
    for axis, newmod in spec["folds"]:
        plm, pA, pB = levels[-1]
        assert plm[axis] == 2 * newmod, \
            f"{name}: fold {axis}/{newmod} is not an index-2 fold of {plm}"
        nlm = tuple(newmod if a == axis else plm[a]
                    for a in range(len(plm)))
        levels.append((nlm, fold_support(pA, axis, newmod),
                       fold_support(pB, axis, newmod)))

    codes: list[TowerCode] = []
    lvl_rows = []
    for i, (glm, tA, tB) in enumerate(levels):
        code = TowerCode(f"{name}/L{i}", glm, tA, tB)
        codes.append(code)
        odd_ok = all(int(kv.sum()) % 2 == 0 for kv in code.kerHZ)
        kappa = code.kappa
        assert kappa == code.ng - code.k // 2, "kappa formula fails"
        lvl_rows.append({"lm": list(glm), "n": code.n, "k": code.k,
                         "wA": len(tA), "wB": len(tB),
                         "parity_scope": bool(len(tA) % 2 and len(tB) % 2),
                         "cycles_all_even": bool(odd_ok), "kappa": kappa})
        log(f"    L{i} {glm}: [[{code.n},{code.k}]] "
            f"|A|={len(tA)} |B|={len(tB)} kappa={kappa} "
            f"even-cycles={odd_ok}")

    W_use = W_eff if W_eff else (2 * (spec.get("d_mid") or 0) - 2 or 16)
    rungs = []
    for i, (axis, newmod) in enumerate(spec["folds"]):
        r = screen_rung(codes[i], codes[i + 1], axis, W_use, rng=rng,
                        do_fibers=spec.get("fibers", True))
        rungs.append(r)
        log(f"    rung {i} ({r['axis']}->{newmod}): twisted={r['twisted']} "
            f"R={r['R_holds']} (k {r['k_cover']}->{r['k_base']}) "
            f"rank p*={r['rank_p']} tau*={r['rank_tau']} "
            f"exact(cov/base)={r['exact_cover']}/{r['exact_base']} "
            f"sigma*=id:{r['sigma_id']} codim_lift={r['codim_lift']}"
            + (f" fibers: {r['fibers']['empty']}/{r['fibers']['tried']} "
               f"empty (rate {r['fibers']['empty_rate']})"
               if "fibers" in r else ""))

    # adjacent-pair two-level lattice (in H1 of the middle code)
    pairs = []
    for i in range(len(rungs) - 1):
        S = rungs[i]["_imP"]           # im p_top*  in H1(level i+1)
        K = rungs[i + 1]["_kerP"]      # ker p_bot* in H1(level i+1)
        mid, bot = codes[i + 1], codes[i + 2]
        if mid.k == 0 or not S:
            pairs.append({"i": i, "note": "H1(mid) trivial"})
            continue
        dS = gf2_rank(S)
        dK = gf2_rank(K)
        dSK = dS + dK - gf2_rank(list(S) + list(K))
        Mb = h1_map(AxisDeck(mid, bot, spec["folds"][i + 1][0])) \
            if bot.k > 0 else None
        if Mb is not None:
            Sb, _ = rref_ints(list(S))
            Wimg = [apply_map(Mb, s) for s in span_points(Sb)]
            dW = gf2_rank(Wimg)
            Wb2, Wp2 = rref_ints([w for w in Wimg if w])
            pre = preimage_basis(Mb, Wb2, Wp2)
            reach = bool(span_eq(pre, Sb)) if dSK == dK else False
        else:
            dW, reach = 0, None
        K_in_S = bool(all(in_span(x, *rref_ints(list(S))) for x in K)) \
            if K else True
        pairs.append({"i": i, "dim_S": dS, "dim_K": dK, "dim_SK": dSK,
                      "K_in_S": K_in_S, "one_branch": bool(dSK == 0),
                      "dim_W": dW, "reach_preimage_eq_S": reach})
        log(f"    pair ({i},{i+1}): dim S={dS} K={dK} S^K={dSK} "
            f"K<=S:{K_in_S} one-branch:{dSK == 0} dim W={dW} "
            f"preimage=S:{reach}")

    # cost gates: one row per certification question (W value)
    W_list = spec.get("W_list") or ([W_eff] if W_eff else [])
    mu = min((r.get("fibers", {}).get("mu_sampled") or 99)
             for r in rungs) if rungs else 6
    mu = min(mu, lvl_rows[-1]["wA"] + lvl_rows[-1]["wB"])
    costs = []
    for W in W_list:
        lvl_nodes = [census_nodes(lv["kappa"], W) for lv in lvl_rows]
        nb, nt = lvl_nodes[-1], lvl_nodes[0]
        cap_max = (W - mu) // 2
        verdict = gate_verdict(nb, cap_max)
        costs.append({
            "W": W, "mu": mu, "cap_max": cap_max, "verdict": verdict,
            "log10_nodes_per_level":
                [round(math.log10(x), 1) if x else None for x in lvl_nodes],
            "win_factor_vs_top": round(nt / nb, 1) if nb else None,
        })
        log(f"    costs @W={W}: per-level log10 nodes = "
            f"{costs[-1]['log10_nodes_per_level']} "
            f"win={costs[-1]['win_factor_vs_top']:.1e}x "
            f"cap_max={cap_max} -> {verdict}")

    for r in rungs:                       # strip non-JSON internals
        r.pop("_imP", None)
        r.pop("_kerP", None)
    return {"name": name, "tag": spec.get("tag", ""),
            "d_top": d_top, "levels": lvl_rows, "rungs": rungs,
            "pairs": pairs, "costs": costs, "notes": spec.get("notes", "")}


# ------------------------------------------------------ tower inventory
def tower_inventory(orders: Sequence[int]) -> dict:
    """C1 + the tower shape for a group given as an order tuple: per-axis
    2-adic depth, the odd part, and the deepest iterated-Z2 fold chain
    (levels can be taken one Z2 at a time in any order — A33/A35)."""
    orders = tuple(int(o) for o in orders)
    v2s = []
    odd = []
    for o in orders:
        v = 0
        q = o
        while q % 2 == 0:
            q //= 2
            v += 1
        v2s.append(v)
        odd.append(q)
    card = math.prod(orders)
    folds: list[tuple[int, int]] = []
    cur = list(orders)
    for a, v in enumerate(v2s):
        for _ in range(v):
            cur[a] //= 2
            folds.append((a, cur[a]))
    return {"orders": list(orders), "card": card,
            "C1_two_divides": card % 2 == 0,
            "v2_per_axis": v2s, "odd_part_per_axis": odd,
            "odd_part": math.prod(odd), "two_part": card // math.prod(odd),
            "depth": sum(v2s), "fold_chain": folds}


# ------------------------------------------------------------- A35 docket
#: The A35 generality-screen docket (banked reference: data/a35/screen.json
#: on branch claude/tower-slice-calculus-generalize-410ed1).  Regression
#: battery for any engine change; order matters (shared RNG stream).
A35_DOCKET: list[dict] = [
    dict(name="bravyi360", tag="VAL-A32",
         top=((30, 6), "x^9 + y + y^2", "y^3 + x^25 + x^26"),
         folds=[(1, 3), (0, 15)], d_top=24, d_mid=12,
         notes="A32 instance: mixed-axis, top rung (R)-fails, d=24 CLOSED"),
    dict(name="ibm288Y", tag="VAL-A33",
         top=((18, 8), "1 + x*y^4 + x^14*y", "1 + x*y^2 + x^2*y^7"),
         folds=[(1, 4), (1, 2)], d_top=20, d_mid=10,
         notes="A33 instance: same-axis, all-(R), d=20 CLOSED"),
    dict(name="gross_xx", tag="NEW",
         top=((12, 6), "x^3 + y + y^2", "y^3 + x + x^2"),
         folds=[(0, 6), (0, 3)], d_top=12, d_mid=6,
         notes="retrospective target: gross over bb72 over (3,6)"),
    dict(name="bb288_yxx", tag="NEW",
         top=((12, 12), "x^3 + y^2 + y^7", "y^3 + x + x^2"),
         folds=[(1, 6), (0, 6), (0, 3)], d_top=18, d_mid=12,
         notes="published record [[288,12,18]]; y-quotient IS gross "
               "(y^7 = y mod 6) -> 4-level tower to n=36"),
    dict(name="c37x_360420", tag="NEW",
         top=((30, 6), "1 + y + x", "y^4 + x + x^11*y^2"),
         folds=[(0, 15), (1, 3)], d_top=20, d_mid=10,
         notes="A30 doubled code 37a70e02:x = [[360,4,20]]; mixed-axis"),
    dict(name="e5e50yy_360420", tag="NEW",
         top=((15, 12), "1 + y + x", "y^4 + x^8*y^2 + x^13"),
         folds=[(1, 6), (1, 3)], d_top=20, d_mid=10,
         notes="A30 doubled code 5e50a976:y = [[360,4,20]]; same-axis"),
    dict(name="c37xx_720", tag="FRONTIER",
         top=((60, 6), "1 + y + x", "y^4 + x + x^11*y^2"),
         folds=[(0, 30), (0, 15), (1, 3)], d_top=None, d_mid=20,
         W_eff=18, W_list=[18, 22, 30, 38],
         notes="rung-2 re-double [[720,4,?]]; freeze-vs-double open "
               "(A14 SS13 / A33 SS6 ranking #1); W=18 certifies d=20 if "
               "frozen, W=38 = full doubling budget"),
    dict(name="a8_336", tag="NEW",
         top=((12, 14), "1 + y + x^3*y^3", "1 + x + x^2*y^7"),
         folds=[(0, 6), (1, 7)], d_top=12, d_mid=6,
         notes="A8/A29 [[336,12,12]] over its [[168,12,6]] base, then y"),
    dict(name="bravyi756", tag="FRONTIER",
         top=((21, 18), "x^3 + y^10 + y^17", "y^5 + x^3 + x^19"),
         folds=[(1, 9)], d_top=34, fibers=False,
         notes="[[756,16,<=34]]: ONE deck only (v2(18)=1, 21 odd); "
               "bottom n=378"),
    dict(name="cover300", tag="DEGENERATE-1LVL",
         top=((5, 30), "1 + y + x", "x*y^6 + x*y^10 + x^2*y^12"),
         folds=[(1, 15)], d_top=16, d_mid=8,
         notes="A15's [[300,8,16]]: v2(30)=1 -> one-level = the A15/A30 "
               "architecture (already closed there)"),
    dict(name="pair72", tag="NEW-TINY",
         top=((6, 6), "x^2 + y + y^3", "1 + x + y^2"),
         folds=[(0, 3), (1, 3)], d_top=8, d_mid=4,
         notes="[[72,4,8]] over [[36,4,4]] over (3,3); smallest 2-level"),
]

#: |G| odd: no free Z2 deck exists at all (layer L1 / wall W1).
A35_NO_DECK: list[tuple[str, tuple[int, int], str]] = [
    ("bb90", (15, 3), "[[90,8,10]] -- A19 survey: 'not a Z2-cover at all'"),
    ("bb98", (7, 7), "[[98,6,12]] (A16 host)"),
    ("f2a6_base", (5, 15), "[[150,8,8]] -- the A24 odd-|G| lane"),
]


# --------------------------------------------------- the falsify-first gate
def _assert_bravyi360(res: dict) -> None:
    r0, r1 = res["rungs"]
    assert [lv["k"] for lv in res["levels"]] == [12, 8, 8]
    assert (not r0["R_holds"]) and r1["R_holds"]
    assert r0["rank_p"] == 6 and r1["rank_p"] == 4 and r1["rank_tau"] == 4
    assert r1["exact_cover"] and r1["exact_base"]
    assert r0["sigma_id"] is False and r0["twisted"] and r1["twisted"]
    p = res["pairs"][0]
    assert p["dim_S"] == 6 and p["dim_K"] == 4 and p["dim_SK"] == 4 \
        and p["K_in_S"] and p["dim_W"] == 2 and p["reach_preimage_eq_S"]


def _assert_ibm288Y(res: dict) -> None:
    r0, r1 = res["rungs"]
    assert [lv["k"] for lv in res["levels"]] == [8, 8, 8]
    assert r0["R_holds"] and r1["R_holds"]
    assert r0["rank_p"] == r0["rank_tau"] == 4 == r1["rank_p"] \
        == r1["rank_tau"]
    assert r0["exact_cover"] and r0["exact_base"] \
        and r1["exact_cover"] and r1["exact_base"]
    assert r0["sigma_id"] and r1["sigma_id"]
    assert r0["twisted"] and r1["twisted"]
    p = res["pairs"][0]
    assert p["dim_S"] == 4 and p["dim_SK"] == 0 and p["one_branch"] \
        and p["dim_W"] == 4


def validate_banked(data_dir, rng: Optional[np.random.Generator] = None,
                    log: Callable[[str], None] = lambda s: None) -> dict:
    """The falsify-first regression gate (charter A38 §6.0): reproduce the
    banked A32/A33 structure, the census node anchors, the A19 deck-survey
    k-verdicts, and — where the banked A32 census files are present under
    ``data_dir`` — the bit-level sector-C fiber layer (397 orbit fibers,
    4,132 lifts, exact m2-histogram), before anything new runs.

    Raises AssertionError on ANY deviation.  Returns the collected
    numbers.  ``data_dir`` is the lab data root (…/experiments/bb_lab/
    data); pass rng seeded 20260811 to reproduce the banked screen's
    sampled fiber sections bit-for-bit."""
    import gzip
    import json
    from pathlib import Path

    data_dir = Path(data_dir)
    if rng is None:
        rng = np.random.default_rng(20260811)
    out: dict = {"anchors": {}}

    # [0] census node formula vs banked anchors
    a32_w = census_nodes(41, 22, 3)   # A32 W-coset <=22: reported 1.9e10
    a32_s = census_nodes(41, 22, 1)   # A32 stab <=22:    reported 6.3e9
    a33_h = census_nodes(68, 18, 3)   # A33 H5 pass:      reported 6.62e10
    for tag, got, want in [("a32_wcoset22", a32_w, 1.9e10),
                           ("a32_stab22", a32_s, 6.3e9),
                           ("a33_h5", a33_h, 6.62e10)]:
        ratio = got / want
        log(f"    {tag}: formula {got:.2e} vs banked {want:.2e} "
            f"(x{ratio:.2f})")
        assert 0.3 < ratio < 3.5, f"node formula off at {tag}"
        out["anchors"][tag] = {"formula": float(got), "banked": want,
                               "ratio": round(ratio, 2)}
    # the exact-anchor identities the A35 note reports (x1.00/x1.01/x3.00):
    assert out["anchors"]["a32_wcoset22"]["ratio"] == 1.00
    assert out["anchors"]["a32_stab22"]["ratio"] == 1.01
    assert out["anchors"]["a33_h5"]["ratio"] == 3.00, \
        "the A33 shared-walk x3.00 accounting datum moved"

    # [0b] fiber enumerator vs banked A32 sector-C fibers (bit-level)
    stab_f = data_dir / "a32" / "gb_census_stab.jsonl"
    if stab_f.exists():
        GB = TowerCode("GB", (15, 3), "x^9 + y + y^2", "1 + x^10 + x^11")
        BY = TowerCode("BY", (30, 3), "x^9 + y + y^2", "1 + x^25 + x^26")
        deck_x = AxisDeck(BY, GB, 0)
        reps: dict[int, np.ndarray] = {}
        for line in stab_f.open():
            r = json.loads(line)
            if r["w"] in (14, 16) and r["canon"] not in reps:
                v = np.zeros(GB.n, dtype=np.uint8)
                v[r["support"]] = 1
                reps[r["canon"]] = v
        assert len(reps) == 397, f"{len(reps)} orbit reps != 64+333"
        tot = empty = 0
        hist: dict[int, int] = {}
        for beta in reps.values():
            cap = (22 - int(beta.sum())) // 2
            lifts = enumerate_lifts(deck_x, beta, cap)
            tot += len(lifts)
            empty += 0 if lifts else 1
            for m2 in lifts.values():
                hist[m2] = hist.get(m2, 0) + 1
        assert tot == 4132, f"sector-C lift total {tot} != banked 4,132"
        assert hist == {0: 416, 1: 554, 2: 561, 3: 1639, 4: 962}, hist
        log(f"    397 orbit fibers (caps 4/3): {tot} lifts == banked "
            f"4,132 EXACTLY, m2-hist matches; {empty} empty fibers")
        out["sectorC_refiber"] = {
            "fibers": 397, "lifts": tot, "empty": empty,
            "m2_hist": {str(k): v for k, v in sorted(hist.items())}}
        band: dict[int, tuple[int, int]] = {}
        for line in gzip.open(data_dir / "a32" / "sectorAC_C18to22.jsonl.gz"):
            r = json.loads(line)
            e, t = band.get(r["wbeta"], (0, 0))
            band[r["wbeta"]] = (e + (0 if r.get("m2_hist") else 1), t + 1)
        want_band = {18: (319, 1733), 20: (5635, 10602),
                     22: (55555, 64619)}
        assert band == want_band, f"banked heavy-band rates moved: {band}"
        out["sectorC_heavy_empty_rates_banked"] = {
            str(w): {"empty": e, "fibers": t, "rate": round(e / t, 3)}
            for w, (e, t) in sorted(band.items())}
    else:
        out["sectorC_refiber"] = "SKIPPED (banked a32 census absent)"
        log("    [0b] SKIPPED: banked a32 census files absent")

    # [1] the two validation towers + the A19 k-verdict towers
    by_name = {s["name"]: s for s in A35_DOCKET}
    for name in ("bravyi360", "ibm288Y", "gross_xx", "bb288_yxx"):
        res = screen_tower(by_name[name], rng=rng, log=log)
        out[name] = res
        if name == "bravyi360":
            _assert_bravyi360(res)
            log("    bravyi360: ALL banked A32 structure REPRODUCED")
        elif name == "ibm288Y":
            _assert_ibm288Y(res)
            log("    ibm288Y: ALL banked A33 structure REPRODUCED")
        elif name == "gross_xx":
            assert [lv["k"] for lv in res["levels"]] == [12, 12, 8], \
                "A19: gross x-deck R-holds (bb72), bb72 decks jump 12->8"
            log("    gross_xx: A19 deck-survey k pattern REPRODUCED")
        elif name == "bb288_yxx":
            assert res["levels"][1]["lm"] == [12, 6]
            assert res["rungs"][0]["R_holds"], "A19: bb288 all-(R) decks"
            log("    bb288_yxx: y-quotient = (12,6) gross frame confirmed")

    # [2] no-deck rows (structural, no compute)
    for nm, lm, note in A35_NO_DECK:
        inv = tower_inventory(lm)
        assert not inv["C1_two_divides"], f"{nm} should be odd-|G|"
        assert inv["depth"] == 0
    out["no_deck"] = [{"name": nm, "lm": list(lm), "note": note}
                      for nm, lm, note in A35_NO_DECK]
    log("    no-deck rows: |G| odd verified (bb90/bb98/f2a6_base)")

    return out
