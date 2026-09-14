"""Lifted-product (LP) CSS codes over abelian group algebras R = F₂[G].

Conventions (agree with bb_lab.checks for the 1×1 case):
  * G = Z_{n₁} × … × Z_{n_r}; elements indexed 0..|G|-1 (row-major).
  * A polynomial is a frozenset of group indices (its support).
  * M(a)[g, h] = 1  iff  (g − h) ∈ supp(a)   ("circulant" of a).
  * An R-matrix is a dict {(i, j): poly}; missing entries are 0.

For A ∈ R^{m_a×n_a}, B ∈ R^{m_b×n_b} the lifted product is
  H_X = [ A ⊗ I_{n_b} | I_{m_a} ⊗ B† ],   H_Z = [ I_{n_a} ⊗ B | A† ⊗ I_{m_b} ],
with † the conjugate transpose (g ↦ −g on supports). Expanding each
R-entry via M(·) gives binary matrices; CSS commutation follows from
M(a)M(b) = M(b)M(a) (R commutative) and M(a)ᵀ = M(a*).

Special case used throughout: LP(A, b) with A square a×a and b 1×1, i.e.
  H_X = [ A | b* I_a ],   H_Z = [ b I_a | A† ]     (BB codes are a = 1).
"""
from __future__ import annotations

from dataclasses import dataclass, field
from functools import reduce as _reduce
from itertools import product
from typing import Iterable, Mapping

import numpy as np

Poly = frozenset  # of group indices


class Group:
    """Finite abelian group Z_{n1} × … × Z_{nr}, elements as indices."""

    def __init__(self, orders: Iterable[int]):
        self.orders = tuple(int(o) for o in orders)
        self.N = int(np.prod(self.orders))
        self._tuples = list(product(*[range(o) for o in self.orders]))
        self._index = {t: i for i, t in enumerate(self._tuples)}
        # precompute add/neg tables lazily for speed
        self._neg = [self._index[tuple((-c) % o for c, o in zip(t, self.orders))]
                     for t in self._tuples]

    def tup(self, i: int) -> tuple[int, ...]:
        return self._tuples[i]

    def idx(self, t: Iterable[int]) -> int:
        return self._index[tuple(c % o for c, o in zip(t, self.orders))]

    def add(self, i: int, j: int) -> int:
        a, b = self._tuples[i], self._tuples[j]
        return self._index[tuple((x + y) % o for x, y, o in zip(a, b, self.orders))]

    def sub(self, i: int, j: int) -> int:
        a, b = self._tuples[i], self._tuples[j]
        return self._index[tuple((x - y) % o for x, y, o in zip(a, b, self.orders))]

    def neg(self, i: int) -> int:
        return self._neg[i]

    def mono(self, *exps: int) -> Poly:
        return frozenset([self.idx(exps)])

    def poly(self, *terms: Iterable[int]) -> Poly:
        return frozenset(self.idx(t) for t in terms)

    def conj(self, p: Poly) -> Poly:
        return frozenset(self._neg[g] for g in p)

    def mul(self, p: Poly, q: Poly) -> Poly:
        acc: dict[int, int] = {}
        for g in p:
            for h in q:
                s = self.add(g, h)
                acc[s] = acc.get(s, 0) ^ 1
        return frozenset(s for s, c in acc.items() if c)

    def padd(self, p: Poly, q: Poly) -> Poly:
        return frozenset(p ^ q)

    def circulant(self, p: Poly) -> np.ndarray:
        """|G|×|G| binary matrix M[g, h] = [ (g − h) ∈ p ]."""
        N = self.N
        M = np.zeros((N, N), dtype=np.uint8)
        for h in range(N):
            for s in p:
                M[self.add(h, s), h] ^= 1
        return M

    def poly_str(self, p: Poly, names: str = "xyzw") -> str:
        if not p:
            return "0"
        terms = []
        for g in sorted(p):
            t = self.tup(g)
            parts = [f"{names[i]}^{e}" if e != 1 else names[i]
                     for i, e in enumerate(t) if e != 0]
            terms.append("*".join(parts) if parts else "1")
        return " + ".join(terms)

    def parse(self, s: str, names: str = "xyzw") -> Poly:
        """Parse '1 + x^3 + x*y^2' into a Poly."""
        s = s.replace(" ", "")
        if s in ("0", ""):
            return frozenset()
        terms = s.split("+")
        out: set[int] = set()
        for t in terms:
            exps = [0] * len(self.orders)
            if t == "1":
                pass
            else:
                for fac in t.split("*"):
                    if "^" in fac:
                        v, e = fac.split("^")
                        e = int(e)
                    else:
                        v, e = fac, 1
                    exps[names.index(v)] += e
            g = self.idx(exps)
            out ^= {g}
        return frozenset(out)


RMatrix = Mapping[tuple[int, int], Poly]


def r_shape(A: RMatrix) -> tuple[int, int]:
    m = 1 + max(i for i, _ in A)
    n = 1 + max(j for _, j in A)
    return m, n


def expand(G: Group, A: RMatrix, shape: tuple[int, int]) -> np.ndarray:
    """Binary expansion of an R-matrix: (m|G|) × (n|G|)."""
    m, n = shape
    N = G.N
    H = np.zeros((m * N, n * N), dtype=np.uint8)
    for (i, j), p in A.items():
        if p:
            H[i * N:(i + 1) * N, j * N:(j + 1) * N] ^= G.circulant(p)
    return H


def r_conj_transpose(G: Group, A: RMatrix) -> dict:
    return {(j, i): G.conj(p) for (i, j), p in A.items()}


def r_kron_left(A: RMatrix, shape_A: tuple[int, int], n: int) -> tuple[dict, tuple[int, int]]:
    """A ⊗ I_n over R: index ((i, s), (j, s))."""
    ma, na = shape_A
    out = {}
    for (i, j), p in A.items():
        for s in range(n):
            out[(i * n + s, j * n + s)] = p
    return out, (ma * n, na * n)


def r_kron_right(m: int, B: RMatrix, shape_B: tuple[int, int]) -> tuple[dict, tuple[int, int]]:
    """I_m ⊗ B over R: index ((s, i), (s, j))."""
    mb, nb = shape_B
    out = {}
    for s in range(m):
        for (i, j), p in B.items():
            out[(s * mb + i, s * nb + j)] = p
    return out, (m * mb, m * nb)


@dataclass
class CSSCode:
    """Duck-types bb_lab.checks.CheckMatrices for bb_lab.sat_distance."""

    H_X: np.ndarray
    H_Z: np.ndarray
    meta: dict = field(default_factory=dict)

    @property
    def num_qubits(self) -> int:
        return int(self.H_X.shape[1])

    @property
    def n(self) -> int:
        return self.num_qubits

    def check_css(self) -> None:
        P = (self.H_X.astype(np.int64) @ self.H_Z.T.astype(np.int64)) & 1
        if P.any():
            raise AssertionError("H_X H_Zᵀ ≠ 0: not a CSS code")

    def weights(self) -> dict:
        return {
            "X_check_weights": sorted(set(self.H_X.sum(1).tolist())),
            "Z_check_weights": sorted(set(self.H_Z.sum(1).tolist())),
            "qubit_degrees": sorted(set((self.H_X.sum(0) + self.H_Z.sum(0)).tolist())),
        }


def lifted_product(G: Group, A: RMatrix, shape_A: tuple[int, int],
                   B: RMatrix, shape_B: tuple[int, int]) -> CSSCode:
    ma, na = shape_A
    mb, nb = shape_B
    Ad = r_conj_transpose(G, A)
    Bd = r_conj_transpose(G, B)
    left_X, shL = r_kron_left(A, (ma, na), nb)          # (ma nb) × (na nb)
    right_X, shR = r_kron_right(ma, Bd, (nb, mb))       # (ma nb) × (ma mb)
    left_Z, shLZ = r_kron_right(na, B, (mb, nb))        # (na mb) × (na nb)
    right_Z, shRZ = r_kron_left(Ad, (na, ma), mb)       # (na mb) × (ma mb)
    HX = np.concatenate([expand(G, left_X, shL), expand(G, right_X, shR)], axis=1)
    HZ = np.concatenate([expand(G, left_Z, shLZ), expand(G, right_Z, shRZ)], axis=1)
    code = CSSCode(HX, HZ, meta={"G": G.orders, "shape_A": shape_A, "shape_B": shape_B,
                                 "n_left": na * nb * G.N, "n_right": ma * mb * G.N})
    return code


def lp_square(G: Group, A: RMatrix, a: int, b: Poly) -> CSSCode:
    """LP(A, b): H_X = [A | b* I_a], H_Z = [b I_a | A†]."""
    code = lifted_product(G, A, (a, a), {(0, 0): b}, (1, 1))
    code.meta.update({"A": {k: sorted(v) for k, v in A.items()}, "b": sorted(b), "a": a})
    return code


def two_block(G: Group, A_bin: np.ndarray, B_bin: np.ndarray) -> CSSCode:
    """Generic two-block code H_X=[A|B], H_Z=[Bᵀ|Aᵀ] from commuting binary A, B."""
    HX = np.concatenate([A_bin, B_bin], axis=1).astype(np.uint8)
    HZ = np.concatenate([B_bin.T, A_bin.T], axis=1).astype(np.uint8)
    return CSSCode(HX, HZ)


def monomial_decomposition(G: Group, A: RMatrix, a: int) -> list[dict]:
    """Split a (3,3)-regular monomial-count matrix A into monomial-permutation
    matrices A₁, A₂, A₃ (each with exactly one monomial per row and column).

    Uses repeated bipartite perfect matchings on the multigraph
    rows–cols with multiplicity = #monomials of the entry.
    """
    import networkx as nx
    remaining = {k: set(v) for k, v in A.items() if v}
    parts = []
    while any(remaining.values()):
        Bg = nx.Graph()
        rows = [("r", i) for i in range(a)]
        cols = [("c", j) for j in range(a)]
        Bg.add_nodes_from(rows, bipartite=0)
        Bg.add_nodes_from(cols, bipartite=1)
        for (i, j), s in remaining.items():
            if s:
                Bg.add_edge(("r", i), ("c", j))
        M = nx.bipartite.maximum_matching(Bg, top_nodes=rows)
        part = {}
        for i in range(a):
            cj = M.get(("r", i))
            if cj is None:
                raise ValueError("no perfect matching: A is not regular")
            j = cj[1]
            g = min(remaining[(i, j)])
            remaining[(i, j)].remove(g)
            part[(i, j)] = frozenset([g])
        parts.append(part)
    return parts


def perm_of_monomial_matrix(G: Group, P: RMatrix, a: int) -> np.ndarray:
    """For a monomial-permutation matrix P (one monomial per row/col), return
    the permutation π on [a]×G such that the binary expansion M(P) has
    M(P)[π(c), c] = 1 (i.e. P maps column c to row π(c))."""
    N = G.N
    pi = np.empty(a * N, dtype=np.int64)
    for (i, j), p in P.items():
        (s,) = tuple(p)
        for h in range(N):
            pi[j * N + h] = i * N + G.add(h, s)
    return pi
