"""Fast bit-packed builders for lifted-product classes over any finite
abelian group G (elements 0..N-1 with an addition table).

Classes:
  * lp_ab  : LP(A, b) with A a×a over R, b 1×1.   n = 2aN.
             H_X = [A | b* I_a],  H_Z = [b I_a | A†]
  * lp_AB  : LP(A, B) with A, B both a×a.          n = 2a²N.
             H_X = [A⊗I | I⊗B†],  H_Z = [I⊗B | A†⊗I]
Polynomials are tuples of group indices.  A is a dict {(i,j): poly}.

Row conventions (binary expansion M(p)[g,h] = [g−h ∈ p]):
  row g of M(p) has ones at h = g − s (s ∈ p);
  row g of M(p)† = M(p*)ᵀ... we build H_Z rows directly:
    (I⊗B) row (j, g): ones at columns (j, g − s), s ∈ B-entry
    (A†⊗I) row (j,g) block i: A†[j,i] = conj(A[i,j]) ⇒ ones at h = g + s, s ∈ A[i,j].
"""
from __future__ import annotations

from itertools import product

from . import f2


class Ab:
    """Abelian group with index arithmetic tables."""

    def __init__(self, orders):
        self.orders = tuple(int(o) for o in orders)
        self.N = 1
        for o in self.orders:
            self.N *= o
        self.tuples = list(product(*[range(o) for o in self.orders]))
        self.index = {t: i for i, t in enumerate(self.tuples)}
        N = self.N
        self.add = [[0] * N for _ in range(N)]
        for i, a in enumerate(self.tuples):
            row = self.add[i]
            for j, b in enumerate(self.tuples):
                row[j] = self.index[tuple((x + y) % o for x, y, o in zip(a, b, self.orders))]
        self.neg = [self.index[tuple((-x) % o for x, o in zip(a, self.orders))] for a in self.tuples]

    def sub(self, i, j):
        return self.add[i][self.neg[j]]

    def units(self):
        """Automorphisms of G as index permutations (only for cyclic or
        products of cyclic groups: we use the full Aut for cyclic, and the
        coordinate-wise units ⋊ coordinate permutations for equal orders)."""
        import math
        perms = []
        if len(self.orders) == 1:
            n = self.orders[0]
            for u in range(1, n):
                if math.gcd(u, n) == 1:
                    perms.append([self.index[((u * t[0]) % n,)] for t in self.tuples])
        else:
            # coordinate-wise units (a subgroup of Aut(G)) — enough for de-duplication purposes
            unit_lists = [[u for u in range(1, o) if math.gcd(u, o) == 1] for o in self.orders]
            for us in product(*unit_lists):
                perms.append([self.index[tuple((u * x) % o for u, x, o in zip(us, t, self.orders))]
                              for t in self.tuples])
            if len(set(self.orders)) == 1:
                base = list(perms)
                r = len(self.orders)
                from itertools import permutations
                for sigma in permutations(range(r)):
                    if sigma == tuple(range(r)):
                        continue
                    for pm in base:
                        perms.append([pm[self.index[tuple(t[sigma[i]] for i in range(r))]] for t in self.tuples])
        return perms

    def poly_str(self, p, names="xyzw"):
        terms = []
        for g in sorted(p):
            t = self.tuples[g]
            parts = [f"{names[i]}^{e}" if e != 1 else names[i] for i, e in enumerate(t) if e]
            terms.append("*".join(parts) if parts else "1")
        return " + ".join(terms) if terms else "0"


def canon_poly(G: Ab, p, perms=None, shifts=True):
    """Canonical form of a poly under shifts (and optional automorphisms)."""
    best = None
    cand = [list(p)]
    if perms:
        cand = [[pm[g] for g in p] for pm in perms]
    for c in cand:
        for sh in (c if shifts else [0]):
            t = tuple(sorted(G.sub(g, sh) for g in c))
            if best is None or t < best:
                best = t
    return best


def _block_row(G: Ab, g: int, poly, sign: int, block: int, N: int) -> int:
    """Bits at columns block*N + h with h = g − s (sign=−1) or g + s (sign=+1)."""
    v = 0
    add = G.add; neg = G.neg
    for s in poly:
        h = add[g][neg[s]] if sign < 0 else add[g][s]
        v |= 1 << (block * N + h)
    return v


def lp_ab_rows(G: Ab, A: dict, a: int, b):
    """Return (rowsX, rowsZ) for LP(A, b). Columns: L(j,h)=j*N+h, R(i,h)=(a+i)*N+h."""
    N = G.N
    rowsX = []
    for i in range(a):
        for g in range(N):
            v = 0
            for j in range(a):
                p = A.get((i, j))
                if p:
                    v |= _block_row(G, g, p, -1, j, N)
            v |= _block_row(G, g, b, +1, a + i, N)     # (b* I)[g,h]=1 iff g−h ∈ −b iff h = g+s
            rowsX.append(v)
    rowsZ = []
    for j in range(a):
        for g in range(N):
            v = _block_row(G, g, b, -1, j, N)           # (b I) row g: h = g − s
            for i in range(a):
                p = A.get((i, j))
                if p:
                    v |= _block_row(G, g, p, +1, a + i, N)  # A†[j,i] = conj(A[i,j]): h = g + s
            rowsZ.append(v)
    return rowsX, rowsZ


def lp_AB_rows(G: Ab, A: dict, B: dict, a: int):
    """LP(A,B), both a×a.  Left qubits (j1, j2, h): col (j1*a + j2)*N + h;
    right qubits (i1, i2, h): col (a² + i1*a + i2)*N + h.
    X-checks (i1, j2, g): A⊗I: ones at L(j1, j2, g − s) for s∈A[i1,j1];
                          I⊗B†: (B†)[j2, i2] = conj(B[i2, j2]) ⇒ R(i1, i2, g + s), s∈B[i2,j2].
    Z-checks (j1, i2, g): I⊗B: L(j1, j2, g − s) for s∈B[i2, j2];
                          A†⊗I: A†[j1, i1] = conj(A[i1,j1]) ⇒ R(i1, i2, g + s), s∈A[i1,j1].
    """
    N = G.N
    rowsX = []
    for i1 in range(a):
        for j2 in range(a):
            for g in range(N):
                v = 0
                for j1 in range(a):
                    p = A.get((i1, j1))
                    if p:
                        v |= _block_row(G, g, p, -1, j1 * a + j2, N)
                for i2 in range(a):
                    p = B.get((i2, j2))
                    if p:
                        v |= _block_row(G, g, p, +1, a * a + i1 * a + i2, N)
                rowsX.append(v)
    rowsZ = []
    for j1 in range(a):
        for i2 in range(a):
            for g in range(N):
                v = 0
                for j2 in range(a):
                    p = B.get((i2, j2))
                    if p:
                        v |= _block_row(G, g, p, -1, j1 * a + j2, N)
                for i1 in range(a):
                    p = A.get((i1, j1))
                    if p:
                        v |= _block_row(G, g, p, +1, a * a + i1 * a + i2, N)
                rowsZ.append(v)
    return rowsX, rowsZ


def k_from_rows(rowsX, rowsZ, n, symmetric=False):
    rX = f2.rank(rowsX)
    if symmetric:
        return n - 2 * rX
    return n - rX - f2.rank(rowsZ)


def connected(rowsX, rowsZ, n):
    parent = list(range(n))
    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]; x = parent[x]
        return x
    for r in rowsX + rowsZ:
        bits = []
        while r:
            lb = r & -r; bits.append(lb.bit_length() - 1); r ^= lb
        r0 = find(bits[0])
        for x in bits[1:]:
            rx = find(x)
            if rx != r0:
                parent[rx] = r0
    return len({find(i) for i in range(n)}) == 1


def logical_reps(rows_same, rows_other, n):
    """Reps of ker(rows_other)/rs(rows_same)... careful: for X-distance we need
    Z-logical reps = ker(H_X)/rs(H_Z): call logical_reps(rowsZ, rowsX)."""
    ker = f2.nullspace(rows_other, n)
    acc = dict(f2.echelon(rows_same))
    reps = []
    for v in ker:
        r = f2.reduce(v, acc)
        if r:
            reps.append(v); acc[r & -r] = r
    return reps


def to_code(rowsX, rowsZ, n, meta=None):
    from .lp import CSSCode
    return CSSCode(f2.matrix_from_rows(rowsX, n), f2.matrix_from_rows(rowsZ, n), meta=meta or {})
