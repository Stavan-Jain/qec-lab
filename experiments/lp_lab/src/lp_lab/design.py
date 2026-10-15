"""Algebraic design of a=2 codes LP([[p,1],[1,q]], b) over an abelian group G.

Facts (verified numerically in tests):
  * k = 2·dim_F2  R/(b*, 1+pq) = 2·dim R/(b, 1+p*q*),  R = F2[G], * = conjugation g↦−g
    — maximal (= 2·dim ker b) iff pq ≡ 1 mod (b*).   [left kernel of H_X = [A | b* I] is
    Ann(b) ∩ Ann(1+p*q*); dim Ann(J) = N − dim J since F2[G] is Frobenius.]
  * d ≤ D_b := min weight of ker(b)\{0} whenever the X-logical (u, 0), u ∈ ker(b)⊕0,
    is nontrivial (generic).  So fom ≤ δ_b·D_b²/(2N) with δ_b = dim ker(b).
The ideal (b) is the F2-span of the translates of b, so everything is linear algebra
on bit-packed ints.
"""
from __future__ import annotations

import itertools
from dataclasses import dataclass

from . import f2
from .fastlp import Ab, canon_poly


def translate(G: Ab, v: int, g: int) -> int:
    """Multiply the ring element v (bitset over G) by the monomial g."""
    out = 0
    add = G.add[g]
    while v:
        lb = v & -v; h = lb.bit_length() - 1; v ^= lb
        out |= 1 << add[h]
    return out


def poly_int(support) -> int:
    v = 0
    for s in support: v |= 1 << s
    return v


def mul(G: Ab, u: int, v: int) -> int:
    out = 0
    while v:
        lb = v & -v; g = lb.bit_length() - 1; v ^= lb
        out ^= translate(G, u, g)
    return out


class Ideal:
    """The principal ideal (b) ⊂ F2[G] as an echelonised F2-subspace, with reduction."""

    def __init__(self, G: Ab, b: int):
        self.G = G
        rows = [translate(G, b, g) for g in range(G.N)]
        self.piv = f2.echelon(rows)          # {lowest bit: row}
        # fully reduce pivots against each other for a canonical reducer
        self.dim = len(self.piv)
        self.codim = G.N - self.dim

    def reduce(self, v: int) -> int:
        return f2.reduce(v, self.piv)

    def contains(self, v: int) -> bool:
        return self.reduce(v) == 0


def kernel_code(G: Ab, b: int) -> list[int]:
    """Basis of ker(b) = {u : b·u = 0} (as bitsets over G)."""
    rows = [translate(G, b, g) for g in range(G.N)]   # matrix of mult-by-b (columns = g) — rows here are b·x^g
    # u ∈ ker iff Σ_g u_g (b x^g) = 0 iff u ⊥ nothing... solve directly: nullspace of the matrix whose g-th column is b·x^g
    # Build matrix M with M[h, g] = [h ∈ supp(b x^g)]  → rows indexed by h.
    N = G.N
    M_rows = [0] * N
    for g, bg in enumerate(rows):
        v = bg
        while v:
            lb = v & -v; h = lb.bit_length() - 1; v ^= lb
            M_rows[h] |= 1 << g
    return f2.nullspace(M_rows, N)


def min_weight_brute(basis: list[int], cap: int = 22) -> int | None:
    """Exact min weight of the span of `basis` (2^dim enumeration) if dim ≤ cap."""
    k = len(basis)
    if k == 0:
        return None
    if k > cap:
        return None
    best = None
    # Gray-code enumeration
    v = 0
    for i in range(1, 1 << k):
        j = (i & -i).bit_length() - 1
        v ^= basis[j]
        w = v.bit_count()
        if best is None or w < best:
            best = w
    return best


@dataclass
class BInfo:
    b: tuple
    delta: int
    D: int | None
    fom_cap: float


def b_table(G: Ab, delta_range=(2, 16), cap=20) -> list[BInfo]:
    perms = G.units()
    seen = set()
    N = G.N
    for b1 in range(1, N):
        for b2 in range(b1 + 1, N):
            seen.add(canon_poly(G, (0, b1, b2), perms))
    out = []
    for b in sorted(seen):
        K = kernel_code(G, poly_int(b))
        delta = len(K)
        if delta < delta_range[0] or delta > delta_range[1]:
            continue
        D = min_weight_brute(K, cap=cap)
        cap_f = (delta * D * D / (2 * N)) if D else float("nan")
        out.append(BInfo(b, delta, D, cap_f))
    out.sort(key=lambda t: -(t.fom_cap if t.fom_cap == t.fom_cap else -1))
    return out


def max_k_designs(G: Ab, b: tuple, p_values=None, modulus: int | None = None):
    """All (p, q) with p = 1 + x^{p2}, q = x^{q1} + x^{q2}, pq ≡ 1 mod (b)  ⇒ k = 2·dim ker(b).
    Yields (p, q) tuples of group indices."""
    N = G.N
    b_conj = tuple(G.neg[g] for g in b)
    mod_int = modulus if modulus is not None else poly_int(b_conj)   # designs need pq ≡ 1 mod (b*) (or a given modulus)
    I = Ideal(G, mod_int)
    mono_res = {}
    for g in range(N):
        mono_res.setdefault(I.reduce(1 << g), []).append(g)
    # multiplication-by-p on R/(b): solve p·y ≡ 1 mod I via linear algebra on a basis of R/(b)
    # Simpler: brute-force y over the coset structure is expensive; instead compute the inverse by
    # solving the linear system  p·y + Σ c_g (b x^g) = 1  over F2 in the unknowns (y, c).
    # We do this with a small Gaussian elimination on N + N unknowns → fine for N ≤ ~150.
    b_int = mod_int
    cols_b = [translate(G, b_int, g) for g in range(N)]
    for p2 in (p_values if p_values is not None else range(1, N)):
        p = (0, p2); p_int = poly_int(p)
        cols_p = [translate(G, p_int, g) for g in range(N)]   # p·x^g
        # Unknown vector z = (y_0..y_{N-1}, c_0..c_{N-1}); equation Σ y_g (p x^g) + Σ c_g (b x^g) = 1.
        # Build augmented system: for each group element h (equation row), bit over 2N unknowns + rhs.
        eq = [0] * N
        for g in range(N):
            v = cols_p[g]
            while v:
                lb = v & -v; h = lb.bit_length() - 1; v ^= lb
                eq[h] |= 1 << g
            v = cols_b[g]
            while v:
                lb = v & -v; h = lb.bit_length() - 1; v ^= lb
                eq[h] |= 1 << (N + g)
        # rhs = 1 at h = 0
        rows = [(eq[h] << 1) | (1 if h == 0 else 0) for h in range(N)]   # bit0 = rhs
        # solve: echelon on unknown bits (bits 1..2N), track rhs in bit 0
        piv = {}
        consistent = True
        for r in rows:
            rr = r
            while (rr >> 1):
                lb = (rr >> 1) & -(rr >> 1)
                bb = piv.get(lb)
                if bb is None:
                    piv[lb] = rr; break
                rr ^= bb
            else:
                if rr & 1:
                    consistent = False; break
        if not consistent:
            continue   # p not a unit mod (b)
        # back-substitute a particular solution: set free vars 0
        sol = 0
        for lb in sorted(piv.keys(), reverse=True):
            rr = piv[lb]
            # value of this pivot var = rhs xor Σ(other set vars in rr)
            val = rr & 1
            others = (rr >> 1) ^ lb
            v = others
            while v:
                l2 = v & -v; v ^= l2
                if sol & l2: val ^= 1
            if val: sol |= lb
        y = sol & ((1 << N) - 1)          # p^{-1} mod (b) representative
        target = I.reduce(y)
        # find q = x^{q1} + x^{q2} with reduce(q) == target
        for q1 in range(N):
            need = I.reduce(target ^ (1 << q1))
            for q2 in mono_res.get(need, []):
                if q2 > q1:
                    yield p, (q1, q2)
