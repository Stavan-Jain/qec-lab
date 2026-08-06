"""A28 — spectral layer for odd abelian groups (semisimple Fourier side).

For odd G = Z_l x Z_m, exponent e | 15 here, characters take values in
GF(16) = F2[w]/(w^4+w+1); mu_15 = GF(16)^x.  chi_(s,t)(x^a y^b) =
g^((15/l)*s*a + (15/m)*t*b mod 15) for a generator g of GF(16)^x.

Facts used downstream (proofs in the A28 note):
  - u ~ uhat is a ring iso onto pointwise functions (semisimple CRT);
    Z(fA) = Z(f) u Z(A) EXACTLY.
  - Z(u) is Galois-closed (2-cyclotomic: (s,t) -> (2s,2t)) for F2 words.
  - dim_F2 {u : Z(u) >= O} = |G| - |O| for Galois-closed O.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from a28_lsc_lib import GroupAlg

# ---- GF(16), poly w^4 + w + 1; ints 0..15, generator g = 2 (=w) ----
_EXP = [1] * 30
_LOG = [0] * 16
_v = 1
for _i in range(15):
    _EXP[_i] = _v
    _LOG[_v] = _i
    _v <<= 1
    if _v & 16:
        _v ^= 0b10011  # w^4 = w + 1
for _i in range(15, 30):
    _EXP[_i] = _EXP[_i - 15]


def gf_mul(a, b):
    if a == 0 or b == 0:
        return 0
    return _EXP[_LOG[a] + _LOG[b]]


class Spectral:
    """Character table + Fourier for F2[Z_l x Z_m], l*m odd, exponent | 15."""

    def __init__(self, G: GroupAlg):
        assert G.l % 2 == 1 and G.m % 2 == 1
        assert 15 % G.l == 0 and 15 % G.m == 0, "exponent must divide 15"
        self.G = G
        self.sl = 15 // G.l
        self.sm = 15 // G.m
        N = G.N
        # chi_table[chi_idx][g_idx] = value in GF(16)
        self.chi = [[0] * N for _ in range(N)]
        for s in range(G.l):
            for t in range(G.m):
                ci = s * G.m + t
                row = self.chi[ci]
                for a in range(G.l):
                    for b in range(G.m):
                        row[a * G.m + b] = _EXP[(self.sl * s * a + self.sm * t * b) % 15]

    def fourier(self, u: int) -> list[int]:
        """uhat over all N characters (index (s,t) -> s*m+t)."""
        G = self.G
        sup = []
        uu = u
        i = 0
        while uu:
            if uu & 1:
                sup.append(i)
            uu >>= 1
            i += 1
        out = [0] * G.N
        for ci in range(G.N):
            row = self.chi[ci]
            acc = 0
            for gi in sup:
                acc ^= row[gi]
            out[ci] = acc
        return out

    def zero_set(self, u: int) -> frozenset[int]:
        return frozenset(ci for ci, v in enumerate(self.fourier(u)) if v == 0)

    # ---- Galois (2-cyclotomic) orbits on the character grid ----
    def galois_orbits(self) -> list[frozenset[int]]:
        G = self.G
        seen = set()
        orbits = []
        for s in range(G.l):
            for t in range(G.m):
                ci = s * G.m + t
                if ci in seen:
                    continue
                orb = set()
                cs, ct = s, t
                while (cs * G.m + ct) not in orb:
                    orb.add(cs * G.m + ct)
                    cs, ct = (2 * cs) % G.l, (2 * ct) % G.m
                orbits.append(frozenset(orb))
                seen |= orb
        return orbits

    def translate_chi(self, ci: int, mu: int) -> int:
        """Index of chi_ci * chi_mu (additive on (s,t))."""
        G = self.G
        s1, t1 = divmod(ci, G.m)
        s2, t2 = divmod(mu, G.m)
        return ((s1 + s2) % G.l) * G.m + ((t1 + t2) % G.m)


def solve_boundary(inst, b_pair):
    """Find one f with del f = b (assumes membership); returns int mask."""
    from a28_lsc_lib import rref
    G = inst.G
    N = G.N
    target = inst.pack(b_pair)
    # rows: for each generator monomial g, (del g | e_g-marker) — augment with
    # identity to read off the combination.
    aug = []
    for gi, (a, b) in enumerate(G.elements()):
        row = inst.pack(inst.boundary(G.monomial(a, b)))
        aug.append(row | (1 << (2 * N + gi)))
    basis, piv = rref(aug, 3 * N, col_order=range(2 * N))  # pivot only in code cols
    v = target
    used = 0
    for brow, c in zip(basis, piv):
        if v & (1 << c):
            v ^= brow & ((1 << 2 * N) - 1)
            used ^= brow >> (2 * N)
    assert v == 0, "not a boundary"
    f = used
    assert inst.boundary(f) == b_pair
    return f
