"""A28 shared library — the Light Stabilizer Classification (LSC) problem.

Problem. For a BB code over an abelian group G = Z_l(x) x Z_m(y) with
polynomials A, B in F2[G], the X-stabilizer space is the "boundary code"

    C(A,B) = { del f := (f*A, f*B) : f in F2[G] }  subset  F2[G]^2 ,

a 1-generator quasi-abelian code of index 2 (dim = |G| - k/2).  The LSC
problem with threshold W asks for the complete list of translation classes
of nonzero boundaries b with |b| = |f*A| + |f*B| <= W.  Instances:

  - gross base   (Z6 x Z6,  W = 11): the (CLASS) input of the doubling
    theorem — hexagons + D-pairs (A4 §6.3, hand-proved).
  - f2a6f17e     (Z5 x Z15, W = 14): 113 classes, SAT 9.6 h + final UNSAT
    34,554 s (data/a17/f2a6_light_classes.jsonl); analytic re-derivation
    via CRT fibering (A22).
  - A17 docket   (Z15 x Z6, W = 18): the two [[180,4,10]] UNKNOWN cells —
    census extrapolated "days-plus" for the SAT lane (A27 §2), open.

Conventions (locked against the f2a6 census file, see a28_v_foundations):
  - group element (a, b) with a in Z_l (x-exponent), b in Z_m (y-exponent);
    linear index i = a*m + b.
  - F2[G] elements are python int bitmasks over N = l*m bits.
  - boundary del f = (f*A, f*B); block 0 = A-block (u), block 1 = B-block (v).
  - pair (u, v) packs into a 2N-bit int as u | (v << N).
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path

DATA = Path(__file__).resolve().parent.parent / "data"


# ----------------------------------------------------------------------
# group algebra F2[Z_l x Z_m] as int bitmasks
# ----------------------------------------------------------------------

@dataclass(frozen=True)
class GroupAlg:
    l: int  # order of x
    m: int  # order of y

    @property
    def N(self) -> int:
        return self.l * self.m

    def idx(self, a: int, b: int) -> int:
        return (a % self.l) * self.m + (b % self.m)

    def coords(self, i: int) -> tuple[int, int]:
        return divmod(i, self.m)

    def monomial(self, a: int, b: int) -> int:
        return 1 << self.idx(a, b)

    def from_support(self, pts) -> int:
        v = 0
        for a, b in pts:
            v ^= self.monomial(a, b)
        return v

    def support(self, u: int) -> list[tuple[int, int]]:
        out = []
        i = 0
        while u:
            if u & 1:
                out.append(self.coords(i))
            u >>= 1
            i += 1
        return out

    def translate(self, u: int, a: int, b: int) -> int:
        """u * x^a y^b (translation by group element (a,b))."""
        out = 0
        for (p, q) in self.support(u):
            out |= self.monomial(p + a, q + b)
        return out

    def mul(self, f: int, g: int) -> int:
        """Convolution product in F2[G]."""
        out = 0
        for (a, b) in self.support(f):
            out ^= self.translate(g, a, b)
        return out

    def weight(self, u: int) -> int:
        return bin(u).count("1")

    def elements(self):
        for a in range(self.l):
            for b in range(self.m):
                yield (a, b)

    def translation_perm(self, a: int, b: int) -> list[int]:
        """Permutation pi with (x^a y^b * u)_bit[pi[i]] = u_bit[i]."""
        return [self.idx(p + a, q + b) for p in range(self.l) for q in range(self.m)]


# ----------------------------------------------------------------------
# the code registry
# ----------------------------------------------------------------------

@dataclass
class LSCInstance:
    name: str
    G: GroupAlg
    A_supp: list
    B_supp: list
    expect_k: int
    d: int            # code distance (context; not recomputed here)
    W: int            # census threshold
    note: str = ""
    A: int = field(init=False)
    B: int = field(init=False)

    def __post_init__(self):
        self.A = self.G.from_support(self.A_supp)
        self.B = self.G.from_support(self.B_supp)

    def boundary(self, f: int) -> tuple[int, int]:
        return self.G.mul(f, self.A), self.G.mul(f, self.B)

    def pair_weight(self, uv: tuple[int, int]) -> int:
        return self.G.weight(uv[0]) + self.G.weight(uv[1])

    def pack(self, uv: tuple[int, int]) -> int:
        return uv[0] | (uv[1] << self.G.N)

    def unpack(self, w: int) -> tuple[int, int]:
        return w & ((1 << self.G.N) - 1), w >> self.G.N

    # --- generator matrix of C(A,B): rows (gA | gB), g in G ---
    def generator_rows(self) -> list[int]:
        return [self.pack(self.boundary(self.G.monomial(a, b)))
                for (a, b) in self.G.elements()]

    def canonical(self, uv: tuple[int, int]) -> int:
        """Min over translations of the packed pair — the class label."""
        u, v = uv
        best = None
        for (a, b) in self.G.elements():
            t = self.G.translate(u, a, b) | (self.G.translate(v, a, b) << self.G.N)
            if best is None or t < best:
                best = t
        return best


REGISTRY: dict[str, LSCInstance] = {}


def _reg(inst: LSCInstance):
    REGISTRY[inst.name] = inst
    return inst


# gross base: [[72,12,6]] over Z6(x) x Z6(y); (CLASS) threshold 11
_reg(LSCInstance(
    name="grossbase", G=GroupAlg(6, 6),
    A_supp=[(3, 0), (0, 1), (0, 2)], B_supp=[(0, 3), (1, 0), (2, 0)],
    expect_k=12, d=6, W=11,
    note="A4 §6.3 (CLASS): 36 hexagons (w6) + 216 D-pairs (w10), nothing else <= 11",
))

# f2a6f17e: [[150,8,8]] over Z5(x) x Z15(y); census threshold 14, 113 classes
_reg(LSCInstance(
    name="f2a6", G=GroupAlg(5, 15),
    A_supp=[(0, 0), (0, 1), (1, 0)], B_supp=[(1, 6), (1, 10), (2, 12)],
    expect_k=8, d=8, W=14,
    note="A22 ground truth: 113 classes (6:1, 10:7, 12:36, 14:69)",
))

# A17 docket UNKNOWN cells: [[180,4,10]] over Z15(x) x Z6(y); threshold 18
_reg(LSCInstance(
    name="docket37", G=GroupAlg(15, 6),
    A_supp=[(0, 0), (0, 1), (1, 0)], B_supp=[(0, 4), (1, 0), (11, 2)],
    expect_k=4, d=10, W=18,
    note="37a70e02e003d1de: census open (SAT lane days-plus)",
))
_reg(LSCInstance(
    name="docket5e", G=GroupAlg(15, 6),
    A_supp=[(0, 0), (0, 1), (1, 0)], B_supp=[(0, 4), (8, 2), (13, 0)],
    expect_k=4, d=10, W=18,
    note="5e50a9765a02eb70: census open (SAT lane days-plus)",
))


# ----------------------------------------------------------------------
# F2 linear algebra on int-bitmask rows
# ----------------------------------------------------------------------

def rref(rows: list[int], ncols: int, col_order=None):
    """Row-reduce; returns (basis_rows, pivot_cols) with pivots chosen in
    col_order (default 0..ncols-1). basis_rows[i] has leading 1 at
    pivot_cols[i] and 0 at every other pivot column."""
    if col_order is None:
        col_order = range(ncols)
    basis: list[int] = []
    pivots: list[int] = []
    work = [r for r in rows if r]
    for c in col_order:
        bit = 1 << c
        # find a row with this bit among remaining work rows
        hit = None
        for i, r in enumerate(work):
            if r & bit:
                hit = i
                break
        if hit is None:
            continue
        piv = work.pop(hit)
        work = [(r ^ piv) if (r & bit) else r for r in work]
        work = [r for r in work if r]
        basis = [(b ^ piv) if (b & bit) else b for b in basis]
        basis.append(piv)
        pivots.append(c)
        if not work:
            break
    return basis, pivots


def rank(rows: list[int], ncols: int) -> int:
    return len(rref(rows, ncols)[0])


def in_span(rows_rref: list[int], pivots: list[int], v: int) -> bool:
    for b, c in zip(rows_rref, pivots):
        if v & (1 << c):
            v ^= b
    return v == 0


def load_f2a6_census() -> list[dict]:
    rows = []
    with open(DATA / "a17" / "f2a6_light_classes.jsonl") as fh:
        for line in fh:
            r = json.loads(line)
            if "b_weight" in r:
                rows.append(r)
    return rows
