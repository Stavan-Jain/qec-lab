"""A27 — A22-fibering feasibility probe for the Z15xZ6 [[180,4,10]] docket cells.

Provenance: written 2026-08-06 during the teaching-PDF session (see
notes/A27_safe_floor_generality.md §3). Pure Python, deliberately
independent of bb_lab primitives (independent implementation = its own
cross-check); no solver, no Lean. Run:

    python3 scripts/a27_fibering_feasibility_z15z6.py     (~seconds)

Checks, for 37a70e02e003d1de (A = 1+y+x, B = y^4+x+x^11 y^2) and
5e50a9765a02eb70 (A = 1+y+x, B = y^4+x^8 y^2+x^13) on G = Z15(x) x Z6(y)
(the three A17-docket UNKNOWN safe-floor-20 cells; polynomials from
data/a17/docket_decision.jsonl):

  P1  code params: k = 4, dim K = dim ker d2 = 2
  P2  fiber split z := x^3 (order 5): is K purely epsilon-sector
      (N*kappa == kappa for the fiber-averaging idempotent N)?
      => delta-kernel dim = dim K - eps-part (expect 0)
  P3  the shared epsilon-quotient code over Z3(w) x Z6(y):
      params (n=36, kbar, dbar by exhaustive cycle enumeration),
      and whether both codes share it verbatim
  P4  spectral checks (semisimplified spectrum): character zeros of
      A_eps over the 9 odd characters (V3-analog: unit iff 0 zeros);
      unit test for Atilde over GF(16); common ss-zeros of
      (Atilde, Btilde) per code (V4-analog)
  P5  sweep arithmetic: sum_{i<=9} C(18,i)
"""
from itertools import product
from math import comb

# ---------- generic F2 linear algebra on int-bitmask rows ----------
def rank_and_kernel(rows, ncols):
    """rows: list of ints (bitmask over ncols). Returns (rank, kernel_basis)
    where kernel basis vectors are ints over ncols (solutions x of Mx=0)."""
    rows = [r for r in rows]
    pivots = {}  # col -> row index in reduced list
    red = []
    for r in rows:
        cur = r
        for c, ri in pivots.items():
            if (cur >> c) & 1:
                cur ^= red[ri]
        if cur:
            c = (cur & -cur).bit_length() - 1
            for i, rr in enumerate(red):
                if (rr >> c) & 1:
                    red[i] = rr ^ cur
            pivots[c] = len(red)
            red.append(cur)
    rank = len(red)
    piv_cols = set(pivots.keys())
    kernel = []
    for free in range(ncols):
        if free in piv_cols:
            continue
        v = 1 << free
        for c, ri in pivots.items():
            if (red[ri] >> free) & 1:
                v |= 1 << c
        kernel.append(v)
    return rank, kernel

def make_rref(rows):
    pivots = {}
    red = []
    for r in rows:
        cur = r
        for c, ri in pivots.items():
            if (cur >> c) & 1:
                cur ^= red[ri]
        if cur:
            c = (cur & -cur).bit_length() - 1
            for i, rr in enumerate(red):
                if (rr >> c) & 1:
                    red[i] = rr ^ cur
            pivots[c] = len(red)
            red.append(cur)
    return red, pivots

def in_span(rref_rows_pivots, v):
    red, pivots = rref_rows_pivots
    cur = v
    for c, ri in pivots.items():
        if (cur >> c) & 1:
            cur ^= red[ri]
    return cur == 0

# ---------- group algebra machinery ----------
class GA:
    """F2 group algebra of Z_l x Z_m; elements = frozenset of (a,b)."""
    def __init__(self, l, m):
        self.l, self.m = l, m
        self.elems = [(a, b) for a in range(l) for b in range(m)]
        self.idx = {g: i for i, g in enumerate(self.elems)}
        self.n = l * m

    def mul(self, P, Q):
        acc = set()
        for (a1, b1) in P:
            for (a2, b2) in Q:
                g = ((a1 + a2) % self.l, (b1 + b2) % self.m)
                if g in acc: acc.remove(g)
                else: acc.add(g)
        return frozenset(acc)

    def mulmat_rows(self, P):
        """Rows of M_P (coeff of P*f at g), as ints over columns h."""
        rows = []
        for g in self.elems:
            r = 0
            for (pa, pb) in P:
                h = ((g[0] - pa) % self.l, (g[1] - pb) % self.m)
                r |= 1 << self.idx[h]
            rows.append(r)
        return rows

    def vec_of(self, P):
        v = 0
        for g in P:
            v |= 1 << self.idx[g]
        return v

def code_data(G, A, B):
    MA = G.mulmat_rows(A); MB = G.mulmat_rows(B)
    n = 2 * G.n
    # C = [M_A | M_B] (syndrome map); cycles = ker C
    Crows = [MA[i] | (MB[i] << G.n) for i in range(G.n)]
    rC, kerC = rank_and_kernel(Crows, n)
    # S: f -> (B f, A f); kernel of S = K = Ann(A) ∩ Ann(B)
    Srows = MB + MA
    rS, kerS = rank_and_kernel(Srows, G.n)
    k = n - 2 * rC
    return dict(MA=MA, MB=MB, Crows=Crows, rC=rC, kerC=kerC,
                Srows=Srows, rS=rS, K=kerS, k=k, n=n)

# ---------- GF(16), modulus x^4 + x + 1, x primitive ----------
def gf16_mul(a, b):
    r = 0
    for i in range(4):
        if (b >> i) & 1:
            r ^= a << i
    for i in range(7, 3, -1):
        if (r >> i) & 1:
            r ^= (0b10011 << (i - 4))
    return r
POW = [1]
for _ in range(1, 15):
    POW.append(gf16_mul(POW[-1], 0b0010))
def gpow(k): return POW[k % 15]
ZETA = gpow(3)   # order 5
OMEGA = gpow(5)  # order 3

# ---------- the two codes ----------
G = GA(15, 6)
A = frozenset({(0,0), (0,1), (1,0)})                    # 1 + y + x
codes = {
    '37a70e02': frozenset({(0,4), (1,0), (11,2)}),      # y^4 + x + x^11 y^2
    '5e50a9':   frozenset({(0,4), (8,2), (13,0)}),      # y^4 + x^8 y^2 + x^13
}

# fiber z = x^3 (order 5); averaging idempotent N = sum_i z^i
N = frozenset({(3*i % 15, 0) for i in range(5)})

# epsilon map on exponents: x^a y^b -> w^{2a mod 3} y^b  (Z3(w) x Z6(y));
# uses x = z^2 w^2 with z = x^3, w = x^5 (so x^k -> w^{2k mod 3} under z->1)
def eps_poly(P):
    acc = set()
    for (a, b) in P:
        g = ((2*a) % 3, b)
        if g in acc: acc.remove(g)
        else: acc.add(g)
    return frozenset(acc)

Gq = GA(3, 6)

def quotient_distance(Gq, Aq, Bq):
    d = code_data(Gq, Aq, Bq)
    stab_rows = []
    for h in Gq.elems:
        dg = frozenset({h})
        Bz = Gq.mul(Bq, dg); Az = Gq.mul(Aq, dg)
        stab_rows.append(Gq.vec_of(Bz) | (Gq.vec_of(Az) << Gq.n))
    rref = make_rref(stab_rows)
    basis = d['kerC']
    dimc = len(basis)
    best = None
    for mask in range(1, 1 << dimc):
        v = 0
        mm = mask
        while mm:
            i = (mm & -mm).bit_length() - 1
            v ^= basis[i]
            mm &= mm - 1
        if not in_span(rref, v):
            w = bin(v).count('1')
            if best is None or w < best:
                best = w
    return d, best

def delta_poly_val(P, i, j):
    """delta-part value at (chi_w = omega^i, chi_ybar = omega^j, z -> zeta),
    semisimplified (y^3 -> 1). x^a -> zeta^{2a} * omega^{i*(2a mod 3)};
    y^b -> omega^{j*b}."""
    v = 0
    for (a, b) in P:
        term = gpow(3*(2*a) % 15)                        # zeta^{2a}
        term = gf16_mul(term, gpow((5*((2*a % 3) * i)) % 15))
        term = gf16_mul(term, gpow((5*(b*j)) % 15))
        v ^= term
    return v

if __name__ == '__main__':
    print("=" * 70)
    Aq = eps_poly(A)
    print(f"shared A_eps = {sorted(Aq)}  (over Z3(w) x Z6(y): 1 + y + w^2)")
    for name, B in codes.items():
        print("=" * 70)
        print(f"[{name}]  B = {sorted(B)}")
        d = code_data(G, A, B)
        dimK = len(d['K'])
        print(f"P1: n={d['n']}  k={d['k']}  dim K = {dimK}")
        pure = 0
        for kv in d['K']:
            P = frozenset(G.elems[i] for i in range(G.n) if (kv >> i) & 1)
            if G.mul(N, P) == P:
                pure += 1
        print(f"P2: K basis vectors with N*kappa == kappa: {pure}/{dimK} "
              f"=> delta-kernel dim {'0 (all eps-sector)' if pure == dimK else '>0 !'}")
        print(f"P3: B_eps = {sorted(eps_poly(B))}")

    Bq37 = eps_poly(codes['37a70e02']); Bq5e = eps_poly(codes['5e50a9'])
    print("=" * 70)
    print(f"P3: quotient codes identical for both: {Bq37 == Bq5e}")
    dq, dbar = quotient_distance(Gq, Aq, Bq37)
    print(f"P3: eps-quotient over Z3xZ6: n={dq['n']} k={dq['k']} "
          f"dim ker d2 = {len(dq['K'])}  d_bar = {dbar}")

    print("=" * 70)
    zeros_Aeps = []
    for i, j in product(range(3), repeat=2):
        val = 1 ^ gpow(5*j) ^ gpow(5*(2*i) % 15)   # 1 + omega^j + omega^{2i}
        if val == 0:
            zeros_Aeps.append((i, j))
    print(f"P4: A_eps ss-character zeros (of 9): {len(zeros_Aeps)} at {zeros_Aeps}"
          f"  => A_eps {'NOT a unit' if zeros_Aeps else 'unit'}")
    zA = [(i,j) for i,j in product(range(3),repeat=2) if delta_poly_val(A,i,j)==0]
    print(f"P4: Atilde ss-zeros (of 9): {len(zA)}  => Atilde "
          f"{'unit (transfer T exists everywhere)' if not zA else 'NOT unit'}")
    for name, B in codes.items():
        zB = [(i,j) for i,j in product(range(3),repeat=2) if delta_poly_val(B,i,j)==0]
        common = [p for p in zA if p in zB]
        print(f"P4: [{name}] Btilde ss-zeros: {len(zB)}; common with Atilde: {len(common)}")

    total = sum(comb(18, i) for i in range(10))
    print("=" * 70)
    print(f"P5: classification budget 2d-2 = 18, active-site cost >= 2 "
          f"=> <= 9 active sites of 18; sweep size sum C(18,<=9) = {total}")
