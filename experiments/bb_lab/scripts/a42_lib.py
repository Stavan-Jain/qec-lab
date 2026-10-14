#!/usr/bin/env python3
"""A42 shared algebra library — the spectral/module lane on the
tour-de-gross column (sibling lane to A40, same branch).

Exact machinery over F_2:

  * F_2[x] arithmetic on ints (bit i = coeff of x^i), factorization,
    Sylvester resultants with F_2[x] entries.
  * F_{2^d} towers (poly-basis ints) and quadratic Artin–Schreier
    extensions F_{2^d}[s]/(s^2+s+gamma) for fibers that leave F_q.
  * The variety V(A~, B^) of the fixed tour-de-gross pair
        A = 1 + y + x^3 y^-1,  B = 1 + x + x^-1 y^-3
    cleared to polynomial form by unit (monomial) multiplication:
        A~ = y*A     = y + y^2 + x^3          (monic in y, deg 2)
        B^ = x*y^3*B = 1 + (x + x^2) y^3     (deg_y 3, lc = x + x^2)
    Unit multiplication changes neither the ideal (A,B) in any group
    algebra nor the torus zero set.
  * The local multiplicity engine: for a point P = (alpha, beta) with
    residue field k(P) and a 2-group truncation, the dimension over
    k(P) of
        k(P)[u,v] / (u^{2^a1}, v^{2^a2}, A_loc, B_loc)
    where x |-> alpha*(1+u)^{ex_u}(1+v)^{ex_v} etc.  With
    (ex_x, ex_y) = ((1,0),(0,1)) and power truncations (N,N) this is
    the plane local intersection multiplicity (once stabilized in N).
  * The frame spectral k-formula: for the BB code on G = Z^2/L,
    L = <(l,0),(d,p)>,
        k = 2 * sum over Frobenius orbits O of V(A,B) with
            alpha^l = 1 and alpha^d beta^p = 1 of
            |O| * dim_{k(P)} k(P)[G_2] / (A_loc, B_loc)
    which is exact: F_2[G] = prod_O F_{q_O}[G_2] (odd/2-part CRT), a
    local factor survives iff the character kills both polynomials to
    a non-unit, and k = 2*dim_{F_2} F_2[G]/(A,B) because F_2[G] is
    self-injective (annihilator duality), so rank H_X = rank H_Z =
    |G| - dim R/(A,B).

Everything is falsify-first: the s0 battery cross-checks the formula
against TowerCode on hundreds of frames before any claim consumes it.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from math import gcd
from typing import Optional

# ---------------------------------------------------------------------------
# The fixed pair (polynomial form; exponents (x-exp, y-exp), all >= 0)
# ---------------------------------------------------------------------------
AT_SUPP = ((0, 1), (0, 2), (3, 0))   # A~ = y + y^2 + x^3
BH_SUPP = ((0, 0), (1, 3), (2, 3))   # B^ = 1 + x y^3 + x^2 y^3
# Laurent originals (for frame transports where only the ideal matters,
# A~/B^ are unit multiples — use these fixed polynomial supports).

# ---------------------------------------------------------------------------
# F_2[x] on ints
# ---------------------------------------------------------------------------


def pdeg(a: int) -> int:
    return a.bit_length() - 1


def pmul(a: int, b: int) -> int:
    r = 0
    while b:
        lb = b & -b
        r ^= a << (lb.bit_length() - 1)
        b ^= lb
    return r


def pdivmod(a: int, b: int) -> tuple[int, int]:
    assert b, "division by zero poly"
    db = pdeg(b)
    q = 0
    while a and pdeg(a) >= db:
        s = pdeg(a) - db
        q ^= 1 << s
        a ^= b << s
    return q, a


def pmod(a: int, b: int) -> int:
    return pdivmod(a, b)[1]


def pgcd(a: int, b: int) -> int:
    while b:
        a, b = b, pmod(a, b)
    return a


def is_irreducible(f: int) -> bool:
    """Rabin test over F_2 (exact, deterministic)."""
    d = pdeg(f)
    if d <= 0:
        return False
    if d == 1:
        return True
    # x^(2^d) == x mod f, and x^(2^(d/q)) != x for prime q | d
    x = 2

    def frob_iter(e: int) -> int:
        # x^(2^e) mod f by repeated squaring of the element
        t = x
        for _ in range(e):
            t = pmod(pmul(t, t), f)
        return t

    if frob_iter(d) != x:
        return False
    dd = d
    primes = set()
    q = 2
    while q * q <= dd:
        while dd % q == 0:
            primes.add(q)
            dd //= q
        q += 1
    if dd > 1:
        primes.add(dd)
    for q in primes:
        if pgcd(frob_iter(d // q) ^ x, f) != 1:
            return False
    return True


def factorize_f2(f: int) -> dict[int, int]:
    """Full factorization of f in F_2[x] into {irreducible: exponent}.
    Trial division by irreducibles in degree order (fine for deg <= 30)."""
    assert f != 0
    out: dict[int, int] = {}
    d = 1
    while pdeg(f) > 0:
        if is_irreducible(f):
            out[f] = out.get(f, 0) + 1
            break
        found = False
        while d <= pdeg(f) // 2:
            for g in range(1 << d, 1 << (d + 1)):
                if not is_irreducible(g):
                    continue
                while pmod(f, g) == 0:
                    out[g] = out.get(g, 0) + 1
                    f = pdivmod(f, g)[0]
                    found = True
            if found:
                break
            d += 1
        if not found:
            # remaining f irreducible of degree > half
            if pdeg(f) > 0:
                out[f] = out.get(f, 0) + 1
            break
    return out


def poly_str(f: int) -> str:
    if f == 0:
        return "0"
    terms = []
    i = 0
    while f:
        if f & 1:
            terms.append("1" if i == 0 else ("x" if i == 1 else f"x^{i}"))
        f >>= 1
        i += 1
    return "+".join(terms)


def factorint(n: int) -> dict[int, int]:
    out: dict[int, int] = {}
    d = 2
    while d * d <= n:
        while n % d == 0:
            out[d] = out.get(d, 0) + 1
            n //= d
        d += 1
    if n > 1:
        out[n] = out.get(n, 0) + 1
    return out


# ---------------------------------------------------------------------------
# Fields
# ---------------------------------------------------------------------------


class F2k:
    """F_{2^d} = F_2[t]/(mod), elements = ints (poly basis)."""

    def __init__(self, mod: int):
        assert is_irreducible(mod), poly_str(mod)
        self.mod = mod
        self.d = pdeg(mod)
        self.q = 1 << self.d
        self.zero = 0
        self.one = 1
        self.t = 2  # the generator/root of mod

    def add(self, a: int, b: int) -> int:
        return a ^ b

    def mul(self, a: int, b: int) -> int:
        return pmod(pmul(a, b), self.mod)

    def pow(self, a: int, e: int) -> int:
        if e < 0:
            return self.pow(self.inv(a), -e)
        r, b = 1, a
        while e:
            if e & 1:
                r = self.mul(r, b)
            b = self.mul(b, b)
            e >>= 1
        return r

    def inv(self, a: int) -> int:
        assert a != 0
        # extended Euclid in F_2[x]
        r0, r1 = self.mod, a
        s0, s1 = 0, 1
        while r1:
            q, r = pdivmod(r0, r1)
            r0, r1 = r1, r
            s0, s1 = s1, s0 ^ pmul(q, s1)
        assert r0 == 1
        return pmod(s0, self.mod)

    def trace(self, a: int) -> int:
        t, s = a, a
        for _ in range(self.d - 1):
            t = self.mul(t, t)
            s ^= t
        assert s in (0, 1)
        return s

    def frob(self, a: int) -> int:
        return self.mul(a, a)

    def embed_f2(self, c: int) -> int:
        return c & 1

    def eq(self, a: int, b: int) -> bool:
        return a == b


class F2kExt:
    """Quadratic Artin–Schreier extension K = F_q[s]/(s^2+s+gamma),
    elements = (a, b) meaning a + b*s; requires Tr_{F_q/F_2}(gamma)=1."""

    def __init__(self, base: F2k, gamma: int):
        assert base.trace(gamma) == 1, "s^2+s+gamma must be irreducible"
        self.base = base
        self.gamma = gamma
        self.d = 2 * base.d
        self.q = 1 << self.d
        self.zero = (0, 0)
        self.one = (1, 0)
        self.s = (0, 1)

    def add(self, x, y):
        return (x[0] ^ y[0], x[1] ^ y[1])

    def mul(self, x, y):
        a, b = x
        c, d = y
        F = self.base
        ac = F.mul(a, c)
        bd = F.mul(b, d)
        ad_bc = F.mul(a, d) ^ F.mul(b, c)
        return (ac ^ F.mul(bd, self.gamma), ad_bc ^ bd)

    def pow(self, x, e: int):
        if e < 0:
            return self.pow(self.inv(x), -e)
        r, b = self.one, x
        while e:
            if e & 1:
                r = self.mul(r, b)
            b = self.mul(b, b)
            e >>= 1
        return r

    def inv(self, x):
        a, b = x
        F = self.base
        # (a+bs)(a+b+bs) = a^2+ab+gamma b^2  (in F_q)
        n = F.mul(a, a) ^ F.mul(a, b) ^ F.mul(self.gamma, F.mul(b, b))
        ninv = F.inv(n)
        return (F.mul(a ^ b, ninv), F.mul(b, ninv))

    def embed(self, a: int):
        """F_q -> K."""
        return (a, 0)

    def embed_f2(self, c: int):
        return (c & 1, 0)

    def eq(self, x, y) -> bool:
        return x == y


def artin_schreier_solve(F: F2k, gamma: int) -> Optional[int]:
    """One solution of b^2 + b = gamma in F_q, or None (Tr = 1)."""
    d = F.d
    # columns of the F_2-linear map b |-> b^2+b in poly basis
    cols = [F.mul(1 << i, 1 << i) ^ (1 << i) for i in range(d)]
    # solve sum b_i cols[i] = gamma by Gaussian elimination on bits
    rows = []  # (vector over F_2 of length d as int over basis, rhs-mask)
    # build augmented system: unknown x (d bits), M x = gamma where
    # M has columns cols. Work with rows = bit positions.
    M = [[(cols[j] >> i) & 1 for j in range(d)] for i in range(d)]
    rhs = [(gamma >> i) & 1 for i in range(d)]
    piv = []
    r = 0
    for c in range(d):
        pr = None
        for i in range(r, d):
            if M[i][c]:
                pr = i
                break
        if pr is None:
            continue
        M[r], M[pr] = M[pr], M[r]
        rhs[r], rhs[pr] = rhs[pr], rhs[r]
        for i in range(d):
            if i != r and M[i][c]:
                for j in range(d):
                    M[i][j] ^= M[r][j]
                rhs[i] ^= rhs[r]
        piv.append(c)
        r += 1
    # consistency
    for i in range(r, d):
        if rhs[i]:
            return None
    x = 0
    for i, c in enumerate(piv):
        if rhs[i]:
            x |= 1 << c
    assert F.mul(x, x) ^ x == gamma
    return x


def mult_order(field, a, group_order: int) -> int:
    """Multiplicative order of a in field^*, group_order = q - 1."""
    o = group_order
    for p, e in factorint(group_order).items():
        for _ in range(e):
            if field.eq(field.pow(a, o // p), field.one):
                o //= p
            else:
                break
    assert field.eq(field.pow(a, o), field.one)
    return o


# ---------------------------------------------------------------------------
# The variety
# ---------------------------------------------------------------------------


@dataclass
class PointOrbit:
    """One Frobenius orbit of V(A,B): representative (alpha, beta)."""
    fx: int            # min poly of alpha over F_2
    dfx: int           # deg fx
    D: int             # [F_2(alpha,beta):F_2] = orbit size
    field: object      # F2k or F2kExt containing both coords
    alpha: object      # element of field
    beta: object       # element of field
    ord_a: int         # multiplicative order of alpha
    ord_b: int         # multiplicative order of beta
    beta_minpoly: int  # min poly of beta over F_2
    plane_mult: int = -1   # local intersection multiplicity (plane)
    tag: str = ""

    def contributes(self, l: int, p: int, d: int) -> bool:
        """Character kills L = <(l,0),(d,p)>: alpha^l = 1, alpha^d beta^p = 1."""
        F = self.field
        if not F.eq(F.pow(self.alpha, l), F.one):
            return False
        v = F.mul(F.pow(self.alpha, d % self.ord_a),
                  F.pow(self.beta, p % self.ord_b))
        return F.eq(v, F.one)


def sylvester_resultant_y() -> int:
    """Res_y(A~, B^) in F_2[x], exact 5x5 determinant."""
    x3 = 0b1000       # x^3
    xx2 = 0b110       # x + x^2
    # A~ = y^2 + y + x^3: coeffs (by y-degree desc) [1, 1, x^3]
    # B^ = (x+x^2) y^3 + 0 y^2 + 0 y + 1: [x+x^2, 0, 0, 1]
    a = [1, 1, x3]
    b = [xx2, 0, 0, 1]
    n = 5
    M = [[0] * n for _ in range(n)]
    for i in range(3):          # deg B^ = 3 rows of A~
        for j, c in enumerate(a):
            M[i][i + j] = c
    for i in range(2):          # deg A~ = 2 rows of B^
        for j, c in enumerate(b):
            M[3 + i][i + j] = c

    def det(mat):
        m = len(mat)
        if m == 1:
            return mat[0][0]
        r = 0
        for j in range(m):
            if mat[0][j] == 0:
                continue
            minor = [[mat[i][jj] for jj in range(m) if jj != j]
                     for i in range(1, m)]
            r ^= pmul(mat[0][j], det(minor))
        return r

    return det(M)


def variety_orbits(verbose: bool = False) -> list[PointOrbit]:
    """All Frobenius orbits of V(A~, B^) on the torus, exactly."""
    res = sylvester_resultant_y()
    fac = factorize_f2(res)
    orbits: list[PointOrbit] = []
    for f, e in sorted(fac.items(), key=lambda kv: (pdeg(kv[0]), kv[0])):
        df = pdeg(f)
        F = F2k(f)
        alpha = F.t
        gamma = F.pow(alpha, 3)                      # alpha^3
        c = F.add(alpha, F.mul(alpha, alpha))        # alpha + alpha^2
        assert c != 0, "lc of B^ vanishes at a resultant root"
        cinv = F.inv(c)
        b0 = artin_schreier_solve(F, gamma)
        found_here = 0
        if b0 is not None:
            for beta in (b0, b0 ^ 1):
                # B^ condition: (alpha+alpha^2) beta^3 = 1
                if F.mul(c, F.pow(beta, 3)) != 1:
                    continue
                orb = _make_orbit(f, df, F, alpha, beta, D=df)
                orbits.append(orb)
                found_here += 1
        else:
            # beta in K = F_q[s]/(s^2+s+gamma); fiber = {s, s+1}, one orbit
            K = F2kExt(F, gamma)
            beta = K.s
            aK = K.embed(alpha)
            cK = K.embed(c)
            if K.eq(K.mul(cK, K.pow(beta, 3)), K.one):
                orb = _make_orbit(f, df, K, aK, beta, D=2 * df, ext=True)
                orbits.append(orb)
                found_here += 1
            else:
                beta2 = K.add(beta, K.one)
                if K.eq(K.mul(cK, K.pow(beta2, 3)), K.one):
                    orb = _make_orbit(f, df, K, aK, beta2, D=2 * df, ext=True)
                    orbits.append(orb)
                    found_here += 1
        assert found_here >= 1, f"resultant root with empty fiber: {poly_str(f)}"
        if verbose:
            print(f"  factor {poly_str(f)} (deg {df}, exp {e}): "
                  f"{found_here} orbit(s)")
    # verify the defining equations end-to-end for every orbit
    for o in orbits:
        F = o.field
        lhs = F.add(F.mul(o.beta, o.beta), o.beta)              # beta^2+beta
        rhs = F.pow(o.alpha, 3)
        assert F.eq(lhs, rhs), "A~ fails"
        c = F.add(o.alpha, F.mul(o.alpha, o.alpha))
        assert F.eq(F.mul(c, F.pow(o.beta, 3)), F.one), "B^ fails"
    return orbits


def _minpoly(field, elt, maxdeg: int) -> int:
    """Min poly of elt over F_2 (field elements as ints or int-pairs)."""
    def to_bits(z):
        if isinstance(z, tuple):
            a, b = z
            return a | (b << field.base.d)
        return z

    pows = [field.one]
    for _ in range(maxdeg):
        pows.append(field.mul(pows[-1], elt))
    basis: list[int] = []
    piv: list[int] = []
    combo: list[int] = []  # combo[i] = mask over powers building basis[i]
    for k, pw in enumerate(pows):
        v = to_bits(pw)
        m = 1 << k
        for i in range(len(basis)):
            if (v >> piv[i]) & 1:
                v ^= basis[i]
                m ^= combo[i]
        if v == 0:
            return m  # bits of the min poly (bit j = coeff of X^j)
        basis.append(v)
        piv.append(v.bit_length() - 1)
        combo.append(m)
    raise AssertionError("no min poly found within maxdeg")


def _make_orbit(f, df, field, alpha, beta, D, ext=False) -> PointOrbit:
    qa = (1 << df) - 1
    ord_a = mult_order(field, alpha, qa)
    qb = (1 << (2 * df if ext else df)) - 1
    ord_b = mult_order(field, beta, qb)
    bm = _minpoly(field, beta, 2 * df + 1)
    return PointOrbit(fx=f, dfx=df, D=D, field=field, alpha=alpha,
                      beta=beta, ord_a=ord_a, ord_b=ord_b, beta_minpoly=bm)


# ---------------------------------------------------------------------------
# The local multiplicity engine
# ---------------------------------------------------------------------------


def local_dim(orbit: PointOrbit, Nu: int, Nv: int,
              ex_x: tuple[int, int], ex_y: tuple[int, int],
              mod_exp: bool) -> int:
    """dim over k(P) of k(P)[u,v]/(u^Nu, v^Nv, A_loc, B_loc) where
    x |-> alpha (1+u)^{ex_x[0]} (1+v)^{ex_x[1]},
    y |-> beta  (1+u)^{ex_y[0]} (1+v)^{ex_y[1]}.
    If mod_exp, exponents of (1+u) are reduced mod Nu (valid when Nu is
    the 2-part group order, since (1+u)^{Nu} = 1 there), likewise v.
    """
    F = orbit.field
    nb = Nu * Nv

    def zeros():
        return [F.zero] * nb

    def add_into(dst, src, scale):
        if F.eq(scale, F.zero):
            return
        for idx in range(nb):
            if not F.eq(src[idx], F.zero):
                dst[idx] = F.add(dst[idx], F.mul(scale, src[idx]))

    def mul_ring(p, q):
        out = zeros()
        for i1 in range(Nu):
            base1 = i1 * Nv
            for j1 in range(Nv):
                c = p[base1 + j1]
                if F.eq(c, F.zero):
                    continue
                for i2 in range(Nu - i1):
                    base2 = (i1 + i2) * Nv
                    row_q = i2 * Nv
                    for j2 in range(Nv - j1):
                        d = q[row_q + j2]
                        if F.eq(d, F.zero):
                            continue
                        out[base2 + j1 + j2] = F.add(
                            out[base2 + j1 + j2], F.mul(c, d))
        return out

    def one_plus_u_pow(e: int, axis: int):
        """(1+u)^e (axis 0) or (1+v)^e (axis 1) as a ring element."""
        N = Nu if axis == 0 else Nv
        if mod_exp:
            e %= N
        assert e >= 0, "negative exponent without group reduction"
        out = zeros()
        out[0] = F.one
        gen = zeros()
        gen[0] = F.one
        if axis == 0:
            if Nu > 1:
                gen[1 * Nv] = F.one
        else:
            if Nv > 1:
                gen[1] = F.one
        # if truncation is 1 along the axis, (1+u) == 1
        b = gen
        ee = e
        while ee:
            if ee & 1:
                out = mul_ring(out, b)
            b = mul_ring(b, b)
            ee >>= 1
        return out

    def transported(supp):
        acc = zeros()
        for (i, j) in supp:
            scale = F.mul(F.pow(orbit.alpha, i), F.pow(orbit.beta, j))
            eu = i * ex_x[0] + j * ex_y[0]
            ev = i * ex_x[1] + j * ex_y[1]
            term = mul_ring(one_plus_u_pow(eu, 0), one_plus_u_pow(ev, 1))
            add_into(acc, term, scale)
        return acc

    A_loc = transported(AT_SUPP)
    B_loc = transported(BH_SUPP)

    # ideal rows: A*u^a v^b and B*u^a v^b for all monomials (shifts)
    def shifts(p):
        rows = []
        for a in range(Nu):
            for b in range(Nv):
                out = zeros()
                any_nz = False
                for i in range(Nu - a):
                    for j in range(Nv - b):
                        c = p[i * Nv + j]
                        if not F.eq(c, F.zero):
                            out[(i + a) * Nv + (j + b)] = c
                            any_nz = True
                if any_nz:
                    rows.append(out)
        return rows

    rows = shifts(A_loc) + shifts(B_loc)
    # Gaussian elimination over F
    piv_of_row: list[int] = []
    basis: list[list] = []
    for row in rows:
        row = row[:]
        for bi, pv in zip(basis, piv_of_row):
            c = row[pv]
            if not F.eq(c, F.zero):
                for idx in range(nb):
                    if not F.eq(bi[idx], F.zero):
                        row[idx] = F.add(row[idx], F.mul(c, bi[idx]))
        pv = next((idx for idx in range(nb)
                   if not F.eq(row[idx], F.zero)), None)
        if pv is None:
            continue
        cinv = F.inv(row[pv])
        row = [F.mul(cinv, c) for c in row]
        basis.append(row)
        piv_of_row.append(pv)
    return nb - len(basis)


def plane_multiplicity(orbit: PointOrbit, N: int = 10) -> int:
    """Local intersection multiplicity of (A~, B^) at the point, via
    truncation at u^N, v^N (stabilized: caller should confirm with a
    second N)."""
    return local_dim(orbit, N, N, ex_x=(1, 0), ex_y=(0, 1), mod_exp=False)


# ---------------------------------------------------------------------------
# Frames: SNF (verbatim convention of a40_s4_phase_triage.py) + k formula
# ---------------------------------------------------------------------------


def snf2(M):
    M = [list(map(int, r)) for r in M]
    U = [[1, 0], [0, 1]]
    V = [[1, 0], [0, 1]]

    def rowop(i, j, q):
        for c in range(2):
            M[i][c] -= q * M[j][c]
            U[i][c] -= q * U[j][c]

    def colop(i, j, q):
        for r in range(2):
            M[r][i] -= q * M[r][j]
            V[r][i] -= q * V[r][j]

    def swap_rows():
        M[0], M[1] = M[1], M[0]
        U[0], U[1] = U[1], U[0]

    def swap_cols():
        for r in range(2):
            M[r][0], M[r][1] = M[r][1], M[r][0]
            V[r][0], V[r][1] = V[r][1], V[r][0]

    for _ in range(200):
        if M[1][0] == 0 and M[0][1] == 0:
            break
        if M[0][0] == 0:
            if M[1][0]:
                swap_rows()
            else:
                swap_cols()
            continue
        if M[1][0] % M[0][0] == 0 and M[0][1] % M[0][0] == 0:
            rowop(1, 0, M[1][0] // M[0][0])
            colop(1, 0, M[0][1] // M[0][0])
            if M[1][0] == 0 and M[0][1] == 0:
                break
            continue
        if abs(M[1][0]) and abs(M[1][0]) < abs(M[0][0]):
            swap_rows()
            continue
        if abs(M[0][1]) and abs(M[0][1]) < abs(M[0][0]):
            swap_cols()
            continue
        if M[1][0]:
            rowop(1, 0, M[1][0] // M[0][0])
        if M[0][1]:
            colop(1, 0, M[0][1] // M[0][0])
    if M[0][0] and M[1][1] % M[0][0]:
        colop(0, 1, -1)
        D2, U2, V2 = snf2(M)
        U3 = [[sum(U2[i][t] * U[t][j] for t in range(2)) for j in range(2)]
              for i in range(2)]
        V3 = [[sum(V[i][t] * V2[t][j] for t in range(2)) for j in range(2)]
              for i in range(2)]
        return D2, U3, V3
    return M, U, V


def v2(n: int) -> int:
    a = 0
    while n % 2 == 0:
        n //= 2
        a += 1
    return a


@dataclass
class FrameData:
    l: int
    p: int
    d: int
    o1: int
    o2: int
    a1: int          # v2(o1)
    a2: int          # v2(o2)
    x_img: tuple[int, int]   # image of x in Z_o1 x Z_o2
    y_img: tuple[int, int]
    x2: tuple[int, int]      # image of x in the 2-part
    y2: tuple[int, int]


def frame_data(l: int, p: int, d: int) -> FrameData:
    M = [[l, 0], [d, p]]
    D, U, V = snf2(M)
    o1, o2 = abs(D[0][0]), abs(D[1][1])
    assert o1 * o2 == l * p
    xi = (V[0][0] % o1, V[0][1] % o2)
    yi = (V[1][0] % o1, V[1][1] % o2)
    a1, a2 = v2(o1), v2(o2)
    x2 = (V[0][0] % (1 << a1), V[0][1] % (1 << a2))
    y2 = (V[1][0] % (1 << a1), V[1][1] % (1 << a2))
    return FrameData(l, p, d, o1, o2, a1, a2, xi, yi, x2, y2)


def capped_local_dim(o: PointOrbit, a1: int, a2: int,
                     x2: tuple[int, int], y2: tuple[int, int]) -> int:
    """local_dim at the 2-part shape (2^a1, 2^a2), with truncation caps.

    Soundness of the caps:
      * If both a_i >= 1 (both 2-axes live), the plane ideal at the
        point is m-primary of colength = plane_mult =: c, so m^c is
        contained in it (Nakayama chain), hence u^{2^a'} and v^{2^a'}
        are redundant relations once 2^{a'} >= c.  Our orbits have
        c <= 2, so capping either axis at a' >= 1 is exact; we cap at
        3 for slack.  (The substitution (u,v) |-> group images is an
        automorphism-composed evaluation precisely when both axes are
        live with unimodular exponent matrix; when it is ramified the
        argument below applies instead.)
      * If one axis is dead (a_i = 0) the local system is univariate
        (a ramified line evaluation) and the ideal is generated by
        u^{m0}·unit with m0 = the contact order, which the colength
        argument does not bound; we cap at 5 (ring dim 32) and VERIFY
        stabilization dim(cap) == dim(cap+1) whenever the cap binds,
        raising if it fails.
    """
    CAP2, CAP1 = 3, 5
    if a1 >= 1 and a2 >= 1:
        e1, e2 = min(a1, CAP2), min(a2, CAP2)
        return local_dim(o, 1 << e1, 1 << e2, ex_x=x2, ex_y=y2,
                         mod_exp=True)
    # univariate (or fully odd) 2-part
    e1, e2 = min(a1, CAP1), min(a2, CAP1)
    dim = local_dim(o, 1 << e1, 1 << e2, ex_x=x2, ex_y=y2, mod_exp=True)
    if (a1 > CAP1) or (a2 > CAP1):
        f1 = min(a1, CAP1 + 1)
        f2 = min(a2, CAP1 + 1)
        dim2 = local_dim(o, 1 << f1, 1 << f2, ex_x=x2, ex_y=y2,
                         mod_exp=True)
        assert dim2 == dim, ("univariate cap not stabilized",
                             a1, a2, dim, dim2)
    return dim


def spectral_k(orbits: list[PointOrbit], l: int, p: int, d: int,
               detail: bool = False):
    """k of the BB code on Z^2/<(l,0),(d,p)> by the spectral formula."""
    fd = frame_data(l, p, d)
    total = 0
    contrib = []
    for o in orbits:
        if not o.contributes(l, p, d):
            continue
        m = capped_local_dim(o, fd.a1, fd.a2, fd.x2, fd.y2)
        if m:
            contrib.append((o, m))
        total += o.D * m
    k = 2 * total
    if detail:
        return k, contrib, fd
    return k


def transported_supports(l: int, p: int, d: int):
    """(suppA, suppB, fd, collided): honest mod-2 transported supports of
    A~ and B^ on Z_o1 x Z_o2."""
    M = [[l, 0], [d, p]]
    D, U, V = snf2(M)
    o1, o2 = abs(D[0][0]), abs(D[1][1])
    fd = frame_data(l, p, d)
    collided = [False]

    def tr(supp):
        cnt: dict[tuple[int, int], int] = {}
        for e in supp:
            key = ((e[0] * V[0][0] + e[1] * V[1][0]) % o1,
                   (e[0] * V[0][1] + e[1] * V[1][1]) % o2)
            c = cnt.get(key, 0) ^ 1
            cnt[key] = c
            if c == 0:
                collided[0] = True
        return frozenset(k for k, v in cnt.items() if v)

    sa = tr(AT_SUPP)
    sb = tr(BH_SUPP)
    return sa, sb, fd, collided[0]
