# A43-S3: the one-sign modification of the IBM pair

2026-09-08. New-file scope only; no SAT/census, old-process control,
commits, sibling edits, or Lean. This note concerns

    A = 1+y+x^3 y^-1,     Bplus = 1+x+x^-1 y^3.

It changes only the last y exponent of original IBM B. The original
minus-sign family and its 127-divisible exception are not redefined.

**Final disposition.** On `(12t,6t)`, k is twelve for odd t and sixteen
for even t, and the proven distance bounds are `2t<=d_t<=12t`.
The full charge translation lattice is 12 Z^2. Direct monomial
separable-product and reciprocal-trinomial color rewrites are excluded
by the local algebra. No general exact-distance theorem or physical
metric isometry is established; section 5 states the remaining lemma.

## 1. Banked algebra: k is not uniformly twelve

Set

    a = y A = x^3+y+y^2,      b = x Bplus = y^3+x+x^2,
    f = x^2+x+1,             H = (x+1)y+x.

The unit-clearing is a physical coordinate permutation: with the lab
X-cycle convention `Bbar U+Abar V=0`, send `(U,V)` to `(xU,yV)`.
This maps to the standard pair `(a,b)`, its boundaries map after a
common generator shift xy, and every physical Hamming weight is
unchanged. This is a monomial shift, not a Clifford distance assumption.
The standard pair has the exact axis symmetry `(x,y,U,V)->(y,x,V,U)`.

Two finite Laurent identities, independently checked in the diagnostic:

    b+(y+1)a = f H,
    f^3 = P a + Q b,

where

    P = x^3+x^2+1+y+y^2+x^-1 y^2,
    Q = y+x^-1 y+x^-1.

For a hand derivation of the second, first verify

    x f^2 = (x+1)^2 a + H^2+(x+1)H.

Multiply by f, replace fH by the first identity, and simplify the
coefficient of a. In particular all torus zeros have x a nontrivial
cube root of unity; swapping x,y gives the same assertion for y.
There is no order-127 component in this modified ideal.

Over F4 write omega^2+omega+1=0. The four common torus zeros are
`(omega,omega)`, its conjugate, `(omega,omega^2)`, and its conjugate.
The two diagonal zeros are simple: the Jacobian determinant is
`x^2 y^2+1`, nonzero there. At `(omega,omega^2)` the local algebra is
`F4[u]/(u^3)`, with

    x=omega+u,       y=omega^2+omega^2 u.

To check exactness, put x=omega+u, y=omega^2+v. The a equation has
nonzero coefficient of v, so it solves v as a formal series in u.
Its first terms are `v=omega^2 u+u^3+O(u^4)`. Substitution in b gives
`omega^2 u^3+O(u^4)`, which has order exactly three. The quotient is
therefore as stated, and its conjugate gives the other length-three
point. The full Laurent quotient has F2 dimension `2+2*3=8`.

This agrees with [Liang et al., Example 3, equations 20--24](https://arxiv.org/html/2503.03827v2#S2.Ex3):
their displayed eight-element quotient basis implies k_max=16, not 12.
The local computation here supplies the finite-period distinctions.

**Theorem 1 (all rectangular dimensions).** If either ell or m is not
divisible by three, k=0. Otherwise, write alpha=v2(ell), beta=v2(m).
Then

    k = 4(1+min(3,2^alpha,2^beta)).

Indeed, at any surviving root, the period equation in each variable
vanishes to order `2^v2(period)`. At a length-three antidiagonal point,
both x-omega and y-omega^2 are nonzero scalar multiples of u, so the
two period equations truncate to length `min(3,2^alpha,2^beta)`.
Each diagonal point retains length one. Sum the F2 dimensions and
double for CSS homology. The latter identity follows directly from
the finite group-algebra Frobenius pairing: the common annihilator of
a,b has dimension `dim R/(a,b)`, so both check ranks equal
`ell*m-dim R/(a,b)`.

In particular on the requested homothetic gross frames `(12t,6t)`:

    t odd:  [[144t^2,12,d_t]],
    t even: [[144t^2,16,d_t]].

This is a proved dimension law, not a distance law. At `(12,12)` the
source explicitly lists `[[288,16,12]]`; it lists original minus-sign
IBM on that rectangle separately as `[[288,12,18]]` in Table 3.
The modified code must not inherit that original d=18 record.

## 2. Current bounded verification

`scripts/a43_s3_codex_sign_variant.py` verifies the two identities and
the ranks at `(3,3),(6,6),(12,6),(12,12),(24,12),(36,18)`; respectively
k is `8,12,12,16,16,12`. All tests pass in about 0.0054 seconds and
18 MB RSS. These are specified small rank tests, not a parameter census.

It also verifies the banked x-local weight-12 gross seed and its
y-transfer stacks at those applicable frames: weights 12,24,24,36,
all nontrivial by exact rowspace checks. The uniform certificate is
given in section 4. No all-t distance equality is claimed.

## 3. Translation periods and exact obstructions to direct rewrites

### 3.1 The full charge module has period lattice 12 Z^2

Here a translation fixes the full charge module precisely when its
monomial is one in the eight-dimensional Laurent quotient. The identity
`f(x)^3=0`, and its y mirror, imply `x^12=y^12=1`: indeed
`x^12+1=(x+1)^4 f(x)^4`. Thus twelve in either axis suffices.

It is also necessary. At the diagonal and antidiagonal reduced points,
`x^a y^b=1` requires respectively `a+b=0 mod 3` and `a+2b=0 mod 3`.
Thus a=3p and b=3q. At the antidiagonal local point, put c=omega^2.
In F4[u]/u^3 one has

    x^(3p) = 1 + (p mod 2)c u
               + (p+binom(p,2) mod 2)c^2 u^2,
    y^(3q) = 1 + (q mod 2)u
               + (q+binom(q,2) mod 2)u^2.

The coefficient of u in their product is `(p mod 2)c+(q mod 2)`.
Since 1,c are F2-independent, both p and q must be even. Their product
then has u^2 coefficient `(p/2 mod 2)omega+(q/2 mod 2)`, which vanishes
only if both p,q are divisible by four. The same reasoning works for
negative integers using their residues modulo four. Hence the full
period lattice is exactly 12 Z^2; a smaller sheared basis cannot retain
all sixteen logical qubits. Individual charges can have shorter periods.

This explains why homothetic gross frames with odd t retain twelve
logical qubits while even t retain sixteen: the shorter direction 6t
does not always respect the full charge period. It is a boundary
condition effect; it is not an order-127 phenomenon.

### 3.2 No two-factor separable presentation under monomial coordinates

**Proposition 2.** No GL2(Z) change of Laurent variables transforms the
ideal `(a,b)` into `(p(s),q(t))` for univariate Laurent polynomials p,q
over F2. Arbitrary invertible changes of the two ideal generators and
monomial unit rescalings do not evade this obstruction.

**Proof.** The four geometric zeros all have order-three coordinates,
and their local lengths over the algebraic closure are `1,1,3,3`.
A GL2(Z) monomial automorphism preserves both properties and is
injective on the set of points. For a separated ideal the geometric
zero set is the Cartesian product of the root sets of p and q.
Each root set is contained in `{1,omega,omega^2}`, so a product of size
four must be 2-by-2. A Frobenius-stable two-element subset of this set
is necessarily `{omega,omega^2}`. The two conjugate roots of p have
equal multiplicity, as do those of q. All four product points would
therefore have the SAME local length, contrary to `1,1,3,3`. QED.

This rules out a direct separable HGP rewrite of precisely the requested
one-sign modification. It does not rule out blocking cells, introducing
ancillas, or applying a non-monomial Clifford circuit; none of those
operations has been shown here to preserve physical distance exactly.

### 3.3 The reciprocal-trinomial color theorem does not apply directly

The applicable color-code source assumes a reciprocal pair
`1+M+N, 1+M^-1+N^-1`, where M,N are independent Laurent monomials.
It is not a theorem about arbitrary pairs of trinomials. See
[Liang and Chen, Theorem 1 and its hypotheses](https://arxiv.org/html/2510.05211v1).

There is a stronger ideal obstruction in this case. For the reciprocal
pair, the quotient in M,N has exactly two simple geometric points,
`(omega,omega^2)` and its conjugate. If the exponent determinant is
Delta, its full Laurent quotient has dimension `2|Delta|`: the ambient
group algebra is free of rank |Delta| over the monomial sublattice
algebra. To match the sign variant's dimension eight requires |Delta|=4.
An exponent map of determinant four has, by Smith normal form, diagonal
invariants (1,4) or (2,2). Over characteristic two its map of algebraic
tori is purely inseparable: each of the two points has one geometric
preimage, of length four. This would give local lengths `4,4`, whereas
the sign variant has `1,1,3,3`. Thus no such reciprocal-trinomial ideal
is obtained by a monomial automorphism and a change of ideal generators.

Again this is an obstruction to that explicit algebraic presentation,
not a theorem ruling out every possible nonlocal permutation or local
Clifford equivalence. A Clifford equivalence to topological layers, if
constructed, still needs an independent physical-weight argument.

## 4. A symbolic weight-12t upper family on homothetic gross frames

Everything in this section uses the requested Laurent pair, before
unit-clearing. Bars mean inversion of both variables. An X-cycle is
`Bplus_bar U+Abar V=0`, and an X-boundary is `(Abar h,Bplus_bar h)`.
A Z-cycle `(Z0,Z1)` satisfies `A Z0+Bplus Z1=0`.

Use the seed with block-coordinate supports

    L0 = {(1,0),(1,2),(2,1),(2,3),(3,1),(3,5)},
    L1 = {(3,1),(3,5),(4,0),(4,1),(4,2),(4,5)}.

It has the exact plane identity

    Bplus_bar L0 + Abar L1 = (1+y^6) C,
    C = 1+x+x^2(y^-3+y^-1)+x^3 y^-2+x^4(y^-2+y^-1).    (4.1)

The sign matters: C differs from the original minus-sign seed's
syndrome polynomial. Hence L is an x-compact X-cycle on `Z x Z_6`.

An explicit transverse dual tile is

    Z0 = {(3,1),(3,2),(4,1),(4,2),(5,4),(5,5)},
    Z1 = {(0,1),(0,3),(1,0),(2,0),(3,0),(3,1),(4,0),(5,0)}.

Its exact plane identity is

    A Z0 + Bplus Z1 = (1+x^6) E,
    E = x^-1(y^4+y^6)+y+x(1+y)+x^2(y^3+y^4).           (4.2)

Thus the tile is a Z-cycle on `Z_6 x Z`, supported in rows 0 through
5. This 14-qubit tile was found by ONE finite binary linear system
with 72 variables, the cylinder syndrome equations, and the condition
of odd pairing with L. The finite identities, not the solve or sampled
frames, are the certificate used for all parameters.

**Theorem 3 (uniform upper bound).** For all ell,m divisible by six,
with ell,m>=6, the modified code has a nontrivial cycle of weight 2m.
In particular, on `(ell,m)=(12t,6t)`, for every integer t>=1,

    d_t <= 12t.

**Proof.** Write m=6t and set

    L_t = sum_{j=0}^{t-1} y^(6j) L,
    Z_ell = sum_{j=0}^{ell/6-1} x^(6j) Z.

The seed bands in L_t are disjoint and each has twelve qubits, so its
weight is 12t=2m. Equation (4.1) telescopes to a multiple of `1+y^m`;
equation (4.2) telescopes to a multiple of `1+x^ell`. Thus L_t is an
X-cycle and Z_ell is a Z-cycle on the specified torus. The dual meets
only seed band zero, and the seed's columns lie in 1 through 4, so
only dual x-tile zero can intersect. Their common qubits are exactly

    (block,x,y) = (0,3,1), (1,3,1), (1,4,0).

Their pairing is odd. No X-boundary can have odd pairing with a
Z-cycle, so L_t is nontrivial. This holds for even t as well as odd t;
it does not assume k=12. QED.

The actual standard-pair axis symmetry supplies the companion upper
bound 2ell, after transporting to and from the standard coordinates.
Thus these six-divisible rectangles satisfy `d<=2 min(ell,m)`.

The first candidate dual made by reflecting/transposing L is a genuine
cylinder Z-cycle but its six x-translates all pair evenly with the
stack. It does not certify this upper family. The separately solved
14-qubit tile does. This is why simply invoking symmetry of a witness
was insufficient even though symmetry of the modified code is valid.

## 5. A rigorous growing floor and the unresolved exact metric

The standard pair has check footprint spread three in both axes. Its
two polynomials are coprime in the plane Laurent ring: any common
factor would divide f(x)^3 by section 1 and f(y)^3 by its axis mirror,
but those univariate polynomials are coprime to each other there.

**Theorem 4 (elementary floor).** Every nontrivial logical of the sign
variant on an ell-by-m rectangle has weight at least
`ceil(min(ell,m)/3)`. Consequently on homothetic gross frames,

    2t <= d_t <= 12t       for EVERY integer t>=1.       (5.1)

**Proof.** Suppose a nontrivial standard-pair X-cycle has weight
w<min(ell,m)/3. If its occupied columns had no run of three empty
columns, there would be at most two empty columns per occupied column,
forcing ell<=3w; contradiction. The same applies to rows. Cut along
the two empty runs. Since each check's coordinate spread is three,
the support lifts injectively to a compact plane cycle, without any
check coupling across a cut. Coprimality in the Laurent UFD implies
that a plane cycle `(u,v)` solving `bbar u+abar v=0` is a boundary:
abar divides u, bbar divides v, and the two quotients agree. Project
its compact trivializer back to the torus, contradicting nontriviality.
The block-swap/inversion map `(u,v)->(vbar,ubar)` exchanges X and Z
cycles and boundaries without changing weight. QED.

This is a proved linear floor for the requested slight modification,
not an exact constant. No claim that d_t=12t for all t is made.
The published `(12,12)` value 12 does not imply an x-compact or
y-compact weight-12 string on larger rectangles: the wrapping condition
can change. No small witness is extended without a finite syndrome
identity and an all-parameter nontriviality certificate.

For the proposed exact law on `(12t,6t)`, the entire remaining lower
bound is the following concrete statement in PHYSICAL coordinates:

> For a=x^3+y+y^2, b=y^3+x+x^2 on `(12t,6t)`, every binary pair
> `(u,v)` with `bbar u+abar v=0` and `wt(u)+wt(v)<12t` belongs to
> `{(abar h,bbar h): h in the torus group algebra}`.

The upper half and dimensions no longer depend on this statement.
Neither the local algebra above nor a toric-layer phase classification
proves it. A sufficient geometric refinement would give a uniformly
bounded-width normal form for minimum representatives, with a
weight-nonincreasing physical realization and a sharp cost for every
allowed winding and junction. No such completeness theorem was found
in this bounded audit; no exponentially wide scan was attempted.

## 6. Final verification scope

Only this note and `scripts/a43_s3_codex_sign_variant.py` are new
outputs. The script uses standard-library set arithmetic and bitset
Gaussian elimination, with a 55-second alarm, 50-second CPU limit, and
peak-RSS assertion below 768 MiB. It does not launch other processes,
modify old scripts or data, import solver engines, or enumerate logical
cosets. The final check, including the finite translation-action
identities, ran in about 0.0076 seconds with 18.5 MB peak RSS.
Aggregate diagnostic computation, including the
initial witness extraction and all replays, remains below one second.

Reproduce from the repository root:

    uv run --no-cache --offline --frozen --no-sync --project experiments/bb_lab python -B experiments/bb_lab/scripts/a43_s3_codex_sign_variant.py

The source role is limited: Liang et al. provide the standard-pair
identification, their k_max=16 basis, and a separate small `(12,12)`
distance record; Liang and Chen specify the hypotheses for reciprocal
trinomial color codes. The dimension, translation-period, obstruction,
upper-family and elementary-floor proofs here are explicit independent
algebra. No new general exact-distance theorem or novelty claim is
inferred from those references or from the small tests.
