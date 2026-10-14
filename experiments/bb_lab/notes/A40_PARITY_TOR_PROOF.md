# All-size cylinder parity and double-cover kernels for the original IBM pair

2026-09-08. Analytic continuation of A40 S11. This note proves the
class-space statement for all rectangular periods divisible by six,
including the order-127 exception. It does **not** prove a distance
floor or a weight-preserving normal form for a logical representative.

## 1. Statement

Work over `R=F2[x^+/-1,y^+/-1]`, with the original unit-cleared pair

```
P=x^3+y+y^2,       Q=1+(x+x^2)y^3.
```

Write `K(P,Q)` for the complex with differentials `(P,Q)` and
`(Q,P)^T`. Clearing the original Laurent polynomials by monomial
units changes the two physical block bases by translations; it gives
the same homology and intertwines period projection and pullback.

Let `ell,m` be positive multiples of six, and let
`Lambda=<(ell,0),(0,m)>`. Define

```
d0=6+7*[127 divides both ell and m].
```

There is a natural period-coordinate identification

```
H_1(K(P,Q) tensor_R F2[Z^2/Lambda]) = M_Lambda direct-sum M_Lambda,
dim_F2 M_Lambda=d0.                                      (1.1)
```

For a primitive lattice direction
`u=(a ell,b m)`, where `gcd(a,b)=1`, the image of compact cylinder
homology for `Z^2/<u>` in (1.1) is exactly

```
W_u={(a z,b z): z in M_Lambda},                          (1.2)
```

with `a,b` reduced modulo two. Consequently there are exactly three
distinct cylinder-image spaces. Any two have zero intersection and
span the torus homology. There are

```
3(2^d0-1)                         nonzero cylinder-liftable classes,
(2^d0-1)(2^d0-2)                  remaining, mixed classes.
```

For `d0=6` these counts are 189 and 3,906. For `d0=13` they are
24,573 and 67,084,290. The mixed classes have no compact cylinder
lift in **any primitive lattice winding**, not merely the three
displayed representatives.

For the index-two cover retaining `u` and doubling a complementary
lattice direction, pullback has kernel exactly `W_u` and rank `d0`.
This is an all-size version of the cover-kernel calculation used in
S11. It applies independently of the previously tested frames.

## 2. A self-contained finite plane-module calculation

Set `L=F2[x^+/-1]` and `M=R/(P,Q)`. Modulo `P`, the elements
`1,y` form a free L-basis. Indeed `y^2=y+x^3`, and y is already
invertible since `y(y+1)=x^3`. Reducing Q in this basis gives

```
Q=a0+b0 y,
a0=1+x^4+x^5,       b0=x+x^2+x^4+x^5.
```

Multiplication by Q is represented by

```
T = [ a0      b0 x^3 ]
    [ b0      a0+b0  ].
```

Direct polynomial identities give, with `g=x^2+x+1`,

```
a0=g(1+x+x^3),       b0=g(x+x^3),
gcd(a0,b0)=g,
det T = a0^2+a0 b0+x^3 b0^2 = g^3 h,
h=x^7+x^6+x^4+x^2+1.                                  (2.1)
```

The gcd assertion is particularly simple: the displayed quotients
sum to one. Since L is a PID, the first Smith invariant is g and
the second is `g^2 h`. Hence, as an L-module,

```
M = L/(g) direct-sum L/(g^2 h).
```

The g-primary and h-primary summands are invariant under y because
y commutes with x. Denote them by `M0` and `W`. Thus

```
M=M0 direct-sum W,
dim_F2 M0=6,       g^2 M0=0,
W=L/(h),           dim_F2 W=7.                         (2.2)
```

The polynomial h is irreducible. A compact exact certificate is
`x^128=x mod h` and `gcd(h,x^2+x)=1`; since its degree seven is
prime, these are the irreducibility criterion. The accompanying
diagnostic verifies both identities with binary polynomial division.
In particular W is F128. On W, x is neither zero nor one, and y is
also neither zero nor one: Q excludes y=0, while P and invertibility
of x exclude y=1. Both x and y therefore have order 127, since
F128's nonzero group has prime order 127.

On M0, `x^6+1=(x+1)^2 g^2` acts as zero. Squaring P then gives
`y^4+y^2+1=0`, so `(y^2+y+1)^2=0` and also `y^6+1=0`.
Consequently both periods annihilate M0 for every frame in §1.
On W a period acts as zero precisely when 127 divides its length,
and otherwise acts invertibly.

No local multiplicity stabilization or growing matrix calculation
is needed for (2.2). Also, det T is nonzero, so multiplication by Q
on the free L-module `R/(P)` is injective. Since R is a domain,
`P,Q` form a regular sequence. Thus `K(P,Q)` is a free resolution
of M, not just an arbitrary complex with the same degree-zero
quotient. This fact is essential for the following Tor argument.

## 3. Computing the torus homology by the other resolution

Put `f=x^ell+1`, `g_period=y^m+1`. These form a regular sequence
in R: after quotienting by the monic x relation, the monic y
relation is still a non-zero-divisor. Therefore the period Koszul
complex is a free resolution of `R/(f,g_period)`.

Using either free resolution to compute the same balanced Tor group,

```
H_1(K(P,Q) tensor R/(f,g_period))
  = Tor_1^R(M,R/(f,g_period))
  = H_1(M tensor K(f,g_period)).                       (3.1)
```

On M0 both differentials in the last complex vanish. On W, either
both vanish or one period is invertible and the Koszul complex is
contractible. Define `M_Lambda` to be M0 together with W exactly
when both periods annihilate W. Equation (1.1) follows. This also
gives `k=2d0`: in particular `k=12` for the whole b=1 column and
`k=12+14*[127 divides r]` for b=0.

These are homological coordinates. The Tor identification is not
claimed to preserve Hamming weight, locality, or support diameter.

## 4. The cylinder image is determined by winding parity

Let `h_u=x^(a ell)y^(b m)+1`. It is nonzero for the primitive
direction u, so `R --h_u--> R` resolves the cylinder algebra.
Its first Tor group is `ann_M(h_u)`. Laurent polynomials represent
finite-support chains on this cylinder; no completion or infinite
series has been introduced.

For any integer n define the Laurent polynomial

```
S_n(X)=(X^n+1)/(X+1).
```

This includes `S_0=0` and negative n. The quotient is a Laurent
polynomial and `S_n(1)=n mod 2`. With `X=x^ell`, `Y=y^m`, write

```
h_u=U f+V g_period,
U=Y^b S_a(X),       V=S_b(Y).                          (4.1)
```

The quotient map from cylinder to torus is lifted to the free
resolutions by the identity in degree zero and column `(U,V)^T`
in degree one. Naturality of Tor therefore sends

```
z in ann_M(h_u)  ->  (U z,V z) in H_1(M tensor K(f,g_period)).
```

On M_Lambda, X=Y=1, so this map is `(a z,b z)`. Every such z
belongs to the source annihilator. Any W contribution outside
M_Lambda maps to the contractible target summand and is zero in
homology. This proves equality, not merely containment, in (1.2).
Primitive `(a,b)` has nonzero parity, so its map is injective on
M_Lambda. The three possible nonzero parities give the two axes
and the diagonal. Their intersection and counting statements are
now elementary F2-linear algebra.

This is a statement about **logical classes**. It does not assert
that a minimum physical representative of a liftable class has an
injective cylinder support lift. Projection may cancel coordinates.
In particular it does not identify a representative's cyclic gaps
with the rank-one/mixed class distinction, and it does not discharge
the geometric C_min lemma from the original-distance program.

## 5. Double-cover pullback and projection

Complete u to a lattice basis `(u,v)` with determinant ±1 in the
original period basis. Put `F=z^u+1`, `G=z^v+1`, where
`z^(c,d)=x^c y^d`. The base algebra is `R/(F,G)`; retaining u
and doubling v gives `R/(F,G^2)` in characteristic two.
Both pairs are regular sequences: their two independent lattice
relations have finite-dimensional quotients in the two-dimensional
Laurent polynomial ring and thus form complete intersections.

On the physical coefficient modules, projection is the quotient
map and pullback is multiplication by G. The latter is the norm
`1+z^v`: it places one copy on each of the two sheets, is well
defined modulo `(F,G^2)`, and exactly doubles physical support.

The two maps have the following lifts between the period Koszul
resolutions, in degrees zero, one, two respectively:

```
pullback:    G,       diag(G,1),       1;
projection:  1,       diag(1,G),       G.                (5.1)
```

For example `(F,G^2) diag(G,1)=G(F,G)`, and
`diag(G,1)(G,F)^T=(G^2,F)^T`; the other two chain identities
are equally immediate. These are chain maps, not just maps on
quotients with presumed homological behavior.

The same M_Lambda survives for the base and cover. The M0 periods
vanish in either. A W character trivial on the cover lattice but
not the base would give a nontrivial homomorphism from the
order-two quotient to the order-127 character group, which is
impossible. Equivalently, squaring cannot turn a nonzero period
scalar in W into zero.

Tensoring (5.1) with M_Lambda sets F=G=0. Pullback on H1 is
therefore `(z1,z2)->(0,z2)`, while projection is
`(z1,z2)->(z1,0)`. The kernel of pullback and the image of
projection are both the u-cylinder space. Each has dimension d0.
This proves the cover assertion of §1 for every primitive direction.

Let the three covers retain `(ell,0)`, `(0,m)`, `(ell,m)`,
respectively. Every mixed class pulls back nontrivially to all three.
For a representative of such a class with weight w,

```
w >= (1/2) max(d_sector(cover_x),
              d_sector(cover_y), d_sector(cover_diagonal)). (5.2)
```

Here the distances are those of the corresponding CSS sector; the
global quantum cover distance is also a valid weaker substitute.
No value is assumed for those distances. In particular (5.2) does
not establish that every mixed class pays 2m, nor remove a distance
constant from the original conjecture.

## 6. Verification and remaining metric question

`scripts/a40_parity_tor_check.py` checks the small polynomial anchor,
the Smith data, order statements, Laurent identity (4.1) including
negative windings, and both chain-map identities. All check sizes
are fixed independently of ell,m. Universal quantifiers come from
the arguments above, not sampled windings or a finite code census.

The gain is a rigorous all-size separation of cylinder-image and
mixed homology, including the exceptional sector. What remains is
the weight problem: a sharp metric certificate in these classes,
or a weight-compatible map to a simpler complex. The three-cover
comparison is a lower-bound transfer, not itself an exact-distance
proof. A reversed dual gives the analogous homology description
for the other CSS sector; it does not supply any missing weight
inequality in that sector.

Independent review: `A42_S7_spatial_equality.md` §10 audits the
Smith calculation, balanced-Tor naturality, all winding images,
the index-two norm/projection maps, and the exceptional summand.
It finds no gap. A separate fixed polynomial replay passes in
0.014 seconds at approximately 21 MB RSS. The review also states
the metric limitations explicitly; it is not a distance certificate.
