# Modified x-cylinder: exact normal form, 48 class floors, and peeling obstructions

**Historical partial; full floor now closed.** The subsequent all-class
proof in `IBM_ONE_CHECK_CYLINDER_TANGENT.md`, assembled with the torus
reduction in `IBM_MODIFIED_EXACT_DISTANCE_THEOREM.md`, proves
`d_cyl(6r)=4r` and the modified square family `[[72r^2,12,4r]]`
for every r>=1. The normal-form theorems and metric counterexamples
below remain valid; they obstruct particular proof steps, not that
distance theorem. The unchanged IBM distance problem remains open.

2026-09-08. Bounded analytic pass. This note does **not** prove or
disprove the full cylinder floor. It proves a nonexpansive class
certificate, gives an exact algebraic normal form, and explicitly
refutes two proposed metric uses of that normal form. No support,
coset, distance, or circumference search is used.

## 1. Hypotheses and proved scope

Let ell=6r, r>=1, and put

```
S=F2[x]/(x^ell+1),       C=S[y,y^-1],
chi=1+x+x^2,       U=y+x,       V=y+x+1,
P=x^3+y+y^2=UV+x chi,       Q=chi V,
D(f,g)=Pf+Qg,       delta(h)=(Qh,Ph).
```

Every vector has finite Laurent y support. Weight counts binary
physical coefficients, with x periodic and y unbounded. Define

```
W=(1+x) sum_(j=0)^(ell/3-1) x^(3j),       d0=2ell/3.
```

The already proved stripe `(0,W)` is a nontrivial cylinder cycle
of weight d0. The parent square-to-cylinder equality is retained;
it does not evaluate the cylinder distance.

This pass proves:

1. Every class has the unique normal representative

```
(chi b, c+U b),       b in ann_S(chi^2), c in ann_S(chi).
```

2. The map

```
L(f,g)=x f(x,x)+g(x,x)
```

is a physical Hamming contraction, is zero on boundaries, and takes
cycles into ann_S(chi). It proves the desired floor for **48 of the
63 nonzero binary classes**, at every ell=6r.
3. Every nonzero g-only cycle has weight at least d0.
4. Every nonzero normal representative has weight at least d0, but
the fixed normal-form map is **not** nonexpansive: it sends a genuine
weight-d0 coset leader to a representative of weight 2d0.
5. Highest-row peeling is also not monotone, even on cycles below
d0. At ell=60 an explicit boundary of weight 38 is sent to weight
40 by the prescribed top-row peel, losing the strict threshold
needed for a contradiction proof.

At the end of this partial pass, the remaining metric obligation
concerned the fifteen nonzero classes in ker L. The limitation of
this normal-form approach is not the existence of the normal form;
it is control of the physical cost of reaching one. The subsequent
shear-and-seam proof closes the floor without that metric operation.

## 2. Elementary ring facts and the exact normal form

Write ann for the annihilator inside S. Since chi^2 divides
x^ell+1,

```
ann(chi)=(W),                 dim_F2 ann(chi)=2,
ann(chi^2)=((x^ell+1)/chi^2),  dim_F2 ann(chi^2)=4.    (2.1)
```

These follow by cancellation in F2[x] before taking the quotient.
Every nonzero element of ann(chi) has exactly d0 sites: its binary
coefficient recurrence is `a_j+a_(j-1)+a_(j-2)=0`, hence repeats one
of the three triples 110,101,011. In particular x^3 W=W in S.

Both P and V are non-zero-divisors in C. For a nonzero finite Laurent
polynomial, multiplication by either has a highest y coefficient
equal to the original highest coefficient, since both are monic.
Also P has unit constant coefficient x^3. Thus

```
C/(P)=S direct-sum S y,
```

with y already invertible in this quotient: `y(y+1)=x^3`.

**Proposition 2.1.** Every cylinder cycle is uniquely a boundary plus
one representative `(chi b,c+U b)` with b,c as in (2.1).

Proof. Reduce g modulo P to `a+y b`, writing `g=Ph+a+y b` with
h a finite Laurent polynomial. Add delta(h) to the cycle. Its new
first block f' satisfies `P f'+Q(a+y b)=0`.

The second summand has y degrees between zero and two. A negative
lowest y degree in f' is impossible because the constant coefficient
of P is the unit x^3. A positive highest y degree in f' is impossible
because P is monic of degree two. Therefore f' is independent of y.
The y^2 coefficient gives `f'=chi b`. Put `c=a+x b`. The remaining
coefficients give

```
chi c=0,       x chi^2 b=0.
```

For example, their full identity before these vanish is

```
D(chi b,c+U b)=x chi^2 b+chi V c.                    (2.2)
```

This proves existence and the converse. For uniqueness, if a normal
representative were a boundary, its second block `a+y b` would equal
Ph. The same lowest/highest degree argument forces h=0. Then a=b=0
and the first block also vanishes. QED.

This gives a direct six-dimensional homology description without a
metric claim. In particular no inverse of chi is assumed to exist.

For completeness, every physical f row of every cycle is divisible
by chi. Evaluating x at a root omega of chi gives
`(1+y+y^2) f(omega,y)=0` in the Laurent domain F4[y,y^-1], hence
`f(omega,y)=0`. This argument is valid because 3 divides ell and
because y is unbounded. It does not say that the f row lies in
ann(chi); divisibility and annihilation are different statements.

## 3. A genuinely nonexpansive homology certificate

Evaluation y=x is defined on C because x is a unit. It sends each
physical monomial x^i y^j to the single coordinate x^(i+j) modulo
ell, including negative j. Consequently parity summation gives

```
|L(f,g)| <= |f|+|g|,
L(f,g)=x f(x,x)+g(x,x).                              (3.1)
```

Here multiplication by x is a physical coordinate permutation.
The needed polynomial specializations are exactly

```
P(x,x)=x chi,       Q(x,x)=chi.
```

Thus a cycle obeys `chi L(f,g)=0`. A boundary obeys
`L(Qh,Ph)=x chi h(x,x)+x chi h(x,x)=0`.

On Proposition 2.1's coordinates the class map is

```
L(chi b,c+U b)=x chi b+c.                            (3.2)
```

This is in ann(chi) and is surjective, since c is arbitrary there.
Its kernel has dimension four. Of the 64 total classes, 48 therefore
have nonzero L; for every representative of any such class, (2.1)
and (3.1) prove

```
|f|+|g| >= |L(f,g)| = d0.                           (3.3)
```

There is no normalization or metric-equivalence assumption in this
proof. The upper stripe has L=W, so its class minimum is exactly d0.

If f=0 in any cycle, `Q g=V(chi g)=0`. Injectivity of multiplication
by V implies chi g=0, row by row. Every occupied g row then has
weight d0, proving the g-only assertion for arbitrary finite y
support. Similarly, a nonzero f-only cycle is impossible because
P is a non-zero-divisor.

## 4. The fixed algebraic normal form can double a true minimum

Every nonzero normal representative has weight at least d0. If
chi b is nonzero, its first-block row is in ann(chi) and already
has that weight. If chi b=0 and b is nonzero, the y-row of the
second block is b in ann(chi), again of weight d0. The remaining
case b=0,c!=0 is equally immediate.

That observation does not prove the full floor. For every ell=6r,
consider the actual cycle

```
z=(0,y^2 W),       |z|=d0.
```

It is nontrivial and is a coset leader by (3.3), since
`L(z)=x^2 W!=0`. Its unique normal representative from §2 is

```
z_normal=(0,(1+y)W),       |z_normal|=2d0.            (4.1)
```

Indeed QW=0 and `PW=(1+y+y^2)W`, using x^3W=W. Therefore
`z+z_normal=delta(W)`. This disproves nonexpansivity of the fixed
normal-form map, even when its input is a nontrivial class-minimum
representative. It does not rule out choosing an adaptive y origin
or a different representative system with a proved metric property.

## 5. Exact top-peeling cost and a below-threshold obstruction

Let the highest nonzero g row of an actual cycle be n, with
`g_n=G`. The cycle equation forces

```
f_j=0 for j>=n,       f_(n-1)=chi G.
```

The usual highest-row peel adds the boundary generated by
`y^(n-2)G`, removing g_n and f_(n-1). Write

```
A=(1+x^3)G,
F=f_(n-2),       H=g_(n-1),       K=g_(n-2),
I=|F intersect A|+|H intersect G|+|K intersect x^3G|.
```

All masks and intersections here are physical x-periodic binary
rows. The exact change in total weight, after minus before, is

```
Delta=|A|+|G|-|chi G|-2I.                            (5.1)
```

To verify (5.1), the three remaining updated rows are F+A, H+G,
and K+x^3G. Their weight changes are summed, then the removed
weights |chi G| and |G| are subtracted. Multiplication by x^3
preserves |G|. Thus a genuinely sufficient and necessary local
rule for this specific peel to be nonincreasing is

```
2I >= |(1+x^3)G|+|G|-|chi G|.                       (5.2)
```

The cycle equation does not automatically supply (5.2).

### 5.1 A fixed compact boundary whose top peel increases weight

For every ell>=12 divisible by six put

```
G=1+x,       H0=1+x+x^2+x^3+x^4,
h=H0+yG,       e=delta(h).
```

All displayed x exponents, including those in the boundaries, lie
between zero and seven, so no accidental x wrap is used. The
two component boundary weights and their support intersection are

```
|delta(H0)|=26,       |delta(yG)|=12,
|delta(H0) intersect delta(yG)|=7.
```

For a direct coefficient check,

```
(1+x^3)H0=1+x+x^2+x^5+x^6+x^7,
chi H0=1+x^2+x^3+x^4+x^6,
(1+x^3)G=1+x+x^3+x^4,       chi G=1+x^3.
```

The overlap is three sites in f row one, two in g row one, and
two in g row two. Therefore

```
|e|=26+12-14=24,
|e+delta(yG)|=26.                                   (5.3)
```

Here n=3 and g_3=G, so this is exactly the prescribed highest-row
peel, not an arbitrary boundary move. The masks in (5.1) have I=1,
and the formula gives Delta=2. For every r>=7, its input weight
24 is strictly below d0=4r. Thus monotonicity fails even on actual
cycles below the proposed distance threshold.

### 5.2 Even strict threshold preservation can fail

Keep the same compact top component and add a separated boundary:

```
J=y^-5(1+x^3),
e_tilde=delta(H0+yG+J),       ell=60,       d0=40.
```

The boundary delta(J) has weight fourteen and has syndrome-qubit
rows disjoint from the component in (5.3): its first block occupies
rows -5,-4 and its second block rows -5,-4,-3. In fact its first
block has weights two and six, and its second block has three rows
of weight two. The same top peel therefore gives

```
|e_tilde|=38<40,
|e_tilde+delta(yG)|=40.                              (5.4)
```

Both vectors in (5.3) and (5.4) are explicitly boundaries. They are
**not** counterexamples to the cylinder floor. They invalidate
unconditional highest-row peeling as either a nonexpansive move
or a move guaranteed to stay strictly below the target. A proof
using peeling must establish a rule such as (5.2), change the
order or grouping of moves, or supply a separate global cost bound.

## 6. Historical residual, subsequent closure, and resource scope

By (3.2), the classes not settled by this partial proof have normal
representatives

```
(chi b, (U+x chi)b),       0!=b in ann(chi^2).       (6.1)
```

There are fifteen such classes: three with chi b=0 and twelve
with chi b!=0. Controlling arbitrary boundary additions to (6.1)
with physical weight would suffice for a normal-form proof; measuring
only these representatives would not. Before the subsequent closure,
any hypothetical counterexample below d0 was restricted to these
classes and had both blocks nonzero. This partial pass supplied
neither that bound nor an original counterexample.

The new checker verifies fixed operator identities and the three
explicit constructions at the single prescribed circumference
ell=60, with y unbounded. It checks the contraction labels, exact
boundary relations, weights, intersections, and the cost identity.
The all-size assertions come from the analytic proofs above, not
from this fixed test. No class or support enumeration is performed.

The single run passed in 0.000490 seconds with peak RSS 17,186,816
bytes. It verified weights 40->80 for fixed normalization, 24->26
for the local peel, and 38->40 for threshold preservation. The
script installs the 55-second wall alarm and 50-second CPU limit
and checks the 768 MiB RSS cap.

Subsequent independent read-only review of
`IBM_ONE_CHECK_CYLINDER_TANGENT.md` passed: its physical shear
realizes L as a column augmentation. Below d0, L=0 forces at least
two occupied sites per occupied column, while all sites in one
column forbid at most three common seams. An unused seam therefore
lifts the entire cycle to an exact plane boundary. This applies
at every height and all r, including the fifteen classes in (6.1).
Together with the reviewed cylinder-to-torus reduction, the complete
result is presented in `IBM_MODIFIED_EXACT_DISTANCE_THEOREM.md`.
No nonexpansive normalization or peeling claim is needed. The
separate one-monomial family is not a substitute or a metric-transfer
argument for this now-closed chi V theorem.
