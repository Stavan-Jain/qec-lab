# A42 S10 — relative filtration is algebraic, but its natural normalization expands weight

2026-09-08. Native Codex sidecar, original omega cylinder only.
The A31 classical-shadow/repeated-root ledger and A42 S6/Spectral
Floor torsion and transfer statements were read as claims to audit.
No sibling note or script was changed.

## 1. Outcome and exact scope

The all-height pure floor is **not proved or disproved**. This note
banks a concrete obstruction to one natural proposed proof:

> For every `q=2^a >= 4`, an explicit compact nontrivial cycle has
> weight `3q+3`, but its unique normalization obtained by removing
> only generator coefficients below the top two u-layers has
> weight `9q/2`. The top-two-layer valuation `q-2` is already
> maximal in this class. Thus this normalization is not
> weight-nonincreasing, even at maximal homological valuation.

This is an actual original-pair vector, not a compressed-check
walk or a finite-window candidate. It is **not a floor
counterexample**: a different representative of the same class,
also given explicitly, has weight `3q`.

A second, simpler family has weights `3q+2` before and `4q` after
the same normalization. An auxiliary two-coordinate complex has
an exactly computed relative distance smaller than its top-layer
representative by a factor two, showing why torsion multiplicity
alone cannot furnish a general relative weight-retention theorem.
That auxiliary complex is not the IBM pair.

Throughout, let `omega^2+omega+1=0`, and put

    u=1+t,        q=2^a,
    K_q=F4[t]/(t^q+1)=F4[u]/u^q,
    P=x^3+omega t+omega^2 t^2,
    Q=1+(x+x^2)t^3,
    D(f,g)=Pf+Qg,       delta(h)=(Qh,Ph).

The chain modules also have the Laurent variable x. All x
supports are finite. Write `Z_q=ker D`, `B_q=im delta`, and
`H_q=Z_q/B_q`. Weight counts nonzero F4 coefficients in the
physical t basis in both blocks, not nonzero u coefficients.
The desired statement is

    every z in Z_q minus B_q has |z| >= 3q.              (1.1)

The banked base minima are three at q=1 and six at q=2. This
review does not rerun their certificates or infer higher minima
from them. Binary HM, the opposite orientation, and the IBM
torus lower bound are outside this sidecar.

## 2. A finite local anchor and the exact relative filtration

Here is a self-contained check of the two local multiplicities
needed below. Work over the DVR `O=F4[t]_(t+1)` and
`S=O[x^+/-1]`. Modulo Q, the basis `1,x` is free over O:

    x^2=x+t^-3,        x(x+1)=t^-3.

In particular adjoining x inverse does not change this quotient.
The reduction of P is `A+B x`, where

    A=t^-3+omega t+omega^2 t^2,       B=1+t^-3.

Its multiplication matrix and determinant are

    [ A       B t^-3 ]
    [ B       A+B    ],
    determinant=A^2+AB+B^2 t^-3.

Modulo `u^4`, using `t=1+u`,

    t^-3=1+u mod u^4,
    A=omega^2 u+omega^2 u^2 mod u^4,
    B=u mod u^4,
    determinant=omega u^3 mod u^4.

Thus the entry gcd has valuation one and the determinant has
valuation three. Smith form over O has invariants u and `u^2`,
up to units. Consequently the localized plane module satisfies

    M=S/(P,Q) = O/(u) direct-sum O/(u^2)                 (2.1)

as an O-module. The nonzero determinant also proves that P is
injective modulo Q. Since S is a domain, the Koszul complex
`K(P,Q)` is a free resolution of M; the order of its two entries
does not affect this conclusion.

The one-term period resolution with differential `u^q` now gives

    H_q = Tor_1^S(M,S/u^q) = M[u^q].                   (2.2)

Therefore `dim_F4 H_1=2`, while `dim_F4 H_q=3` for q>=2.
Localization introduced no new cylinder chains: `O/u^q` is
exactly `F4[u]/u^q`, and x remains a finite Laurent variable.

For q>=2, the chain map from height two to height q is

    tau_q(w)=u^(q-2) w.

Its period-resolution lift is multiplication by `u^(q-2)` in
degree zero and the identity in degree one. Thus it induces the
inclusion `M[u^2] -> M[u^q]`, an isomorphism. In the physical
basis,

    u^(q-2)=1+t^2+t^4+...+t^(q-2),
    |tau_q(w)|=(q/2)|w|.                               (2.3)

This is disjoint replication of the two rows, not a claim about
the weight of a different representative in the same class.

More generally, the image in H_q of cycles divisible by `u^s`
is `M[u^(q-s)]`. For q>=2 its associated homology filtration has
only the following nonzero graded dimensions:

    at s=q-2: dimension one,
    at s=q-1: dimension two.                           (2.4)

These are algebraic statements about classes. They do not say
that every physical nontrivial cycle itself has this valuation.

## 3. The natural normalization exists uniquely, without a metric bound

Every generator has a unique u expansion. For q>=2 let
`L_q` be the generator subspace consisting of polynomials in u
of degree strictly less than `q-2`, with Laurent x coefficients.
For q=2 this subspace is zero.

**Proposition 3.1.** For every `z in Z_q`, there is a unique
`h_< in L_q` for which

    R_q(z)=z+delta(h_<) is divisible by u^(q-2).         (3.1)

This defines a linear projection onto `tau_q(Z_2)`. It induces
the inverse transfer identification on homology after dividing
by `u^(q-2)`. No assertion about weight is included.

Proof. The isomorphism (2.2) gives `z=tau_q(w)+delta(h)` for
some `w in Z_2`. Split h into its coefficients below `q-2` and
its remaining two coefficients. Subtracting the first part gives
(3.1). For uniqueness, write

    Q|_(u=0)=c=1+x+x^2,       P|_(u=0)=p=x^3+1.

If h has leading u coefficient h_s, then the first block of its
boundary has leading coefficient `c h_s`, which is nonzero
in the domain `F4[x^+/-1]`. Thus

    valuation_u(delta(h))=valuation_u(h).               (3.2)

The difference of two possible low-degree generators cannot
have a boundary of valuation at least `q-2`. This proves
uniqueness, and also linearity and the projection assertion.
The remaining part of h is itself a transferred generator, so
the quotient identification is the stated one. QED.

In particular all leading coefficients at levels `s<q-2` of
cycles in `Z_q intersect u^s K_q[x^+/-1]^2` are boundary
leading coefficients: subtract (3.1) and use
(3.2). Nontrivial homology can hide entirely behind such a
boundary leading term. This is the obstruction to detecting
nontriviality at the first actual nonzero coefficient.

## 4. Actual obstruction at maximal class valuation

All claims in this section hold for every `q=2^a>=4`.
They require no distance certificate.

### A six-coordinate height-two cycle

In K_2 define

    h=x+omega,
    r=x^2+omega^2 x+omega^2,
    v_2=(h, 1+t r).                                    (4.1)

This has support size six: two coordinates in the first block
and four in the second. To verify the cycle equation, expand
in u, where `u^2=0`:

    P=p+omega u,
    Q=c+(x+x^2)u,
    1+t r=(1+x)h+u r.

The constant syndrome vanishes because `p=(1+x)c`. The linear
syndrome vanishes by the finite identity

    c r=(omega+x+x^3)(x+omega)
       =[omega+(x+x^2)(1+x)]h.                         (4.2)

Put `v_q=tau_q(v_2)`, so `|v_q|=3q`.
Its first block is supported in exactly two x columns, whereas
the first block Qg of a nonzero compact boundary has width at
least two, meaning at least three possible columns from its
lowest to highest exponent. Indeed Q has constant coefficient
one and leading x coefficient `t^3`, both units. Its two extreme
columns cannot cancel. Thus v_q is nontrivial.

More precisely its class has maximal possible valuation `q-2`.
If it were homologous to a cycle of valuation at least `q-1`,
(3.2) would force the correcting generator to begin at `q-2`.
The first block at that layer would then require

    c h_0=x+omega.

Evaluation at `x=omega^2` rules this out: the left side is zero
and the right side is one. This excludes a higher-valuation
representative in the entire class, not only in a chosen window.

### One local boundary makes normalization expand

Define the actual compact cycle

    z_q=v_q+delta(t^(q-2)).                            (4.3)

It is nontrivial because v_q is. The local boundary's first
block consists of

    (x=0,row=q-2,coefficient=1),
    (x=1,row=1,coefficient=1),
    (x=2,row=1,coefficient=1).

The first meets v_q with coefficient omega and changes it to
omega^2; the other two are new. Hence the first-block weight
increases by two. The second boundary block consists of

    (x=3,row=q-2,coefficient=1),
    (x=0,row=q-1,coefficient=omega),
    (x=0,row=0,coefficient=omega^2).

Only the first is new. The latter two meet coefficients
omega^2 and one, respectively, and do not cancel them.
Consequently

    |z_q|=3q+3.                                       (4.4)

The polynomial `t^(q-2)=(1+u)^(q-2)` has top u coefficient
one and no coefficient at `q-1`. Its low part is exactly

    h_<=t^(q-2)+u^(q-2).

Proposition 3.1 therefore forces

    R_q(z_q)=tau_q(v_2+delta_2(1)).                     (4.5)

The height-two vector in (4.5) is explicitly

    first:  x+omega^2+t(x+x^2),
    second: x^3+omega+t(x^2+omega^2 x+1).

These blocks have weights four and five. Disjoint replication
gives

    |R_q(z_q)|=9q/2 > 3q+3=|z_q|.                     (4.6)

Thus removing only the lower-order generator coefficients is
not weight-nonincreasing. Since this class cannot be moved above
`q-2`, the example also defeats this normalization when stopped
at the class's maximal filtration level.

Allowing an additional top-layer boundary is a different issue:
adding `tau_q(delta_2(1))` to (4.5) restores v_q of weight 3q.
Therefore (4.6) does **not** refute the existence of a suitably
chosen nonexpansive normalization modulo top-layer boundaries.

### A simpler supplementary example

Set `N=u^(q-1)=1+t+...+t^(q-1)` and
`mu_q=(N,(1+x)N)`. It is a cycle, and is nontrivial by the
same extreme-column argument for Q. For

    y_q=mu_q+delta(t^(q-2))

the first block gains one coordinate and the second gains one,
so `|y_q|=3q+2`. Its normalization is the transfer of

    first:  t(1+x+x^2),
    second: x^3+x+omega+t(x+omega^2),

which has height-two weight eight. Thus `|R_q(y_q)|=4q`.
This example lies in the higher, valuation-`q-1` class stratum;
the preceding example is needed to rule out the maximal-level
escape from this obstruction.

## 5. Exactly what classical retention does prove

**Lemma 5.1.** Let a nonzero scalar in K_q have exact u valuation
s. Its physical t weight is at least `2^popcount(s)`. Hence if
a compact vector is `u^s w`,

    |u^s w| >= 2^popcount(s) |w mod u|.                 (5.1)

Proof. Induct on q. Write a scalar of t degree below q as
`L+t^(q/2) H`, with both polynomials of degree below q/2.
Since `t^(q/2)=1+u^(q/2)`, if `s>=q/2` then L=H,
the physical weight doubles, and H has valuation `s-q/2`.
If `s<q/2`, then L+H has valuation s, and
`|L|+|H|>=|L+H|` permits induction without a factor two.
These are precisely the two binary-digit cases for the claimed
factor. Apply this scalar result to each x column and each
block whose coefficient in `w mod u` is nonzero. QED.

This proof never changes the vector by a boundary. In (4.3),
the actual vector has valuation zero and leading coefficient
`(c,p)`, an ordinary boundary. Its nontrivial *class* has
valuation `q-2`. Substituting the latter for the former in
(5.1) is invalid, and normalizing first requires the metric
comparison disproved in (4.6).

There is also an immediate ordinary-code obstruction: the
kernel Z_q contains the six-coordinate boundary `delta(1)`
for every q>=2. Its classical minimum distance therefore
cannot supply the growing relative minimum in (1.1).

### An exact relative counterexample outside the IBM pair

Even a length-one torsion module does not imply retention for
homology classes. Over `F4[u]/u^q`, with q>=2 a power of two,
consider the two-coordinate complex

    d_1=(u^2,u),       d_2=(1,u)^T.

Over the untruncated domain it resolves the quotient by u:
`ker d_1=im d_2` and `coker d_1=F4[u]/u`. After truncation,

    z=(a,b) is a cycle iff ua+b=c N for some c in F4,
    z is a boundary iff c=0,
    N=u^(q-1).

Every nontrivial class has a top-layer representative `(0,c N)`
of weight q. Nevertheless its exact relative minimum is q/2.
Indeed

    (a,b)=(c u^(q-2),0)

represents the same class and has weight q/2. Conversely,

    q=|c N|=|ua+b| <= 2|a|+|b| <= 2|(a,b)|

proves that no lighter representative exists. This is an exact
all-height relative-coset calculation, not an extrapolation.
It is not the original Koszul pair and does not have the IBM
base distances. Its purpose is to disprove a general deduction
from torsion heights and classical retention alone.

## 6. The precise residual and verification scope

For the original pair, every nontrivial class is represented
in `tau_q(Z_2)`, and among those representatives the banked
height-two lower bound proves weight at least 3q. The missing
global metric assertion is still

    |tau_q(w)+delta(h)| >= 3q
    for all w in Z_2 minus B_2 and all compact h.        (6.1)

Equation (6.1) is the remaining floor itself, not a new result.
This pass identifies why the natural proposed route to it fails:
the unique elimination restricted to low generator coefficients
expands weight even for a maximal-filtration nontrivial class.

A stronger sufficient repair would allow optimizing the two
remaining generator coefficients: for each z, require some
`b in B_2` such that

    |R_q(z)+tau_q(b)| <= |z|.                          (6.2)

This would prove **per-class** distance scaling from height two,
not just the global floor. It is stronger than (6.1), remains
unproved, and is not refuted by our examples. It must not be
silently assumed under the name "relative weight retention."
A direct proof of (6.1) could avoid the stronger requirement.

The new standard-library diagnostic is
`scripts/a42_s10_codex_relative_obstruction.py`. It imports
read-only arithmetic helpers from S6, checks the local determinant
to order four, and verifies the displayed original-pair vectors,
cycle equations, normalizations, and weights at the single fixed
period q=4. The universal assertions and nontriviality proofs
are the symbolic arguments above, not the diagnostic's PASS.
There is no period census, width search, distance solver, or
claim of an automated minimum-distance certificate.

The initial and extended fixed checks took 0.000085 and 0.000137
seconds internally; the extended check measured peak RSS
18,399,232 bytes. Each run has a 55-second wall alarm, 50-second
CPU limit, and a 768-MiB measured-RSS check. The final replay
after tightening the leading-coefficient check passed in
0.000135 seconds, with peak RSS 18,432,000 bytes. No old
process, git state, memory file, Lean file, or
sibling research file was touched. Only this new note and its
new diagnostic were written.
