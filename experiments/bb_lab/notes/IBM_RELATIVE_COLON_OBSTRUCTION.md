# Gross calibration: the class-matched colon relaxation is exactly six

2026-09-08. Original IBM pair, fixed gross frame `(ell,m)=(12,6)`.
An algebraic obstruction to a proposed shortcut, not a distance search.

## 1. What is proved

Use the conventions of `A43_S4_block_purity.md`:

    R=F2[x^+/-1,y^+/-1]/(x^12+1,y^6+1),
    P=x^3+y+y^2, Q=1+(x+x^2)y^3,
    D(f,g)=Pf+Qg, delta(h)=(Qh,Ph).

Let K be the common single-block homology image. S4 proves the
canonical identifications

    H/K = (Q:P)/(Q) = (P:Q)/(P).

Write `d_first(eta)` and `d_second(eta)` for the minimum physical
weights in the two component cosets corresponding to the same
nonzero `eta in H/K`. Then

    d_R((Q:P) minus (Q)) = d_R((P:Q) minus (P)) = 3,
    min_(eta != 0) [d_first(eta)+d_second(eta)] = 6.       (1.1)

Thus even the **class-matched**, doubly stabilizer-excluded colon
lower bound is only six at gross, whose banked quantum distance is
twelve. The individual proposed values `6r` and the class-matched
sum `12r` are false already at `r=1`. The class-space theorem itself
is unaffected. Neither an individual numerical witness nor the
colon relaxation is being mistaken for the quantum distance.

## 2. No excluded vector has weight one or two

It suffices to work over F4 with `omega^2+omega+1=0`. Extension of
scalars is faithful for the ideal-membership equalities used here.
All roots of the period polynomials are cube roots of unity.
Their local multiplicities in x and y are four and two.

First consider `R/(Q)`. Its supported geometric points are
`(alpha,beta)`, where `alpha in {omega,omega^2}` and
`beta in {1,omega,omega^2}`. In local coordinates
`x=alpha(1+u), y=beta(1+v)`, Q has the linear expansion
`alpha u+v`; consequently this quotient component is

    F4[v]/(v^2), with u=alpha^2 v.

On it P is a unit when beta=1, is a nonzero multiple of v when
beta=alpha, and is zero when beta=alpha^2. These follow by
substituting into `P=1+beta+beta^2+u+beta v` modulo v^2.

A weight-one polynomial is a unit and cannot be in `(Q:P)`:
at beta=1 both it and P are units. Up to a monomial factor, a
weight-two polynomial is `F=1+x^a y^b`. If `PF=0 mod Q`, then
at beta=1 both the constant and linear coefficients of F vanish:

    alpha^a=1,
    (a mod 2) alpha^2 + (b mod 2)=0.

Hence 3 divides a and both a,b are even. At beta=alpha, P has
valuation one, so the constant coefficient of F must vanish:
`alpha^(a+b)=1`. Therefore 3 also divides b. Both a,b are
divisible by six, which makes F zero in every displayed local
component of `R/(Q)`. Thus F lies in (Q), and is not excluded.

The second block is analogous but needs its own calculation.
The supported points of `R/(P)` have
`alpha in {1,omega,omega^2}`, `beta in {omega,omega^2}`.
Each component is `F4[v]/(v^2)` with `u=beta v`. Q is a unit
at alpha=1, has valuation one at alpha=beta, and is zero at
alpha=beta^2. Vanishing of `F=1+x^a y^b` at the unit component
forces `beta^b=1` and `(a mod 2) beta+(b mod 2)=0`.
Thus 3 divides b and a,b are even. The alpha=beta component
forces 3 to divide a as well. Again F is zero throughout
`R/(P)`, so it is in (P). Weight one is excluded in the same
way as before. Both relative minima are therefore at least three.

This is a fixed algebraic proof using all local components of
the two classical quotients; no support enumeration establishes
the lower bounds.

## 3. Explicit matching weight-three components

Define two complete weight-twelve cycles by their occupied
positions in the bit ordering `block*72+6*x+y`:

    z_A: 0,61,68,75,77,81,134,135,138,141,142,143;
    z_B: 10,12,17,19,26,34,43,44,60,83,84,109.

Their block weights are `(3,9)` and `(9,3)`. Set

    f=1+x^10 y+x^11 y^2,
    g=x y^5+x^2+x^6 y.

These are respectively the first block of z_A and second block
of z_B. Direct binary polynomial multiplication verifies
`D z_A=D z_B=0`. Both complete cycles are nonboundaries.
Their first components are equal modulo (Q), and their second
components are equal modulo (P). Consequently f and g belong
to the two component cosets of the **same** quotient label eta.
That label is nonzero: for instance `f(omega,omega^2)=1`,
where Q vanishes, so f is not in (Q).

All these finite identities, including the two ideal memberships
and nonboundary tests, are replayed by exact binary elimination
in `scripts/ibm_relative_colon_gate.py`. Combined with section 2,
they show both minima for this eta equal three, proving (1.1).

However `Pf+Qg` is nonzero: it has sixteen occupied sites. The
independently inexpensive blocks cannot simply be combined into
a six-site logical. Sharing a quotient homology label does not
guarantee that independently chosen block representatives satisfy
the joint cycle equation. The missing constraint is substantive,
not a removable choice of notation.

## 4. Provenance and bounded verification

The two source cycles were extracted from the already-banked,
census-derived gross population in `data/a43/s1_l3_shape.json`;
no SAT witness or fresh search was used. For these X-cycle
conventions the exact weight-preserving conversion to the cleared
pair is `(U,V) -> (y^-1 Vbar, x^-1 y^-3 Ubar)`. It combines
reflection, exchange of the two blocks, and block translations.
It uses no incorrect axis-mirror equivalence of the original pair.

The diagnostic first passes `validate_banked`, then rechecks all
1,884 stored vectors for weight, zero syndrome, and nonboundary
status. Of them, 84 have both components principal; the others
occupy 51 nonzero quotient labels. Choosing the lightest observed
component separately within each label gives 36 sums of six,
nine sums of eight, and six sums of twelve. Except for (1.1),
these are only witness upper bounds on the relative-coset costs;
no completeness of the observed label minima is inferred.

The initial parent run passed in 0.99 seconds at 47.9 MB peak RSS.
The two displayed source cycles are additionally pinned directly
in the diagnostic, so the exact obstruction has a fixed replay
rather than depending on an unrecorded choice of representatives.
The diagnostic has 55-second wall, 50-second CPU, and 768-MiB
memory guards. No census, growing-parameter scan, or optimization
was run.

The final pinned-source replay passed in 0.92 seconds at 48.1 MB.
An independent read-only review by the S4 agent passed the local
weight-one/two exclusion, the matched-label inference of exact six,
and the distinction from the actual quantum distance. No diagnostic
rerun or new search was needed for that review.

## 5. Remaining scope

S4 gives useful algebraic control of all b=1 class spaces. It
does not by itself produce a growing sharp metric bound. Already
at the base member the separate block relaxation loses a factor
of two, even after matching quotient classes. Any proposed colon
proof of the IBM target must retain enough of the **joint syndrome
compatibility** to close this gap; repeating either excluded
classical bound alone cannot do so. This is not a counterexample
to IBM's conjecture or to a correctly hypothesized literature
theorem.
