# A43 S4 — coincident single-block homology and the doubly excluded colon reduction

2026-09-08. Bounded native Codex sidecar. Original IBM b=1
rectangles, all r>=1. No metric closure or minimum-weight search.

## 1. Theorem and conventions

Let

    ell=6(r+1),       m=6r,       r>=1,
    R=F2[x^+/-1,y^+/-1]/(x^ell+1,y^m+1),
    P=x^3+y+y^2,       Q=1+(x+x^2)y^3,
    D(f,g)=Pf+Qg,       delta(h)=(Qh,Ph),
    H=ker D / im delta.

These are the original Laurent IBM entries cleared by monomial
units. The corresponding changes in the physical chain blocks
are translations. They preserve block support, Hamming weight,
and the principal ideals occurring below. No axis-mirror
equivalence is used.

Define J_first and J_second to be the images in H of cycles
supported in the first and second physical blocks, respectively.
Let K be the image of cycles satisfying both `f in (Q)` and
`g in (P)`. Define the larger-2-adic period direction by

    L=(ell,0) if r is odd,
    L=(0,m)   if r is even.                            (1.1)

This compares valuations at two, not the numerical lengths.

**Theorem.** For every such r,

    J_first = J_second = K = W_L,       dim_F2 K=6,     (1.2)

where W_L is the image of compact cylinder homology wrapping L.
In the period-coordinate identification of A40's parity theorem,
`H=M0 direct-sum M0`, with `dim_F2 M0=6`, it is

    K=M0 direct-sum 0     when r is odd,
    K=0 direct-sum M0     when r is even.               (1.3)

Moreover

    ann_R(PQ)/(ann_R(P)+ann_R(Q))  -> K,
    z -> [(Qz,0)]                                      (1.4)

is an isomorphism. Thus the failure of the semisimple
annihilator identity has exactly dimension six, and fills the
entire common single-block image.

For every cycle, the following are equivalent:

    [(f,g)] in K;       f in (Q);       g in (P).        (1.5)

Consequently a class outside K forces **both** components outside
their respective principal stabilizer ideals, for every
representative of that class. Section 6 proves the resulting
matched relative-coset distance bounds.

There are 63 nonzero classes admitting single-block representatives;
each admits representatives in either block. The other 4,032
nonzero classes do not admit a single-block representative.
Among these are all 3,906 S11 mixed winding classes and the 126
nonzero classes on the other two cylinder parity axes.

Here "single-block" concerns the two physical Koszul coordinates.
It is neither A42 omega-fibre purity nor S11 rank-one/mixed
winding. A representative using both blocks can still belong to
K or another rank-one cylinder class.

## 2. Source scope and plane bookkeeping

The primary source was read directly: [arXiv:2608.27565v1,
Theorem 90](https://arxiv.org/html/2608.27565v1). It assumes a
semisimple ambient ring and uses
`ann(PQ)=ann(P)+ann(Q)` to justify the alternating principal-ideal
exclusion. Its preceding Theorem 88 does not require
semisimplicity. Our result obstructs extending the former
argument to these even rectangles; it does not contradict
either theorem under its stated hypotheses.

For completeness, the finite plane anchor used to discard
irrelevant components is as follows. Modulo P, the basis `1,y`
is free over `F2[x^+/-1]`; y is already invertible because
`y(y+1)=x^3`. In that basis,

    Q=a0+b0 y,
    a0=1+x^4+x^5,       b0=x+x^2+x^4+x^5.

Multiplication by Q has matrix

    [ a0      b0 x^3 ]
    [ b0      a0+b0  ].

With `g0=x^2+x+1` and `h7=x^7+x^6+x^4+x^2+1`, its entry gcd
is g0 and its determinant is `g0^3 h7`. These identities also
prove regularity of the plane pair: the matrix is injective
over the Laurent polynomial domain. The polynomial h7 is
irreducible, as certified by `x^128=x mod h7` and
`gcd(h7,x^2+x)=1`. This is the degree-seven irreducibility test.

The ordinary closed points are the two F4 orbits

    T: (alpha,beta)=(omega,omega),
    G: (alpha,beta)=(omega,omega^2),
    omega^2+omega+1=0.

The remaining orbit W has residue field F128. Its x and y
coordinates both have order 127: reduction of Q gives y as
a field element, and neither coordinate is zero or one.
Hence W can survive a rectangle only if 127 divides both
lengths. It is absent for consecutive r,r+1. We do not discard
it from the plane zero set or claim that the full plane module
has dimension six. Its full dimension is thirteen.

At a component where P or Q is a unit the Koszul complex is
contractible. W is also absent from torus homology when a
period is invertible there. Thus only T and G contribute in
the present family. The established b=0 exception `127|r`,
with k=26, remains untouched; b=0 is outside this theorem.

## 3. Local derivatives and annihilator calculation

Write locally

    x=alpha(1+u),       y=beta(1+v).

At an ordinary point, the actual linear coefficients are

    P=u+beta v+terms of degree at least two,
    Q=alpha u+v+terms of degree at least two.            (3.1)

Equivalently, the ordinary partial derivatives are
`P_x=alpha^2`, `P_y=1`, `Q_x=1`, `Q_y=beta^2`.
All four are nonzero. The normalized coefficients in (3.1),
not the unnormalized partial derivatives, determine the slopes.

Exactly one of ell,m has 2-part two; the other has 2-part
`N>=4`. Let U be the coordinate along the latter period and
V along the former. The local torus algebra is

    A=kappa[U,V]/(U^N,V^2),       kappa=F4.             (3.2)

Odd period factors contribute units, so do not change this
description. In this algebra, smoothness in U implies

    P=unit*(U+a V),       Q=unit*(U+b V).               (3.3)

Indeed the root of either equation modulo `V^2` is its linear
slope times V, and division by the corresponding monic linear
polynomial has unit quotient. The slopes are

| longer 2-part | a for P | b for Q |
|---|---|---|
| x | beta | alpha^2 |
| y | beta^2 | alpha |

They are distinct at T and equal at G. This verifies the
conventions in either orientation without identifying the two
codes by a mirror symmetry.

Here and below a claim about associated local polynomials or
local units is algebraic, not a physical Hamming isometry.

Put `c=a+b`. Direct multiplication in (3.2) gives

    ann(P)=span_kappa{U^(N-1)+a U^(N-2)V, U^(N-1)V},
    ann(Q)=span_kappa{U^(N-1)+b U^(N-2)V, U^(N-1)V}.

Each has dimension two. They are the same space only when
a=b. Their common annihilator has dimension

    d=1 at T,       d=2 at G.                          (3.4)

Since the product is a unit times `U^2+c UV`, its annihilator
has the following four-element basis:

    U^(N-2)+c U^(N-3)V,
    U^(N-1),       U^(N-2)V,       U^(N-1)V.            (3.5)

For example, writing a general element as `A(U)+VB(U)`, the
ordinary part of its product forces `A(U)` to be divisible
by `U^(N-2)`; the V part then fixes just the coefficient at
`U^(N-3)` in B(U). The four displayed free coefficients are
exactly the resulting solutions. Hence

    dim ann(PQ)/(ann(P)+ann(Q))=4-(4-d)=d.              (3.6)

The N>=4 hypothesis matters. In the local model with N=2 and
a!=b, the product is `c UV`, its annihilator has dimension
three, and the quotient in (3.6) is zero instead. No extension
to arbitrary even-by-even shapes is asserted.

### General identification of the hidden quotient

For any finite commutative F2 algebra, the map in (1.4) has
kernel `ann(P)+ann(Q)`. Indeed `(Qz,0)=delta(h)` precisely
when `h in ann(P)` and `z+h in ann(Q)`.

Its image is exactly `J_first intersect J_second`: the same
class is represented by `(0,Pz)`, since their sum is delta(z).
Conversely, equality of a first-only and a second-only class
requires their difference to be delta(z), giving this form.
It is also exactly the image of both-principal cycles:
write one as `(Qr,Ps)` and subtract delta(s), obtaining
`(Q(r+s),0)` with `r+s in ann(PQ)`.

The first-only image is `ann(P)/Q ann(P)`. Its dimension is
the dimension of the kernel of multiplication by Q on ann(P),
namely the common annihilator dimension d. The same holds for
the second-only image. Equation (3.6) therefore proves that
both entire images coincide locally. Inactive components
contribute zero, and each ordinary residue field has F2 degree
two. Summing gives `2*1+2*2=6`, proving the coincidence and
dimension portions of (1.2) and (1.4).

## 4. Tor identifies the exact period axis

It remains to identify this common six-dimensional subspace,
rather than merely count it.

Use the completed local plane ring `O=kappa[[U,V]]`, and let
F be the actual period polynomial in the larger-2-adic
direction and G_period the other actual period polynomial.
Thus

    F=U^N*unit,       G_period=V^2*unit.                (4.1)

The completion is used only for a local finite-length algebra
calculation; it introduces no infinite-support physical chains.

The local plane intersection module `M=O/(P,Q)` has length
d from (3.4), and is killed by V^2. One direct check of this
plane multiplicity, not just its torus truncation, is to take
U=u,V=v first. Solving P=0 gives

    U=beta V+O(V^3).

Substituting in Q gives

    (alpha beta+1)V
      +(alpha^2 beta^2+1+alpha beta)V^2+O(V^3).

At T the linear coefficient is nonzero; at G it vanishes
and the quadratic coefficient is one. Both curve coordinates
are uniformizers by the nonzero slopes, so the same lengths
and the V^2 annihilation hold when the axis roles are reversed.

The period and plane pairs are regular sequences. Balanced Tor
therefore identifies local physical homology with

    H_1(M tensor K(F,G_period))=M direct-sum M.          (4.2)

Both period actions vanish here. Its coordinates are the actual
ordered period coordinates, with the larger 2-part listed first.

The first-block subspace is the image of

    Tor_1^O(O/(P), O/(F,G_period))
      -> Tor_1^O(M, O/(F,G_period)).                    (4.3)

In the plane resolutions this is exactly the chain map sending
the generator for P to the first Koszul coordinate. Thus (4.3)
is genuinely the image of first-only physical cycles.

On the smooth curve `C=O/(P)=kappa[[V]]`, U is V times a unit.
Consequently

    F=G_period h,       h in V^(N-2) C.

Since C is a domain, cycles of `C tensor K(F,G_period)` are
the multiples of `(1,h)`, and its boundaries are multiples
of `G_period(1,h)`. Its H1 is therefore C/(G_period), with
that generator. On passing to M, h vanishes because N-2>=2.
The induced map is

    C/(G_period) -> M direct-sum M,
    z -> (z mod M,0).                                 (4.4)

It is onto the first axis. The quotient to M is surjective
because G_period annihilates M. Repeating the same argument
on the smooth curve Q=0 proves the identical axis statement
for the second block.

Naturality in (4.3)-(4.4) follows by using the two augmentations
of the bounded tensor product of the plane and period free
resolutions. Both augmentations commute with the quotient
`O/(P) -> M` and its plane-resolution lift. Thus balanced Tor
does not replace the physical map by an unrelated isomorphism.

Returning to the fixed order `(x^ell+1,y^m+1)` gives (1.3).
The cylinder map wrapping one of these periods is the inclusion
of its period coordinate: the one-term cylinder resolution
maps to the corresponding column of the period Koszul complex.
This is exactly the W_L identification in (1.2), consistent
with `A40_PARITY_TOR_PROOF.md`. No weight-compatible Tor basis
or support lift of a particular representative is asserted.

## 5. A fully specified gross counterexample to the semisimple shortcut

In the fixed ring with `(ell,m)=(12,6)`, define

    z=(1+y+y^4+y^5)
      +x(y+y^2+y^3+y^4)
      +x^2(1+y^2+y^3+y^5)
      +x^3(y+y^2+y^4+y^5)
      +x^4(1+y+y^3+y^4)
      +x^5(1+y^2+y^3+y^5)
      +(x^6+x^7)(1+y^2).                              (5.1)

Finite polynomial multiplication modulo the two periods gives
`PQz=0`. The pair

    (f,g)=(Q(z+1),P)                                  (5.2)

is therefore a cycle with both components nonzero and both in
their principal stabilizer ideals. Their weights are 13 and
three, respectively. No minimality is claimed.

Here is a short nonboundary certificate independent of a
minimum-distance computation. At T use
`x=omega(1+U), y=omega(1+V)` and truncate by `U^4,V^2`.
The coefficient of `U^2 V^0` in z is one. In (5.1), the x^2
and x^3 rows evaluate to zero at y=omega; the x^6 and x^7
rows contribute omega and omega^2, whose sum is one. These
are exactly the rows selected by the U^2 binomial coefficient.
But at T, `ann(P)+ann(Q)` is spanned by
`U^3,U^2 V,U^3 V`, none of which has this coefficient.
Thus z is outside that annihilator sum. By (1.4), `(Qz,0)`
and hence (5.2) are nontrivial.

This invalidates the unstated extension of Theorem 90's
principal-membership criterion to the original even gross
ring. It is not a counterexample to that theorem's semisimple
statement, nor to the IBM distance conjecture.

## 6. Both exclusions and the class-matched lower bounds

For the torus ring R define

    C_first=(Q:P)/(Q),       C_second=(P:Q)/(P),
    (Q:P)={f: Pf in (Q)},    (P:Q)={g: Qg in (P)}.

There are natural maps

    phi_first:H -> C_first,     [(f,g)] -> f+(Q),
    phi_second:H -> C_second,   [(f,g)] -> g+(P).        (6.1)

They are well defined because a boundary changes the components
by Qh and Ph. They are surjective: for any f in (Q:P), choose
g with Pf=Qg; the resulting pair is a cycle, and conversely
for the second map.

The kernel of phi_first is J_second. Indeed, if f=Qh,
subtract delta(h) to eliminate the first block. Conversely
a second-only representative has zero image. Similarly,
the kernel of phi_second is J_first. Using the theorem,

    H/K = C_first = C_second                            (6.2)

via the respective canonical quotient isomorphisms, and these
spaces have F2 dimension six. Their identification matches
the two component cosets of the same quotient homology class.
It is not an identification of physical weights.

This proves (1.5). In particular, not only does every class
outside K force both exclusions, but every representative of
a class inside K has both components in the principal ideals.
There is no cycle with exactly one component outside its
corresponding principal ideal in this family.

In particular `ann(P) subset (Q)` and `ann(Q) subset (P)`.
Nevertheless `ann(P)/Q ann(P)` and its second-block counterpart
are nonzero and have dimension six. Membership in the full
principal ideal (Q) is not the single-block boundary test:
the correct denominator there is `Q ann(P)`. This distinction
preserves Theorem 88's single-block terms.

For a nonzero `eta in H/K`, let

    d_first(eta)=min{|f|: f+(Q)=phi_first(eta)},
    d_second(eta)=min{|g|: g+(P)=phi_second(eta)},        (6.3)

where (6.1) is understood through its induced maps on H/K.
Every representative of every full homology class lying above
eta satisfies the rigorous bound

    |(f,g)| >= d_first(eta)+d_second(eta).              (6.4)

Thus the minimum weight over all classes outside K is at least

    min_(eta!=0) [d_first(eta)+d_second(eta)]
      >= d_R((Q:P) minus (Q))
           +d_R((P:Q) minus (P)).                     (6.5)

The class-matched minimum can be stronger than independently
minimizing its two summands. Conversely, independently chosen
minimizers need not satisfy the joint cycle equation or lift
to a prescribed full homology class. Neither equality in
(6.4)-(6.5) nor a numerical value for either relative minimum
is claimed. In particular this note does not prove values
6r for the individual terms or 12r for their sum.

All S11 mixed winding classes lie outside K, so (6.4)-(6.5)
apply to them, and also to the other two cylinder axes.
The eliminations used to prove the kernels in (6.1) are only
algebraic: no claim that subtracting a boundary lowers weight
enters the bounds, which apply directly to the original two
components of the representative.

## 7. Bounded verification and remaining scope

The new standard-library diagnostic is
`scripts/a43_s4_codex_block_purity.py`. It constructs only the
fixed gross pair `(12,6)` and its fixed norm cover `(12,12)`
retaining the x period. All matrices are small binary bitset
maps; no period scan, support search, SAT, or distance census
is performed.

The gross calculation gives

    dim H=12,       rank delta=66,
    dim ann(P)=dim ann(Q)=12,
    dim ann(PQ)=24,       dim common ann=6,
    dim(J_first)=dim(J_second)=dim(J_first+J_second)=6,
    dim K=6.

The cover boundary rank is 138. The norm map on gross H has
rank six and vanishes on both single-block images. Since the
cover retains x, this independently checks the predicted x
axis for r=1. The diagnostic also checks all cycle equations,
the exact witness (5.1)-(5.2), its nonboundary status by binary
elimination, and its transverse U^2 coefficient.

The final fixed check passed in **0.002135 seconds**, with
measured peak RSS **16,449,536 bytes**. The first rank pass
took 0.002118 seconds; total compute is far below the allowed
aggregate. Each run has a 55-second wall alarm, a 50-second
CPU limit, and a 768-MiB measured-RSS guard. The universal
theorem comes from §§2–4 and §6, not from these two fixed
matrices.

No hypothesis remains open in the stated class-space theorem.
What remains is metric: establishing useful lower bounds for
the relative cosets in (6.3), and controlling representatives
in K itself. This does not prove footprint connectedness of
class-minimal representatives, a cylinder support lift of a
given minimum representative, or the original IBM distance
law. Only the new note and its new diagnostic were written;
no shared file, git state, Lean file, memory file, or old
process was touched. No additional assignment is initiated.
