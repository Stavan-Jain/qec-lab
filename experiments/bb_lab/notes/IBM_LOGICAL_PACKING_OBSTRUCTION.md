# Why disjoint or fractionally packed logical witnesses cannot prove the IBM target

2026-09-08. A general analytic limitation on a possible topological
proof technique, not a distance counterexample or a literature claim.
No computation is needed for this argument.

## 1. The certificate under consideration

Consider a binary CSS code on n physical qubits. Fix a nonzero logical
class alpha in one sector. Choose logical representatives z_i in the
opposite sector, each pairing to one with alpha. Give them nonnegative
real weights lambda_i satisfying the coordinate capacity constraints

    sum_(i: j in supp(z_i)) lambda_i <= 1, for every qubit j.

Write L_alpha=sum_i lambda_i. For every representative v of alpha,
each intersection supp(v) intersect supp(z_i) has odd cardinality,
in particular at least one. Therefore

    |v| >= sum_i lambda_i |supp(v) intersect supp(z_i)|
        >= L_alpha.

This is a valid class-specific lower-bound certificate. Disjoint
logical witnesses are its special case with all lambda_i=1.

However, if the opposite-sector distance is d_op, every z_i has
weight at least d_op. Summing the coordinate capacities gives

    d_op L_alpha <= sum_i lambda_i |z_i| <= n.          (1)

Thus no such certificate can exceed n/d_op, even with arbitrarily
many witnesses, fractional weights, or a different packing for each
source class. Allowing signed weights would not preserve the displayed
lower-bound argument and is outside this certificate model.

## 2. The two BB sector distances agree without an axis rotation

Let A and B denote the commuting translation matrices of any binary
bivariate bicycle pair on a finite rectangular torus. In the lab
convention,

    H_X=[A | B], H_Z=[B^T | A^T].

Let R reverse both lattice coordinates. It is a permutation involution
with R A R=A^T and R B R=B^T. Define the physical permutation

    Pi(f,g)=(R g,R f).

Then H_Z Pi(f,g)=R H_X(f,g). It maps the corresponding boundaries
(B h,A h) to (A^T R h,B^T R h), so it induces a weight-preserving
isomorphism between the two logical quotient spaces. Consequently

    d_X=d_Z=d.

This exchanges CSS sectors; it is not the unproved assertion that
the original and axis-rotated cylinder pairs have the same metric.
Neither semisimplicity nor a fixed logical dimension is needed.

If a collection of the certificates in section 1 proves a uniform
floor L for every nonzero class in one sector, it proves d>=L.
Equation (1) applied to any one of those packings now gives

    L^2 <= d L <= n.                                  (2)

Hence the raw, uniform packing value cannot exceed sqrt(n).
Integer rounding matters: d>=L implies d>=ceil(L), so the general
integer ceiling from this argument is ceil(sqrt(n)), not sqrt(n).
This is a limitation of the proof mechanism, not a general upper
bound on quantum-code distance.

More precisely, to prove an integer floor D by rounding, every
class must have a packing value strictly greater than D-1.
Then d_op=d>=D, and (1) requires D(D-1)<n. Thus the method cannot
certify D whenever D(D-1)>=n, even with integer rounding.

## 3. Application to the requested parameters

For b=1, n=72r(r+1) and the proposed exact distance is D=12r.
The rounding-aware difference is

    D(D-1)-n=12r(6r-7)>0 for every integer r>=2.

For b=0, n=72r^2 and D=12r-6. Here

    D(D-1)-n=6(12r^2-26r+7)>0 for every integer r>=2.

For the last positivity assertion, the quadratic equals
3+(r-2)(12r-2), which is positive throughout that range.

Thus proving the complete IBM target at these parameters requires
more than a packing of opposite logical witnesses, even a fractional,
class-dependent one. A gross-only packing at r=1 cannot simply be
scaled to furnish the desired all-r certificate.

This rules out only the specified support-packing argument. Joint
syndrome identities, nonlinear intersection constraints, branched
network arguments, or a weight-sensitive chain equivalence could use
additional information absent from it. No claim is made that those
methods succeed; their metric constraints must be proved separately.
