# One changed monomial: an exact growing-distance BB family

2026-09-08. Parent analytic proof. No distance search, width cutoff,
asymptotic estimate, or unproved metric normal form is used.

## 1. The precise modification and theorem

The original IBM Laurent checks in this program are

    A=1+y+x^3 y^-1,       B=1+x+x^-1 y^-3.

Keep A unchanged and replace just the final monomial of B:

    A_flat=1+y+x^3 y^-1,  B_flat=1+x+x^-1.

Thus five of the six original monomials are unchanged. The changed
term is x^-1 y^-3 -> x^-1. Both checks are still trinomials, with
two physical qubits per lattice cell and check weight six.
Independent monomial translations of the two qubit blocks clear
the checks to

    P=y A_flat=x^3+y+y^2,       Q=x B_flat=chi_x,
    chi_x=1+x+x^2,              chi_y=1+y+y^2.

Multiplying the first and second checks by monomials a,b is a
physical equivalence: send a cycle (f,g) to (a^-1 f,b^-1 g).
The boundary (B h,A h) maps to (b B h',a A h') with
h'=(ab)^-1 h. All these changes are coordinate permutations.

**Theorem.** For positive multiples of three ell,m with ell<=m,
the modified BB code on Z_ell x Z_m has exact parameters

    [[2 ell m, 8, 2 ell/3]].                         (1)

Both CSS sector distances equal 2 ell/3. In particular, on the
original square frames ell=m=6r, for every integer r>=1,

    [[72r^2, 8, 4r]].                                (2)

Both periods and the exact distance grow with r. The code is
connected and has weight-six checks throughout the family.

This proves an exact theorem for the explicitly modified checks,
not d=6(2r+b-1) for the unchanged IBM pair. It changes the ordinary
logical dimension from twelve to eight and loses IBM's conjectured
distance coefficient. It does not preserve the doubled tangent
component of the different five-term modification Q=chi_x(y+x+1)
in `IBM_ONE_CHECK_MODIFICATION.md`. No lower bound is transferred
from that five-term code or back to original IBM.

## 2. Commutation, logical dimension, and physical structure

Let T=F2[x^+/-1,y^+/-1]/(x^ell+1,y^m+1). A sector has

    D(f,g)=P f+chi_x g,       delta(h)=(chi_x h,P h).

Commutativity and characteristic two give D delta=0. With
H_X=[M_P | M_chi_x] and H_Z=[M_chi_x^T | M_P^T], these are
the usual BB CSS cycle and boundary spaces.

The exact finite identity

    P=chi_y+(1+x)chi_x                               (3)

implies (P,chi_x)=(chi_y,chi_x). Since ell and m are divisible
by three, both period polynomials vanish modulo these two
quadratics. Consequently

    T/(P,chi_x)=F2[x,y]/(chi_x,chi_y)

has binary dimension four; x and y are units in this quotient.
Writing N=ell m, rank D=N-4. Coordinate inversion takes each
circulant matrix to its transpose; together with block swap
this gives rank delta=rank D=N-4. Thus

    dim H_1 = 2N-(N-4)-(N-4)=8.                     (4)

This argument includes repeated roots when either period is
even; no semisimplicity assumption on the finite torus is used.

The six sites of every check are distinct for ell,m>=3: each
block contains three different monomials, and the blocks are
disjoint physical coordinate sets. The difference between the
x and 1 terms of chi_x connects neighboring syndrome columns;
the difference between the y and y^2 terms of P connects
neighboring syndrome rows. These generate the full torus, and
every qubit is incident to a syndrome site. The incidence graph
is connected, rather than a disjoint union of fixed-size codes.

Inversion of both coordinates followed by block swap maps
(f,g) to (bar(g),bar(f)), and identifies the two CSS cycle and
boundary spaces by a physical permutation. It suffices to
prove the distance for D,delta.

## 3. A classical recurrence used in both directions

For any multiple-of-three length L, the kernel of 1+z+z^2
on F2[z]/(z^L+1) has dimension two. Its coefficient recurrence
forces period three. The three nonzero patterns are repetitions
of 110,101,011, so every nonzero kernel vector has weight 2L/3.

Over F4 the same recurrence gives triples (a,b,a+b). Every
nonzero triple has at least two nonzero entries. Thus every
nonzero F4-valued kernel vector has weight at least 2L/3.
These are elementary recurrence statements for all lengths in
question, not conclusions from finite enumerations.

## 4. Two contracting maps determine every cylinder class

Here x has period ell and y has finite Laurent support, with

    R_ell=F2[x^+/-1]/(x^ell+1),       C=R_ell[y^+/-1].

For s=1,2 define

    E_s(f,g)=x^s f(x,x^s)+g(x,x^s) in R_ell.          (5)

Every input monomial maps to just one output monomial. Therefore

    |E_s(f,g)| <= |f|+|g|                            (6)

in the actual binary physical coefficient bases, even if terms
collide and cancel. The finite identities

    P(x,x)=x chi_x,       P(x,x^2)=x^2 chi_x           (7)

show that a cylinder cycle maps into ker chi_x. They also show
that E_s vanishes on every boundary: the two boundary terms in
(5) are identical and cancel.

It remains essential to prove these maps jointly detect all
cylinder homology, not just selected primary components.
Modulo P, the Laurent ring C has free R_ell basis 1,y. Indeed
P is monic of degree two in y, with unit constant coefficient
x^3; the relation y(y+1)=x^3 makes y invertible in the quotient.
Thus every g has a unique remainder u+yv, u,v in R_ell, and
a finite Laurent decomposition

    g=P h+u+yv.                                     (8)

Reduce the cycle equation modulo P. Since 1,y are a free basis,

    chi_x u=chi_x v=0.                              (9)

The remainder therefore has zero product with chi_x already in
C, not merely modulo P. Substituting (8) into D(f,g)=0 gives
P(f+chi_x h)=0. Multiplication by P is injective on C: for a
nonzero finite Laurent polynomial its highest y coefficient
survives multiplication by P's leading coefficient one.
Hence f=chi_x h. Subtracting the actual boundary delta(h)
puts the class in the unique form

    (0,u+yv),       u,v in ker chi_x.                (10)

Uniqueness follows from the uniqueness of the remainder in
(8). This algebraic reduction is not asserted weight-preserving;
only the contracting maps (5) are used for the distance bound.

On (10), the two invariants are

    E_1=u+xv,       E_2=u+x^2 v.

If both vanish, their sum is (x+x^2)v=v, by chi_x v=0.
Thus v=0 and then u=0. Conversely zero homology has both
invariants zero because they vanish on boundaries. This proves
joint injectivity on the entire cylinder homology.

Every nontrivial cylinder cycle therefore has some nonzero E_s.
Sections 3 and (6) prove

    |(f,g)| >= 2 ell/3.                              (11)

This is an all-class physical lower bound, with no bound on
the Laurent extent in y.

## 5. Low-weight torus cycles have an exact cylinder lift

Suppose z=(f,g) is a torus cycle with

    |z| < min(2 ell/3,2m/3).                         (12)

Choose omega in F4 with omega^2+omega+1=0. Since 3 divides ell,
evaluation x=omega is well-defined on the physical torus.
At that value P=chi_y and chi_x=0, so

    chi_y f(omega,y)=0 modulo y^m+1.

Evaluation cannot create a nonzero row from an empty physical
row, so |f(omega,y)|<=|f|<2m/3. Section 3 forces this evaluated
vector to vanish. Every binary row of f is consequently
divisible by chi_x, and each nonempty f row has at least two
qubits: a single monomial cannot vanish at omega.

If R_f is the number of occupied f rows, then 2R_f<=|f|<m.
An f qubit has syndrome rows j,j+1,j+2 and can cross only the
two row seams j+1,j+2. All f qubits in the same row cross the
same possible seams. A g qubit has all syndrome neighbors in
its own row, since chi_x has no y displacement; it crosses no
row seam. At most 2R_f<m seams are crossed. Cut at an uncrossed
seam and choose the associated m-row interval of the cover.

Every occupied qubit and every syndrome neighbor stays in
that interval without wrap. Distinct physical sites stay
distinct and all syndrome cancellations are unchanged. This
is an injective, weight-preserving lift of the entire cycle
to the x-periodic, y-Laurent cylinder. No connectivity or
minimality assumption about z was needed.

By (11) and (12), the lifted cylinder cycle must be a boundary.
Project its finite Laurent boundary generator back to the
torus. The original z is then a boundary as well. We have proved
for all multiples-of-three periods, without an aspect assumption,

    d_tor(ell,m) >= min(2 ell/3,2m/3).                (13)

## 6. A physical logical attains the floor

The finite binary polynomial

    W_ell=(1+x) sum_(j=0)^(ell/3-1) x^(3j)

has exactly 2 ell/3 distinct sites and satisfies

    chi_x W_ell=x^ell+1.

Therefore z_ell=(0,W_ell), in one physical row, is a torus
cycle. It is nontrivial for every multiple-of-three m, as the
following direct boundary obstruction shows.

If (chi_x h,P h)=(0,W_ell) for a torus polynomial h, then
chi_x h=0 and (3) gives P h=chi_y h. Evaluation y=omega is
well-defined because 3 divides m. It sends chi_y h to zero
and W_ell to the same nonzero polynomial in
F4[x]/(x^ell+1). This is a contradiction. Thus

    d_tor(ell,m) <= |z_ell|=2 ell/3.                 (14)

For ell<=m, (13) and (14) agree and prove (1), in both CSS
sectors. Substituting ell=m=6r proves (2). There are no missing
small r cases or exceptional prime periods.

## 7. Relation to product codes and limits of the claim

Identity (3) gives the invertible linear shear

    (f,g) -> (f,g+(1+x)f),

which takes this complex to the separable pair (chi_y,chi_x)
and takes boundaries to its boundaries. This is an algebraic
description of why the modified module is simpler. It is not
a coordinate permutation and is not used to transfer the
separable code's Hamming distance. Such a transfer would leave
a metric gap. Sections 4-6 instead prove the distance in the
original modified physical basis by explicit contractions,
an exact geometric lift, and a matching logical vector.

This family is consequently a direct one-monomial deformation
of the specified IBM checks, but it is product-code-like rather
than a solution of IBM's nonreduced mixed-class problem. Its
qubit-to-distance-squared ratio is 9/2 on the square frames,
and its logical dimension is eight rather than twelve. No
novelty or optimality claim is made. The unchanged IBM conjecture
remains open. The separate k=12 five-term modification is now
proved in `IBM_MODIFIED_EXACT_DISTANCE_THEOREM.md`; no distance
transfer between those different checks is asserted.

The companion diagnostic checks prescribed polynomial identities
and the prescribed length-six stripe only. The universal claims
are the analytic proofs above, not numerical distance evidence.

Two independent read-only reviews passed the physical block
translations, joint cylinder-class detection, repeated-root
dimension, seam lifting, both-sector equivalence, and nontrivial
stripe. The parent fixed diagnostic passed in 0.000020 seconds
with 15.80 MiB peak RSS. No full distance computation was used.
