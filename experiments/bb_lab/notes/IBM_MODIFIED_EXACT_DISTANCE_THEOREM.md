# Exact distance for an explicit one-check IBM modification

2026-09-08. Self-contained analytic theorem, combining the one-check
algebra, the parent row-seam reduction, and the sharp cylinder proof
in `IBM_ONE_CHECK_CYLINDER_TANGENT.md`. No distance enumeration or
unproved support-normalization lemma is used.

## Theorem and the change from IBM

In the original program's cleared Laurent coordinates the checks are

    P=x^3+y+y^2,       Q_IBM=1+(x+x^2)y^3.

Keep P and replace only the second check by

    chi=1+x+x^2,       Q=chi(y+x+1)
                       =1+x^3+y+xy+x^2y.              (1)

On Z_ell x Z_m, for any positive multiples of six ell,m with ell<=m,
use H_X=[M_P | M_Q], H_Z=[M_Q^T | M_P^T]. This binary CSS code has
the exact parameters

    [[2 ell m, 12, 2 ell/3]].                         (2)

In particular, on the original square frames ell=m=6r,

    [[72r^2, 12, 4r]]     for every integer r>=1.      (3)

Both sector distances have this exact value. The code is connected,
has two physical qubits per cell, and has check weight eight instead
of the original six. All periods, including multiples of 127, are
covered. The original ordinary primary algebra, including its doubled
tangent component, is retained by (1); the exceptional component is
removed. These algebraic facts do not constitute a metric equivalence.

This is a full exact-distance theorem for the stated modification.
It is **not** the unchanged IBM conjecture d=6(2r+b-1), and does not
retain that conjecture's distance coefficient. No claim is made for
ell>m, for the original checks, or for novelty of this construction.

## 1. Plane algebra, dimension, and the two sectors

All polynomials are over F2. Put R=F2[x^+/-1,y^+/-1],
U=y+x, V=y+x+1, and M=R/(P,Q). Exact identities are

    P=UV+x chi,
    chi P+U Q=x chi^2.                               (4)

Modulo P, R is free of rank two over L=F2[x^+/-1], with basis 1,y.
Indeed P is monic in y and y(y+1)=x^3 makes y invertible there.
Multiplication by Q has matrix

    chi [ x+1   x^3 ]
        [  1     x  ].

Its determinant is x chi^3, nonzero in L, so this multiplication is
injective. Its Smith factors are chi and chi^2: the entry gcd is chi
and their product is the determinant up to the unit x. Therefore

    M = L/(chi) direct-sum L/(chi^2) as an L-module,
    dim_F2 M=6.                                      (5)

The injectivity also proves plane exactness directly. If Pf+Qg=0
in R, reduction modulo P first gives g=Ph. Since R is a domain,
the original equation then gives f=Qh. Thus every finite Laurent
plane cycle is a finite Laurent boundary (Qh,Ph).

Both x^6+1 and y^6+1 annihilate M. This follows from (4) and

    x^6+1=(x+1)^2 chi^2,
    (1+y+y^2)^2=P^2+(x^6+1),
    y^6+1=(y+1)^2(1+y+y^2)^2.                        (6)

Consequently F=x^ell+1 and G=y^m+1 annihilate M. On the finite
torus T=R/(F,G), the quotient T/(P,Q) is M and has dimension six.
For N=ell m, the syndrome map D=[P | Q] has rank N-6. Its boundary
map delta=[Q;P] has the same rank: coordinate inversion transposes
each circulant, and block swap then identifies the relevant ranks.
Since P,Q commute, D delta=0. Hence the CSS logical dimension is

    2N-rank D-rank delta=12.                         (7)

All eight check sites are distinct because ell,m>=6 and the two
blocks are disjoint. The support difference between y and y^2 in
P gives a unit vertical move among syndrome sites; the difference
between y and xy in Q gives a unit horizontal move. These generate
the torus, and every qubit is incident to a syndrome site. Thus
the incidence graph is connected, with bounded check weight eight.

The physical permutation (f,g)->(bar(g),bar(f)), where bars reverse
both lattice coordinates, identifies the two sector cycle and
boundary spaces. It therefore suffices to prove the minimum for
D(f,g)=Pf+Qg, with boundaries delta(h)=(Qh,Ph).

For the algebraic comparison, the exact identity

    Q_IBM=(x+x^2)(y+1)P+x^2 Q+chi^2(1+xy)

together with (4) gives (P,Q)=(P,Q_IBM,chi^2). The original Smith
calculation identifies this quotient with its whole ordinary primary
part, not just its reduced points; see `IBM_ONE_CHECK_MODIFICATION.md`
section 2.2. This comparison is not needed for the distance proof.

## 2. The classical recurrence

At every length L divisible by three, the binary cyclic kernel of
1+z+z^2 has dimension two and its nonzero elements repeat one of
110,101,011. Its coefficient recurrence forces period three; every
nonzero word therefore has exactly 2L/3 occupied sites.

Over F4 the same recurrence repeats triples (a,b,a+b). A nonzero
triple has at least two nonzero entries, so every nonzero word has
at least 2L/3 occupied sites. These statements include even L and
do not assume a semisimple ambient group algebra.

## 3. First cut: a light torus cycle lifts to a cylinder

Let z=(f,g) be any torus cycle of weight w<2 ell/3. Fix omega in
F4 with omega^2+omega+1=0. Evaluation x=omega is well-defined,
since 3 divides ell. At that value P=1+y+y^2 and Q=0. Hence

    (1+y+y^2) f(omega,y)=0 modulo y^m+1.

Evaluation cannot create an occupied row from an empty physical
row, so |f(omega,y)|<=|f|<=w<2 ell/3<=2m/3. Section 2 implies
f(omega,y)=0. Every binary row of f thus vanishes at omega, so
each nonempty f row has at least two physical sites.

Let R_f,R_g count occupied rows in the two blocks. Then

    2R_f+R_g <= |f|+|g|=w<m.                         (8)

An f qubit in row j has syndrome neighbors only in rows j,j+1,j+2,
and can cross just the two row seams j+1,j+2. A g qubit has
neighbors only in rows j,j+1, crossing at most seam j+1. All
sites anchored in a given row have those same possible crossed
seams, irrespective of their x positions. By (8), fewer than m
seams are crossed. Cut at an uncrossed seam.

Every occupied site and all its syndrome neighbors then lie in
one m-row interval of the y cover, without wrapping through the
cut. Distinct sites stay distinct and all syndrome cancellations
are preserved. This is an injective, weight-preserving lift to
the x-periodic, y-unbounded cylinder. The entire given vector
lifts; no assumption of connectivity or minimality was used.

## 4. Second cut: a light cylinder cycle lifts to the plane

Apply the lattice shear x=X, y=XY and translate the first block:

    f'=X f(X,XY),       g'=g(X,XY).

It maps an f site (i,j) to (i+j+1,j) and a g site to (i+j,j).
These are bijections on Z_ell x Z. They preserve physical weight,
finite row support, and the sole period X^ell+1. The new checks are

    P'=X^2+Y+XY^2,
    Q'=1+X^3+(X+X^2+X^3)Y.                          (9)

Their syndrome is exactly the sheared original syndrome. The
boundary (Qh,Ph) becomes (Q'H,P'H) with H=Xh(X,XY), so boundaries
are preserved as well. This is a physical permutation, not an
uncontrolled algebraic change of qubit basis.

Both checks in (9) become chi(X) at Y=1. The contracting invariant

    c(X)=f'(X,1)+g'(X,1)

therefore satisfies chi c=0 on the length-ell circle, and |c|<=w.
Since w<2 ell/3, Section 2 forces c=0. Thus each occupied X column
contains an even total number of occupied qubits, counting both
blocks and all Y rows. If K counts these columns, then 2K<=w.

All syndrome offsets of P' in X lie in {0,1,2}; those of Q' lie
in {0,1,2,3}. For every occupied column a, all qubits anchored
there together can cross only the three seams a+1,a+2,a+3.
The number of forbidden column seams is therefore at most

    3K <= 3w/2 < ell.                               (10)

Cut at an uncrossed column seam. Because all offsets in (9) are
nonnegative and at most three, every occupied qubit and every
syndrome neighbor lies in the chosen ell-column interval without
wrap. This gives an injective, weight-preserving finite plane
cycle. There was no Y period on this cylinder to undo.

Plane exactness from section 1, transported through the physical
shear, makes the lifted cycle a finite plane boundary. Project its
generator to the cylinder, undo the shear, and project back to the
original torus. The original z is a boundary. This proves, for all
classes in both sectors,

    d >= 2 ell/3.                                   (11)

The invariant c alone does not detect the doubled component.
Its kernel is handled by even column occupancy and the second
cut; no logical classes have been discarded.

## 5. A matching nontrivial physical stripe

Define the finite polynomial

    W=(1+x) sum_(j=0)^(ell/3-1) x^(3j).

It has exactly 2 ell/3 distinct sites modulo ell and satisfies
chi W=F. Consequently z_*=(0,W) is a torus cycle, with its finite
Laurent lift obeying

    D z_*=V F.                                      (12)

Here is a direct proof that it is not a torus boundary. Otherwise,
for finite Laurent lifts one could write

    z_*=delta(h)+F a+G b

with two-component vectors a,b. Applying D to this equality and
using (12) gives

    F(V-Da)=G Db.

The polynomials F and G are coprime in the Laurent domain R,
since they depend on different independent variables. Thus G
divides V-Da. In M, both Da and G vanish, which would give V=0.
But evaluation x=omega,y=omega annihilates P,Q and sends V to one.
This is a contradiction. Hence z_* is nontrivial and

    d <= |z_*|=2 ell/3.                             (13)

Equations (11) and (13), together with (7), prove (2) and (3).
The matching upper and lower bounds hold at every parameter, not
only along a numerically tested sequence.

## Verification and provenance

The one-check algebra and first seam argument had independent
read-only reviews before the sharp cylinder result. The cylinder
proof and its composition have also passed two independent read-only
reviews, including the actual period lattice, both physical cuts,
the doubled component, repeated-root dimension, and nonboundary
stripe. A separate review passed the self-contained integration.
The parent has replayed the prescribed polynomial/shear/stripe
diagnostic successfully in under 0.001 seconds and 20 MB RSS.
It checks identities, not distances; the universal proof is above.
No Lean formalization or computational distance census is asserted.

All three native agents have completed this proof-and-review pass.
The permitted closely related modified-family objective is achieved;
the original IBM formula remains open. These files are local and
prepared for the PR #27 research publication update. No old paused
headless session was resumed.

Detailed sources within this workbench are
`IBM_ONE_CHECK_MODIFICATION.md`,
`IBM_ONE_CHECK_CYLINDER_REDUCTION.md`, and
`IBM_ONE_CHECK_CYLINDER_TANGENT.md`.
