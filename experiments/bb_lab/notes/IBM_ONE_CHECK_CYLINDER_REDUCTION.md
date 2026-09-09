# The one-check modification: exact square-to-cylinder distance reduction

2026-09-08. Parent analytic continuation. This removes the toroidal
branch for the explicit modified square family; it does not yet
compute the remaining cylinder distance. No search is used.

> **Subsequent closure.** The remaining cylinder floor in this
> reduction is now proved in `IBM_ONE_CHECK_CYLINDER_TANGENT.md`.
> The self-contained full exact theorem, including both cuts and
> the matching stripe, is `IBM_MODIFIED_EXACT_DISTANCE_THEOREM.md`.
> The open-residual language below describes the earlier reduction
> stage, not the current status of this modified family.

## 1. Definitions and hypotheses

Keep the original cleared polynomial P and change only Q:

    P=x^3+y+y^2,
    Q_new=(1+x+x^2)(y+x+1)=x^3+1+y+xy+x^2y.

On a rectangle of periods ell,m, both divisible by six, the cycle
equation is D(f,g)=P f+Q_new g=0 and boundaries are (Q_new h,P h).
All weights below count physical binary coefficients in both blocks.

Write d_tor(ell,m) for the nontrivial minimum on the finite torus.
Write d_cyl(ell) for the nontrivial minimum with only x periodic,
and with finite Laurent support in y. It is the x-periodic cylinder,
not the original A42 pure-omega cylinder.

The proved plane-module calculation in `IBM_ONE_CHECK_MODIFICATION.md`
gives the regular Koszul resolution, a six-dimensional plane module M,
and annihilators x^6+1 and y^6+1. The same note gives a nontrivial
one-row stripe of weight 2 ell/3. Its cylinder nontriviality is
checked in section 4 below.

## 2. An empty row seam below 2m/3

**Lemma 2.1.** Every torus cycle of weight strictly less than 2m/3
has an injective, weight-preserving lift to the x-periodic cylinder.
No connectivity or representative-minimality assumption is needed.

Proof. Choose omega in F4 with omega^2+omega+1=0. Evaluation at
x=omega is well-defined on physical x-periodic polynomials because
3 divides ell. Since omega^3=1,

    P(omega,y)=1+y+y^2,       Q_new(omega,y)=0.

Consequently F(y)=f(omega,y) satisfies

    (1+y+y^2) F(y)=0 modulo y^m+1.                    (1)

The coefficient recurrence in (1) forces F_j=F_(j-3). Each
period-three triple has the form (a,b,a+b). If it is nonzero,
at least two of its entries are nonzero. Because 3 divides m,
every nonzero solution of (1) has F4 Hamming weight at least
2m/3. On the other hand, evaluation can only discard occupied
rows, so |F|<=|f|. The assumed weight bound implies F=0.

Thus each physical row of f evaluates to zero at omega. Its
binary polynomial is divisible by chi=1+x+x^2, or is zero.
In particular each occupied f row has at least two physical
qubits: a single nonzero binary monomial cannot vanish at omega.
Let R_f,R_g count the occupied rows in the two blocks. Then

    2 R_f+R_g <= |f|+|g|.                             (2)

Label the m row seams by the boundaries between consecutive cyclic
syndrome rows. An f qubit in row j has syndrome rows j,j+1,j+2,
and crosses only the two seams j+1,j+2. All f qubits in that row
have the same two possible crossed seams, regardless of x. A g
qubit has syndrome rows j,j+1 and crosses only seam j+1. This
uses the actual nonnegative y offsets of P and Q_new, not just
an abstract locality radius.

The union of crossed seams therefore has cardinality at most
2 R_f+R_g, which by (2) is less than 2m/3 and hence less than m.
There exists a seam crossed by no occupied qubit's syndrome span.
Cut the torus there. Every occupied qubit and every one of its
syndrome neighbors lies in the same chosen m-row interval, with
no wrap through the cut. All incidences and their cancellations
are therefore unchanged on the cylinder. Distinct physical sites
remain distinct; the lift has finite y support and the same weight.
Its syndrome is zero, proving the lemma. QED.

The actual sufficient hypotheses used were |f|<2m/3 and
|f|+|g|<m. The stated symmetric total-weight threshold implies
both. A full-column evaluation was not assumed to preserve weight
or logical information; it was used only to prove the row-support
constraint needed for the geometric cut.

## 3. Exact class-image and metric consequences

The period-Tor description for this modified pair gives

    H_cyl(ell)=M,
    H_tor(ell,m)=M direct-sum M.

Indeed x^ell+1 and y^m+1 annihilate M. The actual projection from
the cylinder to the torus has homology image

    W_x=M direct-sum {0},                             (3)

and is injective there: in a finite plane lift the syndrome of
a cylinder cycle has the form (x^ell+1) A, with no y-period term.
The comparison with the period resolution sends it to ([A],0).
This is a statement about classes, not an isometry of arbitrary
representatives; projection may cancel physical coefficients.

By Lemma 2.1, every class outside W_x has weight at least 2m/3.
There are 4,032 nonzero such binary classes, since the full
logical space has dimension twelve and W_x has dimension six.
In particular all the winding-mixed classes obey this floor
for the modified code. The claim does not concern original Q.

If a torus cycle of weight below both 2m/3 and d_cyl(ell) were
nontrivial, Lemma 2.1 would lift it to a cylinder cycle of the
same weight. That lift cannot be a boundary, since projection
would then make the original vector a torus boundary. This
contradicts the definition of d_cyl. Conversely, the projection
of a nontrivial cylinder minimizer is still nontrivial by (3)
and has no larger weight. Thus

    min(2m/3,d_cyl(ell)) <= d_tor(ell,m) <= d_cyl(ell). (4)

These statements apply to both CSS sectors: inversion of both
lattice coordinates followed by block swap is a physical
permutation taking cycles and boundaries between the sectors.

## 4. A stripe attainer and equality on every square

Let

    W_ell=(1+x) sum_(j=0)^(ell/3-1) x^(3j).

The finite identity chi W_ell=x^ell+1 yields

    D(0,W_ell)=(y+x+1)(x^ell+1).

It is a cylinder cycle of weight 2 ell/3. Its cylinder period
coordinate is [y+x+1] in M. This is one in the simple component
where y=x and chi=0, so it is nonzero. Consequently

    d_cyl(ell) <= 2 ell/3.                            (5)

For ell<=m, (5) and (4) now prove the exact equality

    d_tor(ell,m)=d_cyl(ell).                          (6)

In particular, for every integer r>=1,

    d_tor(6r,6r)=d_cyl(6r).                           (7)

Thus the candidate exact growing theorem d=4r for the modified
square family is equivalent to the single cylinder-floor lemma

    d_cyl(6r)>=4r for every r>=1.                     (8)

No separate toroidal mixed-class estimate remains necessary for
this modified square family. Equality (7) is already proved;
the numerical evaluation (8) is not proved here. This distinction
prevents the reduction from being mistaken for the desired exact
distance theorem.

The existing unconditional lower bound d_cyl(6r)>=2r follows
from the same incidence-lifting argument used for the torus:
a small finite cylinder component with fewer than ell/3 qubits
has no nonzero x-period displacement and hence lifts to an exact
plane boundary. The current interval is therefore still

    2r <= d_cyl(6r)=d_tor(6r,6r) <= 4r.

All assertions in this note retain the specified new check,
physical support, and finite-lift hypotheses. They neither
transfer a lower bound back to original IBM distance nor
replace the remaining cylinder-floor proof with a computation.

Independent read-only review passed the row-divisibility and
seam argument, the actual class-injective projection, and the
equality at the threshold ell=m. The one-check algebra and
stripe certificates also passed that review. The remaining
numerical cylinder floor was open at the reduction stage; the
subsequent sharp cylinder proof now closes it as indicated above.
