# The five-term modification: a sharp all-class cylinder distance

2026-09-08. Bounded native Codex analytic pass. The proof uses a physical
lattice shear, a contracting classical invariant, and an empty seam. It
does not minimize algebraic normal forms or enumerate supports. In
particular, it covers the doubled tangent component, not only the simple
component detected directly by the invariant.

## 1. Theorem and precise scope

Let r be any positive integer, let ell=6r, and set

```
R=F2[x^+/-1,y^+/-1],       F=x^ell+1,
chi=1+x+x^2,              U=y+x,       V=y+x+1,
P=x^3+y+y^2,             Q=chi V.
```

On the cylinder C=R/(F), y has finite Laurent support and no period. Use
the actual physical binary coefficient basis in both blocks, with

```
D(f,g)=Pf+Qg,       delta(h)=(Qh,Ph),
wt(f,g)=|f|+|g|.
```

**Theorem.** Every cylinder cycle of weight less than 2ell/3 is a
boundary. A nontrivial cylinder cycle of weight 2ell/3 exists. Thus

```
d_cyl(6r)=4r                     for every integer r>=1.       (1.1)
```

The lower bound applies to every nonzero homology class, including the
four-dimensional doubled component. There is no dyadic restriction,
bound on the y extent, connectivity assumption, or minimizer assumption.

This is the five-term modification Q=chi(y+x+1), with check weight eight.
It is not the different trinomial modification Q=chi, and it does not
prove a distance law for the unchanged IBM pair.

## 2. Plane regularity and the class space

The identities

```
P=UV+x chi,
chi P+U Q=x chi^2,
V P+x Q=UV^2                                             (2.1)
```

hold in the Laurent domain R. Modulo P, the elements 1,y form a free
L=F2[x^+/-1] basis: P is monic in y, and y(y+1)=x^3 makes y invertible.
Multiplication by Q in this basis has matrix

```
chi [ x+1   x^3 ]
    [  1     x  ],
```

whose determinant is x chi^3, nonzero in L. Hence Q is a
non-zero-divisor modulo P. In particular, every finite Laurent plane
cycle is a finite Laurent boundary. Explicitly, Pf+Qg=0 first implies
g=Ph; then the domain property gives f=Qh.

For completeness, M=R/(P,Q) splits as

```
M = M_simple direct-product M_double,
M_simple = R/(chi,U) = F4,
M_double = R/(chi^2,V+x chi) = F2[x^+/-1]/chi^2.           (2.2)
```

Indeed (2.1) gives chi^2=UV^2=0 in M. Since U+V=1, the ideals
(U) and (V^2) are comaximal and their product vanishes. On the U=0
factor, V=1 and chi=0. On the V^2=0 factor, P gives V=x chi;
conversely chi^2=0 and V=x chi imply both P=Q=0. The resulting
value y=x+1+x chi is a unit, as its reduction modulo chi is nonzero.
This proves (2.2), including the actual tangent relation rather than
the false substitution V=0 on the doubled factor.

Thus dim_F2 M=2+4=6. Since

```
x^6+1=(x+1)^2 chi^2,
```

F annihilates M. A cylinder cycle lifted to a finite Laurent vector
z_tilde has D z_tilde=F A for some A in R. Its class is identified
with [A] in M. Here is a direct verification, without a metric claim:
changing the lift changes A by an element of (P,Q), and boundaries
give [A]=0. Conversely, if A=Pa+Qb, then
z_tilde+F(a,b) is a plane cycle, hence a plane boundary by regularity.
Its cylinder projection is the original class, which is therefore
zero. Finally, F M=0 ensures that every [A] is obtained this way.

## 3. A physical shear with common augmentation

Introduce new Laurent coordinates X,Y by the substitution

```
phi(x)=X,       phi(y)=XY,
f'=X phi(f),   g'=phi(g).
```

In physical site coordinates the two block maps are

```
f site (i,j) -> (i+j+1,j),
g site (i,j) -> (i+j,j).                                (3.1)
```

They are bijections on Z_ell x Z, preserve finite y support, and
preserve the Hamming weight exactly. The x period becomes X^ell+1;
no Y period is introduced. The transformed checks are

```
P'=X^-1 phi(P)=X^2+Y+XY^2,
Q'=phi(Q)=1+X^3+(X+X^2+X^3)Y.                           (3.2)
```

They obey

```
P'f'+Q'g'=phi(Pf+Qg).
```

Boundaries also transform exactly: (Qh,Ph) goes to
(Q'H,P'H), where H=X phi(h). This is a physical equivalence, unlike
a general change of algebraic generators.

Both checks have the same vertical augmentation:

```
P'(X,1)=Q'(X,1)=1+X+X^2=chi(X).                         (3.3)
```

For any actual cycle define

```
c(X)=f'(X,1)+g'(X,1)
    =X f(X,X)+g(X,X)  modulo X^ell+1.                  (3.4)
```

Then chi c=0. Each physical occupied site contributes at most one
output site, so |c|<=wt(f,g).

The binary cyclic kernel of chi at any length divisible by three
consists of zero and the three repetitions of 110,101,011. To prove
this, its coefficient equation is
c_a+c_(a-1)+c_(a-2)=0; the next equation gives c_(a+1)=c_(a-2),
so the sequence has period three. Each nonzero solution has exactly
two ones per triple. Therefore

```
c nonzero  =>  |c|=2ell/3.                              (3.5)
```

## 4. The kernel also has the sharp floor: even occupancy and a seam

Suppose an actual cycle has weight w<2ell/3. By (3.5), c=0. In the
sheared physical coordinates (3.1), this says that the total number
of occupied sites in each column, counting both blocks and all Y
rows, is even. Consequently, if K is the number of occupied columns,

```
2K<=w.                                                  (4.1)
```

Label an X seam by the column immediately to its right. From (3.2),
an f' qubit in column a has syndrome X offsets in {0,1,2} and can
cross only seams a+1,a+2. A g' qubit in column a has offsets in
{0,1,2,3} and can cross only seams a+1,a+2,a+3. Here “cross” may
include an intervening seam even if no syndrome site lies exactly
there; this only overcounts the forbidden cuts.

Crucially, all occupied sites anchored in the same column, regardless
of block or Y row, have their possible crossed seams in that same
set of at most three seams. Hence the number of forbidden seams is
at most

```
3K <= 3w/2 < ell.                                       (4.2)
```

Choose an un-crossed seam and represent the cylinder columns by the
ell consecutive integers beginning just to its right. Every occupied
qubit and all its syndrome neighbors then lie within this interval:
all offsets in (3.2) are nonnegative, and a neighbor beyond the right
endpoint would cross the chosen seam. Thus the cylinder cycle lifts
to a finite Laurent plane cycle with precisely the same occupied
qubits and precisely the same syndrome cancellations. There is no
Y identification to undo. No distinct sites merge during this lift.

The transformed plane complex is regular by Section 2 and the
physical equivalence in Section 3. The lift is therefore a plane
boundary. Project its finite boundary generator back to the cylinder
and invert the shear. The original cycle is a cylinder boundary.
This proves the all-class lower bound in (1.1).

The argument does not require that the support be connected, that a
class-minimal representative be connected, or that subtracting any
boundary reduce weight. It handles the given vector directly.

## 5. Exactly what the invariant detects

The invariant c alone detects the simple component, not the doubled
one. This distinction is essential and is not an obstacle to Section 4.
Let W_ell=(x^ell+1)/chi. For a lift with D z_tilde=F A, evaluating
y=x in the Laurent domain gives

```
chi (x f_tilde(x,x)+g_tilde(x,x))=F A(x,x),
x f_tilde(x,x)+g_tilde(x,x)=W_ell A(x,x).                (5.1)
```

After x-period reduction, the right side depends exactly on
A(x,x) modulo chi. Multiplication by W_ell identifies L/(chi)
with ann_(L/(F))(chi): if F divides W_ell a, cancellation in L
shows that chi divides a, and every element of that annihilator
is a multiple of W_ell by the same cancellation argument.

The residue A(x,x) modulo chi is precisely the projection of [A]
onto M_simple in (2.2). Thus the homological kernel of c is the
entire four-dimensional M_double, including its unit and nilpotent
directions. Nonzero simple projection has the lower bound directly
from (3.5). All classes with zero simple projection satisfy the even
column condition, so any nontrivial such class has the same lower
bound by the empty-seam argument. There is no assertion that c is
injective on the full six-dimensional class space.

## 6. Exact physical representatives

### 6.1 Sharp simple stripe

For ell=6r set

```
W_ell=(1+x) sum_(j=0)^(2r-1) x^(3j),
z_simple=(0,W_ell).
```

The exact unfolded-y identity is

```
D z_simple=V(x^ell+1).                                 (6.1)
```

Its period coefficient V is one on M_simple, so the cylinder
cycle is nontrivial. Its 4r physical sites are distinct and lie
in a single y row. This attains the lower bound, proving (1.1).

### 6.2 A genuine doubled-component family, not a shorter counterexample

An explicit all-r representative of a doubled-component unit is

```
R_ell=(x^ell+1)/chi^2
     =(1+x^2) sum_(j=0)^(r-1) x^(6j),
z_double=(chi R_ell, (y+x^3+x^2) R_ell).                (6.2)
```

Its exact finite Laurent syndrome, with y never folded, is

```
D z_double=xU(x^ell+1).                                (6.3)
```

Indeed y+x^3+x^2=U+x chi, and (2.1) gives
chi P+(U+x chi)Q=x chi^2(1+V)=xU chi^2. On M_simple
the period coefficient xU vanishes. On M_double, U=1+x chi is
a unit, so xU is nonzero and a unit. This certifies nontriviality
in precisely the component not detected by c.

The physical weights can be read without a search. Write
S=sum_(j=0)^(r-1) x^(6j). The f block is
(1+x+x^3+x^4)S, of weight 4r. In the g block the y row is
(1+x^2)S, of weight 2r, and row zero is
(x^2+x^3+x^4+x^5)S, of weight 4r. All displayed sites are
distinct modulo ell. Thus |z_double|=10r. Its invariant c is
identically zero even before period reduction. This is a verified
logical family, but not an improved upper bound or a refutation.

## 7. Consequence and bounded verification record

The already-proved `IBM_ONE_CHECK_CYLINDER_REDUCTION.md` gives
d_tor(ell,m)=d_cyl(ell) when both periods are divisible by six and
ell<=m. Combining that result with (1.1) now gives the exact
five-term-modification family

```
[[2ell m,12,2ell/3]],       6|ell, 6|m, ell<=m,
[[72r^2,12,4r]]            on ell=m=6r.                 (7.1)
```

Both CSS sectors are covered by the physical equivalence already
proved for that modified code. The cylinder proof here does not
assume the torus equality; only this corollary uses it. It does not
assert (7.1) for ell>m or for the original IBM checks. No novelty
claim is made.

The companion script `ibm_one_check_cylinder_tangent.py` checks only
prescribed finite identities: the physical shear and its inverse,
syndrome and boundary conjugacy, both augmentations, the three
period-three patterns, and the two prescribed ell=6 representatives
with their unfolded-y syndromes and primary-component coefficients.
It is not a distance verifier, and the all-r assertions above do
not depend on extrapolating these checks.

The fixed diagnostic passed in 0.000103 seconds with measured peak
RSS 17,367,040 bytes. No support, coset, height, or distance search
was performed. Only this new note and its assigned companion script
were created. The sharp cylinder target for this modification is
closed; the unchanged IBM metric problem remains separate.

## 8. Independent corollary audit: lattice, class injection, and parameters

The final `IBM_ONE_CHECK_CYLINDER_REDUCTION.md` was independently
read against the proof above. Its reduction and the resulting exact
parameters pass the following direct checks. This section records
the certificate without changing that parent note.

### 8.1 The shear is on the correct covering space

The cylinder's sole period vector is (ell,0). The linear part of
(3.1) sends (ell,0) to (ell,0), so it induces an honest bijection
of Z_ell x Z. The first-block translation by (1,0) is also such
a bijection. The inverse maps are

```
f' site (a,b) -> (a-b-1,b),
g' site (a,b) -> (a-b,b).
```

Both preserve finite Laurent support. Thus the ell seams counted
in Section 4 are exactly the seams of the remaining cyclic
coordinate. They are not seams of a different quotient lattice.
All check X offsets lie between zero and three, and ell>=6; the
crossing and lift argument applies also at r=1. The construction
never imposes a Y period or applies the shear to a torus while
silently discarding its transformed second period.

Once an unused seam gives a transformed finite plane cycle,
plane regularity makes it an actual plane boundary. Its finite
generator projects back through the exact chain equivalence in
Section 3. Consequently no nonzero cylinder class, tangent or
otherwise, can survive this plane-boundary conclusion.

### 8.2 A cylinder charge cannot disappear on torus projection

Set G=y^m+1, with 6|m. In M, both F and G vanish. For G, a
self-contained certificate supplementing Section 2 is

```
(y^2+y+1)^2=P^2+(x^6+1),
y^6+1=(y+1)^2(y^2+y+1)^2.
```

Let a cylinder cycle have a finite lift z_tilde with
D z_tilde=F A. Suppose its torus projection were a boundary.
Lifting that equation gives finite Laurent vectors a,b and a
finite Laurent scalar h with

```
z_tilde=delta(h)+F a+G b,
F(A+Da)=G Db.                                           (8.1)
```

The independent period polynomials form a regular sequence in
R; specifically, G is monic with unit leading coefficient in y
over (F2[x^+/-1]/F)[y^+/-1], so it is a non-zero-divisor there.
Reducing (8.1) modulo F therefore gives Db=F c for a finite
Laurent c. Cancellation of F in the domain R then gives
A+Da=G c. Reducing this equality in M yields [A]=0. Section 2
now makes the original cylinder class zero. This proves actual
class injection and rules out a lost nonzero period charge,
without an isometry claim about arbitrary projected vectors.

### 8.3 Metric equality and the complete square parameters

For a torus cycle of weight less than 2m/3, evaluation x=omega
forces f(omega,y)=0: otherwise its chi_y recurrence would occupy
at least 2m/3 y rows. Hence each occupied f row has at least two
qubits. If R_f,R_g count occupied rows in the two blocks, the
number of potentially crossed y seams is at most

```
2R_f+R_g <= |f|+|g| < 2m/3 < m.
```

Here P has y span two and Q has y span one, in the actual
unsheared torus coordinates. The unused seam gives precisely
the finite-y cylinder lift needed in Section 1, with unchanged
weight and syndrome. No statement about connected minimizers
enters this argument.

When 6|ell, 6|m and ell<=m, any torus cycle of weight less than
2ell/3 has such a lift. Section 4 makes the lift a boundary,
and projection makes the torus cycle a boundary. Conversely,
the stripe in Section 6 has weight 2ell/3 and nonzero cylinder
period coefficient [V]; Section 8.2 proves its torus projection
is still nontrivial. Thus both inequalities are attained:

```
d_tor(ell,m)=d_cyl(ell)=2ell/3.                          (8.2)
```

To include redundancies explicitly, write N=ell m. The cokernel
of the torus syndrome map is R/(P,Q,F,G)=M, of binary dimension
six. Hence rank D=N-6. Coordinate inversion and block swap give
rank delta=N-6 and identify the two CSS sectors by a physical
permutation. Therefore

```
k=2N-(N-6)-(N-6)=12.
```

The two blocks contain three and five check sites, respectively.
For ell,m>=6 these are eight distinct physical sites. Substituting
ell=m=6r into (8.2) proves, for every integer r>=1,

```
[[72r^2,12,4r]],       check weight eight,
P=x^3+y+y^2,          Q=(1+x+x^2)(y+x+1).               (8.3)
```

Its ordinary tangent ideal is exactly the one in (2.2), retained
from the original pair as proved in `IBM_ONE_CHECK_MODIFICATION.md`.
Neither the dimension nor the metric conclusion in (8.3) is a claim
for the unchanged IBM checks or for the separate Q=chi family.
This corollary audit used no new computation.
