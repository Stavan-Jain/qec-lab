# An exact k=12 family after changing one original IBM monomial

2026-09-08. Parent-supplied blueprint, independently checked and
written as a self-contained physical-weight proof. No distance,
support, coset, or height search is used. The extension at multiples
of three is included with its different dimension and distance.

## 1. Exact statement and relationship to the original problem

Keep the original Laurent check

```
A=1+y+x^3 y^-1.
```

In the other original check change **only its last monomial**:

```
B_original=1+x+x^-1 y^-3,
B_star=1+x+x^5.
```

Use a square torus of periods `ell=m=7r`, with integer r>=1. These
are new frames, not the original IBM b-bit parameter rectangles.
For a polynomial pair the CSS convention is
`H_X=[A | B_star]`, `H_Z=[B_star^T | A^T]`.

**Theorem 1.1.** This connected, check-weight-six BB CSS family has

```
[[98r^2, 12, 4r]]                   if 3 does not divide r,
[[98r^2, 20, 10r/3]]                if 3 divides r.    (1.1)
```

Both CSS sector distances equal the stated distance. For a constant
k=12 family indexed by every nonnegative integer t, one may set
`r=3t+1`. Then the exact parameters are

```
[[98(3t+1)^2, 12, 4(3t+1)]],       t>=0.             (1.2)
```

This is a theorem for the explicitly modified code, **not** the
original IBM distance conjecture. It does not preserve the original
ordinary primary ideal and is not the earlier Q=chi V modification.
No novelty claim, or claim of inequivalence to known constructions,
is made. The closeness assertion is the literal one-term change
above, not an algebraic identification of physical metrics.

We use the monomial-cleared pair

```
P=y A=x^3+y+y^2,
Q=B_star=1+x+x^5=chi K,
chi=1+x+x^2,       K=1+x^2+x^3.                      (1.3)
```

This clearing preserves physical weight. Explicitly, the block
permutation `(f,g)->(y^-1 f,g)` takes cycles for (A,B_star) to
cycles for (P,Q). It takes `(B_star h,A h)` to
`(Q y^-1 h,P y^-1 h)`, so also takes boundaries to boundaries.
No nonmonomial basis operation is used as a weight isometry.

## 2. The classical x-periodic kernel and its exact weight

Put `S_ell=F2[x]/(x^ell+1)` and

```
C_ell=ann_(S_ell)(Q).
```

The polynomial K is irreducible over F2: its degree is three and
it has no F2 root. A root alpha lies in F8 and has order seven,
since it is nonzero and is not one. In particular K divides
x^7+1. The factor chi is irreducible of degree two, is coprime
to K, and divides x^3+1. Q is squarefree.

If `3 does not divide r`, chi is invertible in S_ell because its
roots have order three, whereas ell=7r. Hence

```
C_ell=ann(K),       dim_F2 C_ell=3.                  (2.1)
```

Every such word is seven-periodic. More explicitly its coefficients
are exactly

```
w_j=Tr_(F8/F2)(c alpha^(-j)),       c in F8.
```

Substitution gives `w_j+w_(j-2)+w_(j-3)=0`, the convolution
recurrence for K. These are all its solutions: the recurrence is
determined by three consecutive coefficients, and the trace words
form a three-dimensional space. To see injectivity in c, for
c!=0 multiplication by c permutes F8 and the trace is nonzero
(`Tr(1)=1`), so the word is not identically zero.

For c!=0 the trace functional takes the value one at four of the
eight field elements, none of which is zero. The seven powers of
alpha run over all nonzero elements. Thus every nonzero word has
four ones per seven sites, and weight 4r on circumference 7r.

If `3 divides r`, the coprime factors split the kernel exactly:

```
C_ell=ann(chi) direct-sum ann(K),
dim_F2 C_ell=2+3=5.                                 (2.2)
```

This follows from Bezout projectors for chi,K on S_ell and does
not require the period polynomial to be squarefree. A nonzero
ann(chi) word has period three and weight two per triple; a
nonzero ann(K) word has period seven and weight four per seven.

For their sum with both components nonzero, the residue pairs
`(j mod 3,j mod 7)` occur exactly once in each block of 21 sites.
Its XOR weight in that block is therefore

```
2*(7-4)+(3-2)*4=10.
```

Pure chi and pure K words have weights fourteen and twelve per
block of 21, respectively. Consequently the exact minimum nonzero
weight in C_ell is

```
d_class(ell)=4r       if 3 does not divide r,
d_class(ell)=10r/3    if 3 divides r.                 (2.3)
```

All arguments include even r. Repeated factors in x^ell+1 do not
enlarge these kernels beyond the stated degrees, nor introduce
longer recurrence solutions.

### Explicit minimum words

The fixed identity

```
K(1+x^2+x^3+x^4)=1+x^7
```

gives the word

```
W_K=(1+x^2+x^3+x^4) sum_(j=0)^(r-1) x^(7j),
K W_K=x^ell+1,       |W_K|=4r.                      (2.4)
```

It attains (2.3) when 3 does not divide r. When 3 divides r, use

```
W_chi=(1+x) sum_(j=0)^(ell/3-1) x^(3j),
w_min=W_chi+W_K.                                    (2.5)
```

The preceding CRT count proves its exact weight 10r/3. Thus these
are prescribed analytic witnesses, not outputs of an optimization.

## 3. A cylinder normal form used only for class detection

Now let only x be periodic and let y have finite Laurent support.
Write `C=S_ell[y,y^-1]`, with cycles `Pf+Qg=0` and boundaries
`(Qh,Ph)`.

**Lemma 3.1.** Every cylinder class has a unique representative

```
(0,u+y v),       u,v in C_ell.                       (3.1)
```

Proof. P is monic of degree two in y with unit constant coefficient
x^3. Hence `C/(P)=S_ell direct-sum S_ell y`; y is already a unit
in that quotient since `y(y+1)=x^3`. Reduce g modulo P by adding
a boundary. Its new form is u+yv. The new first block f' satisfies

```
P f'+Q(u+yv)=0.
```

The second term has y degrees zero and one. If f' had a negative
lowest degree, its contribution there could not vanish because
x^3 is a unit. If f' had nonnegative highest degree, the monic
top term of P would produce degree at least two. Thus f'=0 and
Qu=Qv=0. Conversely these conditions give a cycle. If a vector
(3.1) were a boundary, `Ph=u+yv`; the same lowest/highest-degree
argument forces h=0. This proves uniqueness. QED.

In particular cylinder homology has dimension twice that of C_ell.
The normal-form operation itself is not claimed to preserve or
decrease Hamming weight. Only its algebraic uniqueness is used next.

## 4. Two physical contractions prove the exact cylinder floor

The crucial identities hold in the unbounded Laurent plane:

```
P(x,x^-1)=x^-2 Q(x),
P(x,x^4)=x^3 Q(x).                                  (4.1)
```

Define, on actual cylinder chains,

```
E_minus(f,g)=x^-2 f(x,x^-1)+g(x,x^-1),
E_four(f,g)=x^3 f(x,x^4)+g(x,x^4).                   (4.2)
```

Each input physical monomial maps to one x-periodic coordinate.
Parity summation and the displayed monomial shifts therefore give

```
|E_minus(f,g)|<=|f|+|g|,
|E_four(f,g)|<=|f|+|g|.                             (4.3)
```

This remains true for negative y exponents. For a cycle, (4.1)
gives `Q E_minus=Q E_four=0`. For a boundary, its two terms cancel
in each E. Thus these are genuine nonexpansive homology maps into
the physical classical code C_ell.

On the normal representative (3.1), their values are

```
E_minus=u+x^-1 v,       E_four=u+x^4 v.
```

Since `x^-1 Q=1+x^-1+x^4` and Qv=0,

```
E_minus+E_four=(x^-1+x^4)v=v.                       (4.4)
```

If both maps vanish, v=0 and then u=0. They therefore jointly
detect every nonzero **cylinder** class. At least one output of
every representative of such a class has weight at least d_class,
so (4.3) proves the cylinder floor d_class. Conversely `(0,w_min)`
has that weight and has nonzero E outputs. Consequently

```
d_cyl(ell)=d_class(ell).                             (4.5)
```

This proof uses a normal form to check class detection, not to
estimate the weight of an arbitrary representative. It supplies
the missing metric justification explicitly through (4.3).

## 5. A low-weight square cycle has an uncrossed row seam

Let both periods be ell=7r. Evaluate a torus cycle at x=alpha,
where K(alpha)=0. This evaluation is well-defined because
alpha^ell=1. Set `F(y)=f(alpha,y)`. Then

```
(alpha^3+y+y^2)F(y)=0 modulo y^ell+1.               (5.1)
```

In F8, put beta=alpha^-1 and gamma=alpha^4. From Q(alpha)=0,

```
beta+gamma=1,       beta gamma=alpha^3,
alpha^3+y+y^2=(y+beta)(y+gamma).                    (5.2)
```

The two roots are distinct and both have order seven. With the
coefficient-convolution convention in (5.1), its coefficient
sequences are exactly

```
F_j=a beta^(-j)+b gamma^(-j),       a,b in F8.        (5.3)
```

The inverse powers are important for this convention. The two
initial coefficients determine a,b uniquely, and substitution
verifies the recurrence. Thus this describes all solutions,
including on an even multiple of the seven-row period.

If exactly one of a,b is nonzero, all seven entries are nonzero.
If both are nonzero, exactly one of the seven entries vanishes:
the ratio gamma/beta=alpha^5 runs over F8's nonzero elements.
Every nonzero solution of (5.1) therefore has F8 row weight at
least 6r.

Evaluation cannot create an occupied row, so |F|<=|f|. Suppose
the physical total weight of the torus cycle is less than
`d_class(ell)<=4r`. It follows that F=0. Each nonzero binary
physical f row then vanishes at alpha and contains at least two
qubits; a single binary monomial cannot vanish at alpha.

Let R_f count these occupied f rows. An f qubit has syndrome rows
j,j+1,j+2, so all qubits on a given f row cross at most two row
seams. A g qubit crosses none, since Q is independent of y. Hence
the union of crossed seams has size at most

```
2 R_f <= |f| <= |f|+|g| < d_class(ell) <=4r<ell.     (5.4)
```

There is an uncrossed seam. Cut there to obtain a finite-y cylinder
chain with the same physical weight and zero syndrome. More
explicitly, after rotating the seam to the end of rows 0,...,ell-1,
the last two f rows vanish. All terms of Pf and Qg then stay within
that interval, so torus cancellation is exactly cylinder cancellation.
No condition on the g rows or connectedness of the cycle is needed.

If the torus cycle were nontrivial, this lift could not be a cylinder
boundary, since projecting a boundary gives a boundary. Its weight
would contradict (4.5). This proves the all-class torus lower bound
d_class, including classes not detected by the two E maps on the
torus. Joint injectivity of those maps was asserted only on the
cylinder; it is not an extra torus assumption.

## 6. Nontrivial attaining torus words and both sectors

On the square torus both substitutions in (4.2) are well-defined:
`y^ell` maps to x^-ell or x^(4ell), both equal to one. The E maps
still vanish on torus boundaries. The explicit vector

```
z=(0,w_min)
```

is a cycle since Qw_min=0. Its E_minus value is the nonzero binary
word w_min itself, so z is not a torus boundary. Its physical weight
is exactly d_class by §2. This matches the lower bound from §5.

Let R reverse both torus coordinates. Translation matrices satisfy
`R P R=P^T` and `R Q R=Q^T`. The physical permutation

```
Pi(f,g)=(R g,R f)
```

takes the two CSS cycle/boundary quotients to one another. Thus the
two sector distances are equal, and both equal (2.3).

## 7. Dimension, check weights, connectedness, and even r

Let N=ell^2 and D=[P | Q] on the finite torus ring T. Its cokernel
is `T/(P,Q)`. Before imposing periods, the factorization Q=chi K
and (5.2) give the reduced plane algebra

```
F2[x^+/-1,y^+/-1]/(P,Q) = F4^2 direct-product F8^2.  (7.1)
```

On the chi component x has order three, x^3=1, and
P=y^2+y+1 splits into its two distinct F4 roots. This component
survives both square periods precisely when 3 divides r.
On the K component x=alpha and y is beta or gamma; all three
have order seven, so both F8 components always survive. Therefore

```
dim_F2 coker D = 6+4*[3 divides r].                 (7.2)
```

This quotient is reduced by Q and P themselves. Repeated factors
in the period polynomials when r is even do not add multiplicity.

The boundary matrix is delta=[Q;P]. With Pi and R from §6,

```
delta=Pi D^T R,
```

so its rank equals rank D. Hence the CSS logical dimension is

```
k=2N-rank D-rank delta
 =2 dim coker D=12+8*[3 divides r].                 (7.3)
```

Commutation follows from PQ=QP in the translation ring. Both checks
have three distinct terms for ell>=7, so every CSS check has exactly
six sites. Each qubit is incident to three checks of each CSS type.
For connectedness, the P support difference between y and y^2
generates a one-row move among syndrome vertices, and the Q support
difference between 1 and x generates a one-column move. Thus the
H_X incidence graph, and consequently the full code incidence
graph, is connected.

Equations (2.3), (7.3), and n=2ell^2 complete Theorem 1.1.

## 8. Verification and exact scope boundary

The companion script checks finite Laurent identities, the F8 root
relations, the two contraction identities, the normal-form detector
identity, and the prescribed stripes at ell=7 and ell=21. It checks
their physical weights, zero torus syndromes, and nonzero E_minus
values. It does not enumerate codewords, coefficients in classes,
supports, heights, or distances. The universal minima, seam argument,
dimension count, and parameter exception are the analytic proofs
above, not extrapolations from these two witnesses.

The final fixed diagnostic passed in 0.000080 seconds with measured
peak RSS 17,039,360 bytes. It installs a 55-second wall alarm and
50-second CPU limit and checks the 768 MiB RSS cap. Positive field
powers are multiplied directly, so the order-seven identity is
actually checked rather than imposed by exponent reduction.

The original IBM family remains an unresolved distance problem.
The distinct chi V cylinder floor is now proved independently in
`IBM_ONE_CHECK_CYLINDER_TANGENT.md`; the complete corresponding
modified-code theorem is `IBM_MODIFIED_EXACT_DISTANCE_THEOREM.md`.
That proof passed independent read-only review. No metric transfer
between those codes and this family is asserted or needed. In
particular this family makes no claim to preserve the original
ordinary primary ideal: (7.1) is a different, reduced plane algebra.
What is complete here is the exact theorem for this one-monomial
modification on the stated new square frames, including the
explicitly tracked 3-divisibility case.
