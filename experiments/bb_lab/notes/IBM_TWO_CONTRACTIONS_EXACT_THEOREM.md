# Two physical contractions: an exact primitive-trinomial BB construction

2026-09-08. Bounded Codex theory audit and proof write-up of the
parent's construction. This note is the only new file. No search,
distance computation, new code, or novelty claim is involved.

## 1. Hypotheses and exact theorem

Let

```
q(x)=1+x^s+x^t,       0<s<t,       t>=2,
```

be a **primitive binary polynomial of degree t**: it is irreducible
over F2 and a root alpha generates the multiplicative group of
`K=F_(2^t)`. Put

```
N0=2^t-1,       L=N0*r,       r>=1,
P=x^(s+t)+y+y^2,       Q=q(x),
d0=2^(t-1)*r.
```

No existence assertion for primitive trinomials in every degree is
being made; primitiveness of the specified polynomial is a theorem
hypothesis.

On the physical square `Z_L times Z_L`, use the binary CSS matrices

```
H_X=[M_P | M_Q],       H_Z=[M_Q^T | M_P^T].
```

Here multiplication by a Laurent monomial is the corresponding
translation permutation. The cycle and boundary maps in one sector
are

```
D(f,g)=Pf+Qg,       delta(h)=(Qh,Ph).
```

**Theorem 1.1.** For every positive integer r, including even r,
this is a connected, weight-six-check BB code with exact parameters

```
[[2(2^t-1)^2 r^2, 4t, 2^(t-1) r]].                (1.1)
```

Both CSS sector distances are d0. The proof uses physical Hamming
contractions, an exact empty-seam lift, and a physical nonboundary
attainer. An algebraic normal form is used only to prove joint
injectivity of the contractions, never to transfer Hamming weight.

The family is explicitly modified. This theorem does not establish
the distance conjecture for the original IBM Laurent pair, preserve
its nonreduced ordinary module, or assert optimality or novelty.

## 2. The classical kernel at every repeated period

Write

```
R_L=F2[x,x^-1]/(x^L+1),       C_L=ker(q:R_L->R_L).
```

**Lemma 2.1.** The binary code C_L has parameters

```
[L,t,2^(t-1) r],
```

and every nonzero codeword has exactly the displayed weight. Its
words repeat a length-N0 binary simplex word r times.

**Proof.** Since alpha has order N0, q divides `x^N0+1`, hence
also `F=x^L+1`. Put `W_L=F/q`. Cancellation in the one-variable
polynomial domain shows

```
ker(q mod F) = {W_L a mod F : a mod q}.
```

The parametrization by a modulo q is injective, so the kernel has
dimension t. This argument does not require F to be square-free.

There is also an explicit description of all its words. For c in K,
take the cyclic binary sequence

```
b_j=Tr_(K/F2)(c alpha^(-j)),       0<=j<L.           (2.1)
```

Its q-convolution vanishes because `q(alpha)=0`: the coefficient
at row j is the trace of `c alpha^(-j)(1+alpha^s+alpha^t)`.
The inverse powers alpha^(-j) run through all nonzero field elements
once per N0 positions. For nonzero c, the map `z -> Tr(cz)` is a
nonzero binary linear functional and takes the value one on exactly
`2^(t-1)` elements of K, none of which is zero. Thus (2.1) has
weight `2^(t-1)` in each period and weight d0 in length L.

For completeness, the trace is not identically zero: its polynomial
`z+z^2+...+z^(2^(t-1))` is nonzero of degree below `|K|`. For c
nonzero, multiplication by c is a field bijection. These observations
also prove that the t-dimensional parametrization (2.1) is injective.
It therefore exhausts C_L. QED.

In particular W_L itself is a nonzero word of weight d0. The finite
identity

```
W_L = W_N0 sum_(j=0)^(r-1) x^(j N0)
```

has disjoint blocks since `deg W_N0=N0-t<N0`. It is a literal
repetition identity, not an assumption about large-r behavior.
Repeated roots of F at even r produce no extra vectors in ker q.

## 3. Cylinder normal form and two physical contractions

For this section only x is periodic. Let

```
C=R_L[y,y^-1],
```

and require finite support in the unbounded y coordinate. Define

```
E_s(f,g)=x^s f(x,x^s)+g(x,x^s),
E_t(f,g)=x^t f(x,x^t)+g(x,x^t).                     (3.1)
```

The substitutions are well-defined on Laurent polynomials because
x is a unit. Every input physical site maps to one output physical
site, with possible parity cancellations. Separately for a=s,t,

```
|E_a(f,g)| <= |f|+|g|.                             (3.2)
```

The exact finite identities

```
P(x,x^s)=x^s q(x),
P(x,x^t)=x^t q(x)                                  (3.3)
```

show that E_s and E_t map cycles into C_L and vanish on boundaries.

### 3.1 The normal form is algebraic, not metric

Modulo P, the Laurent ring C has free R_L basis 1,y. Indeed, P
is monic of degree two in y, and its constant coefficient
`x^(s+t)` is a unit. In the quotient,

```
y(y+1)=x^(s+t),
```

so y is invertible. Hence every g admits a unique remainder
`u+yv`, with u,v in R_L, and a finite Laurent decomposition

```
g=Ph+u+yv.                                        (3.4)
```

For a cycle, reduction of `Pf+qg=0` modulo P forces
`qu=qv=0`, because 1,y are a free basis. The equation then becomes
`P(f+qh)=0`. Multiplication by P on C is injective: the highest
y coefficient of a nonzero finite Laurent polynomial survives
multiplication by the leading coefficient one. Therefore f=qh,
and subtraction of the actual boundary delta(h) leaves

```
(0,u+yv),       u,v in C_L.                        (3.5)
```

The remainder makes this representative unique as a homology
normal form. No assertion that the subtraction decreases weight
is used or needed.

### 3.2 Joint injectivity on every cylinder class

On (3.5), the invariants are

```
E_s=u+x^s v,       E_t=u+x^t v.
```

Since qv=0, `(x^s+x^t)v=v`. Consequently

```
v=E_s+E_t,       u=E_s+x^s(E_s+E_t).               (3.6)
```

In particular E_s=E_t=0 if and only if the cylinder class is zero.
This proves joint injectivity on the entire cylinder homology,
including every class at even r.

Every nontrivial cylinder cycle has at least one nonzero invariant.
Using its individual contraction bound (3.2) and Lemma 2.1 gives

```
d_cyl(L)>=d0.                                     (3.7)
```

The stripe `(0,W_L)` is a cylinder cycle and has E_s=W_L nonzero,
so in fact `d_cyl(L)=d0`. We use neither a bound on the sum of
the two output weights nor a metric claim for (3.5).

## 4. The evaluated row recurrence and an empty torus seam

We now put both coordinates on period L. Let

```
beta=alpha^s,       gamma=alpha^t.
```

The primitive-polynomial relation gives `beta+gamma=1`; hence
beta and gamma are distinct nonzero elements of K. At x=alpha,

```
P(alpha,y)=(y+beta)(y+gamma),       Q(alpha)=0.      (4.1)
```

Let h be the multiplicative order of beta/gamma. It is a divisor
of the odd number N0 and is greater than one, so h>=3. Both beta
and gamma, and their ratio, have orders dividing L.

**Lemma 4.1.** Every nonzero cyclic K-valued sequence F annihilated
by `(y+beta)(y+gamma)` has weight at least

```
L(1-1/h) >= 2L/3.                                 (4.2)
```

**Proof.** The coefficient recurrence is second order, with nonzero
constant coefficient and distinct characteristic modes. Its full
solution space consists of

```
F_j=A beta^(-j)+B gamma^(-j),       A,B in K.        (4.3)
```

The two modes are independent because beta and gamma are distinct.
They are both L-periodic, so this describes every cyclic solution,
including when L is even. There is no generalized-root mode.

If precisely one of A,B is nonzero, F has no zero entry. If both
are nonzero, a zero entry requires
`(beta/gamma)^j=A/B`. This has either no solutions or one residue
class modulo h, containing exactly L/h indices modulo L. This
proves (4.2). QED.

Also, for t>=2,

```
d0/L=2^(t-1)/(2^t-1) <= 2/3 < 1.                  (4.4)
```

Consider a torus cycle z=(f,g) with `|z|<d0`. Evaluation x=alpha
is well-defined because alpha^L=1. It sends f to a sequence in
the kernel from Lemma 4.1. Evaluation cannot make an empty
physical row nonzero, so

```
|f(alpha,y)| <= |f| < d0 <= 2L/3.
```

The evaluated sequence must therefore vanish. Every binary row
of f is in the ideal generated by q, the kernel of evaluation
at alpha. In particular a nonempty f row has at least two sites:
a single binary monomial evaluates to a nonzero field element.

Let R_f be the number of occupied f rows. Then

```
2 R_f <= |f| < d0 < L.                            (4.5)
```

An f qubit in row j has syndrome rows j,j+1,j+2 and can cross only
the two row seams j+1,j+2. Every f qubit in that row has the same
two possible crossed seams, independent of its x coordinate.
Every g qubit has all its syndrome sites in its own row, because
Q is independent of y. Thus at most 2R_f<L seams are crossed.

Cut at an uncrossed seam. All occupied qubits and all their syndrome
neighbors have consistent lifts into one L-row interval. None of
their incidences wraps through the cut. Distinct physical sites
remain distinct, and every cancellation in the syndrome is retained.
This gives an injective, weight-preserving lift of the complete
cycle to the x-periodic, y-Laurent cylinder.

By (3.7), this lifted cycle of weight below d0 is a cylinder
boundary. Projecting its finite Laurent boundary generator back
to the torus proves that z was a torus boundary. Consequently

```
d_tor(L,L)>=d0.                                   (4.6)
```

This lift makes no connectedness or minimal-representative
assumption about the particular cycle. Its use of row divisibility
is only for counting crossed physical seams, not for transferring
weight through an algebraic quotient.

## 5. A torus attainer with an actual nonboundary detector

Take the same physical one-row stripe

```
z_L=(0,W_L),       q W_L=x^L+1.
```

It is a torus cycle and has weight d0 by Lemma 2.1. On the square
torus, E_s and E_t are still well-defined: substituting y=x^a
sends `y^L+1` to `x^(aL)+1=0`. Their boundary-vanishing identities
from (3.3) still hold. Since

```
E_s(z_L)=W_L!=0 in R_L,
```

the stripe is not a torus boundary. This argument remains valid
when r is even and W_L vanishes under some field evaluations:
the detector retains W_L in its full physical cyclic ring, not
merely its value at alpha. Thus

```
d_tor(L,L)<=d0.                                   (5.1)
```

Together with (4.6), this proves the exact distance in one sector.
Joint injectivity of E_s,E_t was asserted only on the cylinder.
It is not needed, and is not claimed, for the full torus homology;
the remaining torus winding direction is handled by the seam proof.

## 6. Dimension, both sectors, check weight, and connectedness

### 6.1 Dimension, including all repeated-root periods

Let R be the plane binary Laurent ring. Modulo q, x becomes alpha
in K, and (4.1) gives the exact quotient

```
M=R/(P,q)
 = K[y,y^-1]/((y+beta)(y+gamma))
 = K direct-product K.                            (6.1)
```

Both roots are nonzero, so Laurent inversion in y introduces no
extra relation. The binary dimension is 2t. Since alpha, beta,
and gamma have orders dividing N0, both periods `x^L+1,y^L+1`
already vanish on M. Therefore the finite torus cokernel of D
is exactly M and has dimension 2t.

Writing N=L^2, the rank of D is N-2t. Coordinate inversion
transposes the two circulant matrices; block swap then shows
`rank delta=rank D`. Since D delta=0,

```
k=dim H_1=2N-(N-2t)-(N-2t)=4t.                    (6.2)
```

No semisimplicity assumption on the period rings is used. In
particular even r and arbitrarily repeated period roots introduce
no exception to (6.2).

### 6.2 The other CSS sector

Let bar invert both lattice coordinates. The physical permutation

```
(f,g) -> (bar(g),bar(f))
```

maps the cycle equation `Pf+Qg=0` to
`Q^T bar(g)+P^T bar(f)=0`. It maps boundaries to the other CSS
sector's boundaries as well. It preserves physical weight, so
both sector distances equal d0.

### 6.3 Primitiveness forces the required connectivity condition

**Lemma 6.1.** A primitive polynomial of the specified trinomial
form has `gcd(s,t)=1`.

**Proof.** Suppose d=gcd(s,t)>1 and write t=d k. Then alpha^d is
a root of the degree-k polynomial
`1+z^(s/d)+z^k`. Its minimal polynomial has degree at most k,
so its multiplicative order is at most `2^k-1`. On the other hand,
primitiveness of alpha implies

```
ord(alpha^d)=(2^(d k)-1)/gcd(2^(d k)-1,d)
             >= (2^(d k)-1)/d
             > 2^k-1.
```

The strict inequality follows from
`(2^(d k)-1)/(2^k-1)=1+2^k+...+2^((d-1)k)>d`.
This contradiction proves the lemma. QED.

Thus no extra connectivity hypothesis is required. The Q support
differences generate horizontal steps s and t, hence unit horizontal
displacement by Lemma 6.1. The difference between the y and y^2
terms of P gives unit vertical displacement. Paths through qubits
therefore connect all syndrome sites on the square torus. Every
qubit is incident to a syndrome site, so the full incidence graph
is connected.

Finally, `0<s<t<N0<=L`, so Q has three distinct physical sites.
P has three distinct row offsets 0,1,2 because L>=3, even if its
x exponent wraps at the smallest parameter. The blocks are disjoint,
so every CSS check has exactly six sites. The construction does not
obtain growth by taking disjoint copies of a fixed-size code.

## 7. Explicit examples and honest relation to IBM

### 7.1 The flat k=8 example is included

For `q=1+x+x^2`, take s=1,t=2. The polynomial is irreducible
and primitive of order three. The theorem gives

```
P=x^3+y+y^2, Q=1+x+x^2,
L=3r: [[18r^2,8,2r]].                             (7.1)
```

Restricting to L=6r gives `[[72r^2,8,4r]]`, exactly the previously
proved flat family. Clearing its Laurent checks uses independent
monomial translations of the two blocks: `P=y A_flat` and
`Q=x B_flat`, where

```
A_flat=1+y+x^3 y^-1,       B_flat=1+x+x^-1.
```

This retains five of the original six Laurent monomials.

### 7.2 A fully explicit k=12 example with two changed monomials

Take `q=1+x+x^3`, s=1,t=3. It has no binary root and is therefore
irreducible of degree three. A root is a nonidentity element of
the multiplicative group of F8, whose order is the prime seven;
thus this trinomial is primitive. The theorem gives, for every r>=1,

```
P=x^4+y+y^2,       Q=1+x+x^3,
L=7r: [[98r^2,12,4r]].                            (7.2)
```

In the original Laurent presentation this is

```
A_new=1+y+x^4 y^-1,       B_new=1+x+x^3.
```

Relative to the specified IBM pair
`A=1+y+x^3 y^-1`, `B=1+x+x^-1 y^-3`, exactly two terms change:

```
x^3 y^-1 -> x^4 y^-1 in A,
x^-1 y^-3 -> x^3 in B.
```

Four of six Laurent monomials are unchanged. Clearing A by y
gives precisely P; B is already Q. For completeness, replacing
checks A,B by aA,bB for monomials a,b is a physical equivalence:
map `(f,g)` to `(a^-1 f,b^-1 g)`, and map a boundary generator h
to `(ab)^-1 h`. This sends old cycles and boundaries to the new
ones by coordinate permutations. No nonmonomial shear is involved
in counting the changed terms or transferring physical distance.

The periods and distance genuinely grow, the check weight remains
six, and k=12 is retained numerically. Nevertheless this is not the
original IBM family: the plane module (6.1) is reduced, the periods
are 7r rather than the original square 6r, and the distance coefficient
is different. Its exact square ratio is `n/d^2=49/8`. Keeping k=12
does not preserve the original ordinary tangent structure or solve
the original conjecture, including its exceptional parameter cases.

## 8. Audit disposition and boundaries of the result

The parent's proposed all-r proof passes after making the following
points explicit:

- The primitive-kernel distance is exact for repeated periods too;
  kernel dimension and trace sequences prove this without a scan.
- The two contractions are individually physical contractions and
  jointly injective only on cylinder homology. The normal form is
  not used metrically, and torus winding is handled separately.
- The evaluated row recurrence uses distinct roots and counts at
  most L/h zeros. This is valid at even r, with no extra root modes.
- The stripe detector lives in the full physical x-period ring,
  so nontriviality is not lost at repeated roots.
- Primitiveness itself forces gcd(s,t)=1, settling connectedness
  without an additional assumption.

There is no unresolved metric lemma within Theorem 1.1's hypotheses.
The theorem concerns the specified modified square family only.
No claim about arbitrary nonprimitive trinomials, arbitrary rectangular
aspect ratios, unchanged IBM distance, novelty, or optimality is made.
The separate one-changed-term k=12 construction being investigated
elsewhere is not used or rederived here.

This note was written from exact analytic identities and elementary
finite-field reasoning. No new computation, checker, support search,
coset search, height sweep, large model, old-process action, git change,
shared-file edit, or Lean work was performed.
