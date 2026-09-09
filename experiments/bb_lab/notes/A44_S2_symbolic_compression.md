# A44-S2: exact obstructions to inherited compression and two collective identities

Date: 2026-09-08. New-only continuation of A44-S1 on
`claude/a40-tour-de-gross-column`, PR #27. No earlier note or script was
modified. This note concerns the **original IBM polynomials**.

## 0. Result and stopping point

No weight-14 representative in Q_4, or weight-30 original r=3 attainer,
was constructed. The substantive results are:

* **Exact support tripling for every vector**, not merely the banked
  Frobenius motif. Consequently arbitrary stabilizer dressing at a
  smaller scale cannot evade the tripling cost of an unchanged lift.
* An elementary, self-contained proof that **Q_2 has X and Z distance
  exactly six**. It excludes all nonzero cycles below six using a
  constant-size column-overlap argument, without enumerating logicals.
* Therefore the entire inherited-boundary family
  `w_4+alpha beta q(z^2,tau^2)(alpha,beta)`, for **arbitrary q in Q_2**,
  has minimum weight **18**, not 14. This includes collective inherited
  patches of unrestricted size, not just the three-check obstruction
  from S1.
* Two exact collective factorizations are supplied with their boundary
  multipliers. The partial factorization has weight `8*3^(a-1)` for
  every `r=2^a>=2`, uniformly worse than `2*3^a`. The full scalar
  factorization costs 32 at r=4. Neither is a compression route.

These statements exclude explicit ansatz classes. They **do not exclude
an arbitrary new boundary correction in Q_4**, do not show that the
banked homology class is wrong, and do not refute IBM's distance law.
After these obstructions the tested factorizations are stopped; no SAT,
coset optimization, or further witness search is proposed as an output.

## 1. Ring, weights, motif, and explicit duals

Let `r=2^a`, `a>=0`, and

```
Q_r = F2[z^{+/-1},tau^{+/-1}]/(z^(6r)+1,tau^(2r)+1),
alpha = 1+z tau^-1+z^-4 tau,
beta  = 1+z^-1+z^-2 tau^3.
```

Q_r has `12r^2` group elements and two qubits per group element. A
physical pair `(f,g)` has weight `|f|+|g|`, is an X-cycle when
`beta f+alpha g=0`, and is a boundary when it is `h(alpha,beta)`.
All exponents in displayed finite supports are reduced modulo
`(6r,2r)` unless stated otherwise. The banked motif is

`w_r=(z^r beta^(r-1), z^-r tau^r alpha^(r-1))`.     (1.1)

This quotient is obtained from the original square `(6r,6r)` using
`z=x`, `tau=xy` and the additional period `tau^(2r)=1`. Multiplication
by `1+tau^(2r)+tau^(4r)` pulls a quotient vector back to the original
square, replicating every occupied qubit three times. Pushforward
after pullback is multiplication by three, hence the identity in
characteristic two; pullback is injective on homology. In particular
a quotient weight-14 attainer at r=4 would give an original weight-42
attainer. This is a periodic witness subspace, not a claim that the
original code decomposes into three components.

The two first relevant motifs, with their complete supports, are

```
w_2, on (z,tau) periods (12,4):
  f = {(0,3),(1,0),(2,0)}
  g = {(6,3),(10,2),(11,1)}

w_4, on (z,tau) periods (24,8):
  f = {(0,3),(0,6),(1,0),(2,0),(2,3),
       (3,0),(4,0),(22,1),(23,6)}
  g = {(8,7),(12,6),(13,5),(16,5),(18,3),
       (20,4),(21,3),(22,2),(23,1)}.
```

For explicit dual pairings, let `omega^2+omega+1=0` in F4 and use
`epsilon` as a formal nilpotent. Define `Z(epsilon)` by

```
Z=(1+epsilon)+(1+epsilon)^2 Z^-4,  Z(0)=omega.
```

Iteration from omega determines the series: differences gain a factor
four in epsilon-adic order. It solves
`alpha(Z,1+epsilon)=0`. Modulo the needed powers, its full constants are

```
Z = omega+epsilon+omega^2 epsilon^2                 mod epsilon^4,
Z = omega+epsilon+omega^2 epsilon^2
    +omega epsilon^4+omega epsilon^6                mod epsilon^8.
```

The substitutions `z=Z`, `tau=1+epsilon` respect the periods of Q_r
in `F4[epsilon]/epsilon^(2r)`. Their constant terms have order three
and one, respectively, and `2r` is a power of two. As in S1, define

```
lambda_r(P) = Tr_(F4/F2)(omega^(1-r)
                 [epsilon^(2r-2)] P(Z,1+epsilon)).             (1.2)
```

The physical Z-dual is specified completely by

`U_(r,0)[i,j]=lambda_r(z^i tau^j)`, `U_(r,1)[i,j]=0`.

Since lambda annihilates alpha times any polynomial, U annihilates
every X-check. The expansion

`beta(Z,1+epsilon)=epsilon^2+omega epsilon^3+O(epsilon^4)`

shows that the first component of w_r has leading term
`omega^r epsilon^(2r-2)`. Therefore `<U_r,w_r>=Tr(omega)=1`.
All boundary rewrites in this note keep this pairing equal to one.

In particular the r=2 first component evaluates to
`omega^2 epsilon^2+epsilon^3 mod epsilon^4`, and the r=4 first component
to `omega epsilon^6+omega^2 epsilon^7 mod epsilon^8`. Their respective
functionals are `Tr(omega^2 [epsilon^2])` and `Tr([epsilon^6])`.
For r=1, `w_1=(z,z^-1 tau)` has weight two and lambda is evaluation
at `(omega,1)` followed by trace. These formulas supply the
nontriviality used below without a spectral dimension assumption.

## 2. The lift is an exact Hamming similarity on the whole vector space

For an arbitrary vector in Q_s, define its lift to Q_(2s) by

`D_s(f,g)=(beta f(z^2,tau^2), alpha g(z^2,tau^2))`.  (2.1)

**Theorem 2.1.** For every vector, whether or not a cycle,

`|D_s(f,g)|=3(|f|+|g|)`.                           (2.2)

*Proof.* Squaring exponents embeds each source support injectively in
the even-even coordinates of the doubled frame. The three beta
monomials have different `(z,tau)` parity classes: `00,10,01`.
Their translates of this embedded support are disjoint. Thus the
first lifted block has exactly `3|f|` cells. The three alpha monomials
also have different parity classes, `00,11,01`, giving `3|g|`.
The two physical blocks are disjoint from each other. ∎

The algebraic identities are

```
syndrome(D_s(f,g)) = (beta f+alpha g)(z,tau)^2,
D_s(h alpha,h beta) = alpha beta h(z^2,tau^2)(alpha,beta).
```

Here squaring a source period relation makes a target period relation.
Thus D maps cycles to cycles, boundaries to boundaries, and
`D_s(w_s)=w_(2s)`. Importantly, a boundary inherited through this map
cannot exploit overlaps between the three copies: the copies occupy
distinct parity classes for every possible source representative.

Iterating from Q_1 gives an unconditional obstruction for every a.

**Corollary 2.2.** Among all vectors

```
w_r + (alpha beta)^(r-1) q(z^r,tau^r)(alpha,beta),  q in Q_1,
```

the minimum weight is exactly `2*3^a`. In particular this entire
pre-dressed-seed construction misses `4r-2` for all `a>=2`.

*Proof.* These are exactly the iterated lifts of
`w_1+q(alpha,beta)`. The latter are nontrivial by (1.2). Every cycle
has even weight, since each syndrome column has odd weight three;
therefore these nonzero cycles have weight at least two. The seed
attains two. Repeated application of (2.2) proves the claim. ∎

The next section removes any reliance on a banked distance floor and
allows arbitrary dressing at Q_2, not just at Q_1.

## 3. Exact distance six in Q_2, proved by local column identities

Work on `G=Z_12 x Z_4`. Syndrome columns are translates of

```
P_alpha={0,u,v}, u=(1,-1), v=(-4,1),
P_beta ={0,p,q}, p=(-1,0), q=(-2,3).
```

**Lemma 3.1.** Any two distinct syndrome columns intersect in at most
one syndrome coordinate, including columns from different blocks.

*Proof.* The six nonzero ordered differences within each of P_alpha
and P_beta are distinct in G. The nine differences from P_alpha to
P_beta are also distinct; explicitly the latter are

`{(0,0),(1,0),(2,1),(1,3),(2,3),(3,0),(8,1),(9,1),(10,2)}`.

A two-point intersection between two translates would equate two
different such differences. ∎

**Lemma 3.2.** Three distinct translates of either one of these
trinomials, whose pairwise intersections are all present and distinct,
have sum equal to a translate of the squared trinomial.

*Proof.* For a support `{0,U,V}`, intersection edges between its
translates have differences `+/-U,+/-V,+/-(V-U)`. There are no extra
three-edge relations from the torus. Indeed, if
`aU+bV=0` in G with `|a|,|b|<=3`, then for alpha the conditions are

`a-4b=0 mod 12`, `-a+b=0 mod 4`.

The first gives `a=0 mod 4`; the second gives `b=0 mod 4`.
The bounds imply `a=b=0`. For beta the conditions are

`-a-2b=0 mod 12`, `3b=0 mod 4`,

which likewise force `b=0`, then `a=0`. Thus a triangle of translates
lifts to a triangle in the ordinary lattice with edge directions
`+/-U,+/-V,+/-(V-U)`.

There are two triangle orientations, up to translation. Positions
`{h,h-U,h-V}` share the common syndrome point h, and are excluded by
the distinct-intersection hypothesis. Positions `{h,h+U,h+V}` have
three distinct pair intersections. Their sum is

`z^h (1+z^U+z^V)^2`,

with vector-exponent monomial notation, as claimed. ∎

**Theorem 3.3.** There is no nonzero X-cycle of weight less than six
in Q_2. Its X distance and Z distance both equal six.

*Proof of the lower bound.* Odd weights are impossible by augmentation
of the syndrome equation. Weight two is impossible by Lemma 3.1.
Suppose four distinct columns summed to zero. Each syndrome point
then occurs either twice or four times. Write n_2 and n_4 for their
numbers. Counting incidences gives `2n_2+4n_4=12`, while pairwise
intersections total `n_2+6n_4=6+4n_4`. Lemma 3.1 bounds that total
by six, so `n_4=0`, `n_2=6`. Hence every pair of columns intersects
once, and all six intersections are distinct.

If at least three columns belong to one block, Lemma 3.2 says their
sum is a translate of that block's squared trinomial. All its
z-coordinates have the same parity. Every single column in either
block has both z parities, so it cannot equal the fourth column.
This rules out the 4+0 and 3+1 block distributions and their reversals.

It remains to exclude two columns from each block. In each block
the two columns intersect, so their displacement is, up to sign,
one of the three support differences. Opposite signs only translate
the pair-sum. The complete horizontal-coordinate multisets of the
three alpha pair-sums and three beta pair-sums are:

| pair | displacement | z multiset modulo 12 |
|---|---|---|
| alpha | u | `{0,2,8,9}` |
| alpha | v | `{0,1,4,9}` |
| alpha | v-u | `{0,1,3,7}` |
| beta | p | `{0,9,10,10}` |
| beta | q | `{0,8,9,11}` |
| beta | q-p | `{0,9,11,11}` |

The beta p and q-p rows have a repeated z coordinate; no alpha row
does. The remaining beta row has sorted cyclic-gap multiset
`{1,1,2,8}`. The three alpha rows have respectively
`{1,2,3,6}`, `{1,3,3,5}`, `{1,2,4,5}`. Translation preserves these
invariants, so no alpha pair-sum equals a beta pair-sum. This rules
out 2+2 and completes the lower bound.

The explicitly paired w_2 in §1 attains weight six, proving X distance
six. Reflection of exponents and swapping blocks takes X-cycles and
X-boundaries to Z-cycles and Z-boundaries, preserving weight and
nontriviality: `(f,g)` maps to `(bar(g),bar(f))`. Thus Z distance is
also six. A nontrivial mixed Pauli has a nontrivial X or Z component
whose support is contained in its support, so the ordinary Pauli
distance is six as well. ∎

This is an exact finite quotient result, not a new lower bound on all
logicals of the original `(12,12)` square. Its proof is a six-pattern
identity table and a local triangle argument, not a distance census.

## 4. An unrestricted-size inherited-patch obstruction at the first target

**Theorem 4.1.** In Q_4,

```
min_(q in Q_2) |w_4+alpha beta q(z^2,tau^2)(alpha,beta)| = 18.  (4.1)
```

*Proof.* The vector inside the minimum is exactly
`D_2(w_2+q(alpha,beta))`. Its source is nontrivial by the explicit
Q_2 dual in §1, hence has weight at least six by Theorem 3.3. Theorem
2.1 triples that weight. Taking q=0 attains 18. ∎

Thus any weight-14 solution of C(2) must use a genuinely new correction:
its boundary cannot equal an inherited boundary. More precisely its
multiplier h cannot lie in

`alpha beta {q(z^2,tau^2):q in Q_2} + (ann(alpha) intersect ann(beta))`.

Adding the common annihilator changes no physical boundary. This
qualification matters: a different-looking multiplier can represent
exactly the same inherited boundary.

The obstruction permits **arbitrary q**, not merely low-weight q,
small support, a chosen gauge, or a particular recursive stencil. It
does not exclude the remaining multipliers in the 192-dimensional
coefficient space of Q_4 and does not prove C(2) false.

## 5. Collective factorization identities and their exact physical costs

The following identities are valid polynomials before taking the
quotient; taking the quotient makes them cycles. Put

```
u=z, v=z^-1 tau, uv=tau,
p=u beta = z+1+z^-1 tau^3,
q=v alpha = z^-1 tau+1+z^-5 tau^2,
c=p+q = z+z^-5 tau^2+z^-1 tau+z^-1 tau^3.
```

Then `w_r=(u p^(r-1),v q^(r-1))`. The scalar c also has the useful
period-defect identity

`z^5 c=(1+z^6)+(1+tau^2)(1+z^4 tau)`.             (5.1)

In Q_r both defects have r-th power zero, so `c^r=0`. This explicitly
exhibits the two closure relations consumed by the motif.

### 5.1 Full polarization

Define the boundary multiplier with no division ambiguity:

`H_r=tau sum_(j=0)^(r-2) p^(r-2-j) q^j` for r>=2,
and `H_1=0`. Then

```
w_r+H_r(alpha,beta) = (u c^(r-1),v c^(r-1)).       (5.2)
```

*Proof.* For a power of two r, every binomial coefficient in the
`(r-1)`st row is odd. Consequently

```
p^(r-1)+c^(r-1) = q sum_(j=0)^(r-2) p^(r-2-j)q^j,
q^(r-1)+c^(r-1) = p sum_(j=0)^(r-2) p^(r-2-j)q^j.
```

Multiplying the first equation by u and the second by v, and using
`u q=tau alpha`, `v p=tau beta`, gives (5.2). ∎

At r=4, `H_4=tau(p^2+p q+q^2)` has 13 terms and each rewritten
block has 16 terms. The physical weight is **32**, versus the raw
18 and desired 14. The exact identity is useful for tracking the
closure algebra, but this scalar factorization is not a weight proof.
No infinite weight formula for (5.2) is inferred from its small tests.

### 5.2 Partial polarization is uniformly worse, with an exact formula

Let `r=2s>=2` be dyadic and put

`K_r=tau(p q)^(s-1)=tau^s(alpha beta)^(s-1)`.

Then the exact boundary identity is

```
w_r+K_r(alpha,beta)
   = c^s (z^s beta^(s-1),z^-s tau^s alpha^(s-1)).   (5.3)
```

Here the lower motif on the right is interpreted literally in Q_r;
it is not individually asserted to be a cycle at the larger periods.
For example the first-component difference in scaled coordinates is

`p^(2s-1)+p^(s-1)(p^s+q^s)=p^(s-1)q^s`,

which is exactly the first boundary component after multiplying by u.
The second component is identical with p,q interchanged. This proves
(5.3); the resulting pair is a cycle and retains the dual pairing one.

**Theorem 5.1.** The weight in (5.3) is exactly

`8*3^(a-1)` for `r=2^a>=2`.                        (5.4)

*Proof.* At r=2 the rewritten pair is `(u c,v c)`, with four distinct
terms in each block, hence weight eight. The right side of (5.3) is
exactly its `(a-1)`-fold Frobenius lift: c becomes `c^s` and the
remaining factors become `z^s beta^(s-1)` and
`z^-s tau^s alpha^(s-1)`. Apply the exact weight-tripling theorem.
Equivalently, its four coarse translates each contain `3^(a-1)`
distinct cells per block and cannot overlap. ∎

Thus partial polarization increases the raw `2*3^a` cost by precisely
a factor `4/3`, for every parameter, not only at the first test.

### 5.3 Complete r=4 boundary constants

For direct replay on Q_4, the nine-term multiplier for (5.3) is

```
K_4=tau^2 alpha beta
  = {(0,1),(0,2),(1,1),(18,6),(19,3),
     (20,3),(22,5),(23,2),(23,4)}.
```

The 13-term multiplier for (5.2) is its disjoint union with

`tau c^2 = {(2,1),(14,5),(22,3),(22,7)}`.

These are fixed collective patches, not outputs of a coset search.
Their resulting pairs are exactly

```
partial: (z^2 beta c^2, z^-2 tau^2 alpha c^2),  weights (12,12),
full:    (z c^3, z^-1 tau c^3),                 weights (16,16).
```

All coefficients and shifts are specified above, so the finite support
counts can be checked by simple expansion. Both pair with the explicit
U_4 in §1 to one, since their differences from w_4 are the displayed
boundaries. They are upper-bound witnesses of weights 24 and 32 in
Q_4, respectively, but are worse than the already banked weight 18.

The partial multiplier at this first target is already inherited:
`K_4=alpha beta (tau)^2`. Theorem 4.1 explains structurally why it
cannot compress. Adding the four fresh terms of the full polarization
does not rescue its physical cost.

## 6. Residue, controls, and computational accounting

The surviving first problem remains exactly

`exists h in Q_4: |w_4+h(alpha,beta)|<=14`,

with h outside the inherited-boundary class of §4. The dual is supplied
and stays odd automatically. No new nontriviality lemma is needed.
Neither the full nor partial polarization is a viable weight ansatz at
this parameter. A different fresh-scale patch may still work; its
existence has not been assumed, disproved, or inferred from samples.
No conclusion about the IBM r=3 target 30 follows from this note.

The primary new script is `scripts/a44_s2_codex_identities.py`. It verifies
fixed symbolic identities and the six local proof patterns, checks
their exact support counts, and checks the inherited-boundary identity
on prescribed inputs. There is no solver, representative enumeration,
or optimization objective. A development-only prescribed linear solve
was also tried once; it produced no improving representative, and that
unused routine was removed from the final diagnostic. No kernel or
coset alternatives were enumerated.

The final run passed in 0.001 seconds measured script time, with peak
RSS 16.5 MiB. All S2 diagnostics combined used well under one second
of measured execution, below the 180-second aggregate limit. The
script sets a 55-second alarm and CPU limit and checks peak RSS against
768 MiB; no output files are written. No old process, git state,
existing research note/script, memory, or Lean file was modified.

Reproduce from the repository root:

```
PYTHONDONTWRITEBYTECODE=1 uv run --offline --no-sync --no-cache \
  --project experiments/bb_lab python \
  experiments/bb_lab/scripts/a44_s2_codex_identities.py
```

## 7. Optional one-vector class gate: W4 survives every index-two cover

At the parent's request, the separate script
`scripts/a44_s2_codex_class_gate.py` tests just the banked original
W4 on the square `(24,24)`, not a population. Its support is obtained
from the §1 motif by `z=x,tau=xy` and the three translates
`(0,0),(8,8),(16,16)`; it has weight 54.

The covers use lattice generators `(L,0),(d,M)` and representatives
`0<=x<L,0<=y<M`. Reduction is
`(x,y)->((x-floor(y/M)d) mod L, y mod M)`. Projection to the base
is `(x mod 24,y mod 24)`, so the pullback is specified without any
choice of homology basis.

| cover | (L,M,d) | qubits | pullback weight | rank H_X | nontrivial |
|---|---|---:|---:|---:|---|
| x period doubled | (48,24,0) | 2304 | 108 | 1146 | yes |
| y period doubled | (24,48,0) | 2304 | 108 | 1146 | yes |
| diagonal | (48,24,24) | 2304 | 108 | 1146 | yes |

The checker builds the original fixed-pair syndrome and stabilizer
supports directly, checks that each pullback is a cycle, and produces
an odd Z-dual for each cover. Every displayed dual annihilates every
X-stabilizer generator and pairs with its pullback to one. No parity
law or all-parameter spectral statement is an input to this computation.

Thus the direct gate supplies **no rank-one warning**. Applying the
parent's stated all-parameter kernel/cylinder theorem would classify
W4 as mixed. That general theorem is being written and audited by the
parent; it is neither duplicated nor independently promoted to a new
theorem here. The unconditional conclusion here is the three nonzero
pullback classes. This supplies no weight floor at r=4 and does not
resolve C(2).

The first gate run took 0.0103 seconds and 19,251,200 bytes peak RSS;
a replay to record compact certificates took 0.0106 seconds. Both
remain within the original S2 time and memory budget. The script has
a 55-second alarm/CPU limit and a 768-MiB RSS check and writes no data
files.

### Complete binary dual constants

For a cover (L,M,d), bit `b*L*M+x*M+y` is physical qubit
`(block=b,x,y)` in its canonical representative. Each dual integer
is `core << shift`, with the core and shift below. This is a full
coordinate certificate, not merely a rank or an assertion of
nontriviality. All unused bits are zero.

### x-period doubled

Dual weight 496; shift 1128.

```
0xdb6db6db6db6000000db6db6db6db6000000db6db6db6db6000000db6db6db6db6000000db6db6db6db6000000db6db6db6db6000000db6db6db6db6000000db6db6db6db6000000db6db6db6db6000000db6db6db6db6000000db6db6db6db6000000db6db6db6db6000000db6db6db6db6000000db6db6db6db6000000db6db6db6db6000000000000000000b6db6d
```

### y-period doubled

Dual weight 216; shift 1129.

```
0x26c300a0000014080000000090d060140000825160140000904b6c16800080014014000092082d96d000820805945000904900268a0080080094000092012090d140820020825140904004904b6880000080014092080012082d820800020805b0490000490000080000080086012000012096002000002086c004000004940000000000b2036d00000020836d0000006db6db
```

### diagonal cover

Dual weight 430; shift 1129.

```
0x10000036db6d24db6d020000a69b6d36db6d120800b4d36d16926d000800a2da4d22db4d168004a2db6d200b6d94500030516da2db6d84d140208a2db01005940140320340a08368a6c300a2d36d260b0d025160a6900136da2d12082d34d368a6926d00080024da4db4db4d90400436db6da4d36d020800b6926d36d36d000120b6db4d36db690000008000000000006db6db
```

Reproduce the gate from the repository root:

```
PYTHONDONTWRITEBYTECODE=1 uv run --offline --no-sync --no-cache \
  --project experiments/bb_lab python \
  experiments/bb_lab/scripts/a44_s2_codex_class_gate.py
```

## 8. Conditional analytic review of the all-height Frobenius classification

Reviewed the parent's `A44_FROBENIUS_CLASS_PARITY.md` and the relevant
coordinate conventions in `A40_PARITY_TOR_PROOF.md`. This review is
**conditional on the base parity theorem**, including its regular
sequence, plane-module, natural period-Tor, and cylinder/cover-kernel
statements. It does not replace Parfit's independent review of that
theorem. No diagnostic computation was performed for this review.

**Verdict: no blocking gap in the extension.** Under those inputs, the
conclusion holds for every banked `W_r` with `r=2^a`, `a>=0`. The
following checks explain the chain-map and quantifier issues explicitly.

### 8.1 Flat base change and the quotient comparison map

For `R=F2[x^{+/-1},y^{+/-1}]`, Frobenius extension of scalars is flat:
R is free over its square subring with basis `1,x,y,xy`. Thus it takes
the regular Koszul resolutions for `(beta,alpha)` and `(F,G)` to
those for `(beta^2,alpha^2)` and `(F^2,G^2)`, where
`F=x^ell+1`, `G=y^m+1`. These last periods are exactly the doubled
ones. The coefficient map is Frobenius-semilinear over R and linear
over F2; it must not be interpreted as an R-linear squaring map on
the original modules.

The comparison from the squared plane resolution to the ordinary
plane resolution has components

`1`, `diag(beta,alpha)`, `alpha beta` in degrees 0,1,2.

Both identities are correct in the stated ordering:

```
(beta,alpha) diag(beta,alpha) = (beta^2,alpha^2),
diag(beta,alpha)(alpha^2,beta^2)^T
    = alpha beta (alpha,beta)^T.
```

It lifts the actual coefficient quotient
`R/(beta^2,alpha^2) -> R/(beta,alpha)`. Composing with Frobenius on
physical coefficients gives precisely the banked D, with beta on
the first block and alpha on the second. There is no missing block
swap or reciprocal in this calculation.

Balanced Tor is natural in the two module arguments. Computing it
instead with the doubled period resolution makes the comparison
above a coefficient quotient in every degree. On the surviving
ordinary component M0 the period differentials vanish, so the map
in **consistent, ordered period coordinates** is exactly
`(z1,z2)->(z1^2,z2^2)`.

A direct check of this coordinate formula is also available. Lift a
physical cycle to R and write its syndrome as

`beta f+alpha g=F A+G B`.

The balanced-Tor period coordinates are the ordinary plane classes
of A and B. The lifted D syndrome is
`F^2 A^2+G^2 B^2`, giving precisely their squares. Choices of lifts
and of the period decomposition disappear in the usual Tor quotient.
This check uses the natural coordinates from the base theorem, not
an arbitrary physical homology basis. Inversion and unit changes
used to identify the reciprocal ordinary plane algebra do not
obstruct Frobenius; they commute with squaring as algebra maps.

### 8.2 Reduced square image and the first-step qualification

Write `B=Fr(M0)={z^2:z in M0}` to avoid overloading `M0^2`, which
also denotes two period components elsewhere. The base calculation
gives `(x^2+x+1)^2=(y^2+y+1)^2=0` on M0. Consequently B is generated
by x^2 and y^2, each satisfying `T^2+T+1=0`. Their inverses are
polynomials in the same generators, so Laurent inverses add no
extra generators. There is therefore a surjection

`F4 tensor_(F2) F4 = F4 x F4 -> B`.

The reducedness conclusion is valid because this is a **finite
product of fields**, every quotient of which is another product of
fields. It is not an appeal to the false general claim that every
quotient of a reduced ring is reduced. Frobenius on the finite
reduced algebra B is injective and hence bijective. The same argument
applies after inversion to the reciprocal pair.

Thus, for components already in B, further squaring preserves zero,
nonzero, equality, and inequality. For arbitrary components in M0,
the first square can kill a nonzero nilpotent or a nonzero difference.
The parent's first-step qualification is necessary and correctly
stated. Nothing here proves D injective on the full homology.

### 8.3 Applying the gate and checking the quantifiers

Conditional on the base theorem, rank-one classes are the nonzero
classes on the two period axes or their diagonal. Componentwise
squaring sends each such subspace into itself, possibly to zero;
it cannot send a rank-one class to a mixed class.

The three explicit odd pairings in §7 imply that W4 is mixed under
the base cover-kernel identification. Since `D(W2)=W4` and
`D(W1)=W2`, the preceding observation already forces W2 and W1
to be mixed. W4 lies in B direct-sum B, so injectivity of Frobenius
on B preserves mixedness at every subsequent doubling. Equivalently
one may start stability at W2, which is already the image of W1.
Both versions justify exactly the claimed all-height conclusion.

All square periods here are `6*2^a`; the order-127 sector is absent
at every height. More generally it cannot switch on under a
power-of-two period cover, since 127 is odd. This does not extend
the conclusion to nondyadic W_r, nor assert first-step preservation
for every mixed class of an arbitrary member.

Accordingly the parent may integrate: **assuming the audited base
parity theorem, every W_(2^a) is mixed**. There is no class-type
reason from these three cylinder images to discard this witness
family. Its physical weights remain `6*3^a` on the original square;
no shortening, class minimum, pure-sector floor, or exact IBM
distance follows. The arbitrary fresh-scale correction in §6
remains the metric problem.
