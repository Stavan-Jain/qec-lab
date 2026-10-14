# Exact growing distance: IBM target and proved one-check modifications

Date: 2026-09-08. Branch: `claude/a40-tour-de-gross-column`, PR #27.

## 1. Success criterion and current scope

The primary target remains the fixed IBM Laurent pair

\[
A=1+y+x^3y^{-1},\qquad B=1+x+x^{-1}y^{-3}
\]

on `(ell,m)=(6(r+b),6r)`, `r>=1`, `b in {0,1}`, with exact distance
`6(2r+b-1)`. This passage is present in **v1, §2.7 item 3** of
[Tour de gross](https://arxiv.org/html/2506.03094v1#S2.SS7), checked directly
on 2026-09-08. A42 §1.3's sourcing doubt about its presence is resolved.

For a complete distance claim we require: specified checks for every
parameter, commutation, logical dimension with exceptions, a lower bound
for every nontrivial logical in both CSS sectors, and a nontrivial logical
attaining the same value. Finite censuses and upper-bound witnesses alone
do not settle an infinite family. Neither does a fixed-distance thin strip.

Preference order: original IBM family; an explicit infinite subfamily with
both periods growing; an explicit related BB modification. A modification
must list its changed checks and parameters. Merely changing the name of
an unresolved inequality is not a proof step.

**Current exact result.** The allowed modified-family target now has
a complete analytic proof: `IBM_MODIFIED_EXACT_DISTANCE_THEOREM.md`
proves [[72r^2,12,4r]] for every r>=1 on the original square frames.
It keeps P=x^3+y+y^2 and replaces only the second physical check by
Q=(1+x+x^2)(y+x+1). The ordinary primary algebra, including its
doubled tangent component, is retained; check weight becomes eight
instead of six. This is not the unchanged IBM conjecture, whose
distance coefficient is not retained. The full theorem, explicit
attainer, scope, and independent reviews are recorded in §13.

**Existing dimension exception.** A42 §1.3 gives `k=12` for all `b=1` and
`k=12+14*[127 divides r]` for `b=0`. Thus a uniform `[[n,12,d]]` claim for
the entire square column already needs correction. Additional homology
does not by itself refute the proposed distance; its weights require a
separate argument. The `762 x 762`, `k=26` phenomenon is also explicit in
[Liang et al., Example 4 and §II.2](https://arxiv.org/html/2503.03827v2).

## 2. Three renewed proof tracks

Fresh Codex agents take over the mathematical work; the old headless
processes are not resumed. All new diagnostics have short time and memory
limits. Their write scopes are disjoint:

| Track | Immediate question | New note |
|---|---|---|
| A43-S2 | Audit L3's quantifiers and the distinction between representative geometry and logical class; find a valid toroidal bound | `A43_S2_codex_audit.md` |
| A42-S6 | Prove a substantive consequence of HM / seam doubling, checking both orientations and exceptional periods | `A42_S6_codex_residual.md` |
| A44-S1 | Find a symbolic b=0 attaining family or a precise obstruction; no repeated SAT hunts | `A44_S1_codex_family.md` |

These are work assignments, not claims of success. In particular A43's
all-rectangles L3 is stronger than IBM's conjecture and must survive
adversarial checks before being used as the program's reduction. The
previous comparison/surgery approaches cannot be retired solely because
a stronger unproved claim has been proposed.

**First audit result (verified).** `a43_s2_codex_audit.py` constructs
a weight-48 nontrivial cycle on `(54,6)` with cyclic gaps `(3,1)` and a
weight-66 example on `(78,6)`, disproving the all-rectangles L3 floors
54 and 78 respectively. Nontriviality is checked by both anticommuting
dual witnesses and a rank test. Each example retains an injective
cylinder lift. Thus it is the gap-dense version that is refuted;
neither IBM's b-bit family nor the restricted no-cylinder-lift statement
is refuted by these examples. Parent replay passed in 0.0065 seconds.

The audited sufficient replacement is **C_min**: on the actual IBM
rectangles, a connected, class-minimal logical without either axis-cylinder
lift has weight at least `ell+m-6`. A globally minimum nontrivial logical
is connected, so C_min plus the two applicable cylinder bounds suffices.
This remains an open lemma. Connectivity does not follow merely from
minimality in a prescribed logical class. The A43 audit also provides a
proved standard comparison family `[[72r(r+b),8,4r]]` on the original
rectangles using `A=1+x+x^2, B=1+y+y^2`. Theorem S below instead preserves
`k=12` and changes the aspect ratio.

**Cylinder audit (independently replayed).** The new
`A42_S6_codex_residual.md` proves an exact binary-tower graph identity,
an infinite counterfamily to arbitrary-support syndrome expansion, and
a restricted-support mixed bound. For its explicit omega cycle with
`|S|=3q`, arbitrary binary boundary mixing has physical weight at least
`4q`. A follow-up proposed by the parent and proved in the cylinder
note combines two compressed checks into a fixed 16-state, 128-edge
potential certificate. It proves the exact coset distance
`d(1_S,B_q)=q` for every `3` not dividing `q`. The singleton-fibre
stratum for this specified support then has exact physical minimum
`6q`, with an attaining boundary supplied. Other fibre occupancies
and omega supports remain open. The parent replay passed in 0.1
seconds; the infinite quantifiers follow from the potential inequality
summed on an arbitrary-length closed walk, not finite-period evidence.

The audit also corrects the proof obligations. Original-lane full
cylinder certificates at periods `9,12` do not automatically cover the
opposite orientation: the displayed rotation maps `(A,B)` to
`(B,bar A)`, not to itself. Both lanes were independently checked only
at the earlier periods `3,6`. Strict order-three primary purity and
zero three-sheet binary fold differ when the odd part grows. Finally,
the binary boundary identification needs an order-127 qualification,
and nine-sheet fibres alone do not cover periods divisible by `27`.
These are scope corrections, not conjecture counterexamples. General
HM, the pure tower, the opposite orientation, higher odd-factor mixing,
and the toroidal lower bound remain distinct obligations.

**Square upper-bound audit (independently replayed).**
`A44_S1_codex_family.md` supplies finite polynomial identities proving
the `12r` L12 upper bound on every `(6q,6r)` frame, including `b=0`.
At `(18,18)`, all nonzero one-block X-cycles have weight at least 54,
with a nontrivial one attaining 54: a conjectured weight-30 attainer
must therefore use both blocks. An explicit Frobenius family on
`r=2^a` has nontriviality proved by a formal-series pairing, but its
exact witness weight `6*3^a` misses the desired `12r-6` for `a>=2`.
The proposed correction is an open existence lemma in a specified
class, not a new distance theorem. The parent diagnostic passed in
0.006 seconds. This audit also independently checked Theorem S,
S-opt, and Proposition C below and found no proof gap.

## 3. A complete exact theorem for a related BB family

This section proves a comparison family using standard product-code
ideas. It is **not a proof of the original IBM conjecture, nor a claim of
a new code construction**. It retains 12 logical qubits, two qubits per
lattice cell, fixed trinomials, connected checks, weight-six generators,
and distance proportional to the linear size. Its overhead is worse.

### Theorem S

Let `s,u>=1`. On `Z_(3s) x Z_(7u)` use

\[
A_{\rm sep}=1+x+x^2,\qquad B_{\rm sep}=1+y+y^3,
\]

with the lab convention `H_X=[M_A|M_B]`, `H_Z=[M_B^T|M_A^T]`.
Then the code has exact parameters

\[
\boxed{[[42su,\;12,\;\min(2s,4u)]]}.
\]

Both `d_X` and `d_Z` equal this minimum. In particular `s=2t`, `u=t`
gives, for **every integer `t>=1`**,

\[
\boxed{[[84t^2,\;12,\;4t]]}.
\]

### 3.1 Classical factors

For length `3s`, the kernel of `1+S+S^2`, where `S` is the cyclic shift,
has dimension two. The recurrence is, up to reversing the sequence,
`c_(j+2)=c_j+c_(j+1)`. Hence `c_(j+3)=c_j`. Every nonzero initial pair
gives a period-three word with two ones: `110`, `101`, or `011`.
Thus its parameters are `[3s,2,2s]`.

For length `7u`, the kernel of `1+S+S^3` has dimension three. Its recurrence
is, up to reversal, `c_(j+3)=c_j+c_(j+1)`. The companion polynomial
`z^3+z+1` is irreducible over F2 because it has no root in F2. Its roots
are nonidentity elements of F8*, whose order is seven. Therefore the
companion matrix has order seven and acts transitively on the seven
nonzero three-bit states. Of these states exactly four have first
coordinate one. Every nonzero sequence consequently has four ones per
period seven, giving `[7u,3,4u]`.

Both recurrences are reversible and admit arbitrary initial states, so
there are no missing solutions. They apply unchanged when `s` or `u`
is even: no semisimplicity of the larger cyclic group algebra is assumed.
The transpose of either circulant is conjugate to it by coordinate
reversal. Transpose kernels therefore have the same dimensions and
distances.

### 3.2 Product complex and dimension

For clarity use two square linear maps `P:E->V` and `Q:F->W`, each with
specified coordinate bases. The degree-one space of their tensor complex
is `(E tensor W) direct-sum (V tensor F)`, with

\[
\partial_1(a,b)=(P\otimes I)a+(I\otimes Q)b,
\quad
\partial_2 h=((I\otimes Q)h,(P\otimes I)h).
\]

Over F2, `partial_1 partial_2=0`. Splitting each finite-dimensional
two-term complex into its kernel, cokernel, and an isomorphism summand
gives

\[
H_1\cong(\ker P\otimes\operatorname{coker}Q)
\oplus(\operatorname{coker}P\otimes\ker Q).
\]

Take `P=(1+S_(3s)+S_(3s)^2)^T` and
`Q=(1+S_(7u)+S_(7u)^3)^T`. The BB complex above is this complex after
swapping its two qubit blocks. This is a coordinate permutation, so
physical Hamming weight is unchanged. The two kernel dimensions are two
and three; since the maps are square, the cokernel dimensions match.
Thus `k=2*3+2*3=12` and `n=2*(3s)*(7u)=42su`.

### 3.3 Lower bound without quantum enumeration

We prove that every non-boundary cycle `(a,b)` has weight at least
`min(d(ker P),d(ker Q))`.

For any functional `lambda in W*` annihilating `im Q`, contract `a`:

\[
a_\lambda=(I\otimes\lambda)a\in\ker P.
\]

Membership follows by contracting the cycle equation. Hamming weight
cannot increase under this contraction: each nonzero output coordinate
requires a nonzero input coordinate in a distinct row. If any
`a_lambda` is nonzero, then `|a|>=d(ker P)`.

Likewise any `mu in V*` annihilating `im P` produces

\[
b_\mu=(\mu\otimes I)b\in\ker Q,
\]

of weight at most `|b|`. If any is nonzero, then `|b|>=d(ker Q)`.

It remains to show that if all contractions of both types vanish, the
cycle is a boundary. Vanishing of the first type gives
`a in E tensor im Q`, so choose `h` with `(I tensor Q)h=a`.
Subtract `partial_2 h`, leaving `(0,b')` with
`b' in V tensor ker Q`. Contractions of `b'` against every `mu` still
vanish, because `mu P=0`. Thus
`b' in im P tensor ker Q`. Lift its first factor to obtain
`h' in E tensor ker Q` with `(P tensor I)h'=b'`.
Then `partial_2 h'=(0,b')`, proving the assertion.

Consequently `d_X>=min(2s,4u)`. Apply the same argument to the transposed
complex for `d_Z`. Its two classical kernel distances are unchanged,
so `d_Z>=min(2s,4u)` as well. A general Pauli outside the stabilizer has
a nontrivial X or Z component, whose support is contained in its support;
therefore this also bounds the Pauli distance.

### 3.4 Matching upper bounds and nontriviality

Choose a minimum word `c in ker P`, and a coordinate vector `e_j` whose
class in `coker Q` is nonzero. Such a coordinate exists because the
cokernel is nonzero and coordinate vectors span `W`. Then
`(c tensor e_j,0)` is a cycle of weight `d(ker P)`. A functional
annihilating `im Q` but taking value one on `e_j` contracts this cycle
to `c`, whereas every boundary contracts to zero. Thus it is nontrivial.

Similarly `(0,e_i tensor d)` with `d` a minimum nonzero word of `ker Q`
and `e_i` nonzero in `coker P` is nontrivial of weight `d(ker Q)`.
This proves `d_X<=min(2s,4u)`; the transpose construction proves the same
for `d_Z`. Combined with §3.3, equality follows for every parameter.

### 3.5 Locality, connectedness, and comparison

All three monomials within each block are distinct for these periods,
so each X and Z generator has weight exactly six. The largest displacement
is three cells, independent of the parameters.

The graph joining qubits that share an X check is connected: the constant
terms join the two blocks at every cell; the `1+x` terms connect adjacent
x-cells in one block; the `1+y` terms connect adjacent y-cells in the
other block. Combining these moves reaches every qubit. This is not a
growing number of disconnected copies of a fixed-distance code.

The modification replaces IBM's coupled two-variable polynomials with
separate one-variable factors and changes the periods. At `s=2t,u=t`,
`kd^2/n=16/7`, versus the conjectured IBM `b=1` asymptotic value 24.
At the same asymptotically large distance and `k=12`, its physical-qubit
overhead is therefore 10.5 times IBM's conjectured overhead. The distinction
is material; this theorem establishes an exact growing BB comparison
family, not the desired IBM performance or an original-family subcase.

### 3.6 Validation and provenance

`scripts/ibm_separable_family_check.py` checks the two classical kernels
and transpose kernels for repetition counts `1,2,3,4,8`, and the quantum
ranks, weight-six checks, and explicit anticommuting X/Z witnesses at
`t=1,2,3`. It first runs `validate_banked`, enforces a 60-second time limit
and a 1-GiB RSS guard, and writes
`data/exact_family/separable_checks.json`. These are convention checks;
the universal proof is §§3.1–3.5.

The general product-code distance mechanism is standard:
[Tillich–Zémor, Theorem 15](https://arxiv.org/html/0903.0566v2#S5).
Separable BB codes as products of cyclic codes are discussed explicitly
by [Eberhardt–Pereira–Steffan](https://arxiv.org/html/2412.04181v1).
The proof here specializes that framework and includes its own contraction
argument and attaining representatives. No novelty claim is made.

### 3.7 The separable weight-six, k=12 limitation is exact

**Proposition S-opt.** For any separable BB code built from two binary
univariate trinomials, with square circulant check matrices, nonzero
logical dimension `k=12`, and distance `d`,

\[
 n\ge\frac{21}{4}d^2.
\]

Theorem S with `s=2t,u=t` attains equality. Therefore changing the
univariate trinomials within this class cannot recover IBM's conjectured
asymptotic relation `n~d^2/2` at `k=12`.

**Proof.** Let the two seed lengths be `n_1,n_2`, kernel dimensions
`k_1,k_2`, and kernel distances `d_1,d_2`. Transposition of a circulant
is conjugation by coordinate reversal, so the square product has
`k=2k_1k_2=12` and exact `d=min(d_1,d_2)`.

Neither positive `k_i` can equal one. The kernel dimension is the
degree of the gcd of the check polynomial with `z^(n_i)-1`.
A degree-one gcd over F2, with nonzero constant term, would be `z+1`.
But a binary trinomial takes value one at `z=1`, so it has no such
factor. This remains true after cancellation of coincident cyclic
exponents, because evaluation at one is still one. Hence the only
possibility is `{k_1,k_2}={2,3}`.

For any binary linear `[n_i,k_i,d_i]` code, summing weights over its
nonzero codewords gives at most `2^(k_i-1)n_i`: each active coordinate
is one on exactly half the codewords, while an inactive coordinate
contributes zero. The same sum is at least `(2^k_i-1)d_i`. Thus the
dimension-two seed obeys `n_1>=3d_1/2` and the dimension-three seed
obeys `n_2>=7d_2/4`, after interchanging the seeds if needed. Therefore

\[
n=2n_1n_2\ge\frac{21}{4}d_1d_2\ge\frac{21}{4}d^2.
\]

For the balanced family, `n=84t^2=(21/4)(4t)^2`, proving sharpness.
This is a restriction on the specified separable subclass, not a bound
on general coupled BB codes. It makes preserving the coupled checks
essential if the goal includes IBM-like overhead.

### 3.8 An exact gross-derived benchmark if bounded check weight is dropped

There is a much more direct exact family if one relaxes the fixed-weight
BB requirement rather than changing the bulk trinomials. It is useful
for identifying which part of the goal is difficult, not for declaring
the IBM conjecture solved.

**Proposition C.** Concatenate the banked gross `[[144,12,12]]` code,
whose two CSS distances are both 12, with the square generalized Shor
code on `t` rows and `t` columns. For every integer `t>=1`, this gives

\[
\boxed{[[144t^2,12,12t]]},\qquad d_X=d_Z=12t.
\]

**Inner code proof.** Generate Z checks by adjacent pairs within each
row, and X checks by the union of each adjacent pair of full rows.
They commute, and their independent ranks are `t(t-1)` and `t-1`.
Thus the inner code encodes one qubit. An X operator commuting with
the Z checks is constant along each row; it is nontrivial precisely
when an odd number of full rows is occupied. A Z operator commuting
with the X checks has the same parity in every row; it is nontrivial
precisely when that parity is odd. Minimum representatives are one
full X row and one Z coordinate in every row, each of weight `t`.
For any mixed Pauli, a nontrivial X or Z component also has weight
at least `t` and its support lies in the Pauli's support. Therefore
the inner distance is exactly `t`, including when `t` is even.

**Concatenation proof.** An operator commuting with the concatenated
checks restricts on each inner block to an inner normalizer. Quotient
by the inner stabilizers to obtain an effective Pauli on the 144 gross
qubits. It is a nontrivial outer logical whenever the original operator
is a nontrivial concatenated logical. At least 12 disjoint inner blocks
therefore carry nontrivial inner logicals, and each costs at least `t`
physical qubits. This proves the lower bound `12t`. Encode a banked
weight-12 pure X gross logical using the inner weight-`t` X representative
on each occupied block. It is nontrivial and has weight `12t`. Repeat
for Z. The code space is the outer code embedded by 144 inner encoders,
so it still encodes 12 qubits (code-space dimension `2^12`). This is the standard concatenation
mechanism, not a novelty claim; see, for example,
[Grassl et al.](https://arxiv.org/abs/0901.1319) for its broader setting.

The displayed generating set has inner weights `2,2t` and encoded gross
check weight `6t`. More importantly, this is **intrinsically not an LDPC
family**, not merely a poor choice of generators: any stabilizer whose
effective outer action is nonidentity has a nontrivial inner logical in
some block, hence weight at least `t`. A generating set for the full
stabilizer group must contain such an element, since inner stabilizers
alone do not generate the encoded outer constraints.

At equal large distance and `k=12`, this family uses twice the conjectured
IBM asymptotic number of qubits, compared with Theorem S's factor 10.5.
It starts from the actual gross code at `t=1`, but discards an essential
architectural advantage of BB codes. It is consequently a second
benchmark, **not an accepted replacement for the active fixed-weight
IBM proof target**. Any proposed weight-reduction gadget must separately
prove exact distance preservation, not just a constant-factor bound.

## 4. Structurally different routes for the original pair

### 4.1 Toric normal form with the physical weight retained

[Haah's classification](https://arxiv.org/abs/1812.11193) and
[Liang et al.'s twisted-torus analysis](https://arxiv.org/html/2503.03827v2)
motivate identifying the string charges and translation action first.
This can expose the correct closure arithmetic without dividing every
representative into cylinder and mixed sectors.

However, a local Clifford circuit generally changes Pauli weight. A
topological equivalence, an ideal equality, or a computation of `k` does
not establish the distance constant. On finite tori, coarse-graining and
twists also require explicit compatibility with the periods. In
particular the order-127 species cannot simply be discarded because a
small example does not see it.

Concrete proof target: construct an explicit normal-form map and retain
the induced physical weight function, including stabilizer dressing.
Then prove a weight-nonincreasing reduction of a minimum logical to a
finite list of string/network types and solve their winding problem.
Interactions between strings and cancellations must be charged exactly;
the number of toric layers alone is insufficient. This reduction is
currently open, not a supplied theorem.

There is now a proved obstruction to its simplest version:
`IBM_GRAPHIC_MINOR_OBSTRUCTION.md` exhibits a fixed 42-qubit F7 minor
of the original syndrome column matroid, embedding in every IBM
member with `r>=2`. Thus row-basis changes and qubit permutations
alone cannot turn that matrix into ordinary graph incidence. The
fixed certificate verifies all 128 subset ranks in under a millisecond.
This does not rule out Clifford maps or a metric-aware network normal
form; it shows why their extra weight bookkeeping cannot be skipped.

### 4.2 Exact transfer recurrence after a uniform normal-form theorem

If minimizers can be represented by bounded-width strings and finitely
many junction types, with all logical and winding labels retained, their
cost may be encoded by a finite min-plus automaton independent of `r`.
Under suitable irreducibility and closure hypotheses, powers of a fixed
min-plus matrix have eventual affine-periodic behavior. See
[Kennedy-Cochran-Patrick et al.](https://arxiv.org/abs/2005.05390).

This is a route to an **exact** law, including exceptional residue classes,
by proving the finite transient and exhibiting attaining paths. It is
not justified by a transfer matrix whose width grows with the torus. The
missing step is completeness of a parameter-independent state model;
the earlier exponential races do not supply it. A useful initial target
is one arithmetic subfamily with both periods unbounded, rather than a
uniform all-rectangle inequality.

### 4.3 A weight-compatible chain retraction

Section 3 succeeds because tensor contractions detect every homology
class and never increase support. Seek an analogous finite family of
maps from the original BB complex to small classical complexes, together
with an explicit account of which classes each detects and how their
weight budgets combine. Unlike spectral projection alone, this must
track physical coordinates. A filtration, triangulation, or local basis
change is useful only if its support inequalities are proved.

The A43 cap applies to its particular fractional packing of
anticommuting logical supports. It does not by itself exclude every
linear-programming or chain-map proof. Nor is substituting the desired
distance into that cap an unconditional numerical obstruction. The
contracting-map argument of §3 illustrates a different certificate shape.

## 5. Next milestones

1. Build on the proved 16-state certificate and the spatial uniqueness
   theorem, rather than rerunning growing censuses. For the canonical
   support, let k count holes and b_out count outside error sites.
   The remaining target is `k+3b_out>=3q`. The k=0 stratum is exact
   and unique; k=1,2,3 meet the needed floor at every `3 ∤ q`,
   and §11 extends this to k<=6 at dyadic heights q>=2.
   A six-column generator-strip certificate is also available, but
   width completeness is not proved. Any extension to other omega
   supports must state and prove its additional quantifiers.
2. Prefer a complete theorem for the original `b=1` column or an explicit
   infinite subfamily, using its existing matching upper bound.
3. For `b=0`, keep the order-127 homology exception and the failed band
   ansatz explicit. Require an attaining construction as well as a floor.
4. Use Theorem S as an exact, inexpensive baseline for evaluating new
   proof mechanisms. Improving its parameters or returning to IBM's
   fixed pair is additional work; its existence does not close the IBM
   research goal.

## 6. Earlier continuation checkpoint (2026-09-08)

The preceding goal turn made concrete progress: a verified infinite
counterexample family replaced a false lemma, and Theorem S supplied a
complete exact comparison theorem. The original IBM goal remains open.

The completed A43 audit was followed by a closer modification:
`B_plus=1+x+x^(-1)y^3`, keeping `A=1+y+x^3y^(-1)` and changing only
the sign of one exponent. The unit-cleared pair is the standard
`x^3+y+y^2`, `y^3+x+x^2`. It agrees with the original pair on the
period-six direction. `A43_S3_codex_sign_variant.md` now proves the
dimension law on `(12t,6t)`: 12 encoded qubits for odd `t`, 16 for even
`t`, together with `2t<=d_t<=12t`. Its full translation-period lattice
is `12 Z^2`. The different local multiplicities rule out direct
monomial rewrites as a separable product or a reciprocal-trinomial
color code. The parent replay passed in 0.0078 seconds. Exact distance
for this closer variant is still open, and the original minus-sign
code's distance records cannot be imported.

S-opt closes the prospect of recovering IBM's overhead merely by
optimizing separable weight-six, k=12 factors. Proposition C provides
an exact gross-derived family only by dropping bounded check weight.
The local F7 certificate separately excludes an ordinary graph-incidence
rewrite on the unchanged original qubits. These are specific limits on
shortcuts, not impossibility results for proving IBM's distance.

The three initial agent assignments and the sign-variant follow-up have
returned. The witness agent independently verified the fixed F7 minor,
including its embedding on the `12x12` torus. The cylinder equality
audit found a zero-slack five-cycle that can be padded to length eight:
the compressed graph does not force unique coset leaders. This does
not produce another physical leader, since lifting to the unbounded
syndrome has not been proved. The next weighted-stability argument
must retain additional x-location information; do not reuse the
refuted graph-only uniqueness premise. No old frozen process has
been resumed. Research changes
remain in the working tree; no new commit or push was made in this pass.
The exact fixed-weight IBM goal remains active and incomplete. The
later results in §7 supersede the open physical-uniqueness question
at this checkpoint, but do not revive graph-only uniqueness.

## 7. Spatial, algebraic, and exact-product continuation

### 7.1 The original canonical leader is unique at all heights

`A42_S7_spatial_equality.md` supplies the missing Laurent-coordinate
argument. For every `q>=1` with `3 ∤ q`, a compact solution of

```
(x^3+t+t^2) f + (1+(x+x^2)t^3) g = N_q
```

has weight at least q, and equality forces
`(f,g)=(x^-3 N_q,0)`. The nontrivial equality graph components
are retained. A six-row spatial context around its D-to-C edge
forces one monomial to be both `x^-3` and 1, and the separate
g-only component is excluded by its actual x-extrema. This proves
uniqueness without enumerating heights or assuming compactness of
an arbitrary finite graph walk.

For the fixed omega support from S6, noncanonical errors therefore
have weight at least q+2. The physical floor for k=1,2,3 holes is
respectively `6q+4`, `6q+2`, `6q`; the unrestricted fixed-support
floor improves to `4q+2`. The exact k=0 minimizer is unique.
The parent's independent proof read and fixed-certificate replay
passed; these are infinite-quantifier results, not census evidence.
They do not establish general HM or the pure tower floor.

The parent also gives a separate opposite-orientation argument in
`A42_MIRROR_CANONICAL.md`. It uses the actual rotated pair, clears
it to `P=x^3(1+t)+t^2`, `Q=1+(x+x^2)t^3`, and proves the
corresponding canonical leader `(x^3 N_q,0)` is unique by retaining
one additional syndrome row. It does not identify the original
code with its rotation. Its fixed identities replay successfully;
the independent review in the S7 follow-ups passes as well.
This is only the analogous canonical-support theorem, not the full
opposite-orientation cylinder floor.

### 7.2 Sharp canonical stability remains spatial

`A42_S8_canonical_discharge.md` proves canonical HM for arbitrary
row masks when the boundary generator occupies columns `-3,-2`.
An exact finite integer potential certifies the larger actual strip
`[-3,2]` for **every height**. There is no theorem reducing arbitrary
compact generators to that strip without increasing energy.

S8 also refutes the stronger folded relaxation: after identifying
x modulo six, an explicit boundary at q=512 has inside weight
1207 and outside weight 381, violating `3 outside>=inside` by 64.
Its actual unbounded-cylinder interpretation has positive energy
5466, so this is not an IBM or original-cylinder counterexample.
The finite row table and seam identities, not a failed optimizer,
certify the obstruction. The parent replay passed in 0.0005 seconds.

Spatial monomial-run lemmas bound clean noncanonical excursions,
but their sharp overlapping-defect charge is not proved. Taking
the maximum of the raw-syndrome and full-support triangle bounds
also fails on an explicit alternating-mask family. Neither folded
periodicity nor those two global scalar bounds should be reused as
a claimed proof of full HM.

### 7.3 Original square witnesses need fresh-scale corrections

`A44_S2_symbolic_compression.md` proves the Frobenius lift triples
the weight of **every** input vector, including an arbitrarily
dressed one. A six-pattern local argument proves the finite
quotient Q_2 has exact distance six. Hence every correction at
r=4 inherited through this lift has weight at least 18 in Q_4,
not the desired 14. This excludes arbitrary-sized inherited patches,
not only the previously tested small patches. Other corrections
remain possible; no target attainer or counterexample was found.

The explicit original-square W4 witness survives all three index-two
covers, with an odd dual certificate for each. The parent replay
passed in 0.0105 seconds. Combined with the all-size class theorem
below, this places W4 in the mixed sector, so this construction is
not merely pursuing a rank-one class by mistake. It supplies no
new lower bound for that class. The symbolic identity replay passed
in 0.001 seconds, without a representative census.

`A44_FROBENIUS_CLASS_PARITY.md` upgrades this to every dyadic
witness W_(2^a). An independently reviewed Frobenius/Tor argument
shows that the induced map squares the two period components;
after its first application the square-image algebra is reduced.
A further direct identity identifies the components as
`[C^r],[D0^r]`, with
`C=x^-3 y^2+x^-1 y^4+y^7`, `D0=x+y`. At the ordinary point
`(omega,omega^2)` they evaluate to `(omega^(2r),1)`, nonzero and
unequal for every r=2^a. Thus the whole banked dyadic witness
family is mixed, by an analytic certificate that does not require
the finite W4 gate. Its weights remain `6*3^a`; a matching
conjectured upper bound has not been obtained.

### 7.4 All-size parity from a fixed Smith calculation and Tor

`A40_PARITY_TOR_PROOF.md` gives a self-contained algebraic upgrade
of the S11 class-space calculations. Reducing Q modulo P gives a
two-by-two multiplication matrix with Smith invariants
`g,g^2 h`, where `g=x^2+x+1` and
`h=x^7+x^6+x^4+x^2+1` is irreducible. Thus the plane module splits
into a six-dimensional part killed by both sixth powers minus one
and a seven-dimensional F128 part. This proves the module input
without local-multiplicity stabilization or large member matrices.

Naturality of the period Koszul resolution then shows, for **all**
rectangles with both periods divisible by six, that the image of a
primitive `(a ell,b m)` cylinder is exactly `{(a z,b z)}` in
`M_Lambda^2`. The three parity spaces are therefore the full list,
and the corresponding double-cover pullback kernels are exactly
those spaces. This includes the order-127 exceptional sector.
The ordinary counts remain 189 liftable and 3,906 mixed nonzero
classes; the exception has 24,573 liftable and 67,084,290 mixed.

The fixed algebraic diagnostic passed in 0.0003 seconds. The proof
passed the cylinder agent's independent audit, recorded in S7 §10,
including a separate 0.014-second polynomial check. No Tor basis is claimed
to preserve physical weight. Class liftability does not guarantee
an injective support lift of a minimum representative, so this does
not discharge C_min or prove the desired mixed-class floor.

### 7.5 An exact, gross-derived LDPC family, with a changed architecture

`GROSS_SURFACE_PRODUCT.md` proves, for every integer t>=1,

```
[[144(3t^2-3t+1), 12, 12t]],       d_X=d_Z=12t.
```

It starts with the actual full gross check complex, including all
check dependencies, and tensors it with an open repetition interval
and its reversed dual. The proof tracks every homology class:
one interval multiplies its minimum by t, while augmentation of the
other is a weight-contracting homotopy equivalence. Reversing the
full complex gives the same argument in the other CSS sector.
Matching logical representatives are explicit. The displayed
generators have weight at most eight, bounded qubit participation,
and a connected incidence graph; at t=1 it is the original gross
code. Independent reviews by both cylinder agents found no gap in
the proof, conditional on the banked seed distance 12.

This is an augmented-surface construction, not a newly discovered
general product-distance rule. It avoids Proposition C's growing
check weights but changes the BB architecture and uses asymptotically
six times IBM's qubits at the same k and d. It is an exact fallback,
**not a proof of IBM's conjecture or a silent replacement of the
original fixed-trinomial goal**. The finite product diagnostic checks
commutation, dimension, sparsity, connectivity, and attaining logicals
at t=1,2,3 in under one second; the universal distance proof is the
chain argument, not these finite checks.

The central original-family obstruction remains a sharp metric
inequality for mixed logicals and the remaining pure-cylinder
classes. The all-size algebraic classification, restricted canonical
theorems, and exact fallback do not make that inequality automatic.
No old paused process was resumed, no Lean work was started, and
no changes from this continuation were committed or pushed.

## 8. Metric-proof continuation: what the new reductions do and do not preserve

The previous continuation was substantive progress: the all-size parity
theorem and original/mirror canonical uniqueness are proved, and the
dyadic witnesses have explicit mixed fluxes. This continuation targets
the remaining metric gap rather than treating the exact alternative as
completion of the original goal.

### 8.1 Canonical tails: safe rules with sharp costs

`A42_S9_tail_elimination.md` gives exact four-mask identities for
deleting an extreme boundary-generator column. The worst increases
in canonical energy `E=3 outside-inside` are at most
`6|a intersect t a|` on the left and `3|(t+t^2)a|` on the right;
explicit all-even-height examples attain both constants.

Thus arbitrary-length left tails whose masks have no cyclic adjacent
ones and right tails whose masks are row-constant can be removed
without increasing E. The resulting six-column generator lies in
the already certified strip, so canonical HM holds for this genuinely
unbounded-span subclass. More general context-sensitive deletion
guards are exact, and stop if they fail.

Unrestricted clipping is not valid: one explicit generator has
energy 11q, clipping gives 17q, every single-column move into the
strip also raises energy, and translating the whole generator gives
15q. A correlated replacement preserving every inside bit instead
gives 8q. All these formulas hold for every even q by a proved
two-row repetition identity. The parent fixed replay passes in
0.0001 seconds. Full HM still needs coverage of stalled general
tails by correlated safe moves, or a sharp interface-charge bound.

### 8.2 Pure cycles: normalization is not a free metric step

`A42_S10_relative_filtration.md` proves the localized omega plane
module has DVR invariants u and u^2. Transfer from height two is
an isomorphism on homology at every higher binary height and
multiplies the transferred representative's weight by q/2.
The exact homology filtration has only its last two nonzero layers.

Nevertheless the unique normalization obtained by removing only
lower-layer generator coefficients can increase physical weight.
For every dyadic q>=4, an explicit nontrivial maximal-filtration
class has a representative of weight 3q+3 whose forced normalization
weighs 9q/2. The same class has another representative of weight
3q. This is an obstruction to that normalization rule, **not** a
counterexample to the desired pure floor.

The note proves classical weight retention in the physical basis,
but shows precisely why substituting a class's maximal valuation for
the actual vector's valuation is invalid. It also gives a separate
two-coordinate relative-code counterexample to inferring coset
distance solely from torsion multiplicity. The parent fixed original
q=4 replay passes in 0.00015 seconds. Optimizing the remaining
top-layer boundary, or proving the floor without normalization,
remains necessary.

### 8.3 A linear recovery certificate must have growing spatial width

`IBM_CANONICAL_RECOVERY_OBSTRUCTION.md` examines a different
sufficient mechanism: linearly recover all 3q inside boundary bits
from outside bits, with each input used in at most three outputs.
Such a map would prove HM. If its inputs occupy K block-columns
and its maximum column load is c, the note proves the universal
necessary tradeoff

```
K c >= 3 ceil(q/3)-3.
```

Every systematic recovery cochain must wind around the row period:
a gap of three rows would unwrap it to a plane syndrome-row
combination, contradicted by an explicit period-three cycle supported
entirely inside the canonical strip. Counting the decoder entries
then proves the inequality, without translation-equivariance or a
solver. In particular no fixed-width linear recovery can certify
HM at all heights with load three.

An explicit ten-block-column recovery is algebraically correct but
one of its columns already has load 2q+1 or 2q-1. This separates
exact recovery from a sharp support bound. The theorem passed an
independent review in S9 §9; the parent fixed identities replay in
0.0004 seconds. Growing-width analytic maps and nonlinear discharge
remain possible. This input-window obstruction does not contradict
the generator-tail theorem in §8.1.

### 8.4 Fixed diagonal seams cannot supply the square discount

`A44_S3_linear_seams.md` supplies a six-site plane transporter
and proves that the infinite cylinder of physical period `(2,2)`
has minimum nonzero cycle weight exactly six. Repetition gives
an original-square logical of weight 18r for every integer r>=1,
with an explicit ordinary endpoint charge and nonzero period flux.
It is a cylinder-image class, and its weight does not improve the
known unrestricted 12r upper bound.

More importantly, the entire ansatz is excluded: every nonzero
original-square cycle invariant under translation `(2,2)` has
weight at least 12r for r>=2, strictly above the desired 12r-6.
For every r>=4 its invariant-sector minimum is exactly 18r,
including when restricted to nontrivial cycles. The proof uses
an elementary cylinder lower bound and a bounded-support lifting
argument, not a finite distance search. The parent fixed replay
passes in 0.0001 seconds.

These are **sector minima, not the quantum distance of the whole
code**. An attaining square construction must break that symmetry;
correcting the explicit 18r family to the target changes at least
6r+6 physical sites, so a fixed-size endpoint patch cannot suffice.
General seam networks and growing-period motifs remain open.

### 8.5 Recent primary-source scope check

`IBM_DISTANCE_LITERATURE_CHECK_20260908.md` records a focused
check of recent spectral, cover, divisor-search, and Frobenius
papers. In particular the stabilizer-excluded colon refinement of
[Sabo, Can, Marquis, Theorem 90](https://arxiv.org/html/2608.27565v1)
assumes semisimplicity, whereas the original IBM periods are even.
The unrestricted colon term is capped by the six-site checks.
Neither statement can simply be imported as the missing growing
floor. This motivates a precise block-purity audit, rather than a
new enumeration or a claim that the inspected papers solve the goal.

### 8.6 All-size block-purity theorem and an exact metric obstruction

`A43_S4_block_purity.md` proves a new class-space theorem for the
original b=1 rectangles at every r>=1. The first-only and
second-only physical homology images coincide with the image K
of cycles whose components both lie in their principal stabilizer
ideals. This is exactly the cylinder axis along the period with
the larger **2-adic valuation**, has binary dimension six, and
contains 63 nonzero classes. Both physical components of every
class outside K must escape their respective stabilizer ideals.
The remaining 4,032 nonzero classes include all 3,906 winding-mixed
classes. These assertions do not conflate winding, physical-block,
and omega-fibre purity.

The proof computes the two local ordinary algebras, their
annihilator intersections, and the natural period-Tor maps; it
excludes the order-127 orbit using consecutive r,r+1. It identifies
the semisimple-shortcut obstruction as
`ann(PQ)/(ann(P)+ann(Q))`, of dimension six, rather than assuming
that quotient vanishes. The parent read the proof and replayed
its fixed gross/cover identities in 0.0016 seconds.

The resulting classical quotients both identify with H/K. This
gives a valid **class-matched** lower bound obtained by summing
the minimum weights in the two component cosets of one nonzero
quotient label. However, the additional parent note
`IBM_RELATIVE_COLON_OBSTRUCTION.md` proves that this strengthened
relaxation is **exactly six at gross**, not twelve: both relative
classical minima are exactly three, and one matched label has
three-site representatives in each block. Those independently
chosen blocks have a nonzero joint syndrome and do not form a
logical together.

The lower bound three is an analytic two-term exclusion using
all local components of R/(P) and R/(Q); the upper witnesses
come from two pinned, reverified banked twelve-site cycles. A
bounded pass over the 1,884 stored gross vectors supplies no new
distance census and makes no all-r extrapolation. Thus the
proposed individual values 6r and even the matched sum 12r are
falsified at r=1. A sharp algebraic metric proof must retain joint
syndrome compatibility; principal exclusion and matching six
homology bits alone are insufficient.

The exact-six obstruction passed an independent read-only review
by the S4 agent; its pinned-source diagnostic passed in 0.92 seconds
at 48.1 MB. All three native agents have completed this bounded
pass. No old headless session was resumed, and the new work remains
uncommitted and unpushed. The exact IBM-family distance target is
still open; the comparison-family theorem is not substituted for it.

## 9. Uniform finite-excess structure, not just an additive gap

### 9.1 Canonical stability with the physical coordinates retained

`IBM_CANONICAL_STABILITY.md` proves, for every positive q with
3 not dividing q and every actual finite-x solution of the
original binary canonical syndrome equation,

    k <=17 delta,       |e+e_star| <=55 delta,
    delta=|e|-q,        e_star=(x^-3 N_q,0).

The required sharp coefficient is still 3/2, not 17. Nevertheless
this is a uniform Hamming-stability theorem, not just uniqueness
or an additive gap. The proof marks rows with positive graph
slack or aggregation loss, proves a local endpoint-aware exclusion
of D->C, and controls every remaining clean interval using the
actual Laurent monomial recurrences. It accounts for all transient
edges and the previously problematic recurrent components.

For the fixed canonical mixed stratum this gives the improved
lower bound `max(4q+2,2 ceil(65q/32))`. It is not a new global
IBM-member distance bound. More structurally, each nonzero
connected defect costs at least two units of excess; hence an
excess-two error has exactly one connected defect of at most
110 physical sites.

For dyadic heights, the note proves a complete **bounded-excess**
counterexample reduction. A first support-span proof gave a
height-512 cutoff at excess two. An independently reviewed
canonical-row surgery strengthens it to **q=1,2,4,...,64**, with
the defect in actual columns [-330,324]. It shortens proved
homogeneous D-loop interiors and then pads to a dyadic height,
preserving every actual row equation, excess, and hole count.
For any fixed excess cap Delta there is an analogous explicit
cutoff proportional to Delta. No finite window is asserted for
unrestricted excess, and no enumeration of the large finite gate
has been run or proposed. Its purpose is to constrain an analytic
next-shell classification, not to restart the old distance races.

The fixed graph bookkeeping passes in 0.00016 seconds at 17.9 MB.
Independent review passed both stability constants and the
connected-defect/lifting argument. The sharp next-shell inequality
and the pure omega parity recursion are the two bounded analytic
follow-ups; neither is assumed proved here.

### 9.2 The next-shell hidden-pair branch is eliminated

`A42_S11_next_shell.md` proves that every dyadic original binary
canonical-syndrome solution of weight q+2 has zero aggregation
loss L and total graph slack two. Its proof exhausts the entirely
zero-slack recurrent components when L=2. In the six-state component,
one hidden pair permits only one D->C transition; the actual row
equations restrict it to six local binomials, at most one triangle
excursion and seven D loops. This reduces that branch to two
height-eight words, both excluded algebraically. The D-only and
separate g-only components are excluded independently.

The note also proves that the height-two weight-four shell is
empty, using a two-layer Laurent factorization. The parent audited
the complete arguments and replayed the fixed identities in
0.0004 seconds at 17.6 MB. There is no support or height search.
An independent read-only review also passed the full all-dyadic
hidden-pair elimination, including both remaining height-eight
contradictions and the endpoint qualifications.

The remaining excess-two case has either one slack-two edge or
two slack-one edges, no hidden pairs, and actual monomial locations
still to determine. Zero-slack transient edges can participate in
these charged closed walks and cannot be discarded. The desired
bound k<=3 is **not** proved; in particular the floor for holes
4–6 has not been promoted to a theorem.

### 9.3 Physical parity retains an indispensable joint constraint

`A42_S12_relative_Plotkin.md` derives exact even/odd-row equations
and boundary updates for the pure omega cylinder. Eliminating a
unit gives a lower-height effective Koszul pair together with a
four-term physical cost depending on a shared correction variable.
It does not give a scalar distance recurrence or a weight-preserving
normalization. The note also proves that either relative parity
projection identifies the full homology class via a degree-five
annihilator quotient.

Nevertheless separate parity minimization fails even when matching
the **full** class: at height two its global relaxed value is
exactly four, compared with the banked actual distance six. Two
explicit same-class cycles supply inexpensive opposite halves;
their splice has nonzero syndrome. Transfer gives the same
obstruction at every higher dyadic q, with a matched-parity sum
at most 2q below the desired pure floor 3q. The actual source
cycles weigh 5q and 9q and are not distance counterexamples.

The parent audited the algebra, full-class matching, transfer
hypotheses, and exact cost, then replayed the fixed polynomial
identities in 0.00047 seconds at 19.9 MB. The sharp bound for the
coupled cost is left open; restating it in reduced coordinates
is not counted as a distance proof.

All three native agents completed this bounded pass. The next
sharp local target is the actual realization of the one-slack-two
and two-slack-one cases in §9.2, not a sweep of the complete
finite gate. The main exact-distance goal remains active and
unproved; the new work is still uncommitted and unpushed, and
the old paused headless sessions have not been resumed.

## 10. Charged-path elimination and genuine Cartier descent

### 10.1 The one-charge case is now completely excluded

`A42_S13_single_charge.md` proves that the original canonical
excess-two branch with one slack-two edge and zero aggregation
loss is empty at every dyadic q>=2. The fixed graph has 39 such
edges, of which only 18 admit a zero return. All 18 are treated,
including transient states, both immediate and late D->C returns,
and the wrapped height-four case. Arbitrarily long returns are
excluded by actual Laurent monomial equations, not by truncation.
The height-two input is S11's analytic empty-shell theorem.

The parent replay passed in 0.0002 seconds at 17.5 MB, and an
independent read-only review passed the entire proof, including
its S8/S11 dependencies. This eliminates that whole branch; it
does not alone prove the sharp next-shell hole bound.

### 10.2 Two charges: only eight entry/exit types remain

`A42_S14_double_charge.md` partitions the 92 compatible ordered
pairs of slack-one edges into 46 unordered pairs. Its actual
row arguments exclude 34 of them. The remaining twelve enter
the isolated g-only two-cycle and leave it once. The note also
proves that any actual remaining solution has dyadic q<=16:
the clean g-only interval has at most nine rows, there are at
most two clean triangle circuits, and a D-loop block has at most
seven rows. The resulting pre-dyadic height bound is 29. This
is an all-height analytic statement, not a scan of q<=16.

The parent continuation `IBM_CANONICAL_TERMINAL_CHARGES.md`
excludes four more unordered pairs, namely entries 16 or 19
with exits 9 or 10 in S14's fixed numbering. Their zero return
cannot contain a triangle circuit: its final first-block
monomial would equal both x^-3 and 1. The remaining g-only
word fails its two-step Laurent recurrence, with q=4 and q=8
handled by explicit wrapped equations. Both the S14 author and
an independent second agent passed this additional proof.

Thus the residual types are exactly entries 0,1,4,5 with exits
9 or 10, at q in {4,8,16}, with unrestricted actual Laurent
positions. The parent replayed S14's fixed checker in 0.00056
seconds at 17.3 MiB and the four-pair extension in 0.00011
seconds at 18.9 MB. Independent review also passed S14's complete
34-pair elimination, endpoint-sensitive loop bounds, and height
reduction. No support search was performed. Two bounded
analytic follow-ups split this last set according to whether
the zero return contains D->C; they must still prove their
actual endpoint claims before a complete next-shell theorem
can be asserted.

### 10.3 Cartier maps preserve classes but do not prove doubling

`A42_S15_Cartier_descent.md` proves four genuine semilinear
chain maps from the original pure omega cylinder of height 2n
to height n. They first multiply the two blocks by Q and P,
respectively, then extract x/t parity with coefficient square
roots. The maps respect both cycles and boundaries in the
same omega convention. Local Smith factors u and u^2 and an
explicit height-two basis determine their action at every
dyadic height, without extrapolating a finite calculation.

Every nonzero source class has three or four nonzero outputs;
all outputs belong to the norm-transfer socle. Their exact
ordinary-point support types imply the unconditional pure
cylinder bound

    d_(2^a) >= ceil(3 (4/3)^a).

This is a genuinely growing bound but has exponent
log_2(4/3), not the required linear exponent. It is neither
the pure floor 3q nor an original IBM torus distance theorem.

Two stronger-looking metric steps are now explicitly falsified.
The raw total-output factor three is attained by a compact
boundary of weight 18 whose four outputs total 54. More
importantly, a nontrivial height-two cycle of weight six has
four output classes with exact minima 3,3,3,5. Their sum is
14>12, disproving a universal factor-two bound even after
each output is individually minimized. These are cost
obstructions, not distance counterexamples. Other genuinely
coupled Cartier inequalities remain possible, but none is
assumed here.

The parent audited the local algebra and the physical cost,
and replayed the fixed identities in 0.00107 seconds at
19.3 MB. Independent read-only review passed the all-height
action, type recurrences, base-class minima, and both cost
obstructions. This distinction between valid homology descent
and invalid metric descent is retained in the failure ledger.

The main exact-distance target remains open. All work in this
continuation is local and uncommitted/unpushed; no old headless
process has been resumed and no Lean or assembly agent was
launched.

### 10.4 A different topological shortcut also has a precise ceiling

`IBM_LOGICAL_PACKING_OBSTRUCTION.md` considers disjoint or
fractionally packed opposite-sector logical witnesses. A valid
class-specific packing of raw value L satisfies L<=n/d_opposite.
The physical inversion-and-block-swap isometry equates the two
BB sector distances. Thus an all-class packing cannot certify an
integer floor D when D(D-1)>=n, even after integer rounding.
Both IBM target formulas satisfy that inequality for every r>=2.

Independent review caught an initial rounding imprecision; the
note now distinguishes the raw sqrt(n) ceiling from its integer
rounding and proves the stronger rounding-aware IBM obstruction.
This is not a universal distance upper bound and does not rule
out joint-syndrome, nonlinear-network, or other topological proofs.
It explains why simply scaling disjoint gross logical witnesses
cannot furnish the desired exact growing certificate.

## 11. The entire next shell is closed, with a sharp attaining family

The two analytic follow-ups in §10.2 now eliminate its final
eight pairs completely. `A42_S16_D_excursion.md` rules out every
D->C excursion by the entry-row certificate

    (X F+1)+X(F+1)=1+X !=0.

It covers an immediate charged exit as well as a subsequent
g-only row, and needs no bound on the intervening path lengths.
`A42_S17_no_D_excursion.md` rules out the complementary DC-free
case. Its two endpoint equations force identical delayed g
monomials, and then a nonzero exit monomial would satisfy X E=E.
Its two short height-four words have separate wrapped identities.

Together with S11, S13, S14 and the parent terminal-entry proof,
this closes **all** possible weight-q+2 solutions for every
dyadic q>=2. `IBM_CANONICAL_EXACT_GAP.md` records the precise
new exact theorem:

    min {|e| : D e=N_q, e!=e_star} = q+4,
    e_star=(x^-3 N_q,0).

The lower bound is the complete shell exclusion plus parity and
the earlier unique-leader theorem. The matching all-height
representative is e_star+boundary(x^-3): its boundary has exactly
six sites and meets the leader once, so its weight is q+4.
The special case q=2 has no hidden periodic collision; it and
the prescribed q=4 attainer were directly rechecked. These fixed
checks supplement the explicit all-height support argument.

For the specified canonical mixed-fibre support, this proves
the sharp hiding-mass inequality whenever the hole count k<=6,
extending the earlier k<=3 result. The weight expression
6q+3 delta-2k is at least 6q because any nonleader now costs
delta>=4. The canonical leader attains 6q. Each nonzero connected
defect likewise costs at least four excess units, so its component
count is bounded by delta/4 rather than delta/2.

The parent audited both final branch proofs and replayed their
fixed diagnostics in 0.00021 and 0.00015 seconds, each below
20 MB RSS. Both S16 and S17 passed independent read-only review.
Further reviews passed the complete branch exhaustion, matching
q+4 attainer, boundary-dressing identity, k<=6 corollary, and
component-count improvement. All three native agents completed
this bounded pass. No support, height, or distance search was
used to discharge the universal quantifier.

This is an exact theorem in one affine cylinder problem, **not**
the requested IBM-family exact distance. Arbitrarily many holes,
other omega supports, the sharp pure floor, the opposite cylinder
orientation and the toroidal mixed classes remain outside it;
the square-column attaining family is still separate. The goal
remains active, and this continuation remains uncommitted and
unpushed. No old paused headless session has been resumed.

## 12. Three pair leaders, a transverse obstruction, and a nearby code

`A42_S18_pair_leaders.md` gives a sufficient all-excess route for
the fixed canonical hiding-mass problem. Proving that each of
the three two-row-support unions is a weight-2q coset leader
would imply the desired inequality by summing their three
boundary inequalities. Equivalently these are the stronger
individual hole budgets 2k_A,2k_B,2k_C<=delta. All three hold
unrestricted at q=1,2, for arbitrary row masks in the specified
three-column generator strip, and for single-column generators.
The unrestricted assertion at q>=4 remains open. An explicit
x^3-folded family violates the proposed folded floor for dyadic
q>=8; its nonzero unfolded syndrome prevents it from being a
counterexample to the actual cylinder problem. Independent
review passed the reductions and all stated restricted results.

`A42_S19_transverse_quotient.md` closes another tempting shortcut
negatively. Folding the pure-omega cylinder by x^6+1 is indeed
class-injective, with image the first axis in the period-Tor
coordinates. Nevertheless, for every dyadic q>=4 the explicit
torus vector (u^(q-2) P_omega,0) lies in a nonzero class on that
very axis and has physical weight 3q/2. Its unfolded syndrome is
(x^6+1)u^(q-2), not zero. Thus keeping the exact cylinder-image
classes does not repair the factor-two metric loss of this
quotient. It neither refutes the actual pure floor nor supplies
a cheap cylinder representative. Independent review passed.

`IBM_ONE_CHECK_MODIFICATION.md` studies a genuinely nearby code:
keep P=x^3+y+y^2 and replace only the other physical check by
Q_new=(1+x+x^2)(y+x+1). The ideal is exactly the original ordinary
primary ideal, including its doubled tangent component. On
all multiple-of-six rectangles it has n=2 ell m and k=12, with
check weight eight instead of six. A nontrivial one-row stripe
has weight 2 ell/3, so the original numerical conjecture cannot
hold for this modification. Its exact growing distance remains
an admissible but unproved alternative target.

The parent analytic result `IBM_ONE_CHECK_CYLINDER_REDUCTION.md`
removes the entire toroidal branch for modified squares. Any
cycle of weight below 2m/3 has every occupied f row divisible
by 1+x+x^2 and hence containing at least two sites. Counting
actual crossed row seams then gives a weight-preserving lift
to the x-periodic, y-unbounded cylinder. The projection is
class-injective, and the stripe is already a nontrivial cylinder
cycle. Consequently, for every multiple-of-six ell<=m,

    d_modified_tor(ell,m)=d_modified_cyl(ell).

In particular d_modified(6r,6r)=d_modified_cyl(6r). A single
remaining lemma, d_modified_cyl(6r)>=4r, would match the stripe
and give the exact [[72r^2,12,4r]] family. At present only
2r<=d_modified_cyl(6r)<=4r is proved. This cylinder is not the
original A42 pure-omega cylinder. Independent review passed the
ideal, homology, stripe, row-divisibility, seam, and exact
distance-reduction proofs; no metric normal form is assumed.

The parent replayed all three new fixed diagnostics; each ran
below 0.001 seconds and 20 MB RSS. No distance, support, coset,
or height search was used. Neither the original conjecture nor
the modified family's numerical exact distance is proved.
These results remain local, uncommitted, and unpushed.

## 13. Exact one-check modified family: both cuts close

### 13.1 The main theorem is now unconditional for the modified code

`IBM_ONE_CHECK_CYLINDER_TANGENT.md` closes the single cylinder
floor left in §12. The fully self-contained proof is consolidated
in `IBM_MODIFIED_EXACT_DISTANCE_THEOREM.md`. Its exact statement is

    P=x^3+y+y^2,       Q=(1+x+x^2)(y+x+1),
    [[2ell m,12,2ell/3]],       6|ell,m and ell<=m.

In particular, every integer r>=1 gives the exact family

    [[72r^2,12,4r]] on ell=m=6r.

The result retains the original first check, two qubits per cell,
the square frames, twelve logical qubits, and the whole ordinary
primary module including its doubled component. Only the second
check changes; its weight becomes five and total check weight is
eight rather than six. The exceptional order-127 component is
removed. The distance is genuinely exact and grows with r, but
its coefficient is smaller than IBM's conjectured value.

The proof has no residual metric lemma. A torus cycle below the
claimed floor first lifts weight-preservingly to the x-cylinder
by the row recurrence and occupied-row seam count from §12.
On that actual cylinder, apply the physical shear

    f'=X f(X,XY),       g'=g(X,XY).

The new checks are P'=X^2+Y+XY^2 and
Q'=1+X^3+(X+X^2+X^3)Y. Their Y augmentations are both chi(X).
The contracting invariant c=f'(X,1)+g'(X,1) must vanish below
2ell/3. Therefore each occupied X column contains at least two
qubits, across both blocks. If K is the number of occupied
columns, the union of all possible crossed column seams is
bounded by 3K<=3w/2<ell. An empty column seam lifts the entire
cycle to a finite exact plane boundary. Projecting back proves
the all-class, both-sector lower bound.

Crucially, the contracting invariant need not detect every class.
Its kernel contains the full four-dimensional doubled component;
the second geometric cut controls that kernel directly. There is
no unproved metric normalization, transverse quotient isometry,
or assumption that class-minimal representatives are connected.

The matching physical stripe is (0,W), where
W=(1+x)sum x^(3j) has 2ell/3 sites and chi W=x^ell+1.
Its nontriviality is proved by its finite period syndrome V F:
a hypothetical torus boundary would force V=0 in the plane
quotient, whereas evaluation at x=y=omega sends V to one.
The quotient's dimension six and rank symmetry give k=12,
including all repeated-root and exceptional-prime periods.

### 13.2 Failed normal forms are retained, not silently assumed

The parallel notes `IBM_ONE_CHECK_CYLINDER_NORMAL_FORM.md` and
`IBM_ONE_CHECK_CYLINDER_GEOMETRY.md` independently identify the
same 48 detected classes and fifteen-class kernel. They prove
useful restricted bounds, but also exhibit why standard algebraic
shortcuts were insufficient: the fixed normal form can double
the weight of a true leader, and a prescribed highest-row peel
can take an actual boundary from weight 38 to the threshold 40.
No monomial slope realizes the doubled tangent, and substituting
y=x+1 can amplify physical weight without bound. The final two-cut
proof uses none of those false metric assertions. These failures
remain valid route obstructions even though the target is closed.

### 13.3 Literal one-monomial and general contraction alternatives

Two additional exact results are now proved, with their differing
scope kept explicit:

- `IBM_ONE_MONOMIAL_EXACT_FAMILY.md` keeps original A and changes
  only x^-1 y^-3 to x^-1 in B. It gives [[72r^2,8,4r]] on the
  original square frames for every r>=1. The logical dimension
  is eight, and the module is product-like rather than tangent.
- `IBM_ONE_MONOMIAL_K12_FAMILY.md` keeps original A and changes
  only x^-1 y^-3 to x^5 in B. On new squares of side 7r it gives
  [[98r^2,12,4r]] when 3 does not divide r, and
  [[98r^2,20,10r/3]] when 3 divides r. Setting r=3t+1 gives a
  constant-k=12 infinite family for every t>=0. Its reduced
  plane algebra is different from original IBM's.

Both use two individually contracting monomial evaluations into
an exact classical cyclic code. Their joint cylinder-class
detection is proved by an algebraic normal form, but no physical
weight is assigned to that normal form. A separate low-weight
seam lift handles the torus classes the evaluations miss.

`IBM_TWO_CONTRACTIONS_EXACT_THEOREM.md` establishes the general
analytic mechanism. For a specified primitive binary trinomial
q=1+x^s+x^t of degree t>=2, take
P=x^(s+t)+y+y^2, Q=q, and square side (2^t-1)r. The exact family is

    [[2(2^t-1)^2 r^2, 4t, 2^(t-1)r]] for every r>=1.

Primitiveness is a stated polynomial hypothesis, not an assumed
existence theorem in every degree. The explicit choice s=1,t=3
gives [[98r^2,12,4r]] for all r, by changing two original monomials.
These examples are not claimed novel or distance-equivalent to
the original IBM checks. They are additional exact constructions,
not substitutes for the closer one-check theorem in §13.1.

### 13.4 Verification and completion scope

The parent audited every new proof and replayed all prescribed
new diagnostics, each below 0.001 seconds and 20 MB RSS. The
flat eight-logical-qubit theorem passed two independent reviews.
The main twelve-logical-qubit cylinder proof and its composition
passed two independent reviews, and a separate review passed the
full self-contained theorem, including both cuts, dimension,
and the repeated-root nonboundary argument. No support, coset,
height, or distance search supplied a universal quantifier.

This completes the requested exact growing-distance objective
under its explicitly permitted closely related modification.
It does not prove the original formula 6(2r+b-1), remove its
remaining mixed-sector obligations, or transfer an exact metric
back to the original check pair. No Lean formalization is claimed.
These proofs, diagnostic sources, and scope corrections are included
in the PR #27 research publication update. Earlier local/unpushed
statements record historical checkpoints before this update.
Protected user files and old paused headless sessions have not
been touched.

All three native agents completed the final proof-and-review pass.
