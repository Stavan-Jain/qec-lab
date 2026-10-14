# A42 S9 — exact tail interfaces, sharp clipping obstructions, and safe peeling

2026-09-08. Native Codex sidecar; original Laurent pair only. New note
and new diagnostic only. No growing-width search, periodic x folding,
or linear recovery construction is used.

## 1. Results and limitations

**Unrestricted canonical HM remains open.** This pass proves more
than the general warning that clipping can destroy cancellations:

- At every even height, one explicit original-cylinder generator
  has energy `11q`, but clipping its single left-tail column gives
  `17q`. Moving that column into **any** of the six strip positions
  also strictly raises energy. Translating the entire generator into
  the strip raises energy to `15q`. These are exact counterexamples
  to universal monotonicity of those operations, not to HM.
- The worst increase on deleting the leftmost generator column a is
  at most `6|a intersect t a|`. On deleting the rightmost column a
  it is at most `3|(t+t^2)a|`. Both bounds are attained by explicit
  examples, at every even height.
- Sharper, context-sensitive deletion inequalities give actual safe
  local rules, using only the extreme mask and its three neighbors.
  In particular, left columns with no cyclically adjacent occupied
  rows and right columns constant in the row direction are always
  safely deletable, regardless of the rest of the generator.
- Consequently the six-column theorem extends to generators with
  **arbitrarily long tails of these types**. No restriction on tail
  length or a growing transfer width is involved.
- A correlated replacement repairs the displayed obstruction while
  preserving its entire inside boundary trace. Thus the obstruction
  does not disprove a more flexible strip-completeness theorem.

The remaining work is coverage of stalled, possibly negative-energy
configurations by safe correlated replacements, or a way to account
for their positive deletion charges. That coverage is not proved.

## 2. Conventions and why the extreme-column calculation is local

Let `q>=1`, normally `q=2^a`, and put

```
R=F_2[t]/(t^q+1),       u=t+t^2,       v=t^3,
P=x^3+u,               Q=1+(x+x^2)v,
partial h=(Qh,Ph).
```

The coordinate x is an **unbounded Laurent coordinate**. Every
generator h has finite x support. A row mask a is a binary element
of R; `|a|` is its Hamming weight in the t basis. Intersections in
this note are intersections of binary supports, not polynomial
multiplication. Addition of masks is XOR. Multiplication by t or v
is a permutation, while multiplication by u is not.

Write `h=sum_j x^j h_j`, using zero for missing columns. Then

```
(Qh)_j=h_j+v(h_(j-1)+h_(j-2)),
(Ph)_j=h_(j-3)+u h_j.                                  (2.1)
```

Fix the canonical S: f-column zero and g-columns zero and one,
each at every row. For `s=partial h`, define

```
I(h)=|s intersect S|,
O(h)=|s outside S|,
E(h)=3O(h)-I(h)=3|s|-4I(h).                             (2.2)
```

The canonical mixed physical weight is `6q+E(h)`. Full canonical
HM is `E(h)>=0` for all compact h. The S8 finite certificate proves
this for all h supported in `[-3,2]`, for every height.

If j is the leftmost generator column and `j<=-4`, deleting it
changes f only in `j,j+1,j+2` and g only in `j,j+3`. All of these
positions are outside S. If j is the rightmost column and `j>=3`,
the same conclusion holds. In either case the **inside boundary
trace remains exactly unchanged** under deletion.

For any binary masks A and B,

```
|A+B|=|A|+|B|-2|A intersect B|.                         (2.3)
```

All the following deletion identities follow from (2.1) and (2.3).
They do not assume that a tail is independently supported as a
boundary, nor that a class-minimal representative is connected.

## 3. Leftmost-column deletion: exact identity and sharp bound

**Theorem 3.1.** Let j be the leftmost occupied generator column,
with `j<=-4`. Set

```
a=h_j,       b=h_(j+1),       c=h_(j+2),       d=h_(j+3),
h'=h+x^j a.
```

Thus h' deletes a, without changing the other generator columns.
Put

```
L(a,b,c,d)
  = |v a intersect b|
    + |v a intersect (c+v b)|
    + |a intersect u d|.
```

Then the exact energy change is

```
E(h')-E(h) = 3(2L(a,b,c,d)-4|a|-|u a|).                (3.1)
```

In particular,

```
E(h')-E(h) <= 6|a intersect u d|-3|u a|                 (3.2)
             <= 6|a intersect t a|.                    (3.3)
```

Proof. The boundary of the deleted column has masks a,va,va in
f-columns j,j+1,j+2, and masks ua,a in g-columns j,j+3. Its weight
is `4|a|+|ua|`. The boundary of h' has zero f- and g-masks at j,
and at the other three potentially overlapping positions has

```
(Qh')_(j+1)=b,
(Qh')_(j+2)=c+v b,
(Ph')_(j+3)=u d.
```

Thus its overlap with the deleted boundary is exactly L. All
changed boundary sites are outside S, proving (3.1) by (2.3).
The first two overlaps are at most `2|a|`, proving (3.2).
Finally

```
|u a|=|a+t a|=2|a|-2|a intersect t a|,
```

and `|a intersect u d|<=|a|` proves (3.3). QED.

**Safe local rule.** It is sufficient to check

```
2|a intersect u d| <= |u a|                             (3.4)
```

before deleting a. The exact guard `2L<=4|a|+|ua|` is less
restrictive and is necessary and sufficient for this particular
deletion to be nonincreasing. Both are explicit four-mask tests,
not a test of the unknown global HM inequality.

Useful unconditional cases of (3.4) are:

- `u d=0`, i.e. the third neighboring column is row-constant;
- `a intersect u d` is empty;
- `a intersect t a` is empty, i.e. a has no cyclically adjacent
  occupied rows. In the last case `|ua|=2|a|` and (3.4) follows
  regardless of d.

In particular, at any `q>=2`, a leftmost column of weight one is
always safely deletable. The claim concerns a whole extreme
column: deleting one bit from a multi-bit extreme column does not
have the same exclusive g-mask ua and is not justified by it.

The term `u h_(j+3)` is essential. It is the precise third-neighbor
interface that a two-column clipping argument would miss. Keeping
its actual row mask, not only its support size, gives the sharper
condition (3.4).

## 4. Rightmost-column deletion: a different exact rule

**Theorem 4.1.** Let j be the rightmost occupied column, `j>=3`, and
write

```
a=h_j,       b=h_(j-1),       c=h_(j-2),       d=h_(j-3),
h'=h+x^j a.
```

Set

```
R(a,b,c,d)
  = |a intersect v(b+c)|
    + |a intersect b|
    + |u a intersect d|.
```

Then

```
E(h')-E(h) = 3(2R(a,b,c,d)-4|a|-|u a|)                 (4.1)
             <= 6|u a intersect d|-3|u a|              (4.2)
             <= 3|u a|.                              (4.3)
```

Proof. This time the overlapping f-masks of h' are `v(b+c)` at
j and `v b` at j+1. Its g-mask at j is d. The deleted boundary's
other two positions, f at j+2 and g at j+3, cannot overlap h'.
The total deleted weight is again `4|a|+|ua|`. Formula (2.3) gives
(4.1), and bounding the two f-overlaps by `2|a|` gives the rest.
QED.

**Safe local rule.** It suffices that

```
2|u a intersect d| <= |u a|.                            (4.4)
```

The exact guard is `2R<=4|a|+|ua|`. In particular, if `u a=0`,
then deleting a is always safe, regardless of the remaining
generator. Such a is exactly zero or the all-one mask N. Also,
`d=0` is an unconditional sufficient context for any a.

This is deliberately **not** obtained by declaring the pair
self-mirror: the left rule contains `a intersect u d`, whereas the
right rule contains `u a intersect d`. Their worst-case bounds
also differ. The calculation uses the original P,Q on both sides.

## 5. Explicit all-even-height obstructions and a genuine repair

Let q be any positive even integer, and put

```
N=sum_(r=0)^(q-1) t^r,
J=sum_(r=0)^(q/2-1) t^(2r),
K=N+J=tJ.
```

Then `tJ=K`, `tK=J`, `uJ=uK=N`, `uN=0`, and `vJ=K`.
The masks J and K each have weight q/2.

### 5.1 One left-tail column defeats clipping and every direct move

Take the actual Laurent generator

```
h=N(x^-4+x^-3)+J(x^-1+x).                               (5.1)
```

Its boundary is exactly

```
Qh=N x^-4 + K(x^-1+1) + N x + K(x^2+x^3),
Ph=N(1+x)+J(x^2+x^4).
```

Hence `I=5q/2`, `O=9q/2`, and `E(h)=11q`. There is just one
generator column outside `[-3,2]`: the column -4 with mask N.

| operation on h | resulting energy |
|---|---:|
| leave h unchanged | `11q` |
| delete column -4 | `17q` |
| move its mask N to column -3 | `12q` |
| move its mask N to column -2 | `15q` |
| move its mask N to column -1 | `14q` |
| move its mask N to column 0 | `17q` |
| move its mask N to column 1 | `14q` |
| move its mask N to column 2 | `23q` |
| translate the entire h one column right | `15q` |

A move to column d means `h -> h+x^-4N+x^dN`, with XOR at an
already occupied destination. These are all possible single-column
moves into the assigned strip. Row translations of the moved N
do not change it, so allowing a simultaneous t shift does not
avoid the obstruction.

The original extreme generator columns are -4 and 1. Therefore
the **only** horizontal translation placing the entire generator
inside `[-3,2]` is multiplication by x. Any additional global t
translation leaves all the listed weights unchanged. Thus merely
translating an already narrow generator to fit the canonical strip
is not a universally nonincreasing operation either.

Deleting -4 attains the sharp upper bound in (3.3): here
`a=N`, `b=N`, `c=0`, `d=J`, all three overlaps in (3.1) have weight
q, and the energy increase is `6q=6|N intersect tN|`.

**Verification for every even q.** All displayed boundary masks
and all masks in the table are in `{0,J,K,N}`. The map from
two-row masks to q-row masks that repeats a two-row word q/2
times is injective, linear, commutes with t, and multiplies
Hamming weight by q/2. Applying it coefficientwise in the actual
Laurent x coordinate commutes with P,Q and all the specified
generator moves. Consequently the finite q=2 XOR table certifies
the displayed formulas at every even q. This is a proved
intertwining identity, not extrapolation from selected heights.

### 5.2 Correlated replacements do repair this obstruction

The two-column replacement

```
h_two=h+x^-4N+x^-3N+x^-1N = K x^-1+J x
```

lies in the strip and has `I=3q/2`, `O=7q/2`, hence energy `9q`.
This operation changes the inside trace, but reduces its actual
weighted energy by `2q`; it is a valid repair for (5.1).

There is also an **inside-trace-preserving** repair:

```
h_trace=N(x^-3+x^-2)+K x^-1.                            (5.2)
```

Both (5.1) and (5.2) have inside masks

```
(Qh)_0=K,       (Ph)_0=N,       (Ph)_1=N.
```

For (5.2), the outside weight is `7q/2`, so its energy is `8q`.
Thus it strictly improves the energy by `3q` while preserving all
inside bits. Direct substitution, or the same repeating-word
identity, proves these assertions for every even q.

These repairs explain the next action determined by the obstruction:
a complete elimination argument must permit correlated edits of
several retained columns, rather than choosing only clipping or a
single monomial move of the old tail. No universal validity of the
two displayed repair formulas for arbitrary h is claimed.

### 5.3 The right-hand bound is sharp too

Take

```
h_right=N(1+x)+J(x^2+x^3).
```

Its boundary is

```
Qh_right=N+J x^2+K x^5,
Ph_right=N(x^2+x^4)+J(x^5+x^6).
```

It has `I=q`, `O=4q`, and energy `11q`. Deleting its single
right-tail column `x^3J` leaves inside weight q and raises outside
weight to `5q`, hence raises energy to `14q`.
Here `a=J`, `b=J`, `c=N`, `d=N`. All three overlaps in (4.1)
saturate their possible values, and the increase is exactly
`3q=3|uJ|`, attaining (4.3).

All energies in this section are positive. None is a violation of
original HM or of the IBM distance conjecture. In particular these
examples do not refute an elimination theorem restricted to
negative-energy inputs with more flexible replacement moves.

## 6. An arbitrary-tail-length theorem and a quantitative penalty

Let c be the generator obtained by retaining only columns `[-3,2]`.
Delete left-tail columns in increasing order of their index, and
right-tail columns in decreasing order. At each nonzero deletion,
the removed column is genuinely extreme. Since all other generator
masks are unchanged, (3.3) and (4.3) telescope to

```
E(c)-E(h)
  <= 6 sum_(j<-3) |h_j intersect t h_j|
     +3 sum_(j>2) |u h_j|.                             (6.1)
```

This is a uniform, finite-support inequality. It bounds the
possible cost of clipping in terms of actual row structure in
the discarded columns, with the sharp one-step constants from
§§3–5. It is not asserted to be the optimal many-column bound.
The stronger bounds (3.2) and (4.2), evaluated at each stage, may
replace their respective terms to retain neighbor correlations.

**Theorem 6.1 (safe arbitrary-length tails).** For every positive
height q, suppose h is any finite Laurent generator satisfying

```
h_j intersect t h_j is empty for every j<-3,
u h_j=0 for every j>2.
```

There is no condition on the six retained columns. Then

```
E(h) >= E(c) >= 0.
```

Proof. Every term on the right of (6.1) vanishes. Equivalently,
each left deletion is safe by Theorem 3.1 and each right deletion
is safe by Theorem 4.1. There are finitely many columns, so the
procedure terminates at c without increasing energy. The final
inequality is exactly the all-height S8 strip certificate. QED.

More generally, the same conclusion holds whenever the four-mask
guard in §3 or §4 succeeds at each corresponding deletion. The
guarded algorithm must **stop**, without claiming success, if a
remaining outside extreme column fails its guard. The left
example (5.1) is an explicit such stalled input, so this caveat
is substantive rather than a hypothetical implementation issue.

This theorem adds truly unbounded x spans to the known strip result.
It does not restrict a transfer computation to a larger strip,
and it does not assert that arbitrary masks can be made to satisfy
the hypotheses at no energy cost.

## 7. Exact two-sided interface bookkeeping and the remaining task

For completeness, write `h=ell+c+r`, where ell is supported at
`j<=-4`, c at `[-3,2]`, and r at `j>=3`. The boundaries of ell
and r have disjoint support, both entirely outside S. Their only
possible overlaps with the boundary of c are respectively

```
left:  f-columns -3,-2;     g-columns -3,-2,-1;
right: f-columns 3,4;       g-columns 3,4,5.
```

Define

```
K_left = |partial ell|-2|supp(partial ell) intersect supp(partial c)|,
K_right= |partial r|  -2|supp(partial r)   intersect supp(partial c)|.
```

Then the exact identity is

```
E(h)=E(c)+3K_left+3K_right.                              (7.1)
```

The actual left boundary trace on these five positions depends on
the last three tail masks `a=ell_-4,b=ell_-5,d=ell_-6` and is

```
(v(a+b), va; d,b,a).
```

The right trace depends on `a=r_3,b=r_4,d=r_5` and is

```
(a,b+va; ua,ub,ud).
```

These identities follow immediately from (2.1), and (7.1) from
(2.3). In particular the discarded tail's positive standalone
boundary weight is not its conditional cost after coupling to c:
the overlap subtraction is essential. Example (5.1) has
`K_left=-2q`, despite its tail itself having positive weight.

The S8 endpoint warning has an exact counterpart here: removing
the third-neighbor derivative contribution breaks the calculation.
No periodic x-folding can restore it; the prior q=512 folded
counterexample remains a separate prohibition on that shortcut.

**Precise residual.** For arbitrary h, neither K in (7.1) is
proved nonnegative, and such a claim is false. To close full HM
by tail elimination, one still needs a proved supply of correlated
nonincreasing replacements that remove an outside extreme column
when the safe deletion rules fail, or a sharp inequality paying
the negative conditional costs from retained-strip energy.
In a replacement proof, the remaining outside span must decrease
and no new outer column may be introduced; these give a real
termination measure. Neither the needed coverage of stalled
negative-energy inputs nor the compensating sharp inequality is
established here. They are not being renamed as a proved lemma.

This is distinct from the parent's column-load-at-most-three linear
recovery route. No reconstruction matrix for inside boundary bits
from outside bits has been built or assumed in this sidecar.

## 8. Diagnostic and scope

`a42_s9_codex_tail_certificates.py` uses sparse actual x columns
and exact integer row masks. It verifies the stated q=2 boundary
counts, both sharp deletion charges, every entry in the prescribed
six-destination move table, the unique whole-generator translation,
and both repairs, including equality of all inside trace bits for
(5.2). Every expected count is asserted before a PASS result is
emitted. It performs no optimization or search over strip widths.

The all-even-height extension is the proved row-repetition map
in §5, not a diagnostic extrapolation. The universal deletion
identities and Theorem 6.1 are analytic arguments, not inferred
from the test configuration.

The final diagnostic body took about 0.0001 seconds and used less
than 18 MiB measured peak RSS. Invocations were far below the
55-second wall / 50-second CPU caps and 768 MiB memory limit;
aggregate computation was below one second. No old process,
git state, shared note, or sibling script was modified.

## 9. Brief independent review of the parent's recovery obstruction

After the primary tail work, I read
`IBM_CANONICAL_RECOVERY_OBSTRUCTION.md`. I found no gap in its
width–load argument or the explicit inverse:

- A three-empty-row gap puts the unwrapped orthogonality polynomial
  in `[-3,q-g-1]`, containing at most q integers. Reduction modulo
  q is injective on precisely this interval; there is no endpoint
  off-by-one loss in the plane-lifting argument.
- Coprimality in the Laurent UFD suffices for the plane cochain
  to be a syndrome-row combination; comaximality is not needed.
  Reduction to row period three preserves the single inside bit,
  while a suitable translate of `(1+t,(1+x)(1+t))` is a cycle
  pairing nontrivially with that bit. Its outside support is empty.
  This gives the required contradiction without transferring a
  distance claim between periods.
- Each decoder row has at least `ceil(q/3)-1` outside entries.
  There are 3q such rows and at most Kq input coordinates, so
  counting entries gives exactly `Kc>=3ceil(q/3)-3`. No shift
  equivariance of the decoder was assumed.
- Substitution in the displayed elimination formula for j leaves
  `(1+u)j`, with the i and k terms cancelling. The formulas for
  k and i then follow immediately. The ten input block-columns
  and the inverse-mask weights `2q+1` or `2q-1` are consistent.

This review used no new computation. The theorem restricts a
global decoder's input block-columns, not the generator columns
in Theorem 6.1 or in the S8 strip certificate. It therefore does
not obstruct the nonlinear conditional peeling rules proved here.
