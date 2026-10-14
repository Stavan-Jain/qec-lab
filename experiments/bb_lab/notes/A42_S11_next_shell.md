# A42 S11 — excess two has no aggregation-hidden pair

2026-09-08. Bounded native Codex pass on the original cylinder.
Only this note and `a42_s11_codex_next_shell.py` are new.

## 1. Outcome: a rigorous branch elimination, not the full shell theorem

The requested statement

```
P f+Q g=N_q,       |f|+|g|=q+2
    implies       |supp(f,g) intersect S_*|<=3
```

is **not proved in this pass**. No counterexample to it was found
or asserted. The new results are:

1. **No hidden pair at excess two.** At every binary height
   `q=2^a`, an actual compact solution of weight `q+2` has zero
   weight loss under x-mod-three parity aggregation. Consequently
   its entire excess is graph slack two. This eliminates the other
   possible branch, `L=2, sigma=0`, including the six-state
   component's five-cycles and the separate two-cycle.
2. **The height-two shell is empty.** There is no compact actual
   solution of syndrome `1+t` and weight four. The proof is a
   two-layer Laurent identity, not a search over supports.
3. The remaining problem is exactly the actual spatial realization
   of closed walks with either one slack-two edge or two
   slack-one edges, no aggregation loss, and arbitrarily placed
   Laurent monomials. Zero-slack transient edges must be retained
   in this problem.

At q=1 the desired hole bound is immediate from total weight three.
At q=2 it is vacuous by result 2. The hole bound at general binary
height, and therefore its advertised consequence for holes 4–6,
remain unbanked.

The authoritative S7 uniqueness theorem and S8's monomial-word
lemma are used with their actual spatial hypotheses. No uniqueness
statement is extrapolated to the next shell. The parent's separate
coarse stability theorem is not needed for the proofs here and was
not used to launch a support search.

## 2. Definitions and the exact excess split

Let `q=2^a`, `a>=0`, and work in

```
R_q[x^±],        R_q=F_2[t]/(t^q+1),
P=x^3+t+t^2,    Q=1+(x+x^2)t^3,
N_q=1+t+...+t^(q-1),
S_*={all f sites at x=0; all g sites at x=0,1}.
```

All x supports are finite, with x an unbounded Laurent coordinate.
For an actual solution e=(f,g) of syndrome N_q, write its row
polynomials as f_j(x),g_j(x). The actual equation at every cyclic
row is

```
x^3 f_j+f_(j-1)+f_(j-2)+g_j+(x+x^2)g_(j-3)=1.           (2.1)
```

Aggregate x exponents modulo three by parity separately in each
block and row. The physical weight lost is L. Each aggregated
cell with m occupied physical sites loses `m-(m mod 2)`, so L is
nonnegative and even. With the S7 potential on its 16-state graph,

```
delta=|e|-q=L+sigma,
sigma=sum of nonnegative edge slacks.                   (2.2)
```

Thus delta=2 has just two possibilities:

```
(L,sigma)=(2,0) or (0,2).                               (2.3)
```

When L=2, exactly one block/row/residue cell has an extra pair.
It has either two sites above aggregate parity zero, or three
sites above aggregate parity one. Choose one base site in the
latter case and regard the other two as a pair. In either case
the exception is one homogeneous-residue Laurent binomial p at
one row and in one block. Every other occupied aggregate cell
has exactly one physical site. This description does not discard
the possible three-site cell.

Put `X=x^3`, still an unbounded variable, and let Pi_0 retain
exactly the terms of residue zero, retaining their actual powers
of X. The S7 zero-slack recurrent components are:

- the component with states A,B,C,D,E,F and internal edges
  `AB, BC, BE, BF, CA, DC, DD, ED, FD`;
- the separate two-cycle with empty first block and one g monomial
  per row, alternating between residues one and two.

The other six zero-slack edges are transient. A completely
zero-slack **closed** walk cannot use them. This exclusion will
not apply once sigma=2.

## 3. Two zero-slack components that cannot hide one pair

### 3.1 The D self-loop alone

**Lemma 3.1.** An actual solution with L=2 and sigma=0 cannot have
an aggregated walk consisting solely of D->D.

Proof. If the exceptional pair is in f, then g is physically zero.
Multiplication by P is injective on compact Laurent polynomials:
the highest nonzero x column of an input survives three columns
later with coefficient one. Since `P(x^-3 N_q)=N_q`, this forces
`f=x^-3 N_q`, which has no excess.

If instead the pair is in g, every f monomial has residue zero,
and g is one nonzero binomial p of a single residue at a single
row. The three terms of Qg lie in three different x residues;
their respective residue components are nonzero, since they are
unit multiples of p. But Pf and N_q have only residue zero.
The other two components cannot cancel. QED.

### 3.2 The separate two-cycle

**Lemma 3.2.** An actual solution with L=2 and sigma=0 cannot
realize the separate g-only two-cycle.

Proof. Its length q is even. If the pair is in g, f is physically
zero. The minimum and maximum x columns of Qg are distinct,
separated by at least two, with nonzero unit-multiplied extreme
coefficients. They cannot both be the single column of N_q.

Otherwise f is one binomial p of residue rho at one row, while
every g_j is a monomial. If the residue of g_j is one, that of
g_(j-3) is two, and vice versa.

First suppose rho=0. The nonzero-residue part of (2.1) forces

```
g_j=x^2 g_(j-3) when g_j has residue one,
g_j=x   g_(j-3) when g_j has residue two.
```

Summing the differences of the integer x exponents over all q
rows gives `0=3q/2`, a contradiction. This uses actual exponents,
not their residues.

If rho is one or two, the residue-zero part of (2.1) is unaffected
by Pf. It fixes every g monomial:

```
g_j=x^-2 at residue-one rows,
g_j=x^-1 at residue-two rows.
```

With these values, Qg has residue-zero part N_q, residue-two part
zero, and residue-one part

```
(x^-2+x) J,
```

where J is the alternating mask of the residue-one rows. Thus
rho=2 is impossible immediately. For rho=1, Pf must equal this
last expression. If q>=4, Pf has nonzero values at the three
distinct rows r,r+1,r+2 of the exceptional f pair: respectively
Xp,p,p. These three rows cannot all belong to one alternating
parity class. If q=2, the two values are `(X+1)p` and p, both
nonzero in the Laurent domain, whereas J occupies just one row.
Both cases contradict the required equation. QED.

## 4. One hidden pair allows at most one D->C transition

**Lemma 4.1.** If an actual L=2, sigma=0 solution uses the
six-state component and is not the D loop, its closed walk has
exactly one D->C transition.

Proof. Without D->C, the only recurrent pieces of this component
are the D self-loop and the directed three-cycle A-B-C-A.
The former is excluded by Lemma 3.1 and the latter cannot have
binary length q. Thus a D->C is present.

For a D->C at row j, the S7 spatial contradiction uses the five
rows `j-1,...,j+3`. If those rows contain no exceptional pair,
the equations force the same monomial to equal both X^-1 and 1.
Thus the unique exceptional row must belong to that window.

Consecutive D->C edges in a zero-slack walk are at least five
rows apart: from C the shortest return to D is C-A-B-E-D or
C-A-B-F-D, with four edges. Consequently these five-row windows
are disjoint, including across the cyclic seam. One exceptional
row cannot meet two of them. QED.

Number the single D->C edge as row zero. The whole walk is now
described by two nonnegative integers n,r and one choice E/F:

```
row 0: D->C;    row 1: C->A;
n copies of A->B->C->A;
A->B->E->D or A->B->F->D;
r copies of D->D;
q=5+3n+r.                                               (4.1)
```

This is a graph consequence with arbitrary loop counts, not an
enumeration of physical supports or an assumption of connectivity.

### 4.2 The exceptional binomial has only six possible local forms

Let F_2 be the base residue-zero f monomial on A->B at row two.
Temporarily remove the exceptional pair. The two equations used
by S7, with the pair's corrections restored, have the form

```
X F_2+A=1,       F_2+B=1.                               (4.2)
```

A is the full correction in row two, and B the residue-zero
correction in row three. Hence necessarily

```
A has residue zero,       A+X B=1+X.                    (4.3)
```

The pair must lie in rows -1 through 3 by Lemma 4.1. Directly
substituting its contributions from (2.1) gives precisely the
following possibilities:

| block | row | exceptional pair p | F_2 |
|---|---:|---|---|
| f | 0 | `1+X` | `1` |
| f | 3 | `X^-2+X^-1` | `X^-1` |
| g | 0 | `x^-5+x^-2` | `X^-1` |
| g | 0 | `x^-4+x^-1` | `X^-1` |
| g | 2 | `1+X` | `1` |
| g | 3 | `X^-1+1` | `X^-1` |

Here is the exhaustive algebra behind the table. For an f pair
at rows -1,0,1,2,3, the corrections `(A,B)` are respectively

```
(0,0), (p,0), (p,Pi_0 p), (Xp,Pi_0 p), (0,X Pi_0 p).
```

At row one, residue zero is forced and (4.3) would require p=1,
not a pair. At row two the two terms in (4.3) cancel. Rows zero
and three give the first two table entries; row -1 gives none.

For a g pair at those five rows, the corrections are

```
((x+x^2)p,0), (0,Pi_0((x+x^2)p)), (0,0), (p,0), (0,Pi_0 p).
```

The row -1 case has two nonzero residues in A and cannot satisfy
(4.3). Row zero has the two residue choices one and two, giving
the two displayed g pairs. Rows two and three give the remaining
entries. Row one gives none.

All six permitted pairs lie in cells whose aggregate parity is
zero. Thus even the initially allowed three-site occupied cell
has now been excluded, rather than silently ignored.

## 5. The two arbitrary loop counts become bounded

### 5.1 At most one A-B-C-A copy

**Lemma 5.1.** In (4.1), `n<=1`.

Proof. If n>=2, rows two through seven form the exact word
A-B-C-A-B-C-A. The only possible exceptions from the table are
at rows zero, two, or three. Let p_f3 and p_g2 denote the
exceptional pair in those positions when present, and zero
otherwise. Row seven, an empty C->A row, gives `F_5=1`.
Row five gives

```
X F_5+p_f3+(x+x^2)p_g2=1.
```

Thus `p_f3+(x+x^2)p_g2=1+X`. But p_f3, if present, is
`X^-2+X^-1`, which is not `1+X`. The possible p_g2 is `1+X`,
whose product with x+x^2 has two nonzero residues. There is only
one pair, so both exceptions cannot occur together. None of the
possibilities satisfies the equation. QED.

### 5.2 At most seven D loops

**Lemma 5.2.** In (4.1), `r<=7`.

Proof. If r>=1, the final D-loop row is row -1. The g-mask in
row -2 is zero: it is either another D loop or the entering
E->D/F->D edge. Row one of (2.1), which is C->A, therefore gives

```
f_-1+p_f0=1.
```

The table makes f_-1 equal X when p_f0 is present and 1 otherwise.
In particular it is not the canonical monomial X^-1.

All r D-loop rows are clean and have g=0. Their entering E->D or
F->D row also has g=0. Thus the whole r-entry D word satisfies

```
X M_j+M_(j-1)+M_(j-2)=1
```

starting at its third entry, with no delayed-g endpoint term.
S8 Lemma 5.1, proved there by the exhaustive monomial-successor
identity, says that a word of length at least eight is entirely
canonical. This contradicts the final entry just found. QED.

Combining the lemmas gives `q=5+3n+r<=15`. Since a D->C cycle
has length at least five, the only binary height still possible
is q=8. Its two graph forms are `(n,r)=(1,0)` and `(0,3)`.
This is an analytic reduction from arbitrary binary height; no
height scan or spatial support bound was assumed.

## 6. Both remaining height-eight forms contradict actual positions

### 6.1 One A-B-C-A copy and no D loops

The eight edges are

```
DC, CA, AB, BC, CA, AB, BE/BF, ED/FD.
```

At row four, (2.1) reduces to `F_2+p_f3=1`. In the six-entry
table, the f-row-three case gives `X^-1+(X^-2+X^-1)=X^-2`,
and the g-row-zero or g-row-three cases give X^-1. Thus they
are impossible. Only the f-row-zero pair or the g-row-two pair
remains, and F_2=1.

At row five the equation is

```
X F_5+(x+x^2)p_g2=1.
```

The g-row-two case is impossible by its two nonzero residues.
In the f-row-zero case F_5=X^-1. But the residue-zero projection
of row seven forces F_5=1: its f_6 and f_7 have nonzero residues,
while g_4=g_7=0. This is a contradiction.

### 6.2 No A-B-C-A copies and three D loops

The edges are

```
DC, CA, AB, BE/BF, ED/FD, DD, DD, DD.
```

The residue-zero part of row five gives

```
X F_5+p_f3=1.
```

The f-row-three case would make X F_5 a sum of three distinct
monomials, impossible for the clean D-loop monomial F_5. With
that case excluded, the residue-zero part of row four gives
F_2=1. Again only the f-row-zero or g-row-two pair remains.

Let M_3,M_4 be the clean f monomials on the BE/BF and ED/FD
edges. They have the same nonzero residue. The nonzero-residue
part of row five is

```
M_4+M_3+(x+x^2)p_g2=0.
```

If p_g2 is present, its product has two nonzero residues and
cannot be canceled by M_3+M_4, which has only one. Otherwise
M_4=M_3. But row four's nonzero-residue part is
`X M_4+M_3=0`, forcing `(X+1)M_3=0` in the Laurent domain.
This too is impossible.

**Theorem 6.1.** For every binary height q, every compact actual
solution of syndrome N_q and weight q+2 has `L=0, sigma=2`.

Proof. The case L=2, sigma=0 is exhausted by the recurrent
components. Lemmas 3.1–3.2 exclude the D-only and two-cycle
possibilities. The pure A-B-C-A cycle has nonbinary length.
Lemmas 4–5 reduce the remaining component to the two cases just
excluded. Equation (2.3) leaves exactly L=0, sigma=2. QED.

## 7. Additional exact base-shell result: q=2 has no weight-four solution

Put `epsilon=1+t` at q=2, so `epsilon^2=0`, `N_2=epsilon`, and

```
P=x^3+epsilon,
Q=p+epsilon(x+x^2),       p=1+x+x^2.
```

Write `f=f_0+epsilon f_1`, `g=g_0+epsilon g_1`, where all four
coefficients are finite Laurent polynomials in x. The coefficient
of one in the syndrome equation gives

```
f_0=pH,       g_0=x^3 H
```

for a Laurent H. The epsilon coefficient then gives

```
x^3 f_1+p g_1=1+R H,       R=1+x+x^2+x^4+x^5.           (7.1)
```

If H=0, the physical weight is `2(|f_1|+|g_1|)`. Evaluating
(7.1) at x=1 makes the quantity in parentheses odd, so the
physical weight cannot be four.

If H is nonzero, aggregation of the two rows gives

```
4=|f|+|g| >= |pH|+|H| >=4.
```

The last inequality follows because p is not a Laurent unit,
so pH is never a monomial, and `|pH|+|H|` is even. Equality
requires H to be either `x^c` or `x^c(1+x)`: if H has at least
three terms, the sum is at least six; for a binomial, exactly
two cancellations in pH require its two exponents to differ by
one. This argument allows arbitrary integer c; no support
locations are scanned.

Equality in two-row aggregation also gives

```
supp(f_1) subset supp(pH),
supp(g_1) subset supp(x^3H).                            (7.2)
```

For `H=x^c`, the left side of (7.1) is supported only at
`c+3,c+4,c+5`. The right side has three distinct terms at
`c,c+1,c+2`, at most one of which can cancel with its constant
one. This is impossible.

For `H=x^c(1+x)`, the left side is supported at
`c+3,c+4,c+5,c+6`, whereas

```
RH=x^c(1+x^3+x^4+x^6).
```

Its lowest term must cancel the constant one, forcing c=0.
By (7.2), now f_1 is supported at `{0,3}` and g_1 at `{3,4}`.
The coefficients of x^4 and x^5 in `x^3 f_1+p g_1` are equal.
In the required right side `x^3+x^4+x^6`, they are different.
This final contradiction proves the shell is empty. QED.

## 8. Precise remaining sharp problem

By Theorem 6.1, at binary heights an excess-two solution has at
most one physical site in **every** block/row/x-residue cell.
There are no hidden cancelling pairs or three-site cells left.
Its graph has one of these two charged-edge patterns:

```
one edge of slack two, all other edges of slack zero;
two edges of slack one, all other edges of slack zero.
```

The fixed 128-edge graph contains 39 slack-two edges and 32
slack-one edges. These counts are not a spatial classification.
In particular, its six zero-slack transient edges

```
A -> (0,0,1),
(0,1,0) -> (1,0,zeta), (0,1,0) -> (1,0,zeta^2),
(1,0,1) -> (0,1,0),
E -> (0,1,zeta^2), F -> (0,1,zeta)
```

can now be part of a closed walk completed by charged edges.
Dropping them because they were transient at sigma=0 would be
incorrect. Nor does the remaining graph fix the integer x
positions of its occupied residues.

What remains to prove is that every actual Laurent realization
of these charged closed walks has at most three sites in
f-column zero and g-columns zero/one. This is a strictly reduced
branch of the original target, but the hole bound on that branch
is not established here. The local D-word or g-only burst bounds
do not alone provide its sharp charge, and the parent's coarse
finite-excess bound is not substituted for it.

Accordingly this note does not claim the canonical physical floor
for holes 4–6. It only removes the aggregation-hidden branch
from that remaining proof task and settles the stated base cases.

## 9. Verification and resource discipline

`a42_s11_codex_next_shell.py` checks only the 128 fixed graph
edges, their nonnegative slacks and recurrent components, the
five-edge minimum separation of D->C edges, the components after
deleting D->C, the six displayed Laurent correction identities,
and the scalar identities used at q=2. It does not enumerate
supports, cosets, heights, code families, or charged walks.

All prescribed assertions pass. The diagnostic body took about
0.0002 seconds and peak RSS was below 18 MiB. The script installs
55-second wall and 50-second CPU limits and asserts the 768 MiB
memory cap. No large computation, shared-file edit, old process,
git operation, memory edit, Lean work, or other-agent activity
was initiated. The universal conclusions above rest on the
written Laurent arguments, not extrapolation from the diagnostic.
