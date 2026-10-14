# A42 S13 — the single slack-two branch is empty above height one

2026-09-08. Bounded native Codex proof pass. Original Laurent pair,
actual unbounded x positions, no support search or periodic x folding.

## 1. Theorem and precise scope

**Theorem 1.1.** Let `q=2^a`, `a>=0`, and let f,g have finite
x support in `(F_2[t]/(t^q+1))[x^±]`. Suppose

```
P f+Q g=N_q,       P=x^3+t+t^2,
Q=1+(x+x^2)t^3,    N_q=1+t+...+t^(q-1),
|f|+|g|=q+2.
```

Assume x-mod-three parity aggregation loses no physical weight,
and the S7 16-state walk has exactly one edge of slack two,
with every other edge of slack zero. Then:

- If `q>=2`, there is **no such actual solution**.
- If `q=1`, its number of holes in the canonical support S_* is
  at most three, simply because its total weight is three.

Thus the requested `k<=3` statement is proved for this entire
single-charge branch. This does **not** address the two-slack-one
branch assigned to the parent/S14. Together with S11, it reduces
any remaining dyadic excess-two solution above height one to
that other branch, not to an established full next-shell theorem.

The q=2 case uses the analytic shell-emptiness theorem in S11 §7.
The proof below handles every dyadic q>=4 uniformly, with explicit
short cyclic identities where necessary. The parent's separate
row-compression theorem is not used to enumerate heights or
physical supports.

## 2. Conventions and the complete fixed graph input

Put `X=x^3`, retaining X as an unbounded Laurent variable. The
actual syndrome equation at each cyclic row is

```
X f_j+f_(j-1)+f_(j-2)+g_j+(x+x^2)g_(j-3)=1.             (2.1)
```

Pi_0 means projection onto x exponents divisible by three,
**retaining** those integer exponents. It does not set X=1.

Use numeric field labels `0,1,2,3` for `0,1,zeta,zeta^2`, and
concatenate the three state coordinates. Name

```
A=000, B=100, C=011, D=110, E=112, F=113,
U=001, V=012, W=013.
```

The complete zero-slack adjacency is

```
A -> U,B;              U has no zero outgoing edge;
002 -> 003;            003 -> 002;
010 -> 102,103;        101 -> 010;
C -> A;               B -> C,E,F;
D -> C,D;             E -> W,D;             F -> V,D.
```

States V,W,102,103,111 also have no zero outgoing edge. These
lists include all 17 zero-slack edges and all six transient
ones; none is silently discarded.

At L=0 an occupied residue is exactly one physical monomial,
not several sites with the same parity. The zero edges used in
the spatial proof have these exact residue supports:

| edge | f residues | g residues |
|---|---|---|
| AB | `{0}` | empty |
| BC | empty | `{1,2}` |
| BE | `{2}` | `{1}` |
| BF | `{1}` | `{2}` |
| CA | empty | empty |
| DC | empty | `{0}` |
| DD | `{0}` | empty |
| ED | `{2}` | empty |
| FD | `{1}` | empty |
| AU | empty | `{0}` |
| EW | empty | `{2}` |
| FV | empty | `{1}` |

A single charged edge `source -> target` can close only if its
source is zero-reachable from its target. Of the 39 slack-two
edges, exactly 18 pass this necessary test:

| charged edge | f residues | g residues |
|---|---|---|
| UA | empty | `{0,1,2}` |
| UU | `{1,2}` | `{0}` |
| UB | `{0,1,2}` | empty |
| `002->002` | `{1,2}` | `{1}` |
| `003->003` | `{1,2}` | `{2}` |
| `010->101` | `{0,1,2}` | `{0}` |
| VA | `{0,1}` | empty |
| WA | `{0,2}` | empty |
| BD | `{0}` | `{0,1,2}` |
| `102->010` | `{0,1}` | empty |
| `103->010` | `{0,2}` | empty |
| DE | `{2}` | `{0,2}` |
| DF | `{1}` | `{0,1}` |
| `111->111` | `{0}` | `{1,2}` |
| EC | `{0,2}` | `{0}` |
| EE | `{0}` | `{0,2}` |
| FC | `{0,1}` | `{0}` |
| FF | `{0}` | `{0,1}` |

The diagnostic verifies the exact edge sets and residue table by
the fixed Fourier conversion and S7 potential, not merely the
counts 39 and 18. The remaining 21 edges are excluded by this
zero-reachability test alone.

Throughout the rest of the proof the unique charged row is zero;
all other rows are zero-slack. Row indices are cyclic. The proof
does not identify distinct Laurent positions with the same residue.

## 3. Three spatial facts used below

### 3.1 A clean D->C window is impossible

If the five edges at rows `j-1,...,j+3` all have zero slack and
the edge at j is DC, they necessarily have the form

```
(D/E/F)->D->C->A->B->(C/E/F).
```

Indeed choosing AU instead of AB would leave no zero outgoing
edge for the next step. The row equations at j+2 and j+3 force
the AB monomial to equal X^-1 and 1, respectively, just as in
S7's spatial lemma. Therefore such a five-row window cannot be
entirely zero-slack in an actual L=0 solution.

It follows that every zero DC edge lies, relative to the charged
row zero, at one of

```
j=-3,-2,-1,1 modulo q.                                 (3.1)
```

The position j=0 itself is charged and is not DC in the table.
At q=4 some positions in (3.1) coincide; we do not treat them as
distinct there. In particular short wrapped cases are handled
separately rather than assuming a five-row window avoids charge.

### 3.2 Two consecutive clean ABC circuits are impossible

For the six zero edges A-B-C-A-B-C-A, numbered rows zero through
five locally, (2.1) at row three gives `X f_3=1`, because
`f_1=f_2=g_0=g_3=0`. At row five it gives `f_3=1`, because
`f_4=f_5=g_2=g_5=0`. This is impossible in the Laurent coordinate.

Consequently, if a zero return path has n consecutive complete
ABC circuits, `n<=1`.

### 3.3 A long pure D monomial word is canonical

For monomials M_0,...,M_(r-1) satisfying

```
X M_j+M_(j-1)+M_(j-2)=1       for 2<=j<r,                (3.2)
```

either all entries equal X^-1 or `r<=7`. This is the analytic
S8 monomial-word lemma. Briefly, equal consecutive exponents
force the next exponent to be -1. Unequal exponents admit a
successor only if one is zero, with next exponent the other
minus one. Apart from the fixed pair `(-1,-1)`, the longest
continuation is `2,0,1,0,0,-1,-2`, with seven entries.

When this fact is used below, all delayed g terms in (3.2) are
explicitly checked to vanish. It is not applied blindly at an
interval's first rows.

## 4. Seven edges are disposed of by length or the base shell

For UU and `111->111`, the target has no zero outgoing edge,
so a single-charge closed walk has length one.

For `002->002` or `003->003`, the other q-1 edges alternate
002 and 003. Returning to the starting state requires q-1
even, so q is odd; a dyadic q is therefore one.

For `010->101`, the only zero return is `101->010`.
For `102->010` or `103->010`, the only zero return reaching
the charged source is the corresponding direct edge from 010.
Any other choice ends at a state with no zero outgoing edge.
These three charged cases have length two, whose actual shell
is empty by S11.

This leaves exactly eleven charged edges for q>=4.

## 5. Charges UA, UB, VA, and WA

### 5.1 UA is impossible directly

The preceding edges are CA at row -2 and AU at row -1. Since
q>=4, the edge at row one must be AB: choosing AU would force
the next edge out of U to be charged again. Let its f monomial
be M, of residue zero. Row one of (2.1) gives `X M=1`.

At row two the edge is BC, BE, or BF, with no residue-zero
f or g. Projecting (2.1) onto residue zero gives `M=1`:
f_0 is empty and `(x+x^2)g_-1` has no residue-zero term.
The two equations contradict one another.

### 5.2 UB: all longer return paths collapse to one four-row identity

Again rows -2,-1 are CA,AU. The charged f_0 has one monomial
of each residue; denote its residue-zero monomial by M. The
residue-zero part of row one gives `M=1`.

If DC were at -3, the residue-zero part of row zero would give
`X M=1`, since g_-3 would have only residue zero. Thus DC is
not at -3. The other locations allowed by (3.1) are already
CA, AU, or an edge out of B. Therefore there is no DC anywhere.

With no DC, a zero path from B back to U can only use BC, CA,
ABC circuits, and the final AU. Entering D could never return
to B or U; entering V or W would stop. Hence

```
q=4+3n,
word = (charged UB), BC, CA, (ABC circuit)^n, AU.
```

By §3.2, n<=1. Of the lengths four and seven only four is
dyadic. At q=4 the rows are UB,BC,CA,AU.

Let g_3=M_3 be the AU monomial. Row three gives M_3=1.
Rows two and one give, successively,

```
f_0=1+(x+x^2)M_3=1+x+x^2,
g_1=f_0+1=x+x^2.
```

But row zero would require

```
X(1+x+x^2)+(x+x^2)^2=1.
```

Its left side is `x^2+x^3+x^5`, not 1. Thus UB is impossible.

### 5.3 VA and WA: the transient tails must also be included

For VA, put K=F and eta=1; for WA, put K=E and eta=2.
Write tau=3-eta. The preceding zero edges are

```
AB at -3,       BK at -2,       K(V/W) at -1.
```

The BK row has f residue eta and g residue tau. The final
transient row has empty f and one g of residue eta. The charged
row has f residues `{0,eta}` and empty g. The next edge must
be AB, since AU is a dead end for this charge.

Every candidate DC position in (3.1) is thus already a different
edge. With no DC, the zero return has just n ABC circuits followed
by AB, BK, K(V/W). Hence again `q=4+3n`, and §3.2 and dyadicity
reduce it to q=4.

The four rows are `(V/W)A, AB, BK, K(V/W)`. Let M be the
residue-zero monomial of f_0 and let H=f_1. The residue-zero
parts of rows zero and three give

```
M=X^-1,       H=1.
```

At row one, Pi_0 of `(x+x^2)g_2` is a single monomial T,
because g_2 has one nonzero residue. The same row equation gives

```
X H+M+T=1,       so T=1+X+X^-1.
```

This is not a monomial. Both transient-source charges are excluded.

## 6. Charges DE, DF, EC, FC, EE, and FF

These pairs of cases use only the distinction between residue
zero and the stated common nonzero residue. They do not assume
an axis reflection or a physical isometry swapping E and F.

### 6.1 DE and DF: row two is proved to be a D loop

Let K=E or F and rho be respectively two or one. For a charge
DK, the next edge must be KD: the other zero successor of K
is V or W, with no zero way back to D.

There can be no DC at a position in (3.1). A DC at -3 would
force AB at -1 and hence source B rather than D for the charge.
A DC at -2 would leave source A for the charge, and one at -1
would leave source C. At row one the source is K, not D.
Thus there is no DC anywhere.

The only possible return word is consequently

```
(charged DK), KD, DD, ..., DD.
```

In particular, **DD at row two is a conclusion, not an assumption**;
the possible immediate DC continuation has been excluded using
the actual spatial window.

At q>=4, let the nonzero-residue f monomials at rows zero and
one be M_0,M_1. All earlier g terms in rows one and two vanish.
Their residue-rho equations are

```
X M_1+M_0=0,       M_1+M_0=0.
```

They force `(X+1)M_1=0`, impossible for a nonzero monomial.

### 6.2 EC and FC are locally impossible

Let K=E or F and rho be its nonzero f residue. The preceding
edges are CA at -3, AB at -2, BK at -1; the next edge is CA.
The charge KC has f consisting of a residue-zero monomial A
and a residue-rho monomial M, and g of residue zero. Let T
be the f monomial at row -1, also of residue rho.

Row one gives `A=1` and `M=T`, because its delayed g term is
g_-2=0. Row zero's residue-rho part gives `X M+T=0`, because
g_-3=0 and its other terms have residue zero. Thus
`(X+1)M=0`, again impossible.

### 6.3 EE and FF have no actual nontrivial zero return

For a charge KK with K=E or F, the preceding zero edges are
CA,AB,BK at rows -3,-2,-1. A nonempty zero return must start
with KD at row one; its other successor is a terminal transient
state. All four positions in (3.1) are therefore different
from DC. No DC is possible.

But a zero return starting KD cannot reach K again without
passing through DC. This rules out every q>1 in these two
cases. The q=1 case is already covered by total weight three.

## 7. The remaining charge BD

Here f_0 is one residue-zero monomial A and g_0 is a sum G of
three monomials, one in each residue. The preceding edges are
CA at -2 and AB at -1. A zero return from D to B must use DC.
By (3.1), its only possible positions are

```
DC at 1 or at -3 modulo q.                              (7.1)
```

We first handle q>=8, where these two positions are distinct.
The q=4 wrapped case is treated separately in §7.4.

### 7.1 DC immediately after the charge is excluded, not assumed away

Suppose DC is at row one. The following edges at rows two and
three must be CA,AB, and row four is BC,BE, or BF. AU could
not return to the charged source B.

Row two gives `A=1`. If f_3=H, the residue-zero part of row
four gives `H=1`. Row three consequently gives

```
(x+x^2)G=1+X,
G=x^-1(1+x+x^2)=x^-1+1+x.                              (7.2)
```

The last implication is in the actual Laurent domain; it follows
by factoring `1+x^3=(1+x)(1+x+x^2)`. In particular Pi_0 G=1.

If DC is also at -3, the zero edge into D at -4 has g=0.
Thus row -1 forces its AB monomial J to equal X^-1. On the
other hand, the residue-zero part of row zero, using A=1,
Pi_0 G=1 and the residue-zero-only g_-3, forces J=X.
This is impossible.

If there is no DC at -3, the sole DC is the one at row one.
After CA the zero return to B can only consist of ABC circuits;
entering D again would require a second DC to get back to B.
The word therefore has length `q=4+3n`, with n<=1 by §3.2.
Neither possibility gives a dyadic q>=8. Thus row-one DC is
fully excluded for the heights in this subsection.

### 7.2 A late DC forces one long D block, but not a free endpoint

The only DC is now at -3. Since D's zero successors are only
D and C, the entire word is

```
row 0: charged BD;
rows 1,...,q-4: DD;
row q-3: DC;       row q-2: CA;       row q-1: AB.        (7.3)
```

Write r=q-4 and f_j=F_j on the D rows. Since q>=8, the CA
equation at row q-2 has vanishing delayed g term and gives

```
F_r=1.                                                 (7.4)
```

The D monomial word `F_2,...,F_r` satisfies (3.2): its recurrence
starts at row four, whose delayed term is g_1=0, and every later
delayed g term is also a D-row zero. The potentially nonzero
charged g_0 affects row three, which is **before** the recurrence
is imposed. This verifies the required endpoint condition.

If q>=16, this word has length `q-5>=11`, so §3.3 makes its
last monomial X^-1. This contradicts (7.4). Thus the only
remaining dyadic height in (7.3) is q=8.

### 7.3 The height-eight word is impossible by four monomial equations

Its rows are

```
BD, DD, DD, DD, DD, DC, CA, AB.
```

Let H=f_7. Rows seven and six give `H=X^-1` and `F_4=1`.
Row four then yields

```
F_3+F_2=1+X,       hence F_2 belongs to {1,X}.            (7.5)
```

Row one yields

```
X F_1+A=1+X^-1.
```

Both terms on the left are monomials, so the only choices are

```
(F_1,A)=(X^-1,X^-1) or (X^-2,1).
```

Row two, `X F_2+F_1+A=1`, then gives respectively
`F_2=X^-1` or `F_2=X^-3`. Neither belongs to the set in (7.5).
No spatial support search enters this two-term identity.

### 7.4 The wrapped height-four case

At q=4 the unique word is BD,DC,CA,AB. Let f_0=A, f_3=B,
g_0=G and g_1=C. Here A,B,C are residue-zero monomials and
G has one monomial in each residue. The actual cyclic row
equations, with no unwrapping assumption, give

```
A=1,        B=C,
G=1+X+(1+x+x^2)B,
X B+(x+x^2)G=1.
```

Eliminating G and writing p=1+x+x^2 yields

```
(x+x^3+x^4)B=1+x+x^2+x^4+x^5.                          (7.6)
```

The left side has exactly three distinct terms because B is a
monomial, whereas the right side has five. This contradiction
finishes BD, including the short case where row one and row -3
are the same DC edge.

## 8. Completion of the branch and the remaining program boundary

The fixed zero-reachability test excludes 21 of the 39 possible
charged edges. Section 4 disposes of seven more at dyadic q>=4.
Sections 5–7 exclude all eleven remaining charged edges using
actual Laurent row equations. This proves emptiness for q>=4.
The q=2 shell is already empty by S11, and q=1 has k<=3 by
weight alone. Theorem 1.1 follows.

No new finite-height search, support classification by bounded
enumeration, or larger-width transfer model was used. The finite
graph test only selects which **edge types** admit a zero return;
the proof treats arbitrary return lengths analytically. Powers
of two are used only in the explicit arithmetic length reductions.

The result is limited to L=0 with one slack-two edge. The other
possible excess-two charged pattern, two slack-one edges, remains
outside this note's scope. Consequently this note alone does not
claim full next-shell stability, full canonical HM, or an IBM
distance theorem.

## 9. Certificate and resource record

`a42_s13_codex_single_charge.py` verifies:

- all 128 edge slacks and the exact 17 zero edges;
- all 39 slack-two edges and the exact 18 zero-return candidates;
- their row-residue supports, retaining all transient states;
- the four prescribed length-four edge words used above;
- the finite Laurent identities in §§5 and 7, including the
  three-versus-five-term contradiction (7.6).

The verifier asserts the full expected edge sets, not just their
cardinalities, and emits PASS only after every assertion. It
does not enumerate heights, supports, cosets, or arbitrary words.
Its measured body took about 0.0002 seconds with peak RSS below
18 MiB; aggregate diagnostic time in this pass was below one
second. It installs 55-second wall and 50-second CPU limits and
checks the 768 MiB memory cap. Only the two assigned new files
were written; no shared files, git state, Lean, memory, old
processes, or additional agents were touched.
