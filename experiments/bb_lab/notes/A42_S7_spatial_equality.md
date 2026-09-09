# A42 S7 — restoring spatial location proves canonical-leader uniqueness

2026-09-08. Codex proof-research pass. New scope only; the preceding
S6 note and script were read but not modified.

## 1. Result and limits

**The canonical coset leader is unique at every binary-tower height.**
In fact, the theorem holds for every positive `q` with `3 ∤ q`, not
only `q=2^a`:

    P_0 f + Q_0 g = N_q,       |f|+|g|=q
        implies       (f,g)=(x^-3 N_q,0),

where

    P_0=x^3+t+t^2,    Q_0=1+(x+x^2)t^3,
    N_q=1+t+...+t^(q-1),
    f,g in (F_2[t]/(t^q+1))[x^±].

All chains have finite support in the unbounded `x` direction. The
S6 minimum-weight theorem `min(|f|+|g|)=q` is retained, with its
finite certificate recalled below. The new ingredient is a spatial
obstruction to an edge of its equality graph: two consecutive
syndrome rows require the same first-block monomial to equal both
`x^-3` and `1`.

This does **not** reuse the false small-SCC premise. Both the six-state
component containing five-cycles and the separate two-cycle are
retained and dealt with explicitly. The proof uses no finite-period
enumeration and makes no assumption that a class-minimal chain is
connected.

For the fixed omega support `S_*` from S6, the new consequences are:

- Every noncanonical error in its binary coset has weight at least
  `q+2`.
- With `k=|S_* minus s|`, all binary boundaries with `1<=k<=3` satisfy
  the desired physical floor `6q`; at `k=1,2` the floor is strict by
  four and two, respectively.
- The general fixed-support bound improves from `4q` to `4q+2`.
- The no-pair minimum remains exactly `6q`, and its binary boundary
  attaining that minimum is now unique.

**Still open:** full fixed-support weighted stability, proper-subset
HM beyond the proved strata, general HM for other omega supports, the
pure tower floor, and the remaining IBM cylinder/member assembly.
This note proves no exact-distance theorem for an infinite quantum-code
family.

## 2. The necessary finite graph, including all equality components

Let `q>=1`, `3 ∤ q`. Put `zeta^2+zeta+1=0` and `Phi=1+t+t^2`.
For a solution of `P_0 f+Q_0 g=N_q`, aggregate the `x` coordinates
modulo three by parity, separately in each block and `t` row. This
contracts physical Hamming weight.

The first aggregated binary triple at row `j` has Fourier coordinates
`(a_j,b_j) in F_2 × F_4`; the second has `(c_j,d_j)`. A triple's exact
weight is

    w(0,0)=0,       w(1,0)=3,
    w(0,b)=2,       w(1,b)=1       when b != 0.

The Fourier coordinates here mean augmentation and evaluation at
`x=zeta`. Thus a singleton with field coordinate `1,zeta,zeta^2`
is respectively in residue `0,1,2`. The pair with coordinates
`(0,1)` occupies residues `{1,2}`.

Evaluating the syndrome at `x=1` and `x=zeta` gives

    Phi a+c=N_q,       Phi(b+(1+t)d)=N_q.

Since `3 ∤ q`, `Phi` is invertible modulo `t^q+1`, and
`Phi N_q=N_q`. Therefore the cyclic row constraints are

    c_j=1+a_j+a_(j-1)+a_(j-2),
    b_j=1+d_j+d_(j-1).                               (2.1)

The 16-state graph has state `U=(a_prev,a_prev2,d_prev)`. Choosing
`(a,d) in F_2 × F_4` leads to `V=(a,a_prev,d)`, with cost

    c(U,V)=w(a,1+d+d_prev)+w(1+a+a_prev+a_prev2,d).

The full potential from S6 is as follows; columns give `d_prev`.

| `(a_prev,a_prev2)` | 0 | 1 | zeta | zeta^2 |
|---|---|---|---|---|
| `(0,0)` | 0 | 0 | 1 | 1 |
| `(0,1)` | 0 | 1 | 1 | 1 |
| `(1,0)` | 0 | 1 | 1 | 1 |
| `(1,1)` | 1 | 1 | 1 | 1 |

For all 128 edges,

    slack(U,V)=c(U,V)-1+F(U)-F(V) >= 0.               (2.2)

One short verification is to let
`L={(0,0,0),(0,0,1),(0,1,0),(1,0,0)}`. The only cost-zero edges
are `(0,1,1)->(0,0,0)` and `(1,0,1)->(0,1,0)`; both go from
potential one to zero. The outgoing cost lists from `L`, in choice
order `(0,0),(0,1),(0,zeta),(0,zeta^2),(1,0),(1,1),(1,zeta),(1,zeta^2)`,
are

    (0,0,0): 5,1,3,3,1,5,3,3
    (0,0,1): 3,3,3,3,3,3,3,3
    (0,1,0): 2,2,4,4,4,4,2,2
    (1,0,0): 2,2,4,4,4,4,2,2.

The two cost-one edges in this list stay in `L`. Every edge increasing
potential therefore has cost at least two; this proves (2.2).

Summing around a `q`-step closed walk yields cost at least `q`.
Consequently an actual weight-`q` error has equality in both

    physical weight >= aggregated weight >= q.       (2.3)

It follows that every edge of its walk has zero slack. Equality in
the first inequality also means that every occupied aggregated
residue lifts to **exactly one physical site**, and every empty
aggregated residue lifts to none: for a residue containing `n`
physical sites, the lost weight is `n-(n mod 2)`, a nonnegative
integer that must vanish individually.

### Complete zero-slack graph

Name the six recurrent states

    A=(0,0,0), B=(1,0,0), C=(0,1,1), D=(1,1,0),
    E=(1,1,zeta), F=(1,1,zeta^2).

The internal edges of this six-state component have the following
**exact row supports after equality in (2.3)**. A displayed residue
means one physical site somewhere in that residue, not a prescribed
integer `x` coordinate.

| edge | first-block residues | second-block residues | cost |
|---|---|---|---|
| `A->B` | `{0}` | empty | 1 |
| `B->C` | empty | `{1,2}` | 2 |
| `B->E` | `{2}` | `{1}` | 2 |
| `B->F` | `{1}` | `{2}` | 2 |
| `C->A` | empty | empty | 0 |
| `D->C` | empty | `{0}` | 1 |
| `D->D` | `{0}` | empty | 1 |
| `E->D` | `{2}` | empty | 1 |
| `F->D` | `{1}` | empty | 1 |

There is also the recurrent two-cycle

    (0,0,zeta) <-> (0,0,zeta^2),

with cost one on both edges. Its first block is empty on every row;
its second block has one site per row, alternately in residues two
and one.

The remaining six zero-slack edges are

    A -> (0,0,1),
    (0,1,0) -> (1,0,zeta),
    (0,1,0) -> (1,0,zeta^2),
    (1,0,1) -> (0,1,0),
    E -> (0,1,zeta^2),
    F -> (0,1,zeta).

None is in a directed cycle, as is immediate from the complete list.
These 17 edges exhaust the zero-slack possibilities by direct
substitution in the eight-entry choice table at each of 16 states.
The diagnostic independently checks this finite statement.

The five-cycles `A-B-E-D-C-A` and `A-B-F-D-C-A` remain genuine
zero-slack graph cycles. Adding self-loops at `D` gives length-eight
equality walks. The next lemma, not their lengths, excludes actual
spatial lifts involving them.

## 3. Spatial obstruction: an actual equality walk cannot use D->C

Write `f_j(x),g_j(x)` for the actual Laurent polynomials in row `j`.
The unaggregated syndrome equation reads, for every cyclic row,

    x^3 f_j + f_(j-1) + f_(j-2) + g_j
        + (x+x^2)g_(j-3) = 1.                        (3.1)

Let `Pi_0` extract the terms whose `x` exponents are divisible by
three, retaining their actual exponents. Thus

    Pi_0(x^-3)=x^-3,       Pi_0(1)=1.

This is not aggregation modulo `x^3+1`. Equivalently one may use the
direct-sum decomposition

    F_2[x^±] = F_2[X^±] + x F_2[X^±] + x^2 F_2[X^±],
    X=x^3,

while preserving the Laurent variable `X`.

**Lemma 3.1.** A weight-`q` actual syndrome solution cannot have
transition `D->C` in its zero-slack closed walk.

Proof. Put that transition at row zero. Because the walk is closed,
all its vertices lie in a recurrent component. By the complete
edge list its predecessor and successors must have the form

    U -> D -> C -> A -> B -> V,
    U in {D,E,F},       V in {C,E,F}.

Label the five edges by rows `-1,0,1,2,3`. The residue table supplies
the following physical information:

    g_-1=0;
    f_0=0,        supp(g_0) is in residue 0;
    f_1=g_1=0;
    g_2=0,        supp(f_2) is in residue 0;
    supp(f_3), supp(g_3) are in residues {1,2}.

Empty supports in the last line are permitted. These facts hold
for every choice of `U,V`, every placement of the occupied sites,
and any number of earlier or later loops or excursions.

At row two, (3.1) becomes

    x^3 f_2 = 1,

since `f_1,f_0,g_2,g_-1` vanish. Hence `f_2=x^-3`.

At row three, apply `Pi_0` to (3.1). The terms `x^3 f_3` and `g_3`
have no residue-zero part. The term `f_1` vanishes. As `g_0` has
only residue-zero exponents, `(x+x^2)g_0` also has no residue-zero
part. The only surviving term is `f_2`, giving

    f_2 = 1.

This contradicts `x^-3 != 1` in the unbounded Laurent coordinate.
No bound on the actual positions or a finite-period search enters
the argument. ∎

This subsumes the parent's hand obstruction for the single five-cycle:
it does not need to solve the locations of either `g` site or the
other `f` sites. It excludes every longer equality word containing
`D->C`, including arbitrary self-loop padding and concatenated
three- and five-cycle excursions.

The lemma uses the original spatial row equations only after proving
equality in physical Hamming contraction. Without that equality,
an empty aggregated residue could hide cancelling physical pairs,
so this row-support argument would not apply.

## 4. Excluding the other recurrent component

**Lemma 4.1.** No compact solution with weight `q` can realize the
two-cycle `(0,0,zeta)<->(0,0,zeta^2)`.

Proof. Its first block is empty in the aggregate; equality in (2.3)
forces `f=0` physically. Thus `Q_0 g=N_q`. If `g` is nonzero with
extreme occupied `x` columns `c_min,c_max`, multiplication by `Q_0`
has lowest column `c_min` and highest column `c_max+2`. The
coefficients at those extremes are respectively multiplied by the
units `1,t^3`, so neither disappears, even in the presence of
nilpotents in `F_2[t]/(t^q+1)`. Both cannot be the single column zero
of `N_q`. The zero chain also fails its syndrome equation. ∎

This is the already-known extrema exclusion, now included explicitly
to make the all-height argument complete.

## 5. All-height uniqueness theorem

**Theorem 5.1.** For every positive `q` with `3 ∤ q`, the minimum
physical Hamming weight among compact solutions of
`P_0 f+Q_0 g=N_q` is `q`, and its unique minimizer is

    e_0=(x^-3 N_q,0).

Proof. Equation (2.2) and physical Hamming contraction give the
lower bound `q`. The displayed chain has that weight and that
syndrome because `(t+t^2)N_q=0`.

For any solution attaining the bound, the closed walk lies in the
zero-slack graph. Lemma 4.1 excludes the two-cycle component.
Lemma 3.1 forbids `D->C` in the six-state component. After deleting
that edge, the only recurrent components left there are the
self-loop at `D` and the directed three-cycle `A-B-C-A`. A closed
walk in the latter has length divisible by three and cannot have
length `q`. Thus every row of the equality walk is `D->D`.

Its aggregated first block has one site in residue zero per row;
its second block is empty. Equality in aggregation implies `g=0`
and one physical first-block site per row. It remains to solve
`P_0 f=N_q`. Multiplication by `P_0` is injective on compact Laurent
chains: for any nonzero input, its highest occupied column gives
a nonzero highest output column three steps later, since the
coefficient of `x^3` is one. As `x^-3 N_q` is already a solution,
it is the only solution. ∎

The proof works over the actual finite row ring, including arbitrary
two-adic multiplicity. It needs no global claim about how many plane
zeros exist, no order-127 exclusion, and no homology classification
of the binary complex. In particular it proves the original S7
objective for every `q=2^a`, `a>=0`.

The condition `3 ∤ q` is used in deriving (2.1). No uniqueness claim
for `3 | q` is inferred from this proof.

## 6. Consequences for the canonical omega support

Let

    B_q={ (Q_0 h,P_0 h) : h compact },
    1_(S_*)=(N_q,(1+x)N_q).

Then `|S_*|=3q` and `D_0(1_(S_*))=N_q`. The quaternary cycle
`mu_*=(N_q,(1+x)N_q)` is nontrivial for

    P_w=x^3+zeta t+zeta^2 t^2,      Q_w=Q_0.

Indeed `t N_q=N_q` reduces its cycle equation to
`(x^3+1)+(1+x+x^2)(1+x)=0`; a boundary would require
`Q_w h=N_q`, which the same two-extreme-column argument excludes.
These statements do not assume a minimum-support theorem for all
omega cycles.

For any binary boundary `s`, set

    e=1_(S_*)+s,   E_err=|e|,
    k=|S_* minus s|=|e intersect S_*|,
    e_out=|s minus S_*|,       E_err=k+e_out.

In the period-`3q` physical inverse Fourier lift with omega component
`mu_*` and binary component `s`, the fibre table gives exactly

    |v|=3q+k+3e_out=3q+3 E_err-2k.                    (6.1)

This remains true for any nonzero scalar multiple of `mu_*`.

**Corollary 6.1 (strict gap away from the leader).** Every error of
syndrome `N_q` distinct from `e_0` has weight at least `q+2`.

Proof. Theorem 5.1 excludes weight `q` and lower. Augmentation at
`x=t=1` evaluates both `P_0,Q_0` to one, so every such error has
weight congruent to `q` modulo two. ∎

The leader `e_0` has no coordinates in `S_*`. Consequently `k>0`
forces `e!=e_0`, giving the now-unconditional estimate

    |v| >= 6q+6-2k       whenever k>0.                (6.2)

In particular the canonical-support physical floor is

| holes `k` | proved weight lower bound |
|---|---|
| 0 | `6q` |
| 1 | `6q+4` |
| 2 | `6q+2` |
| 3 | `6q` |

These hold at every `3 ∤ q`, without restricting the boundary's
width, number of components, or other coefficients.

**Corollary 6.2 (all boundaries, fixed support).** For every such `s`,

    |v| >= 4q+2.

Proof. If `e=e_0`, (6.1) gives `|v|=6q>=4q+2`. Otherwise
`E_err>=q+2` and `k<=E_err`, so (6.1) gives
`|v|>=3q+E_err>=4q+2`. ∎

**Corollary 6.3 (unique no-pair minimizer).** If `k=0`, the physical
minimum is exactly `6q`, and the unique binary boundary attaining it is

    s_* = ((1+x^-3)N_q,(1+x)N_q).

Proof. Equation (6.1) with `k=0` attains `6q` precisely when
`E_err=q`. Theorem 5.1 forces `e=e_0`, hence the displayed `s_*`.
It is a boundary: take `h_*=(x^-3+x^-2)N_q`. Multiplication gives
`(Q_0 h_*,P_0 h_*)=s_*`. Its support contains `S_*` and has exactly
`q` extra fibres. ∎

These are statements about a fixed omega component and its allowed
binary mixing. They neither bound every omega cycle nor address
S11 mixed torus homology classes.

## 7. Exact remaining obligation

The full fixed-support HM bound is equivalent to

    k <= (3/2)(E_err-q),                             (7.1)

for every `e=1_(S_*)+s`, `s in B_q`. Uniqueness proves only that
`E_err-q>=2` whenever `k>0`. It closes (7.1) for `k<=3` but does
not control larger `k`. In particular an additive weight gap is
not the required bound proportional to the number of holes.

For the binary tower one may equivalently minimize over all compact
solutions of the syndrome equation, because `P_0` is a unit there
and all binary cycles are boundaries. Explicitly,

    f=x^-3 N_q + C_q g,       C_q=Q_0 P_0^-1.

The remaining target is a quantitative spatial bound for this
actual syndrome coset, not for its two compressed check values alone.
As already audited in S6, the proper subset `(0,(1+x)N_q)` has zero
values under both compressed checks but nonzero actual syndrome
`(1+x^3)N_q`. Thus the checks alone do not settle all-subset HM.

For errors above weight `q`, there are two possible losses that the
new equality argument does not yet charge quantitatively:

- Aggregation can cancel physical pairs in the same residue.
- The closed graph walk can have positive total slack.

Either allows extra terms in the two rows used in Lemma 3.1.
Charging holes against those losses at the constant in (7.1) is a
concrete next lemma; it is not proved here. General HM additionally
ranges over other omega supports and their excess, and the full
cylinder floor still also needs the independent pure-floor work.

## 8. Verification and handoff

The new diagnostic is
`scripts/a42_s7_codex_spatial_certificate.py`. It uses only Python's
standard library. It independently reconstructs the eight Fourier
triples and all 128 graph edges, verifies the Bellman inequalities,
checks all 17 zero-slack edges and recurrent components, and verifies
the nine possible predecessor/successor contexts of Lemma 3.1.
Deleting the forbidden edge leaves exactly the claimed self-loop,
three-cycle, and pure-g two-cycle as recurrent components. The final
Laurent contradiction is the formal inequality `x^-3 != 1`, proved
above; no finite window substitutes for it.

The diagnostic ran once in **0.000163 seconds**, with measured peak
RSS **16,531,456 bytes**. Its built-in alarm is 55 seconds and its
measured-RSS assertion is 768 MiB. No period enumeration, SAT,
census, lifted-vector search, old-process control, Lean work, or git
mutation was performed. Only this new note and its new diagnostic
were written.

Parent handoff: bank Theorem 5.1 as an unconditional all-height
uniqueness theorem, and Corollaries 6.1–6.3 with their fixed-support
scope. Keep (7.1), general HM, the pure tower, and the full IBM
distance law open.

## 9. Independent audit of GROSS_SURFACE_PRODUCT.md, §§2–6

Parent-requested proof audit, 2026-09-08. The construction note and
the seed/check conventions in `scripts/gross_surface_product_check.py`
were read. No diagnostic was executed and neither file was edited.
**Verdict: the stated product construction and its universal distance
proof are sound, given the banked gross seed distances.** This is
the known augmented-surface construction specialized to the retained
gross presentation, not a novelty claim or completion of the original
IBM trinomials problem.

### All-class metric statement, not only existence of one large class

For a nonzero class `alpha in H_1(C)`, let `delta_X(alpha)` denote its
minimum physical Hamming weight. In `D=C tensor R_t`, the degree-two
cycle equations are exactly

    d_1 u_j=0,       u_j+u_(j+1)=d_2 v_j.

Every slice has the same class `alpha`. Hence every representative
of the associated class in D has weight at least
`sum_j |u_j| >= t delta_X(alpha)`. Repeating a minimum representative
of `alpha` in every slice, with all `v_j=0`, attains equality.

The common-class map has zero kernel on homology even with dependent
X generators: when the class is zero, remove the `u_j` using arbitrary
preimages under `d_2`; the remaining `v_j` lie in `ker d_2`. Their
successive partial sums supply a boundary that removes them while
keeping every `u_j` zero. This proves the required class detection
without assuming `d_2` injective.

For `E=D tensor I_t`, vertex augmentation is a homotopy equivalence
and contracts Hamming weight. Its inverse on homology is inclusion
at one vertex, which preserves Hamming weight. Consequently it
preserves the minimum weight of **each individual** degree-two
homology class, not just the least minimum over all classes. Thus

    delta_X^E(alpha)=t delta_X^C(alpha)

under the displayed homology identification. No connectivity or
class-minimal-implies-connected assumption is involved.

### Reversed dual and Pauli distance

The full complex has degrees zero through four, so reversing and
dualizing preserves the qubit degree two. Its factorization is

    E^vee = C^vee tensor I_t tensor R_t,
    (C^vee)_j=C_(2-j)^*.

Indeed `R_t^vee=I_t` and `I_t^vee=R_t`, with the indicated adjacent
maps transposed. The central homology of this reversed complex is
the cohomology determining Z logicals in the parent's convention.
Exchanging the last two factors is a basis permutation, so the
same argument gives per-class multiplication for Z as well. The
fixed gross inputs `d_X=d_Z=12` therefore imply both sector minima
are `12t`. A nontrivial CSS Pauli has a nontrivial X or Z component,
whose support is contained in its total support, establishing full
distance `12t` without an unproved product-distance assertion.

### Dimensions with every redundancy retained

For this BB seed, spatial inversion and block exchange identify
the X and Z matrices up to row and column permutations. Their ranks
are therefore equal; the banked `n=144,k=12` gives rank 66 each.
Thus `dim H_0(C)=dim H_2(C)=6`, with `dim H_1(C)=12`.

Since R has only one-dimensional degree-one homology and I has only
one-dimensional degree-zero homology, the product has

    dim H_0(E)=dim H_4(E)=0,
    dim H_1(E)=6,       dim H_2(E)=12,       dim H_3(E)=6.

In particular the seed's extreme homology does not enter the qubit
homology in degree two. Explicit dimensions and ranks give a second
way to track all dependencies:

    dim E_0 = dim E_4 = 72t(t-1),
    dim E_1 = dim E_3 = 72(2t-1)^2,
    n=dim E_2 = 144(3t^2-3t+1),
    rank d_2 = rank d_3 = 72(3t^2-3t+1)-6.

Hence each displayed CSS check set has `72t(t-1)+6` independent
linear relations, and `n-rank d_2-rank d_3=12`. At `t=1` this
retains exactly six dependencies in each gross check set. No
independent-check hypothesis was silently imported.

### Sparsity and connectedness

Independent degree bookkeeping confirms the three X-source check
weight bounds `8,8,7`, the three Z-target bounds `7,8,8`, and the
four qubit participation lists

    X: 5,4,6,5;       Z: 5,6,4,5.

The tensor terms land in distinct summands, so those counts are
valid upper bounds; missing endpoint terms only reduce them.

The particular seed supports in the diagnostic contain the origin
in both blocks, a unit `y` step in the first block, and a unit `x`
step in the second. Translated X checks therefore connect the two
blocks at every site and generate both translations, proving seed
incidence connectedness. In the product, Z checks of degree type
`(1,0,0)` connect adjacent R positions in the central array, while
X checks of type `(1,1,1)` connect adjacent I positions. The other
three qubit types attach to that array via the checks specified in
the parent note. The argument also handles `t=1`, when those extra
types are absent. Thus the displayed family is connected and LDPC.

No correction to §§2–6 is required by this audit. Its theorem is
the actual gross-derived family `[[144(3t^2-3t+1),12,12t]]`; the
original IBM-cylinder stability and distance objectives remain
separate and open. S8's spatial-stability work was not duplicated.

## 10. Independent audit of A40_PARITY_TOR_PROOF.md

Parent-requested bounded proof audit, 2026-09-08. The actual proof
note and `scripts/a40_parity_tor_check.py` were read. **Verdict: the
universal class-space theorem, cylinder-image equality, and cover
maps pass; no correction is required.** The original-mirror note
was also reviewed separately and passes without a self-mirror
hypothesis. Neither parent note nor diagnostic was edited.

### Plane module and exceptional summand

The reduction of Q modulo `P=y^2+y+x^3` and its multiplication
matrix in the basis `1,y` are correct. In particular

    a0/g=1+x+x^3,       b0/g=x+x^3,

whose sum is one, proves the entry gcd exactly. The determinant
`g^3 h` then gives Smith invariants `g,g^2 h`, not a globally
six-dimensional quotient. The two primary summands are preserved
by y because multiplication by y commutes with x. Thus M0 has
dimension six and is killed by `g^2`, while W has dimension seven.
On W, y acts as a scalar: every L-linear endomorphism of the cyclic
module `L/(h)` is multiplication by an element of `L/(h)`.

The degree-seven irreducibility certificate is sufficient: since
seven is prime, `h | x^128+x` and `gcd(h,x^2+x)=1` exclude every
proper factor degree and repeated factor. Both x and y on W are
nonzero and not one, so their order is 127. On M0,

    x^6+1=(x+1)^2 g^2=0,
    (y^2+y+1)^2=x^6+1=0,
    y^6+1=(y+1)^2 (y^2+y+1)^2=0.

These identities prove period annihilation for all multiples of
six; no local-multiplicity stabilization is needed. W survives a
rectangle precisely when 127 divides both lengths. Consequently
`k=12+14*[127 divides r]` for b=0 and `k=12` for b=1. This
exception is the established A42 result, not a newly discovered
distance obstruction.

### Balanced Tor is natural for the actual physical maps

The nonzero determinant makes Q injective on the free L-module
`R/(P)`. Together with P being nonzero in the domain R, this proves
that `K(P,Q)` resolves M. The rectangular period Koszul complex
also resolves its quotient; monicity in the separate variables
proves regularity despite repeated period factors.

For an explicit naturality justification, let C be `K(P,Q)` and
let E be a free period resolution of an algebra A. The two
augmentations of the bounded double complex `C tensor E` give
quasi-isomorphisms to `C tensor A` and `M tensor E`. For every
chain map between period resolutions, these augmentations commute
with that map. Thus the balanced Tor identification used in the
parent proof intertwines the physical quotient and norm maps with
the displayed maps on `M tensor E`; it is not merely a dimension
comparison. Nothing in this argument makes it a Hamming isometry.

### Cylinder-image equality, including invisible W classes

For all integer windings, including negative ones,

    X^a Y^b+1=Y^b S_a(X)(X+1)+S_b(Y)(Y+1)

is an identity of Laurent polynomials, with `S_n(1)=n mod 2`.
Hence the column `(U,V)^T` really lifts cylinder projection, and
its induced Tor map is `z -> (Uz,Vz)` on `ann_M(h_u)`.

On the surviving summand `M_Lambda`, both X and Y act as one.
Every element of this summand is in the source annihilator, so
the image contains every `(a z,b z)`. Conversely, any source
component in W when W is absent from `M_Lambda` lands in an
acyclic target summand and contributes zero. Primary components
cannot mix under these R-linear maps. This proves equality of
the image, not just an upper bound on it. A primitive winding has
nonzero parity, giving exactly the two axes and the diagonal.
The stated class counts for dimensions six and thirteen follow.

### Cover direction, regularity, and coordinate change

For a lattice basis `(u,v)`, set `F=z^u+1`, `G=z^v+1`.
One elementary justification of both period regularity claims is
that R is free over `F2[(z^u)^+/-1,(z^v)^+/-1]`, with a basis
of lattice-coset monomials. The regular sequences `(F,G)` and
`(F,G^2)` in that subring remain regular under this free extension.

The norm goes from the base to the cover, and is multiplication
by G, not the quotient map. It is well defined since it sends
the base relations into `(F,G^2)`. Its two translates have
disjoint supports on the two sheets, so it doubles physical
Hamming weight exactly. The chain lifts, in degrees zero through
two, are correctly oriented:

    base -> cover: (G, diag(G,1), 1),
    cover -> base: (1, diag(1,G), G).

All four adjacent chain-map identities hold. W cannot become
newly active under this cover: a nonzero G scalar in F128 remains
nonzero when squared. On the common surviving summand, tensoring
the lifts gives `diag(0,1)` and `diag(1,0)` on degree-one homology.

To identify their kernels in the original rectangular coordinates,
write `u=(a ell,b m)` and `v=(c ell,d m)`. The period-basis change
on this homology evaluates to the matrix with columns `(a,b)`
and `(c,d)` modulo two. Its determinant is one because the
integer determinant is plus or minus one. Thus the first axis in
the `(F,G)` coordinates is exactly `W_u` in the original ones.
The pullback kernel and projection image are therefore both W_u,
with dimension `d0`; pullback has rank `d0` as claimed.

### Metric scope and bounded verification

Every mixed class survives each of the three indicated pullbacks.
Applying a cover's sector distance to the doubled representative
proves precisely the factor-one-half comparison in (5.2). It does
not prove a numerical cover distance, an injective cylinder lift
of a minimum representative, or a footprint-connectivity claim.
The class-space label "mixed" here is the S11 winding distinction,
not A42's binary/omega fibre distinction. The original IBM metric
closure and the quantitative S8 task remain untouched.

One independent fixed-size calculation checked the polynomial
quotients, gcd, determinant, seven Frobenius squarings modulo h,
and `g^2 | x^6+1`. It took **0.013543 seconds**, with measured
peak RSS **20,889,600 bytes**, under the 55-second and 768-MiB
limits. No parent diagnostic was rerun. The winding and cover
arguments above were audited algebraically, not inferred from
the diagnostic's finite winding examples. Only this appendix was
written during the review; no searches for vectors or periods,
process control, or git mutation occurred.
