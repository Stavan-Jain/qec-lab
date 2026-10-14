# A43-S2 Codex audit: crossing quantifiers and homology-sensitive bounds

2026-09-08. Scope: this new note and `scripts/a43_s2_codex_audit.py` only.
Research on `claude/a40-tour-de-gross-column`, PR #27. No Lean, integration,
registry changes, process resumptions, commits, or pushes.

## 0. Results and limits

1. **A43 L3 as quantified over every rectangular frame and every gap-dense
   representative is false.** On the original pair at `(ell,m)=(54,6)`, an
   explicitly specified nontrivial X-cycle has weight **48 < 54**, with
   longest empty-column/row runs **(3,1)**. Its nontriviality has both an
   explicit dual certificate and an independent binary rowspace check.
   An infinite family of these counterexamples is proved below.
2. These examples have an injective axis-cylinder lift and a weight-12
   representative in the same class. Thus they **do not refute** L3
   restricted to class-minimal representatives, S11's actual D sector,
   or the IBM bit family. Gap-density, D membership, and homological
   mixedness must remain different predicates.
3. The coordinate rotation used in A40 transports `(A,B)` to `(B,Abar)`,
   not back to `(A,B)`. Its cylinder floor needs its own proof or an
   explicit additional equivalence. The proposed literal axis symmetry
   fails on a weight-12 cycle. The fractional-packing cap is valid for
   the specified LP, but is not a no-go theorem for arbitrary LPs or
   homology-sensitive arguments; its b=0 arithmetic is also incorrect.
4. For the original pair, two explicit F4-valued homology invariants
   give proven, conditional-on-nonzero-invariant weight bounds. They
   are unaffected by stabilizer dressing and do not identify all mixed
   classes. They fall well short of the conjectured distance.
5. The parent has separately banked a connected, k=12 separable
   baseline in `IBM_EXACT_DISTANCE_PROGRAM.md` section 3. No duplicate
   tensor-product proof is needed here. The original L3 correction and
   restricted sufficient target are this note's integration outputs.

**Preserved upstream corrections.** The original IBM passage is the one
the parent reverified in arXiv:2506.03094v1, section 2.7 item 3; its Liang
et al. reference is arXiv:2503.03827. In the ORIGINAL family, b=1 has
k=12, whereas b=0 has k=26 when 127 divides r, as identified by A42
section 1.3. No argument here removes that exception or uses the wrong
global assertion that the plane scheme has dimension six/four zeros.
The F128 orbit matters. Finite common-zero support suffices for the
usual coprimality/plane Koszul argument; a k jump or a new species alone
does not disprove a distance formula.

## 1. Conventions and the quantifier audit

Work over F2, with

    A = 1 + y + x^3 y^-1,       B = 1 + x + x^-1 y^-3,
    Abar = 1 + y^-1 + x^-3 y,  Bbar = 1 + x^-1 + x y^3.

A pair of support polynomials `(U,V)` is an X-cycle precisely when
`Bbar U + Abar V = 0` modulo `x^ell=1`, `y^m=1`. X-boundaries are
`(Abar h,Bbar h)`. A Z-cycle `(Z0,Z1)` obeys `A Z0+B Z1=0` and pairs
with an X-chain by ordinary binary support intersection. These are
exactly A40-S11's READ/BND conventions, not an inferred transpose.

For arbitrary representatives, let G mean both longest cyclic empty
runs are at most three. Let W_x mean existence of an injective compact
cycle lift to `Z x Z_m`; W_y is analogous. Let D be the complement of
`W_x union W_y` in the nontrivial cycles. Finally let I_x and I_y be the
images on homology of the two cylinder projection maps, allowing compact
lifts of any weight. Then

    D is contained in G;
    v in W_x implies [v] in I_x;
    [v] outside I_x implies v outside W_x for EVERY representative.

The reverse implications are not automatic. In particular:

- G is a predicate on a support, not on a homology class.
- Membership in I_x permits a compact lift with cancellations and
  excess weight; it does not assert an injective lift of the given v.
- Being outside both axis images does not exclude an image from a
  twisted cylinder. The stronger S11 "mixed" definition excludes all
  cylinder images. Its finite three-subspace inventory is not assumed
  here, especially on the b=0, 127-divisible exceptional frames.
- Minimizing weight in a class can remove the hypothesis G, so a bound
  for class-minimal gap-dense vectors does not immediately imply the
  same bound for every gap-dense vector in the class.

These distinctions are explicit in A40 section 16.0 and are lost in the
stronger wording of A43 section 2.1. Its two forms can support the same
IBM lower-bound assembly **if both required cylinder floors are proved**.
They are not equivalent assertions over arbitrary rectangles: the
cylinder bound 2m need not reach ell+m-6 when ell>m+6.

## 2. An explicit counterexample, with finite identities

Use the following seed support, whose coordinates are `(block,x,y)`:

    L0 = {(1,0),(1,2),(2,1),(2,3),(3,1),(3,5)},
    L1 = {(3,1),(3,5),(4,0),(4,1),(4,2),(4,5)}.

Equivalently,

    L0(x,y) = x(1+y^2) + x^2(y+y^3) + x^3(y+y^5),
    L1(x,y) = x^3(y+y^5) + x^4(1+y+y^2+y^5).

Their Laurent syndrome has the **exact plane identity**

    Bbar L0 + Abar L1
      = (1+y^6)(1+x+x^3+x^4 y^-1+x^4 y^2).                (2.1)

Thus L is a compact weight-12 cylinder cycle on `Z x Z_6`.

Here is an independent dual tile, periodic in x with period six and
compact in y:

    Z0 = {(0,1),(0,2),(0,3),(1,2),(1,3),(2,2),
          (3,0),(3,2),(4,1)},
    Z1 = {(0,4),(1,2),(2,3)}.

On the plane its exact syndrome is

    A Z0 + B Z1 = (1+x^6)(x^-1 y+y^-1+y+x).               (2.2)

Consequently repeat this tile in x on any ell divisible by six, and
reduce y modulo six. The resulting vector Z is a torus Z-cycle. For
ell>=6, L meets this tiled Z in exactly the qubit `(0,1,2)`. In
particular `<L,Z>=1`, so L is nontrivial on these tori. Equations
(2.1)--(2.2), rather than samples, prove the uniform statement.

The X-boundary at `(c,0)` is the following six-qubit pattern, with rows
reduced modulo six:

    E_c,0 = {(c,0),(c,-1),(c-3,1)},
    E_c,1 = {(c,0),(c-1,0),(c+1,3)}.                     (2.3)

It is `(Abar x^c,Bbar x^c)` and has zero syndrome because the Laurent
polynomials commute.

**Theorem 2.1 (infinite counterexample family).** For every integer
t>=2, on `(ell,m)=(24t+6,6)`, put

    v_t = L + sum_{j=0}^{3t-1} E_{11+8j}.

Then v_t is a nontrivial, gap-dense X-cycle satisfying

    wt(v_t)=12+18t < 24t+6=ell+m-6.

**Proof.** All summands are cycles. They have disjoint occupied column
intervals: L is in columns 1--4, and the jth boundary occupies
`{8+8j,10+8j,11+8j,12+8j}`. Thus their qubit supports are disjoint and
the weight is `12+6(3t)`. Between consecutive intervals there are exactly
three empty columns; the last occupied column is `24t+4`, leaving two
empty columns before column 1 across the seam. Internal column gaps
are smaller. Hence gap_x=3. The occupied rows are exactly
`{0,1,2,3,5}`, so gap_y=1. Adding boundaries preserves the odd pairing
with Z and the class of L. Finally the weight deficit is `6t-6>0`.
Every displayed coordinate also gives a compact cycle on `Z x Z_6`;
columns 1 through `24t+4` project injectively to the torus. Therefore
v_t is in W_x, not D. Its class contains the weight-12 L, so v_t is
not class-minimal. QED.

At t=2 the complete specification is L plus the six generators with
centers `11,19,27,35,43,51`, on `(54,6)`. It has weight 48 instead of
the claimed floor 54. The diagnostic independently obtains rank 318
for the X-boundary rowspace and rank 319 after adjoining v_2, verifies
both syndromes, and verifies the odd pairing. At t=3 it gives weight
66 below 78 on `(78,6)`. These checks illustrate the theorem and do not
prove its quantifier extension.

**What is actually refuted:** the universal representative-level L3,
the identification G=D, and any assertion that adding boundaries
cannot fake two-axis gap-density. **What is not refuted:** a D-sector
floor, any class-minimal floor, or the original IBM distance law; these
rectangles have ell-m>6 and are outside the IBM bit family.

## 3. Axis transport and the LP claim

### 3.1 Axis transport changes the pair

Under `theta(x)=y^-1`, `theta(y)=x`, one has exactly

    theta(A)=B,       theta(B)=Abar.

The coordinate map `(a,b)->(b,-a)` therefore sends the AB code on
`(ell,m)` to the B-Abar code on `(m,ell)`. This is A40 section 13.1's
sound transport, not an AB self-symmetry. For the compact seed L,
theta(L) has zero syndrome for the transported pair on `Z_6 x Z`;
using the ORIGINAL pair instead gives a syndrome of weight **24**.
The diagnostic lists all its check coordinates. This is a concrete
refutation of the literal same-pair coordinate-rotation claim.

A43 section 1.3's one-line "x<->y, A<->B mirror" consequently does not
by itself prove the second cylinder floor. One must either prove
`floor_cyl^(B,Abar)(p)>=2p` or give another genuine weight-preserving
equivalence that supplies it. This audit does not rule out such an
additional equivalence or a separately valid mirrored proof.

The X/Z duality needed for code distance is different and valid:
`(U,V)->(Vbar,Ubar)` preserves weight and exchanges X cycles/boundaries
with Z cycles/boundaries. It does not exchange the two spatial axes.

### 3.2 The fractional-packing theorem has a narrow, sound scope

Fix an X-class c. For every Z-logical z pairing oddly with c, put a
nonnegative variable lambda_z, subject to
`sum_{z containing q} lambda_z<=1` for each physical qubit q. Then

    LP(c) <= n/d_Z,
    LP(c) <= d_X(c),

where d_X(c) is the minimum weight in c. For the first inequality,
sum all qubit capacities and use `wt(z)>=d_Z`. For the second, choose
a minimum representative v of c: every z intersects v at least once,
so `sum lambda_z <= sum_{q in v} load(q) <= wt(v)`. These are proofs.

The cap does not cover adding syndrome constraints, coset inequalities,
integer/parity information, or an exact convex hull of the feasible
representatives. For example, minimize the linear objective
`sum_q p_q` over the convex hull of binary supports of representatives
of c. This is an LP with optimum **exactly d_X(c)**, generally of
impractical size. It is an immediate counterexample to the blanket
"all LP bounds" reading. Nothing in the capacity proof says a sharp
argument must count physical qubits twice.

Substituting the CONJECTURED original distance into n/d_Z yields

    b=1: n/d_Z = 6(r+1),
    b=0: n/d_Z = 12r^2/(2r-1)
                 = 6r+3+3/(2r-1).

A43 has a factor-of-two error in the b=0 expression; the correct value
is **greater** than ell+3. The b=1 cap excludes reaching 12r by this
specific LP for r>=2 if that is the true d_Z; it does not exclude the
r=1 value. Substituting the conjectural d_Z is a conditional analysis
of the method, not an independently established numerical cap at all r.
Likewise on the square toric code, n/d_Z is 2L while its distance is L;
the statement that its distance equals this cap is not correct.

## 4. A homology-sensitive lower bound for the original pair

This section has no cylinder-floor hypothesis. Let both ell and m be
positive multiples of three, and let `omega^2+omega+1=0` in F4. Evaluation
at omega respects either cyclic period. For an X-cycle `(U,V)`, define

    T_x = x^2 U(x,omega) + omega(x+1)V(x,omega),
    T_y = omega y(y+1)U(omega,y) + V(omega,y).             (4.1)

These live in `F4[x]/(x^ell-1)` and `F4[y]/(y^m-1)`, respectively.

**Lemma 4.1.** Both T_x and T_y depend only on the X-homology class.
They obey `(1+x+x^2)T_x=0` and `(1+y+y^2)T_y=0`.

**Proof.** At y=omega,

    Bbar=x^-1(1+x+x^2),
    Abar=omega x^-3(x+1)(1+x+x^2).

Multiplying the cycle equation by x^3 gives the T_x assertion. At
x=omega, one instead has

    Bbar=omega(y+1)(1+y+y^2),
    Abar=y^-1(1+y+y^2),

and multiplication by y gives the T_y assertion. Substitute a boundary
`(Abar h,Bbar h)` into (4.1); the two terms in each expression cancel.
Thus the maps kill all boundaries. QED.

**Lemma 4.2.** Over F4, any nonzero length-s cyclic vector killed by
`1+z+z^2`, for 3 dividing s, has at least `2s/3` nonzero coordinates.

**Proof.** Its coefficients satisfy
`a_i+a_{i-1}+a_{i-2}=0`, hence `a_{i+3}=a_i`. Each period is
`(alpha,beta,alpha+beta)`, up to a cyclic shift. A nonzero such triple
has at least two nonzero entries. QED.

**Theorem 4.3 (class-sensitive bounds).** Write `w0=wt(U)`, `w1=wt(V)`.
Every representative of a class with T_x nonzero satisfies

    w0+2w1 >= 2ell/3.

Every representative of a class with T_y nonzero satisfies

    2w0+w1 >= 2m/3.

If both are nonzero, then

    wt(U,V) >= ceil(max(ell/3, m/3, 2(ell+m)/9)).         (4.2)

**Proof.** Evaluation along rows or columns cannot increase support
size, and multiplication by a binomial increases it by at most a
factor two. Apply Lemma 4.2 to (4.1) to obtain the first two
inequalities. Each weighted sum is at most twice the total weight,
and their sum is three times the total weight. QED.

For a nonvacuous exact control, the six-qubit original cycle on `(6,6)`

    U={(0,4),(2,1),(3,1),(4,4),(4,1),(5,1)}, V=empty

has both T_x and T_y of support size four. The script also verifies
that every padding generator in section 2 has zero invariants and
that L and v_t have identical invariants.

No completeness is claimed: vanishing T_x,T_y need not imply a
boundary, and nonvanishing is not asserted equivalent to S11 mixedness.
In particular these evaluations do not replace the F128 sector. At
near-square IBM periods, the joint bound is only about 8r/3. It is a
sound homology-sensitive estimate, not the missing sharp crossing floor.

## 5. The sound restricted original-IBM reduction

Let `T=ell+m-6` for `ell=6(r+b),m=6r`, r>=1, b in {0,1}. Assume
separately the two **lower** cylinder bounds appropriate to the original
and transported pairs, at least 2m and 2ell respectively. Define:

> **C_min, still open.** Every connected, class-minimal, nontrivial
> X-cycle on this IBM frame which has no injective compact lift to
> either axis cylinder has weight at least T.

Then every nontrivial X-cycle has weight at least T, hence the same
holds for Z by the block-swap/inversion duality.

**Proof of the reduction.** Choose a nontrivial cycle of globally
minimum weight. It is class-minimal. Make a graph on its qubits by
joining any two read by the same Z-check. Each support component is
itself a cycle, because no check meets two distinct components. If all
components were boundaries their sum would be one. Thus a nontrivial
component exists; two or more nonempty components would give a lighter
nontrivial cycle, contradicting global minimality. The minimizer is
therefore connected. An injective cylinder lift is nontrivial since
a boundary upstairs would project to one downstairs. The assumed
cylinder bounds give at least `min(2m,2ell)>=T`. Otherwise C_min applies.
QED.

This connectivity argument uses **global** distance minimality. A40
section 11.2 says "class-minimal" but then discards a component for
being a lighter logical, possibly in another class. That implication
does not follow for a prescribed class. A class-minimal representative
cannot have a nonempty boundary component, but it can in principle
contain several nontrivial components. Certificate claims using the
stronger prune need that distinction audited; this note does not rerun
or invalidate entire certificate populations from that observation.

For any class-minimal v and boundary e, the exact necessary condition
is `2|supp(v) intersect supp(e)|<=wt(e)`, from `wt(v+e)>=wt(v)`.
In particular, each weight-six generator meets it in at most three
qubits. The examples of section 2 meet each added generator in six
qubits and fail even this elementary minimality test. The converse of
this generator-overlap test is not claimed.

C_min is a sufficient restricted replacement for the false universal
L3. It has not been proved, and homology-sector floors alone do not
discharge pure-cylinder classes whose current representatives lie in D.
Restricting only to homologically mixed classes would miss those cases.
The general cylinder statements also remain conditional on the A42
work; the present audit does not duplicate its HM/seam proof attempt.

For the parent's toric-layer direction, the missing bridge is a
**uniform, weight-nonincreasing normal form in physical qubits**, with
explicit bounded network width and coverage of all relevant classes.
Only after such a theorem would a finite-width min-plus recurrence
constitute a complete distance proof. A phase/k classification or a
bounded-depth Clifford map alone supplies neither exact physical cost
nor this completeness theorem. No transfer scan is performed here.

## 6. Verification record and remaining work

Diagnostic: `scripts/a43_s2_codex_audit.py`, standard-library Python,
standalone, no sibling imports, SAT, census, random search, or external
data dependencies. Its two original instances are fixed t=2,3. Four
auxiliary separable-pair checks on `(6,6),(12,6),(12,12),(18,12)` remain
from the initial baseline exploration: commutation, k=8, cycle
membership, non-boundary membership, and witness weights 4,4,8,8.
These are not the parent's k=12 baseline and are not used to make a
distance claim here. The diagnostic does not enumerate distances.

Command from the repository root:

    uv run --no-cache --offline --frozen --no-sync --project experiments/bb_lab python -B experiments/bb_lab/scripts/a43_s2_codex_audit.py

The final diagnostic run reports about 0.0063 seconds of computation
and 18.2 MB peak RSS on macOS. It has a 55-second alarm, a 50-second
CPU cap, and a peak-RSS assertion below 768 MiB. Earlier attempts
stopped immediately on the default uv cache restriction, an unsupported
macOS address-space limit, and an erroneous expected seed/tile
intersection count (the correct count is **one**); no failed assertion
was consumed as evidence. Including those attempts and the direct
plane-syndrome identity check, total diagnostic computation is well
under one second, far below the 180-second aggregate budget. Read-only
inspection and proof writing are not computation sweeps.

Remaining original-IBM task: establish both applicable cylinder bounds
and **C_min**, or a weight-controlled homological/normal-form statement
implying it, with explicit coverage of the 127-divisible b=0 exception.
This note supplies no proof of the original conjectured distance beyond
the previously banked instances. The parent's separate baseline proof
is not reproduced. The two new files are the only outputs for parent
integration.
