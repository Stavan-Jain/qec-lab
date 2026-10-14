# A44-S3: a finite diagonal transporter and an exact seam obstruction

2026-09-08. Codex sidecar, original IBM square column. This note and
`scripts/a44_s3_codex_seams.py` are the only S3 outputs.

## 0. Result and limitations

The target `12r-6` is **not attained here**. Instead, the following
particular linear-cost construction is completely resolved, without
using a conjectural metric floor or a finite-size extrapolation.

* An explicit six-qubit plane transporter has syndrome
  `(1+(xy)^2)E`. Its endpoint charge is nonzero and has no F128
  component. Repeating this transporter gives a nontrivial original
  square logical of **exact weight `18r` for every integer `r>=1`**.
* The infinite cylinder with period vector `(2,2)` has minimum nonzero
  X-cycle weight **exactly six**, also attained nontrivially by this
  transporter. The proof is a Laurent-polynomial calculation.
* On the original `(6r,6r)` square, **every nonzero X-cycle invariant
  under the physical translation `(2,2)` has weight at least `12r`
  for `r>=2`**. For **every `r>=4`**, the minimum among these cycles,
  and separately among the nontrivial ones, is **exactly `18r`**.
  The same conclusions hold on the Z side by inversion and block swap.

Thus fixed diagonal-period-two seams cannot produce IBM's square
discount, even after arbitrary corrections that preserve this symmetry.
This excludes a precisely specified ansatz, not the IBM distance law,
general linear-cost seam networks, or S1's witnesses with growing
diagonal period. The new `18r` family is not a better unrestricted
upper bound than the banked L12 bound `12r`.

The only diagnostic expands fixed Laurent identities, checks the
constants of the elementary lower-bound proof, and verifies one
specified `(24,24)` vector. There is no SAT, coset search, family census,
matrix minimization, or search over seam shapes.

## 1. Coordinates and the finite transporter

Use the original reciprocal pair

```
alpha=1+y^-1+x^-3 y,       beta=1+x^-1+x y^3,
partial_1(f,g)=beta f+alpha g,
partial_2(h)=(alpha h,beta h).
```

Set `z=x`, `tau=xy`. This is an integral unimodular coordinate
change, not a change of the original polynomials or code:

```
alpha=1+z tau^-1+z^-4 tau,
beta=1+z^-1+z^-2 tau^3.                              (1.1)
```

The original square of side `N=6r` has group algebra

```
T_N=F2[z^+/-1,tau^+/-1]/(z^N+1,tau^N+1).
```

Indeed `x^N=z^N` and `y^N=tau^N z^-N`. A physical diagonal
translation `(2,2)` is multiplication by `tau^2`.

Define the six-qubit pair

```
P_0=z^-4+z^-2+tau(z^-1+z^-2+z^-3),
P_1=z^-2,                                          (1.2)
```

or, in original coordinates,

```
P_0=x^-4+x^-2+y+x^-1 y+x^-2 y,
P_1=x^-2.
```

Put

```
L=z^-3+z^-4+z^-5,
E=(1+tau^2)L+tau(z^-4+z^-6)+z^-1 tau^-1.             (1.3)
```

The exact plane identity is

```
beta P_0+alpha P_1=(1+tau^2)E.                     (1.4)
```

For a short hand check, the left side collects as

```
(1+tau^4)L+(tau+tau^3)(z^-4+z^-6)
                     +(tau^-1+tau)z^-1,
```

which factors as stated. There are no period reductions in (1.4).

Let `omega^2+omega+1=0`. At

```
(z,tau)=(omega,1), equivalently (x,y)=(omega,omega^2),
```

both checks vanish and

```
E=omega^-4+omega^-6+omega^-1=1.                    (1.5)
```

In particular E is a nonzero plane charge, not the syndrome of a
finite plane pair.

### 1.1 Why this charge is ordinary, not an illicit global annihilator

No assertion that `x^6+1` or `y^6+1` annihilates the entire plane
module is used. Equation (1.4) explicitly says that this particular
charge is annihilated by `1+tau^2`.

Using the banked plane-module decomposition only to name the charge
sector, the exceptional summand is F128, with z of order 127.
On this summand tau cannot equal one: at tau=1,

```
beta=z^-2(z^2+z+1),
```

and an order-127 z cannot solve `z^2+z+1=0`. Since every nonzero
F128 element has order dividing 127, `tau^2` is also not one.
Consequently `(1+tau^2)[E]=0` forces the F128 component of E to
vanish. Equation (1.5) certifies a surviving ordinary component.

All cycle, weight, and nontriviality proofs below are independent of
that decomposition. In particular they apply also when `127|r`.
The original dimension exception remains `k=12+14*[127|r]`; nothing
here turns those members into k=12 codes or refutes their conjectured
distance because of their extra logical dimensions.

## 2. Exact distance six on the infinite diagonal cylinder

Let

```
C=F2[z^+/-1,tau^+/-1]/(tau^2+1).
```

Only finite Laurent support in z is allowed. Physical weight counts
the coefficients in the basis `z^i tau^j`, `i in Z`, `j in {0,1}`,
separately in each qubit block.

**Theorem 2.1.** Every nonzero pair in `ker_C(beta,alpha)` has
weight at least six. The pair P in (1.2) attains six and is not
in `im_C(alpha,beta)^T`. Thus the cylinder's nontrivial X minimum
is exactly six as well.

**Proof of the lower bound.** Write `eta=1+tau`, so `eta^2=0`, and
use the following univariate Laurent polynomials:

```
p=z^2+z+1,       q=z^3+z+1,
a=z+z^-4,       b=1+z^-1,       c=z^-2,
alpha=a_0+eta a,       a_0=1+a=z^-4 p q,
beta =b_0+eta c,       b_0=b+c=z^-2 p.              (2.1)
```

The finite factor identities needed here are

```
p q=z^5+z^4+1,       c+ab=z^-5 p^3.               (2.2)
```

For a cycle write

```
f=f_0+eta f_1,       g=g_0+eta g_1.
```

Its constant equation `b_0 f_0+a_0 g_0=0` is over the integral
domain `F2[z^+/-1]`. Cancelling p and monomial units gives, for
some finite h,

```
f_0=z^-4 q h,       g_0=z^-2 h.                    (2.3)
```

For every nonzero finite h,

```
|h|+|qh|>=4,
with equality if and only if h is a monomial.      (2.4)
```

Here is the complete small-weight argument for (2.4). If `|h|=1`,
the sum is four. If `|h|=2`, the two translates of the support
`{0,1,3}` of q meet at most once: its six ordered nonzero
differences are `-3,-2,-1,1,2,3`, all distinct. Hence `|qh|>=4`
and the sum is at least six. If `|h|=3`, the product has odd
weight and cannot be a monomial, because the nonunit q cannot
divide a Laurent monomial; thus its weight is at least three.
If `|h|=4`, the product is nonzero and even, hence has weight
at least two. Finally `|h|>=5` makes the sum at least six.

Folding tau to one cannot increase physical weight. If `h!=0`,
(2.3)-(2.4) therefore give a lower bound of four. For a physical
weight-four cycle, equality must hold throughout, and a common
z-translation normalizes h to one. No site can occur in both tau
layers, since those two sites would cancel under folding. Equivalently,

```
supp(f_1) subset supp(f_0)={-1,-3,-4},
g_1 is either 0 or z^-2.                           (2.5)
```

The eta equation and (2.2) give, with this normalization,

```
z^-2 f_1+z^-4 q g_1=z^-5 p.
```

The two possibilities in (2.5) force respectively

```
g_1=0:       f_1=z^-3 p=z^-1+z^-2+z^-3;
g_1=z^-2:    f_1=z^-3 p+z^-4 q=z^-2+z^-4.
```

Both contain `z^-2`, forbidden by (2.5). So weight four is
impossible. If instead `h=0`, the eta equation has the same
form as (2.3):

```
f=eta z^-4 q j,       g=eta z^-2 j.
```

For nonzero j its physical weight is exactly
`2(|qj|+|j|)>=8`. This deals also with the kernel of folding.
Finally every cycle has even physical weight, by augmenting
`beta f+alpha g=0` at `z=tau=1`. Thus no nonzero weight below
six remains possible.

**Attainment and nontriviality.** Equation (1.4) makes P a
cylinder cycle, and its six sites are distinct. Suppose it were
a cylinder boundary. Choosing Laurent lifts would give

```
P=h(alpha,beta)+(1+tau^2)(u,v).
```

Taking the plane syndrome and cancelling the nonzero factor
`1+tau^2` in the Laurent domain yields
`E=beta u+alpha v`, contrary to (1.5). This proves both claims. QED.

## 3. An all-integer linear-cost original-square family

For every integer `r>=1`, set `N=6r` and

```
S_N=sum_(j=0)^(3r-1) tau^(2j),
V_r=S_N(P_0,P_1) in T_N^2.                        (3.1)
```

**Theorem 3.1.** V_r is a nontrivial original-square X-cycle,
of exact weight `18r`. Its displayed period flux at
`(x,y)=(omega,omega^2)` is `(1,1)`.

**Proof.** The plane identity telescopes to

```
partial_1 V_r=(1+tau^N)E
             =(1+x^N)y^N E+(1+y^N)E.             (3.2)
```

Thus it is a cycle. Each of the six sites in (1.2) represents
a different `(z exponent mod N, tau parity, block)` orbit.
The orbit length is `N/2=3r`. They remain distinct even when
`N=6`, so the weight is exactly `6*3r=18r`.

For completeness, the period-flux nontriviality test used in (3.2)
does not require the full parity/Tor theorem. Write `F=x^N+1`,
`G=y^N+1`. If a plane lift of a torus boundary has syndrome
`F C+G D`, it also has a representation with coefficients of F
and G each lying in `(alpha,beta)`. The two representations differ
by a syzygy of the coprime Laurent polynomials F,G, so their
coefficient differences are `(G H,F H)` for some H. Evaluation
at any common zero of `F,G,alpha,beta` therefore makes both C
and D zero. In (3.2) the coefficients are `y^N E,E`, whose
values are `(1,1)` by (1.5), a contradiction. QED.

This is a compact-cylinder, diagonal-winding family, not a mixed
family. It comes from the cylinder with period `(N,N)` through
the exact syndrome `(1+(xy)^N)E`. The equal nonzero fluxes are
consistent with, but do not require assuming, the parent's parity
classification. Its nontriviality is unconditional.

## 4. All diagonal-period-two cycles: the finite-torus obstruction

Let `I_N` be the subspace of physical X-vectors on the original
square satisfying `tau^2 v=v`. Introduce the small-period quotient

```
C_N=F2[z^+/-1,tau^+/-1]/(z^N+1,tau^2+1).
```

Multiplication by S_N gives a linear bijection

```
iota_N:C_N^2 -> I_N,
|iota_N(w)|=(N/2)|w|=3r|w|.                        (4.1)
```

This statement is about coefficient orbits, and is valid even
when `3r` is even. It is not division by the cover degree or
an assertion that multiplication by S_N is injective on all of
T_N. The map intertwines the syndrome maps and is injective
also on syndrome vectors. Consequently a vector in I_N is a
cycle exactly when its unique quotient motif is a C_N-cycle.

### 4.1 The discount is impossible for every r>=2

**Proposition 4.1.** For `N=6r>=12`, C_N has no nonzero
cycle of weight below four. Thus every nonzero invariant
original-square cycle has weight at least `12r`.

**Proof.** Odd weight is excluded by augmentation. For weight
two, two syndrome columns would have to coincide. In the two
tau layers their supports are

```
alpha: tau=0: {0};       tau=1: {1,-4};
beta:  tau=0: {0,-1};    tau=1: {-2}.              (4.2)
```

All these within-layer sites are distinct for `N>=6`. For
two columns of the same type, the layer cardinalities force
the tau displacement to be zero, and the singleton layer
forces the z displacement to be zero. They would be the same
qubit, not a weight-two vector.

For different types, equality forces a tau displacement of
one. Matching the singleton layers then uniquely fixes the
translate of alpha to `z^-2 tau alpha`. Its double layer is
`{-1,-6}`, whereas beta's is `{0,-1}`. They agree if and only
if `N|6`, impossible for `N>=12`. This rules out weight two.
Apply (4.1). QED.

In particular `12r-6 < 12r` lies strictly below this entire
invariant sector for every `r>=2`, without any support-width
assumption or conditional metric floor.

### 4.2 The exact invariant-sector minimum is 18r for all r>=4

**Theorem 4.2.** For every integer `r>=4`, the minimum nonzero
cycle weight in I_N is `18r`. The minimum nontrivial cycle
weight in I_N is also `18r`.

**Proof.** By parity it suffices to exclude quotient weights
two and four. Here is a support-lifting argument excluding both
uniformly once `N>20`.

Suppose a C_N-cycle uses at most four qubits. Make an interaction
graph on those qubits, joining two whenever their syndrome
columns have a common check. Each connected component is itself
a cycle: no check can be incident to columns in different
components, and each check has even total incidence.

Take a nonempty component with s qubits, `s<=4`, and choose a
spanning tree. Each tree edge can be lifted to an actual shared
check on the infinite cylinder C with a z displacement of
absolute value at most five. This follows from (4.2), or directly
from the z-offset ranges `[-4,1]` and `[-2,0]` of alpha and beta.
Lift the qubit anchors successively along the tree. Their total
z span is at most `5(s-1)`; the syndrome of the lifted pair then
has z span at most `5(s-1)+5=5s<=20`.

For `N>20`, reduction modulo `z^N+1` is injective on a Laurent
support contained in an interval of span at most twenty. The
lifted syndrome reduces to the zero syndrome of the chosen
component, so it is already zero on C. The lift has s distinct
qubits. This contradicts Theorem 2.1, which requires at least
six for every nonzero finite-support cylinder cycle.

Thus all nonzero C_N-cycles have weight at least six for
`N=6r>=24`. Equation (4.1) gives the lower bound `18r`.
The nontrivial family V_r attains it by Theorem 3.1. QED.

The threshold is a sufficient analytic threshold, not a claim
that a weight-four motif exists at N=12 or N=18. Their exact
invariant-sector minima are not needed and were not computed.

Inversion of both lattice variables followed by block swap
takes an X-cycle `(f,g)` to the Z-cycle `(bar g,bar f)`, takes
boundaries to boundaries, preserves weight, and preserves
invariance under the subgroup generated by `(2,2)`. Therefore
Proposition 4.1 and Theorem 4.2 also hold for Z-cycles.

## 5. Precise obstruction and residual

The finite-transporter ansatz under consideration is: choose an
arbitrary motif on `Z_N x Z_2` in the `(z,tau)` coordinates and
repeat it along the original diagonal by `(2,2)`. This includes
unbounded transverse motif width, multiple strands, both blocks,
and arbitrary boundary dressing that leaves the final vector
diagonal-period-two invariant. It is larger than repeating the
particular six-site motif P. Theorem 4.2 resolves its minimum
exactly for all `r>=4`; Proposition 4.1 already excludes the IBM
target for every `r>=2`.

A correction that breaks the symmetry is not covered by the
sector theorem. Nevertheless, if a vector of weight at most
`12r-6` is obtained from V_r by adding a physical difference
Delta, the triangle inequality forces

```
|Delta| >= |V_r|-(12r-6)=6r+6.                    (5.1)
```

This bound is on the changed physical qubits, not the number
of stabilizer generators. Thus a fixed-size cap cannot repair
this seam family for an unbounded sequence of r. No claim is
made against a genuinely extensive, symmetry-breaking correction.

**What remains unresolved:** an original-square family at
`12r-6`, even on an infinite arithmetic subsequence. This note
does not leave an untested cancellation lemma inside the fixed
period-two ansatz: that ansatz is excluded. A successful further
construction must at least break this invariance for every
`r>=2`; a weight argument for other seam directions, longer
fixed motif periods, or nonperiodic junction networks is not
supplied. Nor is any unrestricted distance lower bound supplied.

## 6. Diagnostic scope

Run the standalone checker from the repository root with

```
PYTHONDONTWRITEBYTECODE=1 uv run --offline --no-sync --no-cache \
  --project experiments/bb_lab python \
  experiments/bb_lab/scripts/a44_s3_codex_seams.py
```

It asserts (1.4), the ordinary-point evaluations, (2.2), both
explicit failed weight-four assignments in (2.5), and the six
distinct ordered differences of q. It then checks the one
specified original `(24,24)` witness, with weight 72 and the
claimed translation invariance. It does not computationally
certify a minimum; the lower bounds are the proofs above.

The completed diagnostic reported 0.0001 seconds internally,
16.2 MiB peak RSS, and exited successfully. An earlier invocation
stopped before any diagnostic when an unsupported RSS-limit
setting raised `ValueError`; that setting was removed. CPU and
wall alarms are 50 and 55 seconds, and the checker asserts the
768 MiB peak-RSS budget. Both invocations together used well
below one second of diagnostic process wall time. No user or
sibling files, git state, old processes, memory, or Lean files
were touched.
