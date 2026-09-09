# A42 S6 restart — exact residual, two structural identities, and a partial bound

2026-09-08. Codex research pass on `claude/a40-tour-de-gross-column`.
This is a separate note; it does not amend the A42 ledger or assemble a
distance theorem for the IBM family. All chains below have finite support
in the unbounded cylinder coordinate.

## 1. Outcome and scope

The cylinder floor `d_cyl(p) = 2p` is **not proved at all tower heights**.
Neither the IBM distance law nor an exact growing-distance modification
is proved here. The new analytic consequences are:

1. On every binary tower `q = 2^a`, the binary boundary code has an exact
   graph presentation with transfer polynomial `C_q`. Its impulse has
   weight `3^(a+1) + 1`, and `C_q^q` collapses to three monomials. This is
   a finite symbolic identity, not an extrapolation from impulses.
2. Two compressed parity checks, combined in the **16-state graph
   suggested by the parent**, prove the exact all-`q` coset identity
   `d(1_S, B_q) = q` for an explicit nontrivial omega support
   `|S| = 3q`, whenever `3` does not divide `q`. This attains
   the conjectured minimum support size; its global minimality at
   arbitrary tower height is not assumed or proved here.
3. Consequently, for that omega support and **every binary boundary**,
   with no span, size, or periodic-pullback restriction on the boundary,
   the physical mixed representative satisfies

       |v| >= 4q.

   If all occupied omega fibres are singletons, the exact minimum is

       min |v| = 6q.

   For example, at `q = 8`, these are `32` and `48`, respectively;
   the desired unrestricted cylinder floor is `48`.
   These statements cover an unbounded binary tower of **restricted
   omega supports with arbitrary binary mixing**, not all classes or
   all omega cycles.
4. The omega equation implies an exact quadratic equation for its
   support indicator. Thus its binary syndrome is constrained by
   phase correlations, not just support size.

The audit also isolates the missing mirror implication and a distinction
between globally omega-pure chains and chains pure only on the three-sheet
factor. The latter distinction matters when the odd part grows.
Section 4 retains the initial coarse inequalities as intermediate
results; §8 contains the stronger finite certificate and the current
residual. In particular the canonical coset lemma (8.1) is now **closed**.

Established background, **not discoveries of this pass**: A42 §1.3 gives
`k = 26` at `b = 0, 127 | r`, and `k = 12` elsewhere in the IBM family,
including every `b = 1` member. This does not itself refute the distance
law. The parent's checked IBM reference is arXiv:2506.03094v1, §2.7,
item 3; the cited Liang et al. identifier is **2503.03827**, not the
2503.04699 identifier printed in A42's literature paragraph.

## 2. Conventions and the exact HM statement

The original Laurent pair is

    A = 1 + y + x^3 y^-1,
    B = 1 + x + x^-1 y^-3.

Use the equivalent, monomial-cleared Koszul convention

    P = y A = x^3 + y + y^2,
    Q = x y^3 B = 1 + (x + x^2)y^3,
    D(f,g) = P f + Q g,       boundary(h) = (Q h, P h).

Multiplying either row entry by a monomial corresponds to translating
one chain block and multiplying the degree-two coordinate by a unit;
these are chain isomorphisms preserving Hamming weight. The usual BB
X-cycle convention is obtained by also applying the antipode and a block
swap. These operations preserve compactness and cylinder period.

For `p = 3q`, `3 ∤ q`, choose product coordinates
`Z_(3q) = Z_3 × Z_q`, and write `zeta^2 + zeta + 1 = 0`. Fourier transform
only the `Z_3` coordinate. Then `y` becomes `t` in the binary component
and `zeta t` in the quaternary component. Put

    R_q = F_2[t]/(t^q + 1),   K_q = F_4[t]/(t^q + 1),
    P_0 = x^3 + t + t^2,     Q_0 = 1 + (x+x^2)t^3,
    P_w = x^3 + zeta t + zeta^2 t^2,   Q_w = Q_0.

All these rings also have the Laurent variable `x`. Define

    B_q = { (Q_0 h, P_0 h) : h in R_q[x^±] }.

This is the **boundary code**, whether or not its homology vanishes.
If `3 ∤ q` and `127 ∤ q`, the established finite variety calculation
implies all binary cycles are boundaries. Merely assuming `3 ∤ q`
does not justify that conclusion: an order-127 cylinder component
survives when `127 | q`.

Let `mu = (f,g)` be a nontrivial cycle for `(P_w,Q_w)`, with coefficient
support `S`. Let `s` be a binary boundary. In each three-sheet fibre,
the binary fold and quaternary value are `(s_c,mu_c)`. Their weights are
`0,3,2,1` for `(0,0),(1,0),(0,nonzero),(1,nonzero)`, respectively.
Consequently the inverse Fourier lift satisfies

    |v| = 2|S| + 3|s| - 4|S intersect s|.                 (2.1)

It is a nontrivial physical cycle: a boundary would project to a
quaternary boundary. Let `epsilon = |S| - 3q`. For this `mu`, the following
two statements are equivalent:

    4|S intersect s| - 3|s| <= 2 epsilon    for all s in B_q;

    3 d(1_T,B_q) >= |T| - 2 epsilon        for all T subset S.   (HM)

Here distance minimizes over finite-support chains on the cylinder, not
over an arbitrarily selected embedding torus. Its integer minimum exists
because the set of candidate weights is a nonempty subset of the natural
numbers.

Proof of equivalence. In the reverse direction use `T = S intersect s`;
`s outside S` lies in `1_T + B_q`. In the forward direction, if
`e = 1_T + s`, let `u = |T intersect e|` and
`z = |(S minus T) intersect e|`. The first inequality gives

    3|e outside S| >= |T| - u + z - 2 epsilon.

Adding `3|e intersect S| = 3u + 3z` proves
`3|e| >= |T| - 2 epsilon`. Minimize over `e`. This proof audits the
unrestricted coset-leader reformulation in §2.17.3; it does not prove HM.

If the support floor `|S| >= 3q` and HM hold for all the relevant
quaternary cycles, (2.1) gives the full cylinder floor. Conversely that
floor gives both inequalities. Including `T = empty` in HM already
forces `epsilon >= 0`.

## 3. Binary tower: an explicit graph and its fractal impulse

**Proposition 3.1.** Let `q = 2^a`, `a >= 0`. In `R_q[x^±]`, define
`u = t+t^2`. Then

    P_0^-1 = sum_{j=0}^{q-1} x^(-3(j+1)) u^j
           = x^(-3q) P_0^(q-1),
    C_q = Q_0 P_0^-1,
    B_q = { (C_q z,z) : z in R_q[x^±] }.                 (3.1)

In particular every binary cycle is a boundary, uniquely parametrized
by its second block. Moreover

    |P_0^-1| = 3^a,
    |C_q| = 3^(a+1),
    C_q^q = x^(-3q) + x^(-2q) + x^(-q).                (3.2)

Proof. Frobenius and `t^q = 1` give `u^q = 0` and `P_0^q = x^(3q)`.
The geometric sum in (3.1) is therefore an actual finite inverse. A
cycle satisfies `f = C_q g`; choosing `h = P_0^-1 g` proves the graph
description and uniqueness.

For `j < q`, the binomial expansion
`u^j = t^j (1+t)^j` has `2^popcount(j)` terms: each binary digit of `j`
independently chooses the `1` or `t^(2^digit)` term. Their exponent
differences are at most `j < q`, so reduction modulo `t^q+1` causes no
collisions. Distinct `j` lie in different `x` columns. Summing gives
`sum_{j<2^a} 2^popcount(j) = 3^a`.

The three terms of `Q_0` place each summand into three distinct residue
classes of the `x` exponent modulo three. Within each class, different
`j` again have different columns. Hence multiplication triples its
weight. Finally `Q_0^q = 1+x^q+x^(2q)` proves the last identity. ∎

The boundary with second block one coordinate has weight `3^(a+1)+1`;
its first block occupies every column from `-3q` through `-1`. Thus a
single-coordinate second-block prescription requires a unique first
block with a tail whose width grows linearly with `q`. The boundary
code still has its usual local generators `(Q_0 h,P_0 h)`; this statement
concerns the inverse graph, not a lower bound for all nonzero boundaries.

**Exact negative control.** Take `T = (C_q,0)`. Its binary syndrome is
`Q_0`, of weight three, and

    |T| = 3^(a+1),      d(T,B_q) = 1.

The upper bound follows from the codeword `(C_q,1)`; the lower bound
follows from its nonzero syndrome. This is an infinite symbolic family
against arbitrary-support expansion arguments. It is **not** an HM
counterexample: no suitable nontrivial omega support containing this
`T` with small enough excess has been established.

The collapse in (3.2) also prevents treating repeated applications of
`C_q` as a positive, noncancelling expansion. Neither assertion in this
section depends on a distance certificate.

## 4. Two compressed checks and a new all-height mixed bound

Set `Phi(t) = 1+t+t^2`. For any binary chain `e = (e_f,e_g)`, define

    J_1(e) = Phi(t) e_f(1,t) + e_g(1,t),
    J_w(e) = e_f(zeta,t) + (1+t)e_g(zeta,t).             (4.1)

The first takes values in `R_q`, the second in `K_q`; weights count
nonzero coefficients in the `t` basis.

**Proposition 4.1.** If `3 ∤ q`, both checks annihilate `B_q`, and

    |J_1(e)| <= 3|e_f| + |e_g|,
    |J_w(e)| <= |e_f| + 2|e_g|.                        (4.2)

Consequently, for every binary chain `T`,

    d(T,B_q) >= ceil((|J_1(T)| + 2|J_w(T)|)/5).         (4.3)

Proof. At `x = 1`, the cycle row becomes `(Phi,1)`. At `x = zeta`, it
becomes `(Phi,(1+t)Phi)`, since `zeta^3=1`, `zeta+zeta^2=1`, and
`1+t^3=(1+t)Phi`. The roots of `Phi` have order three, so `Phi` is a
unit modulo `t^q+1` when `3 ∤ q`. Thus both checks vanish on cycles,
in particular on boundaries. Evaluation at a nonzero field element
does not increase coefficient support; multiplication by a polynomial
with `k` terms increases weight by at most `k`. This proves (4.2).
For a coset representative `e = T+s`, its checks equal those of `T`.
Adding the first inequality and twice the second gives
`|J_1(T)|+2|J_w(T)| <= 5|e|`. ∎

These are consequences of the actual binary cycle equation, and use
more than the unprocessed syndrome count in A42 §3.14. They do not
assume a pure floor or HM.

Let `N_q = 1+t+...+t^(q-1)` and choose the quaternary cycle

    mu_* = (N_q,(1+x)N_q).                             (4.4)

**Proposition 4.2.** For every positive `q` with `3 ∤ q`, this is a
nontrivial quaternary cylinder cycle, with support

    S_* = {all q fibres at first-block column 0}
          union {all q fibres at second-block columns 0 and 1}.

It has `|S_*| = 3q`, and

    D_0(1_(S_*)) = N_q,
    J_1(1_(S_*)) = J_w(1_(S_*)) = N_q,
    ceil(3q/5) <= d(1_(S_*),B_q) <= q.                (4.5)

Proof. Since `(1+t)N_q=0`, every polynomial in `t` acts on `N_q` by
evaluation at one. Thus on these coefficients

    P_w = x^3+1,    Q_w = 1+x+x^2,
    (x^3+1) + (1+x+x^2)(1+x) = 0.

This proves the cycle equation. If it were a boundary, some compact
`h` would satisfy `Q_w h = N_q`. For nonzero `h` whose extreme columns
are `c_min,c_max`, the constant and leading coefficients of `Q_w`
are units, so the extremes of `Q_w h` are exactly
`c_min,c_max+2`. They cannot both be column zero. Hence no such `h`
exists. This argument works despite nilpotents and does not use a
Tor classification or a finite window.

The same evaluation at `t=1` gives
`P_0 N_q + Q_0(1+x)N_q = N_q`. Direct substitution into (4.1) gives
the two checks in (4.5). Apply (4.3) for the lower bound. For the upper
bound take

    h_* = (x^-3+x^-2)N_q,
    (Q_0 h_*,P_0 h_*) = ((1+x^-3)N_q,(1+x)N_q).

Its difference from `1_(S_*)` is `(x^-3 N_q,0)`, of weight `q`. ∎

In particular the failed raw-syndrome sufficient condition
`|D_0(1_S)| >= |S|-2epsilon` is false on this explicitly nontrivial
family at **every** such `q`: its two sides are `q` and `3q`. This
upgrades the ledger's small examples to an exact infinite counterfamily
to that sufficient condition, without refuting HM itself.

**Corollary 4.3 (new partial floor).** Fix `mu_*` as in (4.4), or any
nonzero scalar multiple, and let `s` be any binary boundary of arbitrary
finite support. In its physical inverse Fourier lift let

    k = |S_* minus s|,       e = |s minus S_*|.

Then

    |v| = 3q + k + 3e,
    k+e >= ceil(3q/5),
    |v| >= 3q + ceil(3q/5).                            (4.6)

If `S_* subset s`, then

    |v| >= 3q + 3 ceil(3q/5).                          (4.7)

Proof. The error `1_(S_*)+s` has weight `k+e`, so (4.5) applies.
Equation (2.1) gives the displayed exact physical weight. For (4.7)
set `k=0`. Every physical cycle has even weight: augmentation of the
cycle equation evaluates both three-term checks to one and forces
the sum of the two block parities to vanish. Thus even rounding is
valid. ∎

This is stronger than the fibre halving bound `3q`, and it holds for
arbitrarily complicated binary mixing, including at unbounded
`a` when `q=2^a`. It does **not** bound all nontrivial `mu` by (4.6).
The finite checks in the accompanying script are sanity checks; the
proofs above establish the quantifiers.

## 5. An exact quadratic constraint on omega supports

For a quaternary coefficient `z`, let `n(z)=z^3` in `{0,1}` and
`Tr(z)=z+z^2`. At a check location `(c,j)`, the six summands of
`P_w f + Q_w g` are

    z_1 = f[c-3,j],       z_2 = zeta f[c,j-1],
    z_3 = zeta^2 f[c,j-2], z_4 = g[c,j],
    z_5 = g[c-1,j-3],     z_6 = g[c-2,j-3].

**Proposition 5.1.** For every omega cycle,

    (P_0 n(f)+Q_0 n(g))[c,j]
      = sum_{i<j'} Tr(z_i z_(j')^2).                   (5.1)

Proof. In `F_4`, expand

    n(sum_i z_i) = (sum_i z_i)(sum_i z_i^2)
                = sum_i n(z_i) + sum_{i<j'} Tr(z_i z_(j')^2).

The left side vanishes for a cycle. Twisting a coefficient by a
nonzero element of `F_4` does not change its norm, so the sum of the
six norms is exactly the binary syndrome on the left of (5.1). ∎

This remains valid when periodicity identifies some source locations;
the six summands retain their multiplicities before cancellation.
It supplies a checkable constraint on allowable `S` and their phases.
Its right side is an `F_2` sum and can cancel. Counting its terms as
positive mass is invalid. For `mu_*`, it is exactly the small source
`N_q` from (4.5). A norm identity by itself therefore does not repair
the failed syndrome proof.

## 6. Which cylinder periods are actually covered?

There are two different questions: which old machine certificates
remain applicable, and which infinite statements follow analytically.
This pass inspected the relevant notes, implementation conventions,
and recorded certificate metadata; it did **not** reexecute their large
enumerations or independently reprove their entire search completeness.
The retained **banked certificate tier** is:

| Object | Covered periods | Qualification |
|---|---|---|
| Original AB full cylinder, equality `2p` | `3,6,9,12` | Small atlas; period-9 register at depth zero; repaired period-12 joint-jet/corridor certificate |
| Rotated `(B,bar A)` full cylinder | `3,6` | Independently run two-lane atlas at periods at most 8; no justified transfer of the later AB-only certificates |
| Either lane, no compact homology | `3 ∤ p` and `127 ∤ p` | Established finite variety plus regularity/Bezout argument; floor is infinity, not `2p` |
| AB globally omega-pure cylinder | `p=3 m 2^a`, `m` odd, `a=0,1,2` | The pure base certificates at `3,6,12` and the exact odd-transfer argument below |
| Rotated globally omega-pure cylinder | `p=3 m 2^a`, `m` odd, `a=0,1` | Pure base minima at `3,6`; no proof of the rotated `a=2` base is supplied here |

The strict pure statements specify the order-three primary component;
they do not include an additional order-127 component. In particular:
full AB period `18` remains open; full AB period `24` still needs the
pure tower bound and binary HM. Period `9` already has a full AB
certificate, so the `HM_9` instance with `q'=1` is covered; period `18`
corresponds to `q'=2` and is not.

**Odd transfer argument.** Put `P=3·2^a` and let `m` be odd. In period
`mP`, the norm `E=1+y^P+...+y^((m-1)P)` is an idempotent, is one modulo
`y^P+1`, and vanishes on its coprime complementary factor. Transfer
from period `P` multiplies a chain by `E`, preserves its primary
content, and makes `m` disjoint identical blocks. Restricted to the
order-three primary component it is an isomorphism of the cycle and
boundary complexes, and multiplies weight by `m`. Thus strict
omega-pure minima scale exactly by `m`. This proves the infinite
period extension in the table, rather than inferring it from the
list `3,6,9,12,15,18,21` of executed pure racers.

**Two meanings of pure.** When `q=2^a`, `K_q` is precisely the
order-three primary factor of the period-`3q` binary ring. When
`q=m 2^a` with odd `m>1` and `3 ∤ m`, `K_q` also contains factors
coming from nontrivial characters of `Z_m`. These factors can have
zero homology while their boundary chains are nonzero. Therefore

    globally omega-pure    implies    binary fold zero,

but the converse is false in this case. The exact three-fibre HM
decomposition still works if `mu` ranges over **all** nontrivial
cycles of `K_q`, but then its prerequisite `|S|>=3q` includes those
additional quaternary boundary directions. The strict pure scaling
theorem does not discharge that stronger prerequisite. One must
either retain these directions explicitly, or include their hiding
in a global-CRT formulation. This is a scope distinction, not a new
counterexample to the conjectured floor.

**127 qualification.** A42 §2.17.3's identification of cycles with
boundaries needs `127 ∤ q` on the cylinder. A torus with transverse
period not divisible by 127 can kill that homology, but it is a
different boundary code. Even when a member has no order-127 class,
an unrolled chain can contain an order-127 contribution along with its
omega contribution; possible weight cancellation cannot simply be
discarded. This is consistent with the already-established dimension
exception, not a claim that new species refute the distance formula.

**Higher powers of three.** The nine-sheet description is exactly for
`p=9q'`, `3 ∤ q'`, not every period divisible by nine. The latter also
includes `27,81,...`, requiring higher fibres or an explicit iteration.
For `p=9q'` and `127 ∤ q'`, homology lies in the quaternary factor of
`F_2 × F_4 × F_64`; the other components are boundaries. The pure
three-character fibre has weight six. With `epsilon=|S|-3q'`, the
remaining inequality is precisely

    sum_(c in S) (6-w_c) <= sum_(c outside S) w_c + 6 epsilon. (HM_9)

At `q'=2` its pure prerequisite is banked, while this inequality is
open. With extra odd factors in `q'`, the same distinction between
strict omega purity and fibre purity applies.

## 7. The mirror audit and the tower residual

Let `theta'(x,y)=(y^-1,x)` on Laurent polynomials. Direct substitution
gives

    theta'(A)=B,        theta'(B)=bar A.

This proves a weight-preserving chain isomorphism from the original
pair to the **transformed** pair, with the cylinder direction exchanged.
It does not identify the two pairs on the same cylinder direction.

**Explicit negative control.** In the unbarred Koszul convention the
original boundary is `(B,A)`. Its image is `(bar A,B)`, which is a cycle
for row `(B,bar A)`. Applied instead to the original row `(A,B)`, its
syndrome is

    A bar A + B^2 != 0.

The coefficient of `x^2 y^0` is one: `A bar A` has only `x` exponents
`0,3,-3`, whereas `B^2` has a term `x^2`. This remains a counterexample
on every cylinder with `x` unbounded. It refutes the proposed direct
self-intertwiner; it does not prove that every other possible isometry
is impossible or that the two minimum weights differ.

The period-9 and period-12 certificate artifacts literally record
`AB (theta'-dual lane by banked duality)`; the corresponding later
engines do not run a separate BAbar lane. Thus the original AB results
remain banked, but transferring those results to the mirror needs
another theorem or an independent certificate. The independently run
two-lane atlas at periods at most eight does not have this problem.
At period three there is also a direct same-direction isometry:
`B=bar B`, so a global antipode followed by a block swap identifies
the transformed pair with the original. This special argument does
not justify an arbitrary-period identification.

For the binary tower set `d_a=d_cyl(3·2^a)`. The established local
intersection data imply transfer is an isomorphism on compact homology
for `a>=2`, and fold is zero there. Equivalently these are the usual
Tor maps for the local modules of lengths one and two. Regularity uses
the **finite complete zero set**, including the order-127 orbit; it
does not use an incorrect claim that the global plane quotient has
only the four F4 zeros or total dimension six. The order-127 points
do not contribute on this binary tower.

With that homological input, a lifted cycle has the exact form

    v = transfer(g) + boundary(h_tilde),
    [g] != 0,

and its fold is a base boundary. Therefore the metric statement still
needed at a rung is

    |transfer(g)+boundary(h_tilde)| >= 2 d_(a-1)
    for all nontrivial base g and all compact h_tilde.       (7.1)

To identify (7.1) with the target `d_a>=6·2^a`, the base hypothesis
`d_(a-1)=6·2^(a-1)` is necessary. This is known at the first unresolved
rung `a=3`; for higher rungs it belongs to the induction hypothesis.
The statement for the pure floor restricts both the cycle and its
boundary changes to the strict pure component. The full statement
allows all boundary changes.

Proof of the lifting form. If the fold of `v` is `boundary(h)`, choose
a chain lift `h_tilde`. Then `v+boundary(h_tilde)` has zero fold. On
chains the zero-fold space is exactly the image of transfer. Transfer
is injective on chains and commutes with the differential, so its
preimage `g` is a cycle. Homological injectivity makes it nontrivial.
The two sheets of `v` separately need not be base cycles: they have
the same seam syndrome. Consequently a base distance bound cannot be
applied separately to their two weights. This is precisely the metric
gap, even though transfer on homology is already known. ∎

## 8. A 16-state certificate closes the canonical coset lemma

The parent proposed combining the compressed checks before taking their
separate Hamming bounds. This succeeds, with a potential taking only
two values. The entire graph has 16 states and 128 edges; no growing
search is used.

**Theorem 8.1 (canonical syndrome, exact).** Let `q>=1` and `3 ∤ q`.
Every compact binary error `e=(e_f,e_g)` satisfying

    P_0 e_f + Q_0 e_g = N_q

has `|e_f|+|e_g| >= q`. Consequently `d(1_(S_*),B_q)=q`. On the binary
tower this also proves the original explicit graph minimization

    min_z ( |x^-3 N_q + C_q z| + |z| ) = q.             (8.1)

Proof. Aggregate each error block's `x` coordinates modulo three,
using parity within each residue. At each row `j` write the Fourier
coordinates of the first resulting binary triple as
`(a_j,b_j) in F_2 × F_4`, and those of the second as `(c_j,d_j)`.
Aggregation cannot increase physical Hamming weight. The exact weight
of a triple with Fourier coordinates `(a,b)` is

| `a` | `b` | `w(a,b)` |
|---|---|---|
| 0 | 0 | 0 |
| 1 | 0 | 3 |
| 0 | nonzero | 2 |
| 1 | nonzero | 1 |

Evaluation of the syndrome at `x=1` gives `Phi a+c=N_q`.
Evaluation at `x=zeta` gives `Phi (b+(1+t)d)=N_q`. Since `Phi` is
invertible and `Phi N_q=N_q`, the necessary constraints are exactly

    c_j = 1+a_j+a_(j-1)+a_(j-2),
    b_j = 1+d_j+d_(j-1).                              (8.2)

All subscripts are cyclic modulo `q`. These equations imply that
the aggregated error defines a closed walk in the graph with state

    U=(a_prev,a_prev2,d_prev) in F_2 × F_2 × F_4.

An edge chooses `(a,d) in F_2 × F_4`, goes to `V=(a,a_prev,d)`, and
has nonnegative integer cost

    c(U,V)=w(a,1+d+d_prev)+w(1+a+a_prev+a_prev2,d).      (8.3)

Its total cost is precisely the aggregated error weight. No claim
that every graph walk lifts to an actual syndrome solution is needed:
the direction from actual errors to walks proves a lower bound.

Here is the **complete Bellman potential** `F`; columns are `d_prev`.

| `(a_prev,a_prev2)` | 0 | 1 | zeta | zeta^2 |
|---|---|---|---|---|
| `(0,0)` | 0 | 0 | 1 | 1 |
| `(0,1)` | 0 | 1 | 1 | 1 |
| `(1,0)` | 0 | 1 | 1 | 1 |
| `(1,1)` | 1 | 1 | 1 | 1 |

For every one of the 128 edges,

    c(U,V)-1 + F(U)-F(V) >= 0.                        (8.4)

For a short verification independent of the program, let
`L={(0,0,0),(0,0,1),(0,1,0),(1,0,0)}`, the zero-potential states.
The only cost-zero edges are

    (0,1,1) -> (0,0,0),    (1,0,1) -> (0,1,0),

both with choice `(a,d)=(0,0)`, so they go from potential one to zero.
The outgoing costs from states in `L` are the following eight-entry
rows, in choice order
`(0,0),(0,1),(0,zeta),(0,zeta^2),(1,0),(1,1),(1,zeta),(1,zeta^2)`:

| source | eight outgoing costs |
|---|---|
| `(0,0,0)` | `5,1,3,3,1,5,3,3` |
| `(0,0,1)` | `3,3,3,3,3,3,3,3` |
| `(0,1,0)` | `2,2,4,4,4,4,2,2` |
| `(1,0,0)` | `2,2,4,4,4,4,2,2` |

The two cost-one edges in this table land in `(0,0,1)` and `(1,0,0)`,
both in `L`. Hence every edge from `L` to its complement costs at
least two. These observations prove (8.4): all other edges have
cost at least one and do not increase potential. Summing (8.4) on
the cyclic walk cancels the potential and gives total cost at least
`q`. Physical error weight is no smaller. The explicit weight-`q`
leader in Proposition 4.2 proves equality for the coset. The self-loop
at `(1,1,0)` with choice `(1,0)` has cost one, so the graph's minimum
cycle mean is exactly one. Finally, on the binary tower `P_0` is a
unit and `ker D_0=B_q`, proving the expression (8.1). ∎

This proves the required infinite statement from a single finite
inequality table. It does not infer anything from finite values of `q`.

**Corollary 8.2 (stronger restricted-support mixed bounds).** For
`mu_*` and arbitrary binary boundary `s`, the notation of Corollary
4.3 now satisfies

    k+e >= q,       |v|=3q+k+3e >= 4q.                 (8.5)

If `S_* subset s` (every occupied omega fibre is a singleton), then
`k=0` and `|v|>=6q`. The boundary produced by `h_*` in Proposition
4.2 contains `S_*` and has exactly `q` extra fibres. It attains `6q`.
Thus this entire no-pair stratum has **exact** minimum `6q` for every
`3 ∤ q`, in particular for an unbounded binary tower. The lower bound
allows all binary boundaries, not just transferred boundaries.

This is an exact infinite statement about a specified set of
representatives. It is **not** an exact-distance theorem for an
infinite family of quantum codes: other omega supports, proper
subsets in HM, and other homology classes remain to be bounded.

### Zero-slack recurrence audit: the proposed all-height uniqueness does not follow

The parent's proposed refinement was to classify the zero-slack edges
of (8.4), force a constant equality walk when `3 ∤ q`, and deduce
uniqueness of the canonical coset leader. Exact inspection of the
**same 16-state graph** finds 17 zero-slack edges and two recurrent
strongly connected components. They are not a separate self-loop
component and a separate three-cycle component.

Write

    A=(0,0,0), B=(1,0,0), C=(0,1,1), D=(1,1,0),
    E=(1,1,zeta), F=(1,1,zeta^2).

The first recurrent component contains all six of `A,B,C,D,E,F`.
Its complete list of internal edges, with original costs, is

| source | target and cost |
|---|---|
| `A` | `B:1` |
| `B` | `C:2`, `E:2`, `F:2` |
| `C` | `A:0` |
| `D` | `C:1`, `D:1` |
| `E` | `D:1` |
| `F` | `D:1` |

It includes the anticipated self-loop and three-cycle, but also the
five-cycles `A-B-E-D-C-A` and `A-B-F-D-C-A`, each of total cost five.
For example,

    A-B-E-D-D-D-D-C-A

is a nonconstant closed walk of length eight and cost eight. This
explicit fixed-size witness already prevents using `3 ∤ q` to force
a constant equality walk at the first unresolved tower height `q=8`.
It is a witness in the compressed graph, **not an actual cylinder
error or a counterexample to uniqueness**.

The other recurrent component is the two-cycle

    (0,0,zeta) -> (0,0,zeta^2) -> (0,0,zeta),

with both costs one. For completeness, the other six zero-slack edges
are transient edges between or outside these recurrent components:

    (0,0,0) -> (0,0,1)       cost 1,
    (0,1,0) -> (1,0,zeta)    cost 2,
    (0,1,0) -> (1,0,zeta^2)  cost 2,
    (1,0,1) -> (0,1,0)       cost 0,
    (1,1,zeta) -> (0,1,zeta^2) cost 1,
    (1,1,zeta^2) -> (0,1,zeta) cost 1.

**One extra component can be excluded analytically for actual errors.**
On the two-cycle the aggregated first block is zero and the second
block has exactly one site per row. If a physical error has weight
`q`, equality in Hamming contraction forces its first block to be
zero as well. Its syndrome equation would then be `Q_0 e_g=N_q`.
This has no nonzero compact solution: `Q_0` has unit constant and
leading coefficients, so a nonzero input with extreme columns
`c_min,c_max` produces extremes `c_min,c_max+2`, which cannot both
be zero. The six-state component, however, still has nonconstant
equality walks; no general obstruction to their actual lifts has
been proved here.

**Limited uniqueness consequence.** The canonical leader is uniquely
`(x^-3 N_q,0)` for `q=1,2,4`. Indeed an equality error yields a
zero-slack closed walk. The two-cycle is excluded as above. In the
six-state component a closed walk avoiding `D` is a repetition of
the three-cycle; a nonconstant closed walk visiting `D` has length
at least five. Thus at these three periods only the self-loop at
`D` remains. It has one first-block site in residue `x=0 mod 3`
per row and no second-block sites. Equality under aggregation
therefore gives `e_g=0`. Finally `P_0` is injective on compact Laurent
chains: a highest nonzero input column produces a highest output
column three steps later, with unchanged nonzero coefficient. Since
`P_0(x^-3 N_q)=N_q`, the first block is uniquely determined. This
small-period argument does not establish all-height uniqueness.

**The proposed small-hole implication is valid conditionally.** Let
`E_err=|1_(S_*)+s|=k+e` for a binary boundary `s`. If the canonical
leader is unique at the given period and `k>0`, this error differs
from the leader: the leader has no coordinates inside `S_*`.
Augmentation of its syndrome equation gives `E_err=q mod 2`, so
uniqueness implies `E_err>=q+2`. Therefore

    |v|=3q+3 E_err-2k >= 6q+6-2k.

This proves the `6q` floor for `1<=k<=3` **under that uniqueness
hypothesis**. It applies unconditionally to the three small periods
just proved, but the SCC audit does not establish it at arbitrary
tower height. No all-height small-hole or uniqueness claim is banked.

The full fixed-support HM target is equivalently the quantitative
stability estimate

    k <= (3/2)(E_err-q).

That stronger estimate remains open. The existing all-height results
`d(1_(S_*),B_q)=q`, `|v|>=4q`, and the exact singleton-stratum minimum
`6q` are unaffected by the additional equality cycles.

### The smallest remaining refinement

Even for `mu_*`, the full desired floor requires

    k+3e >= 3q                                       (8.6)

for arbitrary `s`, whereas the new certificate yields only `k+e>=q`.
The case `k=0` is now closed. In the HM language, `T=S_*` is closed
exactly, but proper `T subset S_*` still require

    3 d(1_T,B_q) >= |T|.

A remaining sufficient local target is therefore this weighted
error inequality for `S_*`; general HM additionally ranges over all
nontrivial omega cycles and their excess. The graph in Theorem 8.1
forgets which physical coordinates of an error lie inside `S_*`, so
its unit-cost certificate does not establish (8.6). The pure tower,
general HM, the mirror, and the toroidal comparison remain separate.

**Explicit proper-subset obstruction to the two-check method.** The
parent supplied the subset

    T = (0,(1+x)N_q) subset S_*,       |T|=2q.

It satisfies `J_1(T)=J_w(T)=0`: at `x=1` the second block cancels;
at `x=zeta`, its contribution is `(1+t)(1+zeta)N_q=0`. Nevertheless,
on the cylinder with `x` unbounded,

    D_0(T) = Q_0(1+x)N_q = (1+x^3)N_q != 0.

Here `t^3 N_q=N_q` and `(1+x+x^2)(1+x)=1+x^3`. Thus `T` is not a
boundary, although both compressed checks vanish. Minimization using
only these two check values admits the zero error and returns zero;
it cannot prove even a positive coset-distance lower bound for this
`T`, much less the HM requirement `d(1_T,B_q)>=ceil(2q/3)`. This is
an analytic obstruction to sufficiency of the two checks for all
subsets, not a counterexample to HM. No additional computation is
needed for these identities.

For clarity, the current general bound for the fixed support `S_*`
and arbitrary binary mixing is **`4q`**, already even; the coarse
`3q+ceil(3q/5)` bound in §4 is superseded. The singleton-fibre
stratum has exact minimum **`6q`**, attained by `h_*`, whose boundary
contains `S_*` with exactly `q` extra fibres. The full cylinder floor
remains open.

### Restricted IBM parameters

For `b=1`, restricting to `r=4^j`, `j>=0`, removes two arithmetic
complications at once: neither `r` nor `r+1` is divisible by three or
127. Indeed `4^j=1 mod 3`; and `2` has order seven modulo 127, so its
powers cannot equal `-1`. The two cylinder periods are then

    6r = 3·2^(2j+1),       6(r+1) = 6(4^j+1).

The first has an unbounded binary tower; the second has binary depth
one and an unbounded odd part. This removes higher-three fibres and
order-127 components, **but does not close the residuals**: the first
still needs the pure tower floor and HM; the second needs a justified
mirror and control of odd-part boundary mixing. Its strictly
omega-pure part alone scales from period six. The toroidal comparison
or class estimates must also be supplied separately. Thus this is
an arithmetically cleaner candidate, not a solved infinite subfamily.

Finally, A42 **mixed** means that a representative has nonzero binary
fold or other non-omega content. S11 **mixed** means a torus homology
class outside the designated transfer subspaces. The mixed chains in
Corollary 4.3 are compact cylinder representatives and have no automatic
membership in the S11 mixed-class sector. Neither a comparison theorem
nor (4.6) supplies the still-missing IBM metric constants by itself.

### Connectivity and contraction audit after the parent's A43 update

Define footprint adjacency on occupied coordinates by joining two
coordinates whenever a syndrome check contains both. Each footprint
component of a cycle is itself a cycle: a check cannot see occupied
coordinates from two different components. Thus a globally
distance-minimal nontrivial cycle has only one component. Otherwise
some component is nontrivial (a sum of boundaries is a boundary),
and that component alone is lighter.

This proves **nothing of the form class-minimal implies connected**.
A minimum representative in a prescribed class can have multiple
nontrivial components whose classes sum to that class; deleting one
changes the class. Only a trivial component can be removed while
preserving the class. A connected-component search can therefore
establish a common lower bound for all nonzero classes if every
connected nontrivial cycle has that bound, but it cannot classify
every class minimum or every minimum representative in a chosen class.

No connectivity hypothesis is used in (2.1), the graph inverse,
Propositions 4.1–5.1, Theorem 8.1, the odd-transfer identity, or the chain-level
seam decomposition. In particular the boundary in Corollary 4.3 is
arbitrary and can connect or disconnect pieces of the resulting
physical support. Theorem 8.1 and the open (7.1) retain their full
compact chain quantifiers; they are not reduced to connected class
minimizers.

In the inspected S1h, S2 jet, S3 corridor, and S3b pure racer code,
no explicit connectivity or connected-class-minimizer filter was
found. The corridor inequality in A42 §2.9.1 is a prefix/suffix weight
identity valid for disconnected supports as well. This specific
audit does not expose a new connectivity dependency in the retained
global floors. It does not promote the old class-profile tables to
complete class minima: A42 §2.2.6 and §3.8 retain their independent
slot/span/gap limitations. Nor does the corridor argument repair the
already-retracted attempt to split omega supports while allowing a
binary boundary to bridge the pieces.

The parent's verified `(54,6)` weight-48 example, with gaps `(3,1)`
but a cylinder lift, is used only as supplied audit input, not as a
new computation here. It reinforces that a large visible gap is a
sufficient unrolling test, not a characterization of all cylinder
lifts. This note uses actual compact lifts when invoking cylinder
statements and makes no converse inference from absent gaps.

Finally physical Hamming contraction under projection supplies a
lower bound only when the projected class is detected as nontrivial.
Repeated norms in the strict pure sector have exact disjoint-support
weight multiplication; that separate identity cannot be extended to
arbitrary mixed chains. The parent's projection-detection comparison
theorem is therefore compatible with this residual analysis and
does not supply the missing factors in (7.1) or HM.

## 9. Diagnostic and handoff

New file `scripts/a42_s6_codex_identities.py` uses only the Python
standard library. It checks the finite polynomial identities at
`q=1,2,4,8,16,32`, the mirror counterexample, the compressed checks,
the canonical syndrome and correction, and the polarization identity
on four explicitly chosen boundary perturbations per period. It is
not a SAT run, atlas, census, or search for minimum distance.
The follow-up adds only the fixed 16-state/128-edge graph: Bellman
relaxation discovers the potential, and every edge inequality is
then checked against the banked table. This is a finite certificate
for Theorem 8.1, not a growing graph search.
The further equality audit computes reachability only among those
same 16 vertices and verifies all 17 zero-slack edges, the two
recurrent components, and the explicit eight-step graph witness.
No new period sweep or search for lifted physical errors was run.

The successful run took **0.103 seconds** and peaked at **17,694,720
bytes RSS**. Two preliminary invocations failed before any checks
because macOS rejected the requested RSS-limit setting; the final
script uses the 55-second alarm, a CPU limit, and measured peak RSS.
No old jobs were resumed or signalled. No existing notes, scripts,
dirty user files, sibling checkout, or git state were edited.

Parent handoff: retain the banked original-lane periods and their
certificate tier; add (3.1)–(5.1) and the exact canonical coset
certificate (8.1)–(8.5). The remaining obligations are the proper-subset
weighted inequality (8.6), general HM, the pure tower, the mirror,
and torus assembly. The successful initial graph run took 0.100
seconds including all existing tiny checks, with peak RSS 18,038,784
bytes. All diagnostic invocations combined remain below one second
of computation and far below the 180-second aggregate budget.
