# A weight-six BB family with proven unbounded figure of merit

2026-09-13. Mathematical theorem and proof. The construction is explicit;
the exact classical minimum `d_s` is a finite, computable parameter. Its
closed formula is not proved here. This note proves an informative
distance lower bound, an explicit upper witness, and a strictly improving
subsequence. The core construction, actual attained distance, weaker bounds,
and unbounded strictly increasing family are now formalized in Lean; the
precise formal scope is recorded in Section 7. Capstone compilation and axiom
audits and all required repository checks have passed.

The proof was independently reviewed by the parallel construction agent.
No novelty claim is made over Yoshida's fractal-code results or the established
hypergraph-product construction. The step developed here after the earlier
candidate note is to average
open fractal histories over cyclic shifts, expressing every history
coordinate as at most `s` physical code coordinates. There is no assumption
that cutting the periodic code produces an initial zero half.

## 1. The theorem

For each integer `s >= 1`, put

```
q_s = 2^s,
N_s = q_s^2 - 1,
p_s(z) = 1 + z + z^(q_s) in F_2[z].
```

Let `C_s` be the binary cyclic code of length `N_s` consisting of sequences

```
c_(j+q_s) = c_j + c_(j+1),        j in Z_(N_s),
```

and let `d_s` be the minimum Hamming weight of a nonzero element of `C_s`.
Thus `d_s` is defined by an explicit finite minimum, with no unproved
distance assertion included in its definition.

On `Z_(N_s) x Z_(N_s)`, use the BB pair

```
A_s = 1 + x + x^(q_s),
B_s = 1 + y + y^(q_s),
H_X = [A_s | B_s],
H_Z = [B_s^T | A_s^T].
```

**Theorem.** This is a connected BB code with exactly weight-six checks,
and its exact quantum parameters are

```
[[ 2*(4^s-1)^2, 2*4^s, d_s ]].                 (1)
```

The analytic proof in this note gives

```
max(2^s, ceil(3^s/(2*s))) <= d_s <= 3^s - 1.    (2)
```

In particular, writing `R_s = k_s*d_s^2/n_s`,

```
(9/4)^s/(4*s^2) < R_s < (16/9)*(9/4)^s.        (3)
```

The values `R_s` are unbounded. Moreover the explicit subsequence

```
s_t = 16*2^t,        t = 0,1,2,...,
```

has strictly increasing **actual** figure of merit:

```
R_(s_(t+1)) > R_(s_t).                         (4)
```

The exact minimum remains `d_s` in (1). Formula (2) does not claim the
conjectural equality `d_s=3^s-1` for every `s`.

## 2. Algebra and the exact upper witness

In this section fix `s` and abbreviate `q=q_s`, `N=N_s`, `p=p_s`.
Characteristic two gives

```
p(z)^q + p(z) = z^(q^2)+z = z*(z^N+1).
```

The constant term of `p^(q-1)+1` vanishes. Hence the polynomial

```
h(z) = (p(z)^(q-1)+1)/z
```

satisfies `p*h=z^N+1`. In particular, `p` divides `z^N+1`.

The kernel of multiplication by `p` modulo `z^N+1` consists exactly of
the products `h*a`, where `deg a<q`. To prove this, the divisibility
`p*v=(z^N+1)*b=p*h*b` implies `v=h*b` by cancellation in `F_2[z]`;
reduction modulo `p` makes the coefficient `b` unique of degree below
`q`. The `q` shifted generators are independent and have degree below
`N`. This cyclic kernel has dimension `q`. Coordinate reversal changes
its convolution equation into the recurrence defining `C_s`, so the
two codes have identical dimensions and Hamming weights.

Since `q-1` has `s` binary digits all equal to one,

```
p(z)^(q-1) = product_(j=0)^(s-1)
                  (1 + z^(2^j) + z^(q*2^j)).
```

The terms are in bijection with pairs `a,b` in `[0,q)` whose binary
supports are disjoint, the exponent being `a+q*b`. Base-`q` expansion
makes these exponents distinct. There are `3^s` such pairs. Deleting
the constant term and shifting down one gives

```
weight(h) = 3^s-1.
```

After reversal this is a nonzero word of `C_s`, proving the upper
bound in (2). There is no cancellation assumption about products
`h*a` with general `a`.

## 3. An open-history weight inequality

For a nonzero polynomial `a(z)` of degree below `q`, define

```
T_s(a) = sum_(t=0)^(q-1) weight(a(z)*(1+z)^t).
```

Weights here count polynomial coefficients before any cyclic reduction.
All products have degree at most `2q-2`.

**Lemma 1.** `T_s(a) >= 3^s` for every such nonzero `a`.

Here is a self-contained proof. Write `pop(u)` for the number of one
bits of a nonnegative integer `u`. First establish the following
weight-retention statement by induction on `m`:

```
f != 0, deg f < 2^m, ord_(z+1)(f)=u
    implies weight(f) >= 2^pop(u).             (5)
```

The case `m=0` is immediate. For the induction step put `b=2^(m-1)`
and write `f=f_0+z^b*f_1`, with both parts of degree below `b`.
Because `(1+z)^b=1+z^b`, the remainder modulo `(1+z)^b` is
`f_0+f_1`.

If this remainder is nonzero, its order at `z+1` equals `u<b`.
The induction hypothesis and the Hamming triangle inequality give

```
weight(f)=weight(f_0)+weight(f_1)
         >= weight(f_0+f_1) >= 2^pop(u).
```

If it is zero, `f_0=f_1!=0` and
`f=(1+z)^b*f_0`. Therefore `u=b+ord_(z+1)(f_0)` and
`weight(f)=2*weight(f_0)`, which gives (5) again.

Now let `u=ord_(z+1)(a)<q`. Applying (5) with `m=s+1` to each
`a(1+z)^t` gives

```
T_s(a) >= sum_(t=0)^(q-1) 2^pop(u+t).
```

For `0<=i<q`, `2^pop(q+i)=2*2^pop(i)`. Consequently the last sum
is

```
sum_(i=0)^(q-1) 2^pop(i) + sum_(i=0)^(u-1) 2^pop(i)
    >= 3^s.
```

The first sum is `3^s` because each of the `s` bits contributes
`1+2`. This proves the lemma.

This open-history inequality is also the binary case of Yoshida's
[*Information storage capacity of discrete spin systems*](https://arxiv.org/pdf/1111.3275),
Theorem 2. His Pascal-row Lemma 2 is a related form of (5). The proof
above supplies the precise polynomial statement used here; the cyclic
boundary step is handled separately below.

## 4. Each history coordinate uses at most s cyclic coordinates

Let

```
V = F_2[z]/(p).
```

This is a `q`-dimensional vector space with basis
`1,z,...,z^(q-1)`. It is not assumed to be a field. The identities

```
z^q=1+z,       z^N=1
```

hold in `V`, so `z` is a unit.

Every linear functional `lambda:V -> F_2` gives a codeword

```
c_j = lambda(z^j),          j in Z_N.           (6)
```

Conversely the first `q` entries of a recurrence word determine
`lambda` on the displayed basis, and the recurrence determines
the rest. This describes all of `C_s` bijectively. The next argument
uses this dual-coordinate parametrization; it does not identify
the coefficient vector of `a` in `h*a` with `lambda`.

Fix a nonzero `lambda`. For each `j in Z_N` form the ordinary
polynomial

```
a_j(z) = sum_(i=0)^(q-1) lambda(z^(j+i))*z^i.
```

It is nonzero: if all its coefficients vanished, `lambda` would
vanish on the basis `z^j,...,z^(j+q-1)`, which is the original
basis multiplied by a unit. Lemma 1 therefore applies to every
`a_j`.

For `0<=t<q` and `0<=r<=2q-2`, define

```
Q_(t,r)(z) = sum_(i=0)^(q-1) binom(t,r-i)*z^i,
```

where binomial coefficients are taken modulo two and are zero
outside `0<=r-i<=t`. The coefficient at `z^r` in the open-history
row `a_j(z)(1+z)^t` is

```
lambda(z^j * Q_(t,r)(z)).                      (7)
```

**Lemma 2.** In `V`, each `Q_(t,r)` is a sum of at most `s`
monomials `z^e`.

**Proof with the truncation indices explicit.** Put `B_t(z)=(1+z)^t`.
Its coefficient at exponent `v` is one precisely when `v` is a
binary submask of `t`, meaning that every one bit of `v` is also
a one bit of `t`. Symmetry of binomial coefficients gives

```
Q_(t,r) = truncate_[0,q)( z^(r-t) * B_t(z) ).   (8)
```

Set `a=r-t`. If `a>=0`, the retained exponents of `B_t` form the
prefix interval `[0,q-a)` in `[0,q)`, possibly empty. If `a<0`,
they form the suffix interval `[-a,q)`, again with the inherent
support restriction `v<=t`. Since `t<q`, at most one side of
`B_t` is cut.

Every prefix interval `[0,L)` in `[0,2^s)` has a partition into
at most `s` aligned dyadic intervals. For `0<L<2^s`, read the
binary expansion of `L` from the largest bit downward and use
one interval for each one bit; the running left endpoint is a
multiple of the next interval's length. There are at most `s`
one bits. The full interval is one dyadic interval, and the empty
interval requires none. Suffix intervals have the same property
by reflection, which preserves dyadic alignment.

Consider one dyadic interval

```
D = [b*2^h, (b+1)*2^h)
```

in this partition. If the fixed high bits `b` are incompatible
with the high bits of `t`, the intersection with the submasks of
`t` is empty. Otherwise let `u = t mod 2^h`; that intersection is
exactly

```
{ b*2^h + v : v is a submask of u }.
```

Its contribution after the shift in (8) is the polynomial

```
z^(a+b*2^h) * (1+z)^u.
```

The exponent `a+b*2^h` is nonnegative, and the degree is below
`q`, because the retained dyadic block lies inside the truncation
window. In `V` this polynomial equals the **single monomial**

```
z^(a+b*2^h+q*u),
```

using `1+z=z^q`. There are at most `s` blocks. This proves the
lemma. Repeated exponents may cancel in `V`; they only improve
the Hamming inequality below. All exponents used as coordinates
are understood modulo `N`.

## 5. Averaging closes the cyclic lower bound

Let `w=weight(c)` for the nonzero codeword (6). For fixed `t,r`,
Lemma 2 writes `Q_(t,r)=sum_(l=1)^K z^(e_l)` in `V`, with
`K<=s`. Equation (7) becomes a sum of at most `s` physical bits
`c_(j+e_l)`.

For binary bits, the indicator that their sum is one is at most
the sum of their individual indicators. Sum this inequality over
all `j in Z_N`. Every translation is a permutation of `Z_N`, so

```
sum_(j in Z_N) [lambda(z^j Q_(t,r))=1]
    <= K*w <= s*w.                            (9)
```

This remains true when `c` has a proper period dividing `N`:
the sum is over all `N` coordinates on both sides.

There are exactly `q*(2q-1)` possible pairs `(t,r)`. Summing (9)
over them and using Lemma 1 for every nonzero `a_j` gives

```
N*3^s <= sum_(j in Z_N) T_s(a_j)
       <= q*(2q-1)*s*w.
```

Since `q>=2`,

```
q*(2q-1) <= 2*(q^2-1) = 2*N.
```

Canceling the positive integer `N` yields the clean inequality

```
3^s <= 2*s*w.                                 (10)
```

It holds for every nonzero word, proving
`d_s>=ceil(3^s/(2s))`.

For completeness the additional bound `d_s>=q` is elementary.
Any `q` consecutive zeros propagate to the whole recurrence
word. Thus every cyclic window of length `q` meets a nonzero
word's support. Counting all `N` windows gives `q*w>=N`, and
`ceil(N/q)=q`. This completes (2).

## 6. BB distance, check weight, and increasing efficiency

The matrices of `A_s` and `B_s` are square circulant seed
matrices on separate coordinates. Their transpose kernels are
coordinate reversals, so they have the same parameters
`[N,q,d_s]`. The square-matrix hypergraph-product distance
theorem gives (1): see Kovalev and Pryadko,
[*Improved quantum hypergraph-product LDPC codes*](https://arxiv.org/pdf/1202.0928),
Section IV-C, equation (12), including its explicit treatment
of cyclic seeds immediately afterward. This uses the true
classical minimum and proves the true quantum minimum equals
it.

The exponents `0,1,q` are distinct modulo `N`. The two physical
blocks in a CSS check are disjoint, hence each check has exactly
six qubits. The support differences include unit steps along
both cyclic factors, so the Tanner graph is connected.

Now

```
R_s = q^2*d_s^2/(q^2-1)^2.
```

The lower distance bound and `q^2-1<q^2` give
`R_s>(9/4)^s/(4s^2)`. The upper distance bound gives
`d_s<3^s`; also `q^2-1 >= (3/4)*q^2` for `q>=2`.
These imply the upper inequality in (3). The lower bound is
unbounded because exponential growth dominates its quadratic
denominator.

Here is an elementary explicit comparison proving (4), without
assuming that consecutive actual distances are given by the
upper witness. From (3),

```
R_(2s) > (9/4)^(2s)/(16*s^2),
R_s    < (16/9)*(9/4)^s.
```

For every integer `s>=16`,

```
(9/4)^s > 2^s >= 256*s^2 > (256/9)*s^2.       (11)
```

The middle inequality is equality at `s=16`. It propagates by
induction because `2*s^2 >= (s+1)^2` for `s>=3`.
Equation (11) makes the lower bound on `R_(2s)` strictly larger
than the upper bound on `R_s`. Taking `s=16*2^t` proves (4).

## 7. Exact finite values and the formal endpoint

The first five exact classical minima, established by complete enumeration
of the explicit kernel, are `2, 8, 26, 80, 242` at `s=1,2,3,4,5`.
The [finite evidence note](BB_GROWING_FOM_FRACTAL_CANDIDATE.md#3-exhaustive-finite-evidence)
and [reproducible checker](../scripts/bb_growing_fom_seed_check.py) record
attaining words and exhaustive counts. These computations are not Lean
certificates and are not premises of the universal Lean theorem.

The compiled formal endpoint uses exactly the concrete BB checks above.
Its distance is the attained minimum among nonzero words of the explicit
recurrence space; it is neither an arbitrary parameter nor a supplied
distance hypothesis. The quantum proof establishes that the minimum weight
of a nontrivial Pauli logical operator equals this classical minimum.

For positive mathematical `s`, the formal statements prove

```
[[2*(4^s-1)^2, 2*4^s, d_s]],
0 < d_s <= 4^s-1,
3^s <= 2*s*d_s <= 2*s*(4^s-1).
```

The rightmost inequality is the immediate consequence of the formal length
bound. The sharper analytic statements `d_s >= 2^s` and `d_s <= 3^s-1`
from Sections 2 and 5 have not been formalized in this endpoint. Likewise,
connectedness is justified analytically here; the exact Pauli weight-six
generating set is formalized.

The formal selected subsequence uses

```
s_t = 16*4^t,
growingCode t : StabilizerCodeWithDistance
  (2*(4^(s_t)-1)^2) (2*4^(s_t)) d_(s_t).

StrictMono growingEfficiency
(bound : Q) < growingEfficiency (bound+1), for every bound : N.
```

Here `growingEfficiency` is proved equal to the actual `k*d^2/n` of the
packaged parameters. The larger index step lets the Lean proof use only
the elementary upper bound `d_s <= N_s`. The sharper `s_t=16*2^t` proved
analytically in Section 6 also uses the unformalized upper witness. These
are two valid subsequences, with distinct proof scopes.

All following declarations are in namespace
`Quantum.Stabilizer.Homological.BB.Fractal`. The local library code is in
the isolated worktree `/Users/stavanjain/Code/QECLean/.Codex/worktrees/bb-growing-fom`.

| Formal component | Declarations and code |
|---|---|
| Actual classical minimum and attainment | `seedDistance`, `exists_seedDistance_word`, `seedDistance_isLeast` in [Distance.lean](https://github.com/Stavan-Jain/QECLean/blob/1e7f44460ad65e9ec6b1e65a198faf79477f8141/QEC/Stabilizer/Codes/BivariateBicycle/Fractal/Distance.lean) |
| Uniform physical weight bound | `functionalWord_weight_lower`, `seedSpace_weight_lower` in [LowerBound.lean](https://github.com/Stavan-Jain/QECLean/blob/1e7f44460ad65e9ec6b1e65a198faf79477f8141/QEC/Stabilizer/Codes/BivariateBicycle/Fractal/LowerBound.lean) |
| Exact Pauli distance and all-code bundle | `chainComplex_pauli_distance`, `stabilizerCodeWithDistance` in [Code.lean](https://github.com/Stavan-Jain/QECLean/blob/1e7f44460ad65e9ec6b1e65a198faf79477f8141/QEC/Stabilizer/Codes/BivariateBicycle/Fractal/Code.lean) |
| Actual weight-six stabilizer generating set | `stabilizerCodeWithDistance_has_weight_six_generating_set` in [CheckWeight.lean](https://github.com/Stavan-Jain/QECLean/blob/1e7f44460ad65e9ec6b1e65a198faf79477f8141/QEC/Stabilizer/Codes/BivariateBicycle/Fractal/CheckWeight.lean) |
| Selected bundle and actual ratio growth | `growingCode`, `growingEfficiency_eq_kd2_div_n`, `growingEfficiency_strictMono`, `growingEfficiency_unbounded` in [Family.lean](https://github.com/Stavan-Jain/QECLean/blob/1e7f44460ad65e9ec6b1e65a198faf79477f8141/QEC/Stabilizer/Codes/BivariateBicycle/Fractal/Family.lean) |

The underlying code definitions use a zero-based Lean index `u`, so
`stabilizerCodeWithDistance u` has mathematical index `s=u+1` and distance
`seedDistance u`. `distanceAt s = seedDistance (s-1)` restores the positive
mathematical notation for `growingCode` and its growth theorems.

The capstones compile and axiom audits report only `propext`,
`Classical.choice`, and `Quot.sound`. All six required repository checks
and the umbrella check have passed, with zero new compiler/linter diagnostics.
See the [validation record](../../../pipeline/attempts/bb_growing_fom/validation.md).

Companion drafts: [QECLean #86](https://github.com/Stavan-Jain/QECLean/pull/86) and [qec-lab #32](https://github.com/Stavan-Jain/qec-lab/pull/32).
