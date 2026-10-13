# Audit of the cyclic fractal distance lower bound

Date: 2026-09-13. Status: **the uniform distance lower bound and the exact-distance BB endpoint now compile in Lean; the exact closed formula remains conjectural**.

The formal endpoint proves `3^s <= 2*s*d_s` and `0<d_s<=4^s-1`, where
`d_s` is an attained minimum, and packages a weight-six family with strictly
increasing unbounded actual `k*d^2/n` on `s=16*4^t`. The sharper analytic
bounds and subsequence `s=16*2^t` remain distinct. See the
[formal endpoint and code links](BB_FRACTAL_UNBOUNDED_FOM_THEOREM.md#7-exact-finite-values-and-the-formal-endpoint)
for exact declarations. Module compilation and axiom audits have passed;
all required repository checks passed and pipeline status is `pr-ready`. No novelty over Yoshida or hypergraph products is claimed.

This note separates the algebraic upper witness from experimental exact-distance evidence, records failed direct proof transfers, and independently audits the successful weaker lower bound. The main result is developed in [BB_FRACTAL_UNBOUNDED_FOM_THEOREM.md](BB_FRACTAL_UNBOUNDED_FOM_THEOREM.md).

## Candidate and the exact missing statement

For `s >= 1`, put

- `q = 2^s`, `N = q^2 - 1`;
- `p(z) = 1 + z + z^q` over `F_2`;
- `C_s = ker(multiplication by p)` in `F_2[z]/(z^N + 1)`.

The separable BB choice `A = p(x)`, `B = p(y)` on `Z_N x Z_N` has six entries in each check row. Its classical factors have dimension `q`, so the usual square hypergraph-product parameter calculation gives quantum length `2 N^2` and dimension `2 q^2`. If the classical distance is `d_s`, its quantum distance is `d_s`.

The proposed exact formula is

`d_s = 3^s - 1`.

It would imply

`k d_s^2 / n = q^2 (3^s - 1)^2 / (q^2 - 1)^2`,

which grows asymptotically as `(9/4)^s`. The unresolved core is the following elementary-looking polynomial inequality:

> For every nonzero `a in F_2[z]` of degree less than `q`,
> `wt(a * ((1 + z + z^q)^(q-1) + 1)) >= 3^s - 1`.

Here `wt` counts nonzero coefficients of an ordinary polynomial; there is no quotient reduction in this inequality. Multiplication by `z` preserves coefficient weight, so this is exactly the lower-bound assertion for all nonzero multiples `a h`, with `h` defined below. The weaker bound `3^s <= 2 s d_s`, proved below, suffices for unbounded `k d_s^2/n`. It does not establish the proposed exact-distance formula.

## Proved: divisibility, dimension, and a weight `3^s - 1` witness

Frobenius gives

`p^q + p = z^(q^2) + z = z (z^N + 1)`.

Since `p(0) = 1`, the polynomials `p` and `z` are coprime. Thus `p` divides `z^N + 1`, and

`h = (z^N + 1)/p = (p^(q-1) + 1)/z`.

In particular, `deg h = N-q`, the ideal `(h)` is exactly `C_s`, and its polynomial basis is `h, z h, ..., z^(q-1) h`. These basis elements have different leading degrees, so the classical dimension is `q`.

Because `q-1` has `s` one-bits,

`p^(q-1) = product_(i=0,...,s-1) (1 + z^(2^i) + z^(q 2^i))`.

Choosing a term from each factor produces an exponent `u + q v`, where `0 <= u,v < q` and the binary digits of `u` and `v` are disjoint. These exponents are all distinct: their low and high `s` binary digits recover the choices. Therefore `p^(q-1)` has exactly `3^s` monomials. Its constant coefficient is one, so adding one removes precisely that monomial. Division by `z` shifts the remaining exponents without changing their count. Consequently

`wt(h) = 3^s - 1`, and hence `d_s <= 3^s - 1`.

This is a rigorous upper bound, not a proof of minimum distance.

## Empirical exact distances, including a complete `s=5` search

Exhaustion of all nonzero message polynomials gives:

| s | q | N | Tested minimum | Number of minimum words |
|---|---|---|---|---|
| 1 | 2 | 3 | 2 | 3 |
| 2 | 4 | 15 | 8 | 15 |
| 3 | 8 | 63 | 26 | 63 |
| 4 | 16 | 255 | 80 | 255 |
| 5 | 32 | 1023 | 242 | 1023 |

The `s=5` check used the ARM NEON Gray-code enumerator now embedded in the
[durable checker](../scripts/bb_growing_fom_seed_check.py), invoked with
`--native-q32`. Its loop variable and bound were explicitly 64-bit: `uint64_t i = 1; i < (1ull << 32); ++i`. Each step toggled the generator numbered `ctz(i)`, covering all messages exactly once.

An audited rerun reported:

```text
min=242
minimum_word_count=1023
total_words=4294967295
total_weight=2196875771904
expected_total_weight=1023 * 2^31 = 2196875771904
```

The total-weight checksum is independently required because each physical coordinate is a nonzero linear functional on the 32-dimensional message space. The SIMD accumulator cannot overflow: each byte collects at most eight byte-popcounts, hence at most 64, followed by a widened horizontal sum. The generator was constructed directly from the distinct exponents `u + 32 v - 1`, with `u & v = 0` and `(u,v) != (0,0)`.

These finite computations support the conjecture. They are not Lean certificates and do not imply the statement for arbitrary `s`.

## Why the immediate Pascal-matrix argument fails

[Yoshida, *Information storage capacity of discrete spin systems*, arXiv:1111.3275](https://arxiv.org/pdf/1111.3275) proves a fractal-code distance theorem with specified boundary restrictions. Its principal-vector lemma controls combinations of rows of a Pascal matrix by the first nonzero coefficient in that basis. Those restrictions do not automatically hold for the cyclic code above.

Writing a cyclic codeword as a `q x q` table using coordinates `f(i,j) = c_(i+q j mod N)` duplicates the corner `c_0`. Interior entries obey a Pascal recurrence, but there are also essential skew boundary equations:

`f(0,j) = f(q-1,j-1) + f(0,j-1)` for `j >= 1`,

`f(i,0) = f(i-1,0) + f(i-1,q-1)` for `i >= 1`.

Interior recurrence and equality of the two opposite corners do not suffice. For `q=4`, the table

```text
1 0 0 1
0 0 0 1
0 0 0 1
0 0 0 1
```

has both opposite corners equal to one and satisfies every interior recurrence, but has weight 5 instead of 9. More generally, the analogous table with its top-left entry and its entire last column set to one has weight `q+1`. It fails the required skew boundary equations.

There is also a precise form of the folding obstruction. Expand a message in the principal basis `a = sum_j c_j (1+z)^j` and write `B(t)=(1+z)^t`. Split `a (p^(q-1)+1)` into `q` blocks of `q` coefficients, labelling a block by `t=q-1-b`, where `b` is its ordinary block index. A principal message `a=B(j)` contributes

- `B(j+t)` if `j+t < q-1`;
- `B(q-1)+1` if `j+t = q-1`;
- `z B(j+t-q)` if `j+t >= q`.

Each principal message separately gives the desired weight. General messages mix the high and low Pascal pieces within a block, allowing cancellations. Applying the ordinary principal-vector lemma independently to these folded blocks does not establish the claimed sum of weights. This is the current concrete missing inequality, rather than an omitted routine boundary check.

## Why the Sidon-gap route is currently insufficient

The roots of `p` in `F_(q^2)` are exactly the affine line

`S = {lambda : lambda^q + lambda = 1}`.

Under a discrete logarithm, this is a Bose–Chowla Sidon set of size `q` in `Z_(q^2-1)`. A long gap in these exponents, or in a unit dilate, would give a BCH lower bound because the complementary exponents are zeros of the cyclic code.

However, the available gap result inspected here only improves the elementary average-gap scale by a constant factor. [Shakan's research note on large gaps in dilates](https://blog.georgeshakan.com/a-large-gap-in-a-dilate-of-a-set/) discusses a bound on the scale `2N/|S|`, hence approximately `2q` here, and describes the hoped-for unbounded improvement as a further problem. A constant multiple of `q` does not imply unbounded `q^2 d_s^2/N^2`. Thus this route supplies no uniform superlinear distance bound in the present audit.

## Proved: a weaker lower bound by averaging shifted open histories

The failed direct folding argument can be bypassed. The following proof was independently checked against the successful construction in the companion theorem note.

### Physical coordinates and the open history

Let `R = F_2[z]/(p)`. Since `p` divides `z^N+1`, multiplication by `z` is invertible and has period dividing `N`. Every linear functional `ell : R -> F_2` determines a periodic recurrence word

`c_j = ell(z^j)`.

The first `q` entries are arbitrary because `1,z,...,z^(q-1)` form a basis. These words satisfy `c_(j+q) = c_j + c_(j+1)` and give the transpose version of the cyclic kernel. Reversal of coordinates identifies it with the original kernel and preserves weight.

For any nonzero initial row `a=(a_0,...,a_(q-1))`, define its ordinary open history by the coefficients of

`a(z)(1+z)^t`, for `0 <= t < q`.

Each row has degree at most `2q-2`. Its total history weight `T(a)` is at least `3^s`. Here is a self-contained binary proof of the needed Pascal estimate. In a row of length `2^m`, expand a vector in the basis `B(j)=(1+z)^j`, and let `j0` be its least nonzero basis coefficient. Induction on `m`, using the two halves `(u+v,v)`, shows its weight is at least `2^popcount(j0)`: if the first half of the coefficients is nonzero, use `wt(u+v)+wt(v)>=wt(u)`; otherwise both output halves equal `v`. The basis is triangular with diagonal one, so the expansion exists uniquely.

Apply this estimate at length `2q` to `a(z)(1+z)^t`. The first principal index is `j0+t`, where `j0<q`. Consequently

`T(a) >= sum_(t=0,...,q-1) 2^popcount(j0+t) >= 3^s`.

The last inequality follows by splitting the interval at `q`: terms that wrap from `q+r` have twice the weight of the omitted terms `r`, while `sum_(r<q)2^popcount(r)=3^s`.

### Explicit dyadic truncation lemma

Write `v subseteq t` when every one-bit of `v` is a one-bit of `t`. Lucas' binary identity gives

`(1+z)^t = sum_(v subseteq t) z^v`.

A prefix or suffix of `[0,q)` is a disjoint union of at most `s` aligned dyadic intervals. For a proper prefix this is the usual binary expansion of its length; a suffix follows by reflection. The full interval is one block and `s>=1`.

Consider one such block `[b 2^h,(b+1)2^h)`. Its intersection with the submasks of `t` is empty unless the fixed high bits `b` are allowed by the high bits of `t`. If compatible, the intersection is exactly

`{b 2^h + v : v subseteq (t mod 2^h)}`.

For an open-history coordinate `(t,r)`, its functional on the initial row is represented by

`Q_(t,r)(z) = sum_(i=0,...,q-1) binom(t,r-i) z^i`.

Equivalently, set `i=r-v` and retain submasks `v` of `t` for which `0<=r-v<q`. There is only one-sided truncation: if `r<q`, it is the prefix `v<=r`; if `r>=q`, the upper restriction is automatic since `v<=t<q`, leaving a suffix. Thus at most `s` dyadic components suffice.

For a compatible component, put `u=t mod 2^h`. Reversing its exponents gives the polynomial

`z^(r-b 2^h-u) (1+z)^u`.

Its first exponent is nonnegative and its degree is below `q`, because all component terms were retained in the original coordinate window. In `R`, the relation `1+z=z^q` turns this polynomial into the single monomial

`z^(r-b 2^h-u+q u)`.

Therefore every open-history coordinate is the XOR of at most `s` physical coordinates of `c`, with fixed offsets independent of the choice of the functional `ell`.

### Averaging and cancellation of the period

For every cyclic shift `j`, use the nonzero initial row `(c_j,...,c_(j+q-1))`. It is nonzero: otherwise its functional vanishes on a basis, contradicting invertibility of multiplication by `z^j` and nonzero `ell`. Hence each shifted open history has weight at least `3^s`.

For a fixed open-history coordinate, the dyadic lemma bounds its indicator of being one by the sum of at most `s` indicators of shifted physical coordinates. Sum over all `N` shifts. Each individual physical-coordinate offset contributes exactly `wt(c)` by periodicity. There are `q(2q-1)` possible history coordinates, so

`N 3^s <= q(2q-1) s wt(c) <= 2N s wt(c)`.

The last inequality is `q>=2`. Canceling positive `N` proves

**`3^s <= 2s wt(c)` for every nonzero recurrence word.**

Thus `3^s <= 2s d_s`. This is a uniform distance proof, independent of the finite searches.

## Strict growth and the formal audit

The [main theorem note](BB_FRACTAL_UNBOUNDED_FOM_THEOREM.md#6-bb-distance-check-weight-and-increasing-efficiency)
gives the analytic inequalities

`(9/4)^s/(4s^2) < R_s < (16/9)(9/4)^s`,

where `R_s=q^2*d_s^2/(q^2-1)^2`. They imply unboundedness and strict increase
on `s_t=16*2^t`. This argument uses the analytic witness `d_s<=3^s-1`.

The Lean proof instead combines the formal floor `3^s<=2*s*d_s` with the
formal length bound `d_s<=N_s`. It proves strict increase on `s_t=16*4^t`
and gives the explicit unboundedness witness `R_(s_(bound+1))>bound` for
every natural threshold. The attained minimum `d_s` is the exact quantum
distance in the packaged `StabilizerCodeWithDistance`.

The independent source review checked the completed `History`, `Truncation`,
`LowerBound`, and `Averaging` modules for the following essential points:

- Tap offsets are independent of the functional and the translation.
- Translation multiplies by a unit, so every translated initial polynomial
  is nonzero.
- The history coefficient identity is exact, including boundary truncation.
- The average includes every physical translation and cancels the positive
  cyclic period; no weight-preserving change of basis is assumed.
- The minimum is attained, and the quantum-distance transfer proves both
  a universal lower bound and an attaining nontrivial Pauli operator.

The audited lower-bound and growth capstones depend only on `propext`,
`Classical.choice`, and `Quot.sound`. The remaining repository-wide
validation is tracked in the [pipeline result](../../../pipeline/attempts/bb_growing_fom/result.md).
The stronger equality `d_s=3^s-1` remains unproved here and is absent from all
formal endpoint assumptions.
