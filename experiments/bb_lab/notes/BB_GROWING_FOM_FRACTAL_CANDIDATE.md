# Cyclic fractal BB codes: finite exact distances and the remaining conjecture

2026-09-13. Research follow-up to qec-lab PRs #27 and #28. No novelty claim
is made over the fractal-code and hypergraph-product literature.

**Current status.** The universal lower bound `3^s <= 2*s*d_s`, the actual
attained quantum distance `d_s`, and an unbounded, strictly increasing family
are now proved in Lean. The packaged subsequence is `s_t=16*4^t`; its checks
have weight six. The capstone modules compile and their audited declarations
use only `propext`, `Classical.choice`, and `Quot.sound`. All required repository checks passed, with no new diagnostics; pipeline
status is `pr-ready`.

See [the theorem and formal endpoint](BB_FRACTAL_UNBOUNDED_FOM_THEOREM.md#7-exact-finite-values-and-the-formal-endpoint)
for precise statements, code links, and the distinction between formal and
analytic bounds. The sharper analytic bounds imply strict increase already
on `s_t=16*2^t`. The exact formula `d_s=3^s-1` remains unproved in this development.

This note retains the construction, reproducible finite computations, and
that remaining exact-distance conjecture. Finite enumeration is separate
from Lean certification and is not used by the universal proof. The family
does not require infinitely many primitive trinomials.

## 1. Construction and unconditional algebra

For an integer `s >= 1`, set

```
q = 2^s,       N = q^2 - 1,
p_s(z) = 1 + z + z^q in F_2[z].
```

On `G = Z_N x Z_N`, define the BB pair

```
A = p_s(x),       B = p_s(y),
H_X = [A | B],    H_Z = [B^T | A^T].
```

Here `A` and `B` denote the corresponding translation matrices. The three
exponents `0,1,q` are distinct modulo `N`, including `s=1`. Each displayed
check consequently has exactly six physical qubits in its support. The support
differences contain the unit step in both cyclic factors, so the full Tanner
graph is connected. This is a separated-variable BB code, also a square-matrix
hypergraph product. The construction does not use a sum of disconnected copies.

In characteristic two,

```
p_s(z)^q + p_s(z) = z^(q^2) + z = z*(z^N + 1).
```

As `p_s(0)=1`, it follows that `p_s` divides `z^N+1`. More explicitly,

```
h_s(z) = (p_s(z)^(q-1) + 1)/z,
p_s(z)*h_s(z) = z^N+1.                          (1)
```

The division by `z` is legitimate: the constant terms cancel. If
`R = F_2[z]/(z^N+1)`, cancellation in `F_2[z]` gives

```
ker(p_s : R -> R) = {h_s(z)*a(z) : deg a < q}.  (2)
```

The `q` shifted polynomials `z^j h_s`, `0 <= j < q`, are linearly independent
because their leading degrees are distinct, and all have degree below `N`.
Thus the classical kernel has dimension exactly `q`.

Reversing the cyclic coordinate identifies the kernel of the transpose check
matrix with this kernel and preserves physical Hamming weight. Write its
classical distance as `d_s`. The BB code therefore has exact parameters

```
[[2*(4^s-1)^2, 2*4^s, d_s]].                    (3)
```

The distance transfer in (3) is the square-matrix hypergraph-product result:
for two square seeds whose transposes have the same parameters, the product
has `[[2*n1*n2,2*k1*k2,min(d1,d2)]]`. See Kovalev and Pryadko,
[*Improved quantum hypergraph-product LDPC codes*](https://arxiv.org/pdf/1202.0928),
Section IV-C, equation (12); the following paragraph addresses full circulant
check matrices explicitly. This is a mathematical transfer theorem, not a
quantum-distance search or a new theorem established by the checking script.

## 2. A closed upper witness at every s

The binary expansion of `q-1` gives

```
p_s(z)^(q-1) = product_(j=0)^(s-1)
                   (1 + z^(2^j) + z^(q*2^j)).  (4)
```

Each selected monomial has exponent `a + q*b`, where `0 <= a,b < q` and
the binary supports of `a` and `b` are disjoint. Conversely every such pair
appears once. Uniqueness of the base-`q` expansion prevents collisions or
cancellations. Each bit has three choices, so (4) has exactly `3^s` terms.
Equation (1) deletes its single constant term and shifts every remaining
term down one place. Consequently

```
weight(h_s) = 3^s - 1,
d_s <= 3^s - 1.                                (5)
```

This is an upper bound only. In particular, the weight of a cyclic generator
does not prove the minimum weight of arbitrary nonzero linear combinations
of its shifts.

A simple unconditional lower bound is `d_s >= q`. Indeed the recurrence
defined by `p_s` propagates any run of `q` zero entries to the whole cyclic
word. Every cyclic interval of length `q` therefore meets a nonzero word's
support. Counting incidences with all `N` intervals gives `q*weight >= N`,
and `ceil(N/q)=q`. This bound does not establish growing figure of merit.

## 3. Exhaustive finite evidence

The standalone script
[`bb_growing_fom_seed_check.py`](../scripts/bb_growing_fom_seed_check.py)
checks the polynomial identity two ways and exhaustively enumerates the
nonzero messages in (2) when `q <= 16` using the Python standard library.
Its optional native mode compiles an embedded C++ ARM NEON loop to enumerate
all `2^32-1` nonzero messages at `s=5`; it requires an ARM64 host and a C++
compiler, with no external libraries.

```
uv run --project experiments/bb_lab python \
  experiments/bb_lab/scripts/bb_growing_fom_seed_check.py

# Full integer polynomials, attaining words, and weight distributions:
uv run --project experiments/bb_lab python \
  experiments/bb_lab/scripts/bb_growing_fom_seed_check.py --json

# Also prove the finite s=5 minimum by a complete native enumeration:
uv run --project experiments/bb_lab python \
  experiments/bb_lab/scripts/bb_growing_fom_seed_check.py --native-q32 --json
```

Each Gray-code step changes exactly one basis vector; the run visits every
nonzero `q`-bit message once. The largest pure-Python enumeration contains
65,535 nonzero words. The native `s=5` run contains 4,294,967,295 nonzero
words, with minimum 242 attained 1,023 times. Its aggregate weight is
2,196,875,771,904, independently matching `1023*2^31`: every coordinate
is a nonzero linear functional of the message and is one in half of all
messages. The Python wrapper also independently reconstructs and checks
the attaining word returned by C++. Completeness follows from (2), not from a heuristic
minimum-distance solver. These are reproducible exact computer calculations;
the script itself is not a kernel-checked Lean proof.

| s | Classical seed | BB code, via (3) | k*d^2/n | Status |
|---|---|---|---:|---|
| 1 | `[3,2,2]` | `[[18,8,2]]` | 1.777778 | exact finite enumeration + HGP theorem |
| 2 | `[15,4,8]` | `[[450,32,8]]` | 4.551111 | exact finite enumeration + HGP theorem |
| 3 | `[63,8,26]` | `[[7938,128,26]]` | 10.900479 | exact finite enumeration + HGP theorem |
| 4 | `[255,16,80]` | `[[130050,512,80]]` | 25.196463 | exact finite enumeration + HGP theorem |
| 5 | `[1023,32,242]` | `[[2093058,2048,242]]` | 57.303272 | exact native enumeration + HGP theorem |
| 6 | `[4095,64,<=728]` | `[[33538050,8192,<=728]]` | unknown | upper witness only |

Checked locally on 2026-09-13, including a second complete native `s=5` run
through the durable Python wrapper. The default invocation does not run
the optional native enumeration and consequently reports only its newly
checked upper bound for `s=5`. The script deliberately prints
`unknown` for the figure of merit of rows whose distance is only upper-bounded.
It does not silently promote the witness weight to a minimum distance.

## 4. The remaining exact-distance conjecture

The candidate statement is

```
For every s >= 1 and every nonzero a in F_2[z] with deg a < 2^s,
    weight(h_s*a) >= 3^s - 1.                   (OPEN)
```

There is no modular reduction in this inequality: all products have degree
at most `N-1`. Together with (5), it would prove

```
d_s = 3^s - 1,
R_s = k*d_s^2/n
    = 4^s*(3^s-1)^2/(4^s-1)^2
    ~ (9/4)^s.                                 (CONDITIONAL)
```

These proposed ratios are strictly increasing. To check this, compare their
positive square roots. Writing `a=3^s`, `b=4^s`, the numerator of the
cross-multiplied difference is

```
2*(3*a-1)*(b-1) - (a-1)*(4*b-1)
    = 2*a*b + 2*b - 5*a + 1 > 0,
```

because `b >= 4`. Unboundedness follows, for example, from
`sqrt(R_s) > (3^s-1)/2^s >= (2/3)*(3/2)^s`.

The exact closed expression for `R_s` and strict increase at every consecutive
positive `s` remain conjectural beyond the finite rows. The universal weight
inequality above would justify a bundle with the explicit closed distance
`3^s-1`. It is not needed for the existing unconditional
`StabilizerCodeWithDistance` bundle at the attained minimum `d_s`, whose
selected actual ratios are already proved strictly increasing and unbounded.

## 5. Relation to the rigorous fractal-code literature

Yoshida's [*Information storage capacity of discrete spin systems*](https://arxiv.org/pdf/1111.3275),
Theorem 2, proves an exact classical fractal-code distance. In its binary
case, `L=2^s`, the parameters are `[2L^2,L,3^s]`. Lemma 2 controls the weight
of a linear combination of Pascal rows by the first nonzero coefficient;
its polynomial form supplies the open-history ingredient used in the
formal lower bound.

The relevant boundary conditions are essential. Yoshida's equations (15)
and (17) specify a finite time history with an initial row whose right half
is zero. Our classical seed is instead a cyclic recurrence
`c_(i+q)=c_i+c_(i+1)`. Viewed as a two-dimensional local rule, its period
lattice is generated by `(q,-1)` and `(-1,q)`, of determinant `q^2-1`.
The Smith normal form is cyclic, explaining its one-variable description.
An arbitrary seed word does not have Yoshida's initial zero half. No
weight-preserving equivalence or suitable shortening argument has been
established here. Thus that published theorem does **not** discharge OPEN.

This distinction also separates the candidate from a hypergraph product
of two ordinary rectangular Newman–Moore codes: their translation group
normally has four cyclic factors, whereas the present BB group has two.

## 6. A limited infinite consequence already available from a fixed seed

For any fixed exact seed `[N,q,d]` above and positive integers `r,t`, use
the same trinomials on `Z_(Nr) x Z_(Nt)`. Formula (2) becomes the literal
repetition of length-`N` words, since

```
(z^(Nr)+1)/p_s = h_s * sum_(j=0)^(r-1) z^(j*N).
```

This identity remains valid for even `r`; extra multiplicity of roots in
the ambient modulus creates no extra kernel dimension. The exact BB
parameters are therefore

```
[[2*N^2*r*t, 2*q^2, d*min(r,t)]].              (6)
```

Taking `t=r+1` gives a strictly increasing figure of merit
`(q^2*d^2/N^2)*r/(r+1)` and growing distance. With the `s=4` finite seed,
this is `[[130050*r*(r+1),512,80*r]]` and the ratio increases towards
25.196463. The `s=5` exact seed similarly gives
`[[2093058*r*(r+1),2048,242*r]]` with ratio increasing towards 57.303272.
These are aspect-ratio improvements with finite limits, retained only as
an earlier fallback. The main family's unbounded improvement follows from
the proved weaker distance bound, independently of this fallback or the
open exact-distance conjecture.
