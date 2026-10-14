# A42 S18: full-column pair leaders — sound reduction, partial certificates

2026-09-08. Bounded Codex analytic pass. Only this note and
`scripts/a42_s18_codex_pair_leaders.py` are new outputs.

## 0. Disposition and honest tiers

The proposed implication from three pair leaders to full canonical
hiding mass is **correct**. This pass does **not** prove the three
unrestricted pair-leader assertions, and finds no counterexample to
them on the original unbounded cylinder.

The following narrower statements are proved:

1. All three unrestricted pair assertions hold at q=1 and q=2.
   The q=2 proof uses the banked canonical exact-gap theorem, not
   a coset search.
2. For every positive height q, all three pair inequalities hold
   for every boundary generator supported in x columns `[-3,-1]`,
   with arbitrary binary row masks. Explicit triangle-inequality
   certificates are given below. This is a generator restriction,
   not a width-completeness theorem.
3. At every height, a generator supported on a single arbitrary
   x column cannot violate any pair inequality.
4. Folding x modulo three cannot prove the unrestricted assertions:
   AB and AC have explicit shorter folded representatives at every
   dyadic q>=8, while BC becomes a folded boundary at every height.
   The AB/AC examples have a displayed nonzero original-syndrome
   residual. They are **not** original-cylinder counterexamples.

The all-excess original pair route therefore remains undecided.
Neither the fixed strip nor the failed quotient relaxation is
presented as a solution of the remaining metric problem.

## 1. Setup and audit of the sufficient reduction

Let q be dyadic unless stated otherwise. Work over

```
R_q[x,x^-1],  R_q=F2[t]/(t^q+1),
u=t+t^2,     v=t^3,     X=x^3,
P=X+u,       Q=1+(x+x^2)v,
N=1+t+...+t^(q-1),     chi=1+x+x^2.
```

The boundary of a compact Laurent generator h is
`s=(Qh,Ph)`. All weights are physical binary Hamming weights,
not numbers of occupied parity residues. Define full columns

```
A=(N,0),       B=(0,N),       C=(0,xN),
S=A+B+C,       T_AB=A+B, T_AC=A+C, T_BC=B+C.
```

The three supports are disjoint, and each pair has weight 2q.
For any boundary s and any pair T,

```
|T+s| = 2q+|s|-2|s intersect T|.                    (1.1)
```

Thus the pair-leader assertion is precisely
`2|s intersect T|<=|s|` for every compact boundary s.
Every coordinate of S belongs to exactly two pairs. Summing
the three inequalities gives

```
4|s intersect S| <= 3|s|
                 = 3|s intersect S|+3|s outside S|,
```

hence full canonical HM, `|s intersect S|<=3|s outside S|`.
This is the same boundary formulation as S8. It implies its
proper-subset formulation as well: for T contained in S, put
J=S intersect supp(s), use `|T|<=|J|+|T symmetric_difference J|`,
and then `|T|<=3|T+s|`. No representative minimality is assumed.

The exact pair syndromes are

```
D T_AB=(x^3+chi)N=(1+x)^3 N,
D T_AC=(x^3+x chi)N=(x+x^2)N,
D T_BC=(1+x^3)N.                                   (1.2)
```

The simplification of AC to the adjacent-column binomial is useful;
none of these actual columns has been identified with another.

### 1.1 Individual hole-column formulation

Put `e_star=(x^-3 N,0)`. The exact identity

```
S+e_star = boundary(x^-3(1+x)N)                     (1.3)
```

also shows that BC is in the same boundary coset as
`((1+x^-3)N,0)`, the two first-block norm lines at x=-3 and x=0.

For dyadic q, `P^q=x^(3q)`, so P is a Laurent unit with inverse
`x^(-3q)P^(q-1)`. Every compact syndrome-zero chain is therefore
a boundary: if `Pf+Qg=0`, take `h=P^-1 g`.
Consequently every compact solution `D e=N` is uniquely of the
form `e=S+s` for a boundary s. Write

```
delta=|e|-q,
k_A=|e intersect A|, k_B=|e intersect B|, k_C=|e intersect C|.
```

For the pair complementary to a column R in `{A,B,C}`,

```
|T_pair+s|=|e+R|=2q+delta-2k_R.                    (1.4)
```

Thus the proposed three assertions are exactly the three separate
estimates `2k_A<=delta`, `2k_B<=delta`, `2k_C<=delta`. Their sum
is the canonical HM estimate. This explains the stronger premise:
it budgets each individual hole column, rather than only their sum.
Equation (1.4) is an equivalence, not a proof of these estimates.

## 2. The two smallest heights, without enumeration

**Proposition 2.1.** Each pair is an unrestricted coset leader at
q=1 and at q=2.

At q=1, every boundary has even total weight: augmentation gives
`Q(1,1)+P(1,1)=1+1=0`. A pair also has even weight. The three
syndromes in (1.2) are nonzero Laurent polynomials, so a member
of a pair coset cannot have weight zero. Its weight is therefore
at least two, attained by the pair itself.

At q=2, suppose a pair representative z had weight less than four.
Its weight is even, again by augmentation, hence at most two.
Let R be the omitted full column. Then `e=z+R` has syndrome N
and weight at most four. The banked theorem in
`IBM_CANONICAL_EXACT_GAP.md` says that the leader has weight two
and every other solution has weight at least six. Therefore
`e=e_star`. But then `z=e_star+R` has weight four: e_star and
each of A, B, C have disjoint supports. This is a contradiction.

The q=2 invocation is exactly within the reviewed gap theorem.
It does not assert a scale-induction step to larger heights.

## 3. Three-column generator theorem at arbitrary height

**Theorem 3.1.** Let q be any positive integer, not necessarily
dyadic. For every a,b,c in R_q, the boundary generated by

```
h=x^-3 a+x^-2 b+x^-1 c                              (3.1)
```

satisfies all three pair inequalities.

**Proof.** The complete block-column expansion is

| x column | first block Qh | second block Ph |
|---:|---|---|
| -3 | a | u a |
| -2 | b+v a | u b |
| -1 | c+v(a+b) | u c |
| 0 | v(b+c) | a |
| 1 | v c | b |
| 2 | empty | c |

The rows in this table are distinct actual x columns. Multiplication
by v is a cyclic row permutation, so it preserves Hamming weight.
Set

```
alpha=|a|, beta=|b|, gamma=|c|,
R1=|b+v a|, R2=|c+v(a+b)|,
Z=|b+c|, V=|a+b|,
U=|u a|+|u b|+|u c|.
```

Direct counting, with all cancellations already contained in the
displayed masks, gives

```
|s|=2alpha+beta+2gamma+R1+R2+Z+U,
I_AB=alpha+Z, I_AC=beta+Z, I_BC=alpha+beta.           (3.2)
```

Here `I_pair=|s intersect T_pair|`. The desired energy differences
have exact decompositions into nonnegative quantities:

```
|s|-2I_AB = (beta+gamma-Z)+gamma+R1+R2+U,
|s|-2I_BC = (Z+gamma-beta)+gamma+R1+R2+U,
|s|-2I_AC = (R1+alpha-beta)
            +(R2+gamma-V)
            +(V+alpha-beta)
            +(beta+gamma-Z)+U.                    (3.3)
```

Every parenthesis is nonnegative by the binary Hamming triangle
inequality and the isometry of v. For example,
`|a+b|<=|c+v(a+b)|+|c|`, and
`|b|<=|b+va|+|a|`. This proves all three assertions. QED.

These constants are sharp within this class. Choosing a=b=N,c=0
gives (1.3), with total boundary weight 4q and intersection 2q
with every pair. Thus all three energies in (3.3) can vanish.

S8 already certifies canonical HM on a larger generator strip by
a different finite potential. The result here is a direct pairwise
inequality on the specified smaller strip. Neither statement
permits clipping an arbitrary generator into that strip.

### 3.1 No single-generator-column counterexample, anywhere

Let `h=x^d H`, `w=|H|`, `z=|uH|`. For every q,

```
|boundary(h)|=4w+z,       0<=z<=2w.
```

The first block has three distinct x columns with masks H,vH,vH;
the second has distinct columns d+3 and d with masks H,uH.
The only nonzero pair intersections are:

| d | I_AB | I_AC | I_BC |
|---:|---:|---:|---:|
| -3 | w | 0 | w |
| -2 | w | 2w | w |
| -1 | w | w | 0 |
| 0 | w+z | w | z |
| 1 | 0 | z | z |

All other d give zero intersection. The inequalities
`2I_pair<=4w+z` follow directly from `z<=2w`. Thus any original
counterexample must use at least two generator x columns, and
Theorem 3.1 says its unique generator cannot be confined to [-3,-1].

## 4. Exact counterfamilies to folding x modulo three

This section concerns the **quotient** `x^3=1`, not the original
cylinder. It explains why a small Fourier transfer graph cannot
simply certify the desired 2q pair minima.

Let q>=2 be dyadic, `Phi=1+t+t^2`, and `J=Phi^-1` in R_q. Explicitly,

```
q=3v+1: J has coefficients 1 at rows 0 or 2 modulo 3;
q=3v+2: J has coefficients 1 at rows 1 or 2 modulo 3.
```

These formulas include the cyclic endpoints. In the interior,
every sum of three successive coefficients is zero. At row zero
the cyclic sum is one; at row one it is zero. This proves
`Phi J=1`. Their weights are, respectively,

```
|J|=(2q+1)/3,       |J|=(2q-1)/3.                  (4.1)
```

In the x-three quotient define the explicit chain

```
z=(J+(N+J)(x+x^2), chi).                           (4.2)
```

Here chi in the second block is in row zero only. Since
`(x+x^2)chi=0` in this quotient and P becomes Phi,

```
D z=(x+x^2)N.
```

Both AB and AC have this folded syndrome. More strongly, z is
in each relevant folded boundary coset, with explicit generators

```
z+T_AB = boundary(N+chi J),
z+T_AC = boundary(xN+chi J).                        (4.3)
```

All three residues in each block remain physically distinct in
the quotient, so its exact weight is

```
|z|=2q-|J|+3
    =(4q+8)/3  if q=1 mod 3,
    =(4q+10)/3 if q=2 mod 3.                       (4.4)
```

For every dyadic q>=8, |J|>3, hence this is strictly below 2q.
At q=8, `J=t+t^2+t^4+t^5+t^7` and the weight is fourteen,
against the proposed pair value sixteen. This is a closed formula,
not an inference from a height sample.

BC fails in the quotient even more directly:

```
T_BC=boundary((1+x)N)       when x^3=1,
```

so its folded coset contains zero.

### 4.1 The original-column residual is not zero

Use the same explicit z from (4.2) without folding x. A direct
Laurent identity gives

```
D z+D T_AC
  =(1+x^3)[chi J+(x+x^2)N+x t^3],
D z+D T_AB
  =(1+x^3)[chi J+(x+x^2)N+x t^3+N].                (4.5)
```

For q>=8, the coefficient at actual x column five is N+J,
which is nonzero. These are therefore not valid representatives
of either original coset. The quotient counterfamilies cannot
be called counterexamples to the parent route or to canonical HM.
They expose precisely the spatial information discarded by this
particular shortcut, in addition to S8's separate x-six obstruction.

## 5. Residual and bounded verification

The remaining question is the original, unrestricted pair-leader
assertion for dyadic q>=4 and multi-column generators outside
the proved strip. No exact original counterfamily is supplied.
The three individual hole-column inequalities in (1.4) are useful
equivalent formulations, but are not newly solved lemmas.

The recovery width-load obstruction rules out assuming a uniformly
bounded local linear decoder. S8 rules out an unrestricted x-six
folding premise and does not establish width completeness. This
note uses neither retired shortcut and performs no excess-four
shell classification.

The companion checker verifies the formal three-column expansion,
the exact coefficient identities in the triangle certificates,
and one prescribed q=8 instance of (4.2)--(4.5), including the
explicit folded boundary generators and nonzero original residual.
It does not search supports, cosets, word layouts, heights, or
families, and it does not label the unrestricted route proved.

Run from the repository root:

```
PYTHONDONTWRITEBYTECODE=1 uv run --offline --no-sync --no-cache \
  --project experiments/bb_lab python \
  experiments/bb_lab/scripts/a42_s18_codex_pair_leaders.py
```

Only these two new assigned files are written. The checker enforces
55-second wall and 50-second CPU limits and checks peak RSS against
768 MiB. No solver, shared-file edit, git operation, old process,
memory edit, Lean work, or additional agent is involved.

At the parent's request to close the bounded pass, the new checker
is supplied without a diagnostic replay in this pass. The results
above are the displayed analytic proofs; no unrun check is reported
as having passed.
