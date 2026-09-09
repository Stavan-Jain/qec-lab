# A42 S8 — canonical discharge: exact local lemmas and a folding obstruction

2026-09-08. Codex bounded audit of the **original** Laurent pair. New
files only. No claim of a full canonical HM proof or an IBM distance
theorem is made here.

## 1. Disposition and quantifiers

The requested unrestricted canonical inequality remains **open**.
This pass banks the following narrower results:

1. A genuine counterexample to the proposed sufficient relaxation that
   identifies `x` modulo six. It occurs at the permitted binary height
   `q=512`: a folded boundary has inside weight `1207` and outside
   weight `381`, so `3 outside - inside = -64`. This is **not** a
   counterexample on the unbounded original cylinder. It rules out a
   discharge proof valid for all boundaries after that folding.
2. An analytic canonical-HM theorem for every row mask when the
   boundary generator has only columns `-3,-2`. A finite exact
   potential also certifies the larger generator strip `[-3,2]`, with
   no restriction on the height. The latter is a computer-assisted
   finite identity, not a width-completeness theorem.
3. Two spatial burst lemmas: a noncanonical monomial D-word has at
   most seven entries and at most three holes; a monomial g-only
   interval has at most nine rows and at most three holes. These
   statements retain the actual Laurent coordinate, not just its
   residue. They do not yet give the sharp charge across defects.
4. An explicit alternating-mask obstruction to taking the maximum of
   the raw syndrome bound and the full-support coset triangle bound.

The all-height leader uniqueness in `A42_S7_spatial_equality.md` is a
**sibling result**, not a discovery of this note. Its proof was read:
the spatial contradiction at `D->C` uses actual Laurent positions and
does not discard the five-cycles. Its conclusions for `k<=3` and the
general `4q+2` floor are retained.

### Definitions

Unless a section explicitly says otherwise, let `q=2^a`, `a>=0`, and
work in `R_q[x^±]`, where `R_q=F_2[t]/(t^q+1)`. All `x` supports are
finite. Set

```
P=x^3+t+t^2,        Q=1+(x+x^2)t^3,
D(f,g)=Pf+Qg,      B_q={(Qh,Ph): h in R_q[x^±]},
N=1+t+...+t^(q-1),
1_S=(N,(1+x)N).
```

Thus S contains the whole f-column zero and the whole g-columns zero
and one. In particular `|S|=3q` and `D(1_S)=N`. The underlying omega
support is the canonical support from S6; we do not vary the omega
cycle in this note.

For a boundary s put

```
I=|s intersect S|,   b=|s outside S|,   k=3q-I,
e=1_S+s,            delta=|e|-q=k+b-q.
```

Full canonical HM asks, equivalently,

```
3 d(1_T,B_q) >= |T|       for EVERY subset T of S;
3b >= I                  for EVERY compact boundary s;
k+3b >= 3q;
2k <= 3 delta.                                             (1.1)
```

Here the binary chain e has k sites in S, not I sites. Its physical
mixed lift has weight

```
3q+k+3b = 6q+3 delta-2k.
```

For completeness, the subset equivalence needs no minimality
assumption. Necessity follows by choosing `T=S intersect s`, since
`1_T+s` has exactly the outside support of s. Conversely, for any
T and s, put `J=S intersect supp(s)`. Then

```
|T| <= |J|+|T symmetric_difference J|
    <= 3|s outside S|+|T symmetric_difference J|
    <= 3|1_T+s|.
```

Minimizing over s proves sufficiency. No condition on the shape or
connectedness of s has entered.

## 2. Why the two global lower bounds do not combine by a maximum

For arbitrary binary row masks i,j,kappa, direct multiplication gives

```
D(i,j+x kappa)
  = ((t+t^2)i+j)
    + x (kappa+t^3 j)
    + x^2 t^3(j+kappa)
    + x^3(i+t^3 kappa).                                  (2.1)
```

Every physical binary coordinate has syndrome weight at most three.
Thus `3|e|>=|D(e)|`. The four columns in (2.1) are **distinct actual
columns** and must not be identified to use their separate contents.

For even q, let i be the alternating mask with ones in the even rows,
and take `T=(i,(1+x)N)`. Then

```
(t+t^2)i=N,        |T|=5q/2,
D(T)=x^3(N+i),     |D(T)|=q/2.                            (2.2)
```

The raw syndrome estimate gives only `d(T,B_q)>=q/6`. The banked
identity `d(1_S,B_q)=q` and `|S minus T|=q/2` give only `q/2` by the
triangle inequality. Their maximum is strictly less than the needed
`5q/6`.

Even the compressed checks from S6 do not rescue this example without
additional spatial information. They have

```
J_1(T)=Phi i=N+i,       J_w(T)=i,       Phi=1+t+t^2.
```

The chain `(i,0)` has these same compressed checks and weight `q/2`.
It is not asserted to have the actual syndrome (2.2); it generally
does not. This separates the missing information precisely: the
identity of the four syndrome columns, rather than only their total
weight or their two Fourier evaluations, still has to be used.

As the parent observed, the simpler mask `(0,(1+x)N)` has invisible
compressed checks but actual syndrome `(1+x^3)N`, of weight `2q`.
The raw bound handles that mask, but (2.2) shows why this does not
close arbitrary mixtures.

## 3. A decisive obstruction to folding x modulo six

### 3.1 The sufficient reduction is sound; its proposed premise is false

Give chain coordinates in S cost one and all other coordinates cost
three, and denote the resulting weighted Hamming norm by W. Then

```
W(1_S+s)=3q+3|s outside S|-|s intersect S|.
```

Folding `x` modulo six maps S injectively to its counterpart in the
quotient. A coordinate outside S can map into S, but its cost then
decreases from three to one. Parity summation cannot increase weighted
Hamming norm. Folding commutes with P and Q and maps boundaries to
boundaries. Therefore a lower bound `W(1_S+s)>=3q` for **all quotient
boundaries** would suffice for the original cylinder.

The universal quotient statement is false, including at binary
height. The following finite identity proves this, without drawing an
infinite law from samples.

### 3.2 Explicit certificate

Encode a polynomial modulo `x^6+1` by a six-bit integer; bit c is the
coefficient of `x^c`. Let `rho` rotate this mask one position to the
left, implementing multiplication by x. If H_j is the mask of the
j-th row of h, then the boundary row masks are exactly

```
F_j=H_j + rho(H_(j-3)) + rho^2(H_(j-3)),
G_j=rho^3(H_j) + H_(j-1) + H_(j-2),                       (3.1)
```

where additions are XOR and the indices are cyclic.

Consider first these seven periodically repeated rows:

| j | H_j | F_j | G_j |
|---|---:|---:|---:|
| 0 | 61 | 1 | 3 |
| 1 | 10 | 0 | 3 |
| 2 | 32 | 1 | 51 |
| 3 | 13 | 1 | 3 |
| 4 | 53 | 9 | 3 |
| 5 | 3 | 0 | 32 |
| 6 | 47 | 1 | 11 |

Each entry follows by substitution in (3.1). The inside bits are bit
zero of F and bits zero and one of G. This table has inside weight
17 and outside weight 5, hence energy `3 outside-inside=-2`.

To obtain binary height, take L whole copies of the H-word followed
by at least three zero rows. The only changes from the repeated
boundary table are the first three rows and the first three padding
rows. Their new `(F,G)` pairs are respectively

```
start: (61,47), (10,44), (32,51);
end:   (60,44), (10,47), (33,0).
```

All later padding rows vanish. Comparing these six rows to the
periodic table gives the exact formulas, for every `L>=1`,

```
inside=17L,       outside=5L+26,
energy=78-2L.                                            (3.2)
```

The proof of (3.2) is locality of (3.1) and the displayed seam table,
not numerical extrapolation. Choose `L=71` and fifteen padding rows:
the height is `7*71+15=512`, and

```
inside=1207,      outside=381,      energy=-64.
```

Equivalently the folded canonical error has `k=329`, `delta=198`,
and violates `2k<=3 delta` by 64. Thus even retaining all six x
residues, all local row constraints, and the power-of-two height is
not sufficient once actual Laurent positions are discarded.

**Non-counterexample warning.** Regard the same masks as polynomials
in actual columns zero through five without x reduction. Their
boundary has inside weight 783 and outside weight 2083, so energy
5466, not a violation. The quotient counterexample proves neither
failure of original canonical HM nor failure of the IBM d-law.
It also does not rule out spatial factor potentials retaining X=x^3.

The tiny certificate script verifies the table, seam, and both
interpretations using integer XOR. The separate discovery script
uses a fixed six-column graph; its negative 42-row cycle reduces to
the seven-row boundary table above. No inference from failure of a
Bellman iteration is used as the counterexample proof.

## 4. Restricted positive results on actual, unfolded columns

### 4.1 An analytic two-column theorem

**Theorem 4.1.** For every positive height q, and arbitrary row masks
`a,b in R_q`, the actual boundary generated by

```
h=x^-3 a+x^-2 b
```

satisfies canonical HM: `3|s outside S|>=|s intersect S|`.
Its energy is zero only for `a=b=0` or `a=b=N`.

Proof. Put `u=t+t^2`. The nonzero possible boundary columns are

```
f_-3=a,  f_-2=b+t^3 a,  f_-1=t^3(a+b),  f_0=t^3 b;
g_-3=u a,  g_-2=u b,    g_0=a,          g_1=b.
```

Because multiplication by t is a permutation of row coordinates,

```
I=|a|+2|b|,
outside=|a|+|b+t^3a|+|a+b|+|u a|+|u b|.
```

The Hamming triangle inequality `|b|<=|a|+|a+b|` gives

```
I <= 3|a|+2|a+b| <= 3 outside.
```

Equality forces `a+b=0`, `b+t^3a=0`, and `u a=u b=0`.
Since t is a unit, the last equations force a and b to be constant
row masks. Together these conditions give precisely the two cases
claimed, and direct substitution verifies both. This proof has no
width assumption on a coset leader: it explicitly restricts h. QED.

### 4.2 A larger, computer-certified strip theorem

The diagnostic `a42_s8_codex_strip.py` certifies the same inequality
for every height and every h supported in actual columns `[-3,2]`.
It does not identify the two ends of this strip. Here is the complete
finite certificate specification and its logical interpretation.

Encode an h row by six bits, now with bit zero at actual column -3.
For states `(b,c,d)` and a new row a, all in `{0,...,63}`, define

```
alpha(a,d)=3 popcount(a XOR (d<<1) XOR (d<<2))
           -4 bit_3(a XOR (d<<1) XOR (d<<2));
beta(a,b,c)=3 popcount((a<<3) XOR b XOR c)
           -4 popcount((((a<<3) XOR b XOR c)>>3) AND 3).
```

There is **no masking of high bits** in these expressions. Thus all
outside syndrome-block columns are charged. The edge
`(b,c,d)->(a,b,c)` has cost `alpha(a,d)+beta(a,b,c)`, exactly the
row contribution to boundary energy.

Define F_0 to be zero on the `64^3` states and, four times, set

```
F_(r+1)(a,b,c)
  = min(F_r(a,b,c),
        beta(a,b,c)+min_d(F_r(b,c,d)+alpha(a,d))).          (4.1)
```

The finite identity checked is `F_5=F_4`; equivalently, the script
checks every one of the `64^4` edge inequalities

```
alpha(a,d)+beta(a,b,c)+F_4(b,c,d)-F_4(a,b,c) >= 0.          (4.2)
```

The minimum slack is zero. The histogram of F_4 is

| value | number of states |
|---:|---:|
| -4 | 1 |
| -3 | 130 |
| -2 | 388 |
| -1 | 391 |
| 0 | 261234 |

These are exact small integer operations, not floating-point
optimization. The intermediate four-step path costs and edge costs
are far within signed 16-bit range. Summing (4.2) over any cyclic
sequence of h rows cancels the potential and proves nonnegative
energy for that strip, for every height, including heights one and
two where row indices coincide. This is a finite computer-assisted
certificate, not a sampling statement about selected heights.

There is **no proof** that an arbitrary compact h can be replaced
by one in this strip while not increasing energy. Merely clipping h
keeps its inside output but can destroy outside cancellations. The
strip result must not be advertised as the unrestricted target.

## 5. Spatial clean-segment lemmas for quantitative stability

The S6/S7 graph gives, for an actual canonical syndrome solution e,

```
delta = L + sigma,
L = physical weight minus x-mod-three aggregated weight,
sigma = sum of the 128-edge graph's nonnegative slacks.    (5.1)
```

Here L is a nonnegative even integer. On a row with no aggregation
loss, each occupied residue really contains just one physical site.
On a zero-slack `D->D` row this means `g_j=0` and
`f_j=X^(n_j)`, where `X=x^3` remains an unbounded Laurent variable.
On the separate zero-slack two-cycle, f is empty and g is a monomial.

These implications are local only on rows with zero aggregation
loss. In particular an empty aggregated residue is not generally
an empty actual residue.

### 5.1 D-only monomial words

**Lemma 5.1.** Let `M_0,...,M_(r-1)` be Laurent monomials in X
satisfying

```
X M_j + M_(j-1) + M_(j-2) = 1      for 2<=j<r.            (5.2)
```

Either all entries are `X^-1`, or `r<=7`. In either case at most
three entries equal 1. For `r<=1`, only the second alternative is
needed.

Proof. Write consecutive exponents as `(a,b)`. If `a=b`, the next
exponent must be -1. If `a!=b`, a next monomial exists only when
exactly one of a,b is zero; its exponent is the other exponent minus
one. This follows because distinct Laurent monomials are linearly
independent.

The pair `(-1,-1)` persists forever. An unequal pair with neither
entry zero has no successor. An equal noncanonical pair has at most
two further entries; the longest such case is `0,0,-1,-2`.
The only longer exceptional continuations are prefixes of

```
0,1,0,0,-1,-2;
1,0,0,-1,-2;
2,0,1,0,0,-1,-2.
```

Indeed, from `(0,b)` a further successor after `b-1` requires `b=1`;
from `(a,0)` further successors beyond `a-1,a-2` require `a=1` or
`a=2`. These cases give exactly the displayed words. Each has at
most three zero exponents and length at most seven. This exhausts
all initial pairs and proves the claim. QED.

In applying the lemma to actual rows, (5.2) requires both `g_j=0`
and `g_(j-3)=0`. A D-only interval supplies the latter automatically
only after its first three rows. The endpoint qualification is
essential; one cannot silently charge endpoint holes as if they
were an interior recurrence.

### 5.2 g-only monomial words

**Lemma 5.2.** Suppose that in a consecutive interval of r actual
rows, f vanishes, g is one Laurent monomial per row, and the full
canonical syndrome equation holds. Then `r<=9`, and at most three
of these g monomials can occupy a site in S.

Proof. On every row j at least three rows into the interval the
syndrome equation reduces to

```
g_j+(x+x^2)g_(j-3)=1.                                    (5.3)
```

If `g_(j-3)=x^m` and `g_j=x^n`, the only possibilities are

```
(m,n)=(-2,-1) or (-1,1).
```

The two shifted monomials on the left are distinct. One must cancel
`x^n` and the other must equal 1, proving this exhaustive two-case
list. Along each row-residue class modulo three, the only possible
nontrivial chain is a subword of

```
x^-2 -> x^-1 -> x.
```

Thus each of the three subwords has length at most three. In a
subword of length at least two, at most its last entry can belong
to S, namely x in g-column one. A subword of length one also has at
most one site in S. Summing gives the two bounds. QED.

This is a local quantitative replacement for the global x-extrema
exclusion of a wholly g-only cycle. It does not rely on graph-only
uniqueness or on the two-cycle's period failing to divide q.

### 5.3 The parent's six-edge obstruction, audited

The exact clean word `A->B->C->A->B->C->A` occupies six rows.
By the S7 residue table, rows zero and three contain only a
residue-zero f monomial; rows one and four have only g in residues
one and two; rows two and five are empty. The actual row equation is

```
x^3 f_j+f_(j-1)+f_(j-2)+g_j+(x+x^2)g_(j-3)=1.
```

At row three it gives `x^3 f_3=1`; at row five it gives `f_3=1`.
They contradict one another in the unbounded coordinate. Thus the
parent's obstruction is valid under the stated **clean actual-row**
hypotheses. It is not valid merely from aggregated row masks.

## 6. The remaining sharp step, without relabeling it as a theorem

The local results explain why a long defect-free excursion cannot
keep producing holes. They do **not** prove that distinct short
excursions can be charged without overlap to L and sigma in (5.1).
Several issues remain: endpoints with nonzero delayed g terms;
zero-slack transient edges; the E/F excursions in the six-state
component; and physical cancelling pairs concealed by aggregation.

The missing quantitative lemma is a spatial packing/discharging
statement assigning all k actual holes a total charge at most
`3(L+sigma)/2`. Its conclusion is precisely the still-open inequality,
not a new proof. A usable local version must specify the windows,
endpoint terms, and capacities so that two overlapping burst
certificates cannot spend the same defect twice.

One very explicit **different sufficient target**, motivated by §4,
is the following strip-completeness statement:

> For each q and every compact h of negative canonical boundary
> energy, there is an h' supported in `[-3,2]` whose energy is no
> larger.

The finite strip certificate would immediately contradict h'. No
such completeness statement is proved or assumed here. The S6
inverse impulse has width growing with q, so an argument imposing
a strip on arbitrary solutions cannot be justified merely by
locality of the original generators.

No further growing-width scan was built. In particular, the
counterexample of §3 is a reason to retain X in a proposed factor
potential, not a reason to increase a transfer width until a sample
passes.

## 7. Verification scope

New scripts in this pass:

- `a42_s8_codex_discharge.py`: the fixed x-mod-six potential audit
  and extraction of a negative cycle; no claimed lower-bound proof.
- `a42_s8_codex_certificates.py`: direct verification of the explicit
  folded counterexample, its seam, and its non-counterexample lift.
- `a42_s8_codex_strip.py`: exact finite strip-potential certificate.

The explicit table certificate takes under a millisecond in its
measured body; the strip identity about 0.03 seconds; the negative
cycle extraction about 0.82 seconds. Every diagnostic invocation
stayed below one second and below 70 MiB measured peak RSS; aggregate
diagnostic time in this pass is below ten seconds. Scripts install
55-second wall alarms and 50-second CPU limits. No SAT solver,
code-family census, or height-dependent transfer width was used.
Earlier failed assertions during construction were not accepted as
certificates; the final counterexample uses q=512, not the initially
insufficient q=128 padding choice.

Run each with the lab environment, for example:

```
uv run --no-cache --offline --frozen --no-sync --project experiments/bb_lab \
  python -B experiments/bb_lab/scripts/a42_s8_codex_certificates.py
```

## 8. Independent review of the parent's augmented-gross proof

After the primary audit, I read §§2–6 of
`GROSS_SURFACE_PRODUCT.md`. Conditional on the stated banked gross
input `d_X=d_Z=12`, I found no gap in the exact-distance argument:

- **Redundancies:** retaining `H_0(C)` and `H_2(C)` causes no extra
  middle homology. R has only degree-one homology and I only
  degree-zero homology, so degree two selects exactly `H_1(C)`.
- **First interval:** the adjacent cycle equations identify all t
  base homology classes. If the common class is zero, subtracting
  the h_j eliminates the u_j. The successive sums of residual
  v_j in `ker d_2` eliminate the remainder without restoring u.
  This handles dependent gross checks, rather than assuming them
  away. Every nonzero class therefore pays t base distances.
- **Second interval:** vertex augmentation is a chain-homotopy
  equivalence and contracts physical Hamming weight, even on chains
  with additional edge components. The first-vertex inclusion gives
  the matching weight-preserving witness. Homology classification
  alone is not being substituted for this metric argument.
- **Dual:** reversing the full length-four tensor complex fixes
  degree two, sends the factors to `C^vee,I,R`, and swapping the
  last two factors is an actual coordinate permutation. Thus the
  same proof applies to Z, without an unqualified quantum-product
  distance formula or a local-Clifford metric claim.
- **LDPC:** the three X check bounds are `8,8,7`; the three Z bounds
  are `7,8,8`. The listed qubit participation bounds follow from
  gross row/column degrees and the path degrees. The connectivity
  argument joins adjacent copies in both interval directions and
  attaches every auxiliary qubit type to the gross-qubit array.

This is a review of the self-contained algebra, not an independent
recertification of the base distance or a literature/novelty claim.
No product diagnostics were rerun. The result remains an augmented
surface construction, not the original scalar-BB goal and not a
resolution of canonical HM.
