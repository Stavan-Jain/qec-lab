# A44-S1 restart: IBM upper witnesses, exact obstructions, and a doubling recurrence

Date: 2026-09-08. Fresh Codex research agent; branch
`claude/a40-tour-de-gross-column`, PR #27. This is a separate research note;
it does not modify or supersede the banked notes.

## 0. Outcome and scope

The IBM equality is **not proved here**, even on an infinite subsequence.
The useful original-pair results are:

1. A40's L12 upper witness and dual satisfy explicit finite Laurent
   identities. In fact they prove `d(6q,6r) <= 12r` for all `q,r>=1`,
   including the square column. They do not supply its six-qubit discount.
2. At the first unresolved square, `(18,18)`, **every nonzero one-block
   X-cycle has weight at least 54**, and a one-block logical attains 54.
   Thus a weight-30 witness must use both blocks. This is an algebraic
   proof of the obstruction, not a repetition of the band census.
3. The uncorrected single helix has syndrome weight **12r** for every
   `r>1`. It is not a stabilizer: it is not a cycle.
4. For `r=2^a` there is an explicit two-block logical `W_r`, a doubling
   recurrence, and a coefficient-extraction dual pairing proving
   nontriviality for every `a`. Its exact weight is `6*3^a`, so the raw
   recurrence misses the desired value from `r=4` onward. The remaining
   candidate is a precisely stated compression lemma in this particular
   homology class, reduced to an odd-degree diagonal quotient (§5).
5. An independent review of the parent's separable comparison theorem
   finds no mathematical gap (§6). That theorem is a standard product
   specialization, substantially less efficient, and does not solve IBM.

Inputs audited: `a44_common.py`, `a44_s0_patterns.py`, both
`data/a44/s0_patterns.{json,log}`, A43 §5, A40 §7.1, the L12/dual JSON
certificates, and A42 §§1.1–1.3. No old process was resumed, signalled, or
relaunched. No SAT, census, Lean, library editing, or git mutation was used.

**Dimension exception retained.** A42 §1.3 gives the original IBM code
`k=12` for all `b=1`, and `k=12+14*[127 divides r]` for `b=0`. This note
does not reprove that spectral calculation. The extra dimensions do not
alone refute the conjectured distance. A43's global four-zero/dimension-six
description is not used: the F128 orbit exists. The corrected Liang et al.
identifier supplied by the parent is `2503.03827`, not `2503.04699`.

## 1. Conventions

Work over

`R_(ell,m)=F2[x^{+/-1},y^{+/-1}]/(x^ell+1,y^m+1)`.

The fixed original checks are

`A=1+y+x^3 y^-1`, `B=1+x+x^-1 y^-3`,

and put

`alpha=bar(A)=1+y^-1+x^-3 y`,
`beta=bar(B)=1+x^-1+x y^3`.

An X-vector is a pair `(f,g)`, with physical weight `|f|+|g|`.
It is a cycle iff `beta f+alpha g=0`; stabilizers are
`h(alpha,beta)`. A Z-vector `(u0,u1)` is a cycle iff
`A u0+B u1=0`. The pairing is coefficientwise binary dot product.
In particular an odd pairing with a Z-cycle proves X-nontriviality.

These follow directly from the lab matrices
`H_X=[M_A|M_B]`, `H_Z=[M_B^T|M_A^T]` and agree with S0's checked
conventions. Squaring a Laurent polynomial means Frobenius squaring,
including doubling its exponents.

The underlying original codes are connected at every positive member
parameter. In the X-check incidence graph, constant terms link the two
qubit blocks at each cell. The `1+y` terms in A and `1+x` terms in B link
neighboring cells in both coordinate directions. Thus every cell and
both qubit blocks lie in one component. The growing patterns below live
in this fixed connected code; no generator-polynomial dilation or
disjoint-copy construction is being substituted.

## 2. Audited L12 theorem, including a square-column upper bound

Write the banked L12 pair as

```
F = x(1+y^2) + x^2(y+y^3) + x^3(y+y^5)
G = x^3(y+y^5) + x^4(1+y+y^2+y^5).
```

Its exact plane identity is

```
beta F + alpha G = (1+y^6) E,
E = 1+x+x^3+x^4(y^-1+y^2).                         (2.1)
```

The banked dual motif is

```
U0 = y+y^2+y^3 + x(y^2+y^3) + x^2 y^2
     + x^3(1+y^2) + x^4 y
U1 = y^4 + x y^2 + x^2 y^3.
```

Its plane identity is

```
A U0 + B U1 = (1+x^6) J,
J = x^-1 y + y^-1 + y + x.                        (2.2)
```

Both identities are verified simply by expanding their finite supports.
The diagnostic independently checks them before any torus calculation.

**Theorem 2.1.** For every `q,r>=1`, on `(ell,m)=(6q,6r)`, let

```
V = (sum_(j=0)^(r-1) y^(6j)) (F,G),
U = (sum_(i=0)^(q-1) x^(6i)) (U0,U1).
```

Then V is a nontrivial X-cycle, `|V|=12r`, U is a Z-cycle, and
`<U,V>=1`. Consequently `d(6q,6r)<=12r`.

*Proof.* Multiplying (2.1) by the y-sum telescopes its syndrome to
`(1+y^(6r))E=0`. Equation (2.2) does the same for U in x. The twelve
cells of `(F,G)` occupy rows 0 through 5 and columns 1 through 4, so
their y-translates are disjoint even at `q=1`. The dual motif occupies
columns 0 through 4 and rows 0 through 4. Only the zero translate in
each direction can meet the other motif. Their sole common physical
qubit is block 0 at `(x,y)=(1,2)`. Thus their pairing is one. ∎

Taking `q=r+1` verifies the banked `UB(r,1)` with no extrapolation from
samples. Taking `q=r` supplies the useful but non-attaining bound `12r`
on the square column. In particular `(18,18)` has `d<=36` independently
of whether its conjectured value 30 is correct.

## 3. A sharp all-one-block obstruction at r=3

Consider first `(f,0)` on `(18,18)`. The equation `beta f=0` separates
the three residue classes of y modulo 3. Within any one class, write

`f=sum_(i=0)^5 y^(c-3i) p_i(x)`.

The coefficient equation gives

`p_(i+1)=T p_i`, where `T=x^-1+x^-2` in `F2[x]/(x^18+1)`.

Closure is equivalent, after multiplication by a unit, to

```
Q p_0=0,
Q=(1+x)^6+x^12
 =(x^2+x+1)^2 (x^4+x^3+1)^2.                      (3.1)
```

**Lemma 3.1.** `gcd(Q,x^18+1)=(x^2+x+1)^2`.

*Proof.* Factor

`x^18+1=(x+1)^2 (x^2+x+1)^2 (x^6+x^3+1)^2`.

Put `e=x^4+x^3+1`. We have `e(1)=1`. Polynomial division gives

```
x^6+x^3+1 = (x^2+x+1)e + (x^2+x),
e = x^2(x^2+x)+1.
```

Thus e is coprime to the other two factors, proving the claim. ∎

The elementary annihilator rule in `F2[x]/M` now gives all seeds:

```
p_0 = (1+x^6+x^12) h,
h=(1+x^2)q,  deg(q)<4.                            (3.2)
```

Indeed `ann(Q)` is the ideal generated by `M/gcd(M,Q)`, and that
generator equals `(1+x^2)(1+x^6+x^12)`. Each h has degree below 6.
Thus every seed is three exact x-copies of a length-six seed. Multiplying
by T respects this embedding, since the copy sum is invariant under
`x^6`. On the length-six seed space, `T^2 h=h`.

Here is the complete length-six cost argument. The space generated by
`1+x^2` consists exactly of words having even parity separately on
positions `{0,2,4}` and `{1,3,5}`. A nonzero h therefore has weight 2
or 4. When it has weight 2, its two positions differ by 2 or 4, and
expanding `T h` gives four distinct positions. When it has weight 4,
`T h` is nonzero because `T^2=1`, so it has weight at least 2.
Consequently

`|h|+|T h| >= 6`, with equality at `h=1+x^2`.

There are three copies in x and three repetitions of this two-row
cycle in y. Every nonzero residue-class band therefore costs at least
`3*3*6=54`. Different y residue classes have disjoint physical support,
so adding them cannot lower this bound. The analogous column recurrence
for `(0,g)`, using alpha, gives the same bound with x and y interchanged
and the appropriate cyclic order of columns.

**Theorem 3.2.** Every nonzero one-block X-cycle on `(18,18)` has
weight at least 54. The minimum among nontrivial one-block X-logicals
is exactly 54.

*Attainment and pairing.* Take

`p_0=(1+x^2)(1+x^6+x^12)`, `p_i=T^i p_0`, `c=1`, and `g=0`.

The constructed cycle has weight 54. Pair it with U from Theorem 2.1
at `q=r=3`. In the dual's row range 0 through 4, the band's occupied
rows are 1 and 4. U0 has no row-4 cells; in row 1 its cells have
x-coordinates `0,4` modulo 6. The band has coordinates `0,2` modulo 6
there. The three x-periods therefore contribute three intersections,
an odd pairing. ∎

This obstruction covers **all** one-block cycles, not only the stored
choice of band phase or minimizer. It does not say the code distance
is 54: Theorem 2.1 already gives 36. It says a proposed weight-30
attainer has to be a genuinely two-block object.

## 4. Why the simple helix is not a cycle

On the square `N=6r` set `tau=xy` and

```
g=x sum_(j=0)^(3r-1) tau^(2j),
f=x^-3 y g.
```

Then `(1+tau^2)g=0`. Direct cancellation gives

`beta f+alpha g=(1+y^-1+x^-4 y+x^-2 y^4)g`.

Using `tau^2 g=g` reduces this exactly to

`beta f+alpha g=(1+y^-1)(1+x^-6)g`.                 (4.1)

For `r=1`, `x^6=1`, and the residual vanishes. For `r>1`, the four
translates on the right are four distinct cosets of the subgroup
`H=<(2,2)>` of `Z_N^2`, each of size `3r`. The translates by `y^-1`
are disjoint from the other two by y-parity. The translates differing
by `x^-6` could coincide only if `(6,0)` belonged to H, which would
require `6=0 mod N`. Hence all four are disjoint for `N>6`, proving

`|syndrome(f,g)|=12r` for every `r>1`.

The S0 script's introductory claim that the helix is a cycle and becomes
trivial contradicts both this proof and its actual JSON/log output.
The executable output is the correct one. Likewise the original
length-six one-block stacking defect is the nonzero polynomial
`(1+x^2+x^4)(1+x^2)=1+x^6` above period six.

## 5. A nontrivial Frobenius recurrence and its exact remaining lemma

### 5.1 An explicit infinite original-pair logical family

Let `r=2^a`, `a>=0`, keep the original square code `(6r,6r)`, and define

```
S_r=1+(xy)^(2r)+(xy)^(4r),
W_r=(S_r x^r beta^(r-1), S_r y^r alpha^(r-1)).       (5.1)
```

At `r=1` this is the banked weight-six two-block helix, with a choice
of translation. The doubling recurrence is completely explicit:

```
D(f,g)=(beta f(x^2,y^2), alpha g(x^2,y^2)),
W_(2r)=D(W_r).                                    (5.2)
```

**Theorem 5.1.** For every `a>=0`, W_r is a nontrivial X-logical of
the original IBM square code, with exact physical weight `6*3^a`.

*Cycle proof.* If `(f,g)` is a cycle at periods `(N,N)`, its syndrome
is in the period ideal. The syndrome of D(f,g) is
`(beta f+alpha g)^2`, which is in the doubled period ideal. Thus D
takes cycles to cycles. The seed is a cycle by direct expansion.
For a boundary it gives

`D(h alpha,h beta)=alpha beta h^2 (alpha,beta)`,

so it also respects logical equivalence.

*Exact weight proof.* Use coordinates `z=x`, `tau=xy` and the quotient

```
Q_r=F2[z^{+/-1},tau^{+/-1}]/(z^(6r)+1,tau^(2r)+1).
alpha=1+z tau^-1+z^-4 tau,
beta=1+z^-1+z^-2 tau^3.                           (5.3)
```

The motif before multiplication by S_r is
`(z^r beta^(r-1), z^-r tau^r alpha^(r-1))` in Q_r.
Since `r-1` has a ones digit in each of its a binary places, a term
of `beta^(r-1)` is indexed by two disjoint-bit integers `u,v<r`:

`z^(-u-2v) tau^(3v)`.

There are `3^a` choices. Distinct v cannot have equal `3v mod 2r`,
since 3 is invertible modulo `2r`; then distinct u cannot alias modulo
`6r`. Thus all terms survive. Terms of `alpha^(r-1)` are

`z^(u-4v) tau^(-u+v)`.

Their z-exponents lie in an interval of width `5(r-1)<6r`, and their
tau-exponents in one of width `2(r-1)<2r`. Equality in Q_r therefore
implies equality as integer exponent pairs, hence the same u,v.
They also contribute exactly `3^a` terms. Pullback to the square
replicates each motif cell three times, proving weight `6*3^a`.

*Uniform nontriviality with a specified dual pairing.* Let
`F4=F2(omega)`, `omega^2+omega+1=0`. There is a unique formal series
`z(v)=omega+v+omega^2 v^2+O(v^3)` solving

`1+z/(1+v)+(1+v)z^-4=0`.                          (5.4)

This assertion requires no root search: equivalently

`z=(1+v)+(1+v)^2 z^-4`.

Iteration starting at omega determines it coefficient by coefficient.
The difference of two right-hand sides has v-adic order four times
that of the difference of the inputs, so any finite truncation is
obtained after finitely many iterations. Alternatively the derivative
with respect to z in (5.4) is the unit `(1+v)^-1`.

Substitute `x=z(v)`, `y=(1+v)/z(v)` into
`T_r=F4[v]/(v^(2r))`. These substitutions respect `x^(6r)=y^(6r)=1`:
their constant terms have order three, and `2r` is a power of two.
They also give `tau^(2r)=1`, `S_r=1`, and `alpha=0`.
Expansion through degree two gives

`beta(z(v),(1+v)/z(v))=v^2+O(v^3)`.                (5.5)

For explicit verification of its coefficient, put `z=omega(1+u)`.
Equation (5.4) gives `u=omega^2 v+omega v^2+O(v^3)`; substituting this
in `1+omega^2(1+u)^-1+omega(1+v)^3(1+u)^-2`
gives (5.5).

Define the F2-linear functional on the physical group ring

```
lambda_r(P)=Tr_(F4/F2)(omega^(1-r) [v^(2r-2)]
                       P(z(v),(1+v)/z(v))).       (5.6)
```

It annihilates `alpha h` for every h. On the first component of W_r,
the extracted leading coefficient is `omega^r`, so its value is
`Tr(omega)=1`. Concretely the Z-vector with coefficients

`U_(r,0)[i,j]=lambda_r(x^i y^j)`, `U_(r,1)[i,j]=0`

annihilates every X-check and pairs oddly with W_r. This proves
nontriviality for every a, rather than assuming a cycle lift preserves
homology. ∎

The raw weights are `6,18,54,162,...` at `r=1,2,4,8,...`; the target
weights are `6,18,42,90,...`. This construction is therefore an exact
infinite **witness-weight** theorem, not an exact-distance theorem.
For `r>=4` the simpler L12 witness already has smaller weight.

### 5.2 The explicit compression candidate; nothing conditional is banked

The proposed remaining lemma is:

> **C(a), currently unproved:** in Q_r of (5.3), for `r=2^a`, there
> exists `h_a` such that
>
> `|(z^r beta^(r-1),z^-r tau^r alpha^(r-1))`
> `       +h_a(alpha,beta)| <= 4r-2`.

If C(a) holds, pullback gives a nontrivial original-pair witness of
weight at most `12r-6`. Nontriviality is already supplied by (5.6),
and boundaries cannot change it. An independent matching lower bound
would still be necessary for exact distance. C(a) concerns a specified
homology class on a specified infinite subsequence, rather than all
logical classes or an unlabelled minimum.

It also has a concrete recursive form. If a corrected motif at r is

`V_r=W_r^motif+h_a(alpha,beta)`,

the next candidate is

```
V_(2r)=(beta V_(r,0)(z^2,tau^2),
        alpha V_(r,1)(z^2,tau^2)) + k_a(alpha,beta),
h_(a+1)=alpha beta h_a(z^2,tau^2)+k_a.             (5.7)
```

Here `h_0=h_1=0` already works. What is missing is a construction of
the correction `k_a` satisfying C(a+1); (5.7) itself is an exact cycle
and homology identity, not a proved weight inequality. There is no
claim that an arbitrary small stencil suffices or that low-r evidence
implies this lemma.

**Smallest new instance.** At `r=4`, work on just
`Q_4=F2[z,tau]/(z^24+1,tau^8+1)` and seek a weight-at-most-14
representative of

`w=(z^4 beta^3, z^-4 tau^4 alpha^3)`.

Its current weight is 18. The desired assertion is exactly

`exists h: |w+h(alpha,beta)|<=14`.                 (5.8)

This quotient has 384 physical qubits, rather than the square's 1152.
The degree-three pullback is injective on homology: pushforward after
pullback multiplies chains by three, which is one over F2. Thus the
quotient is an honest periodic-witness reduction, not disconnected
copies. Failure in this particular class would not refute IBM's law.

There is a useful obstruction even to (5.8). Direct finite support
arithmetic gives maximum overlap two between w and any single check,
and maximum overlap one between distinct checks in Q_4. The script
prints the two nine-term motif supports and verifies these figures by
checking translations only, with no logical enumeration. For a sum
of j distinct checks, j<=3, its support has size at least
`6j-2*binom(j,2)` and its intersection with w has size at most `2j`.
Adding it therefore changes the weight by at least `j(3-j)>=0`.
Consequently **no correction consisting of at most three quotient
checks can improve w at all**. A successful compression must be
collective. This explains why an individual-check reduction stalls;
it does not give a distance lower bound.

## 6. Independent review of the parent's exact comparison theorem

Reviewed `IBM_EXACT_DISTANCE_PROGRAM.md` §3 as supplied during this
restart. No change to that file was made, and its computation was not
rerun. The mathematical claim reviewed is

`A_sep=1+x+x^2`, `B_sep=1+y+y^3` on `(3s,7u)`,
with `[[42su,12,min(2s,4u)]]`, for every `s,u>=1`.

**Contractions.** The lower-bound argument detects every non-boundary
cycle. If all contractions of the first block vanish, that block lies
in `E tensor im Q` and can be removed by a boundary. The remaining
second block is in `V tensor ker Q`; its vanishing contractions put
it in `im P tensor ker Q`, which is again a boundary. A nonzero
contraction has weight at most that of its original block because
each surviving coordinate needs support in a distinct coordinate row.
No unstated weight-preserving change of basis is needed for this step.

**Blocks and transposes.** For X-cycles take `P=M_A^T`, `Q=M_B^T`.
Swapping physical blocks changes `(f,g)` into `(g,f)` and changes the
lab boundary `(P h,Q h)` into `(Q h,P h)`, exactly the displayed
product complex. This is a coordinate permutation. For Z-cycles the
same reasoning uses the untransposed maps. Coordinate reversal
conjugates each cyclic shift to its inverse, proving equality of the
classical transpose distances used in the Z bound.

**Every parameter, including even ones.** The order-two recurrence
has three nonzero initial states, all with period three and two ones.
For the order-three recurrence, irreducibility of `z^3+z+1` identifies
its state update with multiplication by a nonidentity element of F8*.
That group has prime order seven, so all seven nonzero states form
one orbit; four have first bit one. Thus every nonzero word has
period seven and four ones. Reversibility proves these describe all
solutions. Repeating to lengths `3s,7u` introduces no additional
initial states, even when s or u is even. The transpose statement
follows by reversal, without assuming semisimplicity.

**Attainment and connectivity.** A coordinate vector surviving in a
nonzero cokernel always exists. Tensoring it with a minimum kernel
word gives the claimed physical weight, and an annihilating
functional proves nontriviality. The `1+x` and `1+y` check terms,
together with the cross-block constant terms, connect the whole
lattice. Each block has three distinct positions even at `s=u=1`,
so every check has weight six. Both periods grow in the `s=2t,u=t`
subfamily; fixed-distance disconnected copies are not concealed.

**Scope and efficiency.** This is a standard hypergraph-product
specialization with separate cyclic factors. It changes IBM's
two-variable coupling and its periods; it is not an original IBM
subfamily or a novelty claim. At `[[84t^2,12,4t]]`,
`k d^2/n=16/7`; compared with the conjectured IBM b=1 asymptote 24,
the qubit overhead at fixed large distance and k=12 is 10.5 times
larger. For complete precision, any comparison calling 24 simply
the "IBM asymptote" should identify b=1, or exclude the b=0
order-127 subsequence. The contraction proof itself has no such
exception and no gap was found.

**S-opt review: correct within the stated separable subclass.** For
two square circulant factors with positive kernel dimensions p,q, the
same product calculation gives `k=2pq`. Thus `k=12` requires `pq=6`.
A one-dimensional binary cyclic kernel would be fixed by the shift:
its unique nonzero vector must be the all-ones word. An odd-weight
check does not annihilate that word. Equivalently, a degree-one gcd
would require the root 1, whereas a trinomial takes value one there.
Hence neither dimension can be one, leaving `{p,q}={2,3}`.

For a binary `[n_i,p_i,d_i]` kernel, each active coordinate occurs in
exactly `2^(p_i-1)` of its `2^p_i-1` nonzero words; inactive coordinates
only strengthen the resulting inequality. Summing word weights yields
`n_i >= (2^p_i-1)d_i/2^(p_i-1)`. Therefore, after ordering the factors,
`n_1>=3d_1/2` and `n_2>=7d_2/4`. Circulant reversal supplies equal
transpose distances, so the quantum distance is `d=min(d_1,d_2)`.
Consequently `n=2n_1n_2>=21d_1d_2/4>=21d^2/4`. The parent's choice
`n_1=6t,n_2=7t,d_1=d_2=4t` saturates both inequalities. This is a
valid optimality restriction for the specified separable square-map
construction with k=12; it does not bound coupled IBM polynomials,
nonseparable BB codes, or product constructions with different
homological dimensions or nonsquare factor maps. No further
computation is needed for this review.

**Gross/Shor benchmark (§3.8) review: correct, conditional only on the
stated banked gross distances `d_X=d_Z=12`.** The inner Z-check ranks
are `t(t-1)` and the adjacent-full-row X-check rank is `t-1`, including
t=1. Full-row X patterns modulo adjacent row pairs are classified by
the parity of their occupied rows; Z patterns modulo within-row pairs
are classified by their common row parity. These descriptions give
both pure distances t for every t, with no even-t exception. A
concatenated normalizer restricts to an inner normalizer on each
block; a nontrivial outer effective Pauli occupies at least 12 blocks.
The block supports are disjoint, giving the lower bound 12t. The
banked pure outer X and Z attainers lift using inner representatives
of weight t, proving equality in both sectors.

The intrinsic non-LDPC argument is also valid for this stabilizer
code without additional qubits: its stabilizer quotient by the inner
stabilizers contains the nonzero encoded gross constraints. Any
generating set must therefore contain an element with nonidentity
outer image, which has weight at least t in some inner block. Thus
there is no uniformly bounded-weight generating set. The family is
standard concatenation with `n=d^2`, twice the conjectured IBM b=1
asymptotic overhead. It is not a fixed-weight BB solution. Minor
terminology: "dimension remains 12" in the parent's code-space
sentence means 12 encoded qubits; the Hilbert-space dimension is
`2^12`. This review used no new diagnostic computation.

## 7. Validation and handoff

Only the companion script `scripts/a44_s1_codex_checks.py` was added.
It uses the standard library, sparse finite supports, tiny polynomial
series over F4, and binary row reduction only for the two specific
Frobenius witnesses at `r=2,4`. It enforces a 55-second alarm and CPU
limit, checks peak RSS below one GiB, and writes no output files.
It deliberately does not import `a44_s0_patterns.py`, whose top-level
logger would overwrite the banked log.

The final diagnostic run passed in 0.006 seconds of measured script
time with peak RSS 18.4 MiB. Including repeated development runs and
one finite overlap inspection, the original witness investigation
remained below one second of measured execution. The later matroid
audit in §8 adds under three seconds; the total remains far below
the 180-second aggregate cap. An
initial uv cache-access error was resolved using its no-cache mode;
no dependency installation or old process interaction occurred.

Reproduce from the repository root:

```
PYTHONDONTWRITEBYTECODE=1 uv run --offline --no-sync --no-cache \
  --project experiments/bb_lab python \
  experiments/bb_lab/scripts/a44_s1_codex_checks.py
```

Original IBM work still outstanding: construct the b=0 discount, with
(5.8) the first concrete compression test in the supplied recurrence;
prove C(a) if that class is viable; and independently prove a matching
lower bound. Neither the obstruction nor the related-family theorem
may be reported as having completed those tasks.

## 8. Follow-up audit: a local F7 obstruction to graphic syndrome matrices

Audited the parent's `scripts/ibm_graphic_minor_check.py`, then replayed
its frozen-in-memory output using a separately written **row-rank**
elimination and explicit enumeration of the eight possible F2^3 span
vectors. No further witness or distance search was performed.

**Verdict: the certificate and stated geometric implication are valid.**
The parent run found a rank-112 patch with 128 qubit columns. Its
contracted set C contains 35 independent columns, and the seven retained
qubits, in binary-vector order 1 through 7, are

```
(0,0,3), (0,3,4), (1,6,5), (0,3,2),
(1,3,2), (0,0,5), (1,4,2).
```

The independent replay checked that C and these seven qubits are
disjoint, with no repetitions, and verified every one of the 128
rank identities both on the plane and after reduction modulo `(12,12)`.
The parent script took 2.2769 seconds and 17,760,256 bytes peak RSS;
the independent replay took 0.0081 seconds. The literal contracted list
is the parent's certificate output; the parent is preparing its
persistent certificate and dedicated checker. This audit did not modify
that script or create another certificate file.

### 8.1 What the script certifies

Its block-0 columns use beta and its block-1 columns use alpha, with
exactly the Laurent shifts of §1. Thus it represents the original X
syndrome map, not a substituted check convention.

`coordinate_basis` records a basis of original columns and reconstructs
every column from its returned coordinates. The retained three basis
positions first supply all seven distinct nonzero binary projections.
For the seven selected columns, the script keeps the extra coordinate
rows that already belong to the three-dimensional Fano row span and
contracts the basis columns corresponding to the other rows. Contracting
an independent subset of a represented matroid means quotienting the
ambient column space by its span. In these basis coordinates that removes
exactly the associated coordinates.

Most importantly, the final exhaustive subset test is independent of
that discovery heuristic: for every subset J of the retained seven,

`rank(C union J)-rank(C)=rank({binary vectors j+1: j in J})`.

This is precisely the rank-function definition of an F7 restriction
after contracting C and deleting every other column. It also precludes
a retained column lying in C's span. The separately implemented row-rank
replay confirms the identity without using `coordinate_basis` or the
parent's rank function.

**One non-mathematical harness issue:** the current discovery script
prints `no F7 found ...` and exits with status zero if the discovery
loop fails. Consequently an exit-code-only caller must not interpret
success as verification. The output status must be checked, or the
forthcoming fixed-certificate checker should fail explicitly on any
invalid or missing certificate. This does not affect the successful
certificate audited here.

### 8.2 Why the same minor exists on every sufficiently large torus

For the whole two-block qubit patch `0<=x,y<=7`, its syndrome rows lie
in `-3<=x<=8`, `-1<=y<=10`. Each coordinate interval contains exactly
12 integers and has diameter 11. Reduction modulo any `ell,m>=12`
is therefore injective on the entire syndrome-row set and on the
qubit patch. No two of the three terms of a syndrome column collide,
and no distinct syndrome rows merge. Additional torus syndrome rows
are zero on these selected columns and do not affect their ranks.

The plane patch is thus a column-matroid restriction of the actual
torus matrix, with exactly the same rank function. Its F7 minor persists
for every `ell,m>=12`, including all IBM members with `r>=2` and both
b values, without a spectral-dimension assumption. The argument does
not cover the `r=1` frames because at least one period is six.

### 8.3 Exact obstruction and limits

F7 is simple, rank three, and has seven elements. If it were graphic,
its representing graph would have no loops or parallel edges. Write
`r_i=v_i-1` for the rank of each nontrivial connected component. Then
`sum r_i=3`, and the number of edges of a simple graph is at most

`sum r_i(r_i+1)/2 <= ((sum r_i)^2+sum r_i)/2 = 6`.

This contradicts seven elements and handles disconnected graphs as
well as connected ones. Graphic matroids are closed under deletion
and contraction because these are graph edge deletion and contraction
(with any resulting loops or parallel edges allowed). A matrix with
the certified F7 minor therefore cannot have a graphic column matroid.

Invertible row changes, addition or deletion of dependent syndrome
rows, and qubit-column permutations preserve the column-matroid rank
function, up to relabelling. Every binary graph incidence matrix has
a graphic column matroid, including a reduced incidence matrix with
one vertex row omitted per component. Consequently these operations
cannot turn the original X syndrome into graph incidence on any of
the stated tori.

This rules out the specified **syndrome-row/physical-column** graphic
reformulation. It does not rule out CNOT column mixing, general local
Clifford equivalence, additional-qubit constructions, or a different
weighted description of the logical problem. It supplies no distance
bound and does not settle the IBM distance conjecture.
