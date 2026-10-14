# A42 S15 — genuine Cartier chain maps, exact class action, and their metric limitation

2026-09-08. Bounded native Codex continuation of S10/S12. Original
pure omega cylinder only; compact Laurent x support throughout.
No height, support, coset, or distance search was performed.

## 1. Outcome and hypotheses

The proposed four maps are actual semilinear chain maps. Their
coefficient square roots return the **same omega pair**, not its
coefficient conjugate. At every dyadic descent, their homology
images lie in the two-dimensional norm-transfer subspace. Every
nonzero source class has three or four nonzero images, with an
explicit classification below. There is no common kernel.

The resulting class-type recurrences give a genuine but nonsharp
growing lower bound. They do **not** prove `d_(2n)>=2d_n` or the
pure floor `d_q=3q`. Two precise limitations are proved:

1. The total physical cost factor three is sharp even on actual
   cycles, by a compact boundary valid at every dyadic q>=16.
2. Even replacing each output by its individual coset minimum
   does not permit a universal total-cost factor two: a genuine
   nontrivial height-two cycle of weight six has output class
   minima `3,3,3,5`, totaling fourteen.

Thus the obstruction is not just the cost of a poorly chosen
output representative. Different coupled inequalities remain
possible; no impossibility theorem for all Cartier methods is
asserted. S10's non-weight-monotone normalization and S12's
independent parity minima are not used as metric steps.

Fix F4 with `omega^2+omega+1=0`. For q a positive power of two,
write

    R_q=F4[t]/(t^q+1)[x^+/-1],        u=1+t,
    P=x^3+omega t+omega^2 t^2,
    Q=1+(x+x^2)t^3,
    D_q(f,g)=Pf+Qg,        delta_q(h)=(Qh,Ph),
    H_q=ker(D_q)/im(delta_q).

All chains have finite x support. Weight is physical F4 Hamming
support in the x,t basis, summed over both blocks. Let d_q be the
minimum weight of a nonzero class representative, and let
`dist_q(alpha)` be the minimum in a specified class, with
`dist_q(0)=0`. Minima exist because weights are nonnegative
integers and each class has a compact representative.

The exact action and the first growing bound below are proved
without assuming a distance certificate. A separately marked
sharpening uses the previously banked base fact `d_2=6`.
Nothing here closes binary HM or the original IBM torus metric.
The terms one-point and two-point below concern the two ordinary
omega-plane points, not S11 mixed homology classes or A42 mixed
binary fibre cycles.

## 2. Cartier extraction and all three chain degrees

Let q=2n. For i,j in {0,1}, define an additive map
`C_ij:R_(2n)->R_n` on monomials by

    C_ij(lambda x^(2k+i) t^(2ell+j))=lambda^2 x^k t^ell,

and send the other parity classes to zero. Square root in F4 is
squaring. Negative x exponents cause no difficulty. Replacing a
t exponent by that exponent plus 2n preserves parity and changes
the extracted exponent by n. Thus this definition is well-defined
on the periodic ring and preserves compact x support. No x-period
quotient is imposed; this is not an assertion about arbitrary
finite x-periods.

The map is F2-linear and F4-semilinear:

    C_ij(lambda z)=lambda^2 C_ij(z).

Frobenius gives a well-defined embedding `R_n -> R_(2n)`, squaring
both coefficients and monomials. For any a in R_n and any z in
R_(2n),

    C_ij(a^2 z)=a C_ij(z).                              (2.1)

To verify this, first take a monomial a, which shifts both source
parities by even amounts; extend using additivity in characteristic
two. The square of any lift of a is independent of that lift.
In particular (2.1) uses P and Q with their original coefficients
on the target. Omitting coefficient square roots instead would
conjugate the target coefficients.

**Theorem 2.1.** The following three maps form a semilinear chain
map between the height-2n and height-n Koszul complexes:

    degree 2: h       |-> C_ij(h),
    degree 1: (f,g)   |-> Phi_ij(f,g)=(C_ij(Qf),C_ij(Pg)),
    degree 0: s       |-> C_ij(PQ s).

Proof. Applying (2.1) gives

    D_n Phi_ij(f,g)
      =C_ij(P^2 Qf + Q^2 Pg)
      =C_ij(PQ(Pf+Qg)),

    Phi_ij(delta_(2n)(h))
      =(C_ij(Q^2h),C_ij(P^2h))
      =(Q C_ij(h),P C_ij(h))
      =delta_n(C_ij(h)).

These identities hold for arbitrary chains, not just a chosen
normal form of cycles. They prove the cycle and boundary claims.
QED.

The four extractions partition coordinates, so for every chain z,

    sum_(i,j) |Phi_ij(z)| = |Qf|+|Pg| <= 3|z|.          (2.2)

This is an exact cost formula followed by a convolution bound.
All boundary corrections remain admissible on both sides.

## 3. Local homology and a fixed height-two basis

Here is the algebraic anchor, included to make the all-height
quantifier independent of finite diagnostic samples. Let
`O=F4[t]_(t+1)` and `S=O[x^+/-1]`. Modulo Q,

    x^2=x+t^-3,
    P=A+B x,
    A=t^-3+omega t+omega^2 t^2,        B=1+t^-3.

The quotient by Q is free on 1,x over O; x is already invertible
because `x(x+1)=t^-3`. Multiplication by P has matrix

    [ A       B t^-3 ]
    [ B       A+B    ].

Modulo u^4,

    t^-3=1+u,
    A=omega^2 u+omega^2 u^2,
    B=u,
    det= A^2+AB+B^2t^-3 = omega u^3.

Its entry gcd has valuation one and its determinant valuation
three. The Smith invariants are u and u^2, up to units. The
nonzero determinant also makes P injective modulo Q, so the
Koszul complex resolves `M=S/(P,Q)`. Consequently

    M = O/(u) direct-sum O/(u^2) as O-modules,
    H_q = Tor_1^S(M,S/u^q) = M[u^q].                   (3.1)

In particular H_1 has F4 dimension two, and H_q has dimension
three for every q>=2. Localization changes no cylinder chains,
since `O/u^q=F4[u]/u^q`.

The period-resolution lift of multiplication by `u^(q-s)` from
height s to height q is that multiplication in degree zero and
the identity in degree one. Hence its induced map is the
inclusion `M[u^s] -> M[u^q]`. It follows that

    tau_q=u^(q-2): H_2 -> H_q is an isomorphism (q>=2),
    N_q=u^(q-1)=1+t+...+t^(q-1): H_1 -> H_q is injective. (3.2)

Call the latter image the socle subspace; it is `M[u]`. At q=1
it is all of H_1. These are class identifications, not permission
to replace a physical chain without charging its weight.

Put `chi=1+x+x^2`. At height one,

    P_1=(1+x)chi,       Q_1=chi,
    Z_1={(h,(1+x)h): h in F4[x^+/-1]},
    H_1=F4[x^+/-1]/chi.                                (3.3)

At height two use the following three cycle classes as basis:

    e_0=((1+t),(1+x)(1+t)),
    e_1=x e_0,
    e_2=(x+omega, 1+t r),
    r=x^2+omega^2 x+omega^2.                            (3.4)

The first two are norm transfers. For the third, write
`P=p+omega u`, `Q=chi+(x+x^2)u`, where `p=(1+x)chi` and u^2=0.
Its constant syndrome vanishes, and its linear syndrome vanishes
by the finite identity

    chi r=(omega+x+x^3)(x+omega).

For independence, the constant-u first block of any relation
between e_2 and the other two would require chi to divide a
nonzero multiple of `x+omega`, which is impossible. A relation
between e_0,e_1 would then require chi to divide a nonzero
polynomial of degree at most one at the u coefficient. Negative
x powers do not change this divisibility conclusion. Dimension
three in (3.1) completes the basis proof.

At every q>=2, coordinates (a,b,c) mean the class
`tau_q(a e_0+b e_1+c e_2)`. The socle subspace is exactly c=0.

## 4. Exact action at every dyadic height

**Theorem 4.1.** Under (3.2), every Cartier homology map factors

    H_(2n) --tau_(2n)^-1--> H_2 --A_ij--> H_1 --N_n--> H_n.
                                                               (4.1)

In the height-one label h modulo chi, its middle map is

    A_00=(1+x)a^2 + x b^2 + omega^2 c^2,
    A_10=      a^2 +(1+x)b^2 +       c^2,
    A_01=(1+x)a^2 + x b^2 + omega x c^2,
    A_11=      a^2 +(1+x)b^2 +(x+omega^2)c^2.           (4.2)

Proof. The all-height identity is

    tau_(2n)=u^(2n-2)=N_n^2.

Thus Cartier extraction of a transferred cycle factors out N_n
by (2.1). Multiplication by N_n on the target evaluates every
remaining t power at 1. On polynomial lifts, that evaluation
after extraction is exactly the height-two-to-one extraction.
Changing a height-two lift by `(t^2+1)w` changes its extraction
by a multiple of 1+t, which N_n kills. Thus this argument respects
all periodic reductions, including n=1.

It remains to compute only the three basis images. Their first
blocks after multiplying by Q are, at height two,

    Q(e_0)_f=chi(1+t),
    Q(e_1)_f=x chi(1+t),
    Q(e_2)_f=(x+omega)+t(x^3+omega^2 x^2+omega x).

Extracting each parity and taking coefficient square roots gives
the three columns of (4.2). The second block is forced to be
`(1+x)` times the first by (3.3). Semilinearity accounts for the
squares of a,b,c. The chain-map theorem extends this computation
to every representative of every class. QED.

In channel order 00,10,01,11, the four kernels are the lines

    F4*(omega,omega,1),
    F4*(1,0,1),
    F4*(0,omega^2,1),
    F4*(omega^2,1,1),                                  (4.3)

respectively. Each map has rank two: already its a,b columns are
independent. Solving its two coordinate equations gives (4.3).
The four lines are distinct and meet only at zero. Since norm
transfer is injective, a nonzero class has exactly three nonzero
outputs on their union, and four everywhere else. There are twelve
nonzero vectors on those lines and fifty-one other nonzero vectors
in the three-dimensional source. These counts are a consequence
of the displayed linear algebra, not a census.

There is also an exact output relation in H_1:

    A_00+A_01=omega(A_10+A_11)=omega(x+omega)c^2.         (4.4)

Thus the outputs are not four unrelated classes. In particular,
socle inputs have equal j=0 and j=1 output classes for each i.

## 5. Refinement by ordinary-point support

For a socle class, use the two labels `h(omega),h(omega^2)` of
its unique height-one precursor. This gives three nonzero types:
both labels nonzero, only the omega label nonzero, and only the
omega^2 label nonzero.

The full classification, including inputs outside the socle,
is as follows. Set

    A=a^2, B=b^2, C=c^2,
    L=A+omega^2 B+C=(a+omega b+c)^2,
    R=A+omega B.

Evaluation of (4.2), in channel order 00,10,01,11, gives

    at omega:   (omega^2 L, L, omega^2 L, L),
    at omega^2: (omega R+omega^2 C, R+C, omega R+C, R).  (5.1)

If C is nonzero, dividing the even-i entries in the second row
by omega gives the four values `R+omega C,R+C,R+omega^2 C,R`.
Exactly one is zero, since the four multiples of C exhaust F4.
If C=0, each nonzero point component simply remains nonzero in
all four channels. This proves the following table.

| Source type | Number of source classes | Nonzero output types |
| --- | ---: | --- |
| c=0, both point labels nonzero | 9 | four both-point |
| c=0, only omega label nonzero | 3 | four omega-only |
| c=0, only omega^2 label nonzero | 3 | four omega^2-only |
| c!=0, L!=0 | 36 | three both-point, one omega-only |
| c!=0, L=0 | 12 | three omega^2-only, one zero |

Every output is in the target socle, at all heights. In
particular, all further descents of a nonzero output have four
nonzero outputs, preserving its point-support type. The simple
and double u factors are consistent with this orientation: the
coefficient of u in P modulo Q is `omega^2+x`, so the double
factor is at x=omega^2. No y-mirror identification is used.

## 6. Base-class metric and nonsharp analytic recurrences

**Lemma 6.1.** At height one, the nine classes with both point
labels nonzero have exact distance three. The six classes with
exactly one nonzero point label have exact distance five.

Proof. Write every cycle as `(h,(1+x)h)` using (3.3). A nonzero
product `(1+x)h` has at least two terms: evaluation at x=1 is
zero, whereas a nonzero Laurent monomial evaluates nonzero.
Thus every nonzero cycle has weight at least three.

The nine monomial labels `lambda x^k mod chi`, with lambda nonzero
and k=0,1,2, give all nine pairs of nonzero point labels. Indeed
their evaluation ratios exhaust F4's three nonzero elements.
Their cycles have weight three.

For a one-point class, h cannot be a monomial. If h has at least
three terms, the product bound just proved gives total weight at
least five. If h has two terms, multiply by a unit to write it as
`alpha+beta x^m`, with m>=1 and alpha,beta nonzero. For m>=2,
its product with 1+x has four distinct terms. For m=1 it has
three terms unless alpha=beta. In that exception the total
weight is four, but the label is a scalar monomial times 1+x;
it is nonzero at both ordinary points and is not a one-point
class. Finally scalar multiples of `x+omega` and `x+omega^2`
attain weight five and exhaust the six one-point classes. QED.

Define b_a, s_a, and t_a to be the minimum physical weights at
height 2^a among nonzero socle classes of, respectively, both-point,
omega-only, and omega^2-only type. All representatives, including
all possible boundary corrections, enter these minima. Then

    b_0=3,             s_0=t_0=5,
    b_a >= ceil(4 b_(a-1)/3),
    s_a >= ceil(4 s_(a-1)/3),
    t_a >= ceil(4 t_(a-1)/3).                          (6.1)

For a>=1, let v_a and w_a be the minima for source types
`c!=0,L!=0` and `c!=0,L=0`. The remaining inequalities are

    v_a >= ceil((3 b_(a-1)+s_(a-1))/3),
    w_a >= t_(a-1),
    d_(2^a)=min(b_a,s_a,t_a,v_a,w_a).                  (6.2)

Proof of the inequalities. Apply (2.2) to an arbitrary actual
source cycle. The table in §5 specifies the nonzero class type of
each output. Each output's actual weight is bounded below by the
minimum of its specified type. Summing and dividing by three
proves (6.1) and (6.2). No source or target normalization is
assumed to preserve weight. QED.

These are lower inequalities, not sharp distance recurrences.
Nevertheless they imply the entirely analytic growing bound

    d_(2^a) >= ceil(3*(4/3)^a)        for every a>=0.    (6.3)

To see this, (6.1) gives `b_a>=3*(4/3)^a` and
`s_a,t_a>=5*(4/3)^a`. In units of `(4/3)^(a-1)`, the right sides
for b_a,s_a,t_a,v_a,w_a are at least

    4, 20/3, 20/3, 14/3, 5,

respectively. Their minimum is four, proving (6.3) also for the
outside-socle classes. The base case is Lemma 6.1.

**Sharpening using the banked base certificate d_2=6.** The
canonical e_0 has weight six and both point labels nonzero, so
this additional input gives `b_1=6`. Hence

    b_a >= 6*(4/3)^(a-1)              for a>=1,
    d_(2^a) >= ceil(5*(4/3)^(a-1))    for a>=2.         (6.4)

For the second assertion, the five lower coefficients in units
of `(4/3)^(a-1)` are now at least

    6, 20/3, 20/3, 37/6, 5,

using the improved b_(a-1) in v_a. Thus their minimum is five.
The base certificate is not rerun or re-proved here; (6.4) is
explicitly the consequence that uses it. The exponent in these
bounds is `log_2(4/3)`, not one. Neither bound approaches the
sharp target `3*2^a` in relative error.

## 7. Two exact cost obstructions

### 7.1 The factor three is sharp even on cycles

The parent supplied the following finite-plane construction;
its distinctness argument and coefficients were independently
checked here. Take the compact generator h=PQ and its boundary

    z=delta(h)=(P Q^2,P^2 Q).

Each block has nine distinct monomials, so |z|=18. Q^3 has nine
terms, all with t exponents divisible by three. Multiplication
by P separates them into the three distinct t residues 0,1,2,
giving `|P Q^3|=27`. Likewise P^3 has nine terms with x exponents
divisible by three. The three x residues in Q separate them,
giving `|P^3 Q|=27`. All coefficients are nonzero.

For completeness, the nine-term claim follows by multiplying
each three-term polynomial by its Frobenius square. For Q,
write `H=x+x^2` and

    Q^3=1+H t^3+H^2 t^6+H^3 t^9,
    H^3=x^3+x^4+x^5+x^6.

For P, the groups of terms with x degrees 9,6,3,0 have sizes
1,2,2,4, respectively; their t exponents within each group are
distinct. The largest t exponent in P Q^3 is eleven and in
P^3 Q is nine. Therefore no t wrap occurs for any dyadic q>=16.
Equation (2.2) is consequently an equality:

    sum_(i,j)|Phi_ij(z)|=27+27=54=3|z|.                (7.1)

This rules out a smaller universal coefficient based only on
the actual cycle equation. Since z is a boundary, it alone
says nothing about a quotient-aware minimizer inequality.

The suggested factors, with the square-root convention, also
check exactly in the unwrapped plane:

    C_00(PQ)=omega t+omega^2 x t^2,
    C_10(PQ)=x+omega^2 t^2,
    C_01(PQ)=omega^2+x^2 t+omega x t^2,
    C_11(PQ)=x^2 t+omega t^2.                          (7.2)

They generate the four output boundaries in (7.1).

### 7.2 Individual coset minimization still cannot give factor two

Use the actual nontrivial cycle e_2 in (3.4), of weight six.
Its four height-one output labels, in the fixed channel order,
are

    omega^2,       1,       omega x,       x+omega^2.

The first three have both point labels nonzero, while the last
has only its omega label nonzero. By Lemma 6.1 their **exact
coset minima**, already attained by these outputs, are

    3, 3, 3, 5.

In particular all four outputs are nontrivial, but

    sum_(i,j) dist_1([Phi_ij(e_2)])=14 > 12=2|e_2|.    (7.3)

This disproves a universal factor-two bound on the sum of the
four individually minimized output costs, even restricted to
nontrivial source cycles with all four outputs nonzero. Using
the banked d_2=6, e_2 is itself globally distance-minimal, so
restricting the proposed bound to source minimizers does not
repair it. Even without that base certificate, its class has a
minimum at most six and its four output classes are unchanged
by every source boundary correction; the same obstruction holds
with the source class minimum on the right side of (7.3).

A related all-height example is the canonical cycle

    mu_(2n)=(N_(2n),(1+x)N_(2n)),       |mu_(2n)|=6n.

Its raw output weights are `4n,3n,4n,3n`, with sum 14n. This
alone already defeats a raw factor-two estimate. Unlike (7.3),
some of these canonical raw outputs admit shorter representatives:
the label 1+x equals x^2 modulo chi. Thus the two examples
separate representative-level and genuine coset-level losses.

## 8. Exact remaining metric gap and diagnostic scope

The identities (2.1), the full chain map, the universal action
(4.1), the kernels and point-type table, and the lower recurrences
(6.1)-(6.2) are proved. Their infinite quantifiers follow from the
norm-transfer identity and the local Tor calculation, not from
finite checks.

The sharp floor remains open. To improve these recurrences to
doubling, one needs a further weight-sensitive constraint on
simultaneously realizable Cartier outputs and their source,
beyond their class types, relation (4.4), and total-cost factor
three. Simply asserting a factor-two inequality for the four
output coset minima is not such a lemma: it is false by (7.3).
The present note supplies no replacement sharp constraint.

The new diagnostic `scripts/a42_s15_codex_cartier.py` reads only
the previously inspected finite-polynomial helpers in the S6
script, with bytecode writing disabled. It checks:

- one fixed height-four chain identity on arbitrary chains,
  including negative x exponents and periodic representatives;
- the three basis images, four fixed kernel lines, and the
  point-evaluation coefficient identities, without enumerating
  classes;
- one fixed transfer instance from height two to four;
- the fixed height-sixteen boundary in (7.1) and the four factors
  (7.2), with no height loop;
- the actual height-two obstruction (7.3) and the canonical
  height-four raw costs.

One execution passed in 0.001209 seconds, with measured peak RSS
18,825,216 bytes. It enforces a 55-second alarm and 50-second CPU
limit and checks the 768 MiB peak-RSS cap. The script does not
test or assert the all-height sharp floor. Only the two assigned
new S15 files were written; no shared note, git state, old process,
registry, Lean file, or memory file was changed.
