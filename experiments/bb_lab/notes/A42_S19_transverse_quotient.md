# A42 S19 — a nonzero cylinder-image class violates the six-column floor

2026-09-08. Bounded native Codex pure-omega pass. Only this note and
`scripts/a42_s19_codex_transverse.py` are new. No support, distance,
coset, height, SAT, or transfer-state search was performed.

## 1. Outcome and exact scope

The proposed algebraic injection is correct. Folding x modulo six
injects cylinder homology into one winding axis of torus homology.
At every dyadic q>=2 the dimensions over F4 are three and six,
respectively. At q=1 the corresponding dimensions are two and four,
and the injection still holds.

However, the proposed metric premise is **false even with the exact
cylinder-image restriction retained**:

> For every dyadic q>=4, the six-column torus has a nontrivial cycle
> in the cylinder-image axis of physical F4 weight exactly `3q/2`.
> Its two winding coordinates are `(omega^2,0)`.

An explicit witness is

    w_q=(u^(q-2) P,0),       u=1+t.                    (1.1)

Thus its torus class minimum is at most `3q/2`, strictly below 3q.
No assertion that this upper bound is the exact torus minimum is
made. This is a counterexample to the proposed transverse-quotient
premise, **not** a counterexample to the original pure cylinder
floor. The displayed Laurent representative before x folding has
a nonzero seam syndrome.

The height-two precursor has winding coordinates `(omega^2,1)`
and is not in the cylinder image. Norm transfer to every q>=4
kills that second coordinate while preserving the first. This
distinction is essential: the all-height witness cannot be dismissed
as merely an opposite-winding cheap string.

Throughout, fix F4 with `omega^2+omega+1=0` and use exactly

    P=x^3+omega t+omega^2 t^2,
    Q=1+(x+x^2)t^3,
    D(f,g)=Pf+Qg,       delta(h)=(Qh,Ph).

Cylinder chains belong to
`R_q=F4[t]/(t^q+1)[x^+/-1]` and have finite x support. Torus chains
belong to `T_q=R_q/(x^6+1)`. Weight counts nonzero F4 coefficients
in the physical x,t basis, summed across both blocks. No binary
HM weight is used here.

## 2. Local plane module and the transverse annihilator

Let `O=F4[t]_(t+1)`, `S=O[x^+/-1]`, `u=1+t`, and
`M=S/(P,Q)`. Localizing at t=1 is exact for every dyadic cylinder:
`O/u^q=F4[t]/(t^q+1)`. No assertion about other plane points away
from this localization is needed.

For completeness, the S10 local anchor has the following finite
proof. Modulo Q, the ring is free over O on 1,x, with

    x^2=x+t^-3,
    P=A+B x,
    A=t^-3+omega t+omega^2 t^2,       B=1+t^-3.

The equation `x(x+1)=t^-3` makes x invertible already. Multiplication
by P has matrix

    [ A       B t^-3 ]
    [ B       A+B    ].

Modulo u^4,

    t^-3=1+u,
    A=omega^2 u+omega^2 u^2,
    B=u,
    det=A^2+AB+B^2t^-3=omega u^3.

The entry gcd therefore has valuation one and the determinant
valuation three. Smith form over the DVR O gives

    M = O/(u) direct-sum O/(u^2) as O-modules.          (2.1)

In particular u^2 M=0 and `dim_F4 M=3`. The nonzero determinant
makes P injective modulo Q, so the Koszul complex of P,Q is a
free resolution of M. This also justifies the homology calculations
below without a separate local-multiplicity certificate.

Put `v=x^6+1`. There is an exact, unlocalized polynomial identity

    P^2=v+omega^2 u^2+omega u^4.                       (2.2)

Indeed `P^2=x^6+omega^2 t^2+omega t^4`, and substituting t=1+u
gives (2.2). Modulo P,Q, (2.1) and (2.2) imply

    v M=0.                                           (2.3)

Thus x^6+1 annihilates this localized plane module exactly, not
just its reduced points. This is not a claim that x^6+1 annihilates
the entire unlocalized plane module at unrelated t values.

## 3. Projection, both winding coordinates, and exact image membership

### 3.1 Class injection and dimensions

Let `C_q=K(P,Q) tensor_S S/u^q` be the cylinder complex.
Multiplication by v is injective on every chain module: v is
monic in x, and the highest Laurent column of a nonzero finite
chain survives even though its coefficients may be nilpotent.
Consequently

    0 -> C_q --v--> C_q -> C_q/v C_q -> 0              (3.1)

is a short exact sequence of complexes. The cylinder groups are

    H_1(C_q)=M[u^q],       H_0(C_q)=M/u^q M.

By (2.3), v acts as zero on both. The long exact sequence gives

    0 -> H_cyl(q) --pi_*--> H_tor(6,q)
      --kappa--> M/u^q M -> 0.                       (3.2)

The map pi is the actual physical folding of x coefficients,
so it cannot increase Hamming weight and maps boundaries to
boundaries. Equation (3.2) proves its injectivity on homology.
For q>=2 the dimensions in (3.2) are 3,6,3. For q=1 they are
2,4,2, using (2.1).

### 3.2 An explicit comparison map fixes the axes

Here is a concrete coordinate convention, so that image membership
does not depend on an unnamed abstract splitting. For a torus cycle
(f,g), choose lifts F,G in S and write

    P F+Q G=u^q A+v B.                                (3.3)

The periods `(u^q,v)` form a regular sequence: after reduction by
u^q, multiplication by the monic v is still injective. Modulo P,Q,
the pair `(A,B)` is a cycle in the period Koszul complex on M.
It gives the standard comparison isomorphism

    H_tor(6,q) = M[u^q] direct-sum M/u^q M.            (3.4)

The first coordinate is the t-period winding and the second is
the transverse x-period winding. In this convention

    pi_* H_cyl(q)=M[u^q] direct-sum {0}.               (3.5)

Proof and well-definedness. Changing the decomposition (3.3)
changes `(A,B)` by a period boundary `(v H,u^q H)`. Changing lifts
of f,g adds elements of `(P,Q)` to the corresponding coordinates.
An original Koszul boundary has zero syndrome and hence gives
the zero class. Conversely, any period cycle on M lifts to some
(A,B) whose combination `u^q A+v B` lies in `(P,Q)`, and thus
admits F,G as in (3.3).

For injectivity, a zero period class can first be adjusted by
`(v H,u^q H)` so that both A and B lie in `(P,Q)`. Write
`A=P a+Q b`, `B=P c+Q d`. Subtracting the corresponding period
terms from F,G makes their syndrome identically zero in S. Since
the P,Q Koszul complex is exact there, this pair is a boundary;
reducing modulo the periods proves the original torus class was
zero. This establishes the isomorphism directly.

Finally, because v M=0, period cycles are exactly pairs with
first coordinate in M[u^q]; period boundaries change only the
second coordinate by u^q M. This gives (3.4). A lifted cylinder
cycle already has `P F+Q G=u^q A`, so it maps to `(A,0)`, proving
(3.5) and fixing the winding direction. QED.

Equivalently, lift a torus cycle only in x to R_q. Its syndrome
is v times some b. Then

    kappa([f,g])=[b] in R_q/(P,Q)=M/u^q M.            (3.6)

Thus (3.6), or the second coordinate in (3.4), is an exact
cylinder-image test. It is not a test for merely nontrivial torus
homology.

## 4. Explicit nonzero image-class witness of weight 3q/2

For dyadic q>=2 put

    tau_q=u^(q-2)=1+t^2+t^4+...+t^(q-2),
    w_q=(tau_q P,0) in T_q^2.                        (4.1)

Multiplying (2.2) by tau_q gives an exact lift identity in S:

    P(tau_q P)
      =u^q(omega^2+omega u^2)+v u^(q-2).              (4.2)

Thus w_q is an actual torus cycle at every q>=2. Formula (3.3)
gives its exact class coordinates as

    [w_2]=(omega^2,1),
    [w_q]=(omega^2,0) for every dyadic q>=4,           (4.3)

because u^2 M=0 and q-2>=2 in the second line. The scalar
omega^2 is nonzero in M, so these classes are nontrivial. The
second line belongs to the cylinder image, while the first does
not. This is the required winding distinction.

The class omega^2 is a constant unit in M. In particular its
components at both ordinary points x=omega and x=omega^2 are
nonzero. Winding coordinates in (3.4) are not the same distinction
as ordinary-point support, and neither is S11 mixed homology or
binary mixed-fibre support.

In the physical torus ring, `t^2 tau_q=tau_q`, so

    tau_q P=(x^3+omega^2)tau_q+omega t tau_q.           (4.4)

The three displayed supports are disjoint: column three on even
rows, column zero on even rows, and column zero on odd rows. Each
has q/2 sites, with nonzero coefficients. Therefore

    |w_q|=3q/2.                                      (4.5)

For example, at q=4 the first block is exactly

    (x^3+omega^2)(1+t^2)+omega(t+t^3),                 (4.6)

and the second block is zero. It has six sites, whereas the
proposed image-class floor would require twelve.

### An independent nonboundary detector

The nontriviality in (4.3) also has a one-point certificate that
does not require a finite boundary-rank computation. For any torus
cycle and any polynomial lifts F,G, define

    ell_q(f,g)=[u^q] (P(omega,1+u)F(omega,1+u)
                         +Q(omega,1+u)G(omega,1+u)).  (4.7)

Here `[u^q]` means the exact coefficient in the local expansion,
not a numerical evaluation. The expression is divisible by u^q:
the other period v vanishes at x=omega. Changing a t-period lift
adds u^q times P or Q; both vanish at this ordinary point when
u=0, so the displayed coefficient does not change. Changing an
x-period lift has no effect. Every torus boundary has ell_q=0,
by taking lifts of the form `(QH,PH)`.

For (4.1), use `F=tau_q P`, `G=0`. Since

    P(omega,1+u)=omega u+omega^2 u^2,

the expression in (4.7) is

    omega^2 u^q+omega u^(q+2).

Hence `ell_q(w_q)=omega^2!=0`. This independently certifies an
actual nonboundary, including at q=4. Image membership remains
the separate assertion supplied by the second coordinate of (4.3).

## 5. Why the transverse winding disappears under norm transfer

The same phenomenon follows from a useful exact chain map. Norm
transfer from the height-two torus to the height-q torus multiplies
physical chains by tau_q. On the period resolutions ordered as
`(u^2,v)` and `(u^q,v)`, respectively, a lift is

    degree 0: multiplication by tau_q,
    degree 1: diag(1,tau_q),
    degree 2: identity.

Indeed `tau_q u^2=u^q`, and the degree-two boundary vectors
`(v,u^2)` and `(v,u^q)` are compatible with these maps. Therefore
the induced homology map in (3.4) is

    (a,b) |-> (a,u^(q-2)b).                           (5.1)

For q>=4 the entire transverse winding component is killed by
u^2 M=0, while the cylinder winding is unchanged. Applying this
to the actual height-two cycle `(P,0)` with coordinates
`(omega^2,1)` produces exactly (4.3).

This is a construction of new torus cycles, with the precisely
charged physical weight increase q/2 from disjoint row replication.
It is not S10's false weight-nonincreasing normalization of an
arbitrary cycle. Both components and the full joint syndrome are
retained, unlike the S12 independently minimized parity relaxation.

## 6. What the obstruction does and does not prove

Define

    d_image(q)=min{|z| : z is a torus cycle,
                         0!=[z] in pi_* H_cyl(q)}.

The proposed sufficient route was valid: if `d_image(q)>=3q`,
class injection and physical contraction would give the cylinder
floor. But (4.3)-(4.5) prove instead

    d_image(q)<=3q/2       for every dyadic q>=4.       (6.1)

Thus the proposed metric lemma is false, not merely unproved.
Remembering the exact image homology class does not repair it.

Before folding x, the displayed finite Laurent representative
has syndrome

    D_q(tau_q P,0)=u^(q-2)(x^6+1)!=0 in R_q.          (6.2)

The right side is nonzero because u^(q-2) is nonzero and x is
unbounded. Hence this representative is not a cylinder cycle.
Nevertheless its torus class really does come from the cylinder;
indeed (3.6) permits correction of a lift by multiples of x^6+1
to produce an exact compact cylinder cycle with the same folded
vector. Such corrections can add physical support that disappears
under folding. No weight control on those corrections is asserted.
In particular (6.1) does not refute the original pure cylinder
floor, and it does not identify any low-weight actual cylinder
counterexample.

The S8 six-column counterexample concerns binary canonical weighted
HM. It was read as a warning, not imported as evidence here. The
present obstruction is independently an actual F4 torus cycle with
explicit cylinder-image coordinates. The S10, S12, and S15 ledgers
remain consistent with it; no normalization, independent parity
minimum, or false Cartier total-cost improvement was used.

A repair would need additional information about the physical cost
lost at the x seam, not just the image class in (3.5). This pass
does not supply such a repair. The original all-height pure floor
remains open, and this bounded transverse-quotient route stops at
the explicit counterexample (4.1).

## 7. Fixed diagnostic and resource record

The new script uses the inspected, read-only S6 finite-polynomial
helpers, with bytecode writing disabled. It verifies the local
Smith valuations, the exact period identity (2.2), one lifted
height-four winding certificate (4.2), the explicit six-site
height-four torus cycle and its nonzero detector, and the fixed
height-two precursor with transverse coordinate one. A single
torus boundary is also checked to have zero detector as a sanity
check; the proof in §4 establishes this for all boundaries.

Only these fixed identities are evaluated. There is no height loop,
rank-based distance check, support search, or transfer-state model.
All infinite quantifiers follow from the polynomial and chain-map
proofs above, not from the two fixed diagnostic heights.

One execution passed in 0.000204 seconds with measured peak RSS
18,071,552 bytes. The script enforces a 55-second wall alarm and
50-second CPU limit and checks the 768 MiB peak-RSS cap. Only the
two assigned new S19 files were written; no shared note, sibling
script, git state, memory, Lean file, or old process was changed.
