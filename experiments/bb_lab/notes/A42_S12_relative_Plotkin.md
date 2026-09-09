# A42 S12 — exact physical Plotkin equations and an all-height matched-parity obstruction

2026-09-08. Bounded native Codex continuation of S10. Only the
original pure omega cylinder is considered. All Laurent x
supports are finite. No minimum-distance search is performed.

## 1. Outcome

The pure floor `d_q=3q` for all `q=2^a` is **not proved or
disproved**. This pass establishes three exact statements:

1. The physical even/odd split has an explicit two-equation
   syndrome system, with every boundary correction retained.
   Eliminating one unit entry gives a smaller Koszul pair and
   an exact four-term physical cost, not an ordinary Hamming
   distance on that smaller pair.
2. Either parity projection determines an actual cycle uniquely.
   Its relative quotient also determines the **full homology
   class**, through a scalar annihilator modulo a degree-five
   polynomial. The projections are not lower-height cycles.
3. Even after matching the full homology class, separately
   minimizing the two parity halves loses the target constant.
   At q=2, the global matched-parity relaxation is exactly four,
   whereas the banked actual distance is six. Explicit same-class
   witnesses transfer to every height q>=2, giving a nonzero
   class whose separately attained parity minima sum to at most
   `2q`, below the desired `3q`.

The inexpensive splice is not a cycle. Thus these witnesses
obstruct a decoupled Plotkin distance recursion, not a genuinely
joint-syndrome proof. They do not use S10's false weight-monotone
filtration normalization. The unresolved sharp inequality in
§7 is explicitly left open, rather than advertised as progress
merely because it has new coordinates.

## 2. Physical even/odd equations

Let `q=2n`, where n is a positive power of two, and put `s=t^2`.
Use

    R_n=F4[s]/(s^n+1)[x^+/-1],
    omega^2+omega+1=0,
    P=x^3+omega t+omega^2 t^2,
    Q=1+(x+x^2)t^3.

Every physical chain decomposes uniquely as

    f=a(s)+t b(s),       g=c(s)+t d(s).

The map is a partition of coordinates, so

    |(f,g)|=|a|+|b|+|c|+|d|.                          (2.1)

Define the following fixed polynomials in R_n:

    A=x^3+omega^2 s,       B=(x+x^2)s,
    Gamma=A+omega s B,
    C=1+s B^2,
    Lambda=omega+A B.                                 (2.2)

Thus `P=A+omega t` and `Q=1+t B`. The cycle equation is
exactly the pair

    A a+omega s b+c+s B d=0,
    omega a+A b+B c+d=0.                              (2.3)

For a generator `h=h_0+t h_1`, its boundary changes the four
physical parts by

    a -> a+h_0+s B h_1,
    b -> b+B h_0+h_1,
    c -> c+A h_0+omega s h_1,
    d -> d+omega h_0+A h_1.                           (2.4)

Equations (2.1)-(2.4) are identities before any minimization.
Unlike the contiguous-sheet cut in A42 §2.17.4, the split here
is into alternating physical rows. Its coupling is explicitly
present throughout the lower-period system; it is not silently
discarded as a harmless seam term.

## 3. An exact reduced complex with its full physical cost

Set

    b_tilde=b+B a,       d_tilde=d+omega a.             (3.1)

**Proposition 3.1.** All cycles are parametrized, without a
support restriction, by an arbitrary a and an effective cycle

    Gamma b_tilde+C d_tilde=0.                        (3.2)

The original chain is recovered by

    f_even=a,
    f_odd=b_tilde+B a,
    g_even=A a+omega s b_tilde+s B d_tilde,
    g_odd=d_tilde+omega a.                            (3.3)

Under a boundary change, the parameters transform as

    a -> a+h_0+s B h_1,
    b_tilde -> b_tilde+C h_1,
    d_tilde -> d_tilde+Gamma h_1.                      (3.4)

Consequently original homology is isomorphic to the homology
of the effective pair `(Gamma,C)`, with boundaries
`(C h_1,Gamma h_1)`.

Proof. Subtract the boundary with `h_0=a,h_1=0`, only as an
algebraic operation. The new first even component is zero,
and its first odd and second odd components are (3.1).
The first equation of (2.3) gives the new second even part;
substitution into the second gives (3.2). This proves (3.3)
and its converse. Expanding (2.4) gives (3.4). For every fixed
h_1, the free h_0 realizes any desired change in a, so the
stated quotient isomorphism follows. QED.

The exact cost to be minimized is nevertheless

    W_n(a;b_tilde,d_tilde)
      =|a|+|b_tilde+B a|
         +|A a+omega s b_tilde+s B d_tilde|
         +|d_tilde+omega a|.                          (3.5)

In particular

    d_(2n)=min W_n(a;b_tilde,d_tilde),                 (3.6)

over all a and all nonboundary effective cycles. This is an
exact reformulation retaining the boundary cost. It is not a
closed recurrence in d_n, nor permission to set a=0 at no cost.
The effective pair is different from the original pair. One
useful finite identity exhibiting its Frobenius structure is

    Gamma^2=P_(omega^2)(x^2,s)+omega^2 s C,
    P_(omega^2)(x^2,s)=x^6+omega^2 s+omega s^2.

This identity does not establish a Hamming isometry between
the two presentations.

## 4. A full-class invariant of either physical parity

The polynomial Lambda has degree five in x, leading coefficient
s, and constant coefficient omega. Both are units. It is
injective on compact Laurent chains. Moreover

    T_n=R_n/(Lambda)

is a free rank-five module over `F4[s]/(s^n+1)`; x is a unit
there. Let `ann_(T_n)(C)` denote the annihilator of C in T_n.

**Proposition 4.1.** The assignment

    [(f,g)] -> rho=c+A a mod Lambda                    (4.1)

is an isomorphism from the full original H_(2n) to
`ann_(T_n)(C)`. The odd parity of the same chain gives

    sigma=d+A b = B rho mod Lambda.                   (4.2)

Multiplication by B on this annihilator is invertible, with
inverse sB. Thus either relative parity label identifies the
full homology class, not just a quotient of it.

Proof. Eliminating d from (2.3) gives

    s Lambda b=Gamma a+C c.

Since `Gamma+C A=s B Lambda`, this is equivalent to

    s Lambda(b+B a)=C(c+A a).                         (4.3)

It follows that rho is annihilated by C modulo Lambda.
Conversely, for any even pair (a,c) satisfying this divisibility,
(4.3) determines b uniquely, since s Lambda is injective;
the second equation of (2.3) then determines d uniquely.
Both reconstructed parts are compact Laurent polynomials.

Under (2.4), rho changes by `s Lambda h_1`. If rho is zero
modulo Lambda, choose h_1 with `rho=s Lambda h_1` and then
`h_0=a+s B h_1`. The boundary has exactly the prescribed
even pair, and uniqueness forces its odd pair to agree as
well. This proves the kernel assertion. Any annihilator
element is attained by choosing a=0 and c a representative
of that element, so the map is surjective.

Finally the second equation of (2.3) yields

    d+A b=omega a+B c=B(c+A a)+Lambda a.

Equation (4.2) follows. On the annihilator, `C=1+sB^2`
acts as zero, so `sB^2=1`. QED.

The reconstruction assertion also proves that either actual
parity projection is injective on the space of cycles. This
does not make two independently chosen projections compatible.
Their relative labels can agree exactly while their original
joint syndrome fails.

## 5. A fixed same-class obstruction at q=2

All polynomials in this section, except for the displayed
factor t, are in `F4[x]`; here s=1. Define

    B=x+x^2,       A=x^3+omega^2,
    chi=1+x+x^2,
    r=x^2+omega^2 x+omega^2,
    J=x^2+x+omega^2.

Then

    C=chi^2,
    Gamma=(x+omega^2)chi,
    Lambda=Gamma r.

Put

    R=x^5+omega^2=(x+omega)rJ,
    d_0=(x+omega)J=x^3+omega^2 x^2+x+1,
    b_0=(x+omega)^2 J=x^4+x^3+omega^2 x+omega.          (5.1)

The actual compact cycle

    z=(t b_0, R+t d_0)                                (5.2)

has weight ten. Direct finite identities verifying (2.3) are

    R=r d_0,
    omega b_0+B d_0=R,
    A b_0+B R+d_0=0.

Its full homology label is R modulo Lambda. This is nonzero:
at `x=omega^2`, R evaluates to one and Lambda to zero.
Also `C R=Lambda b_0`, as required by Proposition 4.1.

Now define

    kappa=omega x R,
    R_prime=kappa R=omega x^11+omega^2 x,
    z_prime=t kappa z.                                (5.3)

Multiplication commutes with the differential, so z_prime is
a cycle. More importantly, it is in the **same full class**
as z, not merely a matched coarser label. Here is an explicit
boundary certificate. Put

    E=omega x^2 J (x+omega^2)^2,
    h=E+t(B E+b_0).

The polynomial identities

    R_prime+B R=Lambda E,
    C R=Lambda b_0

give

    z+z_prime=delta(h).                               (5.4)

For example, the odd first part of delta(h) is b_0 and
the odd second part is `Lambda E+A b_0=d_0+R_prime`.
For the even first part, use
`(kappa+B)b_0=C E`; it is kappa b_0. The other even
component then follows from the cycle equation, or by expansion.
Thus (5.4) is a finite polynomial identity at the actual pair.

The two parity weights of z are two and eight. Those of
z_prime are sixteen and two. In particular the separately
chosen inexpensive projections are

    even(z)=(0,R),       odd(z_prime)=(0,R_prime).      (5.5)

They each have weight two and come from the same nonzero
full class alpha. Their splice is

    z_splice=(0,R+t R_prime),       |z_splice|=4.

It is not a cycle. Multiplication by Q is injective on compact
chains: its highest x coefficient is the unit `t^3`, so the
highest output column of a nonzero input cannot cancel.
Since the second block of z_splice is nonzero, its syndrome
is nonzero. Direct multiplication gives syndrome weight ten.

### Exactness of the relaxed value four

For a nonzero full class beta at height 2n define

    e_(2n)(beta)=min{|even(z)|: [z]=beta},
    o_(2n)(beta)=min{|odd(z)|: [z]=beta}.                (5.6)

These are minima over projections of actual cycles, not over
arbitrary strings with an unverified syndrome. A nonzero cycle
cannot have a zero parity projection, by Proposition 4.1's
reconstruction argument. Nor can it have a weight-one even
projection: its rho would be a nonzero scalar monomial times
either one or A. Both are units in T_n, since

    Lambda+A B=omega.

But C is nonzero in T_n: it has x degree four, less than
the monic-unit-leading degree five of Lambda. A unit cannot
annihilate it. The same conclusion holds for odd projections,
by shifting the physical chain by t. Thus each quantity in
(5.6) is at least two at every height.

At q=2, (5.5) attains both lower bounds for the same alpha:

    e_2(alpha)=o_2(alpha)=2,
    min_(beta!=0)[e_2(beta)+o_2(beta)]=4.               (5.7)

The banked actual d_2 is six. Equation (5.7) is therefore an
exact obstruction to deriving that distance by independently
optimizing the two halves, even with exact full-class matching.

## 6. The obstruction persists at every tower height

For any `q=2n=2^a>=2`, put

    E_q=1+t^2+...+t^(q-2)=(1+t)^(q-2),
    z_q=E_q z,       z_prime_q=E_q z_prime.

Transfer from period two is injective on chains and, by S10's
local `u,u^2` calculation, an isomorphism on H_q for all
q>=2. Thus both displayed cycles represent the same nonzero
class alpha_q. Alternatively their equality of classes is
already witnessed by `delta(E_q h)` from (5.4); transfer
injectivity on homology supplies nontriviality.

The transfer makes n disjoint copies of each two-row chain
and preserves row parity. Consequently

    |even(z_q)|=2n=q,
    |odd(z_prime_q)|=2n=q,
    e_q(alpha_q)+o_q(alpha_q) <= 2q.                   (6.1)

This is a proved infinite family of upper bounds on the
matched-parity relaxation, strictly below the proposed pure
floor `3q`. No exact value for these higher-height parity
minima is inferred. The actual source-cycle weights are
`5q` and `9q`; neither is a distance counterexample.
Their cheap splice remains a noncycle, by injectivity of Q
or by transfer of its nonzero period-two syndrome.

Thus the defect in a scalar Plotkin recursion is not merely
forgetting the homology class or choosing a nonminimal
filtration representative. Even exact full-class matching
between separately feasible parity projections does not
retain their joint syndrome compatibility.

## 7. Precise residual and bounded verification

The physical formulas (2.3)-(2.4) and cost (3.5) do retain
that compatibility. The unresolved sharp assertion in these
coordinates is

    W_n(a;b_tilde,d_tilde) >= 6n

for every compact a and every nonboundary effective cycle
`Gamma b_tilde+C d_tilde=0`. This is equivalent to the original
floor, not a newly proved lower bound. No closed recurrence
for the numbers d_q has been obtained. In particular, replacing
the joint cost by separate parity minima cannot prove the
target, by (6.1). A successful argument must exploit the
compatibility and the shared correction variables rather than
charge them independently or discard them.

The new diagnostic is `scripts/a42_s12_codex_plotkin.py`. It
checks the fixed q=2 polynomial certificates, the exact
same-class boundary identity, the two inexpensive projections,
and their nonzero splice syndrome. It also checks the general
recursion formulas, the full physical cost, and boundary
updates on a fixed q=4 transfer, one specified correction,
and one specified noncycle. These are identity checks, not
a period enumeration or a minimum-distance verifier.

The first run passed in 0.000477 seconds with measured peak
RSS 18,710,528 bytes. The script has a 55-second wall alarm,
50-second CPU limit, and 768-MiB measured-RSS guard. The final
replay, including the explicit evaluation certificate for the
nonzero homology label, passed in 0.000481 seconds with peak
RSS 18,972,672 bytes. No SAT, coset search, support search,
growing-width matrix, old-process control, extra agent, git
mutation, Lean work, memory edit, or shared-file edit was
performed. Only this new note and its new diagnostic were
written. The bounded pass stops here without claiming the
pure-cylinder or IBM distance objective complete.
