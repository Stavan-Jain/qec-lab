# Exact next-weight gap for the dyadic canonical cylinder

2026-09-08. Consolidation of the completed finite-graph and Laurent
arguments. No new enumeration is used. This is an exact theorem about
one affine cylinder problem, not the full IBM-code distance.

## 1. The theorem

For q=2^a with a>=1, work in

    R_q=(F2[t]/(t^q+1))[x^+/-1],
    P=x^3+t+t^2, Q=1+(x+x^2)t^3, N_q=sum_(j=0)^(q-1) t^j.

All x supports are finite. For a two-block chain e=(f,g), let
D e=Pf+Qg, let |e| be its physical binary Hamming weight, and put

    e_star=(x^-3 N_q,0).

**Theorem 1.1.** Among all solutions of D e=N_q, the unique minimum
is e_star, of weight q. The next distinct weight is exactly q+4:

    min {|e| : D e=N_q, e!=e_star} = q+4.              (1)

In particular the entire weight-q+2 shell is empty for every dyadic
q>=2, not just the part with a large number of canonical holes.
The x coordinates in this assertion are not folded or bounded.

## 2. Complete proof of the empty shell

The original-lane uniqueness theorem in `A42_S7_spatial_equality.md`
gives |e|>=q, with equality only at e_star. Evaluating the syndrome
equation at x=t=1 gives |e|=q modulo two, because P and Q both
evaluate to one. Thus the excess delta=|e|-q is a nonnegative even
integer, and only delta=2 must be excluded to prove the lower bound
in (1).

The S7 graph identity is

    delta=L+sigma,

where L is the nonnegative even aggregation loss and sigma is
the sum of nonnegative integer edge slacks. At delta=2 there
are precisely the following alternatives.

| alternative | complete exclusion |
| --- | --- |
| L=2, sigma=0 | `A42_S11_next_shell.md` |
| L=0, one slack-two edge | `A42_S13_single_charge.md` |
| L=0, two slack-one edges | the four steps below |

For the last alternative, the fixed 16-state graph has 46 unordered
compatible charge pairs. The following disjoint exclusions exhaust
them, with actual Laurent positions retained throughout:

1. `A42_S14_double_charge.md` excludes 34 pairs by spatial pinches,
   impossible return paths, and explicit short cyclic identities.
2. `IBM_CANONICAL_TERMINAL_CHARGES.md` excludes the four pairs
   with entry 16 or 19 and exit 9 or 10.
3. Of the eight remaining pairs, `A42_S16_D_excursion.md` excludes
   every walk containing D->C. Its projected equations are
   X F=1 and F=1, with X=x^3 unbounded.
4. `A42_S17_no_D_excursion.md` excludes those same eight pairs
   when D->C is absent. A triangle-return pinch removes arbitrary
   repetitions; then the exit and CA equations force the same
   monomial E to satisfy X E=E. The two height-four words have
   a separate wrapped pinch.

Steps 3 and 4 are complementary cases, not sixteen additional
pairs: the count is 34+4+8=46. There is no residual edge type,
unclassified actual position, or assumption that a finite-height
test stands for the universal assertion. The q=2 empty shell is
also proved directly in S11 and is included in these exclusions.
This proves that every nonleader has weight at least q+4.

## 3. A matching representative at every height

Take h=x^-3 and its actual boundary

    c=(Qh,Ph)
     =(x^-3+(x^-2+x^-1)t^3, 1+x^-3 t+x^-3 t^2).

For every q>=2, each block has exactly three distinct sites.
At q=2, t^2=1 does not create a collision with the other row-zero
term, since that term has x exponent zero rather than -3. The
first block meets e_star at exactly one site, x^-3 t^0, and the
second block meets none. Therefore

    D(e_star+c)=N_q,
    |e_star+c|=q+|c|-2|c intersect e_star|=q+6-2=q+4.

The syndrome identity follows from D c=PQh+QPh=0. This supplies
an explicit attaining representative, completing equality (1).
The restriction q>=2 is essential: at q=1, P=x^3 and the same
boundary has only four sites, yielding the smaller excess two.

## 4. Consequences for the canonical hiding-mass problem

Let S be f-column zero together with g-columns zero and one, all
q rows, and let k=|supp(e) intersect S|. The canonical mixed-fibre
weight formula is

    W=6q+3 delta-2k.

If e=e_star, then delta=k=0 and W=6q. Otherwise (1) gives delta>=4.
Consequently, for every solution with k<=6,

    3 delta-2k>=0, and hence W>=6q.                   (2)

The minimum over this restricted collection is exactly 6q, since
the leader is admissible. To check its boundary-dressing status
directly, write 1_S=(N_q,(1+x)N_q). Then

    e_star+1_S = boundary(x^-3(1+x)N_q).

Thus (2) closes the holes-four-through-six strata in addition to
the previously closed smaller-hole strata, for this specified S.
It does not assert (2) for larger k or for an arbitrary omega
support. At q=2, S has only six sites, so this particular canonical
model is fully covered; no growing-family conclusion follows.

There is also a stronger connected-defect count. If c=e+e_star
is decomposed into nonzero components of its syndrome-interaction
graph, each component c_i is a cycle. Applying (1) to e_star+c_i
gives

    delta_i=|c_i|-2|c_i intersect e_star|>=4.

The excesses add because the component supports are disjoint.
Hence the number of nonzero components is at most delta/4,
improving the earlier delta/2 bound. This does not prove the
required sharp hole charge for arbitrary components.

## 5. Exact scope boundary

The theorem proves an exact, all-dyadic next-weight value, with a
matching representative. It is not a proof of d=12r for the IBM
family or of its square-column formula. The unrestricted canonical
hiding-mass inequality, arbitrary omega supports, the sharp pure
tower, the opposite cylinder orientation, and winding-mixed torus
classes retain their separate metric obligations. The square-column
attaining family also remains unresolved. No lower bound alone is
promoted to exact IBM distance.

## 6. Review and verification record

The parent audited the cited new branch proofs and replayed their
fixed identity diagnostics. Independent read-only reviews passed
S13, S14, the four terminal-entry pairs, S16 and S17; the older S11
shell exclusion had already passed independent review. A further
review checked the complete branch exhaustion, the q+4 attainer
including q=2, and its use with parity and uniqueness. A separate
consolidation review passed the boundary-dressing identity, the
k<=6 consequence, and the component-count improvement.

The new fixed graph and polynomial checkers each ran in less than
0.002 seconds with peak RSS below 20 MB. The two prescribed
attainer checks at q=2 and q=4 also passed; these were not distance
censuses. The universal proofs do not depend on extrapolating any
finite diagnostic. All three native agents completed this bounded
proof pass. No old headless session was resumed, and no new work
was committed or pushed.
