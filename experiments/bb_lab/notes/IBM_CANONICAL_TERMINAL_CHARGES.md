# Canonical excess two: the terminal-to-alternating-component branches are empty

2026-09-08. Parent analytic continuation of A42 S14. This note removes
four of S14's twelve residual unordered charge pairs without a support,
height, or arbitrary-word search. All x coordinates remain Laurent
integers, not just residues. Independent review is recorded separately
in the main program after it is received.

## 1. Statement and fixed graph notation

Use the original binary cylinder, with dyadic q, finite x support,

    P=x^3+t+t^2, Q=1+(x+x^2)t^3,
    P f+Q g=N_q, |f|+|g|=q+2.

Assume zero aggregation loss and exactly two slack-one edges in the
S7 graph. The S14 edge numbering is fixed by lexicographic source and
choice order. Its relevant four entry/exit types are

| ID | edge | f residues | g residues |
| --- | --- | --- | --- |
| 16 | V=012 -> H_3=003 | empty | {0,1} |
| 19 | W=013 -> H_2=002 | empty | {0,2} |
| 9 | H_2 -> B=100 | {2} | empty |
| 10 | H_3 -> B=100 | {1} | empty |

**Theorem.** No actual solution as above has its two charges drawn
from one of the four unordered pairs

    (16,9), (16,10), (19,9), (19,10).

At q=1 there cannot be two distinct charged rows. The complete q=2
shell is empty by A42 S11 section 7. We prove the assertion for q>=4.
It is not a theorem about arbitrary pairs of slack-one charges or
arbitrary excess. In particular the other eight S14 residual pairs
are not covered by this statement.

Put X=x^3 and c=x+x^2. The exact row equation is

    X f_j+f_(j-1)+f_(j-2)+g_j+c g_(j-3)=1.              (1)

At zero aggregation loss, every occupied residue is one actual
monomial. Pi_0 below retains all powers of X; it never sets X=1.
Write A=000, C=011, D=110, E=112 and F=113, as in S13.

## 2. No D-to-C excursion is possible

None of the 32 slack-one edges has target D. Thus an edge immediately
preceding a zero D->C edge must be zero-slack. A clean five-edge
window around D->C is spatially impossible by S7/S13: the subsequent
A->B monomial would equal both X^-1 and 1. Consequently a charged
edge must occur in one of the next three rows.

If the first following charge is one row later, its source is C;
if two rows later, its source is A; and if three rows later, its
source is U=001 or B. These statements follow from the zero edges
C->A and A->{U,B}. Neither of our charged edges has any of these
sources: their sources are V/W and H_2/H_3. Hence no D->C occurs.
This also handles short cyclic words: q>=4 gives distinct rows up
to the first possible following charge, and no slack-one edge can
serve as the incoming edge to D.

## 3. The return path has no triangle circuit

Let eta=1 for entry 16 and eta=2 for entry 19; put tau=3-eta.
Set K=F for eta=1 and K=E for eta=2. The terminal zero edge is
K->V or K->W. The edge B->K has f residue eta and g residue tau;
the terminal edge has empty f and g residue eta.

Removing D->C, every zero path from the exit's target B back to
the entry's source is precisely

    (B->C->A->B)^n, B->K, K->(V/W),   n>=0.           (2)

Indeed entering D cannot return; the other terminal states cannot
reach the required source. The fixed graph certificate checks this
relevant subgraph, including its only directed circuit.

Place the entry charge at row zero. If n>=1, the last three zero
rows are A->B at -3, B->K at -2, and the terminal edge at -1.
Rows -4 and -5 are C->A and B->C. The predecessor at -6 is either
another A->B or the exit charge. Both have empty g. If H=f_-3,
equation (1) at -3 gives X H=1. At row -1, all terms except H
have nonzero residue eta or vanish, so Pi_0 of (1) gives H=1.
This is impossible. Therefore n=0.

Let h be the number of zero edges in the alternating H_2,H_3
component. The entire cyclic word now consists of

    row 0: entry charge;
    rows 1,...,h: alternating g-only zero edges;
    row h+1: exit charge;
    row h+2: B->K;
    row h+3: terminal edge.                            (3)

Thus q=h+4. Since q>=4 is dyadic, h is even. The exit must therefore
have the same alternating-component source as the entry's target:
16 pairs with 10, and 19 pairs with 9. The other two pairs are
already excluded. In either remaining pair the two f monomials
at h+1,h+2 have residue eta.

## 4. No height at least eight can realize the word

A monomial input M can make 1+c M a monomial only by cancelling
its constant term. The two possibilities are exactly

    M=x^-2 -> 1+c M=x^-1,
    M=x^-1 -> 1+c M=x.                                 (4)

In particular x^-2 is never an output. This is an identity in
the infinite Laurent coordinate, not a bounded exponent test.

For h>=4, the last g-only monomial g_h has residue tau. At row
h+3 the other nonzero terms of (1) all have residue eta. Hence
the residue-zero part of c g_h equals 1. This forces

    g_h=x^-1 if eta=1,    g_h=x^-2 if eta=2.            (5)

Row h is a clean g-only row, with f_h=f_(h-1)=f_(h-2)=0 and
g_(h-3) an actual monomial. Thus g_h=1+c g_(h-3).
For eta=2, this immediately contradicts (4).

For eta=1, equations (4)-(5) force g_(h-3)=x^-2. At q>=16,
h>=12, and the row h-3 equation also has three empty f terms
and a monomial g_(h-6). It would produce x^-2 as an output,
again impossible. The only smaller dyadic possibility q>=8
is q=8, h=4. It requires g_1=x^-2. But row one has empty
f_1,f_0,f_7 and monomial g_6, so g_1=1+c g_6, the same
contradiction. The charged two-monomial g_0 has not been
mistakenly used as a clean monomial in either case.

## 5. The two wrapped height-four identities

At h=0 the four rows are entry, exit, B->K, terminal. Let

    f_1=E_0, f_2=M, g_2=T, g_3=R.

Here E_0,M,R have residue eta and T has residue tau. All are
monomials. The cyclic equations, including every wrapped term,
are

    g_0=1+M,
    X E_0+c T=1,
    X M+E_0+T+c R=1,
    M+E_0+R+c g_0=1.                                  (6)

For eta=1, the second equation's residue-zero and residue-one
parts give T=x^-1 and E_0=x^-2. The third equation gives
R=x^-2 and M=x^-5. Since E_0=R, the last equation reduces to

    (1+x+x^2) M=1+x+x^2.

The Laurent ring is a domain, so M=1, contradicting M=x^-5.

For eta=2, the second equation gives T=x^-2 and E_0=x^-4.
The residue-zero part of the third gives R=x^-1. Its
residue-one part is then x^-2+x=0, impossible. This proves
the theorem for q=4 and completes the all-dyadic exclusion.

## 6. Verification scope

`ibm_canonical_terminal_charges.py` checks the four exact edge
labels, the absence of slack-one edges entering D, the two
relevant no-D->C return subgraphs, the fixed row projections,
and the Laurent substitutions in (4) and (6). These are fixed
identity checks, not a census of possible supports or heights.
The infinite quantifiers are discharged by the arguments above.
