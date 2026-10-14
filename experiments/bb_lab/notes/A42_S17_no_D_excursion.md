# A42 S17 — all eight remaining entry/exit pairs are empty without DC

2026-09-08. Bounded native Codex subproof. Original binary cylinder,
finite Laurent x support, no support or height search. Only this
note and `scripts/a42_s17_codex_no_D.py` are new.

## 1. Theorem and scope

Let q=2^a with a>=0, and work over

    R_q=(F2[t]/(t^q+1))[x^+/-1],
    P=x^3+t+t^2,       Q=1+(x+x^2)t^3,
    Pf+Qg=N_q,        N_q=1+t+...+t^(q-1),
    |f|+|g|=q+2.

Weight is physical binary Hamming support. Assume:

1. Aggregation of x exponents modulo three loses no physical
   weight, separately on every row and block.
2. The S7 cyclic state walk has exactly two slack-one edges, all
   other edges having zero slack.
3. Its two charged edge IDs, in the S14 convention recalled below,
   are one entry from {0,1,4,5} and one exit from {9,10}.
4. The walk contains no zero edge D->C.

**Theorem 1.1. There is no such actual solution, at any dyadic q.**

The eight unordered pair types covered are

    (0,9), (0,10), (1,9), (1,10),
    (4,9), (4,10), (5,9), (5,10),

always subject to the no-DC hypothesis. This hypothesis is not
derived or dropped. The DC branch is assigned elsewhere. Thus
this note does not by itself close the full next shell, canonical
HM, the pure omega floor, or the original IBM distance theorem.

The new local ingredient is stronger than the initially suggested
long g-only recurrence: two equations at the exit and the following
CA row already force the same nonzero monomial E to satisfy XE=E,
where X=x^3 remains an unbounded Laurent variable. It excludes
every sufficiently long alternating excursion, including the
apparently exceptional height-eight words. The two height-four
words have a separate residue-zero contradiction.

## 2. Fixed graph input and physical row convention

Set

    X=x^3,        c=x+x^2.

The exact cyclic row equation is

    X f_j+f_(j-1)+f_(j-2)+g_j+c g_(j-3)=1.            (2.1)

Pi_0 projects onto x exponents divisible by three, retaining their
integer values. In particular X is not set to one. Under assumption
1, every occupied aggregated residue contains exactly one physical
monomial, and every empty residue contains none.

Use the states

    A=000, O=001, B=100, C=011, D=110,
    E=112, F=113, H_2=002, H_3=003.

Digits in the third coordinate denote F4 elements 0,1,omega,omega^2.
For state `(a_prev,a_older,d_prev)` and choice `(a,d)`, the target
is `(a,a_prev,d)`. The f and g Fourier coordinates are respectively

    (a,1+d+d_prev),       (1+a+a_prev+a_older,d).

For a binary triple, the inverse Fourier supports are

    (0,0): empty;        (1,0): {0,1,2};
    (1,1): {0};          (1,2): {1};          (1,3): {2};
    (0,1): {1,2};        (0,2): {0,2};        (0,3): {0,1}.

Slack is the sum of these two support sizes minus one, plus the
source potential minus the target potential. The potential is zero
on {000,001,010,100} and one elsewhere. These fixed formulas give
all 128 edges and the following complete zero-slack adjacency:

    A -> O,B;
    H_2 -> H_3;          H_3 -> H_2;
    010 -> 102,103;      101 -> 010;
    C -> A;             B -> C,E,F;
    D -> C,D;           E -> 013,D;          F -> 012,D.

States O,012,013,102,103,111 have no zero outgoing edge. In
particular all transient edges are retained in this description.

The relevant charged edges, using lexicographic source and choice
order among the slack-one edges, are exactly

| ID | edge | f residues | g residues |
| --- | --- | --- | --- |
| 0 | A -> H_2 | {0,1} | {1} |
| 1 | A -> H_3 | {0,2} | {2} |
| 4 | O -> H_2 | {0,2} | {1} |
| 5 | O -> H_3 | {0,1} | {2} |
| 9 | H_2 -> B | {2} | empty |
| 10 | H_3 -> B | {1} | empty |

The zero edges used spatially below have supports

| edge | f residues | g residues |
| --- | --- | --- |
| AB | {0} | empty |
| BC | empty | {1,2} |
| CA | empty | empty |
| AO | empty | {0} |
| H_2 -> H_3 | empty | {2} |
| H_3 -> H_2 | empty | {1} |

These lists follow by direct substitution into the displayed
Fourier formulas. No claim about actual monomial exponents follows
from the graph alone.

## 3. Complete return path and the triangle pinch

Place the entry charge at row zero. Let h>=0 be the number of zero
edges in the alternating H_2,H_3 component before the unique exit
charge. The zero path from the exit's target B back to the entry's
source is exactly

    (B->C->A->B)^n, B->C, C->A, [A->O],     n>=0,      (3.1)

where the bracketed edge occurs precisely for entries 4 and 5.

Proof of completeness. After deleting DC, D has only its self-loop.
A zero path entering D cannot return to A or O. From B, entering
E or F therefore either gets trapped at D or reaches a different
terminal state. Neither can reach the required charged source A
or O. Consequently every successful return must use BC and CA.
At A it either stops, takes the final AO, or takes AB and repeats.
There is no zero route to the alternating component from this
return path. This proves (3.1), including arbitrary repetitions.
It also shows that a charged closed walk under these hypotheses
cannot visit D at all.

Counting edges before any spatial argument gives

    q=h+4+3n       for entries 0,1,
    q=h+5+3n       for entries 4,5.                    (3.2)

In particular q=1 and q=2 are already impossible here; no base-shell
certificate is needed for this restricted theorem.

**Lemma 3.1. Every actual return has n=0.**

Proof. If n>=1, let row r be the final AB before the final BC,CA.
The six edges at rows r-3 through r+2 are

    (AB or exit), BC, CA, AB, BC, CA.

The predecessor at r-3 has empty g in either case. If the AB
monomial at r is M, equation (2.1) at r is therefore

    X M=1,

since f_(r-1), f_(r-2), g_r, and g_(r-3) all vanish. At r+2,
the same equation is

    M=1,

since the current edge is CA, the preceding edge is BC, and the
delayed g term is on the preceding CA at r-1. The optional AO
comes later and changes neither equation. This is impossible in
the Laurent coordinate. These six rows have no short-period
alias: (3.2) with n>=1 gives q>=7. QED.

Thus the entire word has rows

    0: entry;
    1,...,h: alternating g-only zero edges;
    h+1: exit;
    h+2: BC;
    h+3: CA;
    h+4: AO only for entries 4,5.                     (3.3)

For entries 0,1, q=h+4 is even, so h is even. For entries 4,5,
q=h+5 is even, so h is odd. The alternating component therefore
fixes the only possible exit as follows:

| entry | required exit | h parity | exit f residue rho | final alternating g residue tau |
| --- | --- | --- | --- | --- |
| 0 | 9 | even | 2 | 1 |
| 1 | 10 | even | 1 | 2 |
| 4 | 10 | odd | 1 | 2 |
| 5 | 9 | odd | 2 | 1 |

The final column is used only when h>=1. The other four entry/exit
pairings are excluded by this parity argument. At dyadic q>=8,
h>=4 for A entries and h>=3 for O entries.

## 4. A two-equation spatial endpoint obstruction

**Lemma 4.1.** Let U,V be Laurent monomials with the same nonzero
x residue tau in {1,2}. Let E be a nonzero Laurent monomial of
residue rho=3-tau. The two equations

    X E+c U=1,       E+c V=1                          (4.1)

cannot both hold over F2[x^+/-1].

Proof. The residue-zero parts of XE and E vanish. In c times a
monomial of residue tau, exactly one term has residue zero: the
shift by `3-tau`. Thus Pi_0 applied to the two equations forces

    U=V=x^(tau-3).

Adding the two equations now gives `(X+1)E=0`. The Laurent ring
is a domain, X+1 is nonzero, and E is nonzero, a contradiction.
Explicitly, the incompatible forced values are

| tau | U=V | E from CA equation | E from exit equation |
| --- | --- | --- | --- |
| 1 | x^-2 | x^-1 | x^-4 |
| 2 | x^-1 | x | x^-2 |

No bound or enumeration of Laurent exponents is involved. QED.

**Corollary 4.2. The word (3.3) is impossible whenever h>=3.**

Proof. Let E=f_(h+1) be the exit monomial. At the exit row the
previous two f blocks are empty, and the exit g block is empty.
The actual equation is

    X E+c g_(h-2)=1.                                  (4.2)

At the CA row two rows later, the actual equation is

    E+c g_h=1.                                        (4.3)

Both g_(h-2) and g_h belong to the alternating zero-edge segment:
their indices are at least one. They are actual monomials in the
same nonzero residue because their indices differ by two. The exit
residue is the opposite nonzero residue, by the fixed graph table.
Lemma 4.1 applies. QED.

This uniformly handles every dyadic q>=8, for both A and O entries.
It is unnecessary to propagate a three-step g-only recurrence or
separately bound the length of that recurrence.

### Explicit height-eight endpoints

The supposedly exceptional short words use the following exact
indices, with no unwrapping of cyclic row equations:

| entry/exit | exit row | CA row | U | V |
| --- | --- | --- | --- | --- |
| 0/9 or 1/10 | 5 | 7 | g_2 | g_4 |
| 4/10 or 5/9 | 4 | 6 | g_1 | g_3 |

For A entries the exit equation has empty f_4,f_3 and delayed
monomial g_2. For O entries it has empty f_3,f_2 and delayed
monomial g_1. The CA equations have exactly the matching terms
listed above. No charged g_0 is used in either equation, and no
charged f_0 is incorrectly declared empty.

For O entries, the final AO equation does indeed force g_7=1 at
q=8, and more generally g_(h+4)=1. It is not needed for the
contradiction. The argument also needs no equation at row h-3,
avoiding the earlier recurrence's short-height endpoint issue.

## 5. The two height-four A-entry words

Only A entries can have q=4, and (3.2) and Lemma 3.1 force h=0.
Their entire cyclic words are

    (entry 0), (exit 9), BC, CA,
    (entry 1), (exit 10), BC, CA.

Write the entry f block as `f_0=A_0+M`, where A_0 is its single
residue-zero monomial and M has nonzero residue. The entry g_0
has nonzero residue, and the exit f_1 also has nonzero residue.
The exit g_1 is empty. Row two is BC and row three is CA.

The residue-zero part of row zero of the actual cyclic equation is

    X A_0=1:                                         (5.1)

f_3 and f_2 are empty, g_0 has no zero-residue part, and the wrapped
delayed term is c g_1=0. At row two, the residue-zero equation is

    A_0=1:                                           (5.2)

f_2 is empty, f_1 has nonzero residue, g_2 has only residues 1,2,
and its wrapped delayed term is c g_3=0. Equations (5.1)-(5.2)
contradict X!=1. This covers both residue choices without any
reflection or axis-equivalence assumption.

Together, the minimum graph lengths, the triangle pinch, the
exit-parity restriction, Corollary 4.2, and this wrapped identity
prove Theorem 1.1.

## 6. Audit of the supplied reduction and exact program boundary

The parent's no-DC return paths and n>=1 pinch are correct. Their
period counts and forced exit IDs are also correct. The stated
last-g values follow from (4.3), and the O-entry final AO monomial
is indeed one. Their indicated q=4/q=8 residual word shapes are
valid necessary graph shapes, not evidence of actual realizability.

The new improvement is to compare the exit row directly with CA.
This proves the stronger local obstruction (4.1) and removes every
longer height in one step, including all indicated q=8 cases.
There is no residual within this note's no-DC branch.

What remains outside its scope is an actual two-slack-one walk
containing DC. Neither the absence of charged edges entering D nor
the work here implies that DC is absent for entry sources A or O.
No assertion about that separately assigned branch is made.

## 7. Fixed diagnostic and resource record

The new diagnostic reads the previously inspected graph and
polynomial helpers from S14, with bytecode writing disabled. It
checks the complete 17-edge zero graph and all six relevant charge
IDs and residue supports; the three predecessor cases of the
triangle pinch; the two Laurent endpoint forcings; both prescribed
height-four wrapped words; and the four prescribed height-eight
endpoint contexts. The latter are fixed templates from the proof,
not an enumeration of candidate words or supports.

The universal theorem is proved in §§3-5, not inferred from these
checks. No rational-function elimination, support search, exponent
box, height scan, or growing matrix was needed.

One diagnostic execution passed in 0.000330 seconds with measured
peak RSS 18,546,688 bytes. It installs a 55-second wall alarm and
50-second CPU limit and checks the 768 MiB peak-RSS cap. Only the
two assigned new S17 files were written; no sibling note or script,
shared file, git state, Lean file, memory, or old process was changed.
