# A42 S14: two unit charges force the isolated g-only excursion

2026-09-08. Bounded Codex continuation. Only this note and
`scripts/a42_s14_codex_double_charge.py` are new outputs.

## 0. Result and exact remaining scope

The sharp inequality `k<=3` at excess two is **not proved here**.
This note considers **only** the branch with zero aggregation loss
and exactly two slack-one edges. It does not treat a single
slack-two edge, assigned separately to S13.

The following are proved for every dyadic height `q=2^a`:

1. Of the 46 unordered graph-compatible charged-edge pairs, **34
   cannot occur in any compact actual syndrome solution** in this
   branch. This eliminates 68 of the parent's 92 ordered pairs.
2. The twelve remaining pairs all enter the isolated g-only
   zero-slack two-cycle and subsequently leave it. They are exactly
   an entry from `{A,O,M,N}` of the types specified in §6 and an
   exit to B. No other two-charge solution is possible.
3. Every actual solution in this branch has **q<=16**. This follows
   from bounds on actual clean monomial words, not from a height
   scan, compression assumption, or finite-width search.

The residual sharp problem is explicitly the twelve pair types
in §6 at `q in {4,8,16}`. Their physical Laurent positions have
not been searched or classified, and their hole bound remains
open. A smaller complete height range is **not** a claimed proof
of that bound or a license for an expensive finite search.

## 1. Setup and graph conventions

Use the original canonical syndrome problem

```
R_q=F2[t]/(t^q+1),
P=x^3+t+t^2,       Q=1+(x+x^2)t^3,
N_q=1+t+...+t^(q-1),
P f+Q g=N_q,       |f|+|g|=q+2.
```

All x supports are finite. Put `X=x^3`; X is still an unbounded
Laurent variable. The actual row equation is

```
X f_j+f_(j-1)+f_(j-2)+g_j+(x+x^2)g_(j-3)=1.          (1.1)
```

The hole support S consists of f-column zero and g-columns zero
and one, and `k=|supp(f,g) intersect S|`.

S7 gives `delta=L+sigma`, where L is physical weight lost under
x-mod-three parity aggregation and sigma is total graph slack.
Here **L=0** and the walk has exactly two slack-one edges, with
all other edges zero slack. Thus each occupied residue on every
row is exactly one actual monomial. The monomial exponents are
not fixed by the graph and are not reduced modulo three in (1.1).

S11 excludes `L=2,sigma=0` at dyadic excess two. That result is
context, not a substitution for any two-charge analysis below.
Height one cannot have two charged rows, and the height-two shell
is empty by S11. The argument here concerns all larger dyadic q.

### 1.1 The full sixteen states

Write `omega^2+omega+1=0` and name the states as follows:

```
A=000, B=100, C=011, D=110, E=11omega, F=11omega^2,
G=00omega, H=00omega^2,
O=001, I=010, J=101, K=10omega, L=10omega^2,
M=01omega, N=01omega^2, Z=111.
```

Each three-character expression is the S7 state
`(a_previous,a_earlier,d_previous)`, not a binary integer.
The checker encodes the F4 values `1,omega,omega^2` as `1,2,3`
and uses XOR for addition.

The **complete** zero-slack edge list is

```
AB, BC, BE, BF, CA, DC, DD, ED, FD,
GH, HG,
AO, IK, IL, JI, EN, FM.                              (1.2)
```

The last six transient edges are retained throughout this note.
In particular `AO`, `JI`, `IK/IL`, and `EN/FM` are essential to
the cases below; they cannot be dropped from a charged closed walk.

### 1.2 All thirty-two slack-one edges

An integer in the support columns denotes one occupied x residue,
not a fixed x coordinate. A dash is empty support.

| ID | edge | f residues | g residues |
|---:|---|---|---|
| 0 | AG | 0,1 | 1 |
| 1 | AH | 0,2 | 2 |
| 2 | AK | 2 | 0,2 |
| 3 | AL | 1 | 0,1 |
| 4 | OG | 0,2 | 1 |
| 5 | OH | 0,1 | 2 |
| 6 | OJ | 0 | 1,2 |
| 7 | OK | 1 | 0,2 |
| 8 | OL | 2 | 0,1 |
| 9 | GB | 2 | - |
| 10 | HB | 1 | - |
| 11 | IA | 1,2 | - |
| 12 | IO | - | 1,2 |
| 13 | CJ | 0 | 0 |
| 14 | CK | 1 | 1 |
| 15 | CL | 2 | 2 |
| 16 | MH | - | 0,1 |
| 17 | MJ | 1 | 0 |
| 18 | MK | 0 | 1 |
| 19 | NG | - | 0,2 |
| 20 | NJ | 2 | 0 |
| 21 | NL | 0 | 2 |
| 22 | BI | 1,2 | - |
| 23 | JZ | 0 | 0 |
| 24 | JE | 1 | 1 |
| 25 | JF | 2 | 2 |
| 26 | KN | - | 0,1 |
| 27 | KZ | 1 | 0 |
| 28 | KE | 0 | 1 |
| 29 | LM | - | 0,2 |
| 30 | LZ | 2 | 0 |
| 31 | LF | 0 | 2 |

For a state `u=(a_prev,a_prev2,d_prev)` and choice `(a,d)`,
the target is `(a,a_prev,d)`. The f Fourier data are
`(a,1+d+d_prev)` and the g data are
`(1+a+a_prev+a_prev2,d)`. Inverting the three-bit Fourier map
gives exactly the table. Slack is

```
|f residues|+|g residues|-1+F(u)-F(v),
F(u)=0 on {A,O,I,B}, and 1 elsewhere.
```

These formulas verify all 128 graph edges and both displayed
lists without relying on a graph discovery heuristic.

### 1.3 The finite compatibility ledger

Two charged edges e and h can lie on such a closed walk only if
the source of h is zero-reachable from the target of e, and
conversely. Empty zero paths are allowed. The fixed graph has
92 ordered compatible pairs, or 46 unordered pairs; none uses
the same charged edge twice. The three edges ending at Z have
no compatible partner.

The exact partition used below is:

| proof group | unordered ID pairs |
|---|---|
| forbidden DC, §3.1 | (17,24), (17,28), (20,25), (20,31), (18,28), (21,31) |
| nonbinary length, §3.1 | (11,17), (11,20), (17,25), (17,31), (20,24), (20,28) |
| four-row contradiction, §3.2 | (17,29), (20,26) |
| OJ, §4.1 | (6,11), (6,12), (6,24), (6,25), (6,28), (6,31) |
| KE/LF, §4.2 | (2,28), (3,31), (7,28), (8,31), (13,28), (13,31), (14,28), (15,31), (22,28), (22,31) |
| IA, §4.3 | (11,13), (11,22) |
| opposite-residue bend, §5 | (13,24), (13,25) |
| unresolved, §6 | entry in {0,1,4,5,16,19}, exit in {9,10} |

The groups are disjoint and have sizes `6,6,2,6,10,2,2,12`.
The checker verifies their union is exactly the 46 compatible
pairs. Compatibility is only a necessary graph condition, not
physical realizability or a dyadic-length assertion.

## 2. A DC transition must encounter a specific charged edge immediately

**Lemma 2.1.** Every DC in an actual walk of this branch has a
first subsequent charged edge within its next three rows. That
edge and its preceding zero suffix have one of these forms:

```
DC, (13 or 14 or 15);
DC, CA, (0 or 1 or 2 or 3);
DC, CA, AO, (4 or 5 or 6 or 7 or 8);
DC, CA, AB, 22.                                      (2.1)
```

The final possibility is itself spatially impossible, so ID 22
can be removed from this list.

**Proof.** No slack-one edge ends at D. Its incoming edge is
therefore one of DD, ED, FD, whose actual g block is empty.
If the next three edges were all zero, they would be CA, AB,
and one of BC, BE, BF. The alternative AO cannot be followed
by another zero edge. S7's two spatial equations would then
force the same f monomial to be both `X^-1` and 1.

The state table now gives exactly (2.1) by following up to three
zero edges from C until the first charge. If that first charge
is BI (ID 22), its two f sites have nonzero residues and its g
block is empty. The very same S7 projection still applies:
the AB row has `X F=1`, while the BI row has `F=1`.
Thus BI cannot repair DC. QED.

Let

```
T={0,1,2,3,4,5,6,7,8,13,14,15}.                    (2.2)
```

If neither charged edge belongs to T, no DC is possible.
If only one belongs to T, at most one DC is possible. Indeed
two DC edges whose first subsequent charge is the same would
both lie in the preceding three rows. But zero DC-to-DC returns
take at least five edges. This is a necessary local-context
argument; it does not discard transient states.

For later use, an **f-monomial pinch** means that the residue-zero
parts of two row equations read `X F=1` and `F=1`, with F the
same actual residue-zero f monomial. Such a pinch is impossible
because x is unbounded. Every use below retains the delayed g
row three positions earlier.

## 3. Charged excursions through the other transient states

### 3.1 DC-free graph obstructions

None of IDs `11,17,18,20,21,24,25,26,28,29,31` belongs to T.
For any pair made from these IDs, DC is consequently absent.
The zero graph then has the ABC cycle, the D loop, the isolated
GH cycle, and the displayed terminal branches.

The six pairs in the first ledger row cannot close without DC:
they would require, respectively, a zero path from E to M,
E to M, F to N, F to N, E to M, or F to N. Such paths do not
exist after deleting DC. A D loop does not help, since it then
has no exit.

For `(11,17)` the entire word, up to rotation, is

```
MJ, JI, IA, (AB,BC,CA)^n, AB, BF, FM,
```

where the charged edges are MJ and IA. Its length is `6+3n`.
The pair `(11,20)` has the analogous word using NJ, BE, EN.
Neither has dyadic length.

The remaining DC-free transient words are forced:

| pair | word | length |
|---|---|---:|
| (17,25) | MJ, JF, FM | 3 |
| (17,31) | MJ, JI, IL, LF, FM | 5 |
| (20,24) | NJ, JE, EN | 3 |
| (20,28) | NJ, JI, IK, KE, EN | 5 |
| (17,29) | MJ, JI, IL, LM | 4 |
| (20,26) | NJ, JI, IK, KN | 4 |

The length-three and length-five cases are not dyadic. The
two length-four cases require actual polynomial analysis next.
All these path assertions follow from (1.2); the checker also
verifies the finite path lengths, or their residue modulo three.

### 3.2 The two exceptional four-row words

Number each displayed four-row word starting at row zero.

For `MJ,JI,IL,LM`, write `f_0=U`, `f_2=V`; both have residue
one. Row one reads

```
U+(x+x^2)g_2=1.
```

Here g_2 has residue two. The equation forces `g_2=x^-1`
and `U=x`. Row zero then forces `g_0=1` and `V=XU=x^4`.
But row three's residue-one equation is `V+x g_0=0`, forcing
`V=x`, a contradiction.

For `NJ,JI,IK,KN`, the analogous U,V have residue two and
g_2 has residue one. Row one forces `g_2=x^-2`, `U=x^-1`.
Row zero forces `g_0=1`, `V=x^2`. Row three gives `g_3=1+x`.
The remaining row-two equation would then be

```
x^5+x^-1+x^-2+x+x^3=1,
```

which is false in the Laurent domain. This separately checks
both residue choices; no unproved symmetry exchanging the
actual x-residue geometry is used.

## 4. Monomial pinches exclude eighteen further pairs

### 4.1 Every compatible pair containing OJ is excluded

ID 6 has partners `11,12,24,25,28,31`.

The pair `(6,12)` has the unique word `OJ,JI,IO`, of length
three, and is not dyadic.

For partners 24,25,28,31, the zero return from E or F to O
requires DC. Only ID 6 belongs to T, so Lemma 2.1 forces the
unique prefix `DC,CA,AO,OJ`. The following edge is JI for
partners 28,31, and JE or JF for partners 24,25. On the OJ row
let F be its residue-zero f monomial. That row has `X F=1`
after residue-zero projection. The following row has `F=1`:
JI is empty, while JE/JF have only nonzero residues. The
delayed g rows are respectively DC (residue zero) and CA
(empty), so neither contributes to these projections.

For `(6,11)`, the same argument applies if DC occurs. If no
DC occurs, the only zero return from A to O consists of n ABC
cycles followed by AO. The total length is `4+3n`. Two ABC
cycles already contain the forbidden six-row S8 word, so
`n<=1`. The only dyadic possibility is the word

```
AO, OJ, JI, IA                                      (q=4).
```

In this word the OJ row's residue-zero equation gives `X F=1`
and the JI row gives `F=1`. The delayed g rows, including the
wrapped ones, are empty. This eliminates the last OJ case.

### 4.2 Every remaining pair containing KE or LF is excluded

After §3 and §4.1, the KE partners still needing treatment are
`2,7,13,14,22`; the LF partners are `3,8,13,15,22`.
All corresponding zero returns pass through D and require DC.
Lemma 2.1 therefore forces the following words locally:

| charged pair | forced local word | pinched f row |
|---|---|---|
| (2,28) | DC, CA, AK, KE, ED | KE |
| (3,31) | DC, CA, AL, LF, FD | LF |
| (7,28) | DC, CA, AO, OK, KE, ED | KE |
| (8,31) | DC, CA, AO, OL, LF, FD | LF |
| (14,28) | DC, CK, KE, ED | KE |
| (15,31) | DC, CL, LF, FD | LF |
| (13,28) | CJ, JI, IK, KE, ED | KE |
| (13,31) | CJ, JI, IL, LF, FD | LF |

In every line the pinched row gives `X F=1`, and the following
ED/FD row gives `F=1`, after residue-zero projection of (1.1).
All other listed f sites have nonzero residues or are empty.
The delayed g sites are empty or in residue zero, and hence
their shifts by x and x^2 have no residue-zero component.
For CK/CL the needed incoming-D predecessor is DD, ED, or FD;
its g block is empty. This accounts for the left endpoint too.

For the two pairs containing BI (ID 22), DC would have ID 22
as its first following charge. Lemma 2.1 already excludes this.
Thus these ten pairs are impossible.

### 4.3 The two remaining IA pairs

For `(11,13)`, the local word is

```
CJ, JI, IA, AB, (BC or BE or BF).
```

IA must be followed by AB: AO is a terminal edge and both
charges have already been used. The AB f monomial F is
pinched. Its equation gives `X F=1`, since the delayed g on
CJ has residue zero. The following equation gives `F=1`,
since the delayed g on JI is empty.

For `(11,22)`, BI and IA must be adjacent. The incoming edge
to B before BI is AB. The outgoing A edge after IA is also
AB. Unless the entire word has the nondyadic length three,
the next B edge is BC, BE, or BF. Hence the local word is

```
AB, BI, IA, AB, (BC or BE or BF).
```

Again the second AB f monomial is pinched: the two delayed
g rows are the first AB and BI, both empty.

The checker verifies 25 exact residue-zero projection contexts
covering these pinches, the three possible incoming-D edges,
the three outgoing B edges, and the cyclic q=4 OJ case.
These identities concern actual monomials, not aggregated
equations that would erase the distinction between X and 1.

## 5. The last two nongeneric bends are impossible as well

Consider `(13,24)` or `(13,25)`. The zero return from E or F
to C requires DC, and only ID 13 belongs to T. Lemma 2.1
therefore makes the entire cyclic word

```
DC, CJ, JE, ED, DD^r,
or DC, CJ, JF, FD, DD^r,       q=4+r.                (5.1)
```

Let the DC g monomial be U and the CJ g monomial be A; both
have residue zero. Let the JE/JF f monomial be M of residue
rho, and the ED/FD f monomial be T of residue `sigma=3-rho`.
Thus `(rho,sigma)` is `(1,2)` or `(2,1)`.

The nonzero-residue part of the ED/FD equation is

```
X T+M+(x+x^2)U=0.                                  (5.2)
```

The next row is either DD or the closing DC. Both have only
residue-zero sites, so its nonzero-residue equation is

```
T+M+(x+x^2)A=0.                                    (5.3)
```

The rho components give `M=x^rho U=x^rho A`, hence U=A.
The sigma components then give `X T=T`, impossible for a
nonzero Laurent monomial. This treats r=0 as well: at height
four the next row is the closing DC, whose residue-zero g
site contributes nothing to (5.3). No extra D-loop row or
unavailable endpoint equation is assumed.

Together §§3–5 eliminate all 34 pairs in the first seven ledger
rows. This completes the proof of result 1 in §0.

## 6. Exact residual types and an all-height bound of sixteen

**Theorem 6.1.** Any actual dyadic solution in the two-slack-one,
L=0 branch uses an entry from

```
AG, AH, OG, OH, MH, NG       (IDs 0,1,4,5,16,19)
```

followed, after zero or more alternating GH/HG edges, by an
exit GB or HB (ID 9 or 10). It has q<=16.

**Proof of the type statement.** It is exactly the complement
of the 34 excluded pairs in the finite compatibility ledger.
The zero component {G,H} has no zero exit or entrance. Thus
the entry and exit order, and the alternating intervening path,
are forced. There is only one such visit because the two
charged edges have both been spent.

**Proof of the height bound.** Let m count the zero GH/HG
edges between the two charged rows. These m actual rows have
f=0 and a single g monomial. S8's local g-only lemma applies
to that interval and gives `m<=9`. Its endpoint argument uses
only the recurrence after the first three rows; the charged
entry and exit are not treated as clean monomial rows.

The zero return arc starts at B, the target of the charged
exit. For MH or NG, neither charge belongs to T, so there is
no DC. The return arc must be

```
(BC,CA,AB)^n, BF, FM       for MH,
(BC,CA,AB)^n, BE, EN       for NG.
```

The ABC prefix has `3n` edges. Any eight clean ABC edges
contain S8's forbidden six-edge word starting at A. Therefore
`n<=2`, and these cases have

```
q=4+m+3n<=19.                                      (6.1)
```

For AG/AH or OG/OH there can be at most one DC: the only
eligible first charge is the entry. If DC is absent, the
return arcs are, respectively,

```
(BC,CA,AB)^n, BC, CA;
(BC,CA,AB)^n, BC, CA, AO.
```

Again `n<=2` suffices, giving `q<=19` or `q<=20`.

If DC is present, its first charge must be the entry, at
distance two for AG/AH or distance three for OG/OH. Therefore
the complete zero return arc is

```
(BC,CA,AB)^n, BE, ED, DD^r, DC, CA, [AO],
```

or the version using BF,FD. The bracketed AO occurs precisely
for an O entry. There is no ABC copy after DC: it would put
the first charge outside the three-row window of Lemma 2.1.

Here `n<=2` as before. Also `r<=7`. To see the latter, suppose
r>=1 and number DC as row zero. On the following CA row,
(1.1) reads `f_-1=1`: f_0, f_1, g_1 are empty, and g_-2 is
empty because it is a DD row or the entering ED/FD row. Thus
the final D-loop monomial is noncanonical. The entire r-entry
D word obeys the S8 monomial recurrence from its third entry,
since its entering ED/FD g row is empty. S8's lemma forces any
word of length at least eight to be wholly canonical, a
contradiction. The case r=0 needs no endpoint claim.

Consequently the two possible heights with a DC satisfy

```
q=6+m+3n+r <=28       for AG/AH,
q=7+m+3n+r <=29       for OG/OH.                    (6.2)
```

In all cases q<=29. Since q is a power of two, this implies
q<=16. QED.

This bounds **actual solutions in this branch**, not merely
counterexamples after a compression operation. No physical
positions or configurations at heights 4,8,16 were enumerated.

### 6.1 The precise remaining lemma

For each of the twelve entry/exit pairs above, and every actual
spatial realization of its zero return arc satisfying (1.1),
prove `k<=3`. By Theorem 6.1 it suffices to prove this at
q=4,8,16, with arbitrary finite Laurent x positions and the
required charged-row residues. The entry row can contain two
f monomials or a two-site g block; it must not be included
in the monomial g-only interval when using S8.

That is the unresolved sharp assertion, not a theorem supplied
by the graph partition or the height bound. The single-slack-two
branch remains outside this note's assignment.

## 7. Verification and resource scope

The standalone checker rebuilds the fixed 16-state graph,
checks the 92 ordered/46 unordered pair ledger, verifies the
first-charge prefixes after DC, and checks the relevant fixed
path lengths or lengths modulo three. Its path routines reject
cycles when an exact finite path-length assertion is requested;
they do not truncate walks and call the truncation exhaustive.

It additionally verifies the explicit four-row Laurent constants,
25 formal residue-zero pinch contexts, and both opposite-residue
bends, including their height-four cyclic endpoint. The universal
loop bounds and resulting q<=16 theorem are the analytic proofs
above. The diagnostic explicitly reports twelve unresolved pairs;
its PASS is not a claim to have proved the requested sharp bound.

Run from the repository root with

```
PYTHONDONTWRITEBYTECODE=1 uv run --offline --no-sync --no-cache \
  --project experiments/bb_lab python \
  experiments/bb_lab/scripts/a42_s14_codex_double_charge.py
```

All development invocations are fixed-graph/finite-identity checks,
not support or height searches. They have remained well below one
second in aggregate diagnostic process time and below 20 MiB RSS.
The script enforces 55-second wall and 50-second CPU limits and
asserts the 768 MiB RSS budget. No shared files, git state, memory,
Lean, old processes, or other agents were modified or controlled.
