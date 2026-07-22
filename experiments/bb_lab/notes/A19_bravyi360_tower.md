# A19 — Bravyi [[360,12,≤24]] as a deck-nontrivial tower: instance study and path to certification

**Session 2026-07-20.** Origin: methodology review of IBM's LLM-guided BB-code
discovery paper (arXiv:2606.02418) surfaced Bravyi et al.'s [[360,12,≤24]]
(arXiv:2308.07915, Table 3: `(ℓ,m) = (30,6)`, `A = x⁹+y+y²`,
`B = y³+x²⁵+x²⁶`) as the natural gross-playbook target — highest published
FOM, distance unresolved (IBM's deep MILP: incumbent only; "≤" in Bravyi's
table = upper bound known only). This note records the instance study, the
corrected structural picture, and the staged plan toward certified distance.

Scripts: `scripts/a19_*.py` (run from `experiments/bb_lab` via `uv run`).
Data: `data/a19/` (ISD per-class minima npz, SAT round logs).

## 1. The tower and its exact base distances

```
GB (15,3)  --x-->  BY (30,3)          A = x⁹+y+y² throughout
   |                  |               B descends by exponent reduction:
   y                  y               BX: y³+x¹⁰+x¹¹   BY: 1+x²⁵+x²⁶
   v                  v               GB: 1+x¹⁰+x¹¹
BX (15,6)  --x-->  C  (30,6) = the Bravyi code
```

| code | frame | (n,k) | d | provenance |
|---|---|---|---|---|
| GB | Z₁₅×Z₃ | [[90,8]] | **8 exact** | coset-SAT ladder, 35 orbit reps, 187 s (`a15_coset_distance`); witness support `[12,14,27,29,69,71,84,86]` |
| BY | Z₃₀×Z₃ | [[180,8]] | **12 exact** | floor@11 all-UNSAT (35 reps, 3.1 h, `data/a19/by_floor11.log`) + weight-12 witness (`by_floor13.log`, class 0x10) |
| BX | Z₁₅×Z₆ | [[180,8]] | **12 pending** (≤ 12 witness via ISD; floor@11 round in flight) | `a19_lifts.py` ISD; floor run same protocol as BY |
| C | Z₃₀×Z₆ | [[360,12]] | ≤ 24 constructive; ≥ 6 certified | §3, §4 |

- GB is **not** Bravyi's bb_90 ([[90,8,10]] shares the frame and A; different
  B, different d). BX ≁ BY under the full Aut+monomial orbit
  (`a19_bx_by_equiv.py`: distinct canonical pairs).
- Both lower rungs are k-preserving (8→8) with **deficit 4**:
  d = 12 = 2·8 − 4 on *both* axes — two new natural deficit cells for the A17
  docket, again strictly below the wall value 2d−2, and the first
  same-base/different-axis pair where both axes are deficit-equal.
- All three bases pass `a15_class_certify` in full (D1, D2, (iii), FRM, ANN,
  FLR; one-sided engine floors 20/20/10) ⇒ certified d ≥ 6 each, and by the
  A16 transfer clause **the Bravyi cover inherits certified d ≥ 6** — its
  first certified lower bound by any method.

## 2. Deck-nontriviality (the headline structural verdict)

k jumps 8 → 12 into the cover. Checked exhaustively (`a19_orbit_descents.py`,
`a19_diagonal_descent.py`): all **three** decks (x¹⁵, y³, and the diagonal
x¹⁵y³ via `φ(u,v) = (u+15v mod 30, v mod 3)`), each across the full monomial
orbit (units a ∈ (Z/30)*, b ∈ (Z/6)*) — **no k-preserving descent exists**.
By A12 ((R) ⟺ k-preservation ⟺ Bezout) and A10-L1 (descent verdicts are
code-level), the Bravyi flagship is **genuinely deck-nontrivial on every
deck**: the first natural flagship instance outside the doubling template's
condition 2, after 157/157 historically checked covers satisfied it.

Quantitative structure (`a19_cover_census.py`, `a19_seam_sector.py`, all
kernel-level F₂ linear algebra):

- `dim (1+σ)H₁ = 4 = k̃ − k` for **all three decks** — the A12 open
  quantitative conjecture verified on a natural instance, three ways.
- Pushforward class maps to BX and BY each have rank 6 (kernel 6); the
  doubly-new sector `ker p_x* ∩ ker p_y*` has dim exactly 4, and it **equals
  im(1+σ)₊ for every deck** — the A13 D-module summand, deck-independently.
  Sector census over 155 translation-orbit reps: 138 old / 6 new-x / 6 new-y /
  5 new-xy.

## 3. The weight table is transfer-lift arithmetic (d ≤ 24 constructive)

Deep ISD over ker H_Z (54k iterations, all 4095 classes, vectors saved in
`data/a19/isd_class_minima.npz`):

| sector | min wt seen | witness class | mechanism (verified) |
|---|---|---|---|
| old | 24 | 0x6e2 | section-type: projects to nontrivial base logicals wt 24 (BX) / 22 (BY) |
| new-x | 24 | 0x41d | **τ_x(BX weight-12 logical)** — verified explicitly (0x570) |
| new-y | 24 | 0x002 | **τ_y(BY weight-12 logical)** — verified explicitly (0x380) |
| new-xy | 32 | 0x80c | **τ_x∘τ_y(GB weight-8 logical)** — verified explicitly (0xc04) |

(`a19_lifts.py`: each τ-lift checked to be a cycle, nontrivial, of weight
exactly 2·(base weight), landing in the predicted sector.) The transfer lift
τ = full-preimage projects to **zero** (sheets cancel mod 2) — τ-lifts live in
the *new* sector of their deck, `im τ_* ⊆ ker p_*`. So the cover's light
operators are the tower's arithmetic: 24 = 2·12 (both decks), 32 = 4·8.
**d(C) ≤ 24 is now constructive** (explicit verified operator), independent of
solver incumbents. A budgeted beat-24 SAT gate (w ≤ 22 per-coset over the ten
lightest orbit reps, `a19_beat24.py`) is in flight; 54k ISD iterations never
beat 24. **Epistemics of the gate:** n = 360 SAT is witness-side only — if the
true coset minimum is 24, a w ≤ 22 query is secretly UNSAT and times out
uninformatively (per-coset UNSAT is priced out beyond ~w = 13 at n = 180,
~w = 14 at n = 300; hopeless at n = 360). Hits are decisive, timeouts are NOT
evidence. The *real* beat-24 instrument is at base scale and coincides with
the certification program: a dangerous-sector operator below 24 forces a
base stabilizer with |b| + 2m(b) < 24 (slice identity, (R)-free) — so the
(M)@24 census is simultaneously the lower-bound machine and the localized
sub-24 hunt; old-sector candidates reduce to merges of SAT-enumerable
weight-12–16 base logicals (n = 180, minutes). Blind n = 360 SAT is retired
from the program after this gate's budget expires.

## 4. Gap analysis for certified lower bounds (the corrected picture)

Audit of `docs/gross-distance-proof.md`: **(R) is consumed only by the
safe-sector confinement (Theorem D)**; the dangerous-sector machinery
(slice identity + light-stabilizer inequality (M), Theorem C / Prop 10) never
invokes it. Consequences for this (R)-failing cover:

- **Portable now:** for any cover logical v with `p_y*[v] = 0` (the new-y and
  new-xy sectors), the slice identity over BY applies verbatim:
  `|v| ≥ |b| + 2m(b)`, `b = p_y(v) ∈ Stab_Z(BY) ∪ {0}`, plus the zero-rung
  for b = 0. Symmetrically over BX.
- **Milestone M12 (near-term):** (M)@12 for BY — every nonzero
  `b ∈ Stab_Z(BY)` with `|b| ≤ 11` has `|b| + 2m(b) ≥ 12` — is *gross Prop-10
  depth* (classification to weight 11: hexagons |A|+|B| = 6, D-pairs, …), at
  n = 180. Combined with the projection floor min(d(BX), d(BY)) = 12 on the
  other three sectors, it yields **certified d(C) ≥ 12** — double the A16
  floor, and stronger than any published bound for this code.
- **Full-24 program:** (i) (M)@24 for BY/BX (classification to weight 23 —
  substantially deeper than gross's 11; SAT-assisted census first);
  (ii) zero-rung at 24 = 2·d(base) (τ-tightness: the dangerous floor is
  *attained* by the §3 lifts, exactly the gross pattern); (iii) the genuinely
  new theory: **doubly-old (safe) sector ≥ 24 without (R)** — candidate route
  via A13's D-module decomposition (the seamC ↔ δ₂ transport gap) replacing
  Theorem-D confinement, on F₁₆ components (e = 15; "floors port,
  classifications don't"). A cheap pilot first: doubly-old classes have both
  projections nontrivial (≥ 12 each) — test whether a joint-support argument
  clears 16–18 before building the F₁₆ slot-frame analog.

Record context: the published solver-exact BB record is Bravyi's
[[288,12,18]] (MIP-exact per their Table 3 note; independently E-status in
IBM's Table V). Full certification at [[360,12,24]] would exceed it by 25% in
n and 33% in d; the kernel-verified record (gross, n = 144) by 2.5× in n.

## 5. Falsified-claims ledger (session-internal; doctrine: record them)

- "The 360 code is a doubling tower with d = 4·d(GB)" — dead: (R) fails at
  the top rung, all decks, all presentations.
- "d(BY) = 14" (first ladder witness) — dead: true d(BY) = 12; the witness
  was non-minimal. Same error nearly repeated at BX (MC said 16; ISD found 12
  in 60 s). **Protocol lesson: exhaust cheap witness hunts (ISD with
  early-exit) before buying UNSAT rounds; MC-phase-A of `a15_coset_distance`
  is far weaker than 60 s of ISD.**
- "The cover minimum is an old-sector diagonal lift at 2·d(BY)" — mis-framed:
  τ-lifts project to zero and live in the *new* sectors (§3); the old-sector
  24s are section-type. (Conflated section lift s(u) — weight-preserving,
  p∘s = id, not a cycle in general — with transfer lift τ(u) — weight-2|u|,
  p∘τ = 0, always a cycle.)
- "[[288,12,18]] sits at the wall (its (6,12) base has d = 10)" — dead
  (2026-07-20 addendum): floor@9 + floor@11 all-UNSAT over all 155 orbit
  reps certify d(base) ≥ 12 (§7 verdict).
- "IBM's largest exact is n = 288 at d = 12" — wrong: that is their largest
  *own-discovered* exact; the published record incl. Bravyi's codes is
  [[288,12,18]].

## 6. State and next steps

In flight at session close: BX floor@11 (expect d(BX) = 12 exact); beat-24
SAT gate (3 h budget). Next session, in order:

1. Close BX (floor@11 verdict + witness); log the two deficit cells in the
   A17 docket.
2. **M12**: light-stabilizer classification for BY to weight 11 (Prop-10
   analog; start from the hexagon/D-pair census script pattern) → certified
   d(C) ≥ 12. Then the same at BX for redundancy.
3. Lean staging (engineering track): generated BaseFloors certificates for
   GB/BX/BY (`gen_base_floor_lean.py`; all frames floor-bearing) — kernel
   d ≥ 6 + transfer naming the Bravyi flagship; named-hypothesis packaging of
   d(GB) = 8, d(BY) = 12 on the Z5Z15F2A6 model. Verify BBDoubling's
   dangerous-side rungs are (R)-free (expected from §4) and factor
   `DeckTrivialOnH1` out of them if so.
4. (M)@24 census + the doubly-old pilot (§4) — the moonshot core.

Open question worth naming: both rungs of this tower are deficit-4, yet the
cover recovers 24 = 2·12 via τ-lifts of the *deficit* value while gaining
k̃ − k = 4 logical qubits that sit heavier (32 = 4·8). Whether
"deficit rungs + deck-nontrivial top" is *why* this code's FOM leads the
catalog — i.e., a design principle (hunt covers that break (R) over certified
bases) — is testable against the IBM catalog's other high-FOM codes.

## 7. Same-session addendum: deck survey of the full Bravyi table

`a19_deck_survey.py` (quotients built by coset convolution — the φ-map trick
does not generalize to even m/2): literal-descent k on every order-2 deck of
every published Bravyi code, orbit-swept where a jump appeared.

| code | decks | verdict |
|---|---|---|
| [[72,12,6]] (6,6) | x, y, xy | **JUMP 8→12 on all three — code-level** (full orbit swept: no k-preserving descent exists) |
| [[90,8,10]] (15,3) | — | \|G\| odd: not a Z₂-cover at all |
| [[108,8,10]] (9,6) | y | R-holds (8→8) |
| [[144,12,12]] (12,6) | x / y / xy | x: R-holds (bb72, the known descent); y and xy: JUMP 8→12 (literal) |
| [[288,12,18]] (12,12) | x, y, xy | **R-holds on all three** (12→12) |
| [[360,12,≤24]] (30,6) | x, y, xy | JUMP 8→12 all three, code-level (§2) |
| [[756,16,≤34]] (21,18) | y | R-holds (16→16; base = [[378,16]]) |

Consequences:

- **The strong design-principle hypothesis is falsified**: [[288,12,18]]
  (FOM 13.5, #2 in the table) is deck-trivial on every deck. But the refined
  pattern holds and is sharper: **every k = 12 flagship acquires its k at
  exactly one deck-nontrivial rung and propagates it through (R)-rungs.**
  Gross and 288 both root at bb72 — the (6,6) quarter of [[288,12,18]] is
  literally bb72 (y⁷ ≡ y mod 6) — and *bb72's own k = 12 is deck-born* from
  k = 8 quarter-codes at (3,6). The 360 code's k = 12 is born at its top rung.
- **Wall-attainment prediction (new A17-docket candidate).** [[288,12,18]]
  is an all-(R) iterated bb72-cover; PAR forces its bases' distances even;
  18 = 2·10 − 2. If its (6,12) base `A = x³+y²+y⁷, B = y³+x+x²` (n = 144,
  k = 12) has d = 10, then [[288,12,18]] sits *exactly at the deficit wall
  2d − 2* — which no measured instance has attained (A17 §corrected picture).
  Ladder is cheap at n = 144; queued behind the in-flight SAT rounds.
  **Verdict (2026-07-20, same day): REFUTED.** The ladder ran — floor@9 and
  floor@11 over all 155 orbit reps, both all-UNSAT (191 s / 21.7 min,
  `data/a19/e288_base_floor{9,11}.log`): **d(base) ≥ 12 certified**, so the
  wall value `2d − 2` is out of reach for this cell. **Closed same day:
  floor@13 returned nine weight-12 witnesses (4.0 h) — d(base) = 12
  EXACT, the cell sits at `18 = 2d − 6`, deficit-6 certified 2-for-2**
  (with bb108-y; A17 §8b). The base is [[144,12,12]] — gross parameters,
  the perfect y-double of its bb72 half — so [[288,12,18]] is bb72
  doubled twice: perfect rung, then deficit-6 rung (mirror of A20's
  wall-then-perfect tower). A catalog-wide deck-birth survey
  (`scripts/a19_deck_survey_catalog.py`, IBM supplemental Tables I–III,
  218 k-validated rows) simultaneously confirmed the refined design
  principle for every k = 12 class (literal) and surfaced two follow-up
  targets: **[[288,12,16]]** (IBM class a, A = 1+y+x, B = y³+x⁵+x¹⁰ at
  (12,12)) — a second deck-nontrivial-top flagship, d = 16 MILP-exact —
  and **[[288,8,20]]** (IBM class Y, (18,8), d = 20 MILP-exact), born at
  (18,2) and tower-structured with bases at n = 144/72: the cheapest
  known d ≥ 20 certification target.
- **Lean landing zone verified**: in `BBDoubling.lean`, `dangerous_zero_rung`
  (l. 441) and the shape rungs consume `StrongBaseFloor`/`LogicalFloor` only;
  `DeckTrivialOnH1` first appears in the assembly theorem (l. 1280). The
  dangerous-side machinery is already (R)-free as formalized — the M12/M24
  program needs no refactor to land. m(b): `gross-distance-proof.md` l. 1059
  (off-support minimum), hexagon rung Lemma 11 l. 1477 — the census port
  starts there.
