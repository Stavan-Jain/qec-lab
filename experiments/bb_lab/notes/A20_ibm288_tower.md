# A20 — IBM [[288,8,20]] as an all-(R) tower: the cheapest d ≥ 20 certification target

**Session 2026-07-20.** Origin: the A19 catalog-wide deck-birth survey
(`scripts/a19_deck_survey_catalog.py`, research_log entry
`a19-wall-refutation-and-catalog-survey`) surfaced IBM class Y —
`(ℓ,m) = (18,8)`, `A = 1+xy⁴+x¹⁴y`, `B = 1+xy²+x²y⁷`, `[[288,8,20]]`,
d = 20 **MILP-exact** per arXiv:2606.02418 Table II — as deck-born at
(18,2) with two k-preserving y-rungs above it. A d = 20 code whose whole
tower is SAT-ladder-accessible (bases at n = 144 and n = 72): the cheapest
known concrete target for the program's d ≥ 20 question.

Scripts: `scripts/a20_ibm288_tower.py` (instance study; kernel-level F₂).
Data: `data/a20/` (tower report, Bezout witnesses, ladder logs).

## 1. Tower structure (verified kernel-level, 2026-07-20)

```
Y2 (18,2) [[72,8,d₀]]  --y-->  Y4 (18,4) [[144,8,d₁]]  --y-->  Y8 (18,8) [[288,8,20]]
A: 1+x+x¹⁴y                    1+x+x¹⁴y                        1+xy⁴+x¹⁴y
B: 1+x+x²y                     1+xy²+x²y³                      1+xy²+x²y⁷
```

| rung | frame | (n,k) | deck k's | (R) on y-deck |
|---|---|---|---|---|
| Y8 | Z₁₈×Z₈ | [[288,8]] | x→4, y→8, xy→4 | **HOLDS** — Bezout witness `data/a20/bezout_y_18x8.json` |
| Y4 | Z₁₈×Z₄ | [[144,8]] | x→4, y→8, xy→4 | **HOLDS** — witness `data/a20/bezout_y_18x4.json` |
| Y2 | Z₁₈×Z₂ | [[72,8]] | x→4, y→4, xy→4 | FAILS (k jumps on every deck — the birth rung; A12-consistent) |

**Headline: this is the gross shape, not the Bravyi-360 shape.** k = 8 is
deck-born at (18,2) and propagates through two (R)-rungs; the classical
`BBDoubling` template applies at BOTH rungs, with `DeckTrivialOnH1`
dischargeable from the saved Bezout witnesses via `deckTrivial_of_bezout`
(A12/A13 line). Frames: (18,2) is floor-bearing (4∤18, 4∤2) but the pair is
NOT in the A16 mirrored class (π_y(A) = 1+x+x¹⁴ is not monomial), so base
floors come from SAT ladders, not the class certificate. PAR applies
(trinomials): every distance in the tower is even.

## 2. Distance arithmetic — d₀ = 6 measured; d₁ decides the headline

**d₀ = 6 EXACT** (2026-07-20, 0.5 s: weight-6 witness + all 25 orbit reps
UNSAT@5, `data/a20/y72_ladder.log`). Y2 = [[72,8,6]]. With d₂ = 20 exact
(IBM MILP) and PAR forcing even values, the live scenarios are:

- **d₁ = 10: rung 1 attains the wall** (10 = 2·6 − 2 — the first
  wall-attaining instance ever, immediately after the [[288,12,18]]
  candidate died, A17 §8b) **and rung 2 doubles perfectly** (20 = 2·10).
- **d₁ = 12: rung 1 doubles perfectly** (12 = 2·6) **and rung 2 is a
  deficit-4 cell** (20 = 2·12 − 4, the BX/BY value).
- d₁ ∈ {6, 8}: rung-1 freeze/deficit and rung 2 *exceeds* 2d₁ — novel
  (τ-killed light classes); would contradict no theorem but no precedent.

The n = 144 ladder (in flight, `data/a20/y144_ladder.log`) decides.

## 3. Staged certification plan

1. **Exact d₀** — `a15_coset_distance` full ladder at n = 72 (in flight,
   `data/a20/y72_ladder.log`).
2. **Exact d₁** — same at n = 144 (queued behind the 288-base floor@13
   round; ~30–90 min/round at 8 workers).
3. **Certified d(Y8) ≥ 12-grade floor, cheap**: Theorem-B projection
   transfer needs only min(d₁, µ(Y4)) — ladder + stabilizer floor at
   n = 144. No new theory.
4. **Full 20 via the template at rung 2** (target m = 2d₁ if (d₀,d₁)
   lands 20 = 2d₁): `LogicalFloor d₁` (SAT, n = 144), `DeckTrivialOnH1`
   (witness saved), `DangerousFloorNZ 20` — slice/(M)@20 over Y4:
   light-stabilizer classification to weight 19 at n = 144, SAT-assisted
   census + m-rungs (Prop-10 pattern, deeper: ~3 generator-hexagon sums);
   `SeamCosetFloor 20` — the hard item: UNSAT@19 at n = 144-side is at or
   past the SAT wall, so this is the **first live consumer of the A13
   transport-era safe-sector machinery on an (R)-tower** (odd part Z₉:
   CRT components F₂/F₄/F₆₄ — the slot-frame port question is F₆₄, not
   F₁₆; "floors port, classifications don't").
5. **Lean staging**: BaseFloors-style packaging of d₀/d₁ (named
   hypotheses, Z5Z15F2A6 model), `deckTrivial_of_bezout` instantiation
   from the saved witnesses, `bocksteinVanishes_of_orderFourLift` applies
   (order-4 lifts exist: Z₁₈×Z₁₆ over Y8's y-deck, Z₁₈×Z₈ over Y4's).

## 4. Stakes

d(Y8) = 20 certified — even at mixed grade — would be the first certified
d ≥ 20 for any BB code, on a code IBM's pipeline could only reach by MILP
incumbent. The tower is ~half the scale of Bravyi-360 at every stage
(n = 144 vs n = 180 bases; sector floors at 20 vs 24), making it the
natural pathfinder for the A19 full-24 program: every piece of machinery
built here (census depth, F₆₄/F₄ engine question, template-with-Bezout
Lean wiring) ports upward.

## 5. State

- Instance study: DONE (§1 table; `data/a20/tower_report.json`; Bezout
  witnesses for both y-rungs saved).
- **d₀ = 6 exact** (`data/a20/y72_ladder.log`). d₁ ladder: in flight.
- (M)@20 census, SeamCosetFloor 20, Lean staging: not started.
- Certified floor today: d(Y8) ≥ 6 by Theorem-B projection transfer once
  µ(Y4) ≥ 6 is checked (Y4 stabilizer floor — cheap; queued with d₁).
