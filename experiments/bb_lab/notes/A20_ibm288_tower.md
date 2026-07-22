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
UNSAT@5, `data/a20/y72_ladder.log`). Y2 = [[72,8,6]].

**d₁ = 10 EXACT** (same day, 180 s: weight-10 witness + all 25 orbit reps
UNSAT@9, `data/a20/y144_ladder.log`). Y4 = [[144,8,10]]. Therefore:

- **Rung 1 ATTAINS THE DEFICIT WALL: 10 = 2·6 − 2.** The first measured
  instance at exactly `2d − 2` (A17 §8b entry; supersedes "no measured
  instance attains the wall"), found the same day the [[288,12,18]]
  candidate was refuted. (R) holds on this rung, so the T2 machinery
  applies to the attaining cell directly — the natural specimen for the
  A17 §9 maxSF residue question.
- **Rung 2 doubles perfectly: 20 = 2·10.** τ-lift of the d₁ witness
  verified as a weight-20 nontrivial X-logical of Y8
  (`scripts/a20_tau_lift.py`, `data/a20/y8_weight20_witness.npy`):
  **d(Y8) ≤ 20 constructive**, independent of IBM's MILP, and the
  τ-tightness obligation of the rung-2 template is discharged.

One tower, both extremes: the maximal-deficit cell and a perfect doubling,
stacked. Whatever mechanism sets the deficit acts per-rung, not per-code.

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

## 5. State (end of opening session, 2026-07-20)

Rung-2 template ledger for certified d(Y8) = 20 (target m = 2d₁ = 20):

| obligation | status |
|---|---|
| `LogicalFloor 10` at Y4 | **DONE (solver-grade)** — the d₁ ladder's UNSAT@9 over all 25 orbit reps |
| `DeckTrivialOnH1` | **DONE** — Bezout witness `data/a20/bezout_y_18x8.json`, dischargeable via `deckTrivial_of_bezout` |
| τ-tightness (upper bound) | **DONE** — verified weight-20 witness `data/a20/y8_weight20_witness.npy` |
| `DangerousFloorNZ 20` | **floors DONE, exhaustiveness pending** — census (`a20_m_census.py`): bands 2–16 closed (469 classes: 1×w6, 6×w10, 33×w12, 54×w14, 375×w16; **no weight-8 stabilizers** — gross's Prop-10 gap echoed; µ(Y4) = 6); band 18 at 1,186+ classes, still enumerating. Per-class fiber-pinned UNSAT floors (`a20_m_floors.py`): **all 1,655 censused classes certified in ~17 s total, 0 SAT hits** — incl. the hexagon m ≥ 7 in 9 s (the flagged analytic-tail risk, retired). Remaining: band-18 UNSAT close-out + floors re-run on stragglers |
| `SeamCosetFloor 20` | **FALSE as stated** (2026-07-22): the 15 seam classes collapse to 2 G-orbits (12+3, both δ₂-injective, `a20_seam_floor.py`); class 0x1's coset contains a **weight-18 element** — the seam minimum sits at the wall value 18 = 2d₁ − 2, in the same tower whose rung 1 attained the wall in d. NOT a refutation of d = 20 (the floor is sufficient-not-necessary). Replaced by the **lift-aware seam floor** (`a20_safe_m_floors.py`, running): census every coset element with \|w\| ≤ 19 (all even), and per element the fiber-pinned lift query — every cover cycle over w is *automatically* a nontrivial logical (δ₂-injectivity), so UNSAT at off-support < (20−\|w\|)/2 certifies all its logicals ≥ 20. Class 0x3 probe (\|seamC\| = 84) still running |

Also banked: d₀ = 6, d₁ = 10 exact; the rung-1 wall-attaining cell (A17
§8b); `bocksteinVanishes_of_orderFourLift` applies at both rungs (ZMod
lifts (18,16)/(18,8)). Lean staging (named-hypothesis packaging of d₀/d₁
on the Z5Z15F2A6 model + Bezout instantiation) not started.

## 6. V7 completeness CLOSED — the analytic census; H1 discharged (2026-07-22)

**Headline.** The completeness gap of the boundary census is closed
analytically, at A22-V7 grade (every ingredient a finite algebraic fact —
kernel-checkable in principle, no solver trust): the complete list of
G-translation classes of nonzero b ∈ im ∂₂(Y4) with |b| ≤ 19 is **exactly
the 1,655 SAT classes**. Closed bands reproduced exactly (469); **band 18
is definitively 1,186** — the unterminated SAT enumeration had already
found every class; zero new. Hence §7's H1 is **DISCHARGED** (supersedes
its PENDING marker), and with H2's 1,655/1,655 UNSAT floors the dangerous
sector (`DangerousFloorNZ 20`) is complete outright: **5 of §7's 6
hypotheses now done; only H5 (lift-aware seam floor) remains** for the
d = 20 theorem. No census-resume needed — the stopped enumeration stays
stopped.

Scripts: `a20_v7_lever0.py` (pre-sweep measurements),
`a20_v7_completeness.py` (the engine; rerun to regenerate
`data/a20/v7_complete_classes.jsonl` + `v7_summary.json` +
`v7_completeness.log`), `a20_v7_new_floors.py` (floors for
census-missed classes — vacuously done, 0 rows).

**Lever autopsy (measured first, per the session brief;
`a20_v7_lever0.py`, `data/a20/v7_lever0.json`).**

- **Lever 0 (the A22-P1 analog) — DEAD.** min |f| over the 16
  joint-kernel translates has a heavy band-18 tail: {3: 471, 4: 592,
  5: 101, 6: 7, then 9…23: 15 classes}; min site-support reaches 8/8
  (10 classes). No small-support-f sweep exists (and no converse bound
  would have been available). Bonus structure banked: ker A⋆ = ker B⋆ =
  ker ∂₂ verified directly at Y4 (dims 4/4/4 — the A23 gating fact),
  and the 15 nonzero kernel elements are δ₄-pure with weight exactly
  6·(active sites): weights {36: 12, 48: 3}, site-supports {6: 12, 8: 3}.
- **Lever 5 (A23's flat min-side strata, coordinator-proposed) —
  INFEASIBLE here.** min(|u|,|v|) reaches 9 on 500/1,655 classes, so
  flat strata must run to b = 9: Σ_b C(72,b)/16 ≈ 6.2·10⁹ elements,
  ~31× A23's f2a6 sweep, with only 4 dense parity checks to prune. The
  fibering below is exactly the refinement that makes the same
  min-side sweep per-site-bounded.
- **Levers 1–3 — ADOPTED**, unlocked by two NEW finite facts that
  dissolve the session-1 obstruction ("the inner GF(64)⁸ enumeration
  lacks a finite reduction"):
  - **(P) Excess parity.** W(ε,d4,ξ) ≡ ε (mod 2) for every ξ (fiber
    weight parity = augmentation), so the per-site excess over
    W_min(ε,d4) is EVEN. The dominant B₀ = 18 stratum (54,064 of
    74,528 survivor cells) is excess-0-only; every deviation budget
    halves.
  - **(A) Argmin = 3^cost.** |{ξ : W(ε,d4,ξ) = W_min}| = 3^{W_min}
    exactly (sizes 1/3/9/27). So a **min-cost information set** I of
    the [16, 8] GF(64) graph code {(w, C₆₄w)} — matroid greedy on the
    rows of [I₈; C₆₄] with per-coordinate costs W_min — bounds each
    cell's enumeration by 3^{cost(I)} × even-deviation terms.
    Completeness: excesses are ≥ 0 sitewise, so any |b| ≤ 18 solution
    restricted to I has excess-sum ≤ E := 18 − B₀. Balanced-cell
    fallback (also complete): two one-sided passes at halved budget,
    since min(exc_u, exc_v) ≤ 2⌊E/4⌋ by parity.

**Ingredient list (I1–I7, all finite algebraic facts).** I1 the CRT
fibering + 512-triple weight table (re-derived in-script); I2 the six
component operators extracted from D2 = H_Xᵀ itself by 8 site-delta
applications (fiber value 1 = CRT triple (1,1,1), a unit in every
component), verified against D2 on random f AND on all 1,655 census
rows; I3 im ∂₂ = F₂⁸ × GF(4)⁶ × GF(64)⁸ exactly (ε free, v_ε = C_ε u_ε;
δ₄ = the 6-dim joint image, 4,096 pairs; δ₆₄ free, v₆₄ = C₆₄ u₆₄);
I4 = (P); I5 the outer bound B₀ ≤ |b| (survivor histogram = the
session-1 16.7M pre-filter ÷ 16 kernel redundancy, exact match);
I6 = (A) + the information-set restriction bound; I7 unique-preimage
reconstruction per triple + canonicalization over the 72 translations.

**Run accounting.** 2⁸ × 4⁶ = 1.05M outer cells → 74,528 survivors;
674,288,722 inner candidates (strategies: mixed info-set 73,162 cells,
halved two-pass 1,366); 118,932 light boundary vectors — equal to
Σ orbit sizes over the 1,655 classes EXACTLY, i.e. every G-translate of
every class was produced exactly once (the strongest internal
completeness witness short of Lean). 1,105 s single-core.

**Lean feasibility (assessed, not started).** Structurally A22-L1–L6:
fibering lemmas, finite tables, component-matrix certificates (8×8 over
F₂/GF(4)/GF(64), all F₂-linearizable), B₀ arithmetic, per-cell
info-set rank certificates, glue. The mass is the difference: 74,528
cells / 6.7·10⁸ weight evaluations vs A22's 16,384 subsets / 6·10⁵
checks. Direct decide-scale needs the ÷24 outer symmetry quotient (z³
acts trivially on outer data) and packed-table batching; judged
feasible-with-engineering, one dedicated session to prototype.

## 7. Assembly: certified d(Y8) = 20 — the theorem of record (DRAFT, two verdicts pending)

**Claim.** IBM's class-Y code Y8 = [[288,8,d]] (arXiv:2606.02418 Table II)
has d = 20, certified at solver grade (CryptoMiniSat UNSAT certificates +
kernel-level F₂ linear algebra + explicit verified witnesses), independent
of IBM's MILP.

**Upper bound.** The τ-lift of the d₁ witness is a verified weight-20
nontrivial X-logical (`a20_tau_lift.py`, `data/a20/y8_weight20_witness.npy`)
⟹ d ≤ 20.

**Lower bound — case assembly over any nontrivial cover X-logical v.**
By T2 (`push1_mem_seamCoset_of_deckTrivial`, Lean, under H4) + the L0 iff
(`pull1_mem_boundaries_iff_seamCoset`, Lean), p(v) lies in a seam coset
seamC ζ + im ∂₂ for some ζ ∈ ker ∂₂(Y4):

- **ζ = 0, p(v) = 0** (new sector): v = τ(u) with u a nontrivial base
  logical (τ injective chain map), so |v| = 2|u| ≥ 2·10 = 20 by H3.
- **ζ = 0, p(v) = b ≠ 0** (dangerous): slice identity |v| = |b| + 2·off.
  If |b| ≥ 20, done. Else b is one of the censused classes (H1) and its
  fiber-pinned floor (H2) gives |v| ≥ 20.
- **ζ ≠ 0** (safe): 15 classes = 2 G-orbits (§5). Let w = p(v), an
  element of the coset. If |w| ≥ 20, |v| ≥ |p(v)| ≥ 20. Else w is in the
  lift-aware census (H5) and — every cover cycle over a non-boundary w
  being automatically nontrivial (δ₂-injectivity, §5) — its lift query
  certifies |v| = |w| + 2·off ≥ 20.

Orbit transfer throughout: all certificates are translation-covariant.

**Hypothesis ledger (provenance; solver-grade unless noted).**

| # | statement | status |
|---|---|---|
| H1 | census of nonzero b ∈ im ∂₂(Y4), |b| ≤ 19, is exhaustive | bands 2–16 CLOSED (469); **band 18 PENDING** (1,186 found; V7 moonshot / census-resume) |
| H2 | per-class dangerous floors at 20 | DONE — 1,655/1,655 UNSAT, ~17 s |
| H3 | LogicalFloor 10 at Y4 | DONE — UNSAT@9, all 25 orbit reps |
| H4 | DeckTrivialOnH1 (R) | DONE — Bezout witness `data/a20/bezout_y_18x8.json` (kernel-level; `deckTrivial_of_bezout`-ready) |
| H5 | lift-aware seam floor at 20, both orbits | class 0x1 **RUNNING** (20/20 lift-UNSAT so far); class 0x3 **PENDING** (probe on budget) |
| H6 | weight-20 witness | DONE — verified |

**Grade and staging.** Solver grade throughout; Lean packaging path:
named hypotheses on the Z5Z15F2A6 model — H1/H2 via the A17 (M)-kernel
emitter pipeline, H4 via `deckTrivial_of_bezout`, H5 via the A23 bridge
recipe (kernel cert + seamC additivity + decide facts; kernel equality
ker A⋆ = ker B⋆ = ker ∂₂ VERIFIED at Y4, so the two-trinomial collapse
applies), H3 as a BaseFloors-style certificate. The transport theorems
consumed here are on QECLean main as of 2026-07-20 (PR #60 merged).
