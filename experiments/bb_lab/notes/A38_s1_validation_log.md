# A38 S1 — all-validation session: the F5 engine skeleton + banked reproductions + R4 closure + the burden map + the F2a probe

**Session 2026-08-18** (worktree branch
`claude/distance-descent-theory-plan-e1fd17` continuation). Charter §6
"S1" executed in full: (1) the A35 generality screen promoted into the
library as **`bb_lab.tower`** with a rank-generic group core, gated on
the banked asserts; (2) the three banked closures (A32 bravyi360
d = 24, A33 ibm288Y d = 20, A36 bb288 d = 18) reproduced headlessly,
every count checked against the notes; (3) the **R4-regime methods
closure** (gross_xx lower pair) — A35's pair-regime coverage is now
**4/4**; (4) the **corpus envelope census** — the first empirical
burden map (73 rows priced, no censuses run); (5) the **F2a
odd-Fourier factorization probe** on the banked gross stab census —
verdict: measured absence, with the exact structural obstruction.
**No new distance claims anywhere in this session** (per charter);
every verdict below is a cost verdict or a re-derivation at
certificate tier of an already-known value, labeled as such.

Scripts `a38_*.py`; data `data/a38/`; the copied A35/A36 session
scripts carry provenance headers (source of truth remains the unmerged
branch `claude/tower-slice-calculus-generalize-410ed1` until it
merges — registry renumber-at-merge still owed there, charter §7).

## §1 The F5 engine skeleton: `bb_lab.tower` (+ its gate)

`src/bb_lab/tower.py` — the tower slice calculus as a library, following
the `a30_coset_bz.py → bb_lab.cosetbz` promotion precedent:

- **Rank-generic core**: `TowerCode(name, orders, A, B)` over
  `AbelianGroup` order tuples of ANY rank (never scalar (ℓ,m));
  `AxisDeck(cover, base, axis)` = one-axis index-2 fold with the full
  Lemma-1 constructor asserts (chain map, stabilizer transport,
  twist-invariance im S ⊆ ker E, sections); literal polynomial descent
  (`fold_support`), H1 maps, translation action/orbits/canonical keys,
  the exact two-window node formula (G1), the cap ≤ 4 restricted-MITM
  lift-fiber enumerator, the A35 per-rung/per-pair screen with the cost
  gates, `tower_inventory` (v₂ per axis + odd part + fold chain), and
  `validate_banked()` — the falsify-first regression gate the charter
  §6.0 requires of every engine change.
- **The gate (`scripts/a38_s1_screen_gate.py`, 12.9 s), all PASS**:
  - G1: node anchors ×1.00 (a32 W-coset ≤ 22) / ×1.01 (a32 stab ≤ 22) /
    ×3.00 (a33 H5 — the shared-walk accounting datum, asserted AS
    ×3.00 so future cost models can't silently "fix" it); sector-C
    fiber layer bit-level: **397 orbit fibers, 4,132 lifts, m₂-hist
    {0: 416, 1: 554, 2: 561, 3: 1639, 4: 962} EXACT**; banked
    heavy-band empty rates re-read (319/1733, 5635/10602, 55555/64619).
  - G2: the full 11-tower A35 docket re-screened through the library is
    **field-identical to the banked `screen.json` including the
    RNG-sampled fiber sections** (seed 20260811; the rank-generic code
    consumes the stream exactly as the frozen rank-2 screen did) —
    only `wall_s` exempt. Structure asserts: banked A32 (k 12/8/8,
    ranks 6/4/4, exactness pattern, K ⊆ S, dim W = 2, preimage = S),
    banked A33 (k 8/8/8, all-(R), S ∩ K = 0, dim W₂ = 4), A19
    deck-survey k-verdicts (gross 12/12/8; bb288 y-quotient = the
    (12,6) gross frame).
  - G3: the (3,6) bottom demo (all 2²² cycles: d = 4 exact, stab hist
    {6: 18, 8: 45, 10: 108, 12: 639, 14: 1422, 16: 3411}, 27
    min-weight classes / 9 inside im p\*) == banked; the |A|-even
    parity demo == banked; **rank-3 smoke**: a k = 4 trinomial pair on
    Z₄×Z₃×Z₃ with a twisted k-preserving Z₂ x-fold runs the full rung
    screen (rank law 2 = k/2, cover exactness, codim_lift = k_b − k_c/2)
    — a TYPE-level check that the group core is rank-generic; **no
    trivariate claims** (F4b's sweep stays future work).
- Not yet in the library (deliberate S1 scope): the rung engine
  (`a33_rung_cell.YRungCell` consumed as-is from the frozen script
  layer — duck-compatible with `TowerCode`/`AxisDeck`), the deep
  ordered-split fiber lane (cap > 4), and the BZ census wrappers
  (already library-grade in `bb_lab.cosetbz`). Residue §7.

## §2 Headless reproduction of the three banked closures

All three suites run end-to-end by `scripts/a38_repro_driver.py`
(their internal falsify-first hard asserts ARE the check; any assert
failure = nonzero exit = a finding). Tracked banked outputs diffed
against the reruns afterward; banked files restored, rerun copies
archived under `data/a38/repro/`.

### §2.1 A32 bravyi360 [[360,12,24]] (suite: tower_slice, gb_census, subclosures, sectorAC_full, deep_fibers, dby_floor)

- rc = 0 on all six scripts; **every tracked output file identical to
  the banked version modulo `wall_s` only** (verified field-by-field
  after stripping wall_s: gb_census_summary, subclosure_summary + A14
  jsonl, sectorAC_summary, deep_fibers_summary + jsonl (28 rows),
  dby_floor_summary, tower_validation).
- Key counts reproduced: 8,461 banked BY stabilizers decompose
  (transport + slice asserts); band-16 census = 3 A + 2 B + 1 C
  (A24 §2.5 split); sector A 18–22: 196,557 fibers → 13,109 lifts →
  **13,109 rungs ALL PASS**; sector C 14–16: 397 fibers → 4,132 lifts
  → 2,770 rungs PASS + 765 flat-22; sector C 18–22: 76,954 fibers →
  38,223 lifts → 31,613 rungs PASS; deep fibers 28 → 2,371 lifts →
  2,030 rungs PASS; **flat-22 total 25,492 + 765 + 895 = 27,152 ==
  the note's 27,152 EXACTLY**.
- Wall: 1,281 s vs the banked ~9.4 min (564 s) — contended with the
  S1 gate/R4/envelope runs on the same box; wall time is not a banked
  invariant, counts are.

### §2.2 A33 ibm288Y [[288,8,20]] (suite: tower_cells, validate_banked, h5_close, h5_descent, solver_free)

- rc = 0 on all five scripts. Counts reproduced (all hard-asserted
  in-script): census 1,655/1,655 rows load + v7 == SAT class sets;
  **H2 1,655/1,655 rungs PASS with banked SAT agreement 1,655/1,655**;
  H5 direct: 6.62e10 nodes exact, 118,932 stab vectors ≤ 18 → 1,655
  classes, 12-orbit seam coset **1,680 elements** {14: 6, 16: 84,
  18: 1,590} (280 classes), 3-orbit **EMPTY**, seam rungs
  **1,680/1,680 PASS**; descent: 0x36 census 105,328 / 0xf3 81,108,
  90.6% carry-infeasible, **1,680 SEAM lifts == direct EXACTLY**,
  3-orbit 0; solver-free: d₀ = 6 (84 w6), H3 = 10 (72 w10), Y2 stab
  census 165,517 → 4,605 orbit reps → class-key union == banked
  1,655.
- Output files identical to banked modulo timing fields — the
  first-pass checker flagged `banked_validation.json` + `h2_rungs.jsonl`
  as mismatches, which resolved to per-lane/per-rung timing keys named
  `tot_s`/`secs` rather than `wall_s` (finding: banked artifacts carry
  timing under three key names; the a38 repro checker strips all
  three). No count field differs anywhere. Wall 288 s vs banked
  ~105 s (contended box).

### §2.3 A36 bb288 [[288,12,18]] (suite: tower_cells, direct_close, witness, descent)

- rc = 0 on all four. Counts vs the A36 note, all exact:
  **7.50e9 census nodes** (binomial identity), stab census **33,588
  vectors** {6: 72, 10: 432, 12: 2,268, 14: 3,888, 16: 26,928} →
  **469 orbits** {6: 1, 10: 6, 12: 33, 14: 54, 16: 375}, no weight-8
  stabilizers; seam census **395 elements** {12: 3, 14: 18, 16: 374}
  over the 5 orbit-rep cosets (orbit split 3+3+9+12+36, stab orders
  24/24/8/6/2); V1 4,096/4,096; **469/469 dangerous + 395/395 seam
  rung PASSes**; witness re-verified end-to-end from support AND the
  exhaustive ladder re-run finds it again at the same stratum (weight
  18 over the weight-12 seam element, class 0x40, m₁ = 3, 4 finds all
  ≥ 18); descent lane: bb72 stab census 9,495/272 orbits, cycle census
  1,110 {6: 120, 8: 990} ⟹ d(B72) = 6 census-complete, base-exactness
  1,110/1,110, key-set equality stab 469 == direct / seam 114 ==
  direct.
- **The exhaustive witness ladder re-run regenerates the banked
  witness file byte-identically** (modulo timing): same stratum, same
  support, 4 finds all ≥ 18. Every output file identical to banked
  modulo timing fields. Wall 39.8 s vs banked ~31.8 s.

## §3 The R4-regime methods closure (`scripts/a38_r4_close.py`, 5.7 s)

**Pair-regime coverage now 4/4** (R1 A32, R2 A33, R3 A36, R4 here).
Tower GR (12,6) → B72 (6,6) → B36 (3,6); banked R4 lattice reproduced
(dim S = dim K = 6, dim S∩K = 2, dim W = 4, K ⊄ S). Values NOT new —
d(gross) = 12 is kernel-checked (stronger tier), d(B72) = 6 and
d(B36) = 4 banked; the regime execution is the deliverable.

- **Direct closure** (n = 72 mid-level): d(B72) = 6 census-complete
  re-derived (banked numbers exact); stab census ≤ 10 = 252 vectors
  {6: 36, 10: 216} (== the banked A36 ≤ 16 census truncated) → 7
  orbits → **7/7 dangerous rungs PASS** (V1 4,096/4,096 at k = 12);
  **the 5 seam orbit-rep cosets are census-EMPTY at ≤ 10** ⟹ the seam
  branch closes by emptiness, 0 rungs; τ-diagonal weight-12 witness
  verified end-to-end (all 84 weight-6 B72 logicals lie OUTSIDE
  im p_x\* — measured, cited per-instance only per the A35 ledger).
  Assembly: d ≥ 12 ∧ d ≤ 12 ⟹ **d = 12 re-derived at certificate
  tier**, consistent with the Lean theorem.
- **The regime finding**: this target sits EXACTLY on the G5 ceiling
  (12 = 2·d(B72)), and there the R4 trisection never fires on the
  closure path — seam cosets are empty at W = 2d(mid) − 2 and the
  dangerous sector is 7 orbits. R4's "no automatic pruning" bites in
  the DESCENT direction instead:
- **R4 descent lane** (both n = 72 species from n = 36 data, full
  trisection): the seam species' trisection is CLASS-DETERMINED — 3
  SEAM classes (S∩K ∖ 0) have stab-or-0 shadows, 60 have sector-A
  shadows in W ∖ 0. Sector C: 10 stab-orbit fibers → 11 lifts (7 stab
  + **4 other-class leakage**), 50% empty. Sector B: τ_b family over
  the 54 weight-4 B36 cycles → 54 lifts, ALL other-class; the non-(R)
  signature measured: **ker τ_b\* vs im p_b\* class agreement only
  36/54** (base-exactness fails on the bottom rung, as the lattice
  predicts — the A36 assert "[γ] ∈ im p\* ⟺ τ(γ) stab" provably does
  not port to non-(R) rungs; classify by ker τ\* directly). Sector A
  (the R4 novelty): 102 W-coset element orbit fibers → 84 lifts, **ALL
  84 outside SEAM** — the preimage identity (p_b\*)⁻¹(W) = S of the
  R1 regime fails totally here, quantified. **Key-set equality with
  the direct censuses holds on both species** (stab 7 == 7, seam
  0 == 0): completeness survives the leakage, you just enumerate and
  discard more.

## §4 The corpus envelope census (`scripts/a38_envelope_census.py`, 3.2 s)

The first empirical burden map: `data/a38/envelope_census.{json,md}`.
73 rows priced (pricing ONLY — no censuses run; GREEN/AMBER/RED are
cost verdicts, never distance claims): the 47-group corpus frontier
(best exact-d rep per group + the open d_ub rep where it out-ranks
exact), the 5 instantiable zoo BB entries, and the 11 A35-docket rows
as a calibration battery — **calibration gate 11/11: verdict, node
counts, and caps all equal the banked screen exactly**.

**Verdicts: GREEN 41 / AMBER 3 / RED 6 / no-deck 23.** Headlines:

1. **Every exact-d corpus representative with a deck prices GREEN
   (27/27)** — re-certifying every known best-of-group distance is
   inside the demonstrated envelope wherever 2 | |G|.
2. **The open frontier is CAP-bound, not node-bound**: all 8
   deck-having open d_ub questions (d_ub 26–42, n = 108–168) sit
   WITHIN the 2e11 bottom-node envelope but at caps 9–17 (beyond the
   demonstrated 8) — AMBER at d_ub 26–30, RED at d_ub ≥ 36. **The
   corpus frontier wall is W3 (fiber-cap growth), the F2b/F2c target —
   not W2 (census blowup); [[756]]-style census walls are atypical
   in-corpus.** This sharpens the charter's S3 priority.
3. **W1 (odd |G|) locks out 22 of 57 corpus rows** including the
   highest-d_ub open questions (Z9xZ9 ub 40, Z7xZ9 ub 30) — measured
   demand for F1.
4. The corpus has grown since A18: 47 group shapes / 58,350 rows (41 /
   58,021 at A18 close). Groups whose stale d_ub exceeded a
   later-solved d_exact are NOT frontier questions (checked:
   Z5xZ9/Z3xZ35/Z6xZ15 have zero open rows).

## §5 The F2a probe (`scripts/a38_f2a_probe.py`, 19.0 s): the census walk does NOT factor along the odd part — measured absence, with the obstruction

Target: the cheapest fully-banked census — the gross stabilizer census
at W = 16 (A36 offset-S species). **P1 falsify-first: the census is
reproduced from scratch before any measurement** (33,588 vectors,
banked weight hist, 469 orbits, 7.502e9 nodes exact).

Setup: G = Z₁₂×Z₆, odd part H = Z₃×Z₃; N = Σ_H h central idempotent;
the 4-way odd-CRT idempotent system {N₁N₂, N₁E₂, E₁N₂, E₁E₂}
(partition of unity + orthogonality + idempotency all asserted; the
E₁E₂ block further splits into two F₄ components — recorded, not
consumed).

- **Maschke rank bookkeeping EXACT**: κ = 66 splits 8 + 16 + 16 + 26.
  Side-finding: the invariant sector is the quotient code over
  G/H = Z₄×Z₂ and it is **[[16,0]]** — the gross code's odd-quotient
  has k = 0; the invariant sector carries NO homology, only
  stabilizer.
- **The invariant weight law** (a small theorem, asserted on all
  33,588 vectors): N-components are unions of H-cosets (9 points
  each), and they are stabilizers hence even ⟹ |v_N| ∈ 18ℤ. Measured
  distribution {0: 288, 36: 720, 54: 12,744, 72: 6,120, 90: 13,104,
  108: 612} — **any census at W < 18 is structurally invisible to the
  invariant sector except through cross-sector cancellation.**
- **The measured obstruction is total**: weight additivity
  |v| = Σ_c |v_c| holds on **0/33,588**; all components ≤ W on
  **0/33,588**; pure-sector vectors 0 in every sector; v_N = 0 on only
  288/33,588; per-component maxima (108, 84, 84, 60) ≫ 16. Light
  stabilizers are globally light but sector-wise HEAVY — cancellation
  across the odd-CRT sectors is the rule, not the exception.
- **Verdict: the F2a factorization hypothesis is ABSENT at the census
  bound** — the weight-bounded walk cannot be split into per-sector
  walks; the 2-part stays combinatorial and the odd part enters only
  through cancellation. (The naive "×63.66 gain if it had factored"
  at the measured per-sector maxima is a fiction: it prices the
  per-sector walks but not the recombination product, and the
  factorization premise itself is false.) Charter F2a asked for
  falsify-first with counts reproducing exactly; this is the
  falsification, banked with the joint component-weight distribution
  (`data/a38/f2a_probe.json`) for F1-Q1 to consume: **any odd coupling
  law must be a cancellation law, not a weight-splitting law.**

## §6 Verification map

| claim | check |
|---|---|
| library == frozen screen | gate G2: 11/11 towers field-identical to banked screen.json incl. RNG-sampled fibers (wall_s exempt) |
| node anchors | ×1.00 / ×1.01 / ×3.00 asserted as exact values (shared-walk datum pinned) |
| fiber enumerator bit-level | 397 fibers / 4,132 lifts / m₂-hist EXACT vs banked a32 sector-C layer |
| rank-generic core | rank-3 k = 4 twisted-(R) rung: full screen_rung asserts at rank 3 (type-level) |
| A32 reproduction | 6/6 scripts rc = 0 (internal hard asserts); all tracked outputs identical mod wall_s; flat-22 total 27,152 exact |
| A33 reproduction | 5/5 scripts rc = 0; 1,655 ×(census, rungs, SAT agreement), 1,680 ×(seam census, rungs), descent equality; outputs identical mod wall_s |
| A36 reproduction | 4/4 scripts rc = 0; 7.50e9 nodes, 33,588/469, 395 seam, 469+395 PASSes, witness at w18 over the w12 seam element re-found by the re-run ladder; outputs identical mod wall_s |
| R4 lattice | banked A35 pair numbers hard-asserted (6/6/2/4, K ⊄ S) |
| R4 closure soundness | exact node counts; census cross-asserts vs banked A36 numbers; V1 sector scan 4,096/4,096; covariance 3+3; fold surjectivity; witness re-verified from support |
| R4 descent completeness | key-set equality direct == descent on both species (leakage discarded, counted) |
| envelope calibration | 11/11 docket rows == banked verdict/nodes/cap |
| envelope inputs | corpus read-only (MAIN checkout DuckDB); zoo rows from zoo.yaml presentations |
| F2a probe | census reproduced exactly before any measurement (see §5) |

## §7 Falsified-claims ledger (session-internal)

- **"Base-exactness classifies τ-lifts" as a general rung fact —
  REFUTED in-session on the non-(R) rung** (B72 → B36): ker τ_b\* vs
  im p_b\* agree on only 36/54 of the light-cycle classes. The A36
  1,110/1,110 assert is an (R)-rung fact; classify by ker τ\* directly
  everywhere (the library note carries this).
- **"The R1 preimage identity approximately survives in R4" — NO**:
  100% of sector-A lifts (84/84) fell outside SEAM at this instance;
  reachability is NOT decided one rung below in R4, it is decided by
  enumeration + discard.
- **The F2a factorization hypothesis at the census bound — measured
  absent** (§5): the invariant-sector components of censused vectors
  are forced into 18ℤ by the H-coset-union weight law, so a W < 18
  census cannot see the invariant sector except through cross-sector
  cancellation; weight additivity fails broadly. (Charter F2a asked
  falsify-first; this is the falsification, with the obstruction now a
  theorem-shaped statement rather than a vibe.)
- (Respected, inherited): no SAT anywhere in S1; witness weights never
  reported as floors; light-concentration and Prop-10 cited
  per-instance only; RED/AMBER/GREEN are cost verdicts.

## §8 Residue / next steps

1. **Rung engine promotion** (F5 continuation): `YRungCell` (+ the
   deep ordered-split fiber lane, cap > 4) into `bb_lab.tower` or a
   sibling module, gated the same way (the a33/a36 rung counts are the
   regression battery). S1 consumed it from the frozen script layer.
2. **The burden map's sharpest edge → S3**: the corpus frontier is
   cap-bound (W3), so the F2b ε-recursion / F2c carry-floor session
   should target caps 9–17 at n = 108–168 (concrete instances now
   listed in `envelope_census.md`), with [[720,4]] doubling (cap 16)
   as the charter already planned.
3. **F2a residue → F1/S2**: the probe's structural law (invariant
   components ∈ 18ℤ; the census lives in the twisted sector mod
   cancellation) says the odd-Fourier lane must work per-sector with a
   coupling law, exactly F1-Q1; the measured joint component-weight
   distribution is banked for S2 to consume.
4. **Trivariate**: the rank-generic core is exercised at rank 3
   (type-level only); an actual trivariate screen/sweep is F4b.
5. **A37 coordination**: a8_336 stays deferred to A37 per charter §4
   (its docket row untouched here).
6. Wall-time note: reproduction walls are 1.2–2.3× the banked walls
   when suites share the box; the banked invariants are node/element
   counts, never seconds.

## §9 Session commits (this worktree branch, in order)

- `1d2b296` a38: import A35/A36 session scripts + banked reference data
- `639d0dd` a38: bb_lab.tower — the tower slice calculus as a library (F5)
- `daf57ff` a38: R4-regime methods closure — d(gross) = 12 through the
  lower pair
- `496c1a5` a38: corpus envelope census — the first empirical burden
  map (F5)
- `3a2f6ef` a38: F2a probe (absent-with-obstruction) + headless
  reproduction bank
- (this commit) a38: S1 validation log + charter §4/§6 updates
