# descent_theory_test — pre-registration protocol (Phases 0–1)

**Date frozen:** 2026-08-17 (see MANIFEST.sha256 for the exact timestamp
and file hashes).
**Where:** qec-lab worktree `quantum-distance-presentation-73cb27`
(branch `claude/doubling-argument-bb-codes-461d59`),
`experiments/bb_lab/data/descent_theory_test/`.
**What this is:** the PREDICTION side of a falsification-shaped test of
the BB distance-descent/tower theory (A32 tower slice calculus, A33
port, A35 generality map). Phase 0 built a stratified cohort of target
codes; Phase 1 emitted a frozen prediction record per row. Phases 2–4
(closures, ground truth, scoring) happen LATER, by someone else, against
these frozen artifacts. **This is workbench data, not a research line:
no approach number is claimed and nothing here is a distance claim.**

---

## 1. Theory under test and screen provenance

The screen is the A35 generality screen, **extracted** (not
reimplemented) from git history:

    git show 2f063b0:experiments/bb_lab/scripts/a35_generality_screen.py
      -> tools/a35_generality_screen.py

Only the path-bootstrap block and the output directory were edited
(header comment in the file marks the edit). All dependency modules
(`a30_rung_pass.py`, `a32_tower_slice.py`, `a32_subclosures.py`,
`a33_tower_cells.py`, `src/bb_lab/`) and the banked `data/a32/`,
`data/a33/` assets are **byte-identical** between commit `2f063b0` and
this worktree (`git diff --stat` empty), so no assert had to be
disabled.

**Validation run** (mandatory gate before any cohort work; full log in
`validation/screen_run.log`, output in `validation/screen.json`):

- census node anchors: a32_wcoset22 ×1.00, a32_stab22 ×1.01,
  a33_h5 ×3.00 (the documented shared-walk/3-offset accounting — see
  A35 §6);
- fiber enumerator vs banked A32 sector-C layer: 397 orbit fibers,
  4,132 lifts reproduced EXACTLY, m2-histogram equal;
- banked heavy-band empty rates re-read: 18% / 53% / 86% at caps 2/1/0;
- bravyi360: per-level log10 nodes **17.0 / 13.5 / 9.8**, "ALL banked
  A32 structure REPRODUCED" (hard asserts);
- ibm288Y: **13.7 / 10.8 / 7.8**, "ALL banked A33 structure REPRODUCED";
- gross_xx / bb288_yxx: A19 deck-survey k-patterns reproduced;
  bb288 4-level costs 12.5 / 9.9 / 7.1 / 4.4;
- rank law rank p* = rank τ* = k(cover)/2 on all 22 docket rungs; pair
  regimes match A35 §3/§4 row by row;
- (3,6) bottom demo: d = 4 exact by full enumeration, stabilizer weight
  histogram equal to the note's;
- no-deck and parity anti-instances reproduced.

The per-row driver (`tools/dtt_lib.py`) **imports** the extracted
screen's functions (`census_nodes`, `screen_rung`, `fold_terms`, …); it
does not fork them.

## 2. Contamination rules (operated under, and binding for Phase 2+)

1. NO distance closures, exact SAT/MILP solves, certificate closures, or
   floor certifications were run on any cohort TARGET (level-0) code.
2. Exact distances were computed/consumed ONLY for:
   (a) quotient/base codes with n ≤ 40, by trivial full enumeration
       (`full_enumerate_distance`, guarded by an assert on n and on
       level ≥ 1);
   (b) values already recorded in the corpus DB (opened read-only);
   (c) cheap randomized L1 d_ub sampling — upper-bound witnesses only,
       never floors (20,000 samples for quotient levels, 30,000 for the
       fresh parity-strip targets, seed 20260817).
3. The main-checkout corpus DB was never written; all outputs live under
   this directory; the concurrently-running `data/order144_sweep/` was
   not read or written.
4. Compute: everything `nice`d, no solver processes of any kind.

## 3. Cohort definition (as actually applied)

Pool source: `bb_instances.duckdb` (58,350 rows) in the MAIN checkout,
read-only. Total screened in Phase 0: **183 candidate rows** (budget:
≤ 1,000 rows / ≤ 90 min; used ~20 s). Selected cohort: **130 rows**.

### frontier (96 rows)
Filter: `d_exact IS NULL AND d_ub >= 12 AND n BETWEEN 96 AND 168`,
even |G| only. Group quotas (2-adic depth v₂(|G|) in braces):
Z9xZ6 20 {1}, Z6xZ13 12 {1}; Z6xZ10 8, Z5xZ12 8, Z6xZ14 10,
Z7xZ12 10 {2}; Z12xZ6 14, Z8xZ9 14 {3}.
Candidate order: round-robin across d_ub bands (ascending), then
instance_id; overdraw ×1.6; selection maximizes new pair-regime
coverage, then prefers k = 6 rows, then unseen d_ub values, tie-broken
by instance_id (deterministic, no RNG).

### anchor (14 rows)
Corpus-exact rows, one per (group, d) cell, min(instance_id), except
Z12xZ6/d=12 which preferentially selects the gross code (found:
`6f9e8285eedf94e6`). Cells span d ∈ {6,8,10,12} × depth ∈ {1,2,3,4}
(cell list in `tools/build_cohort.py`; all 14 cells were non-empty).
The screen does NOT consume the known d for its window; d_exact is used
only to price the calibration question W = d_exact − 2.

### scope-control (10 rows)
Odd-|G| rows: Z7xZ7 ×2, Z7xZ11 ×2, Z9xZ9 ×2, Z7xZ9 ×2, Z15xZ3 ×1,
Z5xZ9 ×1 (mix of exact and unfilled rows; min(instance_id) per cell).
Frozen prediction for every row: the deck machinery must REFUSE
(no free Z₂ deck exists; A35 conditions map item 1).

### parity-strip (10 rows)
Freshly sampled weight-3 × weight-4 codes (A odd, B even), k ≥ 2, on
Z6xZ10 (5 rows, 89 tries) and Z12xZ6 (5 rows, 14 tries), RNG seed
20260817, canonicalized with the corpus's Aut×translation×swap
machinery; ids = corpus-style canonical hashes, code_id prefix
`dtt_parity_`. d_ub via fresh L1 sampling (rule 2c). These rows are NOT
in the corpus DB. Frozen prediction: machinery sound (decks, transports,
rank law all apply) but the parity layer is lost — odd cycle weights
possible (measured: present at EVERY level of all 10 rows), W_eff = d−1,
no β = 0 parity kills, odd d values allowed (A35 L3/L3-scope).

### regime coverage
All four A35 §3 pair regimes were occupied without top-up
(rows containing each: R1 7, R2 57, R3 41, R4 3; pair-level counts per
depth in `aggregates.json`). Depth-1 towers have no adjacent-rung pair,
so they cannot contribute regime cells (structural, not a gap). R4 at
depth 2 is unoccupied (both screened depth-2 pair regimes landed
R2/R3); this cell is noted as unfilled-from-pool at freeze time.

### exclusions actually applied
- A18's "40 near-miss SAT timeouts at n = 96/112": **no longer
  prediction-eligible** — every n = 96 and n = 112 corpus row now has
  d_exact (closed post-A18 by `maxsat-tandem@mse23+step2` /
  `sat-cadical`; max d = 12 there). The 40 rows with
  `d_method = maxsat-maxcdcl@mse23` are all n = 168 (incl. four d = 14
  rows). The frontier stratum therefore effectively draws from
  n ∈ [108, 168].
- Odd-|G| rows in the frontier window (Z7xZ9, Z7xZ11, Z9xZ9 shapes,
  v₂ = 0) are excluded from `frontier` and represented in
  `scope-control` instead.
- The legacy distance-less Z9xZ6 block (8,319 eligible rows) is
  included at quota 20 (A18 flagged it "low breadth value"; it is the
  only n = 108 shape).
- Depth-4 frontier rows do not exist in the corpus window (all
  Z4xZ12 n = 96 rows are closed); depth 4 is covered by 2 anchors
  (Z4xZ12, Z12xZ12) and arrives in force via the order-144 fold-in
  (FOLD_IN.md).

## 4. Frozen prediction semantics (per record in predictions.jsonl)

- **Route rule (frozen):** maximal descent ladder, all y-axis Z₂ folds
  first (largest→smallest), then all x-axis folds; truncated when a
  quotient's k hits 0 or construction degenerates (recorded in
  `screen.truncated`). Available decks per level are recorded; alternate
  orders are legitimate Phase-2 routes but predictions price THIS route.
- **Structural laws recorded per rung** (measured at prediction time;
  any Phase-2 violation falsifies the theory layer named):
  rank p* = rank τ* = k(cover)/2 (rank law; 0 violations / 245 rungs),
  exactness at the cover (theorem; 0 violations), codim_lift =
  k(base) − k(cover)/2 (0 violations), (R)-trio co-occurrence
  (exact_base ∧ σ*=id ⟺ k preserved; 0 mismatches), Lemma 1 constructor
  asserts (0 failures), Lemma 2 parity (0 violations among
  odd-polynomial levels).
- **G5 window:** [2, min(d_ub, min_j 2^j·d(L_j))] with per-term
  provenance (corpus-exact 179 / full-enumeration 53 / L1-d_ub 13
  across the cohort); windows with any bounded-only term are flagged
  (`ceiling_is_upper_estimate`) and even-snapped only under parity.
  Chain stalls (d(L_j) > 2·d(L_{j+1}), both exact) are pre-declared
  route-stall predictions.
- **Costs:** two-window census formula Nodes(κ, W), κ = n/2 − k/2,
  r₁+r₂+1 = W (validated ×1.00/×1.01); fiber cap (W − μ)/2 with μ =
  min lightest nonzero stabilizer across levels (exact span when
  κ ≤ 18, else sampled upper bound — provenance recorded); sector
  dispatch 2^k(base) (recorded; > 4096 flagged as beyond the
  demonstrated envelope, not verdict-changing); wall estimate =
  bottom_nodes / 1.1e9 nodes/s (A33 anchor: 6.6e10 ≈ 61 s) +
  20 s/rung overhead (A33/A36 anchors), declared order-of-magnitude
  (×3). The NUMERIC pre-registered claim is the node count, not the
  wall clock.
- **Verdict thresholds (A35 verbatim):** GREEN iff bottom nodes ≤ 2e11
  and cap ≤ 8; AMBER iff ≤ 1e14 and cap ≤ 12; else RED. REFUSED for
  odd |G|; NO-ROUTE if no k>0 fold exists.
- **Operative question:** frontier/parity rows price W at the window
  top (and separately at d_ub when different); anchors price W at
  d_exact − 2 (calibration) plus the blind window question.

Anchors' `g5_containment` self-test (does the route ceiling contain the
known d, computed WITHOUT consuming d_exact in the chain): 11 pass /
0 fail / 3 indeterminate (bounded-only chains).

## 5. Pre-registered pass/fail criteria (Phases 2–4 are scored on these)

(i) **Rank-law violations = 0.** Every Phase-2 engine-measured rung on
    cohort towers must satisfy rank p* = rank τ* = k(cover)/2 and
    exactness at the cover. A single verified counterexample falsifies
    the sheet-SES layer.

(ii) **Node-count formula within ×1.1 of realized counts.** Realized
    census enumeration node counts (per coset base; shared-walk engines
    are normalized by their offset count, the A35 §6 accounting) must
    match Nodes(κ, W) at the executed (κ, W) within a factor 1.1,
    per level, on every executed cohort census.

(iii) **GREEN → closure rate ≥ 90%, RED → non-closure.** A GREEN row
    closes iff Phase 2 produces the certificate-tier floor at the
    priced W within 10× the recorded wall estimate. A RED row closing
    at certificate tier within its stated gate budget falsifies the
    gate calibration. AMBER rows score neither way. Refused rows: any
    successful free-Z₂-deck construction on an odd-|G| cohort row
    falsifies C1.

(iv) **Certificate-vs-SAT agreement = 100%** on the ground-truth
    subsample: wherever Phase 3 obtains a solver-exact d for a cohort
    row (or a quotient), every certificate-tier floor/value produced by
    the descent pipeline must agree (floor ≤ d; exact = d). Any
    certified floor exceeding a verified SAT witness weight falsifies
    soundness.

(v) **Equal-compute head-to-head (Phases 2–4 design).** Paired design
    on the same cohort rows: Arm A = the descent pipeline exactly as
    priced here (frozen routes/W); Arm B = tuned SAT (the corpus's
    production stack: cryptominisat/cadical + tandem-maxsat staging),
    given the SAME total CPU-seconds per row that Arm A actually spent
    (min 60 s), same machine class, both nice'd. Score: per-row closure
    (exact d or floor ≥ window top) within budget; report closure
    rates, wall times, and the rows only one arm closes. Retrodiction
    rows (see FOLD_IN.md) are excluded from foresight scoring and
    reported separately as calibration.

## 6. Adaptations and deviations log

1. Screen output/paths: the two-line path adaptation described in §1
   (no functional change; all hard asserts kept).
2. `TS.BBCode` cannot represent k = 0 codes (empty-pairing crash);
   ladders truncate at the last k > 0 level with the reason recorded.
   No cohort row lost its whole route to this (`depth_used ≥ 1` on all
   120 in-scope rows).
3. μ (lightest stabilizer) is an exact span minimum only when
   κ ≤ 18; otherwise a sampled upper bound (provenance recorded).
   A lighter true μ would RAISE the fiber cap; Phase 4 should re-read
   μ from the executed censuses.
4. The A33 wall anchor (1.1e9 nodes/s) and 20 s/rung overhead are the
   stated conversion constants; wall estimates are ×3
   order-of-magnitude, not pass/fail quantities.
5. Depth-2 R4 pair-regime cell: unoccupied at freeze after 49 screened
   depth-2 towers — every depth-2 pair read dim(S∩K) ∈ {0, 2}, never 1
   (all k = 4/6 rows; R2 or R3 only). Logged as unfilled-from-pool; R4
   itself is covered at depths 3 and 4, so the regime lattice is
   globally occupied. Side observation worth carrying: a depth-2
   k = 4 pair with dim(S∩K) = 1 has not yet been seen; whether that is
   structural or small-sample is open. The order-144 fold-in adds
   depth-4 pair cells, not this one.
6. d_ub values for frontier rows are corpus L1-sampling values
   (A18 `a18_fill_ubs`, 30k–60k samples); they are witnesses, and some
   are loose — hence the window-top operative question rather than
   d_ub. Predictions never treat d_ub as a floor.

## 7. Immutability

`MANIFEST.sha256` lists SHA-256 hashes of every emitted file with an
ISO-8601 UTC timestamp. This directory is expected to be gitignored;
**committing the manifest (or the whole directory) is the user's
option** for stronger immutability — e.g. `git add -f
experiments/bb_lab/data/descent_theory_test/MANIFEST.sha256` in the
main repo, which pins the hashes without pinning the (large) data.
Nothing was committed by the prediction session itself.
