# A39 — Pre-registered falsification test of the descent machinery (order-144 sweep + 164-row frozen-prediction cohort, Phases 0–4)

**Sessions 2026-08-17/18** (worktree `quantum-distance-presentation-73cb27`,
branch `claude/a39-descent-theory-validation`). Mission: test the
A32/A33/A35 tower-descent theory as a *predictive, quantitative,
falsifiable* theory of distance-certification cost and correctness —
pre-register the A35 pricer's per-code predictions, freeze them under a
hash manifest, then shoot at them with real closures, an independent
SAT ground-truth arm, and an equal-compute head-to-head. Companion
survey: a 120 s/code distance sweep over order-144 groups (the depth-4
stratum), folded in under an explicit retrodiction/prediction honesty
rule. **Feeds A38**: criterion (v) below is the F3 complexity ledger's
first controlled measurement; the Phase-1 machinery is F5's
screen-as-engine, validated bit-exact against the banked A35 anchors.

## §0 Verdict

| criterion (pre-registered, PROTOCOL.md §5) | verdict | number |
|---|---|---|
| (i) structural-law violations = 0 | **PASS** | 0 on 381 fresh engine-measured rungs |
| (ii) census node formula within ×1.1 | **PASS** | 152/152 at ratio **1.0000** (assert-verified) |
| (iii) GREEN closure ≥ 90% / RED non-closure | **FAIL as registered** | GREEN 54/138 = 39.1%; RED 0/11 closed ✓; refusals 10/10 ✓ |
| (iv) certificate-vs-SAT agreement = 100% | **PASS** | 17/17 GT solves agree (38/38 incl. anchors + tripwire), 0 disagreements, 3 censored @3600 s |
| (v) equal-compute head-to-head | **descent dominant** | 128/154 vs 47/154 rows; SAT-only uniquely closed **0**; n = 288: 14/34 vs 0/34 |

**Reading.** The theory survives every soundness-shaped criterion and
fails the one calibration-shaped criterion exactly where its window
pricing is loose. The (iii) miss decomposes: 74/138 GREEN rows resolved
with a verified counterexample *below* the priced window top — each an
**exact certificate-tier d** (transport is weight-complete below W), so
the G5 ceilings are loose (median gap 4) but were **never violated**
(0 containment failures); 10/138 hit an **un-modeled memory envelope**
(candidate-set RAM, cap 4·10⁶ vectors — storage, not node count, is the
unpriced axis, confirming A32 §9's remark with data). Post-hoc
companion rate: certificate-tier resolution 128/138 = 92.8% within the
10× scoring window. Wall model conservative: median realized/estimate
0.008, p90 0.354, max 7.2, zero budget kills.

## §1 Design

1. **Order-144 sweep** (companion survey, ran first): 58 sampled
   weight-3 codes, k ≥ 2, canonical dedup, over Z12xZ12 (34) /
   Z18xZ8 / Z24xZ6 / Z36xZ4 / Z16xZ9 / Z48xZ3, hard 120 s/code.
   Outcomes: 24 exact (13 solver CMS-ladder, 11 certificate via the
   doubling front-end — including three [[288,4,16]] and five
   [[288,·,12]]), 34 bounded with honest floors; base-ladder cache 167
   quotient codes (122 exact, 87 free corpus hits). Survey stat that
   did not exist before: **11/58 random order-144 codes certified as
   doublers within 120 s** (a lower bound; certify-lane scope gates
   k ≤ 14, census W ≤ 22, base d_ub ≤ 15). `data/order144_sweep/`.
2. **Phase 0 — stratified cohort (130 rows)**: 96 frontier (corpus
   rows with NO exact d, d_ub ≥ 12, n ∈ [108,168], depths 1–3), 14
   anchors (corpus-exact d ∈ {6..12}, blind calibration), 10 odd-|G|
   scope-controls (frozen prediction: REFUSE), 10 freshly sampled
   wt-3×wt-4 parity-strip rows. All four pair regimes R1–R4 occupied;
   **depth-2×R4 unfillable** from the pool (49 depth-2 towers all read
   dim(S∩K) ∈ {0,2} — open question whether structural at k = 4).
3. **Phase 1 — frozen predictions**: the A35 screen (extracted at
   `2f063b0`, deps byte-identical, zero asserts disabled) driven
   per-row; validation reproduced every banked anchor before
   predicting. Free results at freeze: rank law 0/245 rungs, anchors'
   blind G5 containment 11/11 decidable, parity-strip L3 symptom
   present 10/10. Sealed: `MANIFEST.sha256`, self-hash `abe57971…f55bf`,
   2026-08-17T17:27:38Z. Contamination rule: no exact-distance work on
   any target code; quotient d only by trivial full enumeration
   (n ≤ 40) or corpus lookup.
4. **Fold-in (batch 2)**: sweep rows ingested as the depth-4 stratum —
   24 retrodictions (outcome predates prediction; calibration only,
   excluded from (iii)/(iv)) + 34 genuine predictions;
   `retro_ceiling_fail = 0`; depth-4 regimes R1 11 / R2 113 / R3 40 /
   R4 10. Required one post-freeze amendment (AMENDMENTS.md #1): the
   depth-4 towers hit `TS._preimage`'s degenerate case (descended seam
   image spanning the full bottom H1 — annihilator empty, numpy loses
   2-D shape). Fixed in the test's tools layer only, semantics
   unchanged; **the upstream bug in `scripts/a32_tower_slice.py` is
   still open** (its fix session was deleted).
5. **Phases 2–4**: closures along frozen routes (budget = min(10×
   predicted wall, 30 min) per PROTOCOL deviation D1), GT arm (20
   stratified rows, independent SAT, ≤3600 s caps), equal-compute arm B
   (CMS ladder, per-row max(armA-CPU, 60 s)), scorecard + manifest.
   Engine: `tools_phase2/dtt_close_lib.py` — the A30/A32/A33
   architecture generalized to run any cohort row (coset-BZ censuses
   with exact node asserts, bounded-overflow MITM fibers caps ≤ 8,
   τ-branch, full-enum replacement at cycle-dim ≤ 22, per-rung
   covariance spot-asserts).

## §2 Yields

- **97 new certificate-tier exact distances** (82 counterexample-exact
  + 15 cross-stage sandwich joins, each join's witness re-verified):
  incl. two **[[288,4,18]]** (depth-4 towers matching the published
  record code's distance at n = 288), one [[288,8,12]], [[168,4,16]],
  nine d = 14 (three closed in < 8 s; the corpus's own maxsat campaign
  had left all of these rows open), and the parity-strip
  **[[120,4,14]]** (wt-3×wt-4 — exact without the parity layer,
  demonstrating A35 L3's "simplifier, not precondition" at d = 14).
- 17 certificate floors (d ≥ 12–16), 26 honest UNKNOWNs. By stratum
  (exact/floor/unknown): frontier 79/12/5, parity-strip 9/0/1,
  batch-2 9/5/20.
- **0/88 prior d_ub values were tight** — every L1-sampled upper bound
  in the band overestimated d (30→14, 32→14, 26→12 …). The witness
  side systematically fails here; quantifies why the floor factory was
  the missing half.
- Controls: false floors rejected (gross, pair72 negative controls),
  planted logical FOUND, 26 census-vs-full-enum cross-checks EQUAL,
  10/10 refusals, fiber empty-rate median 0.935 (reproduces the banked
  tower-win mechanism on fresh codes).
- Compute: Arm A 0.56 h CPU total (median 2.6 s/row); GT arm 6.0 h
  wall for 20 rows; Arm B 2.55 h wall; sweep 31.5 min.

## §3 Corpus merge (`scripts/a39_corpus_merge.py`)

Upsert of the A39 outcomes into `data/bb_instances.duckdb`, conflict-
abort semantics, existing `d_exact` never overwritten, solver-lane
sweep values keep solver method strings, certificate values get
`descent-cert@a39*` methods — **the corpus's first non-solver
`d_method` tier** (all 43,614 pre-existing exacts are
maxsat-tandem/cadical/CMS). Dry-run plan, verified clean:
**79 exact updates + 67 inserts (58 sweep + parity rows) + 12 floor
raises + 14 agreements + 0 conflicts**. Apply is one command (writes a
`bb_instances.a39bak.duckdb` backup first):

    cd experiments/bb_lab && uv run python scripts/a39_corpus_merge.py --apply

Status at commit time: dry-run validated; apply pending (session
permission boundary on the main-checkout DB).

## §4 Findings for the theory (beyond the criteria)

1. **G5 ceilings are systematically loose** (median gap 4, never
   violated): the certifiable-d window model needs a calibration term;
   until then, window-top questions overprice what closures actually
   deliver (exact d below top).
2. **Memory is the unpriced cost gate**: 10 GREEN rows censored by the
   4·10⁶-vector candidate-set envelope with node counts exactly as
   predicted. A G6 (storage) gate alongside G1–G5 is owed.
3. **Equal-compute dominance is total on this cohort**: SAT uniquely
   closed zero rows — everything SAT could do, descent did, plus 81
   rows SAT could not touch at the same budget.
4. Depth-2×R4 emptiness; the Z48xZ3 k=4 row (d ∈ [10,32]) still open —
   the cohort's most interesting unresolved instance.
5. Ops lesson: SIGSTOP-pausing solver workers interacts with wall-clock
   kill budgets (two arm-B rows died at SIGCONT after a ~3 h pause;
   re-run, records superseded, accounting basis = getrusage CPU for
   arm A / uninterrupted wall for arm B — incident logged in the
   scorecard).

## §5 Reproduction map

| claim | check |
|---|---|
| frozen predictions + protocol | `data/descent_theory_test/{PROTOCOL.md, predictions*.jsonl, MANIFEST*.sha256}` (batch-1 self-hash `abe57971…f55bf`) |
| screen provenance | `data/descent_theory_test/tools/` + `validation/screen_run.log` (A35 anchors reproduced; extracted at `2f063b0`) |
| Stage-A closures | `data/descent_theory_test/phase2_results.jsonl` (288 records) + `tools_phase2/` engine |
| GT arm / arm B / joins | `data/descent_theory_test/phase3_{groundtruth,armB,joins}.jsonl`, selection in `phase3_gt_selection.json` |
| scorecard + incidents | `data/descent_theory_test/PHASE4_SCORECARD.md`, `MANIFEST_phase2.sha256` (18 files) |
| amendments | `data/descent_theory_test/AMENDMENTS.md` (#1: `_preimage` degenerate-case guard, hashes before/after) |
| order-144 sweep | `data/order144_sweep/{results.jsonl, REPORT.md, *.py}` |
| corpus merge | `scripts/a39_corpus_merge.py` (dry run prints the full plan; `--apply` backs up first) |

## §6 Owed / next

1. Run the corpus-merge apply (one command, §3).
2. Upstream `_preimage` fix + regression test in `a32_tower_slice.py`.
3. G5 ceiling calibration model; G6 memory gate (A38 F3/F5 material).
4. Z48xZ3 [[288,4,d∈[10,32]]] closure (needs the memory-envelope fix
   or a deeper lane).
5. Lean packaging pick: one sandwich join ([[288,4,18]]) as the first
   kernel-checked cross-stage exact.
6. Excluded-from-commit bulk: `tools_phase2/{logs,work}/`,
   `order144_sweep/{certify_runs, sweep.duckdb, smoke_*}` — local only;
   regenerate via the reproduction map.
