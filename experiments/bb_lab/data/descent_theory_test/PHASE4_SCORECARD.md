# descent_theory_test — Phase 4 scorecard

**Status: COMPLETE — all five pre-registered criteria scored.**

## Headline

| criterion | verdict | one line |
|---|---|---|
| (i) structural laws | **PASS** | 0 violations on 381 fresh rungs (rank law, exactness, Lemma 1) |
| (ii) node formula ×1.1 | **PASS** | 152/152 censuses at ratio exactly 1.0000 (2 D2 fallbacks excluded, listed) |
| (iii) GREEN ≥ 90% / RED 0 | **FAIL as registered** | GREEN 54/138 = 39.1%; decomposed: 74 exact-below-loose-ceiling + 10 envelope stops; resolution 128/138 = 92.8%; RED 0/11 closed; refusals 10/10 |
| (iv) SAT agreement | **PASS** | 17/17 GT solves agree (38/38 incl. anchors + tripwire); 0 disagreements; 3 cap-censored |
| (v) equal compute | **Arm A dominant** | 128/154 vs 47/154 closures; Arm B: ~4.5× the compute, 0 unique closures |

The single pre-registered failure, (iii), is a CALIBRATION failure of
the G5 ceiling (window tops loose at median gap 4; 0 containment
failures), plus an un-modeled candidate-set-volume axis in the cost
gates — not a soundness or cost-formula failure.  No RED row closed
(gate calibration held); no certified value was ever contradicted.

- Wall model: median realized/estimate 0.008, p90 0.354, max 7.2 — no
  closure lost to the 10× gate; 0 budget kills across all 262 executed
  questions.
- Net new certificate-tier values on 140 previously-open rows: **97
  exact distances** (82 counterexample-exact + 15 sandwich joins),
  **17 floors** (d ≥ 12–16), 26 UNKNOWN (envelope; both arms fail
  there).

Phases 2–4 executed 2026-08-17 in the same worktree that froze the
predictions (see PROTOCOL.md; MANIFEST.sha256 verified before any run —
sole diff `tools/dtt_lib.py`, the documented AMENDMENTS.md #1 amendment,
hash match confirmed).  All Phase-2+ code lives in `tools_phase2/`
(copies/adaptations only; frozen `tools/` untouched).  All runs nice'd,
caffeinated, ≤ 2 heavy processes, durable state only under this
directory.

## Closure engine (Arm A)

`tools_phase2/dtt_close_lib.py` — the A30/A32/A33 architecture,
generalized to run any cohort row along its FROZEN route:

1. bottom-level coset-BZ census over all 2^k sectors at the priced W
   (`bb_lab.cosetbz`, counting-invariant two-window BZ; exact node-count
   asserts ON, every run);
2. rung-by-rung transport: every cover cycle of weight ≤ W is either a
   bounded-overflow lift of a censused base cycle (complete MITM fiber
   enumerators, caps ≤ 8 = the A35 demonstrated envelope) or τ(γ) with
   γ ∈ ker E (checked ker E = base cycle space per rung at runtime);
   translation-orbit compression with per-rung covariance spot-asserts;
3. tiny levels (cycle-space dim ≤ 22) replaced by direct full
   enumeration; where both lanes ran, census-vs-full-enum orbit sets
   asserted EQUAL;
4. top level: no non-stabilizer cycle ≤ W ⟹ certificate-tier floor
   (d ≥ W+1, +1 under even-cycle parity); otherwise the minimum-weight
   logical found is a verified COUNTEREXAMPLE to the floor at W — and,
   because the transport is weight-complete below W, its weight is the
   EXACT d at certificate tier.

Trust tier of every Arm-A result: certificate (counting-invariant
enumeration + deck transports; not kernel-checked).  Witness weights are
never floors; UNKNOWN stays UNKNOWN.

### Pre-registration deviations (all recorded before scoring)

- D1 **Budget rule**: per-question kill budget = min(10 × predicted
  wall, 30 min), not the tasking's 3×: PROTOCOL (iii) scores closure
  "within 10× the recorded wall estimate", and a 3× kill would censor
  the pre-registered scoring window.  PROTOCOL wins over tasking.
- D2 **Single-window fallback census**: one cohort row's bottom code has
  a degenerate column matroid (no two disjoint information sets); the
  engine falls back to ONE window at r = W (complete unconditionally
  since |c|_I ≤ |c|).  Fallback censuses are excluded from the
  criterion-(ii) formula sample (different r-pair semantics) and listed
  below.
- D3 **r = 0 window skip**: W = 2 questions have pair radius r2 = 0; the
  C kernel over-walks that edge (25 ≠ 0 node assert), so the r = 0
  window's walk is skipped — its only contribution (the coset-base
  elements) is collected directly.  Node accounting unaffected
  (formula's r2 = 0 term is also 0).
- D4 **Memory envelope**: level candidate-sets are capped at 4e6 vectors
  (pre-flight binomial estimate + runtime guards).  The A35 gates bound
  census WALK nodes and fiber caps but not candidate-set volume; rows
  stopped by this guard are recorded BLOWUP_STOP and scored as
  non-closures.  This is an engine-envelope finding, not a soundness
  event.

## Criterion verdicts

### (i) Rank-law / structural violations — **PASS (0 violations)**

Every closure run re-measured its rungs with the frozen screen's own
functions (`h1_map` ranks, exactness spans, `sigma*`, `codim_lift`).
641 rung measurements over 381 unique (row, rung) cells:

- rank p\* = rank τ\* = k(cover)/2: **0 failures / 381**
- exactness at the cover (im τ\* = ker p\*): **0 failures / 381**
- Lemma-1 deck-constructor asserts: 0 failures (every deck built)
- `ker E = base cycle space` (the τ-branch scope check this engine adds):
  **held on every rung of every executed tower** — no row needed the
  general ker-E lane.

The sheet-SES layer survives contact with 154 fresh towers.

### (ii) Node-count formula within ×1.1 — **PASS (152/152 at ×1.0000)**

Every executed bottom census (the criterion's "executed cohort census")
ran the two-window coset-BZ walk with its exact node-count assert ON:
realized nodes = Σ_{s≤r₁} C(κ,s) + Σ_{s≤r₂} C(κ,s), per walk.

- Executed standard censuses: **152**; realized/predicted node ratio
  **1.0000 on all 152** (min = max = 1.0).  The frozen per-question
  `bottom_nodes` were reproduced exactly because the realized bottom
  κ equaled the frozen κ on every row (the engine asserts frame
  equality per level before running).
- 26 bottoms small enough to full-enumerate were ALSO censused and the
  two orbit-sets asserted EQUAL (independent completeness check).
- 2 censuses (row `0079a387…`, both questions) used the D2
  single-window fallback (degenerate column matroid at the bottom code;
  no disjoint information-set pair exists).  Their walks matched their
  own exact node counts but not the two-window formula (×13 / ×53.4) —
  excluded from the formula sample BY THE DEVIATION RULE and listed
  here rather than averaged in.

### (iii) GREEN closure ≥ 90% / RED non-closure — **AS REGISTERED: FAIL**

Scored on operative questions, retrodictions excluded (they were never
run), scope-controls counted only as refusals.

| stratum | GREEN rows | closed@W (floor) | exact-d below top (cex) | engine stop |
|---|---|---|---|---|
| anchor | 14 | 14 | 0 | 0 |
| frontier | 96 | 28 | 63 | 5 |
| parity-strip | 10 | 2 | 7 | 1 |
| order144 (batch 2) | 18 | 10 | 4 | 4 |
| **total** | **138** | **54 (39.1%)** | **74** | **10** |

- **Pre-registered reading: FAIL** (39.1% < 90%).  What failed is the
  G5 pricing of the operative question ("certify d = window top"): on
  74 GREEN rows a verified logical strictly below the window top
  exists.  The windows themselves never failed containment (d ≥ 2 and
  ≤ hi in all 122 resolved rows; gap histogram hi−d:
  0→48, 2→18, 4→27, 6→15, 8→13, 10→1).  Because the transport is
  weight-complete below W, each counterexample is an exact
  certificate-tier d — cross-checked exactly against corpus d on 6/6
  anchor window-questions and against an independent SAT solve on the
  first frontier cex row; full ground-truth check in criterion (iv).
- **Companion (post-hoc, labeled as such): certificate-tier resolution**
  (closed-at-W or exact-d-below-top, within 10× wall): **128/138 =
  92.8%**.
- 10 GREEN rows stopped at the engine's candidate-set memory envelope
  (D4) — an axis the A35 gates (walk nodes ≤ 2e11, cap ≤ 8) do not
  model: 5 frontier census-hit blowups (W=14–18 at n=108–168 bottoms),
  1 parity, 4 batch-2 mid-level cycle-set blowups.  These are honest
  GREEN failures of the COST-GATE layer, distinct from the window
  layer.
- **RED: 0/11 closed** — every RED operative stopped at the envelope
  exactly as gated (caps 13+ / set volumes 1e8+); gate calibration
  holds.  **No RED closed cheaply** (no conservatism event on operative
  questions; 28 RED secondary questions also all stopped).
- AMBER: 0/5 closed (reported separately, scores neither way).
- Refusals: **10/10** odd-|G| scope-control rows mechanically refused
  (both axes odd ⟹ no index-2 subgroup ⟹ no free Z₂ deck).
- Secondary questions (84 batch-1 + 38 batch-2 rows carried one):
  GREEN→cex 26, GREEN→stop 21, AMBER→stop 31, RED→stop 28 — REDs and
  AMBERs behaved as gated; GREEN '=d_ub' secondaries mostly re-found
  the same exact d (cex) or hit the same envelope.

### (iv) Certificate-vs-SAT ground truth agreement — **PASS (17/17, 0 disagreements)**

Independent SAT (CMS ladder, no window hints beyond d_ub, 3600 s caps,
2 nice'd workers) on a 20-row subsample stratified over v2-depth ×
stratum × Stage-A outcome lane (selection in
`phase3_gt_selection.json`; results in `phase3_groundtruth.jsonl`):

- **17 exact SAT solves — every one agrees with the Stage-A
  certificate**: 7 certified floors confirmed (d_sat = floor on the
  floor-attained rows, e.g. `00b630c1…` d=12 after a 1685 s solve;
  `23206eb7…` [[288,8,12]] d=12 after 560 s), 9 counterexample-exact
  rows confirmed (d_sat = cex weight, incl. a [[288,·,12]] batch-2
  row), and the dual-leg parity row `332f839e…` (floor 4 + cex 4)
  confirmed at d = 4.
- **0 disagreements** — no certified floor exceeded any solver value
  anywhere (also re-checked on all 20 cross-stage join candidates and
  on the 6 anchor window questions: nothing beat a floor).
- 3 bounded-only rows timed out at the 3600 s cap (SAT floors
  11/11/12 at d_ub 20/26/36) — the same rows Stage A could not resolve
  (candidate-set envelope); censored, no verdict either way, reported
  here rather than averaged in.
- Earlier plumbing note: the first GT batch attempt errored instantly
  (macOS EPERM from stacked `nice` + `os.nice`); fixed in
  `tools_phase2/` spawn args, errored attempt preserved at
  `tools_phase2/logs/phase3_groundtruth_error_attempt.jsonl`.
- GT arm compute: 21,718 s wall total (2 workers, includes the 3
  full-cap timeouts).

Adding the Stage-A-internal cross-checks (6 anchor window
counterexamples = corpus d; 14 anchor calibration floors = corpus d;
1 frontier SAT tripwire), the certificate machinery is now
solver-confirmed on **38/38 comparable values with zero
disagreements**.

### (v) Equal-compute head-to-head — **Arm A 128/154 vs Arm B 47/154; Arm B uniquely closed 0 rows**

Paired design on the 154 rows Arm A ran (PROTOCOL (v)): Arm B =
SAT-only (the corpus CMS ladder via `bb_lab.sat_distance`, d_ub as the
only hint), per-row budget = max(Arm A's realized CPU on that row,
60 s), enforced at wall × 1.15 + 5 s, cheapest-first scheduling, 2
nice'd workers.  Data: `phase3_armB.jsonl` (2 pause-superseded records
retained + re-run, see Incidents).

| | Arm A (descent/certificate) | Arm B (SAT-only) |
|---|---|---|
| compute | **2,019 s CPU** (getrusage; median 2.6 s/row) | **9,168 s** uninterrupted wall (budget allowance 9,827 s per the 60-s floor) |
| rows resolved/closed | **128/154** (83.1%): 54 floors@W + 74 exact-below-top | **47/154** (30.5%), all exact solves |
| unique closures | **81 rows only Arm A** | **0 rows only Arm B** |
| by stratum (closed/n) | anchor 14/14, frontier 91/96, parity 9/10, batch-2 14/34 | anchor 5/14, frontier 37/96, parity 5/10, batch-2 0/34 |
| floor depth when not closing | 17 floors d ≥ 12–16 at cert tier | timeouts' UNSAT floors median 8, max 11 |

- Headline ratio: **2.7× more rows closed on ~4.5× less compute**
  (≈ 12× per-CPU-second at these budgets); Arm B never closed a row
  Arm A missed.  Both arms fail together on the deep bounded-only rows
  (the 26 Stage-A envelope stops include all 3 GT cap-timeouts).
- Fairness notes, both directions: (a) Arm B's stack here is the CMS
  ladder only — the corpus's tandem-maxsat second stage was not wired
  (recorded limitation; on THESE budgets it would not plausibly flip
  n = 288 rows that need > 10× the allowance, but it caps how far the
  ratio generalizes); (b) the 60-s floor GIFTED Arm B ~4.5× Arm A's
  actual spend, per the pre-registered rule; (c) Arm A row CPU includes
  its structural-measurement overhead.
- Arm-A resolution here = certificate floor at priced W or exact-d
  counterexample (both certificate tier).  On Arm A's 54 floor rows the
  claim is d ≥ W+slack (exactness needs the witness leg — see the
  synthesis section); Arm B "closed" = exact d or floor ≥ window top
  (its pre-registered rule).

## Controls

All engine-family controls PASSED (final states; two early attempts of
the false-floor controls failed on a subprocess-spawn transient and were
re-run — the failed attempts remain in the JSONL as append-only
history):

- **false-floor-must-not-certify** (coset-BZ/tower family): "certify
  d ≥ known d + 2" on the gross code ([[144,12,12]], frozen route
  y,x,x) and pair72 ([[72,4,8]]) both ended COUNTEREXAMPLE at exactly
  the known d (12 / 8), never a floor — PASS.
- **planted-logical-found** (completeness): the weight-12 gross logical
  found above was re-planted and the pipeline's top-level orbit set was
  asserted to contain it — FOUND, PASS.
- **sat-lane-witness-verified** (arm-B family): pair72 x_distance = 8
  with witness verified in a nontrivial coset — PASS.
- Refusal battery: 10/10 odd-|G| rows refused (criterion (iii) row).

## Calibration (anchors + retrodictions)

- **Anchors (blind engine vs corpus-exact)**: 14/14 calibration
  questions (W = d_exact − 2) certified the floor at exactly d_exact;
  6/6 anchor window-questions priced ABOVE d_exact returned
  counterexamples at exactly the corpus d.  The engine independently
  reproduces known distances in both directions (floor side and
  witness side).
- **Batch-2 retrodictions (24 rows, calibration-only, never run)**:
  24/24 blind G5 windows contain the sweep-known d; window_hi − d gap:
  min 0 / median 0 / max 8 — the depth-4 chains' ceilings are tight
  where the quotient distances were freshly exact.
- **Batch-1 frontier windows are the loose ones**: gap histogram above
  (median 4).  All frontier window terms were chain-exact
  (`ceiling_is_upper_estimate` = false on all 96), so the looseness is
  entirely the 2·d(base) doubling bound not being attained — a theory
  calibration finding, not a data-quality artifact.
- **Parity-strip**: machinery ran soundly on all 10 rows (decks,
  transports, rank law all held — included in criterion (i) counts).
  The frozen prediction "odd d values allowed" was not exercised: all
  9 resolved rows returned EVEN d (4,6,8,8,8,10,10,10,14).  Odd-weight
  CYCLES exist at every level (as frozen), but minimum-weight logicals
  landed even anyway.  W_eff = d−1 pricing behaved as specified.

## Incidents

- **Compute pause during Arm B (user-requested).**  The arm-B driver
  and its two in-flight SAT workers were SIGSTOP'd by the main session
  mid-run and SIGCONT'd ≈ 3 h later.  The driver resumed cleanly; the
  two in-flight workers (rows `00dfb7f9a695…`, `00933835c472…`) died at
  resume when suspended wall-clock enforcement fired spuriously (their
  records show wall ≈ 10,600 s against 60 s budgets).  Per instruction:
  those two records are marked `superseded` in `phase3_armB.jsonl`
  (kept, not deleted) and the rows were RE-QUEUED fresh after the
  driver's main pass.  **Accounting basis for criterion (v):** per-row
  budgets are enforced on the wall clock of uninterrupted execution
  (budget × 1.15 + 5 s slack); the solver is single-threaded and
  nice'd, so uninterrupted wall upper-bounds CPU — no row other than
  the two superseded ones spans the pause, arm-A CPU was measured by
  `getrusage` (immune to suspension), and the contaminated batch-level
  wall clock is not used anywhere in scoring.

## Censoring accounting

- Budget kills (hard or soft): **0** across all 260 executed questions
  (largest realized wall 149 s vs cap 1800 s).
- Engine-envelope stops (BLOWUP_STOP / ENVELOPE_STOP): 106 questions
  (operative: 26 — see criterion (iii); secondary: 80).  Every stop
  carries the level and the estimated/realized set size in
  `phase2_results.jsonl`.
- No row was left unattempted; the two Stage-A errors were re-run to
  completion after D2/D3 (below).

## Engine errors

Three question-level errors occurred on first execution, all fixed in
`tools_phase2/` (never in frozen `tools/` or program `scripts/`) and
re-run to green:

1. `0079a387… q0/q1` — `disjoint_info_sets` found no disjoint pair
   (degenerate bottom column matroid; A = 1+y+y² collapses onto the
   folded axis).  Fix = D2 single-window fallback (complete
   unconditionally).  Re-run: floor 4 certified at W=2; W=12 question
   resolved.
2. `00c0a04e… q0` — C-kernel r=0 window edge (`node count 25 != 0`
   assert).  Fix = D3 skip-the-empty-window (its only member, the
   coset base, is collected directly).  Re-run: floor 4 certified;
   the W=16 secondary resolved as a counterexample in 143.6 s
   (wall ratio 7.2, the Stage-A maximum).
3. One control-script spawn transient (controls, not cohort rows) —
   re-run PASSED; both attempts recorded.

## Cross-stage synthesis (post-hoc join of two certificate-grade artifacts)

**Label: SYNTHESIS.**  This section joins Stage-A certified floors with
independently produced minimum-weight witnesses.  No new distance
computation was run beyond witness re-verification (deterministic
re-derivation of the sampler output + an F₂ check that the vector is a
cycle outside the stabilizer rowspace).  Floors were never derived from
witnesses; witnesses never tightened floors.  Full records with both
provenance legs and witness supports: `phase3_joins.jsonl`.

### Family A — order-144 sweep joins (batch-2 rows)

Join rule: Stage-A floor F (counting-invariant tower census) + sweep
`d_ub = F` with the sweep's verified-witness flag, witness re-derived
here with the sweep driver's exact sampler parameters
(`l1_distance_ub(n_samples=100000, seed=3)`) and re-verified from
scratch.  **5/5 candidates joined; 0 alarms:**

| row | code | Stage-A floor leg | witness leg | verdict |
|---|---|---|---|---|
| `23206eb7…` | [[288,8,·]] Z12xZ12 | floor 12 | w=12 re-verified | **d = 12 EXACT** |
| `6637e93c…` | [[288,4,·]] Z12xZ12 | floor 12 | w=12 re-verified | **d = 12 EXACT** |
| `cdf9a988…` | [[288,4,·]] Z12xZ12 | floor 12 | w=12 re-verified | **d = 12 EXACT** |
| `50a94227…` | [[288,4,·]] Z12xZ12 | floor 18 | w=18 re-verified | **d = 18 EXACT** |
| `ecf09d8e…` | [[288,4,·]] Z12xZ12 | floor 18 | w=18 re-verified | **d = 18 EXACT** |

The two d = 18 rows are solver-free floor certifications at n = 288
matching the program's A36-class results, now on fresh sampled codes.
Five further batch-2 floor rows (floor 16 vs sweep d_ub 24–62) do NOT
join and remain honest floors.

### Family B — corpus-d_ub joins (batch-1 rows; separate provenance)

Same sandwich logic, but the witness leg is a FRESH L1 hunt (60k
samples, ≤ 3 seeds) re-finding a weight-F witness, F₂-verified here —
the corpus d_ub value itself is never taken on faith.  **10/10
candidates joined; 0 alarms:** 8 frontier rows at d = 12, one frontier
row at d = 14 (`22aafb05…`), one parity-strip row at d = 10
(`5f81ecbb…`).

Synthesis yield: **15 additional exact distances** on top of Stage A's
counterexample-exact rows; every join lists both legs and the witness
support in `phase3_joins.jsonl`.  A witness below a certified floor
would have been a stop-and-investigate soundness event — none occurred
(the check ran on all 20 candidates).

## Wall-model calibration (informational, ×3 order-of-magnitude claim)

Realized wall / estimate over 262 executed questions: median 0.008,
p90 0.354, max 7.18.  The model (bottom_nodes / 1.1e9 + 20 s/rung) is
strongly conservative on this cohort — bottoms are tiny and the rung
overhead dominates the estimate while the engine's Python transport is
faster than the modeled constant.  No question needed even 8× its
estimate; the pre-registered ×3 order-of-magnitude claim holds
one-sidedly (over-estimation only).

## What this does and does not establish

**Established by this test:**

- The structural core of the descent theory — Lemma-1 transports, the
  rank law rank p\* = rank τ\* = k/2, exactness at the cover — survived
  381 fresh rungs without a single violation (criterion i), exactly as
  the sheet-SES argument predicts.
- The census cost formula is exact, not approximate: 152/152 executed
  censuses landed on Nodes(κ, W) to the node (criterion ii).
- The certificate machinery is SOUND in practice as well as by
  construction: 38/38 solver comparisons agree (criterion iv), no
  witness ever beat a floor, planted logicals are found, false floors
  are rejected.
- Under equal compute the descent lane dominates the SAT lane on this
  cohort: 2.7× the closures on ~4.5× less compute, zero SAT-only
  closures (criterion v).
- The odd-|G| scope boundary is real and mechanical (10/10 refusals).

**Refuted / miscalibrated as pre-registered:**

- Criterion (iii) AS WRITTEN fails badly (39.1% vs ≥ 90%): the frozen
  operative question "certify d = G5 window top" priced the wrong
  TARGET on most frontier rows.  The G5 window is an honest containment
  statement (0 containment failures, anchors 11/11 + retro 24/24 +
  122/122 Stage-A resolutions inside their windows) but its CEILING is
  attained only ~39% of the time — the 2·d(base) chain bound is loose
  at median gap 4 on the Z9xZ6-family frontier.  GREEN should be read
  as "the certificate lane will resolve this row cheaply" (92.8% true),
  not "d equals the window top" (39% true).
- The A35 cost gates miss an axis: 10 GREEN rows (and 26 rows overall)
  stopped on candidate-SET volume, which walk-node and fiber-cap gates
  do not measure.  A binomial set-size gate
  (Σ_{w≤W} C(n_j, w) · 2^{dim_j − n_j} per level) belongs alongside G1/G2
  in any successor screen.

**Not established / out of scope:** nothing here is kernel-checked
Lean; the certificate tier remains "counting-invariant enumeration +
deck transports, asserts on".  The equal-compute ratio is specific to
these budgets, this SAT stack (CMS ladder without the tandem-maxsat
stage), and this cohort's n ≤ 288 window; deep-d rows (the 26
envelope stops) defeated BOTH arms and remain UNKNOWN.  Batch-2
retrodictions never entered criteria (iii)–(v).  Parity-strip rows
exercised the machinery but not the odd-d scope (all resolved d even).
