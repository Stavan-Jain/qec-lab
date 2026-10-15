# A41 S1 — the gates session (all-empirical; no new distance claims)

**Session 2026-08-25** (worktree `distance-descent-theory-plan-e1fd17`,
branch `claude/distance-iterated-z2-covers-9c837c`). Charter:
`A41_closed_form_family_plan.md` §7 S1. Discipline: A38 §6.0 (gate
first, claim tiers exact, witness weights never floors, no SAT on any
floor path — this session runs NO floor-side SAT at all).

Scripts: `a41_g1a_assert.py`, `a41_g3_tierd_screen.py`,
`a41_w7_census.py`, `a41_q6_selfsim_probe.py`; data `data/a41/`.
Parallel: A40 (tour-de-gross) runs in its own worktree; the teaching-doc
update chip runs separately; no file overlap.

## §0 Gate 0 — banked regression

`bb_lab.tower.validate_banked('data', rng=20260811)`: **GREEN** —
banked A32/A33 structure, census node anchors, A19 deck-survey
k-verdicts, and the 11-tower field-identical screen all reproduced.

## §1 G5 — cross-worktree audit (the w = 5 state was ahead of the plan note)

- **The split map lives on unmerged branch
  `worktree-agent-aa0d7676f1aee3b4c`** (A5 log Entries 16–25 +
  durability postscript, 2026-08-17): (1,k≤5), (2,2), (3,3)/P-33,
  (2,4)/P-K4, and all even splits CLOSED; (4,4) profile table proven +
  clean-room adversarial EXACT MATCH; (1,7)/(3,5) profile layers
  census-CLEAN ((1,7) structure-empty on the live population); the
  partition-existence/coverage kill named as the class's leading
  profile-killer (4 profiles across 3 splits).
- **The P-26 / (4,4)-strata production is LOST**: the sharded relaunch
  (20 shard-jobs + K6,K6) was re-parented to launchd-style detached
  processes writing into that session's worktree — the worktree has
  since been deleted, no surviving processes (`ps`/`launchctl` clean),
  and no outputs in the main checkout. What survives in git: the
  engines (`a17_e20_p26_engine.py` sharding + terminal streaming,
  `a17_e21_p26_merge.py`, `a17_e22_44_engine.py`), their validation
  (a13 4-shard byte-exact reproduction), and Entry 23/25's counts
  (a11: ≥ 560 abstract terminals, first-pass all rank-8 FREE).
  **Q1 residue = re-run, not re-derive.** Ops lesson re-learned: the
  A38 checkpoint rule applies to *worktree lifetime* too — long
  detached runs must write outside the worktree or the worktree must
  outlive them.
- **A37's spread lemma is now any-weight** (c953764, unmerged
  `worktree-agent-a7abc8e521d53153e`): Lemma A37.2 re-proven via
  injective cell-to-pair charging for every |A|, |B|, with the
  explicit note that the D1/D2 half ports to w = 5 verbatim (parity
  round-up still needs odd weights). Strengthens charter Q5 route (i).
- **`SidonConvBound` audited reachable**: QECLean local branch
  `a17/sidon-conv-image-bound` @ fa62959 ("first analytic class-level
  lemma in Lean"), unmerged and on no remote — upstreaming owed.

## §2 G1a — the w = 5 census, regenerated and asserted

Protocol archaeology first (falsify-first): Entry 15's numbers came
from a **structural census** (`--sat-cap 0`, live frames = both axes
≥ 5, 4∤ either, |G| ∈ [41, 78] complete + 6 frames "in flight" never
recorded) + a **Tier-1 exact SAT screen at ub = 6** on the 41
small-frame members + the **MC witness hunt** (wmax 9, 800
iters/member). The banked "450 exactly at d = 10" is therefore
**witness-certified d ≤ 10 + the conjectured class floor** — not
per-member solver-exact (one member had UNSAT rounds to d ≥ 9; one
Z₇×Z₉ member has the 7.5 h coset-ladder floor, commit 56e2cf6). The
original sweep jsonl was never banked (data/ is gitignored; it lived
in the main checkout and is gone).

Re-run this session (`data/a41/t42_w5_census.jsonl`,
`t42_census.log`): **every one of the 15 banked Entry-15 frames
reproduced EXACTLY** — all seven fields (poolA/poolB/zdA/zdB/d2/kpos/
members) per frame, including Z₅×Z₁₅'s 42,300/2,460/30,540/1,500/
691,008/67,296/**2,103** — and the banked-scope member total is
**2,144** on the nose (`a41_g1a_assert.py`: GATE GREEN). The [41,60]
Z₂/Z₃-axis frames additionally re-measured empty this session (chunk
runs; Lemma E corroborated). Tail frames (9x9, 6x14, 5x17, 6x15,
9x10, 5x18 — Entry 15's "in flight") are NEW data: results appended
below when the run lands.

MC falsifier gate (wmax 9, 800 iters/member, seed 20260712, re-run
this session over all 2,144 members, 399.5 s at 8 workers): **ZERO
falsifiers** — Entry 15's headline verdict reproduced. Spectrum vs
banked: **450 @ 10 EXACT**, 28 @ 12 exact (= 11+14+2+1 across
frames), 75 @ 14 exact (= 52+4+5+14); the stochastic high buckets
wobble within bucket (684/902/5 vs banked 681/906/4 at 16/18/20) — an
expected seed-order artifact: per-member seeds are index-offset and
this census orders frames differently, and 800 iterations saturate
only the low weights (the orbit-boost argument). Gate reading: the
witness-side data is bit-stable exactly where it is load-bearing
(≤ 14) and bucket-stable above.

## §3 G4 — the k census (banked-scope snapshot)

| frame | members | k |
|---|---|---|
| Z₆×Z₉ | 18 | 4 |
| Z₆×Z₁₀ | 8 | 8 |
| Z₇×Z₉ | 15 | 12 |
| Z₅×Z₁₅ | 2,071 + 32 | 8 (2,071), 16 (32) |

k ∈ {4, 8, 12, 16}: frame-pinned except Z₅×Z₁₅'s 8/16 split. Matches
Entry 15's per-bucket splits (450@10 = 436 k=8 + 14 k=16). Formula
candidate: k = 2·dim Ann(A,B) is the definition; the *frame* fixes the
reachable CRT components (Entry 15's field mechanism), so a closed
form needs the Q3 construction to pin which components vanish — a
design freedom, not an obstacle.

## §4 Q3 — the inhabited shape (new structural finding)

Support-shape taxonomy over the 2,103 Z₅×Z₁₅ members:

- **B's y-row count pattern is (1,2,2) on 2,103/2,103** — uniformly:
  one singleton row + two doubleton rows.
- A's x-column patterns: (1,2,2) × 1,733, (2,3) × 315, (1,4) × 55.

The dominant (and B-side exclusive) shape is **singleton +
(w−1)/2 doubletons** — exactly the gross shape one weight up
(gross A = x³ + y + y², pattern (1,2)). Q3's constructive ansatz is
therefore

    A_w = x^{a₀} + Σ_{i=1}^{(w−1)/2} x^{a_i} (y^{b_i} + y^{c_i}),
    B_w = mirrored,

with D1/D2/(a′)/non-Frobenius imposed on the exponent data and the
Bezout witness (charter §3 lever 1) built into the doubleton
structure (the gross witness (1+x²)B² = 1+x⁶ telescopes through
exactly this shape via Frobenius squaring). Designing and enumerating
fresh members = S2-facing residue; the ansatz is now data-backed
rather than guessed.

**The S2 seed cohort** (joined filters: mc_min = 10 ∧ k = 8 ∧ y-axis
CHEAP-PASS ∧ (1,2,2) A-shape): **296 members**. Canonical first —
the named S2 starting member:

    A = 1 + y + x·y² + x·y¹⁰ + x²·y⁶
    B = 1 + y² + x + x²·y³ + x³·y²

a [[150,8,10]]-candidate whose designed y-rung targets [[300,8,20]] —
note the resonance with the Lean-proven cover300 pair
([[150,8,8]] → [[300,8,16]]): same n and k, two more units of
distance at both levels if Tier D lands.

## §5 G3 — Tier-D screen (cheap tiers, floor 20 = 2·(2w))

Driver `a41_g3_tierd_screen.py` (wraps
`a17_corpus_battery.process_cell`, run_s4=False — k-gate + S0/S1+/S2
only, no SAT). Claim discipline: verdicts are SF@20 statements for the
stored presentation, never "doubles"/"does not double".

Small-frame population (41 members × both axes, **4 s total** — the
cheap tiers are effectively free at n ≤ 252):

| frame | axis | verdicts |
|---|---|---|
| Z₆×Z₉ | x (even) | 13 CHEAP-PASS, 5 CHEAP-REJECT |
| Z₆×Z₉ | y (odd) | 15 CHEAP-PASS, **3 K-GATE-FAIL** |
| Z₆×Z₁₀ | x (even) | 7 CHEAP-PASS, 1 CHEAP-REJECT |
| Z₆×Z₁₀ | y (odd) | 8 CHEAP-PASS |
| Z₇×Z₉ | x, y (both odd) | 15 + 15 CHEAP-PASS |

Signals: (i) the population is largely SF@20-viable at the free tiers
— Tier D is not dead on arrival; (ii) every CHEAP-REJECT so far is on
an even axis, consistent with the charter's odd-axis design rule —
and the reject values are **wall-shaped**: the five Z₆×Z₉-x cells all
sit at cheap_min 18 = 2·10 − 2 (the A17 wall value) and the Z₆×Z₁₀-x
cell at 16 (deficit 4); (iii) three in-the-wild (R)-failures on an
odd axis (Z₆×Z₉ y, k 4 → 8: the k-doubling regime) — the k-gate is
not automatic even off the 2-part, so the Bezout-by-construction
lever is load-bearing, not decorative.

**Z₅×Z₁₅ (2,103 members × both axes, 4,206 cells, 232 s):**

| axis | CHEAP-PASS | CHEAP-REJECT | K-GATE-FAIL |
|---|---|---|---|
| x (double the 5-side) | 2,056 | 45 | 2 |
| y (double the 15-side) | 2,087 | 1 | 15 |

Cross-tab against the MC spectrum — **the 450 tight (mc_min = 10)
members, i.e. the actual d = 2w family candidates**:

| axis | pass | reject | (R)-fail |
|---|---|---|---|
| x | 404 | **44** | 2 |
| y | **442 (98.2%)** | **0** | 8 |

S1 verdict for Tier D: the designed rung is the **y-axis** — zero
cheap-tier SF@20 failures among the tight candidates there, with (R)
failing on only 8/450; every tight-member reject lives on the short
x-axis. The class-wide safe floor is empirically near-uniform on the
right axis — exactly the precondition a family-wide certificate
needs. (Free tiers only: SF@20 survival is S4/BZ-certifiable per
member but not yet certified; SF stays sufficient-not-necessary.)

## §6 G1b — the w = 7 census engine + first frame

`a41_w7_census.py`: the w5 gate stack with the pool enumerator
replaced by a DFS with incremental difference pruning (the sweep's
`combinations` pool is C(|G|−1,6) ≈ 4·10⁸ at |G| = 90 — infeasible),
one DFS feeding both axes' pools, and a bit-packed GF(2) rank for the
zd filter and k-test (the uint8 rref at 3.8M candidates would be
hours). Scope recorded: both axes ≥ 5, |G| ≥ 85 (D1∧D2 bound 84);
small axes excluded by cost — the w5 Lemma-E kill is weight-specific
and NOT re-derived.

Selftest: DFS pools == the reference combinations pools at w = 5 on
Z₆×Z₉ (both axes, exact set equality), packed rank == `rank_f2` on
300 random + padded-hstack cases. **Session catch (falsified claim of
my own)**: the first DFS version admitted supports with an order-2
difference — the ordered difference multiset holds those twice, so
they are never Sidon; the selftest's 25,860 ≠ 19,980 mismatch caught
it before any production use.

**First four frames — all EMPTY, two distinct mechanisms:**

| frame | G | verdict |
|---|---|---|
| Z₅×Z₁₇ | 85 | poolB = 0 — **no mono-y Sidon 7-set exists** on a 5-column frame (B-side death; the B-first ordering skipped the 546k-class A filter) |
| Z₆×Z₁₅ | 90 | poolB = 0 — same 6-column death |
| Z₉×Z₁₀ | 90 | pools rich (155,688 A / 95,544 B classes), zd rich (29,016 / 38,340) — **d2_pairs = 0**: D1∧D2 needs 84 of the 89 nonzero differences, disjointness never packs |
| Z₅×Z₁₇ A-pool datum | — | raw Sidon-7 pool 3.82M supports in 284 s (the DFS is fine; the packed rank removed the filter wall) |

Structural reading (new, and it revises the charter's Tier-F scaling
footing): **at w = 7 the counting bound |G| ≥ 85 is not achievable**
— near-bound frames are D2-empty even when zd-rich, and thin axes
(5, 6 columns) cannot host mono Sidon 7-sets at all. At w = 5 the
members appeared at slack ×1.32–1.83 over the bound (|G| = 54, 75),
so inhabitation at w = 7 is expected (if anywhere) at |G| ≈ 110–160
on frames whose axes carry the good component fields (Entry 15's
mechanism: 5, 7, 9, 10, 15-type axes = GF(8)/GF(16)/GF(64);
11/13/17-type axes were empty at w = 5). Tier F's d = Θ(√n) survives
iff the slack stays a constant factor — now an explicit Q3 question,
not an assumption. Engineering landed for the fat frames: vectorized
D2 pairing (packed-mask block-AND — the Python pair loop dies at
~10⁹ candidates). Good-field probes (7x15, 5x21, |G| = 105) in
flight; results appended below.

## §7 Q6 — the self-similarity probe: NO SIGNAL (park)

`a41_q6_selfsim_probe.py` (banked data only;
`data/a41/q6_selfsim_probe.json`):

- The only orbit-grade level pair banked anywhere is c37xx L2 (n=180)
  vs L1 (n=360), stab censuses to W = 22. Best affine reindex is the
  **shift w_L1 = w_L2 + 6** (not the doubling w ↦ 2w) with ratio
  spread **2.56** over 5 bands — nothing close to proportional
  (threshold 1.5). μ-sequence is flat (6, 6, 6); d-sequence
  10, 10, 20, ≥24.
- **Verdict: park Q6** per the charter's S1 call, with the recorded
  caveat that the evidence base is thin (exactly one comparable pair;
  banked multi-level orbit censuses are scarce). Unbanking more level
  pairs is the cheap way to strengthen this if it ever matters.
- Cross-tower side-find (regularity, NOT self-similarity): bb288-mid
  (6,12) and c37xx-L2 (15,6) — different codes, different k — have
  IDENTICAL stab orbit counts at w = 6 (1), w = 10 (6), w = 14 (54),
  diverging at 12 (33 vs 42) and 16 (375 vs 478). Reads as universal
  low-weight generator-product combinatorics with code-specific
  corrections; recorded as a curiosity for the F2-style census
  theory, no claim.

## §8 Falsified/corrected claims (session ledger)

- **Charter §1's initial "450 members exactly at d = 2w" — grade
  overclaim, corrected** (this session, charter edited): witness-
  certified ≤ 10 + conjectured floor; tightness exclusive to Z₅×Z₁₅.
- **My DFS Sidon enumerator v1 — REFUTED by its own selftest**
  (order-2 differences); fixed before production.
- **"The S1 sweep protocol = the script defaults" — false**
  (protocol archaeology: structural census + ub-6 screen + MC, not
  ub-12 SAT); two mis-protocol chunk runs were discarded, one killed
  mid-run (their frames were all structurally empty, so no wasted SAT
  beyond minutes).
- Inherited ledger respected; no SAT ran on any floor path.

## §9 Residue / next steps

1. Land the MC falsifier pass (all frames incl. tail) + the Z₅×Z₁₅
   G3 screen; fold both into §2/§5 tables.
2. Harvest the w7 5x17 probe; run the remaining w7 frames (5x18,
   6x15, 9x10, then 7x13/7x14/9x11/10x10) — budget by the 5x17
   timing; any member found ⟹ MC pass at wmax 13 (falsifier gate for
   d < 2w = 14).
3. Q3: the seed cohort is cut (296 members; canonical seed named in
   §4). S2 opens there: BZ/S4-certify the seed's y-rung SF@20 + (M)
   (A30 machinery, ≤ 8.4 min/code precedent), extract its Bezout
   witness, then generalize the exponent pattern across the cohort;
   w = 7 fresh members pend the fat-frame hunt (§6).
4. Q1 residue (S2): re-run P-26 shards + (4,4) strata with outputs
   OUTSIDE the worktree (ops lesson §1); merge
   `worktree-agent-aa0d7676f1aee3b4c` to main first.
5. Merge debt for this branch: registry A40/A41 rows will conflict
   trivially with A40's worktree claim — renumber-at-merge precedent
   applies.
6. **A40 session 1 landed while S1 ran** (branch
   `worktree-agent-aafc19cdb2dd71190`, 5 commits, data force-added ✓):
   headline **d([[432,12]]) = 24 EXACT certificate tier** — the
   conjectured tour-de-gross value at (r,b) = (2,1), first proven
   member beyond the BCGMRY pair, no SAT, ~20 min; b is a BIT (no
   [[576,12]] member exists — the charter's §0 mention of a (2,2)
   member is moot); [[432,12]] is a free **Z₃** cover of gross, so the
   family hangs off gross by mixed decks and the ∀r column's missing
   theory is exactly the odd-deck F1 lane (A40's P3 verdict: RED with
   current technology; the family k-row ∀(r,b) is the GREEN sub-lane).
   Charter cross-effects: Q6's vertical framing gets a concrete
   family-specific instance (the tour-de-gross column is
   gross-centered, not a rung chain); A41's Tier-D odd-axis rule and
   A40's odd-deck demand now share the F1 dependency. Verification +
   merge of that branch = owed (this log treats its numbers as
   reported, not independently re-derived).
7. **A40 S2+S3 landed after this log's close** (same branch, PR #27):
   **Theorem UB(r,1)** (d ≤ 12r for the whole b = 1 column; witness =
   stacked A14-freeze-carrier bands) + **Theorem L1** (plane-regular
   span floor: d ≥ ⌈3r/2⌉ for BOTH columns via gap dichotomy + the
   plane-Koszul regular-sequence kill) ⟹ the b = 1 column is pinned
   **Θ(√n) two-sided (1.5r ≤ d ≤ 12r)**; the per-cell cost lemma died
   in both forms on certified data; L2 (matching constant) is one
   named lemma (x-compression) from runnable. Two A41-relevant
   imports: (i) **L1 stamps span floors on every corpus BB pair,
   including this charter's mirrored-Sidon class and odd-|G| frames**
   — check the Q3 constructive members against it (a Tier-F floor
   mechanism independent of the coincidence-table program, possibly
   competitive with the class theorem at large frames); (ii) the
   per-cell refutation is a standing warning for any Tier-D "cost per
   structure" argument — carried into Q4's known-false awareness.
