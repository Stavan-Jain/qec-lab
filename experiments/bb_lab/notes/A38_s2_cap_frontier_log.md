# A38 S2 — the cap frontier: rung/fiber promotion, the [[720,4]] freeze decision, and the ε-recursion chapter

**Session 2026-08-18** (worktree continuation of the A38 thread; S1 =
`A38_s1_validation_log.md`).

**Session identity and the re-rank (charter §6 record).** The charter's
original S2 was F1 session 1 (odd-deck controls). S1's evidence re-ranked
the thread: the corpus envelope census measured the open frontier as
**CAP-bound, not census-bound** — all 8 deck-having open d_ub questions
sit inside the 2e11 bottom-node envelope at fiber caps 9–17 (wall W3),
while W2-style census walls ([[756]]) are atypical in-corpus
(S1 log §4). The highest-value next move is therefore the charter's
F2b/F2c program (original S3) pulled forward, executed here as S2:
promote the rung engine + deep fiber lane (the enabler), execute the
[[720,4]] freeze certification (the first re-double tower decision,
GREEN-priced and unexecuted since A35 §8.2), formulate and validate the
ε-recursion (F2b), and time-box the analytic carry floors (F2c).
**F1 moves to the next session**, carrying S1's F2a finding as its
design constraint: the invariant-sector components of censused vectors
are forced into 18ℤ and weight additivity fails 0/33,588 — any odd
coupling law must be a **cancellation law**, not a weight-splitting law.

Discipline: charter §6.0 throughout — falsify-first hard asserts
against banked numbers before anything new; claim tiers stated exactly;
RED/AMBER/GREEN are cost verdicts, never distance claims; no SAT on any
floor's critical path; witness weights never reported as floors; §5
known-false ledger respected. A39 is claimed by the parallel
validation line (PR #23) — its certificate-tier values are treated as
independent cross-checks where encountered, not duplicated.

Scripts `a38_s2_rung_gate.py`, `a38_c37xx_freeze.py`,
`a38_f2b_recursion.py`, `a38_f2c_carry_floors.py`,
`a38_s2_frontier_probe.py`; data `data/a38/{s2_rung_gate.json,
c37xx/freeze_W{18,22}.json + rungs, f2b_recursion.json,
f2c_carry_floors.json, frontier_z9z6_probe.json}`.

## §1 The promotion: `RungCell` + the deep fiber lane into `bb_lab.tower`

S1's named library residue (its §7.1) landed:

- **`enumerate_lifts_deep`** — the a32 ordered-split MITM (cap ≤ 8),
  verbatim library port.
- **`RungCell`** — the a30/a33 `YRungCell` architecture, rank-generic
  (consumes only the `TowerCode`/`AxisDeck` interface): dangerous rung
  with full 2^k sector dispatch (vacuous / all-trivial / restricted /
  BZ lanes) and the seam feasibility rung. The restricted lane is
  EXTENDED to sizes 7–8 via the a32 ordered-split join; the BZ lane
  builds lazily from `bb_lab.cosetbz` (base.n ≤ 192).
- Rung results now carry `min_overflow` + `ov_hist` computed over ALL
  violations (see the §5 banked-artifact finding).

### §1.1 The gate (`a38_s2_rung_gate.py`) — every banked rung battery

**GATE GREEN, 1,570.5 s (contended box), every banked number EXACT**
(`data/a38/s2_rung_gate.json`):

| battery | banked reference | result |
|---|---|---|
| A1 deep == shallow enumerator | 397 banked a32 sector-C fibers, caps 4/3 | **397/397 bit-level equal** |
| A2 the 28 deep fibers + trichotomy dispatch | `data/a32/deep_fibers.jsonl` (2,371 lifts) + A19 m24 keys + A24 reach16 | **28/28 per-fiber lift hists EXACT; 2,030 rungs ALL PASS through the library RungCell; 341 stab lifts banked-key-checked; 0 new band-16 hits** |
| B1 a33 dangerous | 1,655 rungs, verdict+lane per row == `h2_rungs.jsonl`, M == banked m_req | **1,655/1,655 EXACT** |
| B2 a33 seam | 1,680 rungs, verdict+lane per row == `seam_rungs.jsonl` | **1,680/1,680 EXACT** |
| B3 a33 hexagon + V1 | 256/256 sector-linear; PASS restricted<=6, 240 nontrivial sectors | **EXACT** |
| B4 a33 planted control | see §5 | **deep M=9: 8 distinct ov-8 solutions, all weight 22; BZ M=11: n_viol 540 == banked, min overflow 8, ov≤8 distinct-v0 sets EQUAL across lanes** |
| C a36 direct-close END-TO-END through the library | BZ pass nodes; 33,588/469 stab; 395 seam; rung lane hists; witness | **7,502,279,774 nodes EXACT; censuses + orbit key sets + both lane hists + 469/469 + 395/395 PASS + witness end-to-end + covariance 3+3 — ALL EXACT** |

## §2 The [[720,4]] freeze certification (A35 §8.2 executed)

`a38_c37xx_freeze.py`. Tower (literal polynomial descent at every
level): L0 (60,6) [[720,4,?]] ←x30— L1 (30,6) **[[360,4,20]] = the A30
doubled code 37a70e02:x, d = 20 END-TO-END certificate tier** ←x15— L2
(15,6) **[[180,4,10]], d = 10 same bank** ←y3— L3 (15,3) [[90,4,?]].
Structure asserts: k = 4/4/4/4; pairs (0,1) = R2 (S∩K = 0), (1,2) = R3
(S₁ = ker p₂\*); rungs x-U-(R), x-U-(R), y-T-(R) — the A35 row
reproduced. Architecture: ALL BZ censuses at n ≤ 180; two fiber layers
composed over orbit reps; the kernel-shift lane (§3) for light
shadows; the seam species killed or censused by the A30 certificate.

### §2.1 Q1 (W = 18): **d([[720,4]]) ≥ 20 — certificate tier, 334 s**

- L3-stab census ≤ 18: 219,245 vectors {6: 45, 8: 45, 10: 315,
  12: 1,560, 14: 5,580, 16: 29,115, 18: 182,585}, 4,883 orbit reps,
  μ₃ = 6 (weight-8 stabilizers EXIST at L3 — another Prop-10-gap
  anti-instance, cited per-instance). L3 cosets ≤ 16: 372,624 elements
  over the 15 classes ⟹ **d(L3) = 10 EXACT** (census-complete free
  by-product: (15,3) [[90,4,10]]).
- L2 descent censuses (4,883 stab fibers + 24,846 coset fibers + 90
  τ₃-sources): stab orbit reps ≤ 18: **2,203** {6: 1, 10: 6, 12: 42,
  14: 54, 16: 478, 18: 1,622} — **band-by-band EQUAL to the banked A30
  37a70e02:x dangerous-cell census** (which A30 derived directly at
  n = 360: its 2,203 rung cells, per-lane {bz: 1, r≤4: 6, r≤3: 42,
  r≤2: 54, r≤1: 478, r≤0: 1,622}); S1-coset reps ≤ 18: **EMPTY —
  A30's SeamCosetFloor-20 safe-floor certificate re-derived through
  the composed fibers**; lightest nontrivial L2-cycle seen = 10 =
  d(L2) (tripwire, 111 hits all ≥ 10).
- **Gate (layer 1)**: direct 4-offset BZ census at L2 ≤ 14 == the
  descent slice EXACTLY (103 stab orbits, 0 S1 orbits).
- L1 stab census ≤ 18 (2,203 fibers, all direct): **6,462 orbit reps**
  {6: 1, 10: 6, 12: 88, 14: 54, 16: 1,018, 18: 5,295}, μ₁ = 6;
  L1-logical tripwires 0 (nothing ≤ 18 nontrivial — consistent with
  the A30 d = 20).
- **Gate (layer 2)**: the ≤ 12 L1-stab census re-derived through the
  INDEPENDENT y-quotient (30,3) (different deck, fibers, τ-family) —
  orbit key sets EQUAL (95 orbits).
- **Seam branch: DEAD outright** — seam elements are nontrivial
  L1-logicals, weight ≥ d(L1) = 20 > 18 by the A30 certificate; no
  census needed. b = 0 branch: |v| = 2|u| ≥ 2 d(L1) = 40.
- **6,462 dangerous rungs at target 20: ALL PASS** (124.6 s; lanes
  r≤0: 5,295, r≤1: 1,018, r≤2: 54, r≤3: 88, r≤4: 6, kernel-shift: 1;
  4 lane cross-validations equal; covariance 3/3). The single
  kernel-shift cell (|b| = 6, M = 7) **independently re-derives the
  one bz-lane rung inside A30's banked [[360,4,20]] certificate**
  solver-free and BZ-free.
- Assembly: no nontrivial X-logical ≤ 18 ⟹ **d ≥ 20** (consuming
  d(L1) = 20, d(L2) = 10 at A30 certificate tier); Z side by
  transpose duality. The freeze floor half of A14 §13's question is
  decided: the re-double did not LOSE distance.

### §2.2 Q2 (W = 22): **d([[720,4]]) ≥ 24 — the A14 §13 freeze is REFUTED at certificate tier**

Executed in two phases (the harness kills background tasks at ~1 h
wall — see §9.6): `a38_c37xx_freeze.py 22 --census-only` (2,860 s) +
`a38_c37xx_rungs.py 22` (1,388 s), the checkpoint reloaded with EVERY
vector re-verified from scratch (cycle, stabilizer/seam-class
membership, weight window, the A30 tripwires).

- Census phase (same stack as Q1, W = 22, WC = 18): L3-stab 282,574
  orbit reps; L2 descent: stab 33,689 orbit reps ≤ 22 (the ≤ 18 slice
  again == the banked A30 2,203-cell census band-by-band), S1-cosets
  ≤ 22: {20: 31, 22: 852} — populated exactly where A30's
  SeamCosetFloor-20 said they start; L2 gate ≤ 14 EXACT again; L1-stab
  **109,011 orbit reps** ≤ 22 {…, 20: 12,183, 22: 90,366}; the
  y-quotient L1 gate EXACT again; L1-seam **145 orbit reps
  {20: 1, 22: 144}** (all ≥ 20, tripwired); L1-logical tripwires 18,
  all ≥ 20.
- Rung phase, both directions armed (the complete ≤ 22 sweep IS the
  witness ladder — any weight-20/22 nontrivial logical must appear as
  a rung violation over its own shadow): **dangerous 109,011/109,011
  PASS** at M = (24−|b|)/2 (1,295.8 s; lanes r≤0: 90,366,
  r≤1: 12,183, r≤2: 5,295, r≤3: 1,018, r≤4: 54, **kernel-shift: 95**
  — every deep-cap cell (|b| ∈ {6, 10, 12}, M−1 ∈ 5..8) went through
  the F2b lane instead of a C(354,4)-scale MITM; 6 lane
  cross-validations equal); **seam 145/145 PASS** at M ∈ {2, 1} —
  in particular the single weight-20 seam orbit (the freeze
  mechanism's designated carrier: A14 §13's "undoubled-direction
  logical, lifted" would be a flat lift over exactly such an element)
  admits NO lift of overflow ≤ 1: **the freeze carrier is empty**.
  Covariance 3/3.
- Assembly: b = 0 branch dead (≥ 2 d(L1) = 40); dangerous + seam
  closed + G-transport ⟹ NO nontrivial X-logical of weight ≤ 22:
  **d([[720,4]]) ≥ 24, certificate tier** (consuming d(L1) = 20 and
  d(L2) = 10 at A30 certificate tier); Z side by transpose duality.

**Stated exactly**: the question that closed is "d ≥ 24?" — YES, at
deterministic certificate tier (Q1's "d = 20 if frozen" floor is
subsumed). The A14 §13 freeze prediction (same-axis re-doubles pin at
d(rung-1) = 20) is **refuted for this tower**: the first re-double
that certifiably GAINS distance (≥ 24 > 20). No upper bound on
d([[720,4]]) is claimed — no witness was found ≤ 22, no SAT was run,
and the τ-branch ceiling (2 d(L1) = 40) bounds what this tower's
recursion can certify, not d itself; the full doubling question
(d = 40?, W = 38) remains RED per A35. The A14 §13 battery's five
freeze instances stand — the pattern is now known to be
per-instance, not a law of same-axis re-doubles.

## §3 F2b — the ε-recursion chapter

**The formulation** (per Z₂-deck rung p: C → B, Lemma-1 setting,
parity; `a38_f2b_recursion.py`, all identities asserted numerically on
the a36 tower):

- ε := τ∘p is multiplication by 1+σ on C; ker ε = im ε = im τ (C is a
  free F₂[Z₂]-module). The slice identity in ε-form:
  **|v| = |εv|/2 + 2m** (since |εv| = 2|p v|). The ε-strata of a
  nontrivial logical are exactly the assembly's branches: v ∈ im τ
  (εv = 0 ⟺ b = 0) is the τ-branch — |v| = 2|u| ≥ 2 d(B), the G5
  ceiling as the DEGENERATE case of the recursion (level-(r−1) floor
  consumed as a bare number); v ∉ im τ gives |v| = |b| + 2 m*(b) with
  m*(x) = min overflow over nontrivial lifts of x — the rung content,
  stratified by [b] (stab/seam).
- **(N), the number-only recursion — m*(x) ≥ (t−|x|)/2 from
  level-(r−1) floors alone: REFUTED by banked ground truth** (it is
  the naive SeamCosetFloor). The banked m* table: a36 w12 seam element
  (class 0x40): m* = 3 = the deficit bound (the d-attaining cell
  SATURATES it); a36 w14 stab stratum: m* ≥ 3 > bound 2; a33 w18 seam:
  m* = 2 (lightest 22); **a33 w14 seam: m* = 4 > bound 3** — m*
  depends on the element, not on |x| and d(B) alone.
- **(C), the census-carrying form — the kernel-shift lemma (new,
  validated)**: the carry system's solution set is v0p + ker E, and
  **ker E = Z(B) exactly** (τ is injective on 0-chains, so
  E z = τ₀ HZ_B z = 0 ⟺ z a B-cycle). Hence every candidate with
  overflow ≤ cap is v0p ⊕ z with |z| ≤ |x| + cap + ov(v0p): **the
  level-r rung consumes the level-(r−1) cycle census in a weight
  window — no new enumeration species at level r.** Choosing v0p = the
  row-decomposition lift (x = Σ_{i∈I} rows ⟹ v0p = the same sum of
  cover rows) gives ov(v0p) ≤ (6|I| − |x|)/2, small for light x; and
  when B := |x| + cap + ov(v0p) < d(B), the window holds STABILIZERS
  ONLY — the rung is decided by the level-(r−1) stabilizer census.
- **What the recursion needs per rung, exactly** (the charter's asked
  deliverable): (i) a low-overflow particular lift (row-decomposition
  for stab cells; seam cells have none — solve_E's overflow is
  uncontrolled, so (C) reaches seam cells only through generic
  windows); (ii) the level-(r−1) cycle census complete to B. **Where
  it fails to shrink**: B ≥ d(B) forces the level-(r−1) logical-coset
  censuses into the window (their cost is real — the [[720]] execution
  builds them to ≤ WC for exactly this reason), and the failure is
  forced by the (N)-refutation above: no bare-number form can close
  the light cells, because m* is element-dependent.

**Validation record**: V1 identities 40+60 asserts (a36 tower); V2 the
banked m* table (above); V3 the kernel-shift rung vs the banked a36
dangerous battery — coverage 1/469 cells with a stab-only window
(B < 12; the covered cell IS the expensive cap-5 lane cell:
complementary coverage measured), verdict == banked PASS; at [[720]]
scale: the W = 18 run exercised the lane on the w6 cell (+ 4
cross-validated cells), and the W = 22 run is its production test —
see §2. The A35 §5 cost-gate refinement this forces is in §7
(the cap gate is n-blind).

## §4 F2c — analytic carry floors (time-boxed): NEGATIVE, with the measured gap

`a38_f2c_carry_floors.py` (7.9 s), the structure-free candidate — the
greedy syndrome-weight floor m₂ ≥ min k with the k largest
reduced-column weights summing to ≥ wt(rhs_res); certification =
bound > cap, no enumeration:

- **Battery**: all 425 GB stab-census orbit fibers of the a32 tower
  (caps 3–8 — contains the banked 397-fiber caps-4/3 layer AND the 28
  deep fibers), ground truth recomputed by enumeration per fiber,
  soundness asserted against every enumerated minimum (0 violations).
- **Measured verdict: the bound is vacuous.** Its value is 0 on
  295/425 and 1 on 130/425 fibers — never ≥ 2 — against caps 3–8;
  0 empties exist in the support-carrying ground truth, and 0 would
  have been certified. On the banked heavy bands (empties
  319/1,733 @ cap 2, 5,635/10,602 @ 1, 55,555/64,619 @ 0 —
  reproduced), a cap-0 "certificate" (rhs_res ≠ 0) is identical to
  the cap-0 enumeration itself, and the bound-value distribution
  makes cap-≥ 2 certification implausible (needs bound ≥ 3).
- **The honest structural statement** (recorded, not pursued): by the
  kernel-shift lemma the fiber minimum m₂ is the Hamming distance from
  v0p (off supp b) to the PUNCTURED base cycle code — a structured
  coset-distance problem, i.e. the charter's named self-similarity
  risk verbatim. Any nonvacuous analytic carry floor must therefore be
  a punctured-code distance bound, not a syndrome-weight bound.
- **Stop per the charter's time-box**: the floor is no stronger than
  enumeration on all banked data; negative recorded with the numbers
  (`data/a38/f2c_carry_floors.json`). The production answer to W3 in
  this session is the kernel-shift lane (§3), not an analytic floor.

## §5 Mismatch found and resolved: the banked `found_min_overflow` was a truncation artifact

The gate's first run STOPPED (per §6.0: do not proceed past a mismatch)
on battery B4: the banked `data/a33/rung_validation.json` planted
control reports `found_min_overflow: 10` for the w6-shadow cell at
M = 11 (BZ lane, 500.8 s banked), but the new restricted<=8 lane found
**8 distinct overflow-8 solutions** (weight-22 nontrivial logicals —
no tension with d = 20). Diagnosis, three independent ways:

1. the overflow-8 element re-verified end-to-end (E-system, non-stab
   chain, slice identity; |v0| = 11);
2. a one-offset BZ walk over exactly its sector coset (nodes exact)
   FINDS it — the walk machinery is complete, and the sector IS in the
   nontrivial-sector offset list;
3. the full M = 11 BZ lane re-run reproduces the banked violation
   count n_viol = 540 EXACTLY and reports min overflow 8 once the
   minimum is taken over ALL violations.

Root cause: the frozen `a33_rung_cell.py` control computed its
"found_min_overflow" over the TRUNCATED `violations[:5]` list
(arbitrary append order), not over all 540. **No production verdict
anywhere consumed a truncated-min datum**: the only banked bz-lane
verdicts are (i) this V2 control (whose purpose — find the planted
violation, no sub-20 weight — is intact) and (ii) the A30 rung files
(`rungs_37a70e02_x.json`, `rungs_5e50a976_{x,y}.json`), each containing
exactly one bz-lane rung with verdict PASS (no violation list at all —
completeness of the walk is what a PASS consumes, and item 2 above is a
positive completeness check). The A30 37a70e02:x bz rung (w6 shadow,
M = 7) is additionally re-derived deterministically by this session's
[[720]] stack (§2), which runs the same cell through the restricted<=6
MITM lane — solver-free and BZ-free.

The library `RungCell` now returns `min_overflow`/`ov_hist` over the
full violation set, and the gate asserts the deep-lane and BZ-lane
ov ≤ 8 **distinct-solution sets are equal** on this cell.

## §6 The cap-frontier probe (stretch): the cheapest open row closes at d = 8 — and the burden map's frontier framing corrects

`a38_s2_frontier_probe.py` (8.6 s) on Z9xZ6 [[108,4, d_ub 26]]
(envelope AMBER, cap 9, the cheapest of S1's 8 cap-bound rows):

- The (9,3) bottom census (9.7e5 nodes) measures **d(L1) = 4** and the
  τ-surviving coset minimum **d_τ(L1) = 4** ⟹ τ(u_min) is an explicit
  nontrivial [[108,4]] logical of weight 8, re-verified end-to-end —
  **certified upper bound d ≤ 8** (upper bound only, per discipline).
- The floor half, DIRECT (n = 108 fits the C kernel whole): the
  complete 16-offset L0 coset census ≤ 6 (2.5e4 nodes, exact asserts)
  is EMPTY on every nonzero class ⟹ d ≥ 8.
- **d([[108,4]] Z9xZ6 00f8eb7a) = 8 EXACT, certificate tier, both
  halves this session's own.** Provenance: the parallel A39 line
  (PR #23) had independently closed this row (phase-2
  CERTIFIED_FLOOR 8, closed, W = 6) — the values agree; this probe is
  an independent cross-check by different machinery (τ-witness +
  direct window vs its engine), not a duplication.
- **The interpretive correction to S1's burden map**: the corpus
  d_ub values are SAMPLED upper bounds (this row: d_lb = None,
  d_method = None; A39 measured 0/88 prior d_ub tight across its
  cohort). S1's "the open frontier is CAP-bound at caps 9–17" priced
  the caps against those loose ubs; under certified probing the
  cheapest row's true question sat at d = 8 (cap −1 ⟹ trivially
  GREEN). The REAL cap wall is defined by codes whose TRUE distance
  is large — the [[720]] W = 22 push (§2.2) and [[756]] remain its
  genuine instances; the corpus's 8 "cap-bound" rows need certified
  re-pricing (A39's corpus-merge, once applied, does most of this).

## §7 Verification map

| claim | check |
|---|---|
| library RungCell/deep lane == every banked battery | §1.1 gate: A1 397/397 bit-level; A2 28/28 + 2,030 rungs; B1 1,655/1,655; B2 1,680/1,680; B3; B4 both lanes; C full a36 EXACT |
| the banked-artifact diagnosis | §5: three independent reproductions (deep lane, one-coset walk, full BZ M=11 with n_viol 540 exact + distinct-v0 set equality) |
| [[720]] tower structure | constructor asserts + R2/R3 pair identities + A35 row asserts in-run |
| [[720]] census completeness, layer 1 | direct 4-offset L2 BZ ≤ 14 == descent slice (orbit key sets) |
| [[720]] census completeness, layer 2 | the y-quotient (30,3) independent re-derivation of the ≤ 12 L1 census |
| [[720]] vs banked A30 | 2,203 dangerous cells band-by-band; S1-cosets ≤ 18 empty = SeamCosetFloor-20; the A30 bz-lane rung re-derived solver-free |
| [[720]] rung soundness | in-line E-system/nontriviality/slice asserts per candidate; kernel-shift windows bounded below d(L1); lane cross-validations; covariance |
| d(L1)/d(L2) inputs | A30 certificate tier (banked decide_37a70e02.json), tripwired in-run (0 hits below) |
| F2b identities + m\* table + coverage | `f2b_recursion.json`: V1 100/100, V2 banked artifacts, V3 1/469 coverage verdict == banked |
| F2c bound soundness | asserted ≤ every enumerated minimum on 425 fibers; banked heavy-band rates reproduced |
| frontier probe | node-exact censuses both levels; witness re-verified end-to-end; A39 value agreement |

## §8 Falsified-claims ledger (session-internal)

- **The banked `found_min_overflow: 10` datum (a33 rung_validation.json)
  — REFUTED as reported**; truncation artifact, see §5. The underlying
  walk completeness stands (n_viol 540 reproduced exactly).
- **The A35 cap gate is n-blind — measured**: "cap ≤ 8 demonstrated"
  was demonstrated at n = 90 (C(84,4) ≈ 2e6 half-subsets); the same
  cap-8 ordered-split at n = 360 costs C(354,4) ≈ 6.2e8 (~300×) and at
  n = 180 C(174,4) ≈ 3.7e7 (~19×). G2 should be priced as
  C(n − μ, ⌈cap/2⌉) per fiber, not as a bare cap value. The
  kernel-shift lane (§3) is the session's answer for light shadows.
- **The number-only ε-recursion (N)** — refuted by the banked m\*
  table (§3); only the census-carrying form survives.
- **"The corpus open frontier is cap-bound at caps 9–17" (S1 §4
  headline) — corrected as-interpreted** (§6): the caps were priced
  against sampled d_ub values that certified probing collapses (the
  cheapest row: 26 → 8 exact; A39: 0/88 prior d_ub tight). W3 remains
  real on codes with certified-large d; the corpus rows were not
  evidence for it.
- **F2c syndrome-weight carry floor — vacuous as measured** (§4).
- (respected, inherited): no SAT anywhere; witness weights never
  reported as floors; RED/AMBER/GREEN cost verdicts only; Prop-10 gap
  cited per-instance (new anti-instance at c37xx L3).

## §9 Residue / next steps

1. **Kernel-shift promotion**: `KernelShift` + `row_lift_v0` +
   `min_stab_decompose` live in `a38_c37xx_freeze.py` (imported by the
   F2b script); promote into `bb_lab.tower` next session with the
   c37xx + banked-A36 batteries as the gate (the promotion precedent
   is now twice-established).
2. **The window-population residual** (the honest F2b remainder): the
   kernel-shift converts the cap wall into a window-population wall;
   dense windows (weights near n/2 at small n, heavy shadows at large
   n) are the new named cost. Candidate: orbit-bucketed windows +
   on-the-fly weight filters (constant-factor), or a genuine
   punctured-code distance floor (= the F2c self-similar object —
   blocked, see §4).
3. **Corpus re-pricing**: apply/consume A39's corpus-merge, then
   re-run the envelope census against certified d values — the
   burden map's frontier section should be re-issued (S1 §4's
   headline corrected per §6/§8 here).
4. **F1 (next session)**: unchanged demand (22/57 rows odd-locked);
   carries S1's cancellation-law constraint and this session's
   census-carrying precedent.
5. **Lean**: the freeze certificate's species are the same
   census+rung data-carriage shape as A36 §10.1's named gap — the
   [[720]] Q1 stack adds the composed-fiber and kernel-shift species
   to the design space of the certificate-species work (charter F5).
6. Wall-time notes: the gate's B4ii BZ cell cost 1,473 s against a
   banked 500.8 s — three concurrent suites on one box; node counts,
   never seconds, are the invariants (S1's rule re-confirmed). NEW
   ops rule: **the harness kills background tasks at ~1 h wall** (the
   first Q2 one-shot died at 80,000/109,011 rungs); long stacks must
   checkpoint — the census/rungs two-phase split (`--census-only` +
   `a38_c37xx_rungs.py`, loader re-verifies every vector) is the
   pattern.
7. The Q2 lane histogram is the kernel-shift lane's production
   validation (95 deep-cap cells, 6 cross-validated); its promotion
   (item 1 of this list) should carry that battery.
8. The freeze refutation reopens the value side of the c37xx tower:
   d ∈ {24, …} with the doubling target 40 still RED (cap 16 at
   W = 38) — but the W = 24/26 pushes (caps 9/10 at n = 90 bottoms)
   are now the cheapest way to keep walking d upward; price before
   running (window populations, not caps, are the real cost — §8).

## §10 Session commits (worktree branch `worktree-agent-a785190aa3e0b8975`, in order)

- `7f79f5e` Merge S1 (claude/a38-s1-validation @ 1428632) into the S2
  worktree (no file overlap with main's PRs #20–22)
- `c23b908` a38: promote the rung engine + deep fiber lane into
  bb_lab.tower (S2)
- `093b4a7` a38: [[720,4]] freeze certification, Q1 — d >= 20 at
  certificate tier (334 s)
- `ce08850` a38: F2b epsilon-recursion chapter + F2c carry-floor probe
  (S2)
- `9645e66` a38: cap-frontier probe — the cheapest open row closes at
  d = 8 exact (8.6 s)
- (close-out commit) a38: [[720,4]] Q2 — d >= 24, freeze refuted; S2
  session log + charter updates
