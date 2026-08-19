# A38 — Distance-descent theory, generalized: the thread charter

**Claimed 2026-08-18** (worktree `distance-descent-theory-plan`). This is a
*thread charter*, not a session log: it scopes the research program that
generalizes and extends the distance-descent theory (the tower slice
calculus of A30/A32/A33/A36 and the conditions map of A35, plus the
fibering/safe-floor engine of A29) with one explicit goal:

> **Develop analytic theory that reduces the computational burden of the
> distance computation problem** — move codes down the cost ladder (§0),
> widen the set of codes the descent machinery reaches at all, and make
> the "burden reduced" claim itself a theorem rather than a benchmark.

Sessions executing this charter log per-session notes in the usual way
(A38 session entries, or new numbers claimed for lines that outgrow it)
and update the docket table (§4) and the ledger (§5). Everything here is
plan-tier unless it cites a banked result; falsify-first discipline (§6.0)
applies to every front.

## §0 The cost ladder, and what "reduce the burden" means

Measured anchors from the program so far, worst-to-best per obligation:

| tier | what it is | anchors |
|---|---|---|
| T0 solver | monolithic SAT/enumeration; exponential, opaque, no reuse | n = 288 UNSAT queries historically intractable (A36 §0); d(BY) = 12 was 3.1 h SAT; A33's H5 had ~16.5 h of banked SAT partials |
| T1 certificate | deterministic censuses + bounded fibers + rungs; exact node-count invariants; minutes | A30 ≤ 8.4 min/code; A32 ~9.4 min; A33 ~105 s (~600×); A36 ~31 s; d(BY) re-derived in 0.04 s |
| T2 kernel | Lean, axiom-clean | d(gross) = 12; pair72 d = 8; cover300 d = 16; mitten [[150,30,10]] d = 10; A36's two-tier bb288 packaging (partial) |
| T3 analytic | theorems that delete enumeration entirely | A12 (R) ⟺ k ⟺ Bezout; A13 tower descent; A16 class small-cycle theorem; A17-P3 wall-value parity theorem; μ_e barrenness |

"Reducing the burden" = two directions, both in scope:

1. **Widen T1**: more codes inside the calculus envelope (new deck types,
   group classes, code shapes) — today the envelope is exactly A35 §2.
2. **Raise T3**: replace T1's enumeration layers (censuses, lift fibers)
   with theorems, so the certificate cost stops growing with the gates.

Discipline carried over verbatim: RED/AMBER/GREEN are *cost* verdicts,
never distance claims; SAT witness weights are never reported as floors;
claim tiers are stated exactly on every result.

## §1 Starting state (v1): what the thread inherits

The calculus v1, per A35's layer table (its §1) and conditions map (§2):

- **Hard requirements**: (C1) a free Z₂ deck — 2 | |G|, Z₂ specifically
  (odd-p decks give τ∘p = norm, no F₂ transfer); (C2) two-block
  group-algebra CSS with **central** deck element. Everything else is a
  simplifier ((R), parity, regime lattice) or a quantified cost gate.
- **Cost gates G1–G5** (A35 §5): bottom-census nodes (exact two-window
  binomial formula; demonstrated ≲ 2e11), fiber caps (demonstrated ≤ 8;
  cap = (W_eff − μ)/2 grows linearly in d_target), sector dispatch
  2^k(base) (fine through k = 12), rates-are-a-bonus, and the τ-branch
  ceiling **d_target ≤ 2·d(mid) per rung**.
- **Universal structure (theorem or 22/22-measured)**: sheet SES + transfer
  LES; rank p* = rank τ* = k(cover)/2 (LES half proven; adjointness half
  owed — §3 F3); carry obstruction = the connecting map δ (liftable ⟺
  ∈ im p* = ker δ); the (R)-bound trio (base-exactness ⟺ σ* = id ⟺ k
  preserved, A12); k(base) ≤ k(cover) corollary (literature check owed);
  pair regimes R1–R4, three of four with closed instances (R4 open,
  methods-only).
- **Closed instances**: [[360,12,24]] (A32), [[288,8,20]] (A33),
  [[360,4,20]] ×2 (A30), [[288,12,18]] (A36, kernel-checked mid-level
  consumed), plus solver-free re-derivations (d(BY), d(GB), d(B72) = 6
  census-complete).
- **Adjacent engines**: A29 `bb_lab.fibering` (portable safe floors,
  field-free weight formula, any odd q); A28 certified BZ census lane +
  ε-trisection; A30 §7 doubling front-end (`bb_lab.doubling_certify`);
  `descend.py` (inverse-direction driver); A17-P3 wall taxonomy.
- **Lean layer** is *ahead* of the Python in generality (GB audit,
  2026-07-28): `XDoubleCoverData` accepts any 2-fiber surjection G →+ H
  (univariate type-checks today); `BBTargetFloor` generalizes the doubling
  assembly to m ≤ 2d; BBDeckTower/Bockstein layers are group-free.

## §2 The walls (the enemy list)

Named, quantified obstructions — each front in §3 attacks at least one.

- **W1 no Z₂ deck**: |G| odd ⟹ the calculus does not apply at all.
  Named: bb_90 (15,3), [[98,6,12]] (7,7), the f2a6 base (5,15). The odd
  lanes (A24-phase-1 designs, λ-certificates, Φ-transfer, A29 fibering)
  exist but are not yet a calculus.
- **W2 census blowup**: [[756,16,≤34]] RED at 1e22.6 bottom-census nodes
  (v₂ = 1, so no deeper tower exists) — the gap is ~11 orders "and needs
  new reductions, not more compute" (A35 §8.4). Note (21,18) is
  odd-part-dominated (odd |G|-part 189, 2-part 4): **the [[756]] wall is
  an odd-part wall**, reachable only via W1 progress or census analytics.
- **W3 fiber-cap growth**: caps grow linearly in d_target; [[720,4]]
  doubling needs cap 16 (RED), [[756]] cap 13 — both beyond the
  demonstrated 8.
- **W4 the τ-branch ceiling**: d_target ≤ 2·d(mid) per rung bounds what
  any tower can certify; codes above their tower's ceiling are out even
  when every census is cheap.
- **W5 non-abelian scarcity**: only central involutions port (A26/A35);
  3 of 8 published mitten codes have one; [[150,30,10]]'s center is C₅ —
  excluded. Lemma 1's central-σ re-derivation is flagged-untested.
- **W6 non-CSS**: entirely out of scope of v1 (the sheet SES and the
  X/Z transpose bookkeeping are CSS-shaped).
- **W7 no pre-census feasibility oracle**: carry-infeasibility rates are
  census-population-depth properties (A35 §7 falsification) — planning
  must assume cap work is real.

## §3 Fronts

### F1 — The odd-deck (Maschke) calculus  [attacks W1, W2]

The deepest generalization. For H ≤ G of odd order, char 2 makes
N = Σ_{h∈H} h an idempotent (N² = |H|·N = N), so F₂[G] splits as
N·F₂[G] ⊕ (1+N)·F₂[G] with N·F₂[G] ≅ F₂[G/H] — descent along odd
quotients is a **direct-summand splitting** (Maschke/CRT), not a
filtration. The two-sector geometry is forced: the invariant sector is
the base code verbatim; the twisted sector is a direct sum of isotypic
components (F_{2^r}-linear codes). None of this is new to the program —
A22's CRT fibering, A29's field-free fiber weights, and the A24-phase-1
λ/Φ designs all live here piecemeal. F1's mission is to systematize them
into a second calculus:

- **F1-Q1 (the coupling law)**: the odd analogue of the slice identity.
  Weight does not add across the idempotent splitting; candidate laws:
  per-fiber restriction weights (A29's formula), window averaging
  (teaching doc §7.6), character-support (BCH-style) floors per isotypic
  block. Deliverable: a per-sector floor assembly `d(C) ≥ min(sector
  floors)` with an exact statement of what couples them.
- **F1-Q2 (the H1 split)**: Maschke ⟹ H₁(C) ≅ H₁(base) ⊕ H₁(twisted),
  canonically — the trisection collapses to a *bisection with a section*.
  Verify the rank bookkeeping on controls; write the regime table (the
  odd analogue of R1–R4 should be trivial — check, don't assume).
- **F1-Q3 (certificate species)**: what replaces censuses + rungs? The
  twisted sector's blocks are smaller codes over larger fields; census
  work should shrink by the CRT factorization. Price the gates the A35
  way before running anything.
- **Controls (falsify-first)**: reproduce banked values on bb_90 (15,3),
  [[98,6,12]] (7,7), f2a6 base (5,15) — all have known distances — before
  any new claim. Flagship: re-price [[756]] after CRT reduction (W2).
- **Kill criteria**: if the coupling law on controls is no stronger than
  "max over fibers" (i.e. no better than A29's existing safe floors) and
  the twisted-sector census does not factor, record the negative with the
  measured gap and reduce F1 to "A29 is already the odd theory".

### F2 — Kill the census: analytic replacement of enumeration  [W2, W3]

- **F2a odd-Fourier census factorization**: censuses are weight-bounded
  walks over ideals of F₂[G]; the odd part of G carries characters, the
  2-part doesn't. Hypothesis: census node counts factor along the odd-CRT
  decomposition (per-component walks, recombined by the coupling law of
  F1-Q1) — the 2-part remains combinatorial, matching "univariate-odd is
  the easy regime" (GB audit). Probe falsify-first on a banked census
  (A32's or A36's, counts must reproduce exactly).
- **F2b the ε-recursion chapter** (A29's named residue): ε = 1+σ makes
  F₂[Z_{2^r}] local; the ε-adic strata of a logical are pulled back from
  lower levels. Target theorem: per-stratum floors at level r from level
  r−1 floors with **no new census** — multi-level assembly from one
  censused level. Direct payoff: fiber caps stop growing with depth
  ([[720]] doubling, cap 16 → recursion instead).
  **S2 RESOLUTION (log §3): the bare-number form is FALSE** (banked m\*
  is element-dependent — the naive SeamCosetFloor refutations, restated
  as data); the true statement is CENSUS-CARRYING — the kernel-shift
  lemma: ker E = Z(base) exactly, so level-r rung/fiber candidates =
  particular lift ⊕ level-(r−1) cycle census in the window
  |x| + cap + ov(v0p), stabilizer-only when the window < d(base).
  Validated on banked towers + executed at [[720]] scale. What remains
  of F2b: shrink the WINDOW-population cost (the honest residual — the
  heavy-shadow windows near n/2 are dense), not the cap.
- **F2c analytic carry floors**: fibers enumerate weight-bounded
  solutions of the linear carry system E v₀ = R b; class-level
  solvability is already free (= δ). The 86–93% measured infeasibility is
  a min-weight statement about carry solutions — find Sidon/expansion-style
  floors for it. **Risk stated honestly**: this is itself a structured
  coset-distance problem and may be self-similar in hardness; time-box
  and record.
- **F2d wall theory**: generalize the parity-dead-band classification
  (A32 Thm 5's β = 0 kills, A17-P3's wall-value parity theorem) toward
  "which (stratum, band) pairs are dead for every code of a class" — the
  analytic pruning that shrinks censuses before they run.
- **F2e witness heredity (find-side burden)**: A36 §3 measured the
  d-attaining logical sitting over the τ-diagonal of the τ-diagonal —
  d-attainment projecting onto d-attainment two levels down. Formulate
  precisely; test on every closed tower (the A33 stratum profile did NOT
  transfer, so heredity must be stated class-wise, not stratum-wise). If
  it holds on a class, upper-bound hunts descend too and the witness
  ladder starts at the bottom.
- **Success metric**: [[756]] or [[720]]-doubling moves RED → AMBER
  (≥ 5 orders shaved with the same soundness story), or a theorem that
  deletes a census/fiber species outright on a named class.

### F3 — The complexity ledger: make the speedup a theorem  [meta]

- **Descent width**: package A35 §5's four gate numbers as a computable
  invariant w(C, tower) shipped by the engine (F5), with the exact
  node-count formula as its core.
- **T3.1 (exponent halving)**: prove the cost statement the measurements
  already show (9.4 min → 105 s → 31 s): for fixed check weights, each
  free Z₂ level halves the exponential term's dimension — certified
  distance in time exp(O(n/2^depth)) + poly(n), with constants from G1–G5.
  Corollary shape: fully-enumerable bottoms (κ = O(log n)) give
  quasi-poly certification for bounded-width towers. State hypotheses
  honestly (W4 ceiling is part of the statement).
- **T3.2 (hardness side + novelty)**: literature pass (A31 method) on
  parameterized complexity of minimum distance — what is known for
  structured/quasi-cyclic codes, so the claim is positioned as new where
  it is new and cited where it isn't. Distance for general CSS is
  NP-hard; whether BB-with-tower makes it tractable-in-depth is exactly
  our theorem's territory.
- **Theory residue inherited from A35 §8.5**: (i) prove the adjointness
  half of rank p* = rank τ* = k/2; (ii) settle k(base) ≤ k(cover)
  (literature or standalone remark).

### F4 — Scope: non-abelian, GB/univariate, non-CSS  [W5, W6]

- **F4a mitten central-σ port**: re-derive Lemma 1 with σ a central
  involution in non-abelian G (A26 flagged the abelian hypothesis
  "likely removable — untested"). Run the calculus on the cheapest of the
  3 published mitten codes with central involutions → first non-abelian
  tower slice. Expect the k-doubling regime ((R) fails structurally) and
  no transpose duality (run both sides).
- **F4b rank-generic tooling + the GB lane**: univariate is uniformly the
  easy regime (F₂[Z_n] is a PIR; Bezout = divisibility; A13's chain-ring
  corollary always true) and the A18 corpus already contains cyclic
  groups in CRT coordinates (Z21/Z24/…/Z78 pairs), unframed. But the
  tooling is bivariate-bound (descend.py, DuckDB scalar ℓ/m, ~250
  scripts). Deliverable: rank-generic tower inventory + descent in the
  engine (F5), then a GB/univariate sweep and a first trivariate probe.
- **F4c non-CSS (stretch, gated on F4a)**: symplectic version of the
  sheet SES for a cyclic non-CSS target (T1 [[13,1,5]] is queue #14).
  Scope only; no session until F4a lands.

### F5 — Consolidation: the engine, the Lean species, the paper  [continuous]

- **Engine (`bb_lab.tower` or descend.py v2)**: (G, A, B) in → tower
  inventory (2-part *and* odd part) → gate pricing (A35 screen as a
  library) → auto-run GREEN closures → tier-labelled certificate bundle.
  Corpus pass over A18's 41 groups + the zoo: count the envelope,
  auto-certify the GREENs, rank the REDs as F2 targets. This produces the
  empirical burden-reduction table.
- **Lean certificate species**: the named gap to unconditional bb288
  (A36 §10.1) — census + rung **data-carriage** certificates (KernelCert
  pivot-certificate precedent from A15), designed once for all towers;
  then the general tower theorem over the already-general
  `XDoubleCoverData`/`BBTargetFloor` layer. Co-designed with QECLean
  sessions; budget discipline per code (≤ 5 min added build).
- **Paper hooks**: this thread is the "mechanism + route" story of the
  paper-1 positioning; F3's theorem is its spine; A31's superlative
  guards apply to every claim.

## §4 Validation docket (initial; sessions update in place)

| target | role | status/verdict today |
|---|---|---|
| banked A32/A33/A36 structure asserts | regression floor for any engine change | **S1**: promoted to `bb_lab.tower.validate_banked()`; gate green 12.9 s, 11/11 towers field-identical to the banked screen incl. RNG stream; all three closures re-run headlessly, outputs identical mod timing (A38_s1_validation_log §1–2) |
| R4-regime methods closure (gross_xx lower pair) | completes the 2×2 regime coverage; values already known | **CLOSED S1** (5.7 s): d(gross) = 12 re-derived certificate-tier; coverage 4/4; seam cosets census-EMPTY at the 2·d(mid) ceiling; descent leakage quantified (84/84 sector-A lifts outside SEAM; ker τ\* vs im p\* 36/54) — log §3 |
| a8_336 [[336,12,12]] full closure | third group shape | **check at merge**: the parallel A37 line reports an end-to-end 19.5 s run (commit `bfa0edb`) — do not duplicate |
| [[720,4]] freeze at W = 18–22 | first re-double tower decision (A14 §13 question) | **S2: DECIDED — d([[720,4]]) ≥ 24 certificate tier; the A14 §13 freeze REFUTED on this tower** (Q1 W=18: d ≥ 20, 334 s, 6,462 rungs; Q2 W=22: 109,011 + 145 rungs ALL PASS, both directions armed — the freeze carrier (flat lifts over the weight-20 seam orbit) is EMPTY; kernel-shift lane carried all 95 deep-cap cells; two completeness gates each phase incl. the independent y-quotient re-derivation; banked A30 2,203-cell census + SeamCosetFloor-20 re-derived; d(L3) = 10 exact free; no upper bound claimed — doubling d = 40? stays RED). Log §2 |
| [[720,4]] doubling (cap 16) | F2b/F2c forcing target | RED; S2's kernel-shift lane removes the n-blind cap wall for LIGHT shadows (window below d(mid)); the heavy-shadow window cost is the honest residual |
| [[756,16,≤34]] (cap 13, 1e22.6) | F1/F2a flagship — the odd-part wall | RED, quantified |
| bb_90 / [[98,6,12]] / f2a6 base | F1 controls (banked distances) | no-deck today (W1) |
| cheapest central-σ mitten code | F4a first non-abelian tower | unscoped |
| A18 coprime-pair (cyclic) corpus rows | F4b GB lane, already in-corpus | unframed |
| corpus envelope census (41 groups + zoo) | F5 burden table | **RUN S1** (now 47 groups): 73 rows — GREEN 41 / AMBER 3 / RED 6 / no-deck 23; every exact-d rep with a deck GREEN; calibration 11/11 == banked (log §4, `data/a38/envelope_census.md`). **S2 correction**: the 8 "cap-bound" open rows were priced against SAMPLED d_ub values (d_lb None; A39: 0/88 tight) — the cheapest closes at **d = 8 exact** (S2 §6, == A39 independently); the genuine W3 instances are certified-large-d codes ([[720]] W22, [[756]]), and the corpus rows await A39's corpus-merge for certified re-pricing |
| F2a factorization probe (a36 census) | F2a falsify-first | **ABSENT S1** (19 s): census reproduced exactly, then the odd-CRT split measured — invariant components ∈ 18ℤ (coset-union law), additivity 0/33,588, all-bounded 0/33,588: the walk does not factor at the census bound; any odd coupling law must be a cancellation law (log §5) |

## §5 Known-false ledger (inherited; do not re-propose)

- Light-concentration outside im p* as a lemma (9/27 counterexample at
  (3,6) — A35 §7). Per-instance only.
- Pre-census rate prediction by shallow sampling (A35 §7).
- The Prop-10 weight-8 stabilizer gap as a general lemma (holds at
  gross/bb72/Y4/BY/GB; fails at Y2/(3,6)) — cite per-instance.
- A33's witness-stratum profile as a transferable pattern (A36 §9).
- Universal descent ("every base has a doubling cover") — 13 hard
  negatives (A10).
- C-safe ⇒ doubling weight-agnostically (A11 E3 counterexample); the
  odd-weight form lives on as the (M)-robustness line = **A37, active in
  parallel** — coordinate, don't duplicate.
- D1∧D2 ⇒ floor ≥ 2w (char-2 Frobenius square, [[98,12,4]]).
- Shifting-as-kill-engine (A28, certified gap cell).
- Fine-shard SAT lanes and the Tandem-port/block-parity/floor SAT hybrids
  (negative-with-data; shard/descent-SAT notes).
- "The deficit wall is exactly 2d−2" as a premise (retracted; the wall
  VALUE is the parity theorem — A17-P3).
- [[576,12]] covers of bb288 (closed-negative, A14 §16) — towers go
  *down* from bb288.
- SAT witness weights as floors (standing).

## §6 Sequencing and discipline

**§6.0 Standing discipline**: every engine change re-passes the banked
asserts before touching anything new; every session logs a
falsified-claims section; claim tiers stated exactly; cost verdicts are
not distance claims.

- **S1 (next session)**: F5 engine skeleton — the A35 screen as a
  library + rank-generic group core; corpus envelope census; headless
  reproduction of the three banked closures; R4 methods closure; F2a
  factorization probe against one banked census. All-validation session;
  no new distance claims expected.
  **EXECUTED 2026-08-18** — all five items landed, no new distance
  claims: `bb_lab.tower` + gate green; 3/3 closures reproduced
  exact-mod-timing; R4 coverage 4/4; burden map 73 rows (frontier is
  cap-bound ⟹ sharpens S3); F2a measured ABSENT with the 18ℤ
  obstruction (feeds S2). Log: `A38_s1_validation_log.md`.
- **S2 (RE-RANKED at the S1 checkpoint, executed 2026-08-18)**: the
  original S2 was F1 odd-controls; S1's evidence re-ranked — the
  measured frontier was cap-shaped (W3), so the F2b/F2c program
  (original S3) was pulled forward. **EXECUTED**: (1) rung engine +
  deep fiber lane promoted into `bb_lab.tower` (`RungCell`,
  `enumerate_lifts_deep`), gate green vs EVERY banked rung battery
  (a32 deep fibers + 2,030 rungs, a33 1,655 + 1,680 per-row, a36
  full direct-close through library parts; one banked artifact found
  and resolved — the a33 planted `found_min_overflow` was a
  violations[:5] truncation, no production impact); (2) the [[720,4]]
  freeze question DECIDED: **d ≥ 24 certificate tier — the A14 §13
  freeze REFUTED** (Q1 d ≥ 20 subsumed; the first same-axis re-double
  that certifiably gains distance; the §13 battery's five freeze
  instances stand as per-instance data, not a law); (3) F2b: number-only
  recursion REFUTED on banked m\*, the census-carrying kernel-shift
  lemma (ker E = Z(base)) validated banked + at scale — the
  ε-recursion chapter now has its exact statement and its executable
  form; (4) F2c: time-boxed NEGATIVE (syndrome-weight floor vacuous
  on 425 ground-truth fibers); (5) stretch: the cheapest cap-bound
  corpus row closed **d = 8 exact** (cross-checks A39) and the
  S1 frontier framing corrected (sampled ubs, not real cap demand).
  Log: `A38_s2_cap_frontier_log.md`.
- **S3 (next)**: F1 session 1 — idempotent-splitting formalization,
  F1-Q1/Q2 on the three odd controls; verdict on whether the odd
  calculus exists beyond A29. Carries S1's F2a constraint (coupling
  laws must be cancellation-shaped) and S2's kernel-shift precedent
  (census-carrying, not number-only, recursions). If the [[720]] Q2
  outcome names new F2 residue, triage it here first.
- **S4**: F3 write-up + literature pass; F4a mitten central-σ port.
- **S5+**: F4b GB sweep; Lean certificate species ongoing with QECLean.
- **Checkpoints**: re-rank after S2 and S3. Thread-level fallback: if F1
  and F2 both stall at controls, the thread reduces to F3 + F5
  (consolidation + theorem + engine) — still a complete, publishable arc.

## §7 Relationship to adjacent active lines

- **A37 (M)-robustness** (parallel, unmerged): owns the
  dangerous-sector-never-binds-alone conjecture and has run a8_336
  end-to-end; F2d coordinates with it and this docket defers a8_336 to it.
- **A35(b) high-FOM search / A36(b) constructed doubling pairs**
  (parallel claims of those numbers, unmerged): *consumers* of the
  calculus (discovery-side). This thread supplies the engine; they supply
  demand and fresh instances.
- **A25 d_circ**: different metric, untouched.
- **Registry collisions**: A35 and A36 are each claimed twice on parallel
  branches; renumber-at-merge owed per the A15/A28-A29/A34 precedent —
  the cross-references above will need patching when that lands (this
  note cites the tower-calculus claimants as A35/A36 throughout).
- **A39** (PR #23, parallel): the pre-registered falsification test of
  the descent machinery — 97 certificate-tier exact distances, 0/88
  prior corpus d_ub tight, corpus-merge staged. This thread treats its
  values as independent cross-checks (S2 §6 did exactly that on the
  Z9xZ6 row: values agree by different machinery); the S1 burden map's
  frontier section should be re-issued after its corpus-merge applies.
  Next free number ≥ A40.
