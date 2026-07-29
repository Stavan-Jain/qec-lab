# Descent-strengthened SAT: results (2026-07-29)

**Question.** Can the repo's analytic theory speed up the strengthened
SAT solver (`shard_distance`, the coarse "strengthened monolith" whose
whole cost is the UNSAT round at w = d − 1) — and do the same
improvements transfer to Tandem (the MaxCDCL fork)?

**Answer.** Yes on the plain-SAT side — a new quotient/pushforward
layer (`bb_lab/descent_sat.py`) gives a stable **2.8–4.6× further
CPU speedup** on the refutation round (gross: ~50s → ~12s; on top of
the existing 4–6× over the monolith). No on the Tandem side — the
orchestration-level port loses 10–40% (measured, 3 reps), consistent
with the earlier lesson that hard-constraint scaffolding taxes BnB.
Three cheaper hypotheses were falsified along the way and are recorded
below so they are not re-proposed.

All distances in every run are asserted against published/known
values; every witness is independently re-verified; every structural
premise of the layer is asserted numerically per instance at
construction (`compute_descent`'s battery), falsify-first style.

## 1. The layer

For a free axis involution σ ∈ G with quotient p : G → Ḡ (A10–A15
descent direction), three numerically-asserted facts power a sound
decomposition of every shard:

1. **Intertwining** R·H_Z = H̄_Z·S ⟹ the base parity system on
   pushforward bits w = p₊v is implied — free propagation.
2. **Sector dichotomy** (|⟨σ⟩| = 2): p₊v = 0 ⟺ v = p*v̄ — an
   **n/2-variable instance at weight ⌊w/2⌋** (sector a); else the
   cover instance + ⋁w + scaffolding (sector b).
3. **Class transport** μ̄(p₊v) = Λ·μ(v) (the A14 face): sector (b)
   gets k̄ pinned parities; sector (a) reps must lie in Pμ(ker H̄_Z)
   and have even parity — reps failing this are refuted **with no SAT
   call**. The A14 LES (`ker p₁ = im τ₁`, dim k/2) predicted the
   kill fraction; measured: 144/155 (gross y-deck), 150/155 (x-deck),
   32/35 (bb_108), plus B-R refuting in milliseconds.

## 2. Plain-SAT results (CMS, solver-CPU of the refutation round)

Paired same-run comparisons; baseline = `shard_distance` coarse
default. Machine-load variance is real (gross baseline 27–60s across
the day), so ratios are per-run pairs, not cross-run absolutes.

| code | deck σ | descent | baseline | speedup |
|---|---|---:|---:|---:|
| bb_108 [[108,8,10]] | (0,3) | 1.94s | 5.48s | **2.8×** |
| gross [[144,12,12]] | (0,3) | 12.0/12.9/13.6s | 48.5/49.5/59.9s | **4.1/3.8/4.4×** |
| gross | (6,0) | 49.7/70.7/79.3s | same | 1.0/0.7/0.8× |
| bb_72 [[72,12,6]] | both | 0.02s | (trivial) | — |

**Deck choice is an empirical dial, and the theory-preferred deck
lost.** σ=(6,0) is the doubling deck (base = bb_72, the A14/A8
literal-lift direction) yet it is ≈1× or worse; the y-deck σ=(0,3)
(base Z12×Z3, k̄=8) carries the whole win. Benchmark every deck.

Residual structure: all surviving hardness is the B-L task, i.e. the
moving sector — and within it the Λc = 0 reps (pushforward a nonzero
stabilizer), the same "dangerous sector" shape A15 met analytically.

### Falsified increments (do not re-propose without new mechanism)

- **Certified base-coset floors** (`use_floors=True`): on bb_108 the
  floors {6,10} certify 32/35 reps' sector (b) empty at w ≤ 8 — but
  the surviving 3 reps hold all the hardness (1.82s vs 1.63s), and
  the precompute (2.35s, witness-jumping base SAT) exceeds the gain.
  Kept as a flag (off by default); may matter at larger n where the
  precompute amortizes.
- **Per-block parity rows** (K = G augmentation; par(v_L) = par(v_R),
  each block even): consistently **10–20% slower** on bb_90/bb_108
  (3 reps each). Redundant XOR rows tax CMS Gauss more than they
  prune. This was the last lever applicable to odd-order groups
  (bb_90 has no axis deck), so bb_90 currently gets nothing from the
  layer.
- (Prior art, same shape: strengthened WCNF hurts MaxCDCL ~1.7×.)

### Not yet measured

- gross `floors`/`scaffold`/`full` variants (run paused mid-session);
- odd-order decks via the v = (1+σ)u substitution (2n/3 variables —
  the bb_90 route; designed, not implemented);
- deck stacking beyond one level; per-sector LRAT emission (each
  sector instance is proof-capable — the certificate-mode payoff);
- n ≥ 288 (bb_288 has decks (6,0)/(0,6); the 4× would cut the ~20h
  EPYC-equivalent naive wall to ~5h if it transfers).

## 3. Experiment 1: the Tandem transfer (negative)

`maxsat_distance_descent`: d = min(2·opt_a, opt_b) — sector (a) as a
base-sized MaxSAT, sector (b) as the naive WCNF + fiber-pair XOR
definitions + ⋁w (variants add base rows / the Λ-link, which is
expressible without selectors because the naive a_j literals *are*
μ(v)). Tandem fork binary, `-cost-step=2` auto-gated per side.
Medians of 3 reps:

| code | naive | descent-min | descent-rows | descent-full |
|---|---:|---:|---:|---:|
| gross | **1.98s** | 2.43s (0.8×) | 2.75s (0.7×) | 3.27s (0.6×) |
| l168_d14a [[168,6,14]] | **17.76s** | 19.94s (0.9×) | 16.44s (1.0×) | 20.05s (0.9×) |

Sector (a) itself is essentially free and correct (0.04–0.18s;
on gross the invariant sector *attains* d = 12, opt_a = 12 = opt_b),
but excluding the invariant subspace does not help BnB: its bound
machinery was not spending time there, and the aux Tseitin chains
cost more than the exclusion saves. `rows` at n=168/d=14 is the only
≈neutral variant (1.04–1.08×, within noise).

**Conclusion.** The descent structure helps clause-driven *refutation*
(CMS), not branch-and-bound *optimization* — the two lanes stay
complementary: descent-SAT is the certifiable lane (per-sector LRAT),
Tandem-naive the solve lane. The remaining Tandem idea with a real
mechanism is patch-level per-fiber LB rounding (cost(v) ≥ |partial
p₊v| + class floor remainder, the generalization of the queued
block-parity rounding) — untested, scoped as its own session, value
thesis at n ≥ 288.

## 4. Experiment 2: `-fiber-lb` — the per-fiber LB rounding patch

The patch-level transplant scoped in §3, built the same day. The
descent facts enter Tandem's **bound arithmetic** instead of the
clause database: no aux variables, no clauses, an O(n) fiber scan at
the main search's propagate site.

**Solver side** (patch v4, cumulative 348 lines vs pristine): at each
node, classify fibers (decided-odd / decided-even / half-open-true /
oddable) and compute

    nodeLB = min( moving:    t + max(0, needOdd − oddable0),
                  invariant: max(fiberInvFloor, t + oddable0) )

with `needOdd = max(F_b − odd_now, [odd_now = 0])`, `F_b` = min
moving-sector floor over the classes consistent with the currently
assigned a-vars (memoized table lookup), and the invariant branch
active only while no fiber is decided odd — the split that keeps the
rule sound on the plain naive instance, where σ-invariant optima are
feasible. `nodeLB ≥ UB` fires a soft conflict through the
`UBconflictFlag`/`involvedLits` premise hook with the complete premise
set (assigned-false fiber vars + assigned a-vars; assigned-true fiber
vars are already in `falseLits`), so `analyzeSoftConflict`'s learned
clause is valid given the caller's theorem and stays valid as UB
tightens. The premise-free corner (bound holds with no assumptions)
correctly returns the empty clause at level 0 = "no solution below
UB".

**Caller side**: `emit_fiber_certificate` — fiber pairing from the
verified `DescentData`, per-class table `F(c) = base-coset-floor(Λ·c)`
(budgeted ascending base SAT, `base_coset_floors_budgeted`: sound at
any stopping point; the exact witness-jumping engine stalls on heavy
classes of n=144 bases), invariant floor = 2·(exact invariant-sector
minimum). Emission is one-time per (code, deck) and cached.

**Validation** (all with `-cost-step=2`):

- Toy contract test: valid certificate preserves the optimum; a
  deliberately falsified table (floor 5 on a 1-pair instance) yields a
  wrong verdict — the flag has teeth, caller obligation is real.
- Soundness sweep, 8/8 (code, deck) combinations correct: bb_72 ×2
  decks (d=6), bb_108 (d=10), gross ×2 (d=12), l168_d12 (d=12),
  l168_d14a (d=14).
- Firing telemetry: 77 fiber conflicts on gross during descent, ~0
  during the final refutation — at n=144 the floor table (max ≈ base
  scale) sits below the refutation bound, so the rule is marginal
  there, exactly as pre-registered.
- Timing at n ≤ 168: neutral (bb_108 0.32s vs 0.42s the only mild
  win; gross 1.9s vs 1.7s; l168_d14a 14.6s vs 14.7s). The scan
  overhead is invisible; the information is what's marginal at these
  sizes.

**bb_288 certificates** (budgeted emission, ~11 min both decks):
σ=(6,0): distinct-Λ = 64 = 2^(k/2) (the A14 LES dimension count,
visible at n=288), invFloor = **24 > 18 = d** (the whole σ-invariant
subspace is analytically dead), moving floors {6,7,8} (87.5% of
classes at 8). σ=(0,6): same structure, invFloor 24, floors
{6,8,9,10} — the y-deck again the stronger table.

**bb_288 race** (2026-07-29, 1h box, 3 parallel legs, all
`-cost-step=2`): control vs `-fiber-lb` per deck. All legs reached
incumbent 18 (control t≈1146s; fiber60 t≈1741s; fiber06 t≈1939s —
descent differences are opening-phase noise under 3-way contention)
and none resolved the ≤16 refutation inside the box, consistent with
the prior ~85-CPU-min-unresolved runs. The decisive telemetry: fiber
conflicts fired 90 (fiber60) / 144 (fiber06) in total — **essentially
all before the o-18 lock, ≈0 during the ~30-minute refutation grind**.

**Verdict: v1 fiber-lb is sound, live, and does not bite at bb_288's
refutation bound.** The moving-sector floors (max 8–10, = base-coset
scale) sit far below the effective bound (16), so the class-conditional
term can't reach UB except at depths CDCL rarely visits. This is the
third face of the same finding today (plain-SAT floors increment;
A15's base-floors-are-the-weak-half): *base-coset floors are too small
to carry refutation at 2d̄-scale bounds*. The upgrade with a real
mechanism is analytic, not plumbing: certified **seam/window floors
for the moving sector** (A15's DangerousFloorNZ shape — 2d̄-scale, not
d̄-scale) in the per-class table, which exists today only for codes
where the A15 program has been run. fiber-lb is the delivery vehicle
for exactly that: the solver side is built, validated, and waiting for
stronger certificates.

## 5. Experiment 3: fiber-lb v2 — moving-sector cost floors (the "g" table)

The §4 upgrade, built 2026-07-30. v1's table bounds the *odd-fiber
count* (|p₊v| ≥ base-coset floor — provably stuck at base scale, §4);
v2 adds a second table with different semantics: certified floors on
the **moving completion cost itself**,

    F2(λ) = min{ |v| : H_Z v = 0, p₊v ≠ 0, μ̄(p₊v) = λ,
                 (μ(v) ≠ 0 when λ = 0) },

the quantity the seam/window theory bounds at 2d̄ scale. Patch v5
parses an optional `g` line and folds `mLB = max(t + oddArith, G_b)`
into the moving branch (v1 files unchanged; toy teeth-test for the g
table passes both directions). The λ = 0 fiber — the dangerous sector
— needs the nontriviality pin, else stabilizers with nonzero
pushforward dominate the minimum at check-row weight (measured: gross
λ=0 lifts 6 → 10 with the pin; the A15 DangerousFloorNZ analog).

Certification engine: `moving_cost_floor_budgeted` — ascending
conflict-budgeted CMS refutation of the λ-pinned cover instance,
even-step under the verified parity premise, sound at any stopping
point. This is a *budgeted SAT stand-in* for analytic seam floors.

**bb_288 verdict: the stand-in hits the wall it was meant to bypass.**
Per-λ pinned refutation at w ≤ 10 on n = 288 already exhausts 200k
conflicts (per-λ up to 116 s); the g-table came out {6:63, 8:64,
9:3456, 10:512} — ≤ 1 unit over v1, nowhere near the 16–17 needed.
The pre-registered go/no-go (how many of 64 λ reach ≥ 12/14/17): zero
— the 1 h re-race was skipped as already decided. Class pins do not
break the refutation exponential; only *analytic* certification
(window/seam machinery, per-code work as in the [[300,8,16]] program)
can fill the table at 2d̄ scale. The delivery vehicle — g-table
semantics, solver arithmetic, dangerous-fiber pin — is built,
validated, and waiting for it.

**gross (where per-λ refutation IS tractable): the proof phase
collapses.** The g-table fills to {12, 13} for 87 % of classes
(cap 11); with the unrefined λ=0 entry (6) the run is still a tax
(4.14 s vs 3.11 s — min-over-consistent-classes stays 6 until deep).
After the dangerous-fiber refinement (λ=0 certified 12 in 85 s, 800k
budget) the g-table is ≥ 12 on every nonzero class and invFloor = 12
— the certificate now *contains* a solver-word d ≥ 12 proof — and the
measured run is

    fiber-v2 (refined): 0.25 s  vs  step: 3.12 s   (12.5×, 5 reps,
    3 fiber conflicts total, d = 12 re-verified every run)

The search descends to the weight-12 incumbent and the certificate
kills the entire optimality phase in 3 node conflicts: the
`-init-lb`-style solve-time-is-witness-search-only regime, but with
the floor *self-assembled* from 64 per-λ certified refutations + the
invariant floor instead of externally supplied. Emission cost ~3.4 min
(parallel, one-time per (code, deck)) — the economics favor repeated
solves and, more importantly, demonstrate exactly what an *analytic*
seam-floor table would deliver at bb_288 with no SAT filling at all.
The concrete remaining target is now a finite list: certify
F2(λ) ≥ 17 for the 64 Λ-classes of a bb_288 deck analytically.

- `src/bb_lab/descent_sat.py` — the layer (DescentData battery,
  sector solvers, driver; `use_floors`, `scaffold_rest` flags;
  `invariant_floor`, `base_coset_floors[_budgeted]`)
- `src/bb_lab/maxsat_distance.py` — `maxsat_distance_descent`,
  `emit_fiber_certificate`, `fiber_sigma` plumbing
- `third_party/maxcdcl-qeclab.patch` — Tandem v4 (cost-step,
  prime-vars, init-lb, phase-file, **fiber-lb**); `build_maxcdcl.sh`
  now `make clean`s before the fork rebuild (stale-object segfault)
- `scripts/descent_bench.py`, `scripts/descent_maxsat_ab.py`,
  `scripts/block_parity_probe.py`, `scripts/emit_bb288_fiber.py`,
  `scripts/bb288_fiber_race.py` — the harnesses (all assert d where
  known)
- `tests/test_descent_sat.py` — 7 tests incl. the LES-prediction
  check, floors ⟷ base-distance cross-validation, the falsify-first
  non-descending-input loud failure, and the certificate format
  contract
