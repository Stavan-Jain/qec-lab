# A30 — coset-BZ: the d = 10 doubling cells decided (docket CLOSED)

> **STATUS.** Session 1 (2026-08-06/07, branch
> `claude/light-stabilizer-bounds-d25964`, same session as the A28/A29
> merge). Goal (user brief): using A28 (BZ certified enumeration) and A29
> (general fibering engine) together, find d = 10 base codes where the
> doubling argument provably works — time-boxed to **15 minutes of
> compute per code**. Outcome: **both Z₁₅×Z₆ [[180,4,10]] docket codes
> certify on every axis** (three cells, ≤ 8.4 min/code) — the docket's
> only remaining UNKNOWNs, so the A17 safe-floor docket has **zero
> UNKNOWN cells left** — plus a 0.6 s independent cross-check of the
> CMS-certified Z₂₁×Z₃ [[126,8,8]] cell (and an A29 §5.3 erratum: its
> "4 remaining UNKNOWNs" read a stale docket row; the true count was 3).
> Certificate tier throughout (deterministic C enumeration + counting
> invariants; no SAT, no Lean).

Companions: `A28_light_classification_theory.md` (the BZ census lane this
extends), `A29_general_fibering.md` (seam offsets, kernel orbit reps, and
the engine these cells were beyond), `A17_d7plus_doubling_hunt_plan.md`
(the docket; §6.5 the [[300,8,16]] hypothesis-set precedent),
`A14_safe_floor_criterion_plan.md` §2 (Prop A14.1, consumed for the
one-rep-per-orbit reduction), `A27_safe_floor_generality.md` §2 (the I4
cost anatomy this lane attacks).

---

## 1. The idea: safe floors are coset weights, and coset-BZ decides them

The safe-floor obligation of the doubling template asks: for every safe
class, is the coset minimum ≥ 2d? A29's state-sweep engine decides this
when a linked fiber exists or the unpaired DFS stays small — and its §5.3
diagnosis left the three Z₁₅×Z₆ floor-20 cells beyond reach (no linked
fiber; unpaired tree ≈ C(36, ≤18) ≈ 3·10¹⁰ nodes, "pruning only bites at
the budget boundary").

The A30 move: don't fiber — **enumerate**. The safe class coset is
`t0 + C_AB` with `C_AB = {(A⋆f, B⋆f)}` (dim κ = |G| − dim K) and
`t0 = (seamC_u, seamC_v)` per Prop A14.1(3). The A28 two-window argument
extends verbatim from linear codes to cosets:

- Take disjoint information sets I₁, I₂ for `C_AB` with systematic
  generators G₁, G₂.
- Any coset element c is the *unique* coset element matching its own
  window restriction: `c = c_∅ⱼ ⊕ XOR_{i ∈ supp(c|Iⱼ)} Gⱼ[i]`, where
  `c_∅ⱼ` is the coset element vanishing on Iⱼ. The enumeration is the
  same Gray-style subset walk as A28's census, merely **seeded at c_∅ⱼ
  instead of 0** — a one-line change to the C kernel, plus per-node
  weight checks against several coset base words at once (one shared
  XOR-walk serves every axis and class of a code simultaneously).
- **Asymmetric complete pair** (the production schedule): if
  |c| ≤ W then |c|_{I₁} ≤ r₁ or |c|_{I₂} ≤ r₂ whenever
  r₁ + r₂ + 2 > W. For W = 18: (r₁, r₂) = (9, 8) — the second window
  runs 9× cheaper than the naive (9, 9). See §4 for why this mattered.

The certificate is the A28 species: windows disjoint (listed), systematic
identity blocks, **exact node counts** Σ_{s≤r} C(κ, s) (asserted, thread
counts summed), the coset parity (|A|, |B| odd ⟹ coset weights ≡ |t0|
mod 2, so target 20 decides at W = 18), and the ∅-pattern base words
checked in Python. Completeness is a counting invariant — no UNSAT
anywhere. Refutations would return verified witnesses (weight, coset
membership; at a complete pair the found minimum is exact); these runs
produced none.

The same lane certifies **base floors**: d ≥ d₀ on each CSS side is
"every logical class coset has min ≥ d₀" — 2^k − 1 cosets of the rank-κ
side codes (X: ker D₁/im D₂; Z: ker D₂ᵀ/im D₁ᵀ, D₂ = [[MA],[MB]],
D₁ = [MB, MA]) — at r = ⌊(d₀−1)/2⌋ these run in seconds, and a run at
d₀ + 1 *refutes* with a weight-d₀ logical, pinning exactness and
supplying the lifted 2d cover witness.

## 2. Validation record (10/10 + smoke)

`scripts/a30_coset_bz.py validate` — every recorded coset ground truth in
the program, smallest to largest (`data/a30/validate.json`):

| target | ground truth | coset-BZ | agree |
|---|---|---|---|
| pair72 x@8 / x@10 | A14 exact minima 8/8/8 | certify@8; refute@10 min 8 | ✓ |
| f2a6 y@16 / y@18 | A23/A29: floor 16, d_safe = 16 | certify@16; refute@18 min 16 | ✓ |
| ac46bbea y@16 / y@18 | A29 NEW: d_safe = 16 exact | certify@16; refute@18 min 16 | ✓ |
| 38d3c884 x@16, y@18 | A29 NEW: 16 exact / ≥ 18 | certify both | ✓ |
| a8base x@12 / x@14 | A29 §5.2: minima 12/15/15 | certify@12 3/3; @14: rep0 refute min 12, rep1/2 certify | ✓ |
| f2a6 base d≥8 / d≥10 (smoke) | d = 8 exact | certify@8 both sides; refute@10 at 8 | ✓ |

Both A29-new decisions are hereby **reproduced by an independent
certificate species** (state-sweep certificates vs enumeration counting
invariants — no shared trust base beyond the seam-offset construction).
The smoke test earned its keep: at k = 8 the 255-coset base-floor pass
exposed a C buffer overflow (arrays sized for 16 offsets) that the k = 4
targets could never trigger; fixed and re-validated before production.

## 3. NEW DECISIONS — the three remaining docket UNKNOWNs certified
## (+ one CMS cross-check)

`scripts/a30_coset_bz.py decide` (+ the §5 bonus driver);
`data/a30/decide_*.json`:

| cell | code | safe floor | verdict | certificate | compute |
|---|---|---|---|---|---|
| `37a70e02:x` | [[180,4,10]] Z₁₅×Z₆ | 20 | **CERTIFIED** (NEW) | κ=88, W=18, (9,8) pair, 785·10⁹ nodes, 0 hits | 378 s |
| `5e50a976:x` | [[180,4,10]] Z₁₅×Z₆ | 20 | **CERTIFIED** (NEW) | shared walk, 2 offset lanes, 714·10⁹ nodes, 0 hits | 499 s (both axes) |
| `5e50a976:y` | same code, y-axis | 20 | **CERTIFIED** (NEW) | (same pass) | — |
| `16884e06:y` | [[126,8,8]] Z₂₁×Z₃ | 16 | CERTIFIED (cross-check) | κ=59, W=14, (7,6), 3 orbit reps, 0 hits | **0.6 s** |

Single kernel orbit rep per Z₁₅×Z₆ cell (the G-action orbits all three
nonzero classes; Prop A14.1 makes one coset per axis sufficient); seam
parities all even. Base floors d ≥ 10 (resp. ≥ 8) certified both CSS
sides inside the same runs; exactness pinned by weight-10 (resp. -8)
logicals found at target d+1 in ~1 s (`data/a30/exactness_witnesses.json`).

**Docket accounting + A29 §5.3 erratum.** The docket jsonl carries a
stale first UNKNOWN row for `16884e06:y` beside its later SF-CERTIFIED
row (CMS reps 3× UNSAT; A17 §6.1's "18/21, every floor-16 cell
certified" prose is correct) — A29 §5.3's "4 remaining UNKNOWNs" counted
that stale row. The true remaining set was exactly the three Z₁₅×Z₆
floor-20 cells, all decided above; the Z₂₁×Z₃ run is an independent
second-species confirmation of the CMS certificate (three UNSAT rep
queries redone as one 0.6 s counting certificate). **The A17 docket now has zero UNKNOWN
safe-floor cells**: 14 docket-era CMS certifications (incl. 16884e06),
f2a6 by CMS + A29 + A30, three Z₅×Z₁₅ cells by A29, and the three
Z₁₅×Z₆ cells above by A30.

## 4. The 15-minute box: an honest engineering record

The first production attempt ran the naive symmetric (9,9) schedule and
**hit the budget on both codes** — measured throughput at depth 9 is
1.48·10⁹ nodes/s (8 threads), ~40% under the shallower A28 calibration,
so two full windows needed ~15.2 min of enumeration alone. The budget
enforcement cut them off as designed (base floors banked; the aborted
attempt also surfaced two fixes: partial reports now survive timeouts,
and deadlines use the monotonic clock after a laptop sleep inflated one
wall-clock reading to 26,247 s). The asymmetric (9,8) refinement — a
1.8× cut available by arithmetic alone — brought the runs to 6.3 min and
8.3 min per code at 10 threads. Ladder steps never needed: shard-and-sum
distribution (the certificate composes across machines) and the
ε-recursion/trisection theory route remain documented options for larger
budgets-vs-sizes.

## 5. Consequence: two d = 10 codes where the doubling provably works

Per the [[300,8,16]] precedent (A17 §6.5), the doubling claim
`d([[2n, k, 2d]]) = 2d` rests on exactly
`{LogicalFloor d, LightClassification@2d−2, SeamCosetFloor 2d}` + (R) +
the lifted witness. Scoreboard for the three Z₁₅×Z₆ cells
([[360,4,20]] targets):

| input | 37a70e02:x | 5e50a976:x | 5e50a976:y | tier |
|---|---|---|---|---|
| (R) / k-preservation | ✓ | ✓ | ✓ | A29 run (k = k̃ = 4) + A12 theorem |
| LogicalFloor 10 (both sides) | ✓ | ✓ | ✓ | A30 certificate |
| d_base = 10 exact (witness) | ✓ | ✓ | ✓ | A30 refutation@11 |
| LightClassification@18 | ✓ | ✓ | ✓ | A28 BZ census (2,203 / 2,371 classes) |
| SeamCosetFloor 20 | ✓ | ✓ | ✓ | **A30 certificate (NEW)** |
| cover witness ≤ 20 | ✓ | ✓ | ✓ | lift of the weight-10 base logical |

**These are the program's first d = 10 → 20 doublings**, and the first
codes to clear the deficit wall (A17 P3: every previously examined
orbit ceiling stalled at 2d − 2; bb_90/bb_108 remain refuted-with-
certificates, so the wall was a fact about *those* codes, not the
mechanism). The [[126,8,8]] → [[252,8,16]] bonus cell completes its
template the same way.

**Claim tier, stated precisely**: certified computational data — 
deterministic enumeration with composable counting invariants, the same
tier as A29's engine runs and strictly more auditable than SAT-without-
DRAT. NOT yet kernel-checked: the Lean packaging would follow the A15
pattern ({logicalFloor, lightClassification, seamCosetFloor} instance
layers + re-instantiating the `dangerousFloorNZ_of_lightClassification`
dangerous-sector theorem at (n, d) = (180, 10); the A22/A23-shape
certificate mass is the known route). The dangerous-sector Lean theorem
exists today only for the [[300,8,16]] instance — its architecture is
portable, its Lean artifact is not.

## 5.5. Session 2 (2026-08-07): the rung pass — d = 20 END-TO-END

The §5 scoreboard's one pending arrow (classification ⟹ dangerous floor)
is now discharged: `scripts/a30_rung_pass.py` certifies **(M) at 20** for
every census class of all three cells. Per the teaching doc's dangerous
sector: a nontrivial dangerous logical has sheets (v₀, v₀ + b), shadow
b = p(v) a stabilizer, wt(v) = wt(b) + 2·overflow (slice identity), and
(M) needs overflow ≥ M(b) = d − wt(b)/2 per class. The checker computes
this directly: the cover-cycle condition is the affine system
E·v₀ = rhs(b) with E class-independent (one 180×180 matrix per cell);
solutions decompose into ≤ 2^{k/2} sectors of ker E / im S with
per-sector cover-boundary triviality; violations are hunted by two
complete lanes — restricted-support enumeration with meet-in-the-middle
subset sums over E's reduced columns (M − 1 ≤ 4: the heavy strata), and
the multi-offset coset-BZ over the same κ=88 windows as §3 (light
strata). Every candidate is re-verified (solution check, slice identity,
non-boundary) before it could count as a violation.

**Validation**: f2a6:y — 113/113 PASS in 1.7 s, reproducing the Lean
theorem `dangerousFloorNZ_of_lightClassification`'s instance (A17 §6.5),
with the lane split (1+7+36+69) matching the census histogram exactly;
plus a soundness control that hand-builds a genuine dangerous logical
(τ(u) + S̃(z̃), weight 20, stamp shadow, overflow 7) and confirms the
checker FINDS it at an inflated target.

**Results** (`data/a30/rungs_*.json`):

| cell | classes | verdict | wall | lanes |
|---|---|---|---|---|
| 37a70e02:x | 2,203 | **ALL PASS** | 7.0 s | 1 bz + 6 r≤4 + 42 r≤3 + 54 r≤2 + 478 r≤1 + 1622 r≤0 |
| 5e50a976:x | 2,371 | **ALL PASS** | 7.2 s | 1 + 6 + 43 + 55 + 503 + 1763 |
| 5e50a976:y | 2,371 | **ALL PASS** | 7.0 s | (same strata) |

**Sector-completeness accounting** (why this closes d = 20): every
nonzero cover class's shadow is either (i) a stabilizer — the dangerous
sector: b = 0 via LogicalFloor 10 (§5, certified), light b via census
completeness (A28) + these per-class certificates, heavy b via the slice
identity — or (ii) a base logical whose class lies in im p₁ = im δ₂
(A15-P3/A17), i.e. in one of the three Δ-classes whose full cosets §3
certified ≥ 20, so the projection bound wt(v) ≥ wt(p(v)) ≥ 20 closes it.
With the lifted weight-20 witness (§5): **d([[360,4,20]]) = 20 at
certificate tier, end-to-end, for 37a70e02:x and 5e50a976:x/y** (X-side;
Z by the BB transpose duality, teaching doc §"one side suffices", as in
the [[300,8,16]] packaging). Updated per-code totals: 17.4 min and
23.1 min of wall-clock — the rung pass added ~7 s per cell.

## 6. Deliverables and follow-ons

Shipped: `scripts/a30_coset_bz.py` (validate + decide + the C kernel,
budget-enforced), `data/a30/` (certificates, witnesses, validation).

Follow-ons, in value order: (1) Lean packaging of one Z₁₅×Z₆ cell per
the A15 pattern — the natural first unconditional d = 20 BB theorem of
the program; (2) fold coset-BZ into `bb_lab.fibering`'s
`safe_floor_certify` as the fallback lane for frame-dead codes (A29 §6.2
refusals become decisions); (3) sharded/distributed runs for d = 12-base
hunts (W = 22 ⟹ r-pair (11,10) ≈ 60× today's nodes — cluster
territory, or the theory route); (4) the docket's doubling rows can now
be *executed* — build the [[360,4,20]] instances and SAT-witness the
upper sides as an external sanity row.

## 7. Session 3 (2026-08-10): the front-end — cover code in, d = 20 out,
## composed around Tandem

Follow-on 2 delivered as a library + CLI + UI lane: give it the COVER
code and it does everything this note did, autonomously.

**Library** (`src/bb_lab/cosetbz.py` + `src/bb_lab/doubling_certify.py`;
the a28/a30 session scripts stay frozen as artifacts):
`detect` (axis literal-lifts up to per-block translation normalisation —
a coordinate permutation of the code — plus (R) via k-preservation),
`d_side_exact` (progressive-ladder logical-coset BZ, both CSS sides,
witnesses), `census` (vectorised translation canonicalisation),
`safe_floor`, `RungEngine`/`rung_pass` (the §5.5 checker), witness lift
(τ(u) verified nontrivial — now CHECKED, not assumed), and `certify`
(orchestration, budgets on the monotonic clock, node-estimate guard:
sweeps beyond ~2.5·10¹² nodes are refused with the shard/theory-route
message instead of ground through).  Every verdict carries its claim
tier and a Tandem block.

**The Tandem composition** (per the A27/descent-SAT verdicts: compose
around the solver, never inside it): a CERTIFIED verdict suggests
`-init-lb=<floor> -cost-step=2` with the acknowledgements — the
certified floor deletes Tandem's proof phase, so its run only has to
FIND the weight-2d witness (minutes) and stops at the floor: the
`-fiber-lb` hybrid idea landed soundly.  On detection/(R)/scale failure
or a refutation, the verdict routes to monolithic Tandem as the
fallback lane.  UI: `bb-lab ui` gains a "Certify as doubling cover"
action (`/api/certify`, SSE stage stream, tier-labelled verdict card,
one-click Tandem cross-check with the flags + acknowledgements
prefilled); server-side soundness gating unchanged.

**Validation** (`tests/test_doubling_certify.py`, 5/5 in ~5 s): docket
base detection (both cover shapes); pair72 cover [[72,4,8]] certifies
d = 8 end-to-end in seconds; f2a6:y rung pass 113/113 through the
library path; the by90 rung cover (Z₃₀×Z₃, the [[180,8,12]] Bravyi
rung whose base is the [[90,8,8]] bottom) comes back DOUBLING-REFUTED
(the A14 §13 freeze) — no false certificates; even-weight inputs
refused.

**End-to-end on the [[360,4,?]] covers** (cover code as the ONLY
input): Z₃₀×Z₆ `1+y+x` / `y⁴+x+x¹¹y²` → detect base [[180,4,10]],
d_base = 10 both sides (~1 s), census 2,203 classes with the A28
histogram **bit-identical**, safe floor CERTIFIED, rung pass
2,203/2,203 PASS (lane strata exactly §5.5's), witness lift
established ⟹ **CERTIFIED d = 20 in 1,220 s ≈ 20.3 min** (some CPU
shared with a concurrent session).  Z₁₅×Z₁₂ (the 5e y-cover) run with
the live Tandem cross-check lane — see `data/certify_runs/`.

Tandem itself was built this session (build_maxcdcl.sh, fork verified:
`-cost-step` present, 29 options discovered) and relocated to the
`third_party/maxcdcl/` path the UI probes.

**The Tandem cross-check, measured at both scales.** The second cover
(Z₁₅×Z₁₂, CERTIFIED d = 20 in 946 s) ran the live lane three ways:

| lane | outcome |
|---|---|
| pair72 cover [[72,4,8]], `-init-lb=8` | **d = 8 CONFIRMED in 0.02 s**, `agrees: true`, witness re-verified by the solver path — the composed loop works end-to-end |
| n = 360, `-init-lb=20` unseeded | **no incumbent in 900 CPU-s** (~1.7M conflicts/phase): the monolith cannot even FIND a weight-20 witness in the budget the front-end needs for the entire certificate |
| n = 360, `-phase-file` seeded with the certified witness | no incumbent, **conflict counts identical to unseeded** — MaxCDCL's VSIDS-initialization pass (first 10⁴ conflicts) re-derives polarities and wipes the seed; phase seeding is a null lever in this fork as built (future fork tweak: apply phases after init) |

Honest summary for the UI/story: the cross-check lane is real and
instant at solver-reachable sizes; at n = 360 the certificate's own
in-process-verified witness IS the upper side, and the solver lanes are
fallback-only. Verdicts: `data/certify_runs/*/verdict.json` (three
committed, incl. both n = 360 cross-check records).

## Appendix: verification map

| claim | check |
|---|---|
| validation battery 10/10 | `uv run python scripts/a30_coset_bz.py validate` (~4 min) |
| the three floor-20 certificates | `uv run python scripts/a30_coset_bz.py decide --threads 10` (~15 min) |
| Z₂₁×Z₃ cell + node counts | `data/a30/decide_16884e06.json` (rerun: §5 driver in the session log, 0.6 s) |
| exactness witnesses | `data/a30/exactness_witnesses.json` |
| node-count invariants | asserted in-run (`run_window`: nodes == Σ C(κ,s)) |
| seam offsets / orbit reps | `bb_lab.fibering.seam_offsets`, `kernel_orbit_reps` (A29, battery-exercised) |
