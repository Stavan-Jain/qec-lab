# A27 — Safe-floor deletion: generality tiering, cost anatomy, and a fibering feasibility probe for the Z₁₅×Z₆ [[180,4,10]] docket UNKNOWNs

> **STATUS.** Synthesis note + one new computation, born in the teaching-PDF
> session of 2026-08-06 (`docs/teaching/bb-doubling-theorem.{tex,pdf}` in this
> repo; the Q&A that produced this note is not otherwise recorded). §§1–2 are
> *consolidation* — a tiered answer to "how general are the safe-floor
> deletion techniques?" and a one-place cost anatomy, both assembled from
> existing notes with pointers. §3 is **new**: a numeric feasibility probe of
> the A22 (ε,δ)-fibering ansatz on the three A17-docket UNKNOWN cells
> (`37a70e02e003d1de:x`, `5e50a9765a02eb70:x/y`, Z₁₅×Z₆ `[[180,4,10]]`,
> safe-floor-20 targets). Probe is numeric-only (pure Python, independent of
> `bb_lab`), semisimplified-spectrum only, no Lean, and makes **no claim about
> whether floor 20 is true** — it assesses whether the A22 argument *shape*
> applies. Falsification gates for any follow-on are in §5.

Companions: `A22_analytic_classification.md` (the ansatz this ports),
`A23_analytic_seam_floor.md`, `A17_d7plus_doubling_hunt_plan.md` (the docket),
`A14_safe_floor_criterion_plan.md` (Prop A14.1), `A8` §§4.1/5 (the transfer
verdict), `descent_sat_results.md` (solver-lane acceleration),
`docs/gross-distance-extensibility.md` §7 (the pre-A22 ancestor of §1).

---

## 1. Generality tiering: how far the safe-floor deletion travels

"Getting rid of safe-floor certification" decomposes into three tiers of
decreasing portability. (This supersedes the extensibility doc's §7 framing,
which predates A22/A23.)

**Tier 1 — fully general, proven once (the reduction).** Under (R), for every
free-Z₂ BB cover, with no per-code work: safe classes are exactly
`Δ(Ann(A) ∩ Ann(B) ∖ 0)`, of dimension exactly `k/2`; explicit seam-carry
representatives (`seamC`); coset minima constant on G-translation orbits
(Prop A14.1). Plus: (R) ⟺ k-preservation ⟺ Bezout membership (A12), the
Z_{2^r} tower extension (A13), and general failure fragments (pushforward
monotonicity `d_safe ≤ d̃_safe`; the 2d−2 parity result, A17 deficit wall).
All parametric in Lean (`BBCover`/`BBDoubling`). The *bookkeeping* of
safe-floor certification — which cosets, which representatives, up to which
symmetry — costs nothing on a new code.

**Tier 2 — portable recipe, per-code execution (the fibering).** The move
that *deletes* the certification — CRT/Fourier decomposition + fibering the
weight inequality — has succeeded in two structurally different regimes:
gross (even Z₂² frame, F₄ components, slot-frame walk, A4 §§9–13) and
cover300 (odd base group Z₅×Z₁₅; A22 collapsed the census 9.6 h → 2 min,
A23 collapsed the seam floor to the single inequality
`|A⋆f+e₀| + |B⋆f| ≥ 16`). The Lean scaffolding (pivot/rank certificates,
certified sweeps, transport) is code-agnostic. Three things must be redone
per code, with no meta-theorem bounding their difficulty: the fiber/frame
choice, the layer dictionary (gross's flat 6 became 8/12/14 over Z₃×Z₇), and
the per-component case analysis, whose size is governed by which fields the
CRT produces.

**Tier 3 — provably non-portable (gross's engine core).** The co-point
value-rigidity works because F₄ has exactly three nonzero elements; A8 §2
shows it is simply false over F₈/F₆₄ — which is exactly why the
`[[336,12,12]]` safe floor is A8's declared open core although the doubling
itself holds. A8 §5's verdict stands: *the technique transfers as an
architecture; the quantitative cores are code-specific.*

**Boundaries** (partly per the GB-generalization audit, 2026-07-28, on its
own branch): abelian groups only; univariate gets easier; `k` must stay
modest (2^{k/2} classes, transport divides by |G| only); polynomial weights
odd for the parity arguments (the A11 Entry-3 w=4 counterexample makes that
scope provably necessary); (R) is a genuine prerequisite.

**Scorecard:** two complete analytic deletions (gross, cover300 — both
kernel-checked), one consciously open (Z₆×Z₁₄), one where brute sweeps
sufficed (pair72, dim K = 2).

---

## 2. Cost anatomy, one place

**Worst-case anchors.** I2 (light-stabilizer census) is low-weight codeword
enumeration + completeness — min-distance-type, NP-hard (Vardy 1997). I4
(coset floors) is the coset-weight problem — NP-complete
(Berlekamp–McEliece–van Tilborg 1978). Both certification directions are
coNP-flavored absence proofs; the template halves instance size, collapses
I4's multiplicity to `2^{k/2}−1 → orbit reps` (constant along a doubling
family, since (R) pins k), and swaps the dangerous sector's 2d-threshold
work for a threshold-d statement + O(1) pigeonholes. The exponential core
(n^{O(d)} per absence proof) is concentrated, not dissolved.

**The (150, 14) inversion.** At f2a6f17e's scale, the *solver* costs invert
the usual I4-most-expensive ranking: the I2 census cost 9.6 h bulk (top
rungs dominating), while I4's coset queries ran inside the docket's 1200
s/query CMS budget (proof-grade single-rep re-run: kissat 9,506 s, 6.85 GB
DRAT). Both are the same species at their hardest point — absence-proving
UNSAT at the top weight — and both blow up together two rungs higher: at
(180, 18) the census extrapolates to days-plus and the coset queries
already DNF'd (2 h+, no verdict, the docket's three UNKNOWNs).

**One frame, two floors.** In both fully-analytic instances the frame is
the capital investment and the floors are marginal: gross's CRT layer frame
(A4 §3) powers both the §6.3 classification and the §§9–13 M-im program;
f2a6's (ε,δ)-fibering was built for A22 (census) and A23 rode it a day
later (seam floor). Expect the same shape on any new target.

**Solver-lane accelerations** (details in `descent_sat_results.md` and the
shard/descent memories): strengthened monolith 4–6×; descent/pushforward
lane +2.8–4.6× (y-deck, NOT the doubling deck); CMS for the XOR-heavy coset
instances; Tandem three-part verdict — negative as clause-level drop-in
(0.6–0.9×), the current record-holder for *global* settlement (bb_288
d = 18 exact in 6,924 s), and one live hybrid idea (`-fiber-lb` bound
arithmetic, value thesis n ≥ 288). Constant factors don't move the wall;
frames delete it.

---

## 3. NEW — fibering feasibility probe: Z₁₅×Z₆ [[180,4,10]]

### 3.1 The targets

The three A17-docket UNKNOWN cells (safe-floor-20, w_query = 18, CMS 2 h+
per query, no verdict; rows from `data/a17/docket_decision.jsonl`):

| id | axis | G | A | B |
|---|---|---|---|---|
| `37a70e02e003d1de` | x | Z₁₅(x)×Z₆(y) | `1 + y + x` | `y⁴ + x + x¹¹y²` |
| `5e50a9765a02eb70` | x, y | Z₁₅(x)×Z₆(y) | `1 + y + x` | `y⁴ + x⁸y² + x¹³` |

Both `[[180,4,10]]`, doubling target `[[360,4,20]]`.

### 3.2 Fiber choice and what transfers verbatim

Fiber `z := x³` (order 5; 2 primitive mod 5), sites `Z₃(w)×Z₆(y)` (18
sites), with `x = z²w²` (so `xᵏ ↦ w^{2k mod 3}` under ε). Group-level
ingredients of A22 §0 transfer **verbatim**: the `F₂ × GF(16)` fiber CRT,
the six-entry exact weight table `W(ε,δ)`, and the sweep pattern. Only odd
prime fibers with 2 primitive (p = 3, 5, 11, 13, 19, …) keep the weight
formula exact (two CRT factors ⟹ fiber value determines fiber polynomial);
Z₅ is the only such fiber here that keeps the site count at 18 (Z₃ fibers
give 30 sites and a ~130× bigger sweep).

**Structural novelty vs f2a6:** the site group contains a Z₂ (f2a6's sites
Z₁₅ were odd) ⟹ site-side lemmas need a gross-style radical layer.

### 3.3 Probe results (script: `scripts/a27_fibering_feasibility_z15z6.py`, ~seconds)

| check | result |
|---|---|
| P1 params | both codes n=180, **k=4, dim K = 2** (matches docket k; dim K = k/2 ✓ A14) |
| P2 δ-kernel | K basis 2/2 fixed by the fiber idempotent ⟹ **δ-kernel = 0; ALL safe classes are ε-sector** (forced: any nonzero GF(16)-module has F₂-dim ≥ 4 > 2) |
| P3 ε-quotient | **identical for both codes**: `A_ε = 1+y+w²`, `B_ε = y⁴+wy²+w²` over Z₃×Z₆ — the two codes are *fiber-twists*; quotient = `[[36,4,2]]`, dim ker ∂₂ = 2, **d̄ = 2** (exhaustive over 2²⁰ cycles) |
| P4 ε-spectrum | `A_ε` NOT a unit — 2 character zeros, at (χ_w,χ_ȳ) = (ω,ω), (ω²,ω²) ⟹ **f2a6's free-h collapse fails**; ε-sector is a bona fide mini-BB-code |
| P4 δ-spectrum | `Ã` a **unit** (0 zeros ⟹ transfer T global, no η₀ special point; unit in the full ring by Nakayama); both `B̃` zero-free; **no common zeros** (consistent with δ-kernel 0) |
| P5 sweep | budget 2d−2 = 18, active site ≥ 2 (robust to any h-constraint) ⟹ ≤ 9 active of 18 ⟹ **Σ C(18,≤9) = 155,382** subsets (~9.5× f2a6's 16,384) |

### 3.4 Verdict

**The architecture ports.** δ-side strictly *cleaner* than f2a6 (global
transfer, trivial δ-kernel, strictly overdetermined systems below |S| = 9).
One genuinely new chapter: the **constrained-h taxonomy** — the ε-sector is
the shared `[[36,4,2]]` quotient code with its own kernel (which is exactly
where all `2² − 1 = 3` safe classes live), so the optimal-h accounting must
run over a constrained, kernel-bearing mini-code instead of a free
coordinate. `d̄ = 2` means the ε-side admits light configurations: the
classification tree will be bushier than f2a6's, and any floor-20
conclusion must come from fiber-weight inflation (pure-ε occupied site
costs 5; mixed sites 1–4) plus δ-costs. Because the ε-core is shared, one
ε-chapter serves all three docket cells.

**Outcome symmetry.** The fibering is a decision procedure: it may
*certify* floor 20 or *refute* it with a light safe vector. Either decides
cells the solver could not touch; a refutation feeds the deficit-wall
dataset (P3 of A17).

**Estimates** (calibrated on A21–A23's recorded 2-day arc): decision on the
three cells ≈ 3–7 sessions; full unconditional `[[360,4,20]]` theorem if
floors certify ≈ 8–14 sessions (adds `logicalFloor_10`, rungs, packaging).
Compute is laptop-scale throughout: Phase-0 sweep ~10 min–1.5 h single-core
(embarrassingly parallel), Lean builds ~1–3 h (pivot certs, no
`native_decide` floor monsters), optional SAT spot-checks minutes each.
Contrast: the solver route already spent ~6 CPU-h on these cells for zero
verdicts, and its census at (180, 18) extrapolates to days-plus.

---

## 4. Caveats (read before building on §3)

1. P4 ran on the **semisimplified spectrum**; the Z₂ radical layer in the
   site algebra is unprobed (unit-ness lifts by Nakayama; the sweep's
   linear systems will carry doubled GF(16) dimensions).
2. The probe is an independent pure-Python implementation; its only
   cross-checks against prior data are k = 4 (docket) and dim K = k/2
   (A14). No V6-analog yet: the decomposition identities have **not** been
   checked against actual light vectors of these codes.
3. `d̄ = 2` was computed by exhaustive enumeration of the quotient's 2²⁰
   cycle space — trustworthy but worth one independent re-derivation
   before it becomes load-bearing.
4. No statement here bears on whether floor 20 is *true*.

## 5. Falsification gates + next steps (if pursued)

1. **V6-analog** (first gate): decompose a sample of genuine light
   boundaries of `37a70e02` per the ansatz (constrained SAT or direct
   enumeration at low |b|); any mismatch kills the port cheaply.
2. **ε-chapter**: constrained-h cost taxonomy over the `[[36,4,2]]`
   quotient; gate = a finite, surveyable case table (bushiness budget:
   if the taxonomy exceeds ~10× f2a6's, stop and reassess).
3. **Full sweep** (155,382 subsets, GF(16) pivot certificates) → derived
   census; spot-check against constrained SAT at low rungs.
4. **Seam floor at 20** for the 3 ε-sector classes via the frame
   (A23-analog) — or the refutation witness, whichever the frame yields.
5. If certifying: `logicalFloor_10` (A21-analog, threshold 9 on n=180),
   dangerous rungs (A15 Stage-5 pattern), (R)/Bezout check, Lean packaging
   (`Z15Z6…` module family; watch the naming-collision convention,
   A17 §"Naming note").

## Appendix. Verification map

| claim | check |
|---|---|
| P1–P5 of §3.3 | `python3 scripts/a27_fibering_feasibility_z15z6.py` (~seconds, no deps) |
| docket rows / polynomials | `data/a17/docket_decision.jsonl` (grep the two ids) |
| f2a6 calibration numbers (9.6 h census, 2-min re-derivation, 14-min build) | `A22_analytic_classification.md` §§1, 5 |
| docket run shapes, 1200 s/7200 s budgets, kissat 9,506 s / 6.85 GB DRAT | `A17_d7plus_doubling_hunt_plan.md` §6 |
| solver-acceleration numbers, Tandem verdicts, bb_288 settlement | `descent_sat_results.md` §§1–4, tail |
| Tier-1 theorems | `A14_safe_floor_criterion_plan.md` §2 (Prop A14.1), `A12_deck_homotopy_R.md` §3 |
