# A29 — The general fibering engine: portable safe-floor certification for BB codes

> **STATUS.** Session 1 (2026-08-06, branch `claude/fibering-bb-codes-5c9ff3`).
> **Renumbered A28 → A29 at the 2026-08-06 merge**: this line and the LSC
> line (`A28_light_classification_theory.md`) claimed A28 in parallel
> worktrees; the LSC registry claim was 31 min earlier, so this note,
> its scripts (`a29_*.py`), and its data dir (`data/a29/`) were
> renumbered — the July double-A15 rule.
> Goal (user brief): generalize the A22/A23/A27 fibering techniques so
> safe-floor certification — the most expensive part of discharging the
> doubling theorem on a new code — becomes a portable, per-code-cheap
> procedure, ideally covering the whole corpus. Deliverables: the theory
> (§1–§3, all engine lemmas proven here in elementary form), the engine
> (`src/bb_lab/fibering.py`, code-agnostic), validation against every
> recorded ground truth (§4), and deployment over the A17 docket +
> curated targets (§5). §6 records honest limits and the follow-on
> program. Numeric-only session: no Lean; Lean-feasibility assessed in §7.

Companions: `A22_analytic_classification.md`, `A23_analytic_seam_floor.md`
(the instance-proofs this generalizes), `A27_safe_floor_generality.md`
(the tiering + probe this supersedes in part — see the §3.3 erratum),
`A14_safe_floor_criterion_plan.md` §2 (Prop A14.1, consumed as-is),
`A8_doubling_extension_writeup.md` §4.3 (the open core targeted in §5).

---

## 1. What was actually general in A22/A23 (and what was not)

Dissecting the two shipped instance proofs yields a clean split.

**General with proofs (now engine lemmas, §2):** the fiber partition and
the exact weight formula (for ANY odd fiber order q, not just 2-primitive
primes — see L1; A27 §3.2's restriction is unnecessary); the ε-quotient
ring map; the parity-link lemma whenever `B_ε = x̄^m·A_ε`; the per-site
cost tables including the offset-twist (τ) variant; the cost-invariant
kernel `K₀`; the mask-sweep completeness argument; the ε-fiber structure
(base + joint quotient image) making exact ε-minimization a linear
functional over an affine F₂-space.

**Code-specific in the originals, now automated:** the fiber choice
(A22 §0 chose `z = y³` by hand; the engine enumerates every cyclic odd
subgroup and scores each), the link detection (f2a6's `B̄ = x̄·Ā` was
noticed by hand; the engine searches all site-monomials), and the seam
offsets (A23 §6's P/Q/e₀ machinery is replaced by the always-available
seamC carry chain of Prop A14.1(3), computed by cover-window evaluation).

**Genuinely not general:** the *feasibility* of the certified sweep.
It depends on the (code, fiber, class) triple through three quantities —
the site count S, the budget 2d−2, and the δ-data rank r — and on two
structural bits: whether some fiber carries the ε-monomial link, and how
much ε-twist (τ) the safe classes carry on that fiber. §3 makes this
precise; §5's census measures it over the corpus.

## 2. The engine lemmas (all elementary, all falsify-first-tested)

Setup: G finite abelian, A, B ∈ F₂[G] of odd weight, z ∈ G of odd order
q ≥ 3, sites = G/⟨z⟩ (S of them), fibers = cosets ordered by the section
`rep(s) + k·z`. All statements are proven in-line here; the engine
asserts each on random inputs at frame construction or in validation
scripts.

**L1 (residue coordinates; exact weight formula for every odd q).**
Define `D: F₂^q → F₂^{q−1}`, `D(t)_j = t_j + t_{q−1}`. Then ker D =
{0, N} (N = all-ones), (parity, D) is a bijection, and the two lifts of
a residue r weigh |r| and q − |r|, distinguished by parity (q odd).
Hence for any chain block w: `|w| = Σ_sites W(r_s, ε_s)` with
`W(r, e) = |r|` if `|r| ≡ e (mod 2)` else `q − |r|` — closed form, no
CRT fields, no tables. Corollaries: `w0(r) := min(|r|, q−|r|) ≥ 1` for
r ≠ 0; the ε-flip gap `|q − 2|r|| ≥ 1`. *This kills the A27 §3.2 claim
that only 2-primitive prime fibers keep the weight formula exact: the
field splitting of F₂[Z_q] is scaffolding; the invariant content is
the residue map.*

**L2 (ε is a ring map; parity link; τ).** The fiber-augmentation
ε: F₂[G] → F₂[G/⟨z⟩] is a ring homomorphism, so ε(P⋆f) = P_ε ⋆ ε(f).
If `B_ε = x̄^m · A_ε` for a site-monomial x̄^m (the **link**), then for
every f the fiber parities satisfy `ε(B⋆f)_{s+m} = ε(A⋆f)_s`. For an
affine family (u₀ + A⋆f, v₀ + B⋆f) the twist
`τ(s) := ε(u₀)_s + ε(v₀)_{s+m}` is **f-independent** (both sides pick
up the same correction) — τ is an invariant of the coset, i.e. of the
safe class on that frame.

**L3 (state cost bounds).** At a linked site-pair (s, s+m) with residues
(r₁, r₂):
- τ = 0, pair active ((r₁,r₂) ≠ 0): cost ≥ 2. (Parity match: the even
  lift of a nonzero residue weighs ≥ 2; the odd lift of the zero residue
  weighs q.)
- τ = 1: cost is **odd**; = 1 forces one block's residue to vanish
  (states U0/V0); both-blocks-active forces cost ≥ 3; inactive forces
  cost q. There is **no** even-cost τ=1 site.
- Unlinked (unpaired mode): every δ-active (site, block) position costs
  ≥ w0 ≥ 1.

**L4 (the state sweep; completeness).** Fix a weight target t and let
budget = the largest value ≤ t−1 of the family's weight parity (weights
are ≡ |u₀|+|v₀| mod 2 since |A|,|B| odd). Assign each linked pair-site
one of {OFF(0, both-blocks residue-vanishing), ON(2, no conditions)}
(τ=0) or {U0(1, u-vanishing), V0(1, v-vanishing), ON(3, none)} (τ=1;
OFF(q) is dominated by U0/V0 and dropped), and each unpaired position
{OFF(0, vanishing), ON(1, none)}. Any family member of weight ≤ budget
satisfies the conditions of the state assignment read off from its own
residues, whose state-costs sum ≤ budget. Enumerating all
budget-feasible, domination-maximal assignments and, per assignment,
solving the affine F₂ system in f and checking every solution modulo
K₀ = {f : A⋆f, B⋆f both ε-only} is therefore a **complete** search for
weight-≤-budget members. (Domination: upgrading any conditioned unit to
its active state weakens conditions and raises cost — each feasible
assignment is dominated by a maximal one; K₀-invariance of all costs:
K₀ shifts change no residues and act on ε-parts by the joint quotient
image, which the exact ε-minimization already ranges over.)

**L5 (exact ε-minimization is linear).** For a fixed δ-config, the
reachable ε-patterns are (ε(u), ε(v)) + {(A_ε ḡ, B_ε ḡ) : ḡ ∈ F₂[sites]}
(proof: K₀-elements κ act by (A_ε ε(κ), B_ε ε(κ)), and every ḡ is ε(κ)
for κ = Σ ḡ_s N_s). The weight as a function of the flip pattern is
Σ W(r, e₀ ⊕ flip) = const + Σ ±gap·flip_pos — linear over F₂ — so the
exact coset minimum per δ-config is the minimum of a linear functional
over an affine subspace, enumerable (Gray code) when the ε-image is
small, and relaxable to Σ w0 (sound, ε-free) otherwise.

**L6 (safe classes; offsets).** Under (R) (checked per run as k̃ = k):
safe classes = Δ(ker ∂₂ ∖ 0), δ₂[ζ] = [seamC ζ], coset minima constant
on G-orbits (Prop A14.1). The engine computes seamC ζ directly as the
carry window of the axis-doubled cover boundary of the canonical lift,
and sweeps one offset per G-orbit rep of ker ∂₂ ∖ 0. Frames may be
**mixed across classes** (floors are statements about the code, not the
frame): each class is tried on every feasible frame until decided.
Refutation witnesses are re-verified end-to-end: dressed chain weight,
family membership by construction, and non-boundary (else the safe-class
setup itself is flagged).

## 3. Feasibility calculus (what "portable" quantitatively means)

For a frame with S sites, fiber order q, δ-rank r (= |G| − dim K₀) and
budget b = 2d − 2:

- **paired mode** (link present, τ ≡ 0 classes): sweep size ≈ C(S, b/2),
  overdetermination margin ≈ 2(S − b/2)(q−1) − r. Feasible when
  C(S, b/2) fits the mask cap and the margin is not too negative.
- **τ-twisted classes**: τ1 sites contribute cost-1 states carrying
  one-block conditions; heavy τ classes lose roughly (q−1)·#τ1
  conditions at fixed budget — classes with many τ1 sites on every
  linked frame are the hard residue (38d3c884 below).
- **unpaired mode**: positions 2S, per-position ≥ 1 ⟹ sweep
  ≈ C(2S, b): only small-S or small-b codes.
- **vacuous**: b/2 ≥ S (paired) or b ≥ 2S (unpaired) — the fiber is too
  coarse for the budget. This is why gross (S = 12, b = 22) is out of
  reach of THIS engine on every fiber: consistent with history (gross
  needed the slot-frame walk).

**Erratum for A27 §3.3 P5.** The probe's "active site ≥ 2 (robust to
any h-constraint) ⟹ ≤ 9 active of 18 ⟹ 155,382-subset sweep" for the
Z₁₅×Z₆ cells presumed the ε-parity pairing, but no site-monomial link
exists for the q=5 fiber there (`A_ε = 1+y+w²`, `B_ε = y⁴+wy²+w²` are
not monomial multiples — L3's ≥2 bound is a theorem exactly in the
linked case). A cost-1 site (monomial fiber one block, empty fiber the
other) is not excluded by any h-constraint. The honest sweep on that
frame is the unpaired one (C(36,18) ≈ 9·10⁹ before symmetry), or the
ε-chapter route. The A27 cost *estimates* for those cells are therefore
optimistic; the architecture verdict ("ports, with a new ε-chapter")
stands.

## 4. Validation record (engine vs every recorded ground truth)

| target | ground truth | engine | agree |
|---|---|---|---|
| pair72 base [[36,4,4]] Z₃×Z₆, x-axis | exact safe minima 8/8/8, raw seams 12 (A14 data) | certify@8 (20/20 masks inconsistent); refute@10 with verified weight-8 witnesses; raw seam 12; 1 orbit | ✓ |
| f2a6 [[150,8,8]] Z₅×Z₁₅, y-axis, floor 16 | A23 §9.2: 300 cons / 6135 incons of C(15,7), kdims {0:120, 4:180}, 3000 reps, min 16 | **identical counts**, CERTIFIED, 2.8 s | ✓ bit-for-bit |
| f2a6 floor-18 control | SAT-tight d_safe = 16, split (10,6) | REFUTED with verified weight-16 witnesses (125), min_exact 16 | ✓ |
| f2a6 boundary floor | A22: d(im ∂₂) = 6; the |b|=6 class = one free G-orbit (75 translates) | certify@6; refute@8 with exactly 75 violations | ✓ |
| f2a6 census ≤ 14 | A22: 113 classes, |b|-histogram {6:1, 10:7, 12:36, 14:69} | **113 classes, identical histogram**, 65 s | ✓ exact |
| A22/A27 structural probes | link site = x̄ (f2a6), kernel orbit = 1 rep of weight 40 = e₀, raw seam 18 | same | ✓ |

### 4.1 Census validation

`FiberSweep(FiberFrame(f2a6, z=y³)).census(14)`: the survivor sweep
(6435 homogeneous leaves, 2.59M kernel reps — same order as A22 V7's
1.38M α's) followed by the batched ε-fiber expansion (the flip-cost
functional is linear, so all 2^15 flips of a δ-config are priced by one
matvec) and diagonal-translation canonicalization returns **exactly 113
classes with the A22 histogram** {6:1, 10:7, 12:36, 14:69} in 65 s.
This replaces the enumeration-completeness role of the 9.6 h SAT run
with general-engine machinery (A22's hand-built 2-min rederivation is
subsumed).

## 5. Deployment (docket + curated targets)

### 5.1 NEW DECISIONS — three docket UNKNOWNs certified, two pinned exact

The A17 docket carried **7 distinct UNKNOWN safe-floor cells** (CMS 2 h+
per query, no verdict; A27 §3 discussed only the three Z₁₅×Z₆ ones —
the docket also holds three Z₅×Z₁₅ [[150,8,8]] cells and one Z₂₁×Z₃).
Engine results (`data/a29/docket_engine*.jsonl`):

| cell | code | floor | engine verdict | time | extra |
|---|---|---|---|---|---|
| `ac46bbea:y` | [[150,8,8]] Z₅×Z₁₅, B = 1+y¹¹+xy² | 16 | **CERTIFIED** | 1.4 s | refute@18 ⟹ **d_safe = 16 EXACT** |
| `38d3c884:x` | [[150,8,8]] Z₅×Z₁₅, B = y⁹+y¹²+x²y⁴ | 16 | **CERTIFIED** | ~20 s | refute@18 ⟹ **d_safe = 16 EXACT** |
| `38d3c884:y` | same code, y-axis | 16 | **CERTIFIED** | ~20 s | certify@18 too ⟹ **d_safe ≥ 18** |

Fiber anatomy (verified by direct τ inspection): `ac46bbea:y` carries
the ε-link on the standard y³ fiber with τ ≡ 0 — a straight f2a6-shape
port, hence 1.4 s. `38d3c884` has **no link on the standard fiber**;
the enumeration surfaces two diagonal linked fibers, z = (1,3) with
τ-weight 8 (an ε-twisted seam class — a phenomenon f2a6 never
exhibited, its kernel being pure δ-sector) and z = (1,12) with τ ≡ 0;
the per-class frame fallback decides on the clean frame, while the
τ-state sweep of L3/L4 is what makes twisted frames usable at all.
Hand-picking fibers (the A22 workflow) would have missed all of this.
The exactness pins are solver-free in both directions: sweep
certificates below, verified non-boundary witnesses above.

### 5.2 THE A8 OPEN CORE IS CLOSED

`[[168,12,6]]` over Z₆×Z₁₄ (A = 1+y+x³y³, B = 1+x+x²y⁷), x-doubling to
the SAT-exact `[[336,12,12]]` (A8 §3). A8 §4.3 declared the confined
floor — every base 1-cycle in a nonzero Smith class has weight ≥ 12 —
its open core, expecting a re-derivation of the gross slot-cost engine
with a heterogeneous 8/12/14 layer dictionary over F₈/F₆₄ where the
F₄ co-point rigidity provably fails (A8 §2, the Tier-3 obstruction of
A27 §1). The engine **certifies the floor at 12 in 10.2 s** on the q=7
unpaired frame ((R) re-checked: k = k̃ = 12; 3 kernel orbit classes;
relaxed δ-costs alone suffice — min 12/15/15, zero ε-suspects; leaf
counts 18450/268/268, all kdims ≤ 6). The Tier-3 wall was never on the
path: **the fibering needs no value rigidity at all**, and the "genuinely
undeveloped mathematics" A8 anticipated is unnecessary for this floor.
With condition (1) d(cover) = 12 SAT-exact, (2) k-preservation = Bezout
(A12), and now (3) the safe floor certified, the `[[336,12,12]]`
doubling has all template inputs at engine/certificate level — the
first d = 6 → 12 doubling with a solver-free safe floor.

### 5.3 Validation rows and the engine boundary

Two full-docket passes are recorded. `docket_engine_full.jsonl`
(pre-DFS scorer): f2a6f17e:y — the one Z₅×Z₁₅ SF-CERTIFIED row in
engine reach — re-certifies in 1.3 s, agreeing with its CMS UNSAT@14;
Z₂₁×Z₃ rows report NO-FEASIBLE-FIBER. `docket_engine_final.jsonl`
(current engine; 17 of 21 rows before the run was cut): the three
Z₅×Z₁₅ UNKNOWNs re-certify (0.3 s / ~19 s / ~19 s), 13 Z₂₁×Z₃ rows
return honest UNDECIDED after burning their node caps (~2–10 min
each), and the first Z₁₅×Z₆ floor-20 row spent 3,754 s across its two
unpaired frames before returning UNDECIDED — the measured shape of the
v1 boundary (the remaining three Z₁₅×Z₆ rows were cut as
same-diagnosis stragglers). Every SF-CERTIFIED row the engine can
reach agrees with its SAT verdict; no row anywhere disagrees. The
boundary rows have a precise diagnosis each:

* **Z₂₁×Z₃ (q=3 frames):** no fiber carries the ε-link (all ε-images
  collapse to monomial-vs-trinomial pairs), and in unpaired mode the
  relaxed bound is far too weak (observed min_relaxed = 8 against
  floors of 16) because the ε-penalty carries most of the weight —
  hundreds of thousands of suspects each needing an exact ε-pass, and
  the ε-image (~2²¹) is over the enumeration cap. q=7 frames die the
  other way (δ-kernels ≈ 2²⁶ per leaf).
* **Z₁₅×Z₆ floor-20 (A27's cells):** no linked fiber (the §3 erratum);
  the unpaired q=5 DFS tree is ≈ C(36, ≤18) ≈ 3·10¹⁰ nodes — pruning
  only bites at the depth where the budget boundary sits.

Both failure modes point at the same missing chapter: **charging the
ε-side inside the sweep** (the pure-ε stratum of a safe class costs
q × the quotient-code coset weight of the class's ε-image — the
safe-floor problem *recurses into the ε-quotient BB code*, cf. §6).

### 5.4 Negative control

by90 ([[90,8,8]] Z₃₀×Z₃, the Bravyi-360 tower bottom whose x-rung
freezes at 12 < 16, A14 §13): the engine must NOT certify floor 16.
Disposition: no linked fiber (frame scoreboard in
`data/a29/extra_targets.jsonl`); the unpaired q=5/q=3 DFS attempts
burned caps for ~50 min without any certificate and the run was
killed — the required outcome (**no false certificate on a known-false
floor**) is met by construction and by the run; a witness-producing
refutation at this size awaits the ε-recursion/orbit-reduction
upgrades. (An engine-refuted known-false floor IS on record at small
size: pair72 target-10, §4.)

## 6. Limits, residue, and the follow-on program

1. **Scope (inherited, unchanged):** abelian G; |A|, |B| odd (parity);
   (R) as a checked hypothesis; safe-class semantics via Prop A14.1.
2. **Engine v1 decides** codes with a linked fiber whose paired sweep
   fits (f2a6-class, the new Z₅×Z₁₅ cells), and unlinked codes whose
   unpaired DFS stays small (pair72, the A8 base at q=7). It refuses —
   with named caps, never false claims — when the ε-freedom or the
   node tree blows up.
3. **The ε-recursion chapter (the one genuinely new mathematics
   left):** per L5, a δ-config's exact cost = relaxed + a weighted
   coset-min in the ε-quotient pair map; at δ = 0 this is exactly
   q × (safe-floor problem of the ε-quotient code). Charging partial
   ε-costs during the DFS (per-unit lower bounds conditioned on the
   quotient's own difference-set/small-cycle structure, or recursing
   the engine into the quotient) is the designed route to Z₂₁×Z₃,
   Z₁₅×Z₆-floor-20, and beyond — A27's "constrained-h taxonomy",
   now with a precise interface.
4. **Translation-orbit mask reduction** (÷|G| on leaves/nodes) is
   designed but unimplemented — a constant-factor ~50–90× on the DFS,
   possibly enough for the Z₁₅×Z₆ cells by brute force; the
   ε-recursion is the principled route.
5. **Gross stays out of reach by budget arithmetic** (S = 12 vs
   b = 22 on its only fiber) — consistent with its slot-frame history;
   the engine is a complement, not a replacement, for that machinery.

## 7. Lean-feasibility assessment

The certificate mass of an engine run is exactly the A22/A23 shape the
repo has already shipped twice: per-leaf rank/pivot certificates
(hundreds, cf. A23's 429), the residue/weight lemmas (L1–L3 are
32-case-style decides), seamC dictionary facts (per class, as in
`SeamReduction.lean`), and the K₀/quotient bookkeeping. Three deltas
vs the shipped instances: (i) leaves come from the state-DFS, so the
Lean statement should quantify over state assignments (the completeness
lemma L4 replaces the `table_coverage` Finset enumeration — same
`native_decide`-enumerable shape); (ii) τ-states add a 3-way case
split per twisted site; (iii) the exact-ε pass on suspects becomes a
per-suspect linear-functional minimum certificate (argmin flip + the
Σ-identity — a decide per suspect). For the A8 base the suspects list
is EMPTY (relaxed alone certifies), so its Lean packaging is *simpler*
than A23's: data + sweeps only. Estimated at 1–2 sessions for the A8
instance on the established pattern; the generic engine-emitter (one
generator serving all instances, in the GENERATORS.md discipline) is
the natural S-track follow-on.

## 8. Position in the program narrative (the teaching doc)

`docs/teaching/bb-doubling-theorem.pdf` presents the doubling theorem
as four conditions on the base; read as a recipe, the safe floor is
its cost center, and the doc's honest footnote was "per new code this
is solver-hours with no proof artifact." A29 replaces that footnote:
on engine-reachable codes the safe floor is a uniform seconds-scale
certified procedure. Two narrative corrections flow back into any
future edition: (i) the gross-specific F₄ rigidity the doc must caveat
is *provably incidental for the safe floor* (the A8 §5.2 closure never
touches it — the invariant content of the weight tables is a residue
map and popcount); (ii) the "how general is this?" outlook question
(which A27 answered with a tiering) now has a constructive middle:
Tier 2's "portable recipe, per-code execution" is mechanized, and the
per-code residue is localized in one named object — the ε-quotient
recursion. The loop exposition → gap (A27 P5) → repair (τ-states) →
new theorems (docket cells, A8 core) is itself worth a paragraph in
the teaching doc's closing section.

## Appendix: verification map

| claim | check |
|---|---|
| engine lemmas L1–L5 on live data | frame-construction asserts + `chain_weight_check` in validation scripts |
| pair72 / f2a6 rows of §4 | scratch validation scripts (this session; reproduce via `bb_lab.fibering` calls shown in §4) |
| docket runs | `scripts/a29_docket_engine.py` → `data/a29/docket_engine.jsonl` |
| extra targets | `scripts/a29_extra_targets.py` → `data/a29/extra_targets.jsonl` |
