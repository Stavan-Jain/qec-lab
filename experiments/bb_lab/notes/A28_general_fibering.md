# A28 — The general fibering engine: portable safe-floor certification for BB codes

> **STATUS.** Session 1 (2026-08-06, branch `claude/fibering-bb-codes-5c9ff3`).
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

(filled after the runs)

## 6. Limits, residue, and the follow-on program

(filled at session end)

## 7. Lean-feasibility assessment

(filled at session end)

## Appendix: verification map

| claim | check |
|---|---|
| engine lemmas L1–L5 on live data | frame-construction asserts + `chain_weight_check` in validation scripts |
| pair72 / f2a6 rows of §4 | scratch validation scripts (this session; reproduce via `bb_lab.fibering` calls shown in §4) |
| docket runs | `scripts/a28_docket_engine.py` → `data/a28/docket_engine.jsonl` |
| extra targets | `scripts/a28_extra_targets.py` → `data/a28/extra_targets.jsonl` |
