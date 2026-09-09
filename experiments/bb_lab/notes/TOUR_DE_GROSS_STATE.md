# The tour-de-gross distance conjecture: state of the program

Cross-thread synthesis (A40 sessions 1–9 + A42 sessions 1–2; A42
session 3 in flight), maintained by the orchestrating session.  Every
claim below is banked in `A40_tour_de_gross_column.md` (§ refs) or
`A42_spectral_floor.md` with its tier; this file adds no new claims.
Last updated: 2026-08-29 (post S9).

> **2026-09-08 audit update.** The all-rectangle, gap-dense L3 proposed
> in A43 is refuted: `(54,6)` has a verified nontrivial weight-48 example
> with gaps `(3,1)`, whereas L3 demands 54. It still lifts to a cylinder.
> Thus the gap-defined and no-cylinder-lift sectors must be distinguished;
> the assembly below is a conditional historical proposal and does not
> justify retiring the earlier comparison route. This is not a
> counterexample to an IBM member. Current proof tracks, sources, and an
> exact related comparison family are in `IBM_EXACT_DISTANCE_PROGRAM.md`.
>
> **Cylinder scope correction.** `A42_S6_codex_residual.md` audits the
> later original-lane certificates without transferring them through an
> unproved axis mirror. The rotation produces `(B,bar A)`, a different
> cylinder pair. Strict primary purity also differs from three-sheet
> purity when the odd part grows; order-127 and higher-three fibres need
> their stated qualifications. The historical inventory below must not
> be read as an unconditional all-orientation, all-period theorem.
>
> **Later analytic progress, 2026-09-08.** The original and mirror
> canonical coset leaders are now proved unique for every `3 ∤ q`
> (`A42_S7_spatial_equality.md`, `A42_MIRROR_CANONICAL.md`), closing
> the fixed-support small-hole strata, not full HM. S8 certifies a
> fixed generator strip and refutes its x-periodic relaxation. The
> all-size S11 parity law and three-cover kernels, including the
> order-127 exception, are proved and independently reviewed in
> `A40_PARITY_TOR_PROOF.md`. These classify homology, not its metric.
> `GROSS_SURFACE_PRODUCT.md` separately proves the exact LDPC family
> `[[144(3t^2-3t+1),12,12t]]`; it changes the BB architecture and is
> not an IBM-conjecture proof. The current scope and falsified
> shortcuts are summarized in `IBM_EXACT_DISTANCE_PROGRAM.md` §7.
>
> **Further scope checks, 2026-09-08.** That program's §8 records
> safe but restricted tail elimination, a weight-increasing pure
> normalization counterfamily, an obstruction to fixed-width linear
> HM recovery, and exact diagonal-invariant sector minima. None
> closes the unrestricted metric problem. The new all-b=1
> single-block class theorem is `A43_S4_block_purity.md`; even its
> class-matched colon relaxation is only six at gross, as proved in
> `IBM_RELATIVE_COLON_OBSTRUCTION.md`. The original conjecture
> remains open. These current qualifications supersede any stronger
> reading of the historical assembly below.
>
> **Exact local gap, 2026-09-08.** `IBM_CANONICAL_EXACT_GAP.md`
> proves that the original dyadic canonical affine problem has
> unique minimum q and exact next weight q+4 for every q>=2.
> This closes its fixed-support hole strata through six, not
> unrestricted HM or IBM distance. The original pure-cylinder
> Cartier maps and their metric obstructions are recorded in
> `A42_S15_Cartier_descent.md`. Current details are in
> `IBM_EXACT_DISTANCE_PROGRAM.md` §§10–11.

> **Exact modified-family theorem, 2026-09-08.** Keeping the original
> cleared P and changing only Q to (1+x+x^2)(y+x+1) now gives
> the exact [[72r^2,12,4r]] family for every r>=1 on the original
> square frames. Check weight becomes eight, while the ordinary
> primary algebra, including its tangent component, is retained.
> The full all-class, both-sector proof and matching stripe are in
> `IBM_MODIFIED_EXACT_DISTANCE_THEOREM.md`; the shorter analytic
> route uses two physical seam cuts, not full HM. Original IBM
> distance remains open. Current modified-family results and their
> exact scope supersede older residual language; see the program's
> §13. No claim transfers the new code's metric to the old checks.

## The object

IBM's family (arXiv:2506.03094): BB codes on Z_ℓ × Z_m, ℓ = 6(r+b),
m = 6r, b ∈ {0,1}, fixed pair A = 1+y+x³y⁻¹, B = 1+x+x⁻¹y⁻³;
stated n = 72r(r+b), k = 12; conjectured d = 6(2r+b−1).  At b = 1:
d = 12r = 2m.

## Theorem inventory (with tiers)

| statement | tier | where |
|---|---|---|
| d ≤ 12r, whole b=1 column (UB(r,1)) | closed-form witness + dual, kernel-checked r ≤ 6, hand-grade ∀r | A40 §7.1 |
| d ≥ ⌈3r/2⌉, both columns (L1) | hand proof + machine gates | A40 §8.1 |
| d((ℓ,6)) = 12 ∀ 6\|ℓ ≥ 12 (B6) | certificate, zero residue | A40 §9.6 |
| Exact members: (1,0) 6, (1,1) 12, (2,0) 18, (2,1) 24 | Lean / certificate | repo-wide |
| **d([[864,12]]) ≥ 12** (member (3,1); both sectors 12, independent engines) | certificate (scope-listed) | A40 §13–14 |
| V(A,B) = transverse F₄ ⊕ tangency F₄ (mult 2) ⊕ F₁₂₈ (ord 127) | certificate | A42 S0 |
| k-law: b=1 k = 12 ∀r (127-blocked); b=0 k = 26 at 127 \| r | theorem over the variety table | A42 S0 |
| Parity lemma (all X-cycle weights even) | certificate | A42 |
| Compact-cylinder floors: floor(p) = 2p at p = 3, 6, 9; ≥ 18 at p = 12 (= 24 over the h-DP envelope) | certificate (p=12 corner open) | A42 S1–S2 |
| Theorem H (Tor-purity): fixed 6-class inventory ∀p | theorem + machine checks to p = 24 | A42 |
| **Theorem W (b=1 windowed branch)**: wt ≥ floor_cyl(6r); r=1 ≥ 12 unconditional, r=2 ≥ 18; ∀r = 12r modulo (L-pure, L-band) | assembled, conditionality exact | A42 S2 |
| Momentum-budget machinery: drift calculus, interface tax (τ = 2.0–2.75), pinch, μ-echo Lemma R0, slip caps tight | certificate + 3-line hand lemma | A40 §11–14 |
| Species classification p ≤ 8 (TC63, W7 only; chiral; never b1-close) | certificate | A40 §10 |

## The assembly (2026-09-02, A43 §1: how d(C_{r,b}) = 6(2r+b−1) falls)

The conjectured value is ℓ + m − 6 in BOTH rows (ℓ = 6(r+b), m = 6r).
Lemma K splits any nontrivial logical by its cyclic gaps (threshold 4 =
the check spread):

1. **x-gap ≥ 4** (x-windowed): lifts to the period-m cylinder;
   **Theorem W** (A42) gives ≥ floor_cyl(m) = 2m.
2. **y-gap ≥ 4** (y-windowed): the mirror, ≥ floor_cyl′(ℓ) = 2ℓ.
3. **gap-dense** (doubly spanning, S11's D sector): **Lemma L3**
   (A43 §2, the crossing lemma, b-blind): ≥ ℓ + m − 6.
4. both gaps: trivial (Lemma G/K).

Hence d ≥ min(2m, 2ℓ, ℓ + m − 6) = ℓ + m − 6 exactly when ℓ − m ≤ 6,
i.e. b ∈ {0, 1} — the ONLY place b enters.  Upper bounds: UB(r,1) ≤ 12r
(A40 §5, ∀r); UB(r,0) ≤ 12r − 6 open beyond r ≤ 2 (A44).  Only the
shorter period needs Theorem W sharp; the x-period needs 2p − 6.
The comparison theorem / surgery program (A40 §16) and the −C removal
staging of the momentum lane are superseded by this assembly (A43 §4);
the momentum floor 3.625r (A40 §15) stays as the unconditional ∀r bound.

## The two open lemmas (the whole remaining gap)

- **Theorem W's residual** (windowed): (R1)+(R2) = the seam-doubling
  statement (A42 §2.17.4), (HM₉) when 9 | p (A42 §2.17.5); r ≤ 2
  certified; A42-S6.
- **L3, the crossing lemma** (toroidal): gap-dense nontrivial logicals
  weigh ≥ ℓ + m − 6.  Consistent with every banked population
  (A43 §5: min_D = ℓ + m − 6 at (6,6), (12,12), gross; ≥ 26 at
  (18,12) by A42 §2.17.2); support-level counting refuted, joint
  syndrome-level accounting required; covering/LP arguments capped at
  n/d_Z = ℓ (A43 §3); A43-S2.

## Falsified-claims discipline

Eleven major refutations are first-class in the ledgers (§9.8, §10.8,
§11.7, §12, §13.6, §14.5; A42 §3.6–3.8), including two retractions of
our own prior claims caught by internal audits (the S4 x_winds
classifier; the S1h register blindness).  Nothing outside the ledgers
is consumed by any standing claim.

## Guards

Literature: LLSC 2503.04699 anticipates the k-aggregates (our k-law is
the per-point refinement); per-spectral-point weight floors have no
quantum antecedent found; the A31 novelty pass is owed before any
external claim.  The paper's family passage is double-sourced
internally; a full-text re-extraction remains a minor open item.
