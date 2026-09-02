# The tour-de-gross distance conjecture: state of the program

Cross-thread synthesis (A40 sessions 1–9 + A42 sessions 1–2; A42
session 3 in flight), maintained by the orchestrating session.  Every
claim below is banked in `A40_tour_de_gross_column.md` (§ refs) or
`A42_spectral_floor.md` with its tier; this file adds no new claims.
Last updated: 2026-08-29 (post S9).

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

## The assembly (how d(C_{r,1}) = 12r falls)

Lemma K splits any class-minimal nontrivial logical:

1. **Windowed branch** (an x- or y-gap ≥ 4): lifts to a cylinder;
   **Theorem W** gives ≥ 12r, today modulo **L-band mixed** (+ the
   p = 12 corner, in flight as A42 S3's register-quotient racer).
2. **Doubly-spanning (toroidal) branch**: the momentum budget.
   Per-member DP floors (now 12 at (3,1)); the named ∀r route is the
   **+3q/slab ladder** (J′ = 3q/slab rate conjecture ⟹ analytic
   ~7.5r), with the remaining local piece being the **financed-half
   constant** (net anchor rise ≤ 3 + #K1-inputs; measured+hand).

Then d ≥ 12r − C meets UB(r,1) ≤ 12r, and the −C removal strategy is
staged: discount-event classification killed by the gcd(ℓ,m) = 6
arithmetic (the b-bit, four independent appearances), then parity
(C < 2 suffices) or stratum congruences (C < 6) snap the bound.

## The two open lemmas (the whole remaining gap)

- **L-band mixed half** (windowed): the slot-count → 2p scaling
  mechanism for mixed classes; pure half closed (β-lemma).  Finite
  first instance: p = 12, weights {20, 22}.
- **Financed-half constant + ladder rate** (toroidal): extend μ-echo
  bookkeeping to count K1 inputs; prove J′ = 3q/slab.

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
