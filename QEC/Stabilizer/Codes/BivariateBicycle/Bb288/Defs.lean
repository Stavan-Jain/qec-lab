/-
# The [[288,12,18]] record code as a twisted cover of the gross code — definitions

`bb288 = [[288,12,18]]` over `Z₁₂ × Z₁₂` (`A = x³ + y² + y⁷`,
`B = y³ + x + x²`; arXiv:2308.07915 Table 3, the strongest published
solver-exact BB code) presented as a free ℤ₂ cover of the **gross code**
along `y` (deck `σ = y⁶`): its `y mod 6` descent is literally
`grossA/grossB` (`y⁷ ≡ y`).  Unlike every previous instance of the layer
this lift is **twisted** (the `y⁷` term carries the deck element) — the
`XDoubleCoverData` axioms are pushforward identities, which the twisted
lift satisfies term-by-term, so no new framework is needed (qec-lab A32
Lemma 1 "twist-generic", A35 layer map, A36 §1).

The base bundle entries are the *library's own* `grossA/grossB`, so
`coverData.baseComplex` is definitionally `grossComplex` and the
kernel-checked `d(gross) = 12` discharges the base floor
(`Distance.lean`).  Offline provenance and validation:
`qec-lab:experiments/bb_lab/scripts/a36_gen_lean_data.py` (pushforward
identities re-verified numerically; conventions anchored against the
kernel-checked Z5Z15F2A6 facts) and the A36 note
(`qec-lab:experiments/bb_lab/notes/A36_bb288_d18_certificate.md`).
-/

import QEC.Stabilizer.Codes.BivariateBicycle.Gross.Defs
import QEC.Stabilizer.Framework.Homological.BBDoubling

namespace Quantum
namespace Stabilizer
namespace Homological
namespace BB
namespace Bb288

-- Defeq checks through the `bbChainComplex` structure projections unfold
-- deep `Prod`/`ZMod` instance chains (same as the gross and Z5Z15F2A6
-- instantiations).
set_option maxRecDepth 4096

/-! ## The cover group and polynomials -/

/-- The cover group `Z₁₂ × Z₁₂` (free ℤ₂ cover of `GrossGroup = Z₁₂ × Z₆`
doubling `y`). -/
abbrev G288 : Type := ZMod 12 × ZMod 12

/-- Cover `A = x³ + y² + y⁷` — the **twisted** lift of
`grossA = x³ + y + y²` (the `y⁷` term carries the deck element `y⁶`). -/
def a288 : G288 → ZMod 2 := fun g =>
  if g = (3, 0) ∨ g = (0, 2) ∨ g = (0, 7) then 1 else 0

/-- Cover `B = y³ + x + x²` (literal lift of `grossB`). -/
def b288 : G288 → ZMod 2 := fun g =>
  if g = (0, 3) ∨ g = (1, 0) ∨ g = (2, 0) then 1 else 0

/-! ## The chain complex -/

/-- The `[[288,12,18]]` chain complex. -/
noncomputable def bb288Complex : HomologicalCode := bbChainComplex a288 b288

theorem bb288Complex_numQubits : bb288Complex.numQubits = 288 := by
  change bbNumQubits G288 = 288
  unfold bbNumQubits
  norm_num [Fintype.card_prod, ZMod.card]

/-! ## The covering projection and its finite obligations

All four bundle obligations are kernel `decide`s (no `native_decide` in
this instance): the heartbeat bumps below pay for reducing the
`ZMod.castHom` coercion chain point-by-point in the elaborator and the
kernel. -/

/-- The covering projection `Z₁₂×Z₁₂ →+ Z₁₂×Z₆` (reduce `y` mod 6). -/
def projB288 : G288 →+ GrossGroup :=
  AddMonoidHom.prodMap
    (AddMonoidHom.id (ZMod 12))
    (ZMod.castHom (by norm_num : (6 : ℕ) ∣ 12) (ZMod 6)).toAddMonoidHom

set_option maxHeartbeats 2000000 in
-- 144² fiber pairs, each reducing the `castHom` coercion chain in the
-- kernel (kernel `decide`; ~10× the default budget).
theorem projB288_fiber :
    ∀ g g' : G288, projB288 g' = projB288 g ↔ g' = g ∨ g' = g + (0, 6) := by
  decide

set_option maxHeartbeats 800000 in
-- 72 section points through the same coercion chain.
theorem projB288_sec :
    ∀ p : GrossGroup, projB288 (p.1, (p.2.val : ZMod 12)) = p := by
  decide

set_option maxHeartbeats 2000000 in
-- 72 base points × 144-point fiber sums (the twisted `y⁷` term folds to
-- `y`: the pushforward identity is exactly where the twist dies).
theorem projB288_push_A : fiberSumFn ⇑projB288 a288 = grossA := by
  decide

set_option maxHeartbeats 2000000 in
-- Same shape for `B` (a literal lift).
theorem projB288_push_B : fiberSumFn ⇑projB288 b288 = grossB := by
  decide

/-! ## The cover bundle -/

/-- The parametric cover data: projection `Z₁₂×Z₁₂ →+ Z₁₂×Z₆` (reduce `y`
mod 6), deck `y⁶ = (0,6)`, canonical section, the twisted cover
polynomials, and the library's gross polynomials as the base pair. -/
def coverData : XDoubleCoverData G288 GrossGroup where
  proj := projB288
  deckS := (0, 6)
  sec := fun p => (p.1, (p.2.val : ZMod 12))
  Ac := a288
  Bc := b288
  Ab := grossA
  Bb := grossB
  deckS_ne_zero := by decide
  proj_fiber := projB288_fiber
  proj_sec := projB288_sec
  push_A := projB288_push_A
  push_B := projB288_push_B

lemma coverData_coverComplex : coverData.coverComplex = bb288Complex := rfl

lemma coverData_baseComplex : coverData.baseComplex = grossComplex := rfl

end Bb288
end BB
end Homological
end Stabilizer
end Quantum
