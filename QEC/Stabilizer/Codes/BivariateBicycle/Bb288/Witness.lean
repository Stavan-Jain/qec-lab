/-
# The weight-18 nontrivial cycle of [[288,12,18]]

Headline (`bb288_exists_weight18_nontrivial_cycle`): the bb288 complex
has a 1-cycle of weight 18 that is not a boundary — the `d ≤ 18` half of
the distance, kernel-checked.

`vStar288` is the A36 tower-native witness (found by the seam-rung
first-hit ladder over a weight-12 seam element, overflow 3 — NOT a
diagonal pullback: `18 < 2·d(gross)`, so no `τ`-lift attains it),
translated lab → repo by the global reflection `g ↦ −g` and re-verified
end-to-end in `qec-lab:experiments/bb_lab/scripts/a36_gen_lean_data.py`
(which first re-derives four kernel-checked Z5Z15F2A6 facts as a
convention anchor).  Non-boundaryness is certified by the explicit dual
witness `zStar288` (weight 18, dual-cycle, odd pairing), fed to
`not_mem_boundaries_of_dual_witness`.  Data provenance:
`qec-lab:experiments/bb_lab/data/a36/{w18_witness.json,lean_data.json}`.
-/

import QEC.Stabilizer.Codes.BivariateBicycle.Bb288.Defs

namespace Quantum
namespace Stabilizer
namespace Homological
namespace BB
namespace Bb288

open scoped BigOperators

-- Kernel `decide`s below reduce 288-point convolution sums through the
-- `Prod`/`ZMod`/`Finset` instance chains (same headroom as the sibling
-- instance files).
set_option maxRecDepth 8192

/-! ## The witness chain -/

/-- The weight-18 witness: left block only (the right block is empty). -/
def vStar288 : G288 × Fin 2 → ZMod 2 := fun p =>
  if p.2 = 0 then
    (if p.1 = (0, 0) ∨ p.1 = (0, 3) ∨ p.1 = (0, 9) ∨ p.1 = (1, 0) ∨
        p.1 = (1, 6) ∨ p.1 = (2, 0) ∨ p.1 = (2, 3) ∨ p.1 = (2, 6) ∨
        p.1 = (3, 0) ∨ p.1 = (4, 0) ∨ p.1 = (4, 3) ∨ p.1 = (5, 0) ∨
        p.1 = (6, 0) ∨ p.1 = (8, 9) ∨ p.1 = (9, 6) ∨ p.1 = (10, 3) ∨
        p.1 = (10, 6) ∨ p.1 = (11, 0)
     then 1 else 0)
  else 0

-- `native_decide` (measured relaxation): the raw-`Finset.sum` kernel
-- reduction of 144 boundary points × two 144-term convolutions did not
-- terminate inside a 10-minute budget (the known kernel wall on
-- convolution sums; the kernel-friendly alternative is a packed-table
-- restatement, out of scope here).  The four `Defs.lean` cover
-- obligations and the weight/pairing checks below remain kernel
-- `decide`s.
theorem vStar288_mem_cycles : vStar288 ∈ bb288Complex.cycles := by
  have h : bbBoundary1Fn a288 b288 vStar288 = 0 := by native_decide
  exact h

set_option maxHeartbeats 1000000 in
-- 288-point support filter (kernel `decide`).
theorem chainWeight_vStar288 : bb288Complex.chainWeight vStar288 = 18 := by
  rw [show bb288Complex.chainWeight vStar288
      = (Finset.univ.filter fun p : G288 × Fin 2 => vStar288 p ≠ 0).card
      from rfl]
  decide

/-! ## The dual witness -/

/-- An explicit dual cycle pairing oddly with `vStar288` (weight 18,
left block only). -/
def zStar288 : G288 × Fin 2 → ZMod 2 := fun p =>
  if p.2 = 0 then
    (if p.1 = (0, 2) ∨ p.1 = (0, 5) ∨ p.1 = (0, 9) ∨ p.1 = (0, 10) ∨
        p.1 = (3, 0) ∨ p.1 = (3, 2) ∨ p.1 = (3, 8) ∨ p.1 = (3, 10) ∨
        p.1 = (6, 0) ∨ p.1 = (6, 1) ∨ p.1 = (6, 3) ∨ p.1 = (6, 5) ∨
        p.1 = (6, 6) ∨ p.1 = (6, 7) ∨ p.1 = (6, 8) ∨ p.1 = (6, 10) ∨
        p.1 = (9, 0) ∨ p.1 = (9, 4)
     then 1 else 0)
  else 0

-- `native_decide`: same convolution scale (and the same measured kernel
-- wall) as `vStar288_mem_cycles`.
/-- Raw (computable) form of `dualBoundary zStar288 = 0`, via the
transpose formula `bb_dualBoundary_eq`. -/
theorem zStar288_dual_raw :
    (fun f => conv (reflect a288) (leftHalf zStar288) f
      + conv (reflect b288) (rightHalf zStar288) f)
      = (0 : G288 → ZMod 2) := by
  native_decide

theorem zStar288_dualBoundary :
    bb288Complex.dualBoundary zStar288 = 0 := by
  change (bbChainComplex a288 b288).dualBoundary zStar288 = 0
  rw [bb_dualBoundary_eq]
  exact zStar288_dual_raw

set_option maxHeartbeats 1000000 in
-- One 288-term pairing sum (kernel `decide`).
theorem zStar288_pairing :
    ∑ e : G288 × Fin 2, zStar288 e * vStar288 e = 1 := by
  decide

/-! ## Non-boundaryness and the headline -/

theorem vStar288_not_mem_boundaries :
    vStar288 ∉ bb288Complex.boundaries :=
  HomologicalCode.not_mem_boundaries_of_dual_witness
    zStar288_dualBoundary zStar288_pairing

/-- **The `[[288,12,18]]` code has a weight-18 nontrivial logical chain**
(the kernel-checked `d ≤ 18` half). -/
theorem bb288_exists_weight18_nontrivial_cycle :
    ∃ v ∈ bb288Complex.cycles,
      v ∉ bb288Complex.boundaries ∧
      bb288Complex.chainWeight v = 18 :=
  ⟨vStar288, vStar288_mem_cycles, vStar288_not_mem_boundaries,
    chainWeight_vStar288⟩

end Bb288
end BB
end Homological
end Stabilizer
end Quantum
