/-
# The target-floor assembly for BB covers (deficit rungs: `m ≤ 2·d`)

The tower slice calculus (qec-lab A32/A33/A36) certifies covers whose
distance sits *below* the doubling value — the first such instance is
`[[288,12,18]]` over the gross code (`18 < 2·12`: an (R) rung with
deficit 6; `Codes/BivariateBicycle/Bb288/`).  The `BBDoubling` sector
dichotomy is unchanged: the two sector floors enter at the target `m`,
and the `push v = 0` (τ-diagonal) branch still delivers `2·d`, which
suffices whenever `m ≤ 2·d`.  At `m = 2·d` these specialize to the
`_double_of_logicalFloor` assemblies of `BBDoubling.lean`; the witness
is an arbitrary weight-`m` nontrivial chain rather than a diagonal
pullback (for a strict deficit no `τ`-lift attains `m`).

Lives in its own module (importing `BBDoubling`) so instance builds of
the doubling layer are not invalidated by the extension.
-/

import QEC.Stabilizer.Framework.Homological.BBDoubling

namespace Quantum
namespace Stabilizer
namespace Homological
namespace BB

open scoped BigOperators

namespace XDoubleCoverData

variable {G H : Type}
  [Fintype G] [AddCommGroup G] [DecidableEq G]
  [Fintype H] [AddCommGroup H] [DecidableEq H]
  (D : XDoubleCoverData G H)

/-- **Sector-dichotomy assembly at a target `m ≤ 2·d`**: given the base
logical floor at `d` and the two sector floors at `m`, every nontrivial
cover cycle has weight ≥ `m`. -/
theorem chainWeight_ge_of_sector_floors {d m : ℕ}
    (hbase : D.LogicalFloor d)
    (hM : D.DangerousFloorNZ m) (hS : D.SafeFloor m)
    (hm2d : m ≤ 2 * d) :
    ∀ v : G × Fin 2 → ZMod 2,
      v ∈ D.coverComplex.cycles → v ∉ D.coverComplex.boundaries →
      m ≤ D.coverComplex.chainWeight v := by
  intro v hv hnb
  by_cases hb : D.push1 v ∈ D.baseComplex.boundaries
  · by_cases h0 : D.push1 v = 0
    · exact le_trans hm2d
        (D.dangerous_zero_rung_of_logicalFloor hbase hv hnb h0)
    · exact hM v hv hnb hb h0
  · exact hS v hv hb

/-- **Chain-level `d(cover) = m` at a deficit target**: `m` is attained
(by the explicit witness) and minimal (by the sector floors). -/
theorem chain_distance_eq_of_sector_floors {d m : ℕ}
    (hbase : D.LogicalFloor d)
    (hM : D.DangerousFloorNZ m) (hS : D.SafeFloor m)
    (hm2d : m ≤ 2 * d)
    (vStar : G × Fin 2 → ZMod 2)
    (hv_cyc : vStar ∈ D.coverComplex.cycles)
    (hv_nb : vStar ∉ D.coverComplex.boundaries)
    (hv_w : D.coverComplex.chainWeight vStar = m) :
    IsLeast {w : ℕ | ∃ v : G × Fin 2 → ZMod 2,
      v ∈ D.coverComplex.cycles ∧ v ∉ D.coverComplex.boundaries ∧
      D.coverComplex.chainWeight v = w} m := by
  constructor
  · exact ⟨vStar, hv_cyc, hv_nb, hv_w⟩
  · rintro w ⟨v, hv, hnb, rfl⟩
    exact D.chainWeight_ge_of_sector_floors hbase hM hS hm2d v hv hnb

/-- Dual-side mirror of the target-floor sector bound, via the
chain-level `d_X = d_Z` duality. -/
theorem dual_chainWeight_ge_of_sector_floors {d m : ℕ}
    (hbase : D.LogicalFloor d)
    (hM : D.DangerousFloorNZ m) (hS : D.SafeFloor m)
    (hm2d : m ≤ 2 * d) :
    ∀ c ∈ D.coverComplex.dualCycles, c ∉ D.coverComplex.dualBoundaries →
      m ≤ D.coverComplex.chainWeight c := by
  have hX : ∀ c ∈ (bbChainComplex D.Ac D.Bc).cycles,
      c ∉ (bbChainComplex D.Ac D.Bc).boundaries →
      m ≤ (bbChainComplex D.Ac D.Bc).chainWeight c := fun c hc hnb =>
    D.chainWeight_ge_of_sector_floors hbase hM hS hm2d c hc hnb
  exact (bb_cycle_bound_iff_dual_bound D.Ac D.Bc m).mp hX

/-- Pauli-level lower bound at the target: every nontrivial logical
operator of the cover's homological stabilizer group has weight ≥ `m`. -/
theorem logical_weight_ge_of_sector_floors {d m : ℕ}
    (hbase : D.LogicalFloor d)
    (hM : D.DangerousFloorNZ m) (hS : D.SafeFloor m)
    (hm2d : m ≤ 2 * d)
    (g : NQubitPauliGroupElement D.coverComplex.numQubits)
    (hg : Quantum.StabilizerGroup.IsNontrivialLogicalOperator g
      D.coverComplex.homologicalStabilizerGroup) :
    m ≤ NQubitPauliGroupElement.weight g :=
  HomologicalCode.chainWeight_lower_bound_transfers D.coverComplex m
    (fun c hc hnb =>
      D.chainWeight_ge_of_sector_floors hbase hM hS hm2d c hc hnb)
    (D.dual_chainWeight_ge_of_sector_floors hbase hM hS hm2d) g hg

/-- **Pauli-level `d(cover) = m` at a deficit target**: given the base
logical floor at `d`, the two sector floors at `m ≤ 2·d`, and an
explicit weight-`m` nontrivial witness chain, `m` is the least weight of
a nontrivial logical operator of the cover's homological stabilizer
group. -/
theorem pauli_distance_eq_of_sector_floors {d m : ℕ}
    (hbase : D.LogicalFloor d)
    (hM : D.DangerousFloorNZ m) (hS : D.SafeFloor m)
    (hm2d : m ≤ 2 * d)
    (vStar : G × Fin 2 → ZMod 2)
    (hv_cyc : vStar ∈ D.coverComplex.cycles)
    (hv_nb : vStar ∉ D.coverComplex.boundaries)
    (hv_w : D.coverComplex.chainWeight vStar = m) :
    IsLeast {w : ℕ | ∃ g : NQubitPauliGroupElement D.coverComplex.numQubits,
      Quantum.StabilizerGroup.IsNontrivialLogicalOperator g
        D.coverComplex.homologicalStabilizerGroup ∧
      NQubitPauliGroupElement.weight g = w} m := by
  constructor
  · refine ⟨D.coverComplex.chainXOperator vStar, ?_, ?_⟩
    · exact (HomologicalCode.chainXOperator_isNontrivialLogical_iff
        (X := D.coverComplex) vStar).mpr ⟨hv_cyc, hv_nb⟩
    · rw [HomologicalCode.weight_chainXOperator, hv_w]
  · rintro w ⟨g, hg, rfl⟩
    exact D.logical_weight_ge_of_sector_floors hbase hM hS hm2d g hg

end XDoubleCoverData

end BB
end Homological
end Stabilizer
end Quantum
