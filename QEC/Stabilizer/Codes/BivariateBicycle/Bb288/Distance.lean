/-
# `d([[288,12,18]]) = 18` at the chain and Pauli levels — the two-tier claim

The first instance of the **target-floor assembly** (deficit rung:
`18 < 2·12`; `BBDoubling.lean` §target-floor), and the first tower
instance whose base leg is discharged by a kernel-checked library
theorem rather than a certificate: the base of this cover is the gross
code itself, so `coverData.LogicalFloor 12` follows from the
unconditional `gross_chain_distance_eq_12` (zero new obligations).

Following the `Z5Z15F2A6` precedent (the shipped two-tier shape), the
two sector floors enter as named hypotheses whose provenance is the A36
deterministic certificate layer
(`qec-lab:experiments/bb_lab/notes/A36_bb288_d18_certificate.md`,
data `qec-lab:experiments/bb_lab/data/a36/`) — no SAT anywhere:

* `hM : coverData.DangerousFloorNZ 18` — every nontrivial bb288 logical
  whose gross shadow is a nonzero stabilizer has weight ≥ 18.
  Certificate: the gross stabilizer census ≤ 16 (one 6-offset coset-BZ
  pass, 7.50e9 nodes, exact node-count invariants; 33,588 vectors in
  469 translation orbits) + 469 dangerous rungs ALL PASS at
  `M = (18−|b|)/2` (restricted-MITM lanes, complete by the
  exact-subset-sum argument, every candidate re-verified) + G-transport.
* `hS : coverData.SafeFloor 18` — every bb288 cycle whose gross shadow
  is not a stabilizer has weight ≥ 18.  Certificate: the five seam
  orbit-rep coset censuses ≤ 16 (same BZ pass; 395 elements, minima 12)
  + 395 feasibility seam rungs ALL PASS + G-transport.  (The *naive*
  seam floor at 18 is FALSE — the seam minima sit at `d(gross) = 12` —
  which is exactly why the floor is per-element and certificate-shaped.)

Given those, `bb288_chain_distance_eq_18` / `bb288_pauli_distance_eq_18`
deliver `d = 18`.  Unconditional here: the base floor (library), the
weight-18 witness (`Witness.lean`), the assembly, and hence the
membership halves of both `IsLeast` statements
(`bb288_exists_weight18_nontrivial_cycle`).
-/

import QEC.Stabilizer.Codes.BivariateBicycle.Bb288.Witness
import QEC.Stabilizer.Codes.BivariateBicycle.Gross.LayerInstance
import QEC.Stabilizer.Framework.Homological.BBTargetFloor

namespace Quantum
namespace Stabilizer
namespace Homological
namespace BB
namespace Bb288

-- The base-floor discharge is a defeq pass through the
-- `coverData.baseComplex = grossComplex` bridge (deep `Prod`/`ZMod`
-- instance chains, as in the sibling instance files).
set_option maxRecDepth 8192

/-! ## The base floor, from the library's gross theorem -/

/-- **`LogicalFloor 12` for the gross base, kernel-checked**: the lower
half of the unconditional `gross_chain_distance_eq_12`.  (The bundle's
base complex is definitionally `grossComplex`.) -/
theorem coverData_logicalFloor : coverData.LogicalFloor 12 := by
  intro u hcyc hnb
  exact gross_chain_distance_eq_12.2 ⟨u, hcyc, hnb, rfl⟩

/-! ## The capstones (conditional on the two A36 sector certificates) -/

/-- **Chain-level `d([[288,12,18]]) = 18`**, conditional on the two
certificate-checked sector floors (module docstring for provenance):
18 is the least weight of a nontrivial cycle of the bb288 complex. -/
theorem bb288_chain_distance_eq_18
    (hM : coverData.DangerousFloorNZ 18)
    (hS : coverData.SafeFloor 18) :
    IsLeast {w : ℕ | ∃ v : G288 × Fin 2 → ZMod 2,
      v ∈ coverData.coverComplex.cycles ∧
      v ∉ coverData.coverComplex.boundaries ∧
      coverData.coverComplex.chainWeight v = w} 18 :=
  coverData.chain_distance_eq_of_sector_floors (d := 12)
    coverData_logicalFloor hM hS (by norm_num)
    vStar288 vStar288_mem_cycles vStar288_not_mem_boundaries
    chainWeight_vStar288

/-- **Pauli-level `d([[288,12,18]]) = 18`**, conditional on the two
certificate-checked sector floors: 18 is the least weight of a
nontrivial logical operator of the bb288 homological stabilizer
group. -/
theorem bb288_pauli_distance_eq_18
    (hM : coverData.DangerousFloorNZ 18)
    (hS : coverData.SafeFloor 18) :
    IsLeast {w : ℕ | ∃ g : NQubitPauliGroupElement
        coverData.coverComplex.numQubits,
      Quantum.StabilizerGroup.IsNontrivialLogicalOperator g
        coverData.coverComplex.homologicalStabilizerGroup ∧
      NQubitPauliGroupElement.weight g = w} 18 :=
  coverData.pauli_distance_eq_of_sector_floors (d := 12)
    coverData_logicalFloor hM hS (by norm_num)
    vStar288 vStar288_mem_cycles vStar288_not_mem_boundaries
    chainWeight_vStar288

end Bb288
end BB
end Homological
end Stabilizer
end Quantum
