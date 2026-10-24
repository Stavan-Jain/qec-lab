"""Emit certified orbit metadata with the authored counting-theory tail."""
import os
from pathlib import Path
STAGE=Path(os.environ["QEC360_STAGE"]).resolve()
from pathlib import Path
p = Path(str(STAGE / 'OrbitBatchMetadata.lean'))
old = p.read_text()
allchecked = old[old.index('/-- Every stored batch satisfies'):old.index('\nend OrbitBatchData')]
tail = Path(str(STAGE / 'OrbitMetadataTail.lean.in')).read_text()
s = 'import OrbitBatchData\nimport OrbitBatchCheck\nimport OrbitBatchConcreteSoundness\nimport FiniteSumExceptions\nimport Mathlib.Tactic.FinCases\nset_option maxRecDepth 100000\nset_option maxHeartbeats 0\nnamespace OrbitBatchData\nopen PackedDigits\n\n/-- The fifteen nonzero representatives with a three-element translation stabilizer. -/\ndef fixedThree : Finset Nat := {299,963,971,2941,26949,26959,58329,70679,110969,123726,136237,139734,150152,153236,268422}\n\n/-- The certified stabilizer count at a global census index. -/\ndef fixedAt (i : Nat) : Nat := if i = 0 then 45 else if i = 70587 then 9\n  else if i ∈ fixedThree then 3 else 1\n\n/-- The thirteen batches containing a nontrivial stabilizer. -/\ndef specialBatches : Finset Nat := {0,2,26,56,68,69,108,120,133,136,146,149,262}\n\n/-- All other batches consist entirely of representatives with trivial stabilizers. -/\ntheorem fixedAt_regular (b j : Nat) (hb : b ∉ specialBatches) (hj : j < 1024) :\n    fixedAt (1024*b+j) = 1 := by\n  simp only [specialBatches, Finset.mem_insert, Finset.mem_singleton, not_or] at hb\n  simp only [fixedAt, fixedThree, Finset.mem_insert, Finset.mem_singleton]\n  split_ifs <;> omega\n\n/-- Every valid batch uses the proved 128-bit repetition constant. -/\ntheorem repeat_size_eq (b : Fin 268) :\n    repeat128 (sizeAt b.val) = packF 128 (sizeAt b.val) (fun _ => 1) := by\n  by_cases hb : b.val < 267\n  · simpa [sizeAt,hb] using repeat1024_eq\n  · have he : b.val = 267 := by omega\n    simpa [sizeAt,he] using repeat603_eq\n\n/-- Every stabilizer count fits within one packed lane. -/\ntheorem fixedAt_lt (i : Nat) : fixedAt i < 2^128 := by\n  unfold fixedAt\n  split_ifs <;> decide\n\n'
special = [0, 2, 26, 56, 68, 69, 108, 120, 133, 136, 146, 149, 262]
for b in special:
    s += f'/-- The exceptional fixed counts in batch {b} match their packed certificate. -/\ntheorem expected_special{b} : expected{b} = packF 128 1024 (fun j => fixedAt (1024*{b}+j)) := by\n  decide +kernel\n\n'
s += '/-- Every ordinary expected-count word is the repetition word. -/\ntheorem expected_regular (b : Fin 268) (hb : b.val ∉ specialBatches) :\n    expectedAt b.val = repeat128 (sizeAt b.val) := by\n  fin_cases b\n'
for b in range(268):
    s += '  · exact False.elim (hb (by decide))\n' if b in special else '  · rfl\n'
s += '\n/-- The expected packed stabilizer counts agree with their sparse global description. -/\ntheorem expected_eq (b : Fin 268) :\n    expectedAt b.val = packF 128 (sizeAt b.val) (fun j => fixedAt (1024*b.val+j)) := by\n  by_cases hb : b.val ∈ specialBatches\n  · simp only [specialBatches, Finset.mem_insert, Finset.mem_singleton] at hb\n    rcases hb with ' + ' | '.join(['hb'] * len(special)) + '\n'
for b in special:
    s += f'    · have he : b = {b} := Fin.ext hb\n      subst b\n      exact expected_special{b}\n'
s += "  · rw [expected_regular b hb, repeat_size_eq]\n    apply packF_congr\n    intro j hj\n    have hj' : j < 1024 := by\n      unfold sizeAt at hj\n      split_ifs at hj <;> omega\n    exact (fixedAt_regular b.val j hb hj').symm\n\n"
for b in range(267):
    s += f'/-- Adjacent packed batches {b} and {b + 1} are strictly separated. -/\ntheorem cross_sorted{b} : lane 128 1023 words{b} < lane 128 0 words{b + 1} := by decide +kernel\n'
s += '\n/-- Consecutive batches are strictly separated in natural mask order. -/\ntheorem cross_sorted (b : Fin 267) :\n    lane 128 1023 (wordsAt b.val) < lane 128 0 (wordsAt (b.val+1)) := by\n  fin_cases b\n'
for b in range(267):
    s += f'  · exact cross_sorted{b}\n'
s += '\n' + allchecked
s += '/-- The accepted payload mask bounds every stored word by its allocated lanes. -/\ntheorem words_bounded (b : Fin 268) : wordsAt b.val < 2^(128*sizeAt b.val) := by\n  have h := all_checked b\n  simp only [OrbitBatchCheck.batchCheck, Bool.and_eq_true, beq_iff_eq] at h\n  have he := h.1.1\n  simp only [OrbitBatchCheck.extraCheck, Bool.and_eq_true, beq_iff_eq] at he\n  have hfit := he.1.1\n  rw [repeat_size_eq b, repeat_mul_eq_packF] at hfit\n  calc\n    wordsAt b.val = wordsAt b.val &&& packF 128 (sizeAt b.val) (fun _ => 2^90-1) := hfit.symm\n    _ ≤ packF 128 (sizeAt b.val) (fun _ => 2^90-1) := Nat.and_le_right\n    _ < 2^(128*sizeAt b.val) := packF_lt 128 (sizeAt b.val) _ (by decide)\n      (fun _ _ => show 2^90-1 < 2^128 from by decide)\n\n' + tail
a = s.index('theorem orbit_size_sum :')
b = s.index('\nend OrbitBatchData', a)
s = s[:a] + 'theorem orbit_size_sum : (∑ i : Fin 274011, 45/fixedAt i.val) = 12329961 := by\n  have he := exceptional_counts\n  have hsum := Quantum.Stabilizer.Homological.sum_of_exceptional_values\n    (fun i : Fin 274011 => 45/fixedAt i.val) exceptional 45 274011 17 231\n    (Fintype.card_fin 274011) he.1 he.2\n    (fun i hi => by change 45 / fixedAt i.val = 45; rw [fixedAt_of_not_exceptional i hi])\n  exact hsum.trans (by decide +kernel)\n' + s[b:]
p.write_text(s)
