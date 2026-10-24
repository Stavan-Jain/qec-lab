"""Emit kernel-checkable certificates into QEC360_STAGE."""
import os
from pathlib import Path
STAGE=Path(os.environ["QEC360_STAGE"]).resolve()
import json
r = [json.loads(l) for l in open(str(STAGE / 'qec360-certprobe/gb_orbit_certificates.jsonl'))]
m = json.load(open(str(STAGE / 'qec360-certprobe/gb_orbit_metadata.json')))['masks']
threes = [i for i, z in enumerate(r) if m[z['fixedMask']]['fixedPoints'] == 3]
s = 'import OrbitBatchData\nimport OrbitBatchCheck\nimport Mathlib.Tactic.FinCases\nset_option maxRecDepth 100000\nset_option maxHeartbeats 0\nnamespace OrbitBatchData\nopen PackedDigits\n\n/-- The fifteen nonzero representatives with a three-element translation stabilizer. -/\ndef fixedThree : Finset Nat := {' + ','.join(map(str, threes)) + '}\n\n/-- The certified stabilizer count at a global census index. -/\ndef fixedAt (i : Nat) : Nat := if i = 0 then 45 else if i = 70587 then 9\n  else if i ∈ fixedThree then 3 else 1\n\n/-- Every stored batch word fits in exactly its allocated lanes. -/\ntheorem words_bounded : ∀ b : Fin 268, wordsAt b.val < 2^(128*sizeAt b.val) := by\n  decide +kernel\n\n/-- The expected packed stabilizer counts agree with their sparse global description. -/\ntheorem expected_eq : ∀ b : Fin 268,\n    expectedAt b.val = packF 128 (sizeAt b.val) (fun j => fixedAt (1024*b.val+j)) := by\n  decide +kernel\n\n/-- Every stabilizer count fits within one packed lane. -/\ntheorem fixedAt_lt (i : Nat) : fixedAt i < 2^128 := by\n  unfold fixedAt\n  split_ifs <;> decide\n\n/-- Consecutive batches are strictly separated in natural mask order. -/\ntheorem cross_sorted : ∀ b : Fin 267,\n    lane 128 1023 (wordsAt b.val) < lane 128 0 (wordsAt (b.val+1)) := by\n  decide +kernel\n\n/-- Every stored batch satisfies the same numerical checker. -/\ntheorem all_checked (b : Fin 268) :\n    OrbitBatchCheck.batchCheck (repeat128 (sizeAt b.val)) (wordsAt b.val) (expectedAt b.val) = true := by\n  fin_cases b\n'
for i in range(268):
    s += f'  · exact OrbitBatchCheck.checked{i}\n'
s += '\nend OrbitBatchData\n'
open(str(STAGE / 'OrbitBatchMetadata.lean'), 'w').write(s)
