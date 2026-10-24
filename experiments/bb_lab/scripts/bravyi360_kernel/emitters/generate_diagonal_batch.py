"""Emit checked certificate data into the staging directory selected by QEC360_STAGE."""
import os
from pathlib import Path
STAGE = Path(os.environ["QEC360_STAGE"]).resolve()
from pathlib import Path
import json, sys
sys.path.insert(0, str(STAGE / 'qec360-certprobe'))
from balance_literals import lit
D = json.load(open(str(STAGE / 'diagonal_cap_one_duals.json')))
F = json.load(open(str(STAGE / 'qec360-certprobe/top_frame.json')))
E = F['Ecols']
R = F['Rcols']
cols = [sum(((c >> q & 1) << i for i, c in enumerate(E))) + (sum(((c >> q & 1) << i for i, c in enumerate(R))) << 180) for q in range(180)]
betas = [int(d['beta'], 0) for d in D]
ls = [[int(x, 0) for x in d['lambdas']] for d in D]
assert len(D) == 36
packed = sum((x << 512 * i for i, x in enumerate(betas)))
cs = [sum(((ls[i][j] if j < len(ls[i]) else 0) << 512 * i for i in range(36))) for j in range(max(map(len, ls)))]
rep = sum((1 << 512 * i for i in range(36)))
s = 'import Bravyi360DiagonalData\nimport BatchCap1Raw\nnamespace Quantum.Stabilizer.Homological.BB.Bravyi360.DiagonalBatch\nopen F2Mask FiberData PackedDigits\nset_option maxRecDepth 200000\nset_option maxHeartbeats 0\nset_option exponentiation.threshold 1000000\n'
s += '/-- Full dual rows for the diagonal fibers, combining both carry maps. -/\ndef rows : List (Nat × Nat) := [' + ','.join((f'({i},{hex(c)})' for i, c in enumerate(cols))) + ']\n'
s += '/-- Identity coordinates for the full dual rows. -/\ndef lifts : List (Nat × Nat) := [' + ','.join((f'({i},{hex(1 << i)})' for i in range(180))) + ']\n'
s += f'/-- The thirty-six diagonal shadows in parallel lanes. -/\ndef betas : Nat := {lit(packed)}\n'
s += '/-- Dual coefficient tables excluding every singleton outside support. -/\ndef coeffs : List Nat := [' + ','.join((lit(x) for x in cs)) + ']\n'
s += f'/-- Lane selector for the diagonal batch. -/\ndef rep : Nat := {lit(rep)}\n'
s += '/-- The lane selector is the generic repetition mask. -/\ntheorem rep_eq : rep = packF 512 36 (fun _ => 1) := by decide +kernel\n/-- The full combined dual columns fit the packed lane width. -/\ntheorem rows_bounded : ∀ p ∈ rows, p.1 < 512 ∧ p.2 < 2^512 := by decide +kernel\n/-- The low half of the dual rows is the actual transposed E map. -/\ntheorem e_columns :\n    lifts.map (fun p => (p.1,reconX topETransposeColumns p.2)) =\n      rows.map (fun p => (p.1,p.2 &&& (2^180-1))) := by decide +kernel\n/-- The high half of the dual rows is the actual transposed R map. -/\ntheorem r_columns :\n    lifts.map (fun p => (p.1,reconX topRTransposeColumns p.2)) =\n      rows.map (fun p => (p.1,p.2 >>> 180)) := by decide +kernel\n/-- The combined dual map computes the physical E dual for every coefficient mask. -/\ntheorem dual_E (c : Nat) :\n    reconX topETransposeColumns (reconX lifts c) = reconX rows c &&& (2^180-1) := by\n  rw [reconX_comp,e_columns,← reconX_and]\n/-- The combined dual map computes the physical R dual for every coefficient mask. -/\ntheorem dual_R (c : Nat) :\n    reconX topRTransposeColumns (reconX lifts c) = reconX rows c >>> 180 := by\n  rw [reconX_comp,r_columns,← reconX_shiftRight]\n/-- Every full diagonal dual constrains the actual outside support. -/\ntheorem constraint (coeff beta word : Nat) (hw : word < 2^180)\n    (hann : reconX rows coeff &&& beta = 0)\n    (he : reconX topEColumns word = reconX topRColumns beta) :\n    parityFold 8 (reconX rows coeff &&& outsideMask word beta) =\n      parityFold 8 ((reconX rows coeff >>> 180) &&& beta) :=\n  compressed_dual_constraint 180 8 _ _ _ _ _ _ topETranspose_parityFold\n    topRTranspose_parityFold dual_E dual_R coeff beta word hw hann he\n/-- Every lane is its stated physical diagonal shadow. -/\ntheorem betas_checked : ∀ i : Fin 36, lane 512 i betas = DiagonalData.emptyBeta i := by decide +kernel\n/-- Parallel dual checks exclude zero and singleton outside supports. -/\ntheorem checked :\n    betas < 2^(512*36) ∧ betas &&& (rep*(2^180-1)) = betas ∧\n    coeffs.all (fun c => Nat.blt c (2^(512*36)) && Nat.beq (batchRecon rep c rows &&& betas) 0) = true ∧\n    batchAllowed1 180 8 rep betas rows coeffs = 0 ∧\n    batchZero1 180 8 rep betas rows coeffs = 0 := by decide +kernel\n/-- The diagonal fibers admit no sheet with at most one outside coordinate. -/\ntheorem impossible (i : Fin 36) (word : Nat) (hw : word < 2^180)\n    (he : reconX topEColumns word = reconX topRColumns (DiagonalData.emptyBeta i))\n    (hoff : bitCount 180 (outsideMask word (DiagonalData.emptyBeta i)) ≤ 1) : False := by\n  have h := checked\n  rw [rep_eq] at h\n  have hc := h.2.2.1\n  simp only [List.all_eq_true,Bool.and_eq_true,Nat.blt_eq,Nat.beq_eq] at hc\n  apply cap1_raw_impossible 512 36 180 8 betas coeffs (by decide) (by decide) (by decide) (by decide)\n    h.1 h.2.1 rows topEColumns topRColumns rows_bounded\n    (fun c hm => (hc c hm).1) (fun c hm => (hc c hm).2) h.2.2.2.1 h.2.2.2.2\n    constraint i i.isLt word hw\n  · rwa [betas_checked]\n  · rwa [betas_checked]\nend Quantum.Stabilizer.Homological.BB.Bravyi360.DiagonalBatch\n'
Path(str(STAGE / 'Bravyi360DiagonalBatch.lean')).write_text(s)
print('Generated full diagonal dual batch', len(s), 'bytes')

# Certificate constants are reduced by the kernel and do not need runtime code.
import re
for generated in ['Bravyi360DiagonalBatch.lean']:
    path = STAGE / generated
    text = path.read_text()
    text = re.sub(r"(?m)^def\b", "noncomputable def", text)
    text = re.sub(r"(?m)^instance\b", "noncomputable instance", text)
    path.write_text(text)
