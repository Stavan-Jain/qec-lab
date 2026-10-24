import os
from pathlib import Path as _StagePath
_STAGE = _StagePath(os.environ["QEC360_STAGE"])
import pathlib,json,sys,io
sys.path.insert(0,str(_STAGE / 'qec360-certprobe'))
from balance_literals import lit
P=pathlib.Path(str(_STAGE / 'qec360-certprobe'))
modules=[]
for prefix,n,frame in [('bottom',90,json.load(open(P/'gb_frame.json'))),('top',180,json.load(open(P/'top_frame.json')))]:
 for letter in['E','R']:
  name=prefix+letter;cols=frame[letter+'cols'];m=max(n,max(c.bit_length()for c in cols));assert m==n
  module='Bravyi360'+name[0].upper()+name[1:]+'Adjoint';modules.append(module)
  f=io.StringIO();f.write('import Bravyi360FiberData\nimport MaskAdjoint\nimport MaskDual\nimport ReconWindow\nset_option maxHeartbeats 0\nset_option maxRecDepth 100000\nnamespace Quantum.Stabilizer.Homological.BB.Bravyi360\nopen F2Mask FiberData\n')
  rows=[sum(((c>>j)&1)<<i for i,c in enumerate(cols))for j in range(m)]
  rt=sum(c<<(n*i)for i,c in enumerate(rows));ct=sum(c<<(m*i)for i,c in enumerate(cols))
  f.write(f'/-- Full transposed {name} matrix columns. -/\ndef {name}TransposeColumns : List (Nat × Nat) := ['+','.join('('+str(i)+','+hex(c)+')'for i,c in enumerate(rows))+']\n')
  f.write(f'/-- Packed rows for fixed-matrix certification. -/\ndef {name}RowTable : Nat := {lit(rt)}\n')
  f.write(f'/-- Packed columns for fixed-matrix certification. -/\ndef {name}ColumnTable : Nat := {lit(ct)}\n')
  f.write(f'/-- Read one transposed matrix column. -/\ndef {name}Row (j : Nat) : Nat := Nat.land (Nat.shiftRight {name}RowTable (Nat.mul {n} j)) (2^{n}-1)\n')
  f.write(f'/-- Read one forward matrix column. -/\ndef {name}Column (i : Nat) : Nat := Nat.land (Nat.shiftRight {name}ColumnTable (Nat.mul {m} i)) (2^{m}-1)\n')
  f.write(f'/-- Transposed columns reconstruct their certified point masses. -/\ntheorem {name}Transpose_point : ∀ j : Fin {m}, reconX {name}TransposeColumns (2^j.val) = {name}Row j := by decide +kernel\n')
  f.write(f'/-- Forward columns reconstruct their certified point masses. -/\ntheorem {name}Forward_point : ∀ i : Fin {n}, reconX {name}Columns (2^i.val) = {name}Column i := by decide +kernel\n')
  f.write(f'/-- The two fixed lookup tables are transposes. -/\ntheorem {name}Transpose_bits : ∀ (i : Fin {n}) (j : Fin {m}), ({name}Row j).testBit i = ({name}Column i).testBit j := by decide +kernel\n')
  f.write(f'/-- Transposition gives the adjoint identity on point masses. -/\ntheorem {name}Transpose_single : ∀ (i : Fin {n}) (j : Fin {m}),\n    reconLinear (Equiv.refl (Fin {m})) (Equiv.refl (Fin {n})) {name}TransposeColumns (Pi.single j 1) i =\n    reconLinear (Equiv.refl (Fin {n})) (Equiv.refl (Fin {m})) {name}Columns (Pi.single i 1) j := by\n  intro i j\n  rw [reconLinear_single, reconLinear_single, {name}Transpose_point, {name}Forward_point]\n  simp only [chainOfMask, Equiv.refl_apply, {name}Transpose_bits]\n  rfl\n')
  f.write(f'/-- Both pivot sets and the transposed outputs fit their coordinate ranges. -/\ntheorem {name}Adjoint_bounds : (∀ p ∈ {name}Columns, p.1 < {n}) ∧ (∀ p ∈ {name}TransposeColumns, p.1 < {m} ∧ p.2 < 2^{n}) := by decide +kernel\n')
  f.write(f'/-- The packed transpose computes the full-matrix dual parity. -/\ntheorem {name}Transpose_parity (lam v : Nat) :\n    dualOfMask {m} lam (chainOfMask (Equiv.refl (Fin {m})) (reconX {name}Columns v)) =\n      (bitCount {n} (reconX {name}TransposeColumns lam &&& v) : ZMod 2) := by\n  have h := dualOfMask_apply_of_adjoint\n    (reconLinear (Equiv.refl (Fin {n})) (Equiv.refl (Fin {m})) {name}Columns)\n    (reconLinear (Equiv.refl (Fin {m})) (Equiv.refl (Fin {n})) {name}TransposeColumns)\n    {name}Transpose_single lam v (reconX {name}TransposeColumns lam)\n    (reconLinear_chainOfMask _ _ _ (fun p hp => ({name}Adjoint_bounds.2 p hp).1) lam)\n  rwa [reconLinear_chainOfMask _ _ _ {name}Adjoint_bounds.1 v] at h\n')
  f.write(f'/-- The raw transpose identity applies without input-mask bounds. -/\ntheorem {name}Transpose_parityFold (lam v : Nat) :\n    parityFold 8 (reconX {name}TransposeColumns lam &&& v) =\n      parityFold 8 (lam &&& reconX {name}Columns v) := by\n  have h := {name}Transpose_parity lam v\n  rw [dualOfMask_chainOfMask] at h\n  exact parityFold_eq_of_cast_eq {n} 8 _ _ (by decide)\n    (lt_of_le_of_lt Nat.and_le_left (reconX_lt _ _ _ (fun p hp => ({name}Adjoint_bounds.2 p hp).2)))\n    (lt_of_le_of_lt Nat.and_le_right (reconX_lt _ _ _ {name}Columns_bounded)) h.symm\n')
  f.write(f'#print axioms {name}Transpose_parityFold\n')
  f.write('end Quantum.Stabilizer.Homological.BB.Bravyi360\n')
  (P/(module+'.lean')).write_text(f.getvalue())
(P/'Bravyi360Adjoints.lean').write_text(''.join('import '+m+'\n'for m in modules))

# Data are kernel certificates, not runtime programs.
for _certificate_output in ['Bravyi360BottomEAdjoint.lean', 'Bravyi360BottomRAdjoint.lean', 'Bravyi360TopEAdjoint.lean', 'Bravyi360TopRAdjoint.lean', 'Bravyi360Adjoints.lean']:
 _certificate_path = P / _certificate_output
 _certificate_path.write_text(_certificate_path.read_text().replace("\ndef ", "\nnoncomputable def "))
