import os
from pathlib import Path as _StagePath
_STAGE = _StagePath(os.environ["QEC360_STAGE"])
import pathlib,json,sys
sys.path.insert(0,str(_STAGE / 'qec360-certprobe'))
P=pathlib.Path(str(_STAGE / 'qec360-certprobe'));f=open(P/'Bravyi360DualData.lean','w')
f.write('''import Bravyi360Adjoints
import ReconComposition
set_option maxHeartbeats 0
set_option maxRecDepth 100000
namespace Quantum.Stabilizer.Homological.BB.Bravyi360
open F2Mask FiberData
''')
for name,n,k,width in [('bottom',90,45,256),('top',180,90,512)]:
 frame=json.load(open(P/('gb_frame.json'if name=='bottom'else'top_frame90.json')))
 selected=list(range(k))if name=='bottom'else frame['selectedRows']
 rows=[sum(((c>>i)&1)<<j for j,c in enumerate(frame['Ecols']))|(sum(((c>>i)&1)<<j for j,c in enumerate(frame['Rcols']))<<n)for i in range(k)]
 f.write(f'/-- Combined E and R adjoints for the compressed {name} syndrome. -/\ndef {name}DualColumns : List (Nat × Nat) := ['+','.join(f'({i},{hex(v)})'for i,v in enumerate(rows))+']\n')
 f.write(f'/-- Embed a compressed {name} dual into the original syndrome coordinates. -/\ndef {name}DualLiftColumns : List (Nat × Nat) := ['+','.join(f'({i},{hex(1<<q)})'for i,q in enumerate(selected))+']\n')
 f.write(f'/-- Combined dual entries fit their parallel arithmetic lanes. -/\ntheorem {name}DualColumns_bounded : ∀ p ∈ {name}DualColumns, p.1 < {width} ∧ p.2 < 2^{width} := by decide +kernel\n')
 for letter in['E','R']:
  op=f'p.2 &&& (2^{n}-1)'if letter=='E'else f'p.2 >>> {n}'
  f.write(f'/-- The compressed {letter} dual columns agree with the full matrix. -/\ntheorem {name}Dual{letter}_columns :\n    {name}DualLiftColumns.map (fun p => (p.1,reconX {name}{letter}TransposeColumns p.2)) =\n      {name}DualColumns.map (fun p => (p.1,{op})) := by decide +kernel\n')
  rhs=f'reconX {name}DualColumns coeff &&& (2^{n}-1)'if letter=='E'else f'reconX {name}DualColumns coeff >>> {n}'
  rw='reconX_and'if letter=='E'else'reconX_shiftRight'
  f.write(f'/-- Every compressed {letter} dual evaluates as its original full-matrix dual. -/\ntheorem {name}Dual_{letter} (coeff : Nat) :\n    reconX {name}{letter}TransposeColumns (reconX {name}DualLiftColumns coeff) = {rhs} := by\n  rw [reconX_comp,{name}Dual{letter}_columns,← {rw}]\n')
 if name=='bottom':
  for letter in['E','R']:
   f.write(f'/-- The 90-bit {letter} adjoint also uses the 128-bit parity fold. -/\ntheorem bottom{letter}Transpose_parityFold7 (lam word : Nat) :\n    parityFold 7 (reconX bottom{letter}TransposeColumns lam &&& word) =\n      parityFold 7 (lam &&& reconX bottom{letter}Columns word) := by\n  rw [parityFold_change_width 90 7 8 _ (by decide) (by decide)\n    (lt_of_le_of_lt Nat.and_le_left (reconX_lt _ _ _ (fun p hp => (bottom{letter}Adjoint_bounds.2 p hp).2))),\n    parityFold_change_width 90 7 8 _ (by decide) (by decide)\n    (lt_of_le_of_lt Nat.and_le_right (reconX_lt _ _ _ bottom{letter}Columns_bounded))]\n  exact bottom{letter}Transpose_parityFold lam word\n')
 kfold=7 if name=='bottom'else 8;suffix='7'if name=='bottom'else''
 f.write(f'/-- A compressed dual obstruction excludes every on-shadow solution\nfor the actual {name} matrices. -/\ntheorem {name}Dual_no_supported (coeff beta word : Nat)\n    (hz : reconX {name}DualColumns coeff &&& beta = 0)\n    (ho : parityFold {kfold} ((reconX {name}DualColumns coeff >>> {n}) &&& beta) = 1)\n    (hs : word &&& beta = word) : reconX {name}EColumns word ≠ reconX {name}RColumns beta := by\n  apply no_supported_solution_of_raw_dual {kfold} {name}EColumns {name}RColumns\n    {name}ETransposeColumns {name}RTransposeColumns {name}ETranspose_parityFold{suffix} {name}RTranspose_parityFold{suffix}\n    (reconX {name}DualLiftColumns coeff) beta word\n  · rw [{name}Dual_E]\n    calc\n      (reconX {name}DualColumns coeff &&& (2^{n}-1)) &&& beta =\n          (reconX {name}DualColumns coeff &&& beta) &&& (2^{n}-1) := by ac_rfl\n      _ = 0 := by simp [hz]\n  · rw [{name}Dual_R]\n    exact ho\n  · exact hs\n')
f.write('#print axioms bottomDual_no_supported\n#print axioms topDual_no_supported\nend Quantum.Stabilizer.Homological.BB.Bravyi360\n');f.close()

# Data are kernel certificates, not runtime programs.
for _certificate_output in ['Bravyi360DualData.lean']:
 _certificate_path = P / _certificate_output
 _certificate_path.write_text(_certificate_path.read_text().replace("\ndef ", "\nnoncomputable def "))
