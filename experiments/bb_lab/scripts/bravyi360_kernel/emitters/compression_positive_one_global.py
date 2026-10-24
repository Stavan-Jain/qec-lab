import os
from pathlib import Path as _StagePath
_STAGE = _StagePath(os.environ["QEC360_STAGE"])
import pathlib
P=pathlib.Path(str(_STAGE / 'qec360-certprobe'));C=112
f=open(P/'Bravyi360BottomPositiveOne.lean','w');f.write('import Bravyi360BottomForest\n'+''.join(f'import BottomPositiveOnePart{i}\n'for i in range(14))+'set_option maxRecDepth 100000\nset_option maxHeartbeats 0\nnamespace Quantum.Stabilizer.Homological.BB.Bravyi360.BottomForest\nopen F2Mask FiberData\n')
def branches(lo,hi,indent,kind):
 if hi-lo==1:
  k=lo*250;f.write(indent+f'have heq : b = {lo} := by omega\n'+indent+'subst b\n')
  if kind=='one':f.write(indent+f'exact oneGroup_complete records{k} oneCovers{k} onePacked{k} checked{k} oneChecked{k} r hr hweight word hw he hoff\n')
  if kind=='zero':f.write(indent+f'exact zeroGroup_complete records{k} zeroChecked{k} r hr hweight word hw hoff\n')
  if kind=='weight':f.write(indent+f'have h := List.all_eq_true.mp weightChecked{k} r hr\n'+indent+'simpa using h\n')
 else:
  mid=(lo+hi)//2;f.write(indent+f'by_cases h{mid} : b < {mid}\n'+indent+'·\n');branches(lo,mid,indent+'  ',kind);f.write(indent+'·\n');branches(mid,hi,indent+'  ',kind)
f.write('''/-- Every cap-one forest chunk has complete outside-support seeds. -/
theorem one_records_complete (b : Nat) (hb : b < 112) (r : Fiber)
    (hr : r ∈ recordsAt b) (hweight : r.weight = 20)
    (word : Nat) (hw : word < 2^90)
    (he : reconX bottomEColumns word = reconX bottomRColumns r.beta)
    (hoff : bitCount 90 (outsideMask word r.beta) ≤ 1) :
    outsideMask word r.beta ∈ r.seeds.map Seed.off := by
''');branches(0,C,'  ','one')
f.write('''/-- Every cap-zero forest chunk lists the unique empty outside support. -/
theorem zero_records_complete (b : Nat) (hb : b < 112) (r : Fiber)
    (hr : r ∈ recordsAt b) (hweight : r.weight = 22)
    (word : Nat) (hw : word < 2^90)
    (hoff : bitCount 90 (outsideMask word r.beta) ≤ 0) :
    outsideMask word r.beta ∈ r.seeds.map Seed.off := by
''');branches(0,C,'  ','zero')
f.write('''/-- The five shallow positive shadow weights exhaust every forest chunk. -/
theorem records_weights (b : Nat) (hb : b < 112) (r : Fiber)
    (hr : r ∈ recordsAt b) : r.weight ∈ [14,16,18,20,22] := by
''');branches(0,C,'  ','weight')
f.write('''/-- The canonical positive forest family has only the five shallow weights. -/
theorem fiberAt_weights (i : Fin 27867) : (fiberAt i.val).weight ∈ [14,16,18,20,22] :=
  records_weights _ (by omega) _ (fiberAt_mem i)
/-- Complete outside-support coverage for every positive cap-one fiber. -/
theorem oneFiberAt_complete (i : Fin 27867) (hweight : (fiberAt i.val).weight = 20)
    (word : Nat) (hw : word < 2^90)
    (he : reconX bottomEColumns word = reconX bottomRColumns (fiberAt i.val).beta)
    (hoff : bitCount 90 (outsideMask word (fiberAt i.val).beta) ≤ 1) :
    outsideMask word (fiberAt i.val).beta ∈ (fiberAt i.val).seeds.map Seed.off :=
  one_records_complete _ (by omega) _ (fiberAt_mem i) hweight word hw he hoff
/-- Complete outside-support coverage for every positive cap-zero fiber. -/
theorem zeroFiberAt_complete (i : Fin 27867) (hweight : (fiberAt i.val).weight = 22)
    (word : Nat) (hw : word < 2^90)
    (hoff : bitCount 90 (outsideMask word (fiberAt i.val).beta) ≤ 0) :
    outsideMask word (fiberAt i.val).beta ∈ (fiberAt i.val).seeds.map Seed.off :=
  zero_records_complete _ (by omega) _ (fiberAt_mem i) hweight word hw hoff
#print axioms oneFiberAt_complete
end Quantum.Stabilizer.Homological.BB.Bravyi360.BottomForest
''');f.close()

# Data are kernel certificates, not runtime programs.
for _certificate_output in ['Bravyi360BottomPositiveOne.lean']:
 _certificate_path = P / _certificate_output
 _certificate_path.write_text(_certificate_path.read_text().replace("\ndef ", "\nnoncomputable def "))
