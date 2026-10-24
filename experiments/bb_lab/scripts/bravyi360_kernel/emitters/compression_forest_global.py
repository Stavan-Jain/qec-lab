import os
from pathlib import Path as _StagePath
_STAGE = _StagePath(os.environ["QEC360_STAGE"])
import pathlib
P=pathlib.Path(str(_STAGE / 'qec360-certprobe'));N=27867;B=250;C=(N+B-1)//B
f=open(P/'Bravyi360BottomForest.lean','w')
f.write('import BottomForestGroup\n'+''.join(f'import BottomForestFastPart{i}\n'for i in range(14))+'set_option maxRecDepth 100000\nset_option maxHeartbeats 0\nnamespace Quantum.Stabilizer.Homological.BB.Bravyi360.BottomForest\nopen F2Mask FiberData\n')
def dispatch(names,lo=0):
 if len(names)==1:return names[0]
 mid=len(names)//2;return f'(if Nat.blt b {lo+mid} then {dispatch(names[:mid],lo)} else {dispatch(names[mid:],lo+mid)})'
f.write('/-- Balanced lookup of canonical shallow-fiber record chunks. -/\nnoncomputable def recordsAt (b : Nat) : List Fiber := '+dispatch([f'records{i*B}'for i in range(C)])+'\n')
f.write('/-- Canonical shallow-fiber record at its fixed family position. -/\nnoncomputable def fiberAt (i : Nat) : Fiber := (recordsAt (i/250)).getD (i%250) default\n')
for i in range(C):
 size=min(B,N-i*B);f.write(f'/-- Number of records in forest chunk {i}. -/\ntheorem length{i} : records{i*B}.length = {size} := by decide +kernel\n')
f.write(f'''/-- Every canonical record chunk passes its full checker. -/
theorem recordsAt_checked (b : Nat) (hb : b < {C}) : groupCheck (recordsAt b) = true := by
''')
def branches(lo,hi,ind,kind):
 if hi-lo==1:f.write(ind+f'have heq : b = {lo} := by omega\n'+ind+'subst b\n'+ind+('exact checked'+str(lo*B) if kind=='checked'else'exact length'+str(lo))+'\n')
 else:
  mid=(lo+hi)//2;f.write(ind+f'by_cases h{mid} : b < {mid}\n'+ind+'·\n');branches(lo,mid,ind+'  ',kind);f.write(ind+'·\n');branches(mid,hi,ind+'  ',kind)
branches(0,C,'  ','checked')
f.write(f'''/-- The last forest record chunk has its precise smaller size. -/
theorem recordsAt_length (b : Nat) (hb : b < {C}) :
    (recordsAt b).length = if b < {C-1} then 250 else {N%B} := by
''');branches(0,C,'  ','length')
f.write(f'''/-- A valid global fiber index selects a real record of its chunk. -/
theorem fiberAt_mem (i : Fin {N}) : fiberAt i.val ∈ recordsAt (i.val/250) := by
  have hm := Nat.mod_lt i.val (by decide : 0 < 250)
  have he := Nat.mod_add_div i.val 250
  have hb : i.val/250 < {C} := by omega
  have hi : i.val%250 < (recordsAt (i.val/250)).length := by
    rw [recordsAt_length _ hb]
    split_ifs <;> omega
  unfold fiberAt
  rw [List.getD_eq_getElem _ _ hi]
  exact List.getElem_mem hi
/-- Every valid global forest record passes the structural checker. -/
theorem fiberAt_checked (i : Fin {N}) : fiberCheck (fiberAt i.val) = true :=
  groupCheck_member _ (recordsAt_checked _ (by omega)) _ (fiberAt_mem i)
/-- Every shallow forest solution above a listed outside support reaches the common top bank. -/
theorem fiberAt_complete (i : Fin {N}) (word : Nat)
    (he : reconX bottomEColumns word = reconX bottomRColumns (fiberAt i.val).beta)
    (hoff : outsideMask word (fiberAt i.val).beta ∈ (fiberAt i.val).seeds.map Seed.off) :
    ∃ j : Fin 58130, bottomJoin word (fiberAt i.val).beta = topBankMask j ∨
      bottomJoin (word ^^^ (fiberAt i.val).beta) (fiberAt i.val).beta = topBankMask j :=
  group_solution_bank _ (recordsAt_checked _ (by omega)) _ (fiberAt_mem i) word he hoff
#print axioms fiberAt_complete
end Quantum.Stabilizer.Homological.BB.Bravyi360.BottomForest
''');f.close()

# Data are kernel certificates, not runtime programs.
for _certificate_output in ['Bravyi360BottomForest.lean']:
 _certificate_path = P / _certificate_output
 _certificate_path.write_text(_certificate_path.read_text().replace("\ndef ", "\nnoncomputable def "))
