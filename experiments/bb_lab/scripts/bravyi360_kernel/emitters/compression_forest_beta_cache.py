import os
from pathlib import Path as _StagePath
_STAGE = _StagePath(os.environ["QEC360_STAGE"])
import pathlib
P=pathlib.Path(str(_STAGE / 'qec360-certprobe'));C=112
f=open(P/'Bravyi360BottomBetaCache.lean','w');f.write('''import Bravyi360BottomForest
import PackedListLookup
set_option maxRecDepth 100000
set_option maxHeartbeats 0
namespace Quantum.Stabilizer.Homological.BB.Bravyi360.BottomForest
open F2Mask PackedDigits
/-- The existing packed matrix input is also an exact cache of the shadow field. -/
theorem packedGroupCheck_beta (rs : List Fiber) (p : PackedGroup)
    (h : packedGroupCheck rs p = true) : packList 256 (rs.map Fiber.beta) = p.rInput := by
  simp only [packedGroupCheck,Bool.and_eq_true,Nat.beq_eq,and_assoc] at h
  rcases h with ⟨_,_,_,_,_,_,hrI,_,_,_,_,_,_,_⟩
  simpa only [rPairs,List.map_map,Function.comp_def] using hrI
''')
def dispatch(lo,hi):
 if hi-lo==1:return f'packed{lo*250}'
 mid=(lo+hi)//2;return f'(if Nat.blt b {mid} then {dispatch(lo,mid)} else {dispatch(mid,hi)})'
f.write('/-- Existing checked packed tables, addressed by their short forest chunk. -/\nnoncomputable def packedAt (b : Nat) : PackedGroup := '+dispatch(0,C)+'\n')
f.write('/-- Every packed forest chunk retains its existing numerical certificate. -/\ntheorem packedAt_checked (b : Nat) (hb : b < 112) : packedGroupCheck (recordsAt b) (packedAt b) = true := by\n')
def branch(lo,hi,ind):
 if hi-lo==1:f.write(ind+f'have heq : b = {lo} := by omega\n'+ind+'subst b\n'+ind+f'exact packed_checked{lo*250}\n')
 else:
  m=(lo+hi)//2;f.write(ind+f'by_cases h{m} : b < {m}\n'+ind+'·\n');branch(lo,m,ind+'  ');f.write(ind+'·\n');branch(m,hi,ind+'  ')
branch(0,C,'  ')
f.write('''/-- Read a shadow directly from its already checked packed cache. -/
noncomputable def betaAtFast (i : Nat) : Nat := lane 256 (i%250) (packedAt (i/250)).rInput
/-- The packed shadow accessor is exactly the field of the canonical forest record. -/
theorem betaAtFast_eq (i : Fin 27867) : betaAtFast i.val = (fiberAt i.val).beta := by
  have hm := Nat.mod_lt i.val (by decide : 0 < 250)
  have he := Nat.mod_add_div i.val 250
  have hb : i.val/250 < 112 := by omega
  have hi : i.val%250 < (recordsAt (i.val/250)).length := by
    rw [recordsAt_length _ hb]
    split_ifs <;> omega
  symm
  apply getD_field_eq_lane 256 (recordsAt (i.val/250)) Fiber.beta default
    (packedAt (i.val/250)).rInput
  · intro r hr
    exact lt_trans (fiberCheck_beta r (groupCheck_member _ (recordsAt_checked _ hb) r hr))
      (by decide : 2^90 < 2^256)
  · exact packedGroupCheck_beta _ _ (packedAt_checked _ hb)
  · exact hi
#print axioms betaAtFast_eq
end Quantum.Stabilizer.Homological.BB.Bravyi360.BottomForest
''');f.close()
