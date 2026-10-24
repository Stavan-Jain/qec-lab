import os
from pathlib import Path as _StagePath
_STAGE = _StagePath(os.environ["QEC360_STAGE"])
import pathlib
P=pathlib.Path(str(_STAGE / 'qec360-certprobe'));f=open(P/'Bravyi360BottomEmptyWeights.lean','w')
f.write('''import Bravyi360BottomEmptyData
import BatchWidePopcount
set_option maxHeartbeats 0
set_option maxRecDepth 100000
namespace Quantum.Stabilizer.Homological.BB.Bravyi360.BottomEmpty
open F2Mask PackedDigits
''')
for cap,N in enumerate([224352,20904,859]):
 weight=22-2*cap;chunks=(N+1023)//1024;last=N%1024
 for ix,off in enumerate(range(0,N,1024)):
  size=min(1024,N-off);b=f'beta{cap}Words{ix}';rep=f'(rep256 {size})'
  f.write(f'''/-- Packed exact shadow weights in cap{cap} batch {ix}. -/
theorem weights{cap}_{ix} : {b} < 2^(256*{size}) ∧
    {b} &&& ({rep}*(2^90-1)) = {b} ∧
    popcount256 {rep} {b} = {rep}*{weight} := by decide +kernel
/-- Every shadow in cap{cap} batch {ix} has the exact claimed weight. -/
theorem weight{cap}_{ix} (i : Nat) (hi : i < {size}) : bitCount 90 (lane 256 i {b}) = {weight} := by
  have h := weights{cap}_{ix}
  rw [rep{size}_eq] at h
  have hp := h.2.2
  rw [repeat_mul_eq_packF] at hp
  exact popcount256_weight {size} 90 {b} {weight} (by decide) h.1 h.2.1 hp (by decide) i hi
''')
 f.write(f'''/-- Exact weights for all cap{cap} chunks. -/
theorem weight{cap}_batch (b : Nat) (hb : b < {chunks}) (i : Nat)
    (hi : i < if b < {chunks-1} then 1024 else {last}) :
    bitCount 90 (lane 256 i (beta{cap}Words b)) = {weight} := by
''')
 def branch(lo,hi,ind):
  if hi-lo==1:f.write(ind+f'have heq : b = {lo} := by omega\n'+ind+'subst b\n'+ind+f'exact weight{cap}_{lo} i hi\n')
  else:
   mid=(lo+hi)//2;f.write(ind+f'by_cases h{mid} : b < {mid}\n'+ind+'·\n');branch(lo,mid,ind+'  ');f.write(ind+'·\n');branch(mid,hi,ind+'  ')
 branch(0,chunks,'  ')
 f.write(f'''/-- Exact shadow weight throughout the cap{cap} family. -/
theorem beta{cap}_weight (i : Fin {N}) : bitCount 90 (beta{cap} i.val) = {weight} := by
  have hm := Nat.mod_lt i.val (by decide : 0 < 1024)
  have heq := Nat.mod_add_div i.val 1024
  have hb : i.val / 1024 < {chunks} := by omega
  have hl : i.val % 1024 < (if i.val / 1024 < {chunks-1} then 1024 else {last}) := by
    split_ifs <;> omega
  exact weight{cap}_batch _ hb _ hl
''')
f.write('#print axioms beta0_weight\nend Quantum.Stabilizer.Homological.BB.Bravyi360.BottomEmpty\n');f.close()

# Data are kernel certificates, not runtime programs.
for _certificate_output in ['Bravyi360BottomEmptyWeights.lean']:
 _certificate_path = P / _certificate_output
 _certificate_path.write_text(_certificate_path.read_text().replace("\ndef ", "\nnoncomputable def "))
