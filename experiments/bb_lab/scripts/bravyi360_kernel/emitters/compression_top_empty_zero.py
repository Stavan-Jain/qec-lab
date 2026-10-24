import os
from pathlib import Path as _StagePath
_STAGE = _StagePath(os.environ["QEC360_STAGE"])
import pathlib,json,sys
sys.path.insert(0,str(_STAGE / 'qec360-certprobe'))
from balance_literals import lit
P=pathlib.Path(str(_STAGE / 'qec360-certprobe'));rs=[json.loads(l)for l in open(P/'top_empty_dual90_certificates.jsonl')if json.loads(l)['cap']==0];f=open(P/'Bravyi360TopEmpty0.lean','w')
f.write('''import Bravyi360TopEmptyData
import Bravyi360DualData
import BatchDualRaw
import OutsideMaskBounds
set_option maxHeartbeats 0
set_option maxRecDepth 100000
namespace Quantum.Stabilizer.Homological.BB.Bravyi360.TopEmpty
open F2Mask FiberData PackedDigits
/-- Repetition mask used by the 512-bit parallel certificates. -/
def rep512 (size : Nat) : Nat := ((1 <<< (512*size))-1)/((1 <<< 512)-1)
''')
for size in[1024,len(rs)%1024]:
 f.write(f'/-- Kernel-checked repetition mask for batches of size{size}. -/\ntheorem rep{size}_eq : rep512 {size} = packF 512 {size} (fun _ => 1) := by decide +kernel\n')
for k,start in enumerate(range(0,len(rs),1024)):
 vals=rs[start:start+1024];size=len(vals);coeff=sum(int(r['duals'][0],16)<<(512*i)for i,r in enumerate(vals));b=f'beta0Words{k}';c=f'coeff0Words{k}';rep=f'(rep512 {size})'
 f.write(f'/-- Dual coefficients for empty top fibers in batch{k}. -/\ndef {c} : Nat := {lit(coeff)}\n')
 f.write(f'/-- Parallel checks of the empty top fibers in batch{k}. -/\ntheorem checked0_{k} :\n    {b} < 2^(512*{size}) ∧ {c} < 2^(512*{size}) ∧\n    {b} &&& ({rep}*(2^180-1)) = {b} ∧\n    batchRecon {rep} {c} topDualColumns &&& {b} = 0 ∧\n    rawParityFold 8 ((batchRecon {rep} {c} topDualColumns >>> 180) &&& {b}) &&& {rep} = {rep} := by decide +kernel\n')
 f.write(f'/-- Every lane of top empty batch{k} has an impossible supported fiber. -/\ntheorem impossible0_{k} (i : Nat) (hi : i < {size}) (word : Nat)\n    (hs : word &&& lane 512 i {b} = word) :\n    reconX topEColumns word ≠ reconX topRColumns (lane 512 i {b}) := by\n  have h := checked0_{k}\n  rw [rep{size}_eq] at h\n  have hd := batchDual0_raw_sound 512 {size} 180 8 {b} {c} (by decide) (by decide) (by decide)\n    h.1 h.2.1 h.2.2.1 topDualColumns topDualColumns_bounded h.2.2.2.1 h.2.2.2.2 i hi\n  exact topDual_no_supported _ _ _ hd.1 hd.2 hs\n')
f.write('''/-- Size of one cap-zero batch. -/
def size0 (b : Nat) : Nat := if b < 24 then 1024 else 301
/-- Every cap-zero batch has an impossible supported fiber. -/
theorem cap0_batch_impossible (b : Nat) (hb : b < 25) (i : Nat) (hi : i < size0 b) (word : Nat)
    (hs : word &&& lane 512 i (beta0Words b) = word) :
    reconX topEColumns word ≠ reconX topRColumns (lane 512 i (beta0Words b)) := by
''')
def branch(lo,hi,indent):
 if hi-lo==1:
  f.write(indent+f'have he : b = {lo} := by omega\n'+indent+'subst b\n'+indent+f'exact impossible0_{lo} i hi word hs\n')
 else:
  mid=(lo+hi)//2;f.write(indent+f'by_cases h{mid} : b < {mid}\n'+indent+'·\n');branch(lo,mid,indent+'  ');f.write(indent+'·\n');branch(mid,hi,indent+'  ')
branch(0,25,'  ')
f.write('''/-- Every listed cap-zero top fiber is impossible on its shadow support. -/
theorem cap0_no_supported (i : Fin 24877) (word : Nat) (hs : word &&& beta0 i.val = word) :
    reconX topEColumns word ≠ reconX topRColumns (beta0 i.val) := by
  have hm := Nat.mod_lt i.val (by decide : 0 < 1024)
  have he := Nat.mod_add_div i.val 1024
  have hb : i.val / 1024 < 25 := by omega
  have hl : i.val % 1024 < size0 (i.val / 1024) := by
    unfold size0
    split_ifs <;> omega
  exact cap0_batch_impossible _ hb _ hl word hs
/-- Every cap-zero top fiber has no bounded solution with zero outside weight. -/
theorem cap0_impossible (i : Fin 24877) (word : Nat) (hw : word < 2^180)
    (he : reconX topEColumns word = reconX topRColumns (beta0 i.val))
    (hoff : bitCount 180 (outsideMask word (beta0 i.val)) ≤ 0) : False :=
  cap0_no_supported i word (supported_of_outside_weight_zero 180 word _ hw hoff) he
#print axioms cap0_impossible
end Quantum.Stabilizer.Homological.BB.Bravyi360.TopEmpty
''');f.close();print('topzero',len(rs),'records')

# Data are kernel certificates, not runtime programs.
for _certificate_output in ['Bravyi360TopEmpty0.lean']:
 _certificate_path = P / _certificate_output
 _certificate_path.write_text(_certificate_path.read_text().replace("\ndef ", "\nnoncomputable def "))
