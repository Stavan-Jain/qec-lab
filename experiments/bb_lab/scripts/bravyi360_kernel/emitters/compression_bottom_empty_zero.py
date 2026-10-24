import os
from pathlib import Path as _StagePath
_STAGE = _StagePath(os.environ["QEC360_STAGE"])
import pathlib,json,sys
sys.path.insert(0,str(_STAGE / 'qec360-certprobe'))
from balance_literals import lit
from compact_literals import packlit
P=pathlib.Path(str(_STAGE / 'qec360-certprobe'));rs=[dict(beta=a[0],duals=a[2:])for l in open(P/'emptycap45_short5.tsv')if (a:=l.split())];f=open(P/'Bravyi360BottomEmpty0.lean','w')
f.write('''import Bravyi360BottomEmptyData
import Bravyi360DualData
import BatchDualRaw
import OutsideMaskBounds
set_option maxHeartbeats 0
set_option maxRecDepth 100000
namespace Quantum.Stabilizer.Homological.BB.Bravyi360.BottomEmpty
open F2Mask FiberData PackedDigits
''')
for k,start in enumerate(range(0,len(rs),1024)):
 vals=rs[start:start+1024];size=len(vals);coeff=sum(int(r['duals'][0],16)<<(256*i)for i,r in enumerate(vals));b=f'beta0Words{k}';c=f'coeff0Words{k}';rep=f'(rep256 {size})'
 f.write(f'/-- Dual coefficients for empty bottom fibers in batch{k}. -/\ndef {c} : Nat := {packlit([int(r['duals'][0],16)for r in vals],256)}\n')
 f.write(f'/-- Parallel checks of the empty bottom fibers in batch{k}. -/\ntheorem checked0_{k} :\n    {b} < 2^(256*{size}) ∧ {c} < 2^(256*{size}) ∧\n    {b} &&& ({rep}*(2^90-1)) = {b} ∧\n    batchRecon {rep} {c} bottomDualColumns &&& {b} = 0 ∧\n    rawParityFold 7 ((batchRecon {rep} {c} bottomDualColumns >>> 90) &&& {b}) &&& {rep} = {rep} := by decide +kernel\n')
 f.write(f'/-- Every lane of bottom empty batch{k} has an impossible supported fiber. -/\ntheorem impossible0_{k} (i : Nat) (hi : i < {size}) (word : Nat)\n    (hs : word &&& lane 256 i {b} = word) :\n    reconX bottomEColumns word ≠ reconX bottomRColumns (lane 256 i {b}) := by\n  have h := checked0_{k}\n  rw [rep{size}_eq] at h\n  have hd := batchDual0_raw_sound 256 {size} 90 7 {b} {c} (by decide) (by decide) (by decide)\n    h.1 h.2.1 h.2.2.1 bottomDualColumns bottomDualColumns_bounded h.2.2.2.1 h.2.2.2.2 i hi\n  exact bottomDual_no_supported _ _ _ hd.1 hd.2 hs\n')
f.write('''/-- Size of one cap-zero batch. -/
def size0 (b : Nat) : Nat := if b < 219 then 1024 else 96
/-- Every cap-zero batch has an impossible supported fiber. -/
theorem cap0_batch_impossible (b : Nat) (hb : b < 220) (i : Nat) (hi : i < size0 b) (word : Nat)
    (hs : word &&& lane 256 i (beta0Words b) = word) :
    reconX bottomEColumns word ≠ reconX bottomRColumns (lane 256 i (beta0Words b)) := by
''')
def branch(lo,hi,indent):
 if hi-lo==1:
  f.write(indent+f'have he : b = {lo} := by omega\n'+indent+'subst b\n'+indent+f'exact impossible0_{lo} i hi word hs\n')
 else:
  mid=(lo+hi)//2;f.write(indent+f'by_cases h{mid} : b < {mid}\n'+indent+'·\n');branch(lo,mid,indent+'  ');f.write(indent+'·\n');branch(mid,hi,indent+'  ')
branch(0,220,'  ')
f.write('''/-- Every listed cap-zero bottom fiber is impossible on its shadow support. -/
theorem cap0_no_supported (i : Fin 224352) (word : Nat) (hs : word &&& beta0 i.val = word) :
    reconX bottomEColumns word ≠ reconX bottomRColumns (beta0 i.val) := by
  have hm := Nat.mod_lt i.val (by decide : 0 < 1024)
  have he := Nat.mod_add_div i.val 1024
  have hb : i.val / 1024 < 220 := by omega
  have hl : i.val % 1024 < size0 (i.val / 1024) := by
    unfold size0
    split_ifs <;> omega
  exact cap0_batch_impossible _ hb _ hl word hs
/-- Every cap-zero bottom fiber has no bounded solution with zero outside weight. -/
theorem cap0_impossible (i : Fin 224352) (word : Nat) (hw : word < 2^90)
    (he : reconX bottomEColumns word = reconX bottomRColumns (beta0 i.val))
    (hoff : bitCount 90 (outsideMask word (beta0 i.val)) ≤ 0) : False :=
  cap0_no_supported i word (supported_of_outside_weight_zero 90 word _ hw hoff) he
#print axioms cap0_impossible
end Quantum.Stabilizer.Homological.BB.Bravyi360.BottomEmpty
''');f.close();print('bottomzero',len(rs),'records')

# Data are kernel certificates, not runtime programs.
for _certificate_output in ['Bravyi360BottomEmpty0.lean']:
 _certificate_path = P / _certificate_output
 _certificate_path.write_text(_certificate_path.read_text().replace("\ndef ", "\nnoncomputable def "))
