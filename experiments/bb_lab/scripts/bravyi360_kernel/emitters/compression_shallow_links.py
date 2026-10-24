import os, pathlib
P=pathlib.Path(os.environ['QEC360_STAGE'])/'qec360-certprobe';N=274011;B=1024;C=268
header='''import Bravyi360BottomCaseData
import Bravyi360BottomBetaCache
import Bravyi360BottomEmptyData
namespace Quantum.Stabilizer.Homological.BB.Bravyi360.BottomShallow
open F2Mask PackedDigits
set_option maxRecDepth 100000
set_option maxHeartbeats 0
'''
f=open(P/'BottomShallowLinkCheck.lean','w');f.write(header+'/-- The forest association uses the already checked packed shadow cache. -/\nnoncomputable def forestLinkCheck (beta index : Nat) : Bool :=\n  Nat.blt index 27867 && Nat.beq beta (BottomForest.betaAtFast index)\n/-- A packed forest association has the original canonical-record meaning. -/\ntheorem forestLinkCheck_eq_true (beta index : Nat) :\n    forestLinkCheck beta index = true ↔\n      index < 27867 ∧ beta = (BottomForest.fiberAt index).beta := by\n  simp only [forestLinkCheck,Bool.and_eq_true,Nat.blt_eq,Nat.beq_eq]\n  constructor\n  · rintro ⟨hi,hb⟩\n    exact ⟨hi,by simpa only [BottomForest.betaAtFast_eq ⟨index,hi⟩] using hb⟩\n  · rintro ⟨hi,hb⟩\n    exact ⟨hi,by simpa only [BottomForest.betaAtFast_eq ⟨index,hi⟩] using hb⟩\n'+'''/-- A bottom census representative is identified with its precise non-deep certificate family. -/
def linkCheck (i : Nat) : Bool :=
  let b := OrbitBatchData.repAt i
  let j := BottomCase.localIndex i
  match BottomCase.tag i with
  | 0 => Nat.blt j 224352 && Nat.beq b (BottomEmpty.beta0 j)
  | 1 => Nat.blt j 20904 && Nat.beq b (BottomEmpty.beta1 j)
  | 2 => Nat.blt j 859 && Nat.beq b (BottomEmpty.beta2 j)
  | 3 => forestLinkCheck b j
  | 4 => Nat.beq b 0
  | 5 => true
  | _ => false

/-- Direct interpretation of a checked shallow-family association. -/
theorem linkCheck_sound (i : Nat) (h : linkCheck i = true) :
    (BottomCase.tag i = 0 ∧ BottomCase.localIndex i < 224352 ∧
      OrbitBatchData.repAt i = BottomEmpty.beta0 (BottomCase.localIndex i)) ∨
    (BottomCase.tag i = 1 ∧ BottomCase.localIndex i < 20904 ∧
      OrbitBatchData.repAt i = BottomEmpty.beta1 (BottomCase.localIndex i)) ∨
    (BottomCase.tag i = 2 ∧ BottomCase.localIndex i < 859 ∧
      OrbitBatchData.repAt i = BottomEmpty.beta2 (BottomCase.localIndex i)) ∨
    (BottomCase.tag i = 3 ∧ BottomCase.localIndex i < 27867 ∧
      OrbitBatchData.repAt i = (BottomForest.fiberAt (BottomCase.localIndex i)).beta) ∨
    (BottomCase.tag i = 4 ∧ OrbitBatchData.repAt i = 0) ∨ BottomCase.tag i = 5 := by
  unfold linkCheck at h
  split at h <;> simp_all [forestLinkCheck_eq_true,Bool.and_eq_true,Nat.blt_eq,Nat.beq_eq]
end Quantum.Stabilizer.Homological.BB.Bravyi360.BottomShallow
''');f.close()
import json
rows = [json.loads(line) for line in open(P/'bottom_partition.jsonl')]
for part,start in enumerate(range(0,C,32)):
 f=open(P/f'BottomShallowLinkPart{part}.lean','w');f.write('import BottomShallowCacheCheck\nnamespace Quantum.Stabilizer.Homological.BB.Bravyi360.BottomShallow\nset_option maxRecDepth 100000\nset_option maxHeartbeats 0\n')
 for b in range(start,min(start+32,C)):
  size=min(B,N-b*B)
  for family,chunk,count in [(0,1024,220),(1,1024,21),(2,1024,1),(3,250,112)]:
   indices=[r['familyIndex'] for r in rows[B*b:B*b+size] if r['family']==family]
   base=min(indices)//chunk if indices else 0
   assert not indices or max(indices)//chunk<=base+1
   if family<3: name=lambda x:f'BottomEmpty.beta{family}Words{min(x,count-1)}'
   else: name=lambda x:f'BottomForest.packed{250*min(x,count-1)}.rInput'
   f.write(f'/-- Adjacent family {family} data blocks for census chunk {b}. -/\nnoncomputable def cache{family}_{b} : LinkCache := ⟨{base},{name(base)},{name(base+1)}⟩\n')
  f.write(f'/-- Checked exact family associations using fixed data blocks for chunk {b}. -/\ntheorem cachedLinkChecked{b} : ((List.range {size}).all fun j => linkCachedCheck OrbitBatchData.words{b} BottomCase.words{b} {B*b} cache0_{b} cache1_{b} cache2_{b} cache3_{b} j) = true := by decide +kernel\n')
  f.write(f'/-- Exact certificate-family associations for census chunk {b}. -/\ntheorem linkChecked{b} : ((List.range {size}).all fun j => linkCheck ({B*b}+j)) = true := by\n  apply List.all_eq_true.mpr\n  intro row hrow\n  have h := List.all_eq_true.mp cachedLinkChecked{b} row hrow\n  have hc := linkCachedCheck_sound OrbitBatchData.words{b} BottomCase.words{b} {B*b} row\n    cache0_{b} cache1_{b} cache2_{b} cache3_{b} (by rfl) (by rfl) (by rfl) (by rfl)\n    (by rfl) (by rfl) (by rfl) (by rfl) h\n  have hcFixed : linkChunkCheck (OrbitBatchData.wordsAt {b}) (BottomCase.wordsAt {b}) ({B}*{b}) row = true := hc\n  rw [linkChunkCheck_eq {b} row (by have := List.mem_range.mp hrow; omega)] at hcFixed\n  exact hcFixed\n')
 f.write(f'#print axioms linkChecked{start}\nend Quantum.Stabilizer.Homological.BB.Bravyi360.BottomShallow\n');f.close()
f=open(P/'BottomShallowLinks.lean','w');f.write(''.join(f'import BottomShallowLinkPart{k}\n'for k in range(9))+'namespace Quantum.Stabilizer.Homological.BB.Bravyi360.BottomShallow\nset_option maxRecDepth 100000\nset_option maxHeartbeats 0\n')
f.write('''/-- Every valid census entry has its exact checked shallow-family association. -/
theorem link_checked (i : Fin 274011) : linkCheck i.val = true := by
  have hm := Nat.mod_lt i.val (by decide : 0 < 1024)
  have he := Nat.mod_add_div i.val 1024
  have hb : i.val / 1024 < 268 := by omega
  generalize hn : i.val / 1024 = b at hb
''')
def branch(lo,hi,indent):
 if hi-lo==1:
  size=min(B,N-lo*B)
  f.write(indent+f'have heq : b = {lo} := by omega\n'+indent+'subst b\n'+indent+f'have hj : i.val%1024 ∈ List.range {size} := List.mem_range.mpr (by omega)\n'+indent+f'have h := List.all_eq_true.mp linkChecked{lo} _ hj\n'+indent+f'have heq : {B*lo}+i.val%1024 = i.val := by omega\n'+indent+'simpa only [heq] using h\n')
 else:
  mid=(lo+hi)//2;f.write(indent+f'by_cases h{mid} : b < {mid}\n'+indent+'·\n');branch(lo,mid,indent+'  ');f.write(indent+'·\n');branch(mid,hi,indent+'  ')
branch(0,C,'  ')
f.write('#print axioms link_checked\nend Quantum.Stabilizer.Homological.BB.Bravyi360.BottomShallow\n');f.close()

# Data are kernel certificates, not runtime programs.
for _certificate_output in ['BottomShallowLinkCheck.lean', 'BottomShallowLinkPart0.lean', 'BottomShallowLinkPart1.lean', 'BottomShallowLinkPart2.lean', 'BottomShallowLinkPart3.lean', 'BottomShallowLinkPart4.lean', 'BottomShallowLinkPart5.lean', 'BottomShallowLinkPart6.lean', 'BottomShallowLinkPart7.lean', 'BottomShallowLinkPart8.lean', 'BottomShallowLinks.lean']:
 _certificate_path = P / _certificate_output
 _certificate_path.write_text(_certificate_path.read_text().replace("\ndef ", "\nnoncomputable def "))
