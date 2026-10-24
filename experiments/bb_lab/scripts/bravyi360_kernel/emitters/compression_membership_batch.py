import os
from pathlib import Path as _StagePath
_STAGE = _StagePath(os.environ["QEC360_STAGE"])
import json,pathlib,sys
sys.path.insert(0,str(_STAGE / 'qec360-certprobe'))
from balance_literals import lit
P=pathlib.Path(str(_STAGE / 'qec360-certprobe'));rs=[json.loads(l)for l in open(P/'gb_orbit_membership.jsonl')];src=[json.loads(l)['mask']for l in open(P/'gb_orbit_certificates.jsonl')];G=[int(x,16) if isinstance(x,str) else x for x in json.load(open(str(_STAGE / 'tower_direct43_cert.json')))['G43Rows']];assert [r['mask']for r in rs]==src
W=128;B=1024
prefix='''import OrbitBatchData
import Bravyi360ReachabilityData
import BatchLinear
set_option maxRecDepth 100000
set_option maxHeartbeats 0
namespace Quantum.Stabilizer.Homological.BB.Bravyi360.Membership
open F2Mask PackedDigits ReachabilityData
/-- The common generators and coefficient positions fit in the packed membership lanes. -/
theorem generator_bound : ∀ p ∈ reachableGenerators, p.1 < 128 ∧ p.2 < 2^128 := by decide +kernel
'''
f=open(P/'Bravyi360Membership.lean','w');f.write(prefix)
for ix,start in enumerate(range(0,len(rs),B)):
 cs=rs[start:start+B];size=len(cs);coeff=sum(r['coeff43']<<(W*i)for i,r in enumerate(cs));beta=sum(r['mask']<<(W*i)for i,r in enumerate(cs));rep=sum(1<<(W*i)for i in range(size));v=0
 for j,col in enumerate(G):v^=((coeff>>j)&rep)*col
 assert v==beta
 f.write(f'/-- Reconstruction coefficients for census batch {ix}. -/\ndef coeff{ix} : Nat := {lit(coeff)}\n')
 f.write(f'/-- Coefficients occupy only this batch. -/\ntheorem coeff{ix}_bound : coeff{ix} < 2^(128*{size}) := by decide +kernel\n')
 f.write(f'/-- The packed coefficient table reconstructs this shared census batch. -/\ntheorem checked{ix} : batchRecon (OrbitBatchData.repeat128 {size}) coeff{ix} reachableGenerators = OrbitBatchData.words{ix} := by decide +kernel\n')
 # REP equality supplied later: same helper used by orbit.
 f.write(f'/-- The packed reconstruction certificate certifies every lane in this batch. -/\ntheorem member{ix} :\n    ∀ i < {size}, XorSpan (reachableGenerators.map Prod.snd) (lane 128 i OrbitBatchData.words{ix}) := by\n  apply batchRecon_mem_span 128 {size} coeff{ix} OrbitBatchData.words{ix} (by decide) coeff{ix}_bound reachableGenerators generator_bound\n  rw [← OrbitBatchData.repeat{size}_eq]\n  exact checked{ix}\n')
f.write('/-- Every batch has a checked reconstruction in the common reachable span. -/\ntheorem batch_mem (b : Nat) (hb : b < 268) : ∀ i < OrbitBatchData.sizeAt b, XorSpan (reachableGenerators.map Prod.snd) (lane 128 i (OrbitBatchData.wordsAt b)) := by\n')
def branch(lo,hi,indent):
 if hi-lo==1:
  f.write(indent+f'have he : b = {lo} := by omega\n'+indent+'subst b\n'+indent+f'exact member{lo}\n')
 else:
  mid=(lo+hi)//2
  f.write(indent+f'by_cases h{mid} : b < {mid}\n')
  f.write(indent+'·\n'); branch(lo,mid,indent+'  ')
  f.write(indent+'·\n'); branch(mid,hi,indent+'  ')
branch(0,268,'  ')
f.write('/-- Every listed orbit representative belongs to the reachable43-generator span. -/\ntheorem repAt_mem (i : Fin 274011) : XorSpan (reachableGenerators.map Prod.snd) (OrbitBatchData.repAt i.val) := by\n  have hm := Nat.mod_lt i.val (by decide : 0 < 1024)\n  have he := Nat.mod_add_div i.val 1024\n  have hb : i.val / 1024 < 268 := by omega\n  have hl : i.val % 1024 < OrbitBatchData.sizeAt (i.val / 1024) := by\n    unfold OrbitBatchData.sizeAt\n    split_ifs <;> omega\n  exact batch_mem (i.val / 1024) hb (i.val % 1024) hl\n')
f.write('#print axioms member0\nend Quantum.Stabilizer.Homological.BB.Bravyi360.Membership\n');f.close();print(len(rs),'members',len(range(0,len(rs),B)),'batches')

# Data are kernel certificates, not runtime programs.
for _certificate_output in ['Bravyi360Membership.lean']:
 _certificate_path = P / _certificate_output
 _certificate_path.write_text(_certificate_path.read_text().replace("\ndef ", "\nnoncomputable def "))
