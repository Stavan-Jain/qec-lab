import os
from pathlib import Path as _StagePath
_STAGE = _StagePath(os.environ["QEC360_STAGE"])
import pathlib,json,sys
sys.path.insert(0,str(_STAGE / 'qec360-certprobe'))
from balance_literals import lit
P=pathlib.Path(str(_STAGE / 'qec360-certprobe'));rs=[int(json.loads(l)['b'],16)for l in open(str(_STAGE / 'all_top_rung_certificates.jsonl'))];B=1024;W=180
assert len(rs)==58130 and len(set(rs))==len(rs)
f=open(P/'Bravyi360TopBank.lean','w');f.write('import Mathlib.Data.Nat.Bitwise\nnamespace Quantum.Stabilizer.Homological.BB.Bravyi360\n')
chunks=[]
for ix,start in enumerate(range(0,len(rs),B)):
 xs=rs[start:start+B];v=sum(x<<(W*i)for i,x in enumerate(xs));name=f'topBankChunk{ix}';chunks.append(name)
 f.write(f'/-- Consecutive middle shadows, starting at index {start}. -/\ndef {name} : Nat := {lit(v)}\n')
def tree(lo,hi):
 if hi-lo==1:return chunks[lo]
 mid=(lo+hi)//2
 return f'(if Nat.blt i {mid} then {tree(lo,mid)} else {tree(mid,hi)})'
f.write('/-- Read one chunk of the common middle-shadow bank in logarithmic dispatch depth. -/\ndef topBankChunk (i : Nat) : Nat :=\n  if Nat.blt i '+str(len(chunks))+' then '+tree(0,len(chunks))+' else 0\n')
f.write('''/-- The common bank has one entry for every classified top rung. -/
def topBankSize : Nat := 58130
/-- Read a middle shadow from the common top-rung bank. -/
def topBankMask (i : Nat) : Nat :=
  Nat.land (Nat.shiftRight (topBankChunk (Nat.div i 1024)) (Nat.mul 180 (Nat.mod i 1024)))
    (2^180-1)
/-- Every common-bank entry fits in the middle code's coordinate space. -/
theorem topBankMask_lt (i : Nat) : topBankMask i < 2^180 := by
  apply lt_of_le_of_lt Nat.and_le_right
  omega
end Quantum.Stabilizer.Homological.BB.Bravyi360
''');f.close();print(len(rs),'records',len(chunks),'chunks')

# Data are kernel certificates, not runtime programs.
for _certificate_output in ['Bravyi360TopBank.lean']:
 _certificate_path = P / _certificate_output
 _certificate_path.write_text(_certificate_path.read_text().replace("\ndef ", "\nnoncomputable def "))
