import os
from pathlib import Path as _StagePath
_STAGE = _StagePath(os.environ["QEC360_STAGE"])
import pathlib,json,sys
sys.path.insert(0,str(_STAGE / 'qec360-certprobe'))
from balance_literals import lit
from compact_literals import packlit
P=pathlib.Path(str(_STAGE / 'qec360-certprobe'))
f=open(P/'Bravyi360BottomEmptyData.lean','w');f.write('import PackedDigits\nimport PackedDataExpand\nset_option maxRecDepth 100000\nset_option maxHeartbeats 0\nnamespace Quantum.Stabilizer.Homological.BB.Bravyi360.BottomEmpty\n')
def dispatch(names,lo=0):
 if len(names)==1:return names[0]
 mid=len(names)//2;return f'(if Nat.blt i {lo+mid} then {dispatch(names[:mid],lo)} else {dispatch(names[mid:],lo+mid)})'
for cap,fn in enumerate(['emptycap45_short5.tsv','emptycap45_1.tsv','emptycap45_2.tsv']):
 rs=[int(l.split()[0],16)for l in open(P/fn)];names=[]
 for ix,off in enumerate(range(0,len(rs),1024)):
  b=sum(x<<(256*j)for j,x in enumerate(rs[off:off+1024]));name=f'beta{cap}Words{ix}';names.append(name)
  f.write(f'/-- Shared shadow data for empty cap{cap}, chunk{ix}. -/\ndef {name} : Nat := {packlit(rs[off:off+1024],256)}\n')
 f.write(f'/-- Balanced lookup of empty cap{cap} shadow chunks. -/\ndef beta{cap}Words (i : Nat) : Nat := {dispatch(names)}\n')
 f.write(f'/-- Shadow at position i in the certified empty cap{cap} family. -/\ndef beta{cap} (i : Nat) : Nat := PackedDigits.lane 256 (i%1024) (beta{cap}Words (i/1024))\n')
 f.write(f'/-- Number of empty cap{cap} certificates. -/\ndef count{cap} : Nat := {len(rs)}\n')
f.write('/-- Repetition mask for the 256-bit bottom certificate records. -/\ndef rep256 (size : Nat) : Nat := ((1 <<< (256*size))-1)/((1 <<< 256)-1)\n')
for size in [1024,96,424,859]:
 f.write(f'/-- Kernel-checked repetition mask for {size} records. -/\ntheorem rep{size}_eq : rep256 {size} = PackedDigits.packF 256 {size} (fun _ => 1) := by decide +kernel\n')
f.write('end Quantum.Stabilizer.Homological.BB.Bravyi360.BottomEmpty\n');f.close()

# Data are kernel certificates, not runtime programs.
for _certificate_output in ['Bravyi360BottomEmptyData.lean']:
 _certificate_path = P / _certificate_output
 _certificate_path.write_text(_certificate_path.read_text().replace("\ndef ", "\nnoncomputable def "))
