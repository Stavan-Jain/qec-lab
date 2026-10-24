import os
from pathlib import Path as _StagePath
_STAGE = _StagePath(os.environ["QEC360_STAGE"])
import pathlib,json,sys,collections
sys.path.insert(0,str(_STAGE / 'qec360-certprobe'))
from balance_literals import lit
P=pathlib.Path(str(_STAGE / 'qec360-certprobe'));records=[json.loads(l)for l in open(P/'top_empty_dual90_certificates.jsonl')]
rs0=[r for r in records if r['cap']==0];groups=collections.defaultdict(list)
for r in records:
 if r['cap']==1:groups[len(r['duals'])].append(r)
rs1=[r for group in groups.values()for r in group];rs2=[r for r in records if r['cap']==2]
f=open(P/'Bravyi360TopEmptyData.lean','w');f.write('''import PackedDigits
namespace Quantum.Stabilizer.Homological.BB.Bravyi360.TopEmpty
''')
def dispatch(names,lo=0):
 if len(names)==1:return names[0]
 mid=len(names)//2;return f'(if Nat.blt i {lo+mid} then {dispatch(names[:mid],lo)} else {dispatch(names[mid:],lo+mid)})'
for cap,rs in enumerate([rs0,rs1,rs2]):
 names=[]
 for ix,offset in enumerate(range(0,len(rs),1024)):
  b=sum(int(r['beta'],16)<<(512*j)for j,r in enumerate(rs[offset:offset+1024]));name=f'beta{cap}Words{ix}';names.append(name)
  f.write(f'/-- Shared shadow data for empty cap{cap}, chunk{ix}. -/\ndef {name} : Nat := {lit(b)}\n')
 f.write(f'/-- Balanced lookup of the empty cap{cap} shadow batches. -/\ndef beta{cap}Words (i : Nat) : Nat := {dispatch(names)}\n')
 f.write(f'/-- Shadow of empty cap{cap} certificate i, in the certified family order. -/\ndef beta{cap} (i : Nat) : Nat := PackedDigits.lane 512 (i%1024) (beta{cap}Words (i/1024))\n')
 f.write(f'/-- Number of empty cap{cap} certificates. -/\ndef count{cap} : Nat := {len(rs)}\n')
f.write('end Quantum.Stabilizer.Homological.BB.Bravyi360.TopEmpty\n');f.close();print('topdata sizes',len(rs0),len(rs1),len(rs2))

# Data are kernel certificates, not runtime programs.
for _certificate_output in ['Bravyi360TopEmptyData.lean']:
 _certificate_path = P / _certificate_output
 _certificate_path.write_text(_certificate_path.read_text().replace("\ndef ", "\nnoncomputable def "))
