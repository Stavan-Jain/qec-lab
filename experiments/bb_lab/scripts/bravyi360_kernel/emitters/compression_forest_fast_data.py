import os
from pathlib import Path as _StagePath
_STAGE = _StagePath(os.environ["QEC360_STAGE"])
import pathlib,json,runpy,sys
P=pathlib.Path(str(_STAGE / 'qec360-certprobe'));sys.path.insert(0,str(P));from balance_literals import lit
from compact_literals import packlit
rs=runpy.run_path(str(P/'final_forest_data.py'))['rs']
bank=[int(json.loads(l)['b'],16)for l in open(str(_STAGE / 'all_top_rung_certificates.jsonl'))]
def pack(xs):return packlit(xs,256)
def rep(n):return f'(Nat.div (Nat.sub (Nat.shiftLeft 1 {256*n}) 1) (Nat.sub (Nat.shiftLeft 1 256) 1))'
for part,start in enumerate(range(0,len(rs),2000)):
 s=(P/f'BottomForestPart{part}.lean').read_text().replace('import BottomForestCheck','import BottomForestPacked\nimport PackedDataExpand')
 for i in range(start,min(start+2000,len(rs)),250):
  eps=[];rps=[];js=[]
  for b,w,rhs,cs,ss in rs[i:min(i+250,len(rs))]:
   eps.append((b,0));eps.extend((word,image)for root,word,image,steps in cs);eps.extend((seed,rhs)for off,seed,ts in ss);rps.append((b,rhs));js.extend((v,v^b,bank[j])for off,seed,ts in ss for v,j in ts)
  vals=[rep(len(eps)),pack([a for a,b in eps]),pack([b for a,b in eps]),rep(len(rps)),pack([a for a,b in rps]),pack([b for a,b in rps]),rep(len(js)),pack([a for a,b,c in js]),pack([b for a,b,c in js]),pack([c for a,b,c in js])]
  old=f'theorem checked{i} : groupCheck records{i} = true := by decide +kernel'
  new=f'def packed{i} : PackedGroup := ⟨'+','.join(vals)+f'⟩\n/-- All cached packings and full matrix equations are kernel-checked. -/\ntheorem packed_checked{i} : packedGroupCheck records{i} packed{i} = true := by decide +kernel\n/-- The cached checker certifies the original complete forest predicate. -/\ntheorem checked{i} : groupCheck records{i} = true := packedGroupCheck_sound _ _ packed_checked{i}'
  assert old in s;s=s.replace(old,new)
 (P/f'BottomForestFastPart{part}.lean').write_text(s.replace('\ndef ','\nnoncomputable def '))
s=(P/'BottomForestFastPart0.lean').read_text();s=s[:s.index('/-- Canonical shallow fibers 250 onward.')]+ '#print axioms checked0\nend Quantum.Stabilizer.Homological.BB.Bravyi360.BottomForest\n';(P/'BottomForestFastFirst250.lean').write_text(s)
print('fast forest generated')

# Data are kernel certificates, not runtime programs.
for _certificate_output in ['BottomForestFastPart0.lean', 'BottomForestFastPart1.lean', 'BottomForestFastPart2.lean', 'BottomForestFastPart3.lean', 'BottomForestFastPart4.lean', 'BottomForestFastPart5.lean', 'BottomForestFastPart6.lean', 'BottomForestFastPart7.lean', 'BottomForestFastPart8.lean', 'BottomForestFastPart9.lean', 'BottomForestFastPart10.lean', 'BottomForestFastPart11.lean', 'BottomForestFastPart12.lean', 'BottomForestFastPart13.lean']:
 _certificate_path = P / _certificate_output
 _certificate_path.write_text(_certificate_path.read_text().replace("\ndef ", "\nnoncomputable def "))
