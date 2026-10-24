import os
from pathlib import Path as _StagePath
_STAGE = _StagePath(os.environ["QEC360_STAGE"])
import pathlib,json,runpy,sys
P=pathlib.Path(str(_STAGE / 'qec360-certprobe'));sys.path.insert(0,str(P));from balance_literals import lit
rs=runpy.run_path(str(P/'final_forest_data.py'))['rs'];F=json.load(open(P/'gb_frame.json'));E=F['Ecols'];R=F['Rcols'];cols=[sum(((x>>i)&1)<<j for j,x in enumerate(E))|(sum(((x>>i)&1)<<j for j,x in enumerate(R))<<90)for i in range(45)]
D={int(r['beta'],16):r for l in open(P/'nonempty_dual_covers.jsonl')if (r:=json.loads(l))['cap']==1}
def recon(x):
 z=0
 while x:l=x&-x;z^=cols[l.bit_length()-1];x^=l
 return z
def pack(xs):return sum(x<<(256*i)for i,x in enumerate(xs))
count=0
for part,start in enumerate(range(0,len(rs),2000)):
 f=open(P/f'BottomPositiveOnePart{part}.lean','w');f.write(f'import BottomCap1Cover\nimport BottomForestFastPart{part}\nset_option maxHeartbeats 0\nset_option maxRecDepth 100000\nnamespace Quantum.Stabilizer.Homological.BB.Bravyi360.BottomForest\n')
 for i in range(start,min(start+2000,len(rs)),250):
  covers=[]
  for b,w,rhs,cs,ss in rs[i:min(i+250,len(rs))]:
   ds=[]
   if w==20:
    count+=1;d=D[b];ds=[(int(v,16),recon(int(v,16)))for v in d['duals']];allowed=((1<<90)-1)^b;zero=1
    for c,z in ds:
     assert z&b==0;target=((z>>90)&b).bit_count()%2;allowed&=z if target else~z;zero&=1-target
    assert all(o.bit_count()<=1 for o,s,t in ss)
    assert allowed==sum(o for o,s,t in ss)
    assert zero==int(any(o==0 for o,s,t in ss))
   if w==22:assert any(o==0 for o,s,t in ss)
   covers.append(ds)
  terms=['['+','.join('('+hex(c)+','+hex(z)+')'for c,z in ds)+']'for ds in covers];pairs=[p for ds in covers for p in ds];vals=[pack([1]*len(pairs)),pack([c for c,z in pairs]),pack([z for c,z in pairs])]
  f.write(f'/-- Cap-one dual metadata aligned with forest chunk {i}. -/\nnoncomputable def oneCovers{i} : List (List (Nat × Nat)) := ['+','.join(terms)+']\n')
  f.write(f'/-- Cached dual image packings for chunk {i}. -/\nnoncomputable def onePacked{i} : DualPacking := ⟨'+','.join(lit(x)for x in vals)+'⟩\n')
  f.write(f'/-- Exact cap-one outside coverage and dual images are kernel-checked. -/\ntheorem oneChecked{i} : oneGroupCheck records{i} oneCovers{i} onePacked{i} = true := by decide +kernel\n')
  f.write(f'/-- Every shallow positive shadow has one of the five certified weights. -/\ntheorem weightChecked{i} : (records{i}.all fun r => [14,16,18,20,22].contains r.weight) = true := by decide +kernel\n')
  f.write(f'/-- Every cap-zero record contains the unique empty outside support. -/\ntheorem zeroChecked{i} : (records{i}.all fun r => if Nat.beq r.weight 22 then (r.seeds.map Seed.off).contains 0 else true) = true := by decide +kernel\n')
 f.write(f'#print axioms oneChecked{start}\nend Quantum.Stabilizer.Homological.BB.Bravyi360.BottomForest\n');f.close()
print('positive cap1',count)

# Data are kernel certificates, not runtime programs.
for _certificate_output in ['BottomPositiveOnePart0.lean', 'BottomPositiveOnePart1.lean', 'BottomPositiveOnePart2.lean', 'BottomPositiveOnePart3.lean', 'BottomPositiveOnePart4.lean', 'BottomPositiveOnePart5.lean', 'BottomPositiveOnePart6.lean', 'BottomPositiveOnePart7.lean', 'BottomPositiveOnePart8.lean', 'BottomPositiveOnePart9.lean', 'BottomPositiveOnePart10.lean', 'BottomPositiveOnePart11.lean', 'BottomPositiveOnePart12.lean', 'BottomPositiveOnePart13.lean']:
 _certificate_path = P / _certificate_output
 _certificate_path.write_text(_certificate_path.read_text().replace("\ndef ", "\nnoncomputable def "))
