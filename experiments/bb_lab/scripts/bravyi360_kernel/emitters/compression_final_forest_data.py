import os
from pathlib import Path as _StagePath
_STAGE = _StagePath(os.environ["QEC360_STAGE"])
import pathlib,json,collections,sys
P=pathlib.Path(str(_STAGE / 'qec360-certprobe'));F=json.load(open(P/'gb_frame.json'));E=F['Ecols'];R=F['Rcols'];rows=[sum(((c>>i)&1)<<j for j,c in enumerate(E))for i in range(90)]
def recon(cols,x):
 z=0
 while x:
  l=x&-x;z^=cols[l.bit_length()-1];x^=l
 return z
links={}
for l in open(P/'bottom_top_links.jsonl'):
 r=json.loads(l);links[(int(r['beta'],16),int(r['v0'],16))]=r['topIndex']
bank=collections.defaultdict(lambda:collections.defaultdict(list))
for l in open(P/'fibers_le4.tsv'):
 a,v,_=l.split();b=int(a,16);v=int(v,16);bank[b][v^(v&b)].append(v)
rs=[]
for l in open(P/'nonempty_peel_certificates.jsonl'):
 r=json.loads(l);b=int(r['beta'],16);adj=collections.defaultdict(list);cs=[]
 for q,p,i in r['edges']:adj[q].append((p,i));adj[p].append((q,i))
 for c in r['components']:
  c=int(c,16);root=(c&-c).bit_length()-1;known=1<<root;todo=[root];steps=[]
  for q in todo:
   for p,i in adj[q]:
    if known>>p&1:continue
    known|=rows[i]&b;steps.append(i);todo.append(p)
  assert known==c;cs.append((root,c,recon(E,c),steps))
 ss=[]
 for off,seed in r['seeds']:
  off=int(off,16);seed=int(seed,16);targets=sorted(bank[b][off]);ss.append((off,seed,[(v,links[(b,v)])for v in targets]))
 rs.append((b,b.bit_count(),recon(R,b),cs,ss))
for part,start in enumerate(range(0,len(rs),2000)):
 f=open(P/f'BottomForestPart{part}.lean','w');f.write('import BottomForestCheck\nset_option maxHeartbeats 0\nset_option maxRecDepth 100000\nnamespace Quantum.Stabilizer.Homological.BB.Bravyi360.BottomForest\n')
 for i in range(start,min(start+2000,len(rs)),250):
  f.write(f'/-- Canonical shallow fibers {i} onward. -/\ndef records{i} : List Fiber := [')
  terms=[]
  for b,w,rhs,cs,ss in rs[i:min(i+250,start+2000,len(rs))]:
   comp='['+','.join('⟨'+str(q)+','+hex(c)+','+hex(im)+',['+','.join(map(str,es))+']⟩'for q,c,im,es in cs)+']'
   seeds='['+','.join('⟨'+hex(o)+','+hex(s)+',['+','.join('('+hex(v)+','+str(j)+')'for v,j in ts)+']⟩'for o,s,ts in ss)+']'
   terms.append('⟨'+hex(b)+','+str(w)+','+hex(rhs)+','+comp+','+seeds+'⟩')
  f.write(','.join(terms)+']\n')
  f.write(f'/-- All equations, forest coverage and top-bank joins in this chunk are checked. -/\ntheorem checked{i} : groupCheck records{i} = true := by decide +kernel\n')
 f.write(f'#print axioms checked{start}\nend Quantum.Stabilizer.Homological.BB.Bravyi360.BottomForest\n');f.close()
print('finalforestdata',len(rs),'fibers')
