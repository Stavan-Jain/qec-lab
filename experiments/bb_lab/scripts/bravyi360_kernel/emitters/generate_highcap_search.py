#!/usr/bin/env python3
"""Generate complete bottom high-cap parity-search certificates and canonical forest bridges."""
import json,re,sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parent))
from balanced_nat_literal import literal
import argparse
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--bottom-frame',type=Path,required=True)
parser.add_argument('--forest-certificates',type=Path,required=True)
parser.add_argument('--out',type=Path,required=True)
args=parser.parse_args();OUT=args.out;OUT.mkdir(parents=True,exist_ok=True)
frame=json.loads(args.bottom_frame.read_text());E,R=frame['Ecols'],frame['Rcols']
records=[json.loads(l)for l in args.forest_certificates.open()]
indices={int(r['beta'],16):i for i,r in enumerate(records)}
cases=[dict(beta=r['beta'],cap=(22-int(r['beta'],16).bit_count())//2) for r in records
       if int(r['beta'],16).bit_count() in (14,16)]
assert len(cases)==471
selected=range(len(cases))
def app(cols,x):
 y=0
 while x:
  b=x&-x;x-=b;y^=cols[b.bit_length()-1]
 return y
T=[sum(((v>>q)&1)<<j for j,v in enumerate(E)) for q in range(90)]
header="""import BottomHighSearchSupport
import BottomForestFastPart{part}
set_option maxHeartbeats 0
set_option maxRecDepth 100000
set_option synthInstance.maxSize 100000
namespace Quantum.Stabilizer.Homological.BB.Bravyi360.BottomForest.HighSearch
open F2Mask PackedDigits FiberData

"""
parts=[header.format(part=i)for i in range(14)]
chunks=[[]for _ in range(112)]
metadata=[]
for ix in selected:
 body=''
 r=cases[ix];beta=int(r['beta'],16);cap=r['cap'];index=indices[beta];peel=records[index]
 expected=[int(s[0],16)for s in peel['seeds']]
 basis={}
 for i in range(90):
  if not(beta>>i&1):continue
  v=E[i]
  while v:
   p=v.bit_length()-1
   if p in basis:v^=basis[p]
   else:basis[p]=v;break
 pivots=sorted(basis,reverse=True)
 def proj(x):
  for p in pivots:
   if x>>p&1:x^=basis[p]
  return x
 cols=[proj(c)for c in E];rhs=app(R,beta);target=proj(rhs)
 rowmasks=[sum(((c>>q)&1)<<j for j,c in enumerate(cols))for q in range(90)]
 nodes=0;sols=set();usedrows=set()
 def search(syn,avail,left,selected):
  global nodes
  nodes+=1
  if syn:
   rr=[];x=syn
   while x:
    bit=x&-x;x-=bit;q=bit.bit_length()-1;adj=rowmasks[q]&avail
    if not adj:usedrows.add(q);return f'(.zero {q})'
    rr.append((adj.bit_count(),adj,q))
   rr.sort();used=0;pack=[]
   for _,adj,q in rr:
    if not used&adj:
     pack.append(q);used|=adj
     if len(pack)>left:
      usedrows.update(pack);return '(.pack ['+','.join(map(str,pack))+'])'
   choices=rr[0][1];last=f'(.zero {rr[0][2]})';usedrows.add(rr[0][2])
  else:
   sols.add(selected)
   if left==0:return '.done'
   choices=avail;last='.done'
  branches=[]
  while choices:
   bit=choices&-choices;choices-=bit;avail^=bit;j=bit.bit_length()-1
   branches.append((j,search(syn^cols[j],avail,left-1,selected|bit)))
  for j,yes in reversed(branches):last=f'(.split {j} {yes} {last})'
  return last
 tree=search(target,((1<<90)-1)^beta,cap,0)
 assert sols==set(expected)
 rows=sorted(usedrows);remap={q:i for i,q in enumerate(rows)}
 tree=re.sub(r'\(\.zero (\d+)\)',lambda m:f'(.zero {remap[int(m[1])]})',tree)
 tree=re.sub(r'\(\.pack \[([0-9,]*)\]\)',lambda m:'(.pack ['+','.join(str(remap[int(v)])for v in m[1].split(','))+'])',tree)
 lp=[]
 for q in rows:
  lam=sum(((proj(1<<j)>>q)&1)<<j for j in range(90))
  assert app(T,lam)==rowmasks[q] and app(T,lam)&beta==0
  lp.append(lam)
 newtarget=sum(((target>>q)&1)<<i for i,q in enumerate(rows))
 assert newtarget==sum(((lam&rhs).bit_count()%2)<<i for i,lam in enumerate(lp))
 rp=sum(rowmasks[q]<<(256*i)for i,q in enumerate(rows));lp=sum(l<<(256*i)for i,l in enumerate(lp))
 for name,value,desc in [('rows',rp,'Packed projected parity rows'),('lambdas',lp,'Dual witnesses for the projected rows')]:
  body+=f'/-- {desc} for high-cap case {ix}. -/\nnoncomputable def {name}{ix} : Nat := {literal(value,4096)}\n\n'
 body+=f'/-- Exact seed outside supports for high-cap case {ix}. -/\nnoncomputable def allowed{ix} : List Nat := ['+','.join(map(hex,expected))+']\n\n'
 body+=f'/-- Exhaustive parity-search certificate for high-cap case {ix}. -/\nnoncomputable def tree{ix} : ParitySearchTree := {tree}\n\n'
 body+=f'/-- Every search branch for high-cap case {ix} is kernel-checked. -/\ntheorem checked{ix} : checkParitySearch 8 (fun q => lane 256 q rows{ix}) {hex(newtarget)} (allowed{ix}.contains) tree{ix} {hex(((1<<90)-1)^beta)} 0 {cap} = true := by decide +kernel\n\n'
 body+=f'/-- Initial projected duals for high-cap case {ix} satisfy the full bottom equations. -/\ntheorem initial{ix} : checkSearchInitial 8 90 {len(rows)} rows{ix} lambdas{ix} {hex(newtarget)} {hex(beta)} {hex(rhs)} bottomETransposeColumns = true := by decide +kernel\n\n'
 body+=f'/-- The actual bottom carry gives the certified right-hand side for case {ix}. -/\ntheorem rhs{ix} : reconX bottomRColumns {hex(beta)} = {hex(rhs)} := by decide +kernel\n\n'
 body+=f'''/-- Every sufficiently light solution in high-cap case {ix} has a listed outside support. -/
theorem complete{ix} {{word : Nat}} (hw : word < 2^90)
    (he : reconX bottomEColumns word = reconX bottomRColumns {hex(beta)})
    (hoff : bitCount 90 (outsideMask word {hex(beta)}) ≤ {cap}) :
    outsideMask word {hex(beta)} ∈ allowed{ix} := by
  exact outside_mem_of_search 8 90 {len(rows)} rows{ix} lambdas{ix} {hex(newtarget)} {hex(beta)} {hex(rhs)} {cap}
    bottomEColumns bottomETransposeColumns tree{ix} allowed{ix} (by decide)
    (fun lam v _ _ => bottomETranspose_parityFold lam v) transpose_bounded initial{ix} checked{ix}
    hw (he.trans rhs{ix}) hoff

'''
 parts[index//2000]+=body
 chunks[index//250].append((ix,beta,cap))
 metadata.append(dict(index=ix,fiber_index=index,beta=hex(beta),cap=cap,nodes=nodes,used_rows=len(rows),solutions=len(expected)))

for chunk,entries in enumerate(chunks):
 start=chunk*250;part=start//2000
 conds=[f'(r.beta = {hex(beta)} ∧ outsideBudget r = {cap} ∧ r.seeds.map Seed.off = allowed{ix})'for ix,beta,cap in entries]
 condition=' ∨ '.join(conds)or'False'
 parts[part]+=f"""/-- Every high-cap record in canonical chunk {chunk} matches a certified search. -/
theorem matches{start} : records{start}.all (fun r => decide (
    (r.weight = 16 ∨ r.weight = 14) → {condition})) = true := by decide +kernel

/-- The checked searches give the exact seed-support list for canonical chunk {chunk}. -/
theorem chunk_complete{start} (r : Fiber) (hr : r ∈ records{start})
    (hweight : r.weight = 16 ∨ r.weight = 14) (word : Nat) (hw : word < 2^90)
    (he : reconX bottomEColumns word = reconX bottomRColumns r.beta)
    (hoff : bitCount 90 (outsideMask word r.beta) ≤ outsideBudget r) :
    outsideMask word r.beta ∈ r.seeds.map Seed.off := by
  have h := (of_decide_eq_true ((List.all_eq_true.mp matches{start}) r hr)) hweight
"""
 if len(entries)==1:
  ix,beta,cap=entries[0]
  parts[part]+=f'  exact search_census_transfer r {hex(beta)} {cap} allowed{ix} complete{ix} h word hw he hoff\n'
 elif entries:
  parts[part]+='  rcases h with '+' | '.join('h'+str(j)for j in range(len(entries)))+'\n'
  for j,(ix,beta,cap)in enumerate(entries):
   parts[part]+=f'  · exact search_census_transfer r {hex(beta)} {cap} allowed{ix} complete{ix} h{j} word hw he hoff\n'
 else:parts[part]+='  exact False.elim h\n'
 parts[part]+='\n'
for part,body in enumerate(parts):
 body+=f'#print axioms chunk_complete{part*2000}\nend Quantum.Stabilizer.Homological.BB.Bravyi360.BottomForest.HighSearch\n'
 # Empty canonical chunks close by contradiction; their solution premises are unused.
 def clean_empty(m):
  block=m.group(0)
  return re.sub(r'\b(hw|he|hoff)\b',r'_\1',block) if 'False.elim h' in block else block
 body=re.sub(r'theorem chunk_complete\d+\b.*?(?=\n\n)',clean_empty,body,flags=re.S)
 (OUT/f'BottomHighSearchPart{part}.lean').write_text(body)
body='import Bravyi360BottomForest\n'+''.join(f'import BottomHighSearchPart{i}\n'for i in range(14))
body+='set_option maxHeartbeats 0\nnamespace Quantum.Stabilizer.Homological.BB.Bravyi360.BottomForest\nopen F2Mask FiberData\n/-- Every high-cap record in the canonical forest has a checked outside-support census. -/\ntheorem highRecords_complete (b : Fin 112) (r : Fiber) (hr : r ∈ recordsAt b.val)\n    (hweight : r.weight = 16 ∨ r.weight = 14) (word : Nat) (hw : word < 2^90)\n    (he : reconX bottomEColumns word = reconX bottomRColumns r.beta)\n    (hoff : bitCount 90 (outsideMask word r.beta) ≤ outsideBudget r) :\n    outsideMask word r.beta ∈ r.seeds.map Seed.off := by\n  fin_cases b\n'
for chunk in range(112):
 start=chunk*250
 body+=f'  · change r ∈ records{start} at hr\n    exact HighSearch.chunk_complete{start} r hr hweight word hw he hoff\n'
body+='/-- Every canonical high-cap fiber contains the outside support of each actual solution. -/\ntheorem highFiberAt_complete (j : Fin 27867)\n    (hweight : (fiberAt j.val).weight = 16 ∨ (fiberAt j.val).weight = 14)\n    (word : Nat) (hw : word < 2^90)\n    (he : reconX bottomEColumns word = reconX bottomRColumns (fiberAt j.val).beta)\n    (hoff : bitCount 90 (outsideMask word (fiberAt j.val).beta) ≤ outsideBudget (fiberAt j.val)) :\n    outsideMask word (fiberAt j.val).beta ∈ (fiberAt j.val).seeds.map Seed.off :=\n  highRecords_complete ⟨j.val/250,by omega⟩ (fiberAt j.val) (fiberAt_mem j)\n    hweight word hw he hoff\n#print axioms highFiberAt_complete\nend Quantum.Stabilizer.Homological.BB.Bravyi360.BottomForest\n'
(OUT/'Bravyi360BottomHighSearchComplete.lean').write_text(body)
(OUT/'high_search_metadata.json').write_text(json.dumps(metadata,indent=2)+'\n')
print(json.dumps(dict(cases=len(metadata),nodes=sum(r['nodes']for r in metadata),parts=14)))
