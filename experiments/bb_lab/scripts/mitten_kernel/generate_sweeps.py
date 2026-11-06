#!/usr/bin/env python3
"""Generate complete split certificates checked by Lean's packed arithmetic.

Regenerate all classifier modules from the exported FloorData tables with::

    python3 experiments/bb_lab/scripts/mitten_kernel/generate_sweeps.py \\
        --root /path/to/QECLean --all --force

The Python computations only propose certificate data.  Lean checks the exact
mask coverage, linear maps, packed weight flags, and every surviving candidate.
"""
import argparse, ast, hashlib, itertools, json, re, time
from pathlib import Path

NS='Quantum.Stabilizer.Homological.LP.M150'
MOD='QEC.Stabilizer.Codes.Mitten.M150.'
BANNER='''/-
GENERATED FILE — DO NOT HAND-EDIT
Generator: qec-lab:experiments/bb_lab/scripts/mitten_kernel/generate_sweeps.py
Data source: QEC/Stabilizer/Codes/Mitten/M150/FloorData.lean
Regenerate: python3 experiments/bb_lab/scripts/mitten_kernel/generate_sweeps.py
  --root /path/to/QECLean --all --force
All candidate coverage and packed computations are checked by Lean's kernel.
-/
'''
def header(imports,packed=True):
 opens='PackedDigits F2Mask' if packed else 'F2Mask'
 return BANNER+'\n'+''.join('import '+MOD+x+'\n' for x in imports)+'\nnamespace '+NS+'\nopen '+opens+'\nset_option maxRecDepth 100000\n\n'
def decl(name,typ,value,doc='Generated certificate data.'):
 return 'set_option maxHeartbeats 0 in\n-- Elaborating this generated certificate can exceed the default heartbeat budget.\n'+f'/-- {doc} -/\ndef {name} : {typ} := {value}\n\n'
def theorem(name,typ,proof):
 return 'set_option maxHeartbeats 0 in\n-- Elaborating this generated certificate can exceed the default heartbeat budget.\n'+f'/-- Kernel-checked certificate for {name}. -/\ntheorem {name} : {typ} := by\n{proof}\n\n'
def natlist(vals):
 if not vals:return '[]'
 size=5 if isinstance(vals[0],tuple) else 6
 return '[\n'+''.join('  '+', '.join(str(v) for v in vals[i:i+size])+(',' if i+size<len(vals) else '')+'\n' for i in range(0,len(vals),size))+']'
def cases(count):
 return ' |\n    '.join(' | '.join('rfl' for _ in range(min(10,count-i)))
   for i in range(0,count,10))
def xf(tbl,m):
 v=0
 while m:
  j=(m&-m).bit_length()-1
  if j<len(tbl):v^=tbl[j]
  m&=m-1
 return v

def weights(split):
 p,q,r,m=split
 return (q if m==2 else p), (r if m in [0,2] else q), [q,r,p,r][m]

def chunk_caps(outer_weight,inner_weight,outer_cap,inner_cap):
 return (32768 if inner_weight==0 else 4096 if inner_weight==1 else outer_cap), \
   (32768 if outer_weight==0 else 4096 if outer_weight==1 else inner_cap)

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--root',type=Path,required=True);ap.add_argument('--instance',default='X0');ap.add_argument('--split',type=int,default=69);ap.add_argument('--chunk',type=int,default=512);ap.add_argument('--outer-chunk',type=int,default=256);ap.add_argument('--all',action='store_true');ap.add_argument('--force',action='store_true');args=ap.parse_args()
 out=args.root/'QEC/Stabilizer/Codes/Mitten/M150';s=(out/'FloorData.lean').read_text()
 if not args.force and any((out/name).exists() for name in ['MaskChunks.lean','SweepMaps.lean']):
  ap.error('generated data already exists; pass --force to regenerate')
 lists={k:ast.literal_eval(v) for k,v in re.findall(r'def (\w+) : List Nat := (\[[\s\S]*?\])',s)}
 for key,parts in re.findall(r'def (\w+) : List Nat := (\w+(?: \+\+ \w+)+)',s):
  lists[key]=[v for part in parts.split(' ++ ') for v in lists[part]]
 tables={};splits={}
 for n in ['X0','X1','Z0','Z1']:
  b=re.search(r'def m4T'+n+r' : SweepTables where\n([\s\S]*?)\n/--',s).group(1)
  fields={}
  for k,v in re.findall(r'  (\w+) := (\w+|\[\])',b):fields[k]=lists[v] if v in lists else int(v) if v.isdigit() else []
  tables[n]=fields
  splits[n]=ast.literal_eval(re.search(r'def splits'+n+r' : List \(Nat × Nat × Nat × Nat\) := (\[[\s\S]*?\])',s).group(1))
 def images(n,m):
  t=tables[n]
  if m==0:return t['tUW'],t['tTW'],[0]
  if m==1:return t['tUT'],t['tWT'],[0]
  if m==2:return [xf(t['tP'],c) for c in t['tWC']],[xf(t['tP'],c) for c in t['tTC']],t['ann']
  return [xf(t['tPT'],c) for c in t['tUC']],[xf(t['tPT'],c) for c in t['tWC2']],t['annT']
 masks={k:sorted(sum(1<<i for i in c) for c in itertools.combinations(range(30),k)) for k in range(5)}
 outer_needed=set();inner_needed=set()
 for family in splits.values():
  for split in family:
   ow,iw,_=weights(split);ocap,icap=chunk_caps(ow,iw,args.outer_chunk,args.chunk)
   outer_needed.add((ow,ocap));inner_needed.add((iw,icap))
 chunks={(k,cap):[masks[k][i:i+cap] for i in range(0,len(masks[k]),cap)]
   for k,cap in sorted(inner_needed)}
 outer_chunks={(k,cap):[masks[k][i:i+cap] for i in range(0,len(masks[k]),cap)]
   for k,cap in sorted(outer_needed)}
 md=header(['SweepAdapters'],packed=False)
 base_for_weight={};family_for_cap={}
 for prefix,chunk_map in [('masks',chunks),('outerMasks',outer_chunks)]:
  for (k,cap),cs in chunk_map.items():
   family=f'{prefix}{k}S{cap}'
   previous=family_for_cap.get((k,cap));base=base_for_weight.get(k)
   for i,vs in enumerate(cs):
    name=f'{family}Chunk{i}'
    if previous:
     value=f'{previous}Chunk{i}'
    elif cap==32768 and base:
     value=f'{base}Chunks.flatten'
    else:
     value=natlist(vs)
    md+=decl(name,'List Nat',value,'An ascending block of fixed-weight thirty-bit masks.')
   names=',\n'.join('  '+', '.join(f'{family}Chunk{i}' for i in range(j,min(j+3,len(cs))))
     for j in range(0,len(cs),3))
   md+=decl(f'{family}Chunks','List (List Nat)','[\n'+names+'\n]',
     'Blocks covering all masks of the specified weight.')
   proof='  decide +kernel' if base is None else f'''  calc
    {family}Chunks.flatten = {base}Chunks.flatten := by decide +kernel
    _ = masksOfWt 30 {k} := {base}_cover'''
   md+=theorem(f'{family}_cover',f'{family}Chunks.flatten = masksOfWt 30 {k}',proof)
   for i in range(len(cs)):
    name=f'{family}Chunk{i}'
    member='List.mem_cons_self'
    for _ in range(i):member='List.mem_cons_of_mem _ (\n      '+member+')'
    md+=theorem(name+'_bound',f'∀ x ∈ {name}, x < 2 ^ 30',f'''  intro x hx
  apply masksOfWt_bound 30 {k} x
  rw [← {family}_cover]
  exact List.mem_flatten.mpr ⟨{name},
    {member}, hx⟩''')
   base_for_weight.setdefault(k,family);family_for_cap.setdefault((k,cap),family)
 (out/'MaskChunks.lean').write_text(md+'end '+NS+'\n')
 sm=header(['FloorData','PackedClassification'])
 for n in tables:
  for m in sorted(set(t[3] for t in splits[n])):
   l,r,ann=images(n,m)
   for side,cols in [('Left',l),('Right',r)]:
    name=f'{n}M{m}{side}'
    sm+=decl(name+'Columns','List (Nat × Nat)','[\n'+',\n'.join(f'  ({i}, {v})' for i,v in enumerate(cols))+'\n]','Columns of the split reconstruction map.')
    sm+=theorem(name+'_columns_bound',f'∀ p ∈ {name}Columns, p.1 < 32 ∧ p.2 < 2 ^ 30','  decide +kernel')
    sm+=theorem(name+'_eq',f'∀ x < 2 ^ 30, reconX {name}Columns x = sweep{side} m4T{n} {m} x',f'  apply sweep{side}_eq_recon\n  decide +kernel')
    sm+=theorem(name+'_bound',f'∀ x < 2 ^ 30, sweep{side} m4T{n} {m} x < 2 ^ 30',f'  intro x hx\n  rw [← {name}_eq x hx]\n  exact reconX_lt _ _ _ (fun p hp => ({name}_columns_bound p hp).2)')
   sm+=theorem(f'{n}M{m}_ann_bound',f'∀ a ∈ sweepAnn m4T{n} {m}, a < 2 ^ 30','  decide +kernel')
 (out/'SweepMaps.lean').write_text(sm+'end '+NS+'\n')
 targets=[(n,i) for n in tables for i in range(len(splits[n]))] if args.all else [(args.instance,args.split)]
 imports=['MaskChunks','SweepMaps','PackedGridSweep']
 manifest=[]; instance_code={n:header(imports) for n in tables}; started=time.monotonic()
 for n,si in targets:
  p,q,r,m=splits[n][si];ow,iw,tw=weights(splits[n][si])
  ocap,icap=chunk_caps(ow,iw,args.outer_chunk,args.chunk)
  ocs=outer_chunks[ow,ocap];ics=chunks[iw,icap]
  xp=f'outerMasks{ow}S{ocap}';yp=f'masks{iw}S{icap}'
  left,right,ann=images(n,m);nm=f'{n}Split{si}'
  code=''
  cls_index={v:i for i,v in enumerate(tables[n]['cls'])}
  def candidate_index(x,y,a,l,v):
   output=l^v^a
   if m==0:u,w,t=x,output,y
   elif m==2:u,w,t=output,x,y
   else:u,w,t=x,y,output
   if m==2:
    c=xf(tables[n]['tWC'],x)^xf(tables[n]['tTC'],y)
    if any((tables[n][key]&c).bit_count()%2 for key in ['ln0','ln1']):return 0
   elif m>=3:
    c=xf(tables[n]['tUC'],x)^xf(tables[n]['tWC2'],y)
    if any((tables[n][key]&c).bit_count()%2 for key in ['lnT0','lnT1']):return 0
   if not (u|w|t):return 0
   return cls_index[u|(w<<30)|(t<<60)]
  for ai,a in enumerate(ann):
   for xi,xs in enumerate(ocs):
    xname=f'{xp}Chunk{xi}'
    lx=[xf(left,x) for x in xs]
    for yi,ys in enumerate(ics):
     tag=f'{nm}A{ai}X{xi}Y{yi}';yname=f'{yp}Chunk{yi}'
     ly=[xf(right,y) for y in ys]
     cert=[(i*len(ys)+j,candidate_index(x,ys[j],a,l,v))
       for i,(x,l) in enumerate(zip(xs,lx)) for j,v in enumerate(ly)
       if (l^v^a).bit_count()==tw]
     code+=decl(tag+'Certs','List (Nat × Nat)',natlist(cert),
       'Grid-lane indices and census witnesses for every matching weight.')
     words=f'gridWords {xname} {yname}\n      {n}M{m}LeftColumns {n}M{m}RightColumns'
     code+=theorem(tag+'_check',f'checkGrid m4T{n} {m} {a} {tw}\n    {xname} {yname}\n    ({words}) {tag}Certs = true','  decide +kernel')
     code+=theorem(tag+'_pairs',f'∀ x ∈ {xname}, ∀ y ∈ {yname},\n    popCntGo 30 (sweepOutput m4T{n} {m} x y {a}) = {tw} →\n      sweepAccept m4T{n} {m} x y {a} = true',f'''  apply checkGrid_sound m4T{n} {m} {a} {tw}
    {xname} {yname}
    ({words}) {tag}Certs
  · decide
  · decide
  · decide
  · intro x hx
    exact {n}M{m}Left_bound x ({xname}_bound x hx)
  · intro y hy
    exact {n}M{m}Right_bound y ({yname}_bound y hy)
  · apply gridWords_eq m4T{n} {m} {xname} {yname}
      {n}M{m}LeftColumns {n}M{m}RightColumns
    · decide
    · exact {xname}_bound
    · exact {yname}_bound
    · intro z hz
      exact ⟨({n}M{m}Left_columns_bound z hz).1,
        lt_of_lt_of_le ({n}M{m}Left_columns_bound z hz).2 (by decide)⟩
    · intro z hz
      exact ⟨({n}M{m}Right_columns_bound z hz).1,
        lt_of_lt_of_le ({n}M{m}Right_columns_bound z hz).2 (by decide)⟩
    · exact {n}M{m}Left_eq
    · exact {n}M{m}Right_eq
    · intro x hx
      exact {n}M{m}Left_bound x ({xname}_bound x hx)
    · intro y hy
      exact {n}M{m}Right_bound y ({yname}_bound y hy)
  · exact {tag}_check''')
    code+=theorem(f'{nm}A{ai}X{xi}_pairs',f'∀ x ∈ {xname}, ∀ y ∈ masksOfWt 30 {iw},\n    popCntGo 30 (sweepOutput m4T{n} {m} x y {a}) = {tw} →\n      sweepAccept m4T{n} {m} x y {a} = true',f'''  intro x hx y hy hw
  rw [← {yp}_cover] at hy
  obtain ⟨ys, hys, hy⟩ := List.mem_flatten.mp hy
  simp only [{yp}Chunks, List.mem_cons, List.not_mem_nil, or_false] at hys
  rcases hys with '''+cases(len(ics))+''.join(
      f'\n  · exact {nm}A{ai}X{xi}Y{yi}_pairs x hx y hy hw' for yi in range(len(ics))))
   code+=theorem(f'{nm}A{ai}_pairs',f'∀ x ∈ masksOfWt 30 {ow}, ∀ y ∈ masksOfWt 30 {iw},\n    popCntGo 30 (sweepOutput m4T{n} {m} x y {a}) = {tw} →\n      sweepAccept m4T{n} {m} x y {a} = true',f'''  intro x hx y hy hw
  rw [← {xp}_cover] at hx
  obtain ⟨xs, hxs, hx⟩ := List.mem_flatten.mp hx
  simp only [{xp}Chunks, List.mem_cons, List.not_mem_nil, or_false] at hxs
  rcases hxs with '''+cases(len(ocs))+''.join(
     f'\n  · exact {nm}A{ai}X{xi}_pairs x hx y hy hw' for xi in range(len(ocs))))
  code+=theorem(nm+'_check',f'checkSplit m4T{n} ({p}, {q}, {r}, {m}) = true',f'''  apply checkSplit_of_pairs
  · intro x hx y hy a ha
    apply sweepOutput_lt
    · exact {n}M{m}Left_bound x (masksOfWt_bound 30 {ow} x hx)
    · exact {n}M{m}Right_bound y (masksOfWt_bound 30 {iw} y hy)
    · exact {n}M{m}_ann_bound a ha
  · intro x hx y hy a ha hw
    have hann : a ∈ ({natlist(ann)} : List Nat) := ha
    simp only [List.mem_cons, List.not_mem_nil, or_false] at hann
    rcases hann with '''+' | '.join('rfl' for _ in ann)+''.join(f'\n    · exact {nm}A{ai}_pairs x hx y hy hw' for ai in range(len(ann))))
  if args.all:
   instance_code[n]+=code
  else:
   (out/(nm+'.lean')).write_text(header(imports)+code+'end '+NS+'\n')
  manifest.append(nm)
  print(nm,(p,q,r,m),'outer',len(masks[ow]),'inner',len(masks[iw]),'ann',len(ann),flush=True)
 if args.all:
  for n,code in instance_code.items():
   names=[f'{n}Split{i}_check' for i in range(len(splits[n]))]
   proof='  simp only [checkAll, splits'+n+', List.all_cons, List.all_nil,\n    '+',\n    '.join(', '.join(names[i:i+4]) for i in range(0,len(names),4))+', Bool.and_true]'
   code+=theorem('checkAll_'+n,f'checkAll m4T{n} splits{n} = true',proof)
   code+='#print axioms checkAll_'+n+'\n\n'
   (out/('Classifier'+n+'.lean')).write_text(code+'end '+NS+'\n')
  for side in ['X','Z']:
   (out/('FloorSweep'+side+'.lean')).write_text(BANNER+'\n'+''.join('import '+MOD+name+'\n' for name in ['Classifier'+side+'0','Classifier'+side+'1','Join'+side]))
 modules=['Classifier'+n for n in tables] if args.all else manifest
 generated=['MaskChunks','SweepMaps']+modules+(['FloorSweepX','FloorSweepZ'] if args.all else [])
 metadata={'modules':modules,'splits':manifest,'chunk_size':args.chunk,
   'outer_chunk_size':args.outer_chunk,
   'source_sha256':hashlib.sha256((out/'FloorData.lean').read_bytes()).hexdigest(),
   'generated_sha256':{name+'.lean':hashlib.sha256((out/(name+'.lean')).read_bytes()).hexdigest()
     for name in generated}}
 (out/'sweep-manifest.json').write_text(json.dumps(metadata,indent=2)+'\n')
 print('Generated',len(manifest),'splits in',round(time.monotonic()-started,2),'seconds',flush=True)
if __name__=='__main__': main()
