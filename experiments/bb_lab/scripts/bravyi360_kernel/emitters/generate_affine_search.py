"""Emit complete parity-search and exceptional deep-fiber certificates."""
import argparse
from pathlib import Path
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--top-certificates',type=Path,required=True)
parser.add_argument('--deep-certificates',type=Path,required=True)
parser.add_argument('--diagonal-certificates',type=Path,required=True)
parser.add_argument('--top-frame',type=Path,required=True)
parser.add_argument('--bottom-frame',type=Path,required=True)
parser.add_argument('--out',type=Path,required=True)
args=parser.parse_args();args.out.mkdir(parents=True,exist_ok=True)


def emit_search(args):
    import json,time,collections,sys,re
    Ftop=json.load(open(args.top_frame));Fbottom=json.load(open(args.bottom_frame))
    rows=[]
    for line in open(args.top_certificates):
     r=json.loads(line)
     if r['type']!='flat' and (r.get('hits',0)>0 or r['cap']>=3):rows.append(('top',r,Ftop))
    for r in json.load(open(args.deep_certificates))['fibers']:
     if r['type']!='flat':r['b']=r['beta'];rows.append(('deep',r,Fbottom))
    for r in json.load(open(args.diagonal_certificates))['logical']:
     if r['cap']==3:r['b']=hex(sum(1<<j for j in r['support']));rows.append(('diag',r,Ftop))

    def apply(cols,x):
     out=0
     while x:
      bit=x&-x;x-=bit;out^=cols[bit.bit_length()-1]
     return out

    import pathlib
    from balanced_nat_literal import literal
    base=args.out
    summary=[];t=time.time();totalnodes=0
    out='import SearchTree\n\nset_option maxRecDepth 1000000\nset_option maxHeartbeats 0\nnamespace Quantum.Stabilizer.Homological.F2Mask.AffineSearchData\n\n'
    for ix,(tag,r,F)in enumerate(rows):
     E=F['Ecols'];R=F['Rcols'];n=len(E);beta=int(r['b'],16);supp=r['support'];prows=[int(x,16)for x in r['leftMapRows']]
     def px(y):return sum(((p&y).bit_count()%2)<<j for p,j in zip(prows,supp))
     def proj(y):return y^apply(E,px(y))
     cols=[proj(x)for x in E];target=proj(apply(R,beta));rowMasks=[sum(((c>>q)&1)<<j for j,c in enumerate(cols))for q in range(n)]
     cap=r['cap'];expected={sum(1<<j for j in js)for js in r['off_subsets']};nodes=[0];solutions=set();usedRows=set()
     def search(syn,avail,left,selected):
      nodes[0]+=1
      if syn:
       rr=[];x=syn
       while x:
        bit=x&-x;x-=bit;q=bit.bit_length()-1;adj=rowMasks[q]&avail
        if not adj:usedRows.add(q);return f'(.zero {q})'
        rr.append((adj.bit_count(),adj,q))
       rr.sort();used=0;pack=[]
       for _,adj,q in rr:
        if not used&adj:
         pack.append(q);used|=adj
         if len(pack)>left:
          usedRows.update(pack);return '(.pack ['+','.join(map(str,pack))+'])'
       choices=rr[0][1];last=f'(.zero {rr[0][2]})';usedRows.add(rr[0][2])
      else:
       solutions.add(selected)
       if left==0:return '.done'
       choices=avail;last='.done'
      branches=[]
      while choices:
       bit=choices&-choices;choices-=bit;avail^=bit;j=bit.bit_length()-1
       branches.append((j,search(syn^cols[j],avail,left-1,selected|bit)))
      for j,yes in reversed(branches):last=f'(.split {j} {yes} {last})'
      return last
     tree=search(target,((1<<n)-1)^beta,cap,0)
     assert solutions==expected
     totalnodes+=nodes[0]
     selectedRows=sorted(usedRows); remap={q:i for i,q in enumerate(selectedRows)}
     tree=re.sub(r"\(\.zero (\d+)\)",lambda m:f"(.zero {remap[int(m[1])]})",tree)
     tree=re.sub(r"\(\.pack \[([0-9,]*)\]\)",lambda m:"(.pack ["+",".join(str(remap[int(v)])for v in m[1].split(","))+"])",tree)
     target=sum(((target>>q)&1)<<i for i,q in enumerate(selectedRows))
     packed=sum(rowMasks[q]<<(256*i)for i,q in enumerate(selectedRows))
     out+=f'/-- Packed parity neighborhoods for affine fiber {ix}. -/\ndef rows{ix} : Nat := {literal(packed,4096)}\n\n'
     out+=f'/-- Candidate outside supports for affine fiber {ix}. -/\ndef allowed{ix} : List Nat := ['+','.join(hex(v)for v in sorted(expected))+']\n\n'
     out+=f'/-- Exhaustive parity-search certificate for affine fiber {ix}. -/\ndef tree{ix} : ParitySearchTree := {tree}\n\n'
     out+=f'/-- The kernel verifies every branch of affine fiber {ix}. -/\ntheorem checked{ix} : checkParitySearch 8 (fun q => (rows{ix} >>> (256*q)) &&& (2^256-1)) {hex(target)} (allowed{ix}.contains) tree{ix} {hex(((1<<n)-1)^beta)} 0 {cap} = true := by decide +kernel\n\n'
     summary.append(dict(index=ix,tag=tag,beta=r['b'],cap=cap,target=hex(target),usedRows=sorted(usedRows),nodes=nodes[0],solutions=sorted(expected)))
     if ix==0:base.joinpath('AffineSearchFirst.lean').write_text(out+'#print axioms checked0\nend Quantum.Stabilizer.Homological.F2Mask.AffineSearchData\n')
     if ix%100==99:print('generated',ix+1,flush=True)
    out+='#print axioms checked0\n#print axioms checked656\nend Quantum.Stabilizer.Homological.F2Mask.AffineSearchData\n'
    base.joinpath('AffineSearchAll.lean').write_text(out)
    json.dump(summary,open(base/'affine_search_tree_metadata.json','w'))
    print('generated',len(rows),'cases',totalnodes,'nodes',len(out),'bytes',time.time()-t,'sec',flush=True)


def emit_initial(args):
    import json,sys
    from pathlib import Path
    from balanced_nat_literal import literal
    base=args.out;metadata=json.load(open(base/'affine_search_tree_metadata.json'))
    Ftop=json.load(open(args.top_frame));Fbottom=json.load(open(args.bottom_frame))
    records=[]
    for line in open(args.top_certificates):
     r=json.loads(line)
     if r['type']!='flat' and (r.get('hits',0)>0 or r['cap']>=3):records.append(('top',r,Ftop))
    for r in json.load(open(args.deep_certificates))['fibers']:
     if r['type']!='flat':r['b']=r['beta'];records.append(('bottom',r,Fbottom))
    for r in json.load(open(args.diagonal_certificates))['logical']:
     if r['cap']==3:r['b']=hex(sum(1<<j for j in r['support']));records.append(('top',r,Ftop))
    def app(cols,x):
     y=0
     while x:
      b=x&-x;x-=b;y^=cols[b.bit_length()-1]
     return y
    out='import SearchCertificate\nimport AffineSearchAll\nimport Bravyi360Adjoints\n\nset_option maxRecDepth 1000000\nset_option maxHeartbeats 0\nnamespace Quantum.Stabilizer.Homological.F2Mask.AffineSearchData\nopen PackedDigits BB.Bravyi360 BB.Bravyi360.FiberData\n\n'
    for label in ['top','bottom']:
     out+=f'/-- The {label} transpose fits in all256-bit search lanes. -/\ntheorem {label}Transpose_bounded : ∀ p ∈ {label}ETransposeColumns, p.1 < 256 ∧ p.2 < 2^256 := by decide +kernel\n\n'
    for ix,((label,r,F),m)in enumerate(zip(records,metadata)):
     E=F['Ecols'];R=F['Rcols'];n=len(E);beta=int(r['b'],16);supp=r['support'];prows=[int(x,16)for x in r['leftMapRows']]
     T=[sum(((v>>q)&1)<<j for j,v in enumerate(E))for q in range(n)]
     ls=[];outrows=[];target=0;rhs=app(R,beta)
     for i,q in enumerate(m['usedRows']):
      lam=1<<q
      for j,p in zip(supp,prows):
       if E[j]>>q&1:lam^=p
      row=app(T,lam);assert row&beta==0
      target|=((lam&rhs).bit_count()%2)<<i;ls.append(lam);outrows.append(row)
     assert target==int(m['target'],16)
     lp=sum(v<<(256*i)for i,v in enumerate(ls));count=len(ls)
     out+=f'/-- Dual row witnesses for affine search fiber {ix}. -/\ndef lambdas{ix} : Nat := {literal(lp,4096)}\n\n'
     out+=f'/-- The exact right-hand side of affine search fiber {ix}. -/\ndef rhs{ix} : Nat := {hex(rhs)}\n\n'
     out+=f'/-- The canonical carry map gives the stated right-hand side. -/\ntheorem rhsChecked{ix} : reconX {label}RColumns {hex(beta)} = rhs{ix} := by decide +kernel\n\n'
     out+=f'/-- The kernel verifies every initial dual row of affine search fiber {ix}. -/\ntheorem initialChecked{ix} : checkSearchInitial 8 {n} {count} rows{ix} lambdas{ix} {hex(target)} {hex(beta)} rhs{ix} {label}ETransposeColumns = true := by decide +kernel\n\n'
     out+=f'/-- Every light solution of affine search fiber {ix} has a listed outside support. -/\ntheorem complete{ix} {{word : Nat}} (hw : word < 2^{n})\n    (he : reconX {label}EColumns word = reconX {label}RColumns {hex(beta)})\n    (hweight : bitCount {n} (outsideMask word {hex(beta)}) ≤ {r["cap"]}) :\n    outsideMask word {hex(beta)} ∈ allowed{ix} := by\n  exact outside_mem_of_search 8 {n} {count} rows{ix} lambdas{ix} {hex(target)} {hex(beta)} rhs{ix} {r["cap"]}\n    {label}EColumns {label}ETransposeColumns tree{ix} allowed{ix} (by decide)\n    (fun lam v _ _ => {label}ETranspose_parityFold lam v) {label}Transpose_bounded initialChecked{ix} checked{ix}\n    hw (he.trans rhsChecked{ix}) hweight\n\n'
     if ix==0:base.joinpath('AffineInitialFirst.lean').write_text(out.replace('import AffineSearchAll','import AffineSearchFirst')+'#print axioms complete0\nend Quantum.Stabilizer.Homological.F2Mask.AffineSearchData\n')
    out+='#print axioms complete0\n#print axioms complete656\nend Quantum.Stabilizer.Homological.F2Mask.AffineSearchData\n'
    base.joinpath('AffineInitialAll.lean').write_text(out);print(len(records),len(out),'bytes')


def emit_deep(args):
    import json
    from pathlib import Path
    base=args.out;rows=[r for r in json.load(open(args.deep_certificates))['fibers']if r['type']!='flat']
    F=json.load(open(args.bottom_frame));E=F['Ecols'];R=F['Rcols']
    def app(cols,x):
     y=0
     while x:
      b=x&-x;x-=b;y^=cols[b.bit_length()-1]
     return y
    def pairs(xs):return '['+','.join(f'({i},{hex(x)})'for i,x in enumerate(xs)if x)+']'
    out='import AffineMaskCheck\nimport Bravyi360FiberData\n\nset_option maxRecDepth 1000000\nset_option maxHeartbeats 0\nnamespace Quantum.Stabilizer.Homological.BB.Bravyi360.DeepCensusData\nopen F2Mask\n\n'
    for ix,r in enumerate(rows):
     beta=int(r['beta'],16);supp=r['support'];prows=[int(x,16)for x in r['leftMapRows']]
     P=[sum(((row>>q)&1)<<j for row,j in zip(prows,supp))for q in range(90)]
     K=[sum(((int(d,16)>>i)&1)<<j for i,j in enumerate(supp))for d in r['deps']]
     D=[0]*90
     for i,f in enumerate(r['freeIndices']):D[supp[f]]=1<<i
     census=sorted(int(x['v0'],16)for x in r['lifts'])
     seeds=[]
     for js in r['off_subsets']:
      off=sum(1<<j for j in js)
      for k in range(1<<len(K)):
       v=off^app(P,app(R,beta)^app(E,off))^app(K,k);assert v in census or v^beta in census
      seeds.append(off)
     seeds.sort()
     out+=f'/-- The deck-normalized middle lifts of deep base fiber {ix}. -/\ndef lifts{ix} : List Nat := ['+','.join(map(hex,census))+']\n\n'
     out+=f'/-- The complete affine reconstruction certificate for deep base fiber {ix}. -/\ndef cert{ix} : AffineMaskCertificate := {{\n beta := {hex(beta)}\n sectionColumns := {pairs(P)}\n kernelColumns := {pairs(K)}\n decoderColumns := {pairs(D)}\n dimension := {len(K)}\n seeds := ['+','.join('('+hex(off)+',[])'for off in seeds)+'] }\n\n'
     out+=f'/-- Every affine candidate in deep base fiber {ix} is a listed lift up to deck action. -/\ntheorem checked{ix} : (cert{ix}).checkCensus 90 FiberData.bottomEColumns FiberData.bottomRColumns lifts{ix} = true := by decide +kernel\n\n'
    out+='#print axioms checked0\n#print axioms checked1\nend Quantum.Stabilizer.Homological.BB.Bravyi360.DeepCensusData\n'
    base.joinpath('DeepCensusData.lean').write_text(out);print(len(out),'bytes')


def emit_split(args):
    from pathlib import Path
    import re
    base=args.out;ex={4,14,18,19,653,654,655,656}
    for stem,marker,footer,common in [
     ('AffineSearch','Packed parity neighborhoods for affine fiber','end Quantum.Stabilizer.Homological.F2Mask.AffineSearchData',None),
     ('AffineInitial','Dual row witnesses for affine search fiber','end Quantum.Stabilizer.Homological.F2Mask.AffineSearchData','AffineInitialCommon')]:
     p=base/(stem+'All.lean');text=p.read_text();pattern=r'/-- '+re.escape(marker)+r' (\d+)\. -/'
     starts=list(re.finditer(pattern,text));header=text[:starts[0].start()];blocks={}
     for j,m in enumerate(starts):blocks[int(m[1])]=text[m.start():starts[j+1].start() if j+1<len(starts) else text.index('#print axioms',m.start())]
     chunks=[('Exceptions',[i for i in blocks if i in ex])]
     rest=[i for i in blocks if i not in ex]
     chunks += [(f'Part{j//100}',rest[j:j+100]) for j in range(0,len(rest),100)]
     names=[]
     if common:
      ch=header.replace('import AffineSearchAll\n','')+footer+'\n';base.joinpath(common+'.lean').write_text(ch)
      header='import AffineInitialCommon\nIMPORTSEARCH\n\nset_option maxRecDepth 1000000\nset_option maxHeartbeats 0\nnamespace Quantum.Stabilizer.Homological.F2Mask.AffineSearchData\nopen PackedDigits BB.Bravyi360 BB.Bravyi360.FiberData\n\n'
     for suffix,ids in chunks:
      if not ids:continue
      name=stem+suffix;names.append(name)
      hd=header.replace('IMPORTSEARCH','import AffineSearch'+suffix)
      content=hd+''.join(blocks[i] for i in ids)+footer+'\n'
      base.joinpath(name+'.lean').write_text(content)
     base.joinpath(stem+'All.lean').write_text(''.join('import '+name+'\n'for name in names))
     print(stem,[(suffix,len(ids))for suffix,ids in chunks])


def emit_exceptions(args):
    from pathlib import Path
    import json
    base=args.out;m=json.load(open(base/'affine_search_tree_metadata.json'))
    out='import AffineInitialExceptions\nimport DeepCensusData\n\nnamespace Quantum.Stabilizer.Homological.BB.Bravyi360.AffineExceptions\nopen F2Mask F2Mask.AffineSearchData FiberData\n\n'
    for ix in [4,14,18,19,655,656]:
     beta=m[ix]['beta']
     out+=f'/-- No solution exists in this exceptional top fiber with at most three outside bits. -/\ntheorem empty{ix} {{word : Nat}} (hw : word < 2^180)\n    (he : reconX topEColumns word = reconX topRColumns {beta})\n    (hweight : bitCount 180 (outsideMask word {beta}) ≤ 3) : False := by\n  have h := complete{ix} hw he hweight\n  simpa only [allowed{ix}, List.not_mem_nil] using h\n\n'
    for j,ix in enumerate([653,654]):
     beta=m[ix]['beta']
     out+=f'/-- The outside-support census in the deep affine certificate is the checked search census. -/\ntheorem deepOutside{j} : (DeepCensusData.cert{j}).seeds.map Prod.fst = allowed{ix} := by decide +kernel\n\n'
     out+=f'/-- Every sheet in the deep base fiber is a listed normalized lift or its deck translate. -/\ntheorem deepComplete{j} {{word : Nat}} (hw : word < 2^90)\n    (he : reconX bottomEColumns word = reconX bottomRColumns {beta})\n    (hweight : bitCount 90 (outsideMask word {beta}) ≤ 5) :\n    word ∈ DeepCensusData.lifts{j} ∨ (word ^^^ {beta}) ∈ DeepCensusData.lifts{j} := by\n  apply AffineMaskCertificate.census_of_check (DeepCensusData.cert{j}) 90 bottomEColumns bottomRColumns\n    (DeepCensusData.lifts{j}) (DeepCensusData.checked{j}) hw he\n  rw [deepOutside{j}]\n  exact complete{ix} hw he hweight\n\n'
    out+='#print axioms empty655\n#print axioms deepComplete0\nend Quantum.Stabilizer.Homological.BB.Bravyi360.AffineExceptions\n'
    base.joinpath('AffineExceptions.lean').write_text(out)


def emit_bank(args):
    import re,json,pathlib
    P=args.out
    s=(P/'DeepCensusData.lean').read_text()
    bank=[int(json.loads(l)['b'],16) for l in open(args.top_certificates)]
    lookup={b:i for i,b in enumerate(bank)}
    def join(a,b):
     z=0
     for i in range(90):
      j=(i//45)*90+i%45
      z|=((a>>i)&1)<<j
      z|=(((a^b)>>i)&1)<<(j+45)
     return z
    out='''import DeepCensusData
    import AffineExceptions
    import Bravyi360Join
    import Bravyi360TopBank

    namespace Quantum.Stabilizer.Homological.BB.Bravyi360.DeepBank
    open F2Mask
    set_option maxHeartbeats 0
    set_option maxRecDepth 200000

    /-- Verify the bank index attached to each explicitly listed normalized sheet. -/
    def check (beta : Nat) (entries : List (Nat × Nat)) : Bool :=
      entries.all (fun p => decide (p.2 < 58130 ∧ bottomJoin p.1 beta = topBankMask p.2))

    /-- A checked lookup table puts every listed sheet into the common top bank. -/
    theorem mem_bank {beta a : Nat} {entries : List (Nat × Nat)}
        (hc : check beta entries = true) (ha : a ∈ entries.map Prod.fst) :
        ∃ j : Fin 58130, bottomJoin a beta = topBankMask j := by
      obtain ⟨p, hp, hpa⟩ := List.mem_map.mp ha
      have h := of_decide_eq_true ((List.all_eq_true.mp hc) p hp)
      exact ⟨⟨p.2, h.1⟩, hpa ▸ h.2⟩
    '''
    for k,b in enumerate([0x7008000c00001200200049,0x100020004018007000e001]):
     vals=[int(x,16) for x in re.search(r'def lifts'+str(k)+r' : List Nat := \[(.*?)\]',s).group(1).split(',')]
     pairs=[(a,lookup[join(a,b)]) for a in vals]
     out+=f'''\n/-- Bank indices of the normalized cap-five lifts in exceptional fiber {k}. -/
    def entries{k} : List (Nat × Nat) := [{','.join(f'({hex(a)},{j})' for a,j in pairs)}]

    /-- The indexed sheet table is exactly the certified affine census. -/
    theorem entries{k}_sheets : entries{k}.map Prod.fst = DeepCensusData.lifts{k} := by decide +kernel

    /-- Each normalized exceptional lift belongs to the common top bank. -/
    theorem checked{k} : check {hex(b)} entries{k} = true := by decide +kernel

    /-- Every light solution of exceptional deep fiber {k} reaches the top bank up to deck action. -/
    theorem complete{k} {{a : Nat}} (ha : a < 2^90)
        (he : reconX FiberData.bottomEColumns a = reconX FiberData.bottomRColumns {hex(b)})
        (hw : bitCount 90 (outsideMask a {hex(b)}) ≤ 5) :
        ∃ j : Fin 58130, bottomJoin a {hex(b)} = topBankMask j ∨
          bottomJoin (a ^^^ {hex(b)}) {hex(b)} = topBankMask j := by
      rcases AffineExceptions.deepComplete{k} ha he hw with h | h
      · obtain ⟨j, hj⟩ := mem_bank checked{k} (entries{k}_sheets.symm ▸ h)
        exact ⟨j, Or.inl hj⟩
      · obtain ⟨j, hj⟩ := mem_bank checked{k} (entries{k}_sheets.symm ▸ h)
        exact ⟨j, Or.inr hj⟩
    '''
    out+='\n#print axioms complete0\n#print axioms complete1\nend Quantum.Stabilizer.Homological.BB.Bravyi360.DeepBank\n'
    (P/'DeepBank.lean').write_text(''.join(line[4:] if line.startswith('    ') else line for line in out.splitlines(keepends=True)))
    print('generated DeepBank',len(out))


def emit_deep_order(args):
    import json
    fibers=json.loads(args.deep_certificates.read_text())['fibers']
    ordered=[r for r in fibers if r['type']=='flat']+[r for r in fibers if r['type']!='flat']
    assert len(ordered)==28
    packed=sum(int(r['beta'],16)<<(90*i) for i,r in enumerate(ordered))
    source=f'''import PackedDigits

namespace Quantum.Stabilizer.Homological.BB.Bravyi360.DeepFiberData
/-- The exact certificate shadows: twenty-six flat fibers, followed by two cap-five fibers. -/
def betaPacked : Nat := {hex(packed)}

/-- Read one deep shadow in its original affine-certificate coordinate frame. -/
def betaAt (i : Nat) : Nat := PackedDigits.lane 90 i betaPacked

end Quantum.Stabilizer.Homological.BB.Bravyi360.DeepFiberData
'''
    args.out.joinpath('DeepFiberData.lean').write_text(source)

emit_search(args)
emit_initial(args)
emit_deep(args)
emit_split(args)
emit_exceptions(args)
emit_bank(args)
emit_deep_order(args)

# Certificate data is evaluated by the kernel; no runtime duplicate is needed.
from kernel_only_data import mark_kernel_data
mark_kernel_data(args.out / name for name in ['AffineSearchExceptions.lean', 'AffineSearchPart0.lean', 'AffineSearchPart1.lean', 'AffineSearchPart2.lean', 'AffineSearchPart3.lean', 'AffineSearchPart4.lean', 'AffineSearchPart5.lean', 'AffineSearchPart6.lean', 'AffineInitialCommon.lean', 'AffineInitialExceptions.lean', 'AffineInitialPart0.lean', 'AffineInitialPart1.lean', 'AffineInitialPart2.lean', 'AffineInitialPart3.lean', 'AffineInitialPart4.lean', 'AffineInitialPart5.lean', 'AffineInitialPart6.lean', 'AffineExceptions.lean', 'DeepCensusData.lean', 'DeepBank.lean', 'DeepFiberData.lean'])
