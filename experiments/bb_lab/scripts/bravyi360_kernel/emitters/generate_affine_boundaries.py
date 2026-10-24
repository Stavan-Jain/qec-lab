import argparse,json
from pathlib import Path
from balanced_nat_literal import literal
parser=argparse.ArgumentParser(description='Generate kernel-checked nonempty affine fiber certificates.')
parser.add_argument('--top-certificates',type=Path,required=True)
parser.add_argument('--top-frame',type=Path,required=True)
parser.add_argument('--out',type=Path,required=True)
args=parser.parse_args();base=args.out;base.mkdir(parents=True,exist_ok=True);rs=[]
for line in args.top_certificates.open():
 r=json.loads(line)
 if r['type']!='flat' and (r.get('hits',0)>0 or r['cap']>=3):rs.append(r)
rs=[(i,r)for i,r in enumerate(rs)if r['hits']]
F=json.loads(args.top_frame.read_text());E=F['Ecols'];R=F['Rcols']
def app(cols,x):
 y=0
 while x:
  b=x&-x;x-=b;y^=cols[b.bit_length()-1]
 return y
def pairs(xs):return '['+','.join(f'({i},{hex(x)})'for i,x in enumerate(xs)if x)+']'
for part in range(7):
 group=rs[100*part:100*part+100];coeffs=[];sheets=[];shadows=[];blocks=[]
 for ix,r in group:
  beta=int(r['b'],16);supp=r['support'];prows=[int(x,16)for x in r['leftMapRows']]
  P=[sum(((row>>q)&1)<<j for row,j in zip(prows,supp))for q in range(180)]
  K=[sum(((int(d,16)>>i)&1)<<j for i,j in enumerate(supp))for d in r['deps']]
  D=[0]*180
  for i,f in enumerate(r['freeIndices']):D[supp[f]]=1<<i
  vs=[int(v,16)for v in r['solutions']];idx={v:i for i,v in enumerate(vs)};count=len(vs);start=len(coeffs)
  coeffs.extend(int(c,16)for c in r['coverGeneratorCoeffs']);sheets.extend(vs);shadows.extend([beta]*count)
  seeds=[]
  for js in r['off_subsets']:
   off=sum(1<<j for j in js);is_=[]
   for k in range(1<<len(K)):
    v=off^app(P,app(R,beta)^app(E,off))^app(K,k);is_.append(idx[v])
   seeds.append((off,is_))
  seeds.sort();sheet=sum(v<<(512*i)for i,v in enumerate(vs));shadow=sum(beta<<(512*i)for i in range(count))
  s=f'''/-- The local candidate sheets of top fiber {ix}. -/
def sheets{ix} : Nat := {literal(sheet,4096)}
/-- The shadow accompanying each candidate. -/
def shadows{ix} : Nat := {literal(shadow,4096)}
/-- The local candidates form the specified slice of their checked boundary group. -/
theorem slices{ix} : sheets{ix} = boundedSlice 512 {start} {count} sheetTable ∧
    shadows{ix} = boundedSlice 512 {start} {count} shadowTable := by decide +kernel
/-- Every local candidate inherits a physical boundary witness. -/
theorem boundaryAt{ix} (i : Nat) (hi : i < {count}) : ∃ coefficient,
    reconX boundaryColumns coefficient = topJoin (lane 512 i sheets{ix}) (lane 512 i shadows{ix}) :=
  predicate_of_slice (fun a b => ∃ c, reconX boundaryColumns c = topJoin a b)
    512 count {start} {count} sheetTable shadowTable _ _ (by decide) (by decide +kernel)
    boundaryAt slices{ix}.1 slices{ix}.2 i hi

'''
  s+=f'/-- Affine reconstruction data for top fiber {ix}. -/\ndef cert{ix} : AffineMaskCertificate := {{\n beta := {hex(beta)}\n sectionColumns := {pairs(P)}\n kernelColumns := {pairs(K)}\n decoderColumns := {pairs(D)}\n dimension := {len(K)}\n seeds := ['+','.join('('+hex(off)+',['+','.join(map(str,inds))+'])'for off,inds in seeds)+'] }\n\n'
  rhs=hex(app(R,beta))
  s+=f'''/-- Every reconstructed candidate occurs in the local boundary table. -/
theorem checked{ix} : (cert{ix}).checkIndexedBoundary 180 512 {count} sheets{ix} shadows{ix} FiberData.topEColumns FiberData.topRColumns = true := by
  have hr : {rhs} = reconX FiberData.topRColumns (cert{ix}).beta := AffineSearchData.rhsChecked{ix}.symm
  have h : (cert{ix}).checkIndexedBoundaryCached 180 512 {count} sheets{ix} shadows{ix} {rhs} DiagonalData.eColumn FiberData.topEColumns FiberData.topRColumns = true := by decide +kernel
  rwa [AffineMaskCertificate.checkIndexedBoundaryCached_eq _ _ _ _ _ _ _ _ _ _ topEColumn_eq hr] at h

'''
  blocks.append((ix,s,count))
 n=len(coeffs);out=f'''import AffineBoundaryMaps

namespace Quantum.Stabilizer.Homological.BB.Bravyi360.AffineBoundaryGroups.Group{part}
open F2Mask PackedDigits
set_option maxHeartbeats 0
set_option maxRecDepth 200000
/-- Number of explicit candidates in this group. -/
def count : Nat := {n}
'''
 for name,vs in [('coefficientTable',coeffs),('sheetTable',sheets),('shadowTable',shadows)]:
  out+=f'/-- Packed {name} for this boundary group. -/\ndef {name} : Nat := {literal(sum(v<<(512*i)for i,v in enumerate(vs)),4096)}\n/-- The table occupies the stated lane count. -/\ntheorem {name}_bounded : {name} < 2^(512*count) := by decide +kernel\n'
 out+='''/-- One low bit in each candidate lane. -/
def repeatedOne : Nat := (2^(512*count)-1)/(2^512-1)
/-- The repeated-bit formula equals the standard packed constant. -/
theorem repeatedOne_eq : packF 512 count (fun _ => 1) = repeatedOne := by decide +kernel
/-- Verify every boundary reconstruction and sheet join in this group in parallel. -/
theorem checked : batchRecon (packF 512 count (fun _ => 1)) coefficientTable topBoundaryFast =
    batchRecon (packF 512 count (fun _ => 1)) sheetTable topLeftFast ^^^
      batchRecon (packF 512 count (fun _ => 1)) (sheetTable ^^^ shadowTable) topRightFast := by
  rw [repeatedOne_eq]
  decide +kernel
/-- Every candidate in this group is a physical boundary. -/
theorem boundaryAt (i : Nat) (hi : i < count) : ∃ coefficient, reconX boundaryColumns coefficient =
    topJoin (lane 512 i sheetTable) (lane 512 i shadowTable) := by
  refine ⟨lane 512 i coefficientTable, ?_⟩
  have h := boundary_join_of_batch 512 count coefficientTable sheetTable shadowTable
    topBoundaryFast topLeftFast topRightFast (by decide) coefficientTable_bounded
    sheetTable_bounded shadowTable_bounded affineBoundaryMaps_bounded.1
    affineBoundaryMaps_bounded.2.1 affineBoundaryMaps_bounded.2.2 checked hi
  simpa only [topBoundaryFast_eq, topLeftFast_eq, topRightFast_eq, topJoin] using h
'''+f'end Quantum.Stabilizer.Homological.BB.Bravyi360.AffineBoundaryGroups.Group{part}\n'
 base.joinpath(f'AffineBoundaryGroup{part}.lean').write_text(out)
 header=f'import AffineColumnCache\nimport AffineInitialPart{part}\nimport AffineBoundaryGroup{part}\nimport Bravyi360FiberData\n\nset_option maxHeartbeats 0\nset_option maxRecDepth 1000000\nnamespace Quantum.Stabilizer.Homological.BB.Bravyi360.AffineBoundaryCached\nopen F2Mask PackedDigits AffineBoundaryGroups.Group{part}\n\n'
 base.joinpath(f'AffineBoundaryCachedPart{part}.lean').write_text(header+''.join(s for _,s,_ in blocks)+'end Quantum.Stabilizer.Homological.BB.Bravyi360.AffineBoundaryCached\n')
 sem=f'import AffineInitialPart{part}\nimport AffineBoundaryCachedPart{part}\nimport TopFiberBoundary\n\nnamespace Quantum.Stabilizer.Homological.BB.Bravyi360.TopNonempty\nopen F2Mask FiberData\n\n'
 for ix,r in group:
  b=r['b'];cap=r['cap'];count=len(r['solutions'])
  sem+=f'''/-- The affine witness and search certificates use the same outside-support census. -/
theorem census{ix} : (AffineBoundaryCached.cert{ix}).seeds.map Prod.fst = AffineSearchData.allowed{ix} := by decide +kernel

/-- Every short top lift in nonempty affine fiber {ix} is a boundary. -/
theorem boundary{ix} : TopFiberBoundary {b} := by
  intro word hw he hweight
  have hcap : (22-bitCount 180 {b})/2 = {cap} := by decide +kernel
  rw [hcap] at hweight
  have hoff := AffineSearchData.complete{ix} hw he hweight
  have hseed : outsideMask word (AffineBoundaryCached.cert{ix}).beta ∈
      (AffineBoundaryCached.cert{ix}).seeds.map Prod.fst := by
    rw [census{ix}]
    exact hoff
  exact AffineMaskCertificate.indexedBoundary_of_check (AffineBoundaryCached.cert{ix}) 180 512
    {count} AffineBoundaryCached.sheets{ix} AffineBoundaryCached.shadows{ix}
    topEColumns topRColumns boundaryColumns topJoin AffineBoundaryCached.boundaryAt{ix}
    (AffineBoundaryCached.checked{ix}) hw he hseed

'''
 sem+='end Quantum.Stabilizer.Homological.BB.Bravyi360.TopNonempty\n'
 base.joinpath(f'TopNonemptyPart{part}Cached.lean').write_text(sem)
 print(part,n,len(out))
items=sorted((int(r['b'],16),ix)for ix,r in rs)
out=''.join(f'import TopNonemptyPart{p}Cached\n'for p in range(7))+'\nset_option maxRecDepth 200000\nset_option maxHeartbeats 0\n\nnamespace Quantum.Stabilizer.Homological.BB.Bravyi360.TopNonempty\n\n'
out+='/-- All 649 nonempty affine top fibers, each with its boundary theorem. -/\ndef certifiedCensus : List {beta : Nat // TopFiberBoundary beta} := [\n'+',\n'.join(f'  ⟨{hex(b)}, boundary{i}⟩'for b,i in items)+']\n\n'
out+='''/-- The sorted shadow census for the nonempty affine top fibers. -/
def shadows : List Nat := certifiedCensus.map Subtype.val

/-- Every shadow in the nonempty top census has only boundary lifts of weight at most 22. -/
theorem boundary_of_mem {beta : Nat} (h : beta ∈ shadows) : TopFiberBoundary beta := by
  obtain ⟨entry, _, he⟩ := List.mem_map.mp h
  simpa only [he] using entry.property

#print axioms boundary_of_mem
end Quantum.Stabilizer.Homological.BB.Bravyi360.TopNonempty
'''
base.joinpath('TopNonempty.lean').write_text(out)

# Certificate data is evaluated by the kernel; no runtime duplicate is needed.
from kernel_only_data import mark_kernel_data
mark_kernel_data(args.out / name for name in ['AffineBoundaryGroup0.lean', 'AffineBoundaryGroup1.lean', 'AffineBoundaryGroup2.lean', 'AffineBoundaryGroup3.lean', 'AffineBoundaryGroup4.lean', 'AffineBoundaryGroup5.lean', 'AffineBoundaryGroup6.lean', 'AffineBoundaryCachedPart0.lean', 'AffineBoundaryCachedPart1.lean', 'AffineBoundaryCachedPart2.lean', 'AffineBoundaryCachedPart3.lean', 'AffineBoundaryCachedPart4.lean', 'AffineBoundaryCachedPart5.lean', 'AffineBoundaryCachedPart6.lean', 'TopNonemptyPart0Cached.lean', 'TopNonemptyPart1Cached.lean', 'TopNonemptyPart2Cached.lean', 'TopNonemptyPart3Cached.lean', 'TopNonemptyPart4Cached.lean', 'TopNonemptyPart5Cached.lean', 'TopNonemptyPart6Cached.lean', 'TopNonempty.lean'])
