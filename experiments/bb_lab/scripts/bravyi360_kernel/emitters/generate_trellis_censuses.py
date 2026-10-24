import argparse,json
from pathlib import Path
parser=argparse.ArgumentParser(description='Generate certified packed trellis programs and five censuses.')
parser.add_argument('--trellis-data',type=Path,required=True)
parser.add_argument('--extra-cosets',type=Path,required=True)
parser.add_argument('--out',type=Path,required=True)
args=parser.parse_args();base=args.out;base.mkdir(parents=True,exist_ok=True)
D=json.loads(args.trellis_data.read_text());E=json.loads(args.extra_cosets.read_text())
ns='Quantum.Stabilizer.Homological.F2Mask.Trellis.Bravyi90'
src='import FastTrellis\n\nset_option maxRecDepth 1000000\nset_option maxHeartbeats 0\n\nnamespace '+ns+'\n\n'
for tag,key,slot in [('stab','spanRows',42),('cycle','cycleSpanRows',50)]:
 rows=[int(x,16) for x in D[key]]
 starts={(r&-r).bit_length()-1:i for i,r in enumerate(rows)};ends={r.bit_length()-1:i for i,r in enumerate(rows)}
 active=[];steps=[];readOrder=[]
 for j in range(90):
  initial=len(active);read=j in starts
  if read:active.append(starts[j]);readOrder.append(starts[j])
  w=len(active);ps=[i for i,r in enumerate(active) if rows[r]>>j&1];drop=active.index(ends[j]) if j in ends else w
  if drop<w:active[drop],active[-1]=active[-1],active[drop];active.pop()
  steps.append((j,initial,w,len(active),read,ps,drop))
 src+=f'/-- Terminal stage of the certified encoder. -/\nprivate def {tag}Raw90 (_offset degree : Nat) : PackedProgram {slot} degree 0 0 := PackedProgram.done 0\n\n'
 src+=f'/-- Terminal metadata suffix. -/\nprivate def {tag}Steps90 (_offset : Nat) : List PackedStep := []\n\n'
 src+=f'/-- The terminal stage leaves its table unchanged. -/\nprivate theorem {tag}_numeric90 (offset degree : Nat) : ({tag}Raw90 offset degree).numeric = (fun table => ({tag}Steps90 offset).foldl (packedStep {slot} degree) table) := by rfl\n\n'
 for j,initial,w,final,read,ps,drop in reversed(steps):
  tail=f'{tag}Raw{j+1} offset degree'
  if drop<w:tail=f'PackedProgram.forgetBitFast {drop} ({tail})'
  off=f'(offset.testBit {j} ^^ {str(bool(len(ps)%2)).lower()})'
  tail=f'PackedProgram.emit {ps} {off} (by decide) ({tail})'
  if read:tail=f'PackedProgram.read (width := {initial}) ({tail})'
  step=f'({w}, {str(read).lower()}, {ps}, {drop}, {off})'
  src+=f'/-- Certified suffix beginning at physical coordinate {j}. -/\nprivate def {tag}Raw{j} (offset degree : Nat) : PackedProgram {slot} degree {initial} 0 := {tail}\n\n'
  src+=f'/-- Coordinate metadata followed by the remaining numerical updates. -/\nprivate def {tag}Steps{j} (offset : Nat) : List PackedStep := {step} :: {tag}Steps{j+1} offset\n\n'
  src+=f'''/-- One coordinate preserves agreement between the certified and flattened evaluators. -/
private theorem {tag}_numeric{j} (offset degree : Nat) :
    ({tag}Raw{j} offset degree).numeric =
      (fun table => ({tag}Steps{j} offset).foldl (packedStep {slot} degree) table) := by
  funext table
  change ({tag}Raw{j+1} offset degree).numeric (packedStep {slot} degree table {step}) = _
  rw [{tag}_numeric{j+1}]
  exact List.foldl_cons.symm

'''
 src+=f'/-- Coordinate metadata for the entire fast recurrence. -/\ndef {tag}Steps (offset : Nat) : List PackedStep := {tag}Steps0 offset\n\n'
 instructions=[]
 for j,initial,w,final,read,ps,drop in steps:
  if read:instructions.append(f'.read {initial}')
  instructions.append(f'.emit {ps} (offset.testBit {j} ^^ {str(bool(len(ps)%2)).lower()})')
  if drop<w:instructions.append(f'.forget {drop} {w}')
 src+=f'/-- Flattened affine instruction list for the certified encoder. -/\ndef {tag}Instructions (offset : Nat) : List Instruction := [\n  '+',\n  '.join(instructions)+']\n\n'
 src+=f'/-- Structural equality of the certified and efficient numeric functions. -/\nprivate theorem {tag}_raw_numeric (offset degree : Nat) : ({tag}Raw0 offset degree).numeric = (fun table => ({tag}Steps offset).foldl (fastPackedStep {slot} degree) table) := by\n  funext table\n  rw [fastPackedFold_eq]\n  simpa only [{tag}Steps] using congrFun ({tag}_numeric0 offset degree) table\n\n'
 src+=f'/-- Structural equality of the certified and flattened instruction lists. -/\nprivate theorem {tag}_raw_program (offset degree : Nat) : ({tag}Raw0 offset degree).certificate.program = {tag}Instructions offset := by rfl\n\n'
 src+=f'/-- The certified encoder sums each coefficient once. -/\nprivate theorem {tag}_raw_merges (offset degree : Nat) : ({tag}Raw0 offset degree).merges = {len(rows)} := by rfl\n\n'
 src+=f'/-- Certified encoder with a flattened implementation for efficient kernel evaluation. -/\ndef {tag}Program (offset degree : Nat) : PackedProgram {slot} degree 0 0 :=\n  PackedProgram.flatten ({tag}Raw0 offset degree) ({tag}Instructions offset)\n    (fun table => ({tag}Steps offset).foldl (fastPackedStep {slot} degree) table) {len(rows)}\n    ({tag}_raw_program offset degree) ({tag}_raw_numeric offset degree) ({tag}_raw_merges offset degree)\n\n'
 src+=f'/-- Certified generators paired with their coefficient inputs. -/\ndef {tag}Inputs : List (Nat × Nat) := [\n'+',\n'.join(f'  ({row}, {1<<readOrder.index(r)})' for r,row in enumerate(rows))+']\n\n'
 src+=f'/-- Generators encoded by the trellis. -/\ndef {tag}Generators : List Nat := {tag}Inputs.map Prod.fst\n\n'
 for name,field,n in [('reads','reads',len(rows)),('outputs','outputs',90)]:
  src+=f'/-- The encoder has the required number of {name}. -/\ntheorem {tag}_{name} (offset : Nat) : {field} ({tag}Instructions offset) = {n} := by rfl\n\n'
 src+=f'/-- The encoder sums each coefficient once. -/\nprivate theorem {tag}_merges (offset degree : Nat) : ({tag}Program offset degree).merges = {len(rows)} := by rfl\n\n'
src+='end '+ns+'\n'
base.joinpath('ConcreteTrellis.lean').write_text(src)
specs=[('stab',0,22,D['stabs']['hist']),('W121',int(D['coset'],16),22,D['oneW']['hist']),('W151',int(E[0]['offset'],16),22,E[0]['hist']),('W238',int(E[1]['offset'],16),22,E[1]['hist']),('cycle',0,10,D['cycles']['hist'])]
ns='Quantum.Stabilizer.Homological.F2Mask.Trellis.Bravyi90'
for name,offset,deg,hist in specs:
 tag='cycle' if name.startswith('cycle') else 'stab';slot=50 if tag=='cycle' else 42
 expected=sum(int(c)<<(slot*int(w)) for w,c in hist.items());count=sum(int(c) for c in hist.values())
 out=f'import ConcreteTrellis\n\nset_option maxRecDepth 1000000\nset_option maxHeartbeats 0\n\nnamespace {ns}\n\n'
 out+=f'/-- The affine offset of the {name} census. -/\ndef {name}Offset : Nat := {offset}\n\n'
 out+=f'/-- The packed exact truncated weight enumerator of the {name} census. -/\ndef {name}Answer : Nat := {expected}\n\n'
 out+=f'/-- Kernel evaluation of the entire {name} packed recurrence. -/\ntheorem {name}_computed : ({tag}Steps {name}Offset).foldl (fastPackedStep {slot} {deg}) 1 = {name}Answer := by decide +kernel\n#print axioms {name}_computed\n\n'
 out+=f'/-- The {name} encoder starts at its specified affine offset. -/\ntheorem {name}_zero : runWord ({tag}Instructions {name}Offset) 0 0 = {name}Offset := by decide +kernel\n#print axioms {name}_zero\n\n'
 out+=f'/-- Every stored generator is realized by its certified input coefficient. -/\ntheorem {name}_input_checks : {tag}Inputs.all (fun gi => decide (gi.2 < 2 ^ reads ({tag}Instructions {name}Offset)) && decide (runWord ({tag}Instructions {name}Offset) gi.2 0 = gi.1 ^^^ {name}Offset)) = true := by decide +kernel\n\n'
 out+=f'/-- Reducing the packed weight enumerator counts all retained words. -/\ntheorem {name}_answer_count : {name}Answer % (2 ^ {slot} - 1) = {count} := by decide +kernel\n\n'
 out+=f'#print axioms {name}_computed\n#print axioms {name}_input_checks\nend {ns}\n'
 base.joinpath(f'Census{name}.lean').write_text(out)
print('Generated generic programs and five numeric census modules.')
bounds='import TrellisSpanCensus\n' + ''.join(f'import Census{name}\n' for name,_,_,_ in specs if name != 'cycle7')
bounds+='\nset_option maxRecDepth 1000000\nset_option maxHeartbeats 0\n\nnamespace '+ns+'\n\n'
for name,offset,deg,hist in specs:
 if name == 'cycle7':continue
 tag='cycle' if name.startswith('cycle') else 'stab';slot=50 if tag=='cycle' else 42
 count=sum(int(c) for c in hist.values())
 span='word' if offset==0 else f'(word ^^^ {name}Offset)'
 bounds+=f'''/-- Every light finite subset of this certified affine span obeys the exact census bound. -/
theorem {name}_card_le (target : Finset Nat)
    (hspan : ∀ word ∈ target, XorSpan {tag}Generators {span})
    (hlight : ∀ word ∈ target, bitCount 90 word ≤ {deg}) :
    target.card ≤ {count} := by
  have h := PackedProgram.card_le_of_span_compute ({tag}Program {name}Offset {deg})
    {tag}Inputs {name}Offset target {name}Answer (by decide)
    (by decide) (by decide)
    {name}_zero {name}_input_checks
    (fun word hw => by simpa only [{tag}Generators, {name}Offset, Nat.xor_zero] using hspan word hw)
    (fun word hw => by change bitCount 90 word ≤ {deg}; exact hlight word hw)
    {name}_computed (by decide +kernel)
  simpa only [{name}_answer_count] using h

#print axioms {name}_card_le

'''
bounds+='end '+ns+'\n'
base.joinpath('ConcreteCensusBounds.lean').write_text(bounds)

# Certificate data is evaluated by the kernel; no runtime duplicate is needed.
from kernel_only_data import mark_kernel_data
mark_kernel_data(args.out / name for name in ['ConcreteTrellis.lean', 'Censusstab.lean', 'CensusW121.lean', 'CensusW151.lean', 'CensusW238.lean', 'Censuscycle.lean', 'ConcreteCensusBounds.lean'])
