import os
from pathlib import Path as _StagePath
_STAGE = _StagePath(os.environ["QEC360_STAGE"])
import pathlib,json,sys,collections
sys.path.insert(0,str(_STAGE / 'qec360-certprobe'))
from balance_literals import lit
from compact_literals import packlit
P=pathlib.Path(str(_STAGE / 'qec360-certprobe'))
def gen(top):
 if top:
  groups=collections.defaultdict(list)
  for l in open(P/'top_empty_dual90_certificates.jsonl'):
   r=json.loads(l)
   if r['cap']==1:groups[len(r['duals'])].append(r)
  rs=[r for g in groups.values()for r in g]
 else:rs=[dict(beta=a[0],duals=a[2:])for l in open(P/'emptycap45_1.tsv')if(a:=l.split())]
 prefix='Top'if top else'Bottom';short=prefix.lower();w=512 if top else 256;n=180 if top else 90;k=8 if top else 7;N=len(rs);chunks=(N+1023)//1024;last=N%1024
 f=open(P/f'Bravyi360{prefix}Empty1.lean','w');f.write(f'''import Bravyi360{prefix}Empty0
import BatchCap1Raw
set_option maxHeartbeats 0
set_option maxRecDepth 100000
namespace Quantum.Stabilizer.Homological.BB.Bravyi360.{prefix}Empty
open F2Mask FiberData PackedDigits
''')
 if top and last not in [1024,301]:f.write(f'/-- Repetition mask for the last cap-one batch. -/\ntheorem rep{last}_eq : rep{w} {last} = packF {w} {last} (fun _ => 1) := by decide +kernel\n')
 for ix,off in enumerate(range(0,N,1024)):
  sub=rs[off:off+1024];size=len(sub);num=max(len(r['duals'])for r in sub);ls=[]
  for j in range(num):
   vals=[int(r['duals'][j],16)if j<len(r['duals'])else 0 for r in sub]
   ls.append(lit(sum(v<<(w*i)for i,v in enumerate(vals)))if top else packlit(vals,w))
  b=f'beta1Words{ix}';cs=f'coeff1Words{ix}';rep=f'(rep{w} {size})';rows=f'{short}DualColumns'
  f.write(f'/-- Dual coefficient tables for cap-one batch {ix}. -/\ndef {cs} : List Nat := ['+','.join(ls)+']\n')
  f.write(f'''/-- Parallel completeness checks for cap-one batch {ix}. -/
theorem checked1_{ix} :
    {b} < 2^({w}*{size}) ∧ {b} &&& ({rep}*(2^{n}-1)) = {b} ∧
    {cs}.all (fun c => Nat.blt c (2^({w}*{size})) && Nat.beq (batchRecon {rep} c {rows} &&& {b}) 0) = true ∧
    batchAllowed1 {n} {k} {rep} {b} {rows} {cs} = 0 ∧
    batchZero1 {n} {k} {rep} {b} {rows} {cs} = 0 := by decide +kernel
/-- Every lane in cap-one batch {ix} excludes a low-overflow solution. -/
theorem impossible1_{ix} (i : Nat) (hi : i < {size}) (word : Nat) (hw : word < 2^{n})
    (he : reconX {short}EColumns word = reconX {short}RColumns (lane {w} i {b}))
    (hoff : bitCount {n} (outsideMask word (lane {w} i {b})) ≤ 1) : False := by
  have h := checked1_{ix}
  rw [rep{size}_eq] at h
  have hc := h.2.2.1
  simp only [List.all_eq_true,Bool.and_eq_true,Nat.blt_eq,Nat.beq_eq] at hc
  exact cap1_raw_impossible {w} {size} {n} {k} {b} {cs} (by decide) (by decide) (by decide) (by decide)
    h.1 h.2.1 {rows} {short}EColumns {short}RColumns {rows}_bounded
    (fun c hm => (hc c hm).1) (fun c hm => (hc c hm).2) h.2.2.2.1 h.2.2.2.2
    {short}Dual_constraint i hi word hw he hoff
''')
 f.write(f'/-- Size of one cap-one batch. -/\ndef size1 (b : Nat) : Nat := if b < {chunks-1} then 1024 else {last}\n')
 f.write(f'''/-- Cap-one exclusion for every shadow batch. -/
theorem cap1_batch_impossible (b : Nat) (hb : b < {chunks}) (i : Nat) (hi : i < size1 b)
    (word : Nat) (hw : word < 2^{n})
    (he : reconX {short}EColumns word = reconX {short}RColumns (lane {w} i (beta1Words b)))
    (hoff : bitCount {n} (outsideMask word (lane {w} i (beta1Words b))) ≤ 1) : False := by
''')
 def branch(lo,hi,ind):
  if hi-lo==1:f.write(ind+f'have heq : b = {lo} := by omega\n'+ind+'subst b\n'+ind+f'exact impossible1_{lo} i hi word hw he hoff\n')
  else:
   mid=(lo+hi)//2;f.write(ind+f'by_cases h{mid} : b < {mid}\n'+ind+'·\n');branch(lo,mid,ind+'  ');f.write(ind+'·\n');branch(mid,hi,ind+'  ')
 branch(0,chunks,'  ')
 f.write(f'''/-- Every listed cap-one fiber is impossible with at most one outside coordinate. -/
theorem cap1_impossible (i : Fin {N}) (word : Nat) (hw : word < 2^{n})
    (he : reconX {short}EColumns word = reconX {short}RColumns (beta1 i.val))
    (hoff : bitCount {n} (outsideMask word (beta1 i.val)) ≤ 1) : False := by
  have hm := Nat.mod_lt i.val (by decide : 0 < 1024)
  have heq := Nat.mod_add_div i.val 1024
  have hb : i.val / 1024 < {chunks} := by omega
  have hl : i.val % 1024 < size1 (i.val / 1024) := by
    unfold size1
    split_ifs <;> omega
  exact cap1_batch_impossible _ hb _ hl word hw he hoff
#print axioms cap1_impossible
end Quantum.Stabilizer.Homological.BB.Bravyi360.{prefix}Empty
''');f.close();print(prefix,N)
gen(True);gen(False)

# Data are kernel certificates, not runtime programs.
for _certificate_output in ['Bravyi360TopEmpty1.lean', 'Bravyi360BottomEmpty1.lean']:
 _certificate_path = P / _certificate_output
 _certificate_path.write_text(_certificate_path.read_text().replace("\ndef ", "\nnoncomputable def "))
