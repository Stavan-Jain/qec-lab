import argparse
from pathlib import Path
parser=argparse.ArgumentParser(description='Emit literal top boundary and sheet embedding columns.')
parser.add_argument('--out',type=Path,required=True)
args=parser.parse_args();args.out.mkdir(parents=True,exist_ok=True)
A=[(9,0),(0,1),(0,2)];B=[(0,3),(25,0),(26,0)]
G=[]
for i in range(180):
 x,y=divmod(i,6)
 G.append(sum(1<<(block*180+((x-a)%30)*6+(y-b)%6) for block,polynomial in enumerate([A,B]) for a,b in polynomial))
L=[1<<((i//90)*180+(i%90//3)*6+i%3)for i in range(180)];R=[x<<3 for x in L]
def pairs(xs):return '['+','.join(f'({i},{hex(x)})'for i,x in enumerate(xs))+']'
out='import Bravyi360Join\nimport Bravyi360Packaging\n\nset_option maxRecDepth 200000\nset_option maxHeartbeats 0\nnamespace Quantum.Stabilizer.Homological.BB.Bravyi360\nopen F2Mask\n\n'
for name,xs,old in [('topBoundaryFast',G,'boundaryColumns'),('topLeftFast',L,'topLeft'),('topRightFast',R,'topRight')]:
 out+=f'/-- Literal columns for the shared {old} map. -/\ndef {name} : List (Nat × Nat) := {pairs(xs)}\n\n/-- The literal columns equal the canonical map. -/\ntheorem {name}_eq : {name} = {old} := by decide +kernel\n\n'
out+='/-- Literal-column evaluator of the top join. -/\ndef topJoinFast (a beta : Nat) : Nat := reconX topLeftFast a ^^^ reconX topRightFast (a ^^^ beta)\n\n/-- The literal evaluator is the canonical top join. -/\ntheorem topJoinFast_eq : topJoinFast = topJoin := by\n  funext a beta\n  simp only [topJoinFast, topJoin, topLeftFast_eq, topRightFast_eq]\n\nend Quantum.Stabilizer.Homological.BB.Bravyi360\n'
args.out.joinpath('AffineFastMaps.lean').write_text(out)

# Certificate data is evaluated by the kernel; no runtime duplicate is needed.
from kernel_only_data import mark_kernel_data
mark_kernel_data(args.out / name for name in ['AffineFastMaps.lean'])
