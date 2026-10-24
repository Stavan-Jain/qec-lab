"""Emit kernel-checkable certificates into QEC360_STAGE."""
import os
from pathlib import Path
STAGE = Path(os.environ['QEC360_STAGE']).resolve()
import pathlib
p = pathlib.Path(str(STAGE))
ns = 'Quantum.Stabilizer.Homological.BB.Bravyi360.BottomForest'
f = open(p / 'Bravyi360BottomPositiveTwo.lean', 'w')
f.write('import Bravyi360BottomForest\n' + ''.join((f'import BottomPositiveTwoPart{i}\n' for i in range(14))) + f'set_option maxHeartbeats 0\nnamespace {ns}\nopen F2Mask FiberData PackedDigits\n')
f.write('/-- Every cap-two record in any canonical chunk has a complete outside-support list. -/\ntheorem twoRecords_complete (b : Fin 112) (r : Fiber) (hr : r ∈ recordsAt b.val)\n    (hweight : r.weight = 18) (word : Nat) (hw : word < 2^90)\n    (he : reconX bottomEColumns word = reconX bottomRColumns r.beta)\n    (hoff : bitCount 90 (outsideMask word r.beta) ≤ 2) :\n    outsideMask word r.beta ∈ r.seeds.map Seed.off := by\n  fin_cases b\n')
for i in range(112):
    k = i * 250
    f.write(f'  · change r ∈ records{k} at hr\n    exact twoGroup_complete bottomPairSpreadTable bottomPairRep bottomPair_lookup bottomPair_rep\n      records{k} twoCovers{k} twoPacked{k} checked{k} twoChecked{k} r hr hweight word hw he hoff\n')
f.write('/-- Every canonical cap-two fiber contains the outside support of each actual solution. -/\ntheorem twoFiberAt_complete (j : Fin 27867) (hweight : (fiberAt j.val).weight = 18)\n    (word : Nat) (hw : word < 2^90)\n    (he : reconX bottomEColumns word = reconX bottomRColumns (fiberAt j.val).beta)\n    (hoff : bitCount 90 (outsideMask word (fiberAt j.val).beta) ≤ 2) :\n    outsideMask word (fiberAt j.val).beta ∈ (fiberAt j.val).seeds.map Seed.off :=\n  twoRecords_complete ⟨j.val/250,by omega⟩ (fiberAt j.val) (fiberAt_mem j)\n    hweight word hw he hoff\n#print axioms twoFiberAt_complete\n' + f'end {ns}\n')
f.close()
