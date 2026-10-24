"""Emit kernel-checkable certificates into QEC360_STAGE."""
import os
from pathlib import Path
STAGE=Path(os.environ["QEC360_STAGE"]).resolve()
import json, sys
sys.path.insert(0, str(STAGE))
from balanced_nat_literal import literal
D = json.load(open(str(STAGE / 'cycles10_normalization.json')))
R = D['records']
p = lambda xs, w: sum((x << w * i for i, x in enumerate(xs)))
s = 'import GBReflectedAction\nnamespace Quantum.Stabilizer.Homological.BB.Bravyi360.LowCycleNormalization\nopen F2Mask\nset_option maxRecDepth 200000\nset_option maxHeartbeats 0\nset_option exponentiation.threshold 400\n'
for name, n in [('cyclesPacked', p([r['mask'] for r in R], 128)), ('repsPacked', p(D['representatives'], 128)), ('lookupPacked', p([r['repIndex'] + 64 * r['x'] + 1024 * r['y'] for r in R], 16))]:
    s += f'def {name} : Nat := {literal(n, 4096)}\n'
s += '/-- The sorted complete list of light bottom cycles. -/\ndef maskAt (i : Nat) : Nat := (cyclesPacked >>> (128*i)) &&& (2^90-1)\n/-- The45 nonzero representatives, logical classes first then stabilizers. -/\ndef repAt (i : Nat) : Nat := (repsPacked >>> (128*i)) &&& (2^90-1)\n/-- The representative index selected for a nonzero light cycle. -/\ndef repIndex (i : Nat) : Nat := (lookupPacked >>> (16*i)) &&& 63\n/-- The x parameter in the established numerical orbit convention. -/\ndef xAt (i : Nat) : Nat := (lookupPacked >>> (16*i+6)) &&& 15\n/-- The y parameter in the established numerical orbit convention. -/\ndef yAt (i : Nat) : Nat := (lookupPacked >>> (16*i+10)) &&& 3\n/-- Every listed nonzero light cycle has a checked translation to one representative. -/\ntheorem checked : ∀ i : Fin 1954, maskAt i.val = 0 ∨\n    repIndex i.val < 45 ∧ xAt i.val < 15 ∧ yAt i.val < 3 ∧\n    GBOrbit.entry (maskAt i.val % 2^45) (maskAt i.val / 2^45)\n      (15*yAt i.val+xAt i.val) = repAt (repIndex i.val) := by\n  decide +kernel\n\n/-- Every chain with a listed mask is zero or translates to its certified representative. -/\ntheorem normalize (i : Fin 1954) (v : BaseGroup × Fin 2 → ZMod 2)\n    (hm : maskOf baseQubitEnum v = maskAt i.val) :\n    v = 0 ∨ ∃ (j : Fin 45) (g : BaseGroup),\n      maskOf baseQubitEnum (translate1 g v) = repAt j.val := by\n  rcases checked i with hz | ⟨hr,hx,hy,ht⟩\n  · left\n    apply maskOf_injective baseQubitEnum\n    rw [hm, hz, maskOf_zero]\n  · right\n    let x : ZMod 15 := ⟨xAt i.val,hx⟩\n    let y : ZMod 3 := ⟨yAt i.val,hy⟩\n    refine ⟨⟨repIndex i.val,hr⟩,(-x,y),?_⟩\n    rw [← entry_reflected_mask v x y, hm]\n    exact ht\n\nend Quantum.Stabilizer.Homological.BB.Bravyi360.LowCycleNormalization\n'
open(str(STAGE / 'LowCycleNormalization.lean'), 'w').write(s)
print('Generated1954 lookupcertificate', len(s), 'bytes')
