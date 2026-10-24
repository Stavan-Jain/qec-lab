"""Emit checked certificate data into the staging directory selected by QEC360_STAGE."""
import os
from pathlib import Path
STAGE = Path(os.environ["QEC360_STAGE"]).resolve()
import json
from pathlib import Path
rows = [json.loads(l) for l in open(str(STAGE / 'qec360-certprobe/bottom_partition.jsonl'))]
deep = json.load(open(str(STAGE / 'qec360-certprobe/bottom_deep_order.json')))
assert [r['orbitIndex'] for r in rows] == list(range(274011))
for r in rows:
    if r['family'] == 5:
        assert r['beta'] == deep[r['familyIndex']]

def literal(n):
    if n.bit_length() <= 4096:
        return hex(n)
    bits = (n.bit_length() + 8191) // 8192 * 4096
    return f'(Nat.lor {literal(n % 2 ** bits)} (Nat.shiftLeft {literal(n >> bits)} {bits}))'

def pack(xs, w):
    return sum((x << w * i for i, x in enumerate(xs)))

def dispatch(names, off=0):
    if len(names) == 1:
        return names[0]
    half = len(names) // 2
    return f'(if i < {off + half} then {dispatch(names[:half], off)} else {dispatch(names[half:], off + half)})'
ns = 'Quantum.Stabilizer.Homological.BB.Bravyi360.BottomCase'
s = 'import OrbitBatchData\nnamespace ' + ns + '\nopen F2Mask PackedDigits\nset_option maxRecDepth 200000\nset_option maxHeartbeats 0\n'
for i in range(268):
    s += f"/-- Family tags for bottom-orbit chunk {i}. -/\ndef words{i} : Nat := {literal(pack([r['family'] + 8 * r['familyIndex'] for r in rows[1024 * i:1024 * (i + 1)]], 24))}\n"
s += '/-- Select one packed bottom-orbit tag chunk. -/\ndef wordsAt (i : Nat) : Nat := ' + dispatch([f'words{i}' for i in range(268)]) + '\n'
s += '/-- Read the family and local index of a bottom representative. -/\ndef entry (i : Nat) : Nat := lane 24 (i % 1024) (wordsAt (i / 1024))\n'
s += '/-- Select the local index within a certified bottom-fiber family. -/\ndef localIndex (i : Nat) : Nat := entry i / 8\n'
s += '/-- The twenty-eight deep bottom shadows in shared certificate order. -/\ndef deepWords : Nat := ' + literal(pack(deep, 90)) + '\n'
s += '/-- Read one deep bottom shadow. -/\ndef deepBeta (i : Nat) : Nat := lane 90 i deepWords\n'
deep_indices = [r['orbitIndex'] for r in sorted((r for r in rows if r['family'] == 5), key=lambda r: r['familyIndex'])]
assert len(deep_indices) == 28
for r in rows:
    raw = r['family']
    j = r['familyIndex']
    i = r['orbitIndex']
    guarded = (5 if j < 28 and i == deep_indices[j] else 4) if raw == 5 else min(raw, 4)
    assert guarded == raw
s += '/-- Census positions of the twenty-eight deep representatives. -/\ndef deepIndexWords : Nat := ' + literal(pack(deep_indices, 19)) + '\n'
s += '/-- Read one deep representative position. -/\ndef deepOrbitIndex (i : Nat) : Nat := lane 19 i deepIndexWords\n'
s += '/-- Select a certified family; deep tags additionally carry their sparse-index guard. -/\ndef tag (i : Nat) : Nat :=\n  if entry i % 8 = 5 then\n    if localIndex i < 28 ∧ i = deepOrbitIndex (localIndex i) then 5 else 4\n  else min (entry i % 8) 4\n'
s += 'end ' + ns + '\n'
Path(str(STAGE / 'Bravyi360BottomCaseData.lean')).write_text(s)
s = 'import Bravyi360BottomCaseData\nnamespace ' + ns + '\nset_option maxRecDepth 200000\nset_option maxHeartbeats 0\n'
s += '/-- The sparse deep positions contain their stated canonical shadows. -/\ntheorem deep_sparse_checked : ∀ j : Fin 28,\n    OrbitBatchData.repAt (deepOrbitIndex j) = deepBeta j := by decide +kernel\n/-- Every guarded deep tag points to its actual deep shadow. -/\ntheorem deep_checked (i : Fin 274011) (hd : tag i = 5) :\n    localIndex i < 28 ∧ OrbitBatchData.repAt i = deepBeta (localIndex i) := by\n  have hg : localIndex i < 28 ∧ i.val = deepOrbitIndex (localIndex i) := by\n    unfold tag at hd\n    split_ifs at hd with hraw hguard\n    · exact hguard\n    · contradiction\n    · have hm := Nat.min_le_right (entry i % 8) 4\n      omega\n  refine ⟨hg.1, ?_⟩\n  exact (congrArg OrbitBatchData.repAt hg.2).trans\n    (deep_sparse_checked ⟨localIndex i, hg.1⟩)\n'
s += 'end ' + ns + '\n'
Path(str(STAGE / 'Bravyi360BottomCaseChecks.lean')).write_text(s)
print('Generated bottom partition and28deep shadow linkage.')

# Certificate constants are reduced by the kernel and do not need runtime code.
import re
for generated in ['Bravyi360BottomCaseData.lean', 'Bravyi360BottomCaseChecks.lean']:
    path = STAGE / generated
    text = path.read_text()
    text = re.sub(r"(?m)^def\b", "noncomputable def", text)
    text = re.sub(r"(?m)^instance\b", "noncomputable instance", text)
    path.write_text(text)
