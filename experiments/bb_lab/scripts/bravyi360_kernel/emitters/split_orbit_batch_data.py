"""Emit kernel-checkable certificates into QEC360_STAGE."""
import os
from pathlib import Path
STAGE=Path(os.environ["QEC360_STAGE"]).resolve()
import re, json
s = open(str(STAGE / 'OrbitBatchComplete1024.lean')).read()
prefix = s[:s.index('def reps0 :')].replace('namespace OrbitBatchComplete', 'namespace OrbitBatchCheck')
lines = s[s.index('def reps0 :'):].splitlines()
data = []
proof = []
for line in lines:
    if line.startswith(('def reps', 'def words', 'def expected')):
        data.append(line)
    elif line.startswith('theorem checked'):
        proof.append(line)
head = 'import PackedDigits\nset_option maxRecDepth 400000\nnamespace OrbitBatchData\n/-- The geometric128-bit lane repetition word. -/\ndef repeat128 (size : Nat) : Nat := ((1 <<< (128*size))-1)/((1 <<< 128)-1)\n'
getters = []

def dispatch(lo, hi):
    if hi - lo == 1:
        return f'words{lo}'
    mid = (lo + hi) // 2
    return f'(if Nat.blt i {mid} then {dispatch(lo, mid)} else {dispatch(mid, hi)})'
for stem in ['expected']:
    getters.append(f'/-- Packed {stem} for each consecutive1024-element census batch. -/\ndef {stem}At : Nat → Nat\n' + ''.join((f'  | {i} => {stem}{i}\n' for i in range(268))) + '  | _ => 0\n')
getters.insert(0, '/-- Read a packed word batch by balanced index dispatch. -/\ndef wordsAt (i : Nat) : Nat := if Nat.blt i 268 then ' + dispatch(0, 268) + ' else 0\n')
getters.append('/-- The final batch has603 entries; all earlier batches have1024. -/\ndef sizeAt (i : Nat) : Nat := if i < 267 then 1024 else if i = 267 then 603 else 0\n')
getters.append('/-- The sorted canonical representative at a global zero-based census index. -/\ndef repAt (i : Nat) : Nat := PackedDigits.lane 128 (i%1024) (wordsAt (i/1024))\n')
open(str(STAGE / 'OrbitBatchData.lean'), 'w').write(head + '\n'.join(data) + '\n' + '\n'.join(getters) + '/-- The repeated lane constant has the expected packed semantics. -/\ntheorem repeat1024_eq : repeat128 1024 = PackedDigits.packF 128 1024 (fun _ => 1) := by decide +kernel\n/-- The final partial batch uses the same packed semantics. -/\ntheorem repeat603_eq : repeat128 603 = PackedDigits.packF 128 603 (fun _ => 1) := by decide +kernel\nend OrbitBatchData\n')
prefix = prefix.replace('import Mathlib.Data.Nat.Bitwise', 'import OrbitBatchData').replace('namespace OrbitBatchCheck', 'namespace OrbitBatchCheck\nopen OrbitBatchData')
prefix = prefix[:prefix.index('def repeat128')]
open(str(STAGE / 'OrbitBatchCheck.lean'), 'w').write(prefix + '\n'.join(('/-- This packed batch passes all finite orbit checks through the kernel. -/\n' + p for p in proof)) + '\n#print axioms checked267\nend OrbitBatchCheck\n')
print(len(data), len(proof), len(open(str(STAGE / 'OrbitBatchData.lean')).read()))
