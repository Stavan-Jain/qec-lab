"""Emit checked certificate data into the staging directory selected by QEC360_STAGE."""
import os
from pathlib import Path
STAGE = Path(os.environ["QEC360_STAGE"]).resolve()
import json
from pathlib import Path
p = json.load(open(str(STAGE / 'top_case_partition.json')))

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
ns = 'Quantum.Stabilizer.Homological.BB.Bravyi360.TopCaseData'
s = 'import Bravyi360TopFlatChecks\nimport Bravyi360TopEmptyData\nimport Bravyi360TopBank\n\nnamespace ' + ns + '\nopen F2Mask PackedDigits\nset_option maxRecDepth 200000\nset_option maxHeartbeats 0\n'
for i in range(57):
    s += f"/-- Family tags for top-bank chunk {i}. -/\ndef tags{i} : Nat := {literal(pack(p['tags'][1024 * i:1024 * (i + 1)], 18))}\n"
s += '/-- Select a packed chunk of family tags. -/\ndef tagsChunk (i : Nat) : Nat := ' + dispatch([f'tags{i}' for i in range(57)]) + '\n'
s += '/-- Read the family and local index of a top-bank entry. -/\ndef tagAt (i : Nat) : Nat := lane 18 (i % 1024) (tagsChunk (i / 1024))\n'
s += '/-- The sorted nonempty affine shadows. -/\ndef nonemptyWords : Nat := ' + literal(pack(p['nonempty'], 180)) + '\n'
s += '/-- Read a nonempty affine shadow. -/\ndef nonemptyBeta (i : Nat) : Nat := lane 180 i nonemptyWords\n'
s += '/-- The four exceptional empty affine shadows. -/\ndef exceptionalWords : Nat := ' + literal(pack(p['exceptions'], 180)) + '\n'
s += '/-- Read an exceptional empty affine shadow. -/\ndef exceptionalBeta (i : Nat) : Nat := lane 180 i exceptionalWords\n'
s += '/-- Number of shadows in each certified family. -/\ndef familySize : Nat → Nat\n  | 0 => 30338\n  | 1 => 649\n  | 2 => 4\n  | 3 => 24877\n  | 4 => 2119\n  | 5 => 143\n  | _ => 0\n/-- Read a shadow from its certified family. -/\ndef familyBeta (family i : Nat) : Nat :=\n  match family with\n  | 0 => TopFlatData.betaAt i\n  | 1 => nonemptyBeta i\n  | 2 => exceptionalBeta i\n  | 3 => TopEmpty.beta0 i\n  | 4 => TopEmpty.beta1 i\n  | 5 => TopEmpty.beta2 i\n  | _ => 0\n/-- Overflow allowance for each family of empty fibers. -/\ndef familyCap : Nat → Nat\n  | 2 => 3\n  | 3 => 0\n  | 4 => 1\n  | 5 => 2\n  | _ => 0\n/-- Byte-table Hamming weight for a 180-bit word. -/\ndef fast180 (m : Nat) : Nat := Nat.add (bitCountFast 22 m) (bitCount 4 (Nat.shiftRight m 176))\n/-- The byte-table implementation computes the specified Hamming weight. -/\ntheorem fast180_eq (m : Nat) : fast180 m = bitCount 180 m := by\n  exact (congrArg₂ Nat.add (bitCountFast_eq 22 m)\n    (congrArg (bitCount 4) (Nat.shiftRight_eq_div_pow m 176))).trans\n      (bitCount_add 176 4 m).symm\n/-- A top-bank tag identifies a valid shadow and, when needed, a sufficient overflow bound. -/\ndef Valid (i : Nat) : Prop :=\n  let family := tagAt i % 8\n  let j := tagAt i / 8\n  family < 6 ∧ j < familySize family ∧ topBankMask i = familyBeta family j ∧\n    (2 ≤ family → (22 - fast180 (topBankMask i)) / 2 ≤ familyCap family)\ninstance (i : Nat) : Decidable (Valid i) := inferInstanceAs (Decidable (_ ∧ _ ∧ _ ∧ _))\n'
s += 'end ' + ns + '\n'
Path(str(STAGE / 'Bravyi360TopCaseData.lean')).write_text(s)
s = 'import Bravyi360TopCaseData\nnamespace ' + ns + '\nset_option maxRecDepth 200000\nset_option maxHeartbeats 0\n'
s += '/-- A direct Boolean version of the family-association predicate. -/\ndef entryCheck (i : Nat) : Bool :=\n  let family := tagAt i % 8\n  let j := tagAt i / 8\n  let beta := topBankMask i\n  Nat.blt family 6 && Nat.blt j (familySize family) &&\n    Nat.beq beta (familyBeta family j) &&\n    (Nat.blt family 2 || Nat.ble ((22 - fast180 beta) / 2) (familyCap family))\n/-- Boolean family validation implies the mathematical family predicate. -/\ntheorem entryCheck_sound (i : Nat) (h : entryCheck i = true) : Valid i := by\n  simp only [entryCheck, Bool.and_eq_true, Bool.or_eq_true, Nat.blt_eq,\n    Nat.ble_eq, Nat.beq_eq] at h\n  exact ⟨h.1.1.1, h.1.1.2, h.1.2, fun hf => h.2.resolve_left (by omega)⟩\n'
s += 'end ' + ns + '\n'
Path(str(STAGE / 'Bravyi360TopCaseCheck.lean')).write_text(s)
for part in range(6):
    s = 'import Bravyi360TopCaseCheck\nnamespace ' + ns + '\nset_option maxRecDepth 200000\nset_option maxHeartbeats 0\n'
    for c in range(part * 10, min(57, (part + 1) * 10)):
        n = min(1024, 58130 - 1024 * c)
        s += f'/-- The {n} family associations in bank chunk {c} are kernel-checked. -/\ntheorem checked{c} : (List.range {n}).all (fun j => entryCheck ({1024 * c} + j)) = true := by\n  decide +kernel\n'
    s += 'end ' + ns + '\n'
    Path(f'{STAGE}/Bravyi360TopCasePart{part}.lean').write_text(s)
s = ''.join((f'import Bravyi360TopCasePart{j}\n' for j in range(6))) + 'namespace ' + ns + '\nset_option maxRecDepth 200000\nset_option maxHeartbeats 0\n'
s += '/-- Every packed chunk passes the family-association check. -/\ntheorem all_chunks (c : Fin 57) :\n    (List.range (min 1024 (58130 - 1024 * c.val))).all\n      (fun j => entryCheck (1024 * c.val + j)) = true := by\n  fin_cases c\n'
s += ''.join((f'  · exact checked{c}\n' for c in range(57)))
s += '/-- Every top-bank entry belongs to its certified family with the required overflow bound. -/\ntheorem all_checked (i : Fin 58130) : Valid i := by\n  have hc : i.val / 1024 < 57 := by omega\n  have hj : i.val % 1024 < min 1024 (58130 - 1024 * (i.val / 1024)) := by omega\n  have h := (List.all_eq_true.mp (all_chunks ⟨i.val / 1024, hc⟩))\n    (i.val % 1024) (List.mem_range.mpr hj)\n  have heq : 1024 * (i.val / 1024) + i.val % 1024 = i.val := by omega\n  change entryCheck (1024 * (i.val / 1024) + i.val % 1024) = true at h\n  rw [heq] at h\n  exact entryCheck_sound _ h\n'
s += 'end ' + ns + '\n'
Path(str(STAGE / 'Bravyi360TopCaseChecks.lean')).write_text(s)
print('Generated 58,130 top-bank partition tags and 57 numeric chunks.')

# Certificate constants are reduced by the kernel and do not need runtime code.
import re
for generated in ['Bravyi360TopCaseData.lean', 'Bravyi360TopCaseCheck.lean', 'Bravyi360TopCaseChecks.lean', 'Bravyi360TopCasePart0.lean', 'Bravyi360TopCasePart1.lean', 'Bravyi360TopCasePart2.lean', 'Bravyi360TopCasePart3.lean', 'Bravyi360TopCasePart4.lean', 'Bravyi360TopCasePart5.lean']:
    path = STAGE / generated
    text = path.read_text()
    text = re.sub(r"(?m)^def\b", "noncomputable def", text)
    text = re.sub(r"(?m)^instance\b", "noncomputable instance", text)
    path.write_text(text)
