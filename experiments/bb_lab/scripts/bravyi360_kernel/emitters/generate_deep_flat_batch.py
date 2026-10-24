"""Emit checked certificate data into the staging directory selected by QEC360_STAGE."""
import os
from pathlib import Path
STAGE = Path(os.environ["QEC360_STAGE"]).resolve()
import json
from pathlib import Path
num = lambda x: int(x, 0) if isinstance(x, str) else int(x)
D = json.load(open(str(STAGE / 'deep_flat_bank_certificates.json')))
cycles = [num(r['mask']) for r in json.load(open(str(STAGE / 'cycles10_normalization.json')))['records']]
bank = [num(json.loads(l)['b']) for l in open(str(STAGE / 'all_top_rung_certificates.jsonl'))]
light = []
refs = []
for r in D:
    a = num(r['v0Flat'])
    b = num(r['beta'])
    rr = []
    for j, g in enumerate(cycles):
        f = a ^ g
        t = r['tags'][j]
        if b.bit_count() + 2 * (f & ~b).bit_count() > 22:
            rr.append(0)
            continue
        rr.append(len(light))
        light.append((f ^ b if t % 2 else f, b, t // 2, bank[t // 2]))
    refs.append(rr)
assert len(light) == 2549

def literal(n):
    if n.bit_length() <= 4096:
        return hex(n)
    bits = (n.bit_length() + 8191) // 8192 * 4096
    return f'(Nat.lor {literal(n % 2 ** bits)} (Nat.shiftLeft {literal(n >> bits)} {bits}))'

def pack(xs, w):
    return sum((x << w * i for i, x in enumerate(xs)))
ns = 'Quantum.Stabilizer.Homological.BB.Bravyi360.DeepFlatBatch'
s = 'import Bravyi360DeepFlatData\nimport BatchReconLane\nnamespace ' + ns + '\nopen F2Mask PackedDigits\nset_option maxRecDepth 200000\nset_option maxHeartbeats 0\nset_option exponentiation.threshold 1000000\n'
for key, col, w in [('first', 0, 180), ('beta', 1, 180), ('index', 2, 16), ('bank', 3, 180)]:
    for c in range(3):
        s += f'/-- Packed {key} fields for deep light-outcome batch {c}. -/\ndef {key}{c} : Nat := {literal(pack([r[col] for r in light[1024 * c:1024 * (c + 1)]], w))}\n'
    s += f'/-- Select a deep light-outcome {key} batch. -/\ndef {key}Chunk (i : Nat) : Nat := if i = 0 then {key}0 else if i = 1 then {key}1 else {key}2\n'
    s += f'/-- Read a deep light-outcome {key}. -/\ndef {key}At (i : Nat) : Nat := lane {w} (i % 1024) ({key}Chunk (i / 1024))\n'
for i in range(26):
    s += f'/-- References from the base census to deep light outcomes for case {i}. -/\ndef refs{i} : Nat := {literal(pack(refs[i], 12))}\n'
    s += f'/-- Read a deep light-outcome reference for case {i}. -/\ndef ref{i} (j : Nat) : Nat := lane 12 j refs{i}\n'
s += '/-- Both bottom embedding matrices fit the common packed lane width. -/\ntheorem columns_bounded :\n    (∀ p ∈ bottomLeft, p.1 < 180 ∧ p.2 < 2 ^ 180) ∧\n    (∀ p ∈ bottomRight, p.1 < 180 ∧ p.2 < 2 ^ 180) := by decide +kernel\n/-- Packed deep input lanes fit their three blocks. -/\ntheorem inputs_bounded : ∀ i : Fin 3,\n    firstChunk i < 2 ^ (180 * 1024) ∧ betaChunk i < 2 ^ (180 * 1024) := by decide +kernel\n/-- All light deep outcomes reconstruct their indicated top masks in parallel. -/\ntheorem joins_checked : ∀ i : Fin 3,\n    batchRecon (repeatBlock 1 180 10) (firstChunk i) bottomLeft ^^^\n      batchRecon (repeatBlock 1 180 10) (firstChunk i ^^^ betaChunk i) bottomRight = bankChunk i := by\n  decide +kernel\n/-- Each light deep outcome points to its exact common top-bank row. -/\ntheorem bank_checked : ∀ i : Fin 2549,\n    indexAt i < topBankSize ∧ bankAt i = topBankMask (indexAt i) := by decide +kernel\n/-- The packed matrix equality proves the physical join for each light deep outcome. -/\ntheorem result_checked (i : Fin 2549) :\n    indexAt i < topBankSize ∧ bottomJoin (firstAt i) (betaAt i) = topBankMask (indexAt i) := by\n  let b : Fin 3 := ⟨i.val / 1024, by omega⟩\n  have ha := (inputs_bounded b).1\n  have hb := (inputs_bounded b).2\n  have he := congrArg (lane 180 (i.val % 1024)) (joins_checked b)\n  have hrep : repeatBlock 1 180 10 = packF 180 1024 (fun _ => 1) :=\n    repeatBlock_eq_packF 1 180 10 (by decide) (by decide)\n  rw [lane_xor, hrep,\n    batchRecon_lane 180 1024 (firstChunk b) (by decide) ha bottomLeft columns_bounded.1\n      (i.val % 1024) (Nat.mod_lt _ (by decide)),\n    batchRecon_lane 180 1024 (firstChunk b ^^^ betaChunk b) (by decide)\n      (Nat.xor_lt_two_pow ha hb) bottomRight columns_bounded.2\n      (i.val % 1024) (Nat.mod_lt _ (by decide)), lane_xor] at he\n  exact ⟨(bank_checked i).1, he.trans (bank_checked i).2⟩\n'
s += 'end ' + ns + '\n'
Path(str(STAGE / 'Bravyi360DeepFlatBatch.lean')).write_text(s)
for part in range(5):
    s = 'import Bravyi360DeepFlatBatch\nnamespace Quantum.Stabilizer.Homological.BB.Bravyi360.DeepFlatData\nopen F2Mask\nset_option maxRecDepth 200000\nset_option maxHeartbeats 0\n'
    for i in range(part * 6, min(part * 6 + 6, 26)):
        s += f'/-- The base census selects either a heavy lift or a verified packed light outcome. -/\ntheorem selection{i} : ∀ j : Fin 1954,\n    let a := Nat.xor (firstAt {i}) (LowCycleNormalization.maskAt j)\n    let b := betaAt {i}\n    let t := tag{i} j\n    let k := DeepFlatBatch.ref{i} j\n    22 < fast90 b + 2 * fast90 (Nat.land a (Nat.xor (2^90-1) b)) ∨\n      (k < 2549 ∧ DeepFlatBatch.firstAt k = (if t % 2 = 1 then a ^^^ b else a) ∧\n        DeepFlatBatch.betaAt k = b ∧ DeepFlatBatch.indexAt k = t / 2) := by decide +kernel\n/-- The small-cycle census accounts for every light lift of deep flat shadow {i}. -/\ntheorem checked{i} (j : Fin 1954) :\n    let a := Nat.xor (firstAt {i}) (LowCycleNormalization.maskAt j)\n    let b := betaAt {i}\n    let t := tag{i} j\n    22 < fast90 b + 2 * fast90 (Nat.land a (Nat.xor (2^90-1) b)) ∨\n      (t / 2 < topBankSize ∧\n        bottomJoin (if t % 2 = 1 then a ^^^ b else a) b = topBankMask (t / 2)) := by\n  rcases selection{i} j with h | ⟨hk, ha, hb, ht⟩\n  · exact Or.inl h\n  · right\n    have h := DeepFlatBatch.result_checked ⟨DeepFlatBatch.ref{i} j, hk⟩\n    simpa only [ha, hb, ht] using h\n'
    s += 'end Quantum.Stabilizer.Homological.BB.Bravyi360.DeepFlatData\n'
    Path(f'{STAGE}/Bravyi360DeepFlatPart{part}.lean').write_text(s)
print('Generated packed verification of2549light joins plus26light/heavy selection checks.')

# Certificate constants are reduced by the kernel and do not need runtime code.
import re
for generated in ['Bravyi360DeepFlatBatch.lean', 'Bravyi360DeepFlatPart0.lean', 'Bravyi360DeepFlatPart1.lean', 'Bravyi360DeepFlatPart2.lean', 'Bravyi360DeepFlatPart3.lean', 'Bravyi360DeepFlatPart4.lean']:
    path = STAGE / generated
    text = path.read_text()
    text = re.sub(r"(?m)^def\b", "noncomputable def", text)
    text = re.sub(r"(?m)^instance\b", "noncomputable instance", text)
    path.write_text(text)
