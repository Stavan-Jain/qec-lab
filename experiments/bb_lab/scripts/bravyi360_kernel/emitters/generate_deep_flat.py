"""Emit checked certificate data into the staging directory selected by QEC360_STAGE."""
import os
from pathlib import Path
STAGE = Path(os.environ["QEC360_STAGE"]).resolve()
import json
from pathlib import Path
num = lambda x: int(x, 0) if isinstance(x, str) else int(x)
D = json.load(open(str(STAGE / 'deep_bottom_certificates.json')))['fibers']
D = [r for r in D if r['type'] == 'flat']
assert len(D) == 26
C = json.load(open(str(STAGE / 'cycles10_normalization.json')))
print('normalization shape', list(C)[:10])
cycles = [num(r['mask']) for r in C['records']]
bank = [num(json.loads(l)['b']) for l in open(str(STAGE / 'all_top_rung_certificates.jsonl'))]
lookup = {b: i for i, b in enumerate(bank)}
A = [(9, 0), (0, 1), (0, 2)]
B = [(0, 0), (25, 0), (26, 0)]

def row(i):
    x, y = divmod(i, 3)
    return sum((1 << bl * 90 + (x - a) % 30 * 3 + (y - b) % 3 for bl, P in enumerate([A, B]) for a, b in P))

def sheet(v, k):
    return sum(((v >> bl * 90 + x * 3 + y + 45 * k & 1) << bl * 45 + x * 3 + y for bl in range(2) for x in range(15) for y in range(3)))
rows = [row(i) for i in range(90)]
s0 = [sheet(v, 0) for v in rows]
shadow = [sheet(v, 0) ^ sheet(v, 1) for v in rows]

def recon(cs, m):
    z = 0
    while m:
        low = m & -m
        z ^= cs[low.bit_length() - 1]
        m ^= low
    return z

def join(a, b):
    x = a ^ b
    low = (1 << 45) - 1
    return a & low ^ a >> 45 << 90 ^ (x & low) << 45 ^ x >> 45 << 135

def literal(n):
    if n.bit_length() <= 4096:
        return hex(n)
    bits = (n.bit_length() + 8191) // 8192 * 4096
    return f'(Nat.lor {literal(n % 2 ** bits)} (Nat.shiftLeft {literal(n >> bits)} {bits}))'
ns = 'Quantum.Stabilizer.Homological.BB.Bravyi360.DeepFlatData'
s = 'import Bravyi360MiddleSparse\n\nnamespace ' + ns + '\n'
for name, vals in [('sheetColumns', s0), ('shadowColumns', shadow)]:
    s += '/-- Packed columns for the middle boundary restriction. -/\ndef ' + name + ' : List (Nat × Nat) :=\n [' + ',\n  '.join((f'({i},0x{v:x})' for i, v in enumerate(vals))) + ']\n\n'
s += 'end ' + ns + '\n'
Path(str(STAGE / 'Bravyi360BottomBoundarySheetsData.lean')).write_text(s)
checked = 0
for i, r in enumerate(D):
    b = num(r['beta'])
    a = num(r['v0Flat'])
    c = num(r['coverGeneratorCoeff'])
    assert recon(s0, c) == a and recon(shadow, c) == b and (not a & ~b)
    tags = []
    for g in cycles:
        first = a ^ g
        if b.bit_count() + 2 * (first & ~b).bit_count() > 22:
            tags.append(0)
            continue
        w = join(first, b)
        if w in lookup:
            tag = 2 * lookup[w]
        else:
            w = join(first ^ b, b)
            assert w in lookup, (i, hex(g))
            tag = 2 * lookup[w] + 1
        tags.append(tag)
        checked += 1
    r['tags'] = tags
json.dump(D, open(str(STAGE / 'deep_flat_bank_certificates.json'), 'w'))
s = 'import Bravyi360BottomBoundarySheetsData\nimport Bravyi360Join\nimport Bravyi360TopBank\nimport LowCycleNormalization\nnamespace Quantum.Stabilizer.Homological.BB.Bravyi360.DeepFlatData\nopen F2Mask\nset_option maxRecDepth 200000\nset_option maxHeartbeats 0\n/-- Byte-table weight of a90-bit word, including the final two bits. -/\ndef fast90 (m : Nat) : Nat := Nat.add (bitCountFast 11 m) (bitCount 2 (Nat.shiftRight m 88))\n/-- The byte-table weight agrees with the90-bit Hamming-weight specification. -/\ntheorem fast90_eq (m : Nat) : fast90 m = bitCount 90 m := by\n  exact (congrArg₂ Nat.add (bitCountFast_eq 11 m)\n    (congrArg (bitCount 2) (Nat.shiftRight_eq_div_pow m 88))).trans\n      (bitCount_add 88 2 m).symm\n'
for name, key in [('beta', 'beta'), ('first', 'v0Flat'), ('coeff', 'coverGeneratorCoeff')]:
    packed = sum((num(r[key]) << 90 * i for i, r in enumerate(D)))
    s += f'/-- Packed {name} data for the26 deep flat shadows. -/\ndef {name}Packed : Nat := {literal(packed)}\n'
    s += f'/-- Read one deep-flat {name}. -/\ndef {name}At (i : Nat) : Nat := ({name}Packed >>> (90*i)) % 2^90\n'
for i, r in enumerate(D):
    packed = sum((t << 17 * j for j, t in enumerate(r['tags'])))
    s += f'/-- The top-bank targets of all1954 light cycles for deep shadow {i}. -/\ndef tags{i} : Nat := {literal(packed)}\n'
    s += f'/-- Read one target and deck choice for deep shadow {i}. -/\ndef tag{i} (j : Nat) : Nat := (tags{i} >>> (17*j)) % 2^17\n'
    s += f'/-- The small-cycle census accounts for every light lift of deep flat shadow {i}. -/\ntheorem checked{i} : ∀ j : Fin 1954,\n    let a := Nat.xor (firstAt {i}) (LowCycleNormalization.maskAt j)\n    let b := betaAt {i}\n    let t := tag{i} j\n    22 < fast90 b + 2 * fast90 (Nat.land a (Nat.xor (2^90-1) b)) ∨\n      (t / 2 < topBankSize ∧\n        bottomJoin (if t % 2 = 1 then a ^^^ b else a) b = topBankMask (t / 2)) := by\n  decide +kernel\n'
s += 'end ' + ns + '\n'
Path(str(STAGE / 'Bravyi360DeepFlatChecks.lean')).write_text(s)
print('Generated26deepflatcases,', checked, 'light-cycle/index checks.')

# Certificate constants are reduced by the kernel and do not need runtime code.
import re
for generated in ['Bravyi360BottomBoundarySheetsData.lean']:
    path = STAGE / generated
    text = path.read_text()
    text = re.sub(r"(?m)^def\b", "noncomputable def", text)
    text = re.sub(r"(?m)^instance\b", "noncomputable instance", text)
    path.write_text(text)
