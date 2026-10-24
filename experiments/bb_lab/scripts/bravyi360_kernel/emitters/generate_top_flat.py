"""Emit checked certificate data into the staging directory selected by QEC360_STAGE."""
import os
from pathlib import Path
STAGE = Path(os.environ["QEC360_STAGE"]).resolve()
import json
from pathlib import Path

def num(x):
    return int(x, 0) if isinstance(x, str) else int(x)

def recon(cols, m):
    z = 0
    while m:
        low = m & -m
        z ^= cols[low.bit_length() - 1]
        m ^= low
    return z
A = [(9, 0), (0, 1), (0, 2)]
B = [(0, 3), (25, 0), (26, 0)]

def row(i):
    x, y = divmod(i, 6)
    return sum((1 << bl * 180 + (x - a) % 30 * 6 + (y - b) % 6 for bl, P in enumerate([A, B]) for a, b in P))

def sheet(v, k):
    return sum(((v >> bl * 180 + x * 6 + y + 3 * k & 1) << bl * 90 + x * 3 + y for bl in range(2) for x in range(30) for y in range(3)))
rows = [row(i) for i in range(180)]
s0 = [sheet(v, 0) for v in rows]
shadow = [sheet(v, 0) ^ sheet(v, 1) for v in rows]
flat = []
for i, l in enumerate(open(str(STAGE / 'all_top_rung_certificates.jsonl'))):
    r = json.loads(l)
    if r['type'] != 'flat':
        continue
    b = num(r['b'])
    c = num(r['coverGenerators'])
    a = recon(s0, c)
    assert recon(shadow, c) == b and (not a & ~b)
    assert a == num(r['v0'])
    flat.append({'index': i, 'beta': b, 'coeff': c})
assert len(flat) == 30338
json.dump(flat, open(str(STAGE / 'top_flat_boundary_certificates.json'), 'w'))
ns = 'Quantum.Stabilizer.Homological.BB.Bravyi360.TopFlatData'
s = 'import Bravyi360Defs\n\nnamespace ' + ns + '\n'
for name, vals in [('sheetColumns', s0), ('shadowColumns', shadow)]:
    s += '/-- Sparse packed columns of the top boundary restriction. -/\ndef ' + name + ' : List (Nat × Nat) :=\n  [' + ',\n   '.join((f'({i},0x{v:x})' for i, v in enumerate(vals))) + ']\n\n'
s += 'end ' + ns + '\n'
Path(str(STAGE / 'Bravyi360TopBoundarySheetsData.lean')).write_text(s)

def literal(n):
    if n.bit_length() <= 4096:
        return hex(n)
    bits = (n.bit_length() + 8191) // 8192 * 4096
    return f'(Nat.lor {literal(n % 2 ** bits)} (Nat.shiftLeft {literal(n >> bits)} {bits}))'

def pack(xs, width=180):
    return sum((x << width * i for i, x in enumerate(xs)))
s = 'import Bravyi360TopBoundarySheetsData\nimport BatchLinear\nnamespace Quantum.Stabilizer.Homological.BB.Bravyi360.TopFlatData\nopen F2Mask PackedDigits\nset_option maxRecDepth 200000\nset_option maxHeartbeats 0\n/-- Simultaneously check a batch of flat boundary witnesses. -/\ndef checkFlatBatch (coeff beta : Nat) : Bool :=\n  let rep := repeatBlock 1 180 10\n  let first := batchRecon rep coeff sheetColumns\n  let shadow := batchRecon rep coeff shadowColumns\n  Nat.beq shadow beta && Nat.beq (Nat.land first (Nat.xor (2 ^ (180 * 1024) - 1) beta)) 0\n'
for c in range(30):
    sub = flat[1024 * c:1024 * (c + 1)]
    cp = pack([r['coeff'] for r in sub])
    bp = pack([r['beta'] for r in sub])
    s += f'/-- Packed coefficients for flat boundary batch {c}. -/\ndef coeff{c} : Nat := {literal(cp)}\n'
    s += f'/-- Packed shadows for flat boundary batch {c}. -/\ndef beta{c} : Nat := {literal(bp)}\n'
    s += f'/-- Every witness in flat boundary batch {c} has its specified shadow and support. -/\ntheorem checked{c} : checkFlatBatch coeff{c} beta{c} = true := by decide +kernel\n'

def dispatch(names, off=0):
    if len(names) == 1:
        return names[0]
    half = len(names) // 2
    return f'(if i < {off + half} then {dispatch(names[:half], off)} else {dispatch(names[half:], off + half)})'
for c in range(30):
    sub = flat[1024 * c:1024 * (c + 1)]
    s += f"/-- Original top-bank indices for flat batch {c}. -/\ndef index{c} : Nat := {literal(pack([r['index'] for r in sub], 16))}\n"
for name, width in [('coeff', 180), ('beta', 180), ('index', 16)]:
    s += f'/-- Select one packed {name} batch. -/\ndef {name}Chunk (i : Nat) : Nat := ' + dispatch([f'{name}{c}' for c in range(30)]) + '\n'
    s += f'/-- Read one flat-witness {name}. -/\ndef {name}At (i : Nat) : Nat := lane {width} (i % 1024) ({name}Chunk (i / 1024))\n'
s += '/-- Every flat witness batch fits its1024 lanes. -/\ntheorem coeffChunks_bounded : ∀ i : Fin 30, coeffChunk i < 2 ^ (180 * 1024) := by decide +kernel\n'
s += '/-- Every packed flat witness batch passes the two matrix checks. -/\ntheorem all_checked (i : Fin 30) : checkFlatBatch (coeffChunk i) (betaChunk i) = true := by\n  fin_cases i\n'
for c in range(30):
    s += f'  · exact checked{c}\n'
s += 'end ' + ns + '\n'
Path(str(STAGE / 'Bravyi360TopFlatChecks.lean')).write_text(s)
print('Generated', len(flat), 'flat boundary witnesses and30packed batches.')

# Certificate constants are reduced by the kernel and do not need runtime code.
import re
for generated in ['Bravyi360TopBoundarySheetsData.lean', 'Bravyi360TopFlatChecks.lean']:
    path = STAGE / generated
    text = path.read_text()
    text = re.sub(r"(?m)^def\b", "noncomputable def", text)
    text = re.sub(r"(?m)^instance\b", "noncomputable instance", text)
    path.write_text(text)
