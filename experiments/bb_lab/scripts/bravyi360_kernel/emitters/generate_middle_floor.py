"""Emit checked certificate data into the staging directory selected by QEC360_STAGE."""
import os
from pathlib import Path
STAGE = Path(os.environ["QEC360_STAGE"]).resolve()
from pathlib import Path
import json
D = json.load(open(str(STAGE / 'basefloor_certificates.json')))
F = json.load(open(str(STAGE / 'qec360-certprobe/gb_frame.json')))

def asint(x):
    return int(x, 0) if isinstance(x, str) else int(x)

def mask(v):
    return sum((1 << i for i in v))

def pack(v, w):
    return sum((asint(x) << w * i for i, x in enumerate(v)))

def emit(name, v, w):
    return f'/-- Packed {name} certificate data. -/\ndef {name}Packed : Nat := 0x{pack(v, w):x}\n/-- Read one {name} certificate. -/\ndef {name} (i : Nat) : Nat := ({name}Packed >>> ({w} * i)) % 2 ^ {w}\n\n'
A = [(9, 0), (0, 1), (0, 2)]
B = [(0, 0), (25, 0), (26, 0)]

def row(i):
    x, y = divmod(i, 3)
    return sum((1 << bl * 90 + (x - a) % 30 * 3 + (y - b) % 3 for bl, P in enumerate([A, B]) for a, b in P))

def recon(cols, m):
    r = 0
    for i, c in enumerate(cols):
        if m >> i & 1:
            r ^= asint(c)
    return r
E = list(map(asint, F['Ecols']))
R = list(map(asint, F['Rcols']))
empty = [mask(r['support']) for r in D['empty']]
flat = [mask(r['support']) for r in D['flat']]
for r, b in zip(D['empty'], empty):
    if r['cap'] == 0:
        continue
    rhs = recon(R, b)
    old = list(map(asint, r['lambdas']))
    choices = [recon(old, m) for m in range(1 << len(old))]
    choices = [v for v in choices if (v & rhs).bit_count() % 2]
    off = set(range(90)) - set(r['support'])
    selected = []
    while off:
        v = max(choices, key=lambda v: sum((not (v & E[q]).bit_count() % 2 for q in off)))
        covered = {q for q in off if not (v & E[q]).bit_count() % 2}
        assert covered
        selected.append(v)
        off -= covered
    r['lambdas'] = selected
assert all((len(r['lambdas']) <= 6 for r in D['empty']))
lam = [(list(map(asint, r['lambdas'])) + [0] * 6)[:6] for r in D['empty']]
s = 'import Bravyi360FiberData\n\nnamespace Quantum.Stabilizer.Homological.BB.Bravyi360.MiddleFloorData\nopen F2Mask\n'
s += emit('eColumn', E, 90) + emit('emptyBeta', empty, 90) + emit('emptyRhs', [recon(R, b) for b in empty], 90)
s += emit('emptyLambda', [x for row in lam for x in row], 90)
s += emit('emptyCap', [r['cap'] for r in D['empty']], 2)
s += emit('flatBeta', flat, 90) + emit('flatFirst', [asint(r['v0']) for r in D['flat']], 90)
s += emit('flatCoeff', [asint(r['coverGeneratorCoeff']) for r in D['flat']], 90)
rows = [row(i) for i in range(90)]
s += emit('boundaryColumn', rows, 180)
s += emit('flatMask', [recon(rows, asint(r['coverGeneratorCoeff'])) for r in D['flat']], 180)
s += 'end Quantum.Stabilizer.Homological.BB.Bravyi360.MiddleFloorData\n'
Path(str(STAGE / 'Bravyi360MiddleFloorData.lean')).write_text(s)

# Certificate constants are reduced by the kernel and do not need runtime code.
import re
for generated in ['Bravyi360MiddleFloorData.lean']:
    path = STAGE / generated
    text = path.read_text()
    text = re.sub(r"(?m)^def\b", "noncomputable def", text)
    text = re.sub(r"(?m)^instance\b", "noncomputable instance", text)
    path.write_text(text)
