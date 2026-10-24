"""Emit checked certificate data into the staging directory selected by QEC360_STAGE."""
import os
from pathlib import Path
STAGE = Path(os.environ["QEC360_STAGE"]).resolve()
from pathlib import Path
import json
D = json.load(open(str(STAGE / 'diagonal_certificates.json')))
C = json.load(open(str(STAGE / 'diagonal_cap_one_duals.json')))
F = json.load(open(str(STAGE / 'qec360-certprobe/top_frame.json')))

def num(x):
    return int(x, 0) if isinstance(x, str) else int(x)

def mask(v):
    return sum((1 << i for i in v))

def pack(v, w):
    return sum((num(x) << w * i for i, x in enumerate(v)))

def emit(name, v, w):
    return f'/-- Packed {name} certificates. -/\ndef {name}Packed : Nat := 0x{pack(v, w):x}\n/-- Read one {name} certificate. -/\ndef {name} (i : Nat) : Nat := ({name}Packed >>> ({w} * i)) % 2 ^ {w}\n\n'

def recon(cols, m):
    z = 0
    for i, c in enumerate(cols):
        if m >> i & 1:
            z ^= num(c)
    return z
E = list(map(num, F['Ecols']))
R = list(map(num, F['Rcols']))
beta = [num(r['beta']) for r in C]
gam = [num(r['gamma']) for r in C]
lam = [(list(map(num, r['lambdas'])) + [0] * 21)[:21] for r in C]
s = 'import Bravyi360FiberData\n\nnamespace Quantum.Stabilizer.Homological.BB.Bravyi360.DiagonalData\nopen F2Mask\n'
s += emit('eColumn', E, 180) + emit('emptyBeta', beta, 180) + emit('emptyRhs', [recon(R, b) for b in beta], 180) + emit('emptyGamma', gam, 90) + emit('emptyLambda', [x for r in lam for x in r], 180)
flat = D['flatStabilizers']
s += emit('flatGamma', [mask(r['gammaSupport']) for r in flat], 90) + emit('flatBeta', [mask(r['support']) for r in flat], 180) + emit('flatFirst', [num(r['v0']) for r in flat], 180) + emit('flatCoeff', [num(r['coverGeneratorCoeff']) for r in flat], 180)
A = [(9, 0), (0, 1), (0, 2)]
B = [(0, 3), (25, 0), (26, 0)]

def row(i):
    x, y = divmod(i, 6)
    return sum((1 << bl * 180 + (x - a) % 30 * 6 + (y - b) % 6 for bl, P in enumerate([A, B]) for a, b in P))
rows = [row(i) for i in range(180)]
s += emit('flatMask', [recon(rows, num(r['coverGeneratorCoeff'])) for r in flat], 360)
s += 'end Quantum.Stabilizer.Homological.BB.Bravyi360.DiagonalData\n'
Path(str(STAGE / 'Bravyi360DiagonalData.lean')).write_text(s)

# Certificate constants are reduced by the kernel and do not need runtime code.
import re
for generated in ['Bravyi360DiagonalData.lean']:
    path = STAGE / generated
    text = path.read_text()
    text = re.sub(r"(?m)^def\b", "noncomputable def", text)
    text = re.sub(r"(?m)^instance\b", "noncomputable instance", text)
    path.write_text(text)
