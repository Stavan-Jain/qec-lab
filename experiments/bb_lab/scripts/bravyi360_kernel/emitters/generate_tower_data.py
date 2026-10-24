"""Emit checked certificate data into the staging directory selected by QEC360_STAGE."""
import os
from pathlib import Path
STAGE = Path(os.environ["QEC360_STAGE"]).resolve()
from pathlib import Path
import json

def num(x):
    return int(x, 0) if isinstance(x, str) else x

def columns(name, values, width, description):
    out = f'/-- {description} -/\ndef {name} : List (Nat × Nat) := [\n'
    out += ',\n'.join((f'  ({i}, {num(v)})' for i, v in enumerate(values))) + '\n]\n\n'
    description = 'Every column fits the corresponding cover syndrome space.' if name.startswith(('bottom', 'top')) else 'Every output column fits its coordinate space.'
    out += f'/-- {description} -/\ntheorem {name}_bounded : ∀ p ∈ {name}, p.2 < 2 ^ {width} := by decide +kernel\n\n'
    return out

def transpose(rows, n):
    return [sum(((num(r) >> j & 1) << i for i, r in enumerate(rows))) for j in range(n)]
s = 'import QEC.Stabilizer.Framework.Homological.F2Mask\n\nnamespace Quantum.Stabilizer.Homological.BB.Bravyi360.FiberData\n\n'
for level, width, frame in [('bottom', 90, 'gb_frame'), ('top', 180, 'top_frame')]:
    d = json.load(open(str(STAGE / 'qec360-certprobe/') + '/' + frame + '.json'))
    for label in ['E', 'R']:
        s += columns(level + label + 'Columns', d[label + 'cols'], width, f'Packed columns of the {level} {label} lift equation.')
s += 'end Quantum.Stabilizer.Homological.BB.Bravyi360.FiberData\n'
Path(str(STAGE / 'Bravyi360FiberData.lean')).write_text(s)
d = json.load(open(str(STAGE / 'tower_direct43_cert.json')))
s = 'import QEC.Stabilizer.Framework.Homological.BZCertificate\n\nnamespace Quantum.Stabilizer.Homological.BB.Bravyi360.ReachabilityData\n\n'
for name, vals, w in [('reachableGenerators', d['G43Rows'], 90), ('reachableDecoder', transpose(d['decoderRows'], 90), 43), ('reachableCorrection', transpose(d['coverCheckCoefficients'], 180), 90)]:
    s += columns(name, vals, w, f'Packed columns of the {name} linear map.')
s += 'end Quantum.Stabilizer.Homological.BB.Bravyi360.ReachabilityData\n'
Path(str(STAGE / 'Bravyi360ReachabilityData.lean')).write_text(s)

# Certificate constants are reduced by the kernel and do not need runtime code.
import re
for generated in ['Bravyi360FiberData.lean', 'Bravyi360ReachabilityData.lean']:
    path = STAGE / generated
    text = path.read_text()
    text = re.sub(r"(?m)^def\b", "noncomputable def", text)
    text = re.sub(r"(?m)^instance\b", "noncomputable instance", text)
    path.write_text(text)
