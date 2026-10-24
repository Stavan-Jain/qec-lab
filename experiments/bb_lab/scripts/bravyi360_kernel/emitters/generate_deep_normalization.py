"""Emit checked certificate data into the staging directory selected by QEC360_STAGE."""
import os
from pathlib import Path
STAGE = Path(os.environ["QEC360_STAGE"]).resolve()
import json
from pathlib import Path
can = json.load(open(str(STAGE / 'qec360-certprobe/bottom_deep_order.json')))
cert = [int(r['beta'], 0) for r in json.load(open(str(STAGE / 'deep_flat_bank_certificates.json')))] + [135437471593392630495445065, 19343403427685214168604673]

def translate(m, x, y):
    return sum(((m >> b * 45 + a * 3 + c & 1) << b * 45 + (a + x) % 15 * 3 + (c + y) % 3 for b in range(2) for a in range(15) for c in range(3)))
rows = []
for b in can:
    matches = [(cert.index(z), x, y) for x in range(15) for y in range(3) if (z := translate(b, x, y)) in cert]
    assert matches
    rows.append(matches[0])
assert len({r[0] for r in rows}) == 28
json.dump(rows, open(str(STAGE / 'deep_normalization.json'), 'w'))
ns = 'Quantum.Stabilizer.Homological.BB.Bravyi360.DeepNormalization'
s = 'import Bravyi360BottomCaseData\nimport DeepFiberData\nimport TranslationNormalization\n\nnamespace ' + ns + '\nopen F2Mask\nset_option maxRecDepth 200000\nset_option maxHeartbeats 0\n'
for name, idx, w in [('indexAt', 0, 5), ('xAt', 1, 4), ('yAt', 2, 2)]:
    packed = sum((r[idx] << w * i for i, r in enumerate(rows)))
    s += f'/-- Packed {name} choices for the deep representative normalization. -/\ndef {name} (i : Nat) : Nat := (0x{packed:x} >>> ({w} * i)) % 2 ^ {w}\n'
s += '/-- Translate a canonical deep representative to the checked affine representative. -/\ndef displacement (i : Nat) : BaseGroup := (xAt i, yAt i)\n'
s += '/-- Each deep representative has a checked translation into the affine certificate order. -/\ntheorem checked : ∀ i : Fin 28, indexAt i < 28 ∧\n    maskOf baseQubitEnum (translate1 (displacement i)\n      (chainOfMask baseQubitEnum (BottomCase.deepBeta i))) = DeepFiberData.betaAt (indexAt i) := by\n  decide +kernel\n'
s += 'end ' + ns + '\n'
Path(str(STAGE / 'Bravyi360DeepNormalization.lean')).write_text(s)
print('Verified28canonical-to-affine translations; nonzero translations:', sum((x != 0 or y != 0 for _, x, y in rows)))

# Certificate constants are reduced by the kernel and do not need runtime code.
import re
for generated in ['Bravyi360DeepNormalization.lean']:
    path = STAGE / generated
    text = path.read_text()
    text = re.sub(r"(?m)^def\b", "noncomputable def", text)
    text = re.sub(r"(?m)^instance\b", "noncomputable instance", text)
    path.write_text(text)
