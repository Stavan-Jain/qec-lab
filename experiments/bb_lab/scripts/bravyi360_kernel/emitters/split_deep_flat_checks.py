"""Emit checked certificate data into the staging directory selected by QEC360_STAGE."""
import os
from pathlib import Path
STAGE = Path(os.environ["QEC360_STAGE"]).resolve()
from pathlib import Path
import re
p = Path(str(STAGE / 'Bravyi360DeepFlatChecks.lean'))
s = p.read_text()
if not 'theorem checked0' in s:
    s = Path(str(STAGE / 'Bravyi360DeepFlatChecksUnsplit.txt')).read_text()
else:
    Path(str(STAGE / 'Bravyi360DeepFlatChecksUnsplit.txt')).write_text(s)
blocks = re.findall('/-- The small-cycle census accounts.*?  decide \\+kernel\\n', s, re.S)
assert len(blocks) == 26
for b in blocks:
    s = s.replace(b, '')
Path(str(STAGE / 'Bravyi360DeepFlatData.lean')).write_text(s)
ns = 'Quantum.Stabilizer.Homological.BB.Bravyi360.DeepFlatData'
imports = []
for k in range(5):
    name = f'Bravyi360DeepFlatPart{k}'
    imports.append('import ' + name)
    t = 'import Bravyi360DeepFlatData\n\nnamespace ' + ns + '\nset_option maxRecDepth 200000\nset_option maxHeartbeats 0\n'
    t += '\n'.join(blocks[k * 6:(k + 1) * 6])
    t += '\nend ' + ns + '\n'
    Path(str(STAGE / '') + '/' + name + '.lean').write_text(t)
p.write_text('\n'.join(imports) + '\n')
print('Split26 deep-flat checks into5 modules, shareddata separate.')

# Certificate constants are reduced by the kernel and do not need runtime code.
import re
for generated in ['Bravyi360DeepFlatData.lean', 'Bravyi360DeepFlatChecks.lean']:
    path = STAGE / generated
    text = path.read_text()
    text = re.sub(r"(?m)^def\b", "noncomputable def", text)
    text = re.sub(r"(?m)^instance\b", "noncomputable instance", text)
    path.write_text(text)
