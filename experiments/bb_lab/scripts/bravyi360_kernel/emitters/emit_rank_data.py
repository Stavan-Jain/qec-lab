"""Emit kernel-checkable certificates into QEC360_STAGE."""
import os
from pathlib import Path
STAGE=Path(os.environ["QEC360_STAGE"]).resolve()
"""Emit packed Lean rank witnesses; all emitted facts are checked by Lean separately."""
import json
D = json.load(open(str(STAGE / 'rank174_reflected_cert.json')))
with open(str(STAGE / 'Bravyi360RankData.lean'), 'w') as o:
    o.write('import QEC.Stabilizer.Framework.Homological.BZCertificate\n\nnamespace Quantum.Stabilizer.Homological.BB.Bravyi360.RankData\n\n')
    for key in ['X', 'Z']:
        d = D[key]
        basis = [int(x, 16) for x in d['basisRows']]
        dec = [int(x, 16) for x in d['decoderColumns']]
        cs = [int(x, 16) for x in d['allRowsCoefficients']]

        def reconstruct(c):
            z = 0
            for i, v in enumerate(basis):
                if c >> i & 1:
                    z ^= v
            return z
        for name, ls, wd in [('basis', basis, 360), ('decoder', dec, 174), ('coeff', cs, 174), ('allRows', list(map(reconstruct, cs)), 360)]:
            packed = sum((v << wd * i for i, v in enumerate(ls)))
            o.write(f'/-- Packed {key}-check {name} certificate. -/\ndef {name}{key} : Nat := {hex(packed)}\n\n')
        for name, wd in [('basis', 360), ('decoder', 174), ('coeff', 174), ('allRows', 360)]:
            o.write(f'/-- Read one {key}-check {name} certificate entry. -/\ndef {name}{key}At (i : Nat) : Nat := ({name}{key} >>> ({wd} * i)) % 2 ^ {wd}\n\n')
        for name, n in [('basis', 174), ('decoder', 360)]:
            o.write(f'/-- Column list used by the {key}-check {name} map. -/\ndef {name}{key}List : List (Nat × Nat) := List.ofFn (fun i : Fin {n} => (i.val, {name}{key}At i.val))\n\n')
    w = D['witness24']
    o.write(f"/-- An explicit weight24 nontrivial cycle in certificate coordinates. -/\ndef witness : Nat := {w['mask']}\n\n/-- Dual certificate pairing oddly with the weight24 cycle. -/\ndef witnessDual : Nat := {w['dualSeparatingMask']}\n\nend Quantum.Stabilizer.Homological.BB.Bravyi360.RankData\n")
