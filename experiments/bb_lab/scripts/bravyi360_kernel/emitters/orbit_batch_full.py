"""Emit kernel-checkable certificates into QEC360_STAGE."""
import os
from pathlib import Path
STAGE=Path(os.environ["QEC360_STAGE"]).resolve()
import json, sys, time
sys.path.insert(0, str(STAGE))
from balanced_nat_literal import literal
rows = [json.loads(l) for l in open(str(STAGE / 'qec360-certprobe/gb_orbit_certificates.jsonl'))]
meta = json.load(open(str(STAGE / 'qec360-certprobe/gb_orbit_metadata.json')))['masks']
B = int(sys.argv[1]) if len(sys.argv) > 1 else 1024
prefix = open(str(STAGE / 'OrbitBatchProbe.lean')).read().split('def reps :')[0].replace('namespace OrbitBatchProbe', 'namespace OrbitBatchFull')
prefix = prefix.replace('((v &&& (r*hi1)) >>> 2)', '((v >>> 2) &&& (r*(hi1/4)))')
prefix = prefix.replace('((2^(45-s)-1)*(1+2^45)*r)', '(r*((2^(45-s)-1)*(1+2^45)))')
prefix = prefix.replace('(((2^s-1) <<< (45-s))*(1+2^45)*r)', '(r*(((2^s-1) <<< (45-s))*(1+2^45)))')
prefix = prefix.replace('if i < 15 then v else if i < 30 then yrot r v', 'if i/15 = 0 then v else if i/15 = 1 then yrot r v')
o = open(f'{STAGE}/OrbitBatchFull{B}.lean', 'w')
o.write(prefix)
o.write('def repeat128 (size : Nat) : Nat := ((1 <<< (128*size))-1)/((1 <<< 128)-1)\n')
for batch, offset in enumerate(range(0, len(rows), B)):
    group = rows[offset:offset + B]
    p = lambda xs: sum((x << 128 * i for i, x in enumerate(xs)))
    o.write(f'def reps{batch} : Nat := repeat128 {len(group)}\n')
    o.write(f"def words{batch} : Nat := {literal(p([r['mask'] for r in group]), 4096)}\n")
    extras = [f"(Nat.shiftLeft {meta[r['fixedMask']]['fixedPoints'] - 1} {128 * i})" for i, r in enumerate(group) if meta[r['fixedMask']]['fixedPoints'] > 1]
    expected = f'reps{batch}' + ''.join((' + ' + x for x in extras))
    o.write(f'def expected{batch} : Nat := {expected}\n')
    o.write(f'/-- This batch passes canonicality and stabilizer-size checks through the kernel. -/\ntheorem checked{batch} : batchCheck reps{batch} words{batch} expected{batch} = true := by decide +kernel\n')
o.write('end OrbitBatchFull\n')
o.close()
print(f'Generated {batch + 1} batches, {len(rows)} records, {STAGE}/OrbitBatchFull{B}.lean')
