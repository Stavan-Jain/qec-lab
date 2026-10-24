"""Emit kernel-checkable certificates into QEC360_STAGE."""
import os
from pathlib import Path
STAGE=Path(os.environ["QEC360_STAGE"]).resolve()
import re, json, sys
p = str(STAGE / 'OrbitBatchFull1024.lean')
s = open(p).read()
M = [sum(((1 << (1 << k)) - 1 << j * (1 << k + 1) for j in range(128 // (1 << k + 1)))) for k in range(7)]
extra = 'def pcMask (k : Nat) : Nat := ' + (' if k = 0 then ' + str(M[0]) + ''.join((' else if k = ' + str(i) + ' then ' + str(M[i]) for i in range(1, 6))) + ' else ' + str(M[6])) + '\n'
extra += 'def popcount128 (r v : Nat) : Nat := (List.range 7).foldl (fun v k => (v &&& (r*pcMask k)) + ((v >>> (2^k)) &&& (r*pcMask k))) v\n'
extra += 'def extraCheck (r v : Nat) : Bool :=\n  let mask := r*(2^90-1)\n  let weights := popcount128 r v\n  let ps := r <<< 7\n  let prev := r >>> 128\n  let s := prev <<< 90\n  ((v &&& mask) == v) && ((((r*22+ps-weights) &&& ps) == ps)) &&\n    (((((v >>> 128)+s-(v &&& (prev*(2^90-1)))-prev) &&& s) == s))\n'
s = s.replace('def batchCheck', extra + 'def batchCheck')
s = s.replace('  ((List.range 45).all', '  extraCheck r v && ((List.range 45).all')
s = s.replace('namespace OrbitBatchFull', 'namespace OrbitBatchComplete').replace('end OrbitBatchFull', 'end OrbitBatchComplete')
open(str(STAGE / 'OrbitBatchComplete1024.lean'), 'w').write(s)
rows = [json.loads(l) for l in open(str(STAGE / 'qec360-certprobe/gb_orbit_certificates.jsonl'))]
for off in range(0, len(rows), 1024):
    g = rows[off:off + 1024]
    r = sum((1 << 128 * i for i in range(len(g))))
    v = sum((z['mask'] << 128 * i for i, z in enumerate(g)))
    q = v
    for k, m in enumerate(M):
        q = (q & r * m) + (q >> (1 << k) & r * m)
    assert q == sum((z['mask'].bit_count() << 128 * i for i, z in enumerate(g)))
    assert 22 * r + (r << 7) - q & r << 7 == r << 7
    rp = r >> 128
    sp = rp << 90
    assert (v >> 128) + sp - (v & rp * ((1 << 90) - 1)) - rp & sp == sp
print('Full274011 popcount bound and strict order batchchecks verified')
