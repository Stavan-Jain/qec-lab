"""Emit kernel-checkable certificates into QEC360_STAGE."""
import os
from pathlib import Path
STAGE=Path(os.environ["QEC360_STAGE"]).resolve()
import pathlib, json, sys
sys.path.insert(0, str(STAGE / 'qec360-certprobe'))
from balance_literals import lit
P = pathlib.Path(str(STAGE / 'qec360-certprobe'))
rows = [json.loads(l) for l in open(P / 'cap2_pairplanes_top.jsonl')]
w = 256
count = 182
full = (1 << count) - 1
rep = sum((1 << w * i for i in range(count)))
stride = w * 7
table = sum((sum(((x >> i & 1) << w * i for i in range(7))) << stride * x for x in range(128)))

def spread(x):
    return sum(((x >> i & 1) << w * i for i in range(count)))
f = open(str(STAGE / 'Bravyi360TopEmpty2Data.lean'), 'w')
f.write('import CapTwoFast\nimport DualPlaneCertificate\nimport Bravyi360TopEmptyData\nset_option maxHeartbeats 0\nset_option maxRecDepth 100000\nnamespace Quantum.Stabilizer.Homological.BB.Bravyi360.TopEmpty\nopen F2Mask PackedDigits FiberData\n')
f.write('/-- Shared seven-bit spread lookup used by the cap-two certificate. -/\ndef cap2SpreadTable : Nat := ' + lit(table) + '\n')
f.write('/-- Repeated unit lanes for all padded cap-two coordinates. -/\ndef cap2Rep : Nat := ((1 <<< (256*182))-1)/(2^256-1)\n')
f.write('/-- The finite lookup table spreads every seven-bit input exactly. -/\ntheorem cap2_lookup : ∀ x : Fin 128, dilationLookup (256*7) cap2SpreadTable x = dilate 256 7 x := by decide +kernel\n')
f.write('/-- The repeated lane expression has its semantic packed form. -/\ntheorem cap2_rep : cap2Rep = packF 256 182 (fun _ => 1) := by decide +kernel\n')
for i, r in enumerate(rows):
    b = int(r['beta'], 16)
    ps = [int(x, 16) for x in r['planes']]
    cs = [int(x, 16) for x in r['duals']]
    t = r['target']
    d = spread(full ^ b) * full & (full ^ b) * rep
    for j, p in enumerate(ps):
        z = d & (spread(p) * full ^ p * rep)
        d = z if t >> j & 1 else d ^ z
    assert d == 0
    f.write(f'/-- Compressed dual coefficients and parity planes for empty cap-two fiber {i}. -/\ndef planes2_{i} : List (Nat × Nat) := [' + ','.join((f'({hex(c)},{hex(p)})' for c, p in zip(cs, ps))) + ']\n')
    f.write(f'/-- The cap-two duals are valid for the actual top matrices in fiber {i}. -/\ntheorem duals2_{i} : dualPlaneCertificate 180 8 topDualColumns (beta2 {i}) planes2_{i} {t} := by decide +kernel\n')
    f.write(f'/-- Every pair of outside coordinates fails some certified parity in fiber {i}. -/\ntheorem empty2_{i} : filterPlanes ((planes2_{i}.map Prod.snd).map (fastPairPlane 256 7 26 cap2SpreadTable cap2Rep)) {t} (fastPairDomain 256 7 26 cap2SpreadTable cap2Rep (beta2 {i})) = 0 := by decide +kernel\n')
f.write('end Quantum.Stabilizer.Homological.BB.Bravyi360.TopEmpty\n')
f.close()
p = pathlib.Path(str(STAGE / 'Bravyi360TopEmpty2.lean'))
s = p.read_text()
s = s[:s.index('/-- Every listed top cap-two fiber')]
s += '/-- Every listed top cap-two fiber is impossible with at most two outside coordinates. -/\ntheorem cap2_impossible (i : Fin 143) (word : Nat) (hw : word < 2^180)\n    (he : reconX topEColumns word = reconX topRColumns (beta2 i.val))\n    (hoff : bitCount 180 (outsideMask word (beta2 i.val)) ≤ 2) : False := by\n  fin_cases i\n'
for i, r in enumerate(rows):
    s += f"  · exact cap2_impossible_of_data (beta2 {i}) {r['target']} word planes2_{i} (by decide +kernel) hw duals2_{i} empty2_{i} he hoff\n"
s += '\n#print axioms cap2_impossible\nend Quantum.Stabilizer.Homological.BB.Bravyi360.TopEmpty\n'
p.write_text(s)
print('generated', len(rows), 'top cap-two records')
