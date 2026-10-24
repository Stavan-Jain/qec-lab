"""Emit kernel-checkable certificates into QEC360_STAGE."""
import os
from pathlib import Path
STAGE = Path(os.environ['QEC360_STAGE']).resolve()
import pathlib, json, sys
P = pathlib.Path(str(STAGE / 'qec360-certprobe'))
sys.path.insert(0, str(P))
from balance_literals import lit
from compact_literals import packlit
rs = [json.loads(l) for l in open(P / 'nonempty_peel_certificates.jsonl')]
F = json.load(open(P / 'gb_frame.json'))
E = F['Ecols']
R = F['Rcols']
cols = [sum(((x >> i & 1) << j for j, x in enumerate(E))) | sum(((x >> i & 1) << j for j, x in enumerate(R))) << 90 for i in range(45)]
D = {int(r['beta'], 16): r for l in open(P / 'cap2_pairplanes_bottom.jsonl') if (r := json.loads(l))['allowed']}

def recon(x):
    z = 0
    while x:
        l = x & -x
        z ^= cols[l.bit_length() - 1]
        x ^= l
    return z

def pack(xs):
    return sum((x << 256 * i for i, x in enumerate(xs)))

def spread(x):
    return sum(((x >> i & 1) << 128 * i for i in range(91)))
full = (1 << 91) - 1
rep = sum((1 << 128 * i for i in range(91)))
table = sum((sum(((x >> i & 1) << 128 * i for i in range(7))) << 896 * x for x in range(128)))
with open(str(STAGE / 'Bravyi360PairTables.lean'), 'w') as f:
    f.write('import CapTwoFast\nset_option maxRecDepth 100000\nset_option maxHeartbeats 0\nnamespace Quantum.Stabilizer.Homological.BB.Bravyi360\nopen PackedDigits\n')
    f.write('/-- Shared lookup for spreading seven bits into 128-bit lanes. -/\ndef bottomPairSpreadTable : Nat := ' + lit(table) + '\n')
    f.write('/-- Unit lanes for the 91 sentinel-pair coordinates. -/\ndef bottomPairRep : Nat := ((1 <<< (128*91))-1)/(2^128-1)\n')
    f.write('/-- Every small lookup entry has its exact spreading semantics. -/\ntheorem bottomPair_lookup : ∀ x : Fin 128, dilationLookup (128*7) bottomPairSpreadTable x = dilate 128 7 x := by decide +kernel\n')
    f.write('/-- The repeated-lane constant has its packed interpretation. -/\ntheorem bottomPair_rep : bottomPairRep = packF 128 91 (fun _ => 1) := by decide +kernel\n')
    f.write('end Quantum.Stabilizer.Homological.BB.Bravyi360\n')
count = 0
numcodes = 0
for part, start in enumerate(range(0, len(rs), 2000)):
    f = open(P / f'BottomPositiveTwoPart{part}.lean', 'w')
    f.write(f'import BottomCap2Cover\nimport PackedDataExpand\nimport Bravyi360PairTables\nimport BottomForestFastPart{part}\nset_option maxHeartbeats 0\nset_option maxRecDepth 100000\nnamespace Quantum.Stabilizer.Homological.BB.Bravyi360.BottomForest\n')
    for i in range(start, min(start + 2000, len(rs)), 250):
        covers = []
        for r in rs[i:min(i + 250, len(rs))]:
            b = int(r['beta'], 16)
            w = b.bit_count()
            ds = []
            codes = []
            target = 0
            if w == 18:
                count += 1
                d = D[b]
                ds = [(int(v, 16), recon(int(v, 16))) for v in d['duals']]
                target = d['target']
                allowed = {int(o, 16) for o, s in r['seeds']}
                assert allowed == {int(x, 16) for x in d['allowed']}
                for o in sorted(allowed):
                    qs = [q for q in range(90) if o >> q & 1]
                    assert len(qs) <= 2
                    if not qs:
                        codes += [128 * q + q for q in range(91) if not b >> q & 1]
                    else:
                        q, t = (qs + [90])[:2]
                        codes += [128 * q + t, 128 * t + q]
                codes = sorted(set(codes))
                numcodes += len(codes)
                got = spread(full ^ b) * full & (full ^ b) * rep
                for j, (c, z) in enumerate(ds):
                    p = z & (1 << 90) - 1
                    assert z & b == 0 and (z >> 90 & b).bit_count() % 2 == target >> j & 1
                    v = got & (spread(p) * full ^ p * rep)
                    got = v if target >> j & 1 else got ^ v
                assert got == sum((1 << c for c in codes))
                assert all(((1 << c // 128 ^ 1 << c % 128) & (1 << 90) - 1 in allowed for c in codes))
            covers.append((ds, target, codes))
        terms = ['⟨[' + ','.join(('(' + hex(c) + ',' + hex(z) + ')' for c, z in ds)) + '],' + str(target) + ',[' + ','.join(map(str, codes)) + ']⟩' for ds, target, codes in covers]
        pairs = [p for ds, _, _ in covers for p in ds]
        vals = [f'((1 <<< (256*{len(pairs)}))-1)/(2^256-1)', packlit([c for c, z in pairs], 256), packlit([z for c, z in pairs], 256)]
        f.write(f'/-- Cap-two outside censuses aligned with forest chunk {i}. -/\ndef twoCovers{i} : List TwoCover := [' + ','.join(terms) + ']\n')
        f.write(f'/-- Cached dual image packings for cap-two chunk {i}. -/\ndef twoPacked{i} : DualPacking := ⟨' + ','.join((x for x in vals)) + '⟩\n')
        f.write(f'/-- Every cap-two outside census and actual dual image is kernel-checked. -/\ntheorem twoChecked{i} : twoGroupCheck bottomPairSpreadTable bottomPairRep records{i} twoCovers{i} twoPacked{i} = true := by decide +kernel\n')
    f.write(f'#print axioms twoChecked{start}\nend Quantum.Stabilizer.Homological.BB.Bravyi360.BottomForest\n')
    f.close()
print('positive cap2', count, 'codes', numcodes)
