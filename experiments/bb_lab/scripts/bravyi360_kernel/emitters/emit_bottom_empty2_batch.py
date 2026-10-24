"""Emit bottom cap-two checks with batched actual-matrix dual images."""
import os
from pathlib import Path
STAGE=Path(os.environ["QEC360_STAGE"]).resolve()
import pathlib, json, sys
P = pathlib.Path(str(STAGE / 'qec360-certprobe'))
sys.path.insert(0, str(P))
from compact_literals import packlit
rows = [json.loads(l) for l in open(P / 'cap2_pairplanes_bottom.jsonl')]
rows = [r for r in rows if not r['allowed']]
assert len(rows) == 859
F = json.load(open(P / 'gb_frame.json'))
E = F['Ecols']
R = F['Rcols']
cols = [sum(((x >> i & 1) << j for j, x in enumerate(E))) | sum(((x >> i & 1) << j for j, x in enumerate(R))) << 90 for i in range(45)]

def recon(c):
    z = 0
    while c:
        t = c & -c
        z ^= cols[t.bit_length() - 1]
        c ^= t
    return z
ns = 'Quantum.Stabilizer.Homological.BB.Bravyi360.BottomEmpty'
for part, start in enumerate(range(0, len(rows), 100)):
    ds = []
    recs = []
    for i, r in enumerate(rows[start:start + 100], start):
        rr = [(int(x, 16), recon(int(x, 16))) for x in r['duals']]
        assert [z & (1 << 90) - 1 for c, z in rr] == [int(x, 16) for x in r['planes']]
        ds.extend(rr)
        recs.append(f'⟨beta2 {i},[' + ','.join((f'({hex(c)},{hex(z)})' for c, z in rr)) + f"],{r['target']}⟩")
    s = 'import BottomEmptyTwoBatch\nimport Bravyi360PairTables\nimport Bravyi360BottomEmptyData\nimport PackedDataExpand\nset_option maxHeartbeats 0\nset_option maxRecDepth 100000\nnamespace ' + ns + '\nopen F2Mask FiberData PackedDigits BottomForest\n'
    s += f'/-- Actual-matrix cap-two records in bounded batch {part}. -/\ndef batchTwoRecords{part} : List TwoRecord := [' + ','.join(recs) + ']\n'
    s += f'/-- Packed dual coefficients and full matrix images for batch {part}. -/\ndef batchTwoPacking{part} : DualPacking := ⟨((1 <<< (256*{len(ds)}))-1)/(2^256-1),' + packlit([c for c, z in ds], 256) + ',' + packlit([z for c, z in ds], 256) + '⟩\n'
    s += f'/-- All matrix images and all empty pair censuses in this batch pass the kernel. -/\ntheorem batchTwoChecked{part} : twoBatchCheck bottomPairSpreadTable bottomPairRep batchTwoRecords{part} batchTwoPacking{part} = true := by decide +kernel\n'
    s += f'/-- The stored batch has its stated cardinality. -/\ntheorem batchTwoLength{part} : batchTwoRecords{part}.length = {len(recs)} := by decide +kernel\n'
    s += '#print axioms batchTwoChecked' + str(part) + '\nend ' + ns + '\n'
    pathlib.Path(f'{STAGE}/Bravyi360BottomEmpty2BatchData{part}.lean').write_text(s)
s = ''.join((f'import Bravyi360BottomEmpty2BatchData{k}\n' for k in range(9))) + 'import Mathlib.Tactic.FinCases\nset_option maxRecDepth 100000\nset_option maxHeartbeats 0\nnamespace ' + ns + '\nopen F2Mask FiberData PackedDigits BottomForest\n'
s += '/-- Select one bounded batch of actual cap-two records. -/\ndef batchTwoRecordsAt : Nat → List TwoRecord\n' + ''.join((f'  | {k} => batchTwoRecords{k}\n' for k in range(9))) + '  | _ => []\n'
s += '/-- Select a cap-two record at its canonical empty-fiber index. -/\ndef batchTwoRecordAt (i : Nat) : TwoRecord := (batchTwoRecordsAt (i/100)).getD (i%100) default\n'
s += '/-- Batch records use exactly the shared empty-fiber shadow order. -/\ntheorem batchTwoBeta : ∀ i : Fin 859, (batchTwoRecordAt i.val).beta = beta2 i.val := by\n  decide +kernel\n/-- Every valid batch has the expected number of records. -/\ntheorem batchTwoLength (b : Fin 9) :\n    (batchTwoRecordsAt b.val).length = if b.val < 8 then 100 else 59 := by\n  fin_cases b\n'
for k in range(9):
    s += f'  · exact batchTwoLength{k}\n'
s += '/-- Every canonical empty cap-two index denotes a member of its bounded batch. -/\ntheorem batchTwoRecordAt_mem (i : Fin 859) :\n    batchTwoRecordAt i.val ∈ batchTwoRecordsAt (i.val/100) := by\n  have hb : i.val/100 < 9 := by omega\n  have hi : i.val%100 < (batchTwoRecordsAt (i.val/100)).length := by\n    rw [batchTwoLength ⟨i.val/100,hb⟩]\n    change i.val%100 < if i.val/100 < 8 then 100 else 59\n    have hilim := i.isLt\n    have hdecomp := Nat.mod_add_div i.val 100\n    split_ifs <;> omega\n  unfold batchTwoRecordAt\n  rw [List.getD_eq_getElem _ _ hi]\n  exact List.getElem_mem hi\n/-- Every listed batch excludes all solutions with outside weight at most two. -/\ntheorem batchTwoRecords_impossible (b : Fin 9) (r : TwoRecord)\n    (hr : r ∈ batchTwoRecordsAt b.val) (word : Nat) (hw : word < 2^90)\n    (he : reconX bottomEColumns word = reconX bottomRColumns r.beta)\n    (hoff : bitCount 90 (outsideMask word r.beta) ≤ 2) : False := by\n  fin_cases b\n'
for k in range(9):
    s += f'  · exact twoBatch_impossible bottomPairSpreadTable bottomPairRep bottomPair_lookup bottomPair_rep\n      batchTwoRecords{k} batchTwoPacking{k} batchTwoChecked{k} r hr word hw he hoff\n'
s += '/-- Every listed bottom cap-two fiber is impossible with at most two outside coordinates. -/\ntheorem cap2_impossible (i : Fin 859) (word : Nat) (hw : word < 2^90)\n    (he : reconX bottomEColumns word = reconX bottomRColumns (beta2 i.val))\n    (hoff : bitCount 90 (outsideMask word (beta2 i.val)) ≤ 2) : False :=\n  batchTwoRecords_impossible ⟨i.val/100,by omega⟩ (batchTwoRecordAt i.val)\n    (batchTwoRecordAt_mem i) word hw\n    (by simpa only [batchTwoBeta i] using he) (by simpa only [batchTwoBeta i] using hoff)\n#print axioms cap2_impossible\n' + f'end {ns}\n'
pathlib.Path(str(STAGE / 'Bravyi360BottomEmpty2Batch.lean')).write_text(s)
print('Generated859 batched cap2 records, dual reconstructions', sum((len(r['duals']) for r in rows)))
pathlib.Path(str(STAGE / 'Bravyi360BottomEmpty2.lean')).write_text('import Bravyi360BottomEmpty2Batch\n')
