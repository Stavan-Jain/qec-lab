"""Pick top designs from a design JSONL: max d_ub, then nauty orbits; print non-2BGA ones first."""
import json, sys
from lp_lab.fastlp import Ab, lp_ab_rows, to_code
from lp_lab.structure import tanner_automorphisms
fn = sys.argv[1]; top = int(sys.argv[2]) if len(sys.argv) > 2 else 6
recs = [json.loads(l) for l in open(fn)]
if not recs: print("empty"); sys.exit()
dmax = max(r["d_ub"] for r in recs)
cands = [r for r in recs if r["d_ub"] == dmax]
print(f"{fn}: {len(recs)} designs, max d_ub={dmax} attained by {len(cands)}")
G = Ab(tuple(recs[0]["orders"])); n = recs[0]["n"]
rows = []
for r in cands[:400]:
    A = {(0, 0): tuple(r["p"]), (0, 1): (0,), (1, 0): (0,), (1, 1): tuple(r["q"])}
    rX, rZ = lp_ab_rows(G, A, 2, tuple(r["b"])); code = to_code(rX, rZ, n)
    aut = tanner_automorphisms(code.H_X, code.H_Z)
    rows.append((aut["orbits_X"], aut["order"], r))
rows.sort(key=lambda t: (-t[0], t[1]))
from collections import Counter
print("X-orbit histogram:", Counter(o for o, _, _ in rows))
for o, order, r in rows[:top]:
    print(f"  orbitsX={o} |Aut|={order:.0f} [[{n},{r['k']},≤{r['d_ub']}]] b={G.poly_str(tuple(r['b']))} | p={G.poly_str(tuple(r['p']))} | q={G.poly_str(tuple(r['q']))}  args: --orders {','.join(map(str,G.orders))} --p {','.join(map(str,r['p']))} --q {','.join(map(str,r['q']))} --b {','.join(map(str,r['b']))}")
