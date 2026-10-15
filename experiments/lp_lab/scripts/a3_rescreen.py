"""Re-screen top designs with many ISD trials (no early stop) and report orbits/toric layouts."""
import json, sys
from lp_lab.fastlp import Ab, lp_ab_rows, to_code, logical_reps
from lp_lab.distance import isd_min_weight
from lp_lab.structure import tanner_automorphisms
from lp_lab import f2
fn = sys.argv[1]; dmin = int(sys.argv[2]); trials = int(sys.argv[3]) if len(sys.argv) > 3 else 1000
recs = [json.loads(l) for l in open(fn)]
recs = [r for r in recs if r["d_ub"] >= dmin]
print(f"{fn}: {len(recs)} designs with screened d_ub ≥ {dmin}")
G = Ab(tuple(recs[0]["orders"])); n = recs[0]["n"]
out = []
for r in recs:
    A = {(0, 0): tuple(r["p"]), (0, 1): (0,), (1, 0): (0,), (1, 1): tuple(r["q"])}
    rX, rZ = lp_ab_rows(G, A, 2, tuple(r["b"]))
    LZ = logical_reps(rZ, rX, n)
    res = isd_min_weight(f2.matrix_from_rows(rZ, n), LZ, n, trials=trials, seed=11)
    code = to_code(rX, rZ, n); aut = tanner_automorphisms(code.H_X, code.H_Z)
    out.append((res.weight, aut["orbits_X"], aut["order"], r))
out.sort(key=lambda t: (-t[0], -t[1]))
for w, ox, order, r in out:
    print(f"  d_ub({trials} trials)={w:2d} (screen {r['d_ub']}) k={r['k']} orbitsX={ox} |Aut|={order:.0f}  --orders {','.join(map(str,G.orders))} --p {','.join(map(str,r['p']))} --q {','.join(map(str,r['q']))} --b {','.join(map(str,r['b']))}   b={G.poly_str(tuple(r['b']))} p={G.poly_str(tuple(r['p']))} q={G.poly_str(tuple(r['q']))}")
