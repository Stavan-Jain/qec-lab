"""Deep re-screen (many ISD trials) + orbits + toric layouts for designs with screened d_ub ≥ dmin."""
import json, sys
from lp_lab.lp import Group, lp_square
from lp_lab.fastlp import Ab, lp_ab_rows, logical_reps
from lp_lab.distance import isd_min_weight
from lp_lab.structure import tanner_automorphisms, two_block_perms, toric_layouts
from lp_lab import f2
files = sys.argv[1:-2]; dmin = int(sys.argv[-2]); trials = int(sys.argv[-1])
rows = []
for fn in files:
    recs = [json.loads(l) for l in open(fn) if l.strip()]
    recs = [r for r in recs if r["d_ub"] >= dmin]
    if not recs: continue
    orders = tuple(recs[0]["orders"]); G = Ab(orders); Gg = Group(orders); one = Gg.mono(*([0] * len(orders))); n = recs[0]["n"]
    for r in recs:
        A = {(0, 0): tuple(r["p"]), (0, 1): (0,), (1, 0): (0,), (1, 1): tuple(r["q"])}
        rX, rZ = lp_ab_rows(G, A, 2, tuple(r["b"])); LZ = logical_reps(rZ, rX, n)
        res = isd_min_weight(f2.matrix_from_rows(rZ, n), LZ, n, trials=trials, seed=17)
        Ag = {(0, 0): frozenset(r["p"]), (0, 1): one, (1, 0): one, (1, 1): frozenset(r["q"])}; b = frozenset(r["b"])
        code = lp_square(Gg, Ag, 2, b); P = two_block_perms(Gg, Ag, 2, b, code)
        lay = toric_layouts(P); aut = tanner_automorphisms(code.H_X, code.H_Z)
        rows.append((res.weight, len(lay), aut["orbits_X"], aut["order"], orders, r, lay[0]["torus"] if lay else None))
        print(f"  d_ub({trials})={res.weight:2d} screen={r['d_ub']:2d} k={r['k']} toric={len(lay)}{(' '+str(lay[0]['torus'])) if lay else ''} orbitsX={aut['orbits_X']} |Aut|={aut['order']:.0f} --orders {','.join(map(str,orders))} --p {','.join(map(str,r['p']))} --q {','.join(map(str,r['q']))} --b {','.join(map(str,r['b']))}  b={Gg.poly_str(b)} p={Gg.poly_str(Ag[(0,0)])} q={Gg.poly_str(Ag[(1,1)])}", flush=True)
rows.sort(key=lambda t: (-t[0], -t[1], -t[2]))
print("=== ranked ===")
for w, nl, ox, order, orders, r, torus in rows[:25]:
    print(f"  d_ub={w} toric={nl} {torus} orbitsX={ox} |Aut|={order:.0f} --orders {','.join(map(str,orders))} --p {','.join(map(str,r['p']))} --q {','.join(map(str,r['q']))} --b {','.join(map(str,r['b']))}")
