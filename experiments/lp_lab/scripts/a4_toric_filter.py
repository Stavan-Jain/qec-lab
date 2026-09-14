"""Among designs attaining the max screened d_ub, list those with a toric layout (generalised Lemma 4)."""
import json, sys
from lp_lab.lp import Group, lp_square
from lp_lab.structure import two_block_perms, toric_layouts, tanner_automorphisms
fn = sys.argv[1]; dmin = int(sys.argv[2])
recs = [json.loads(l) for l in open(fn)]
recs = [r for r in recs if r["d_ub"] >= dmin]
orders = tuple(recs[0]["orders"]); G = Group(orders); one = G.mono(*([0] * len(orders)))
print(f"{fn}: {len(recs)} designs with d_ub ≥ {dmin}")
hits = 0
for r in recs:
    A = {(0, 0): frozenset(r["p"]), (0, 1): one, (1, 0): one, (1, 1): frozenset(r["q"])}; b = frozenset(r["b"])
    code = lp_square(G, A, 2, b); P = two_block_perms(G, A, 2, b, code)
    lay = toric_layouts(P)
    if lay:
        hits += 1
        aut = tanner_automorphisms(code.H_X, code.H_Z)
        print(f"  TORIC {lay[0]['torus']} ({len(lay)} layouts) d_ub={r['d_ub']} k={r['k']} orbitsX={aut['orbits_X']} --orders {','.join(map(str,orders))} --p {','.join(map(str,r['p']))} --q {','.join(map(str,r['q']))} --b {','.join(map(str,r['b']))}  b={G.poly_str(b)} p={G.poly_str(A[(0,0)])} q={G.poly_str(A[(1,1)])}")
print(f"{hits}/{len(recs)} have toric layouts")
