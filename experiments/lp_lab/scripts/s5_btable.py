"""Print b-class tables (δ_b, D_b, fom cap) for larger abelian groups to locate high-cap targets."""
import sys, time
from lp_lab.fastlp import Ab
from lp_lab.design import b_table
for spec in sys.argv[1:]:
    orders = tuple(int(x) for x in spec.split(",")); G = Ab(orders); N = G.N
    t = time.time(); tab = b_table(G, delta_range=(3, 10), cap=16)
    tops = [bi for bi in tab if bi.D is not None][:6]
    print(f"G={spec} n={4*N}: {len(tab)} classes δ∈[3,10] ({time.time()-t:.0f}s)")
    for bi in tops:
        k = 2 * bi.delta
        print(f"   b={G.poly_str(bi.b):26s} δ={bi.delta:2d} D={bi.D:3d} k={k:2d} cap={bi.fom_cap:6.2f}  (d needed to beat 13.5: {((13.5*4*N/k))**0.5:.1f})")
