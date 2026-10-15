"""Exact d_X for an a=2 LP([[p,1],[1,q]], b) code via the symmetry-broken CMS ladder.
Writes data/certs/<tag>_sat.json with witness and certified UNSAT rungs."""
import argparse, json, time
from pathlib import Path
import numpy as np
from lp_lab.lp import Group, lp_square
from lp_lab.distance import logical_basis, distance_upper_bounds
from lp_lab.satdist import x_distance_sym
from lp_lab import f2

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--orders", required=True); ap.add_argument("--p", required=True); ap.add_argument("--q", required=True)
    ap.add_argument("--b", required=True); ap.add_argument("--lower", type=int, default=2); ap.add_argument("--upper", type=int, default=None)
    ap.add_argument("--workers", type=int, default=4); ap.add_argument("--tag", required=True)
    args = ap.parse_args()
    orders = tuple(int(x) for x in args.orders.split(",")); G = Group(orders); N = G.N
    one = G.mono(*([0] * len(orders)))
    p = frozenset(int(x) for x in args.p.split(",")); q = frozenset(int(x) for x in args.q.split(",")); b = frozenset(int(x) for x in args.b.split(","))
    A = {(0, 0): p, (0, 1): one, (1, 0): one, (1, 1): q}
    code = lp_square(G, A, 2, b); code.check_css(); n = code.n
    LZ = f2.matrix_from_rows(logical_basis(code.H_Z, code.H_X, n), n); k = LZ.shape[0]
    upper = args.upper
    if upper is None:
        rX, _ = distance_upper_bounds(code, trials=200, seed=3); upper = rX.weight
    print(f"[[{n},{k},?]] G={orders} p={G.poly_str(p)} q={G.poly_str(q)} b={G.poly_str(b)}; ladder {args.lower}..{upper}", flush=True)
    reps = [0, N, 2 * N, 3 * N]   # one qubit per translation orbit (L0, L1, R0, R1)
    t = time.time()
    res = x_distance_sym(code.H_X, code.H_Z, LZ, reps, args.lower, upper, workers=args.workers)
    out = {"n": n, "k": k, "orders": orders, "p": sorted(p), "q": sorted(q), "b": sorted(b), "d_X": res.distance,
           "A_symmetric_so_dZ_equals_dX": True, "unsat_weights": res.unsat_weights, "log": res.log,
           "witness_support": [int(i) for i in np.flatnonzero(res.witness)], "seconds": time.time() - t,
           "method": "CryptoMiniSat ladder, native XOR rows, seqcounter cardinality, translation-orbit symmetry breaking (4 reps)"}
    Path("data/certs").mkdir(parents=True, exist_ok=True)
    Path(f"data/certs/{args.tag}_sat.json").write_text(json.dumps(out, indent=1))
    print(f"*** [[{n},{k},{res.distance}]]  kd^2/n = {k*res.distance**2/n:.2f}   ({time.time()-t:.0f}s)")

if __name__ == "__main__":
    main()
