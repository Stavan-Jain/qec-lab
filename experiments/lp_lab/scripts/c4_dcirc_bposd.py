"""d_circ upper bounds via BP-OSD (paper's method) for BB72 (validation) or an a=2 candidate."""
import argparse, json, time
from pathlib import Path
from lp_lab.lp import Group, lp_square
from lp_lab.structure import two_block_perms
from lp_lab.circuit import build_memory_circuit
from lp_lab.circdist import circuit_distance_bposd
from lp_lab.distance import logical_basis
from lp_lab import f2

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--code", default="bb72"); ap.add_argument("--orders"); ap.add_argument("--p"); ap.add_argument("--q"); ap.add_argument("--b")
    ap.add_argument("--cycles", type=int, default=3); ap.add_argument("--trials", type=int, default=60); ap.add_argument("--tag", default=None)
    args = ap.parse_args()
    if args.code == "bb72":
        G = Group((6, 6)); A = {(0, 0): G.parse("x^3 + y + y^2")}; B = G.parse("y^3 + x + x^2"); b = G.conj(B); a = 1
    elif args.code == "gross":
        G = Group((12, 6)); A = {(0, 0): G.parse("x^3 + y + y^2")}; B = G.parse("y^3 + x + x^2"); b = G.conj(B); a = 1
    else:
        orders = tuple(int(x) for x in args.orders.split(",")); G = Group(orders); one = G.mono(*([0] * len(orders)))
        p = frozenset(int(x) for x in args.p.split(",")); q = frozenset(int(x) for x in args.q.split(",")); b = frozenset(int(x) for x in args.b.split(","))
        A = {(0, 0): p, (0, 1): one, (1, 0): one, (1, 1): q}; a = 2
    code = lp_square(G, A, a, b); n = code.n; P = two_block_perms(G, A, a, b, code)
    LZ = f2.matrix_from_rows(logical_basis(code.H_Z, code.H_X, n), n); LX = f2.matrix_from_rows(logical_basis(code.H_X, code.H_Z, n), n)
    rep = {"code": args.code, "n": n, "k": int(LZ.shape[0]), "cycles": args.cycles, "trials": args.trials}
    for basis, logs in (("Z", LZ), ("X", LX)):
        circ = build_memory_circuit(P, code.H_X, code.H_Z, logs, basis=basis, cycles=args.cycles, p=1e-3)
        t = time.time(); r = circuit_distance_bposd(circ, trials=args.trials, seed=1, verbose=True)
        print(f"{args.code} basis {basis}: d_circ ≤ {r['dcirc_ub']} ({r['mechanisms']} mechanisms, {r['detectors']} detectors, {time.time()-t:.0f}s)", flush=True)
        rep[basis] = {k: v for k, v in r.items() if k != "witness_mechanisms"}
        rep[basis]["witness"] = [str(m) for m in (r["witness_mechanisms"] or [])]
    Path("data/certs").mkdir(exist_ok=True)
    Path(f"data/certs/{args.tag or args.code}_dcirc_bposd.json").write_text(json.dumps(rep, indent=1))

if __name__ == "__main__":
    main()
