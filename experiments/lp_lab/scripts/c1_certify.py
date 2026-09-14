"""Certify a finalist LP(A,b) code: exact d (SAT ladder), Aut orbits (nauty),
wheel-layer planarity, Stim depth-7 circuit sanity + circuit-distance search.

Usage: c1_certify.py --orders 18,2 --p 0,5 --q 7,20 --b 0,2,4 [--a 2] [--dlow 8] [--no-circuit]
For a=3 all-monomial: --e e11,e12,e21,e22 instead of --p/--q.
"""
from __future__ import annotations
import argparse, json, time
from pathlib import Path
import numpy as np
from lp_lab.lp import Group, lp_square
from lp_lab.distance import code_k, distance_upper_bounds, exact_distance, logical_basis
from lp_lab.structure import tanner_automorphisms, tanner_automorphisms_swap, two_block_perms, thickness2_certificate
from lp_lab import f2


def build(args):
    orders = tuple(int(x) for x in args.orders.split(","))
    G = Group(orders)
    b = frozenset(int(x) for x in args.b.split(","))
    if args.e:
        e = [int(x) for x in args.e.split(",")]
        one = G.mono(*([0] * len(orders)))
        A = {(0, 0): one, (0, 1): one, (0, 2): one, (1, 0): one, (2, 0): one,
             (1, 1): frozenset([e[0]]), (1, 2): frozenset([e[1]]), (2, 1): frozenset([e[2]]), (2, 2): frozenset([e[3]])}
        a = 3
    elif args.r:
        p = frozenset(int(x) for x in args.p.split(",")); q = frozenset(int(x) for x in args.q.split(","))
        r = frozenset(int(x) for x in args.r.split(","))
        one = G.mono(*([0] * len(orders)))
        A = {(0, 0): p, (0, 1): one, (1, 1): q, (1, 2): one, (2, 2): r, (2, 0): one}
        a = 3
    else:
        p = frozenset(int(x) for x in args.p.split(",")); q = frozenset(int(x) for x in args.q.split(","))
        one = G.mono(*([0] * len(orders)))
        A = {(0, 0): p, (0, 1): one, (1, 0): one, (1, 1): q}
        a = 2
    return G, A, a, b


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--orders", required=True); ap.add_argument("--p"); ap.add_argument("--q"); ap.add_argument("--e"); ap.add_argument("--r", help="a=3 3-cycle pattern: A=[[p,1,0],[0,q,1],[1,0,r]]")
    ap.add_argument("--b", required=True); ap.add_argument("--dlow", type=int, default=1)
    ap.add_argument("--no-circuit", action="store_true"); ap.add_argument("--cycles", type=int, default=None)
    ap.add_argument("--circ-set-size", type=int, default=4)
    ap.add_argument("--tag", default=None); ap.add_argument("--skip-exact", action="store_true")
    args = ap.parse_args()
    G, A, a, b = build(args)
    code = lp_square(G, A, a, b); code.check_css()
    n = code.n; k = code_k(code)
    desc = {"orders": G.orders, "a": a, "A": {f"{i},{j}": G.poly_str(v) for (i, j), v in A.items()}, "b": G.poly_str(b)}
    print("code:", json.dumps(desc)); print(f"n={n} k={k} weights={code.weights()}")
    rep = {"desc": desc, "n": n, "k": k}
    t = time.time(); rX, rZ = distance_upper_bounds(code, trials=300, seed=7)
    print(f"ISD (300 trials): d_X<= {rX.weight}, d_Z<= {rZ.weight}  ({time.time()-t:.1f}s)")
    rep["isd"] = {"dX_ub": rX.weight, "dZ_ub": rZ.weight}
    # structure
    aut = tanner_automorphisms(code.H_X, code.H_Z); aut.pop("orbit_labels")
    aut_swap = tanner_automorphisms_swap(code.H_X, code.H_Z)
    print("Aut(Tanner):", aut, "| with X<->Z swap:", aut_swap)
    rep["aut"] = aut; rep["aut_swap"] = aut_swap
    rep["not_2bga"] = aut["orbits_X"] >= 2
    P = two_block_perms(G, A, a, b, code)
    cert = thickness2_certificate(P)
    good = [c for c in cert if c["planar_A"] and c["planar_B"]]
    print(f"thickness-2 wheel layers: {len(good)}/9 layer choices planar; e.g. {good[0] if good else None}")
    rep["thickness2"] = {"planar_choices": len(good), "example": good[0] if good else None}
    from lp_lab.structure import toric_layouts
    lay = toric_layouts(P)
    print(f"toric layouts (generalised Lemma 4): {len(lay)} verified; e.g. {lay[0] if lay else None}")
    rep["toric_layouts"] = {"count": len(lay), "examples": lay[:4]}
    # exact distance
    if not args.skip_exact:
        symmetric = all(A.get((i, j)) == A.get((j, i)) for i in range(a) for j in range(a))
        t = time.time()
        rx = exact_distance(code, "X", lower=args.dlow, upper=rX.weight, verbose=True, code_id="finalist")
        dX = rx.distance; print(f"exact d_X = {dX} ({time.time()-t:.1f}s)")
        if symmetric:
            dZ = dX; print("A symmetric ⇒ d_Z = d_X")
        else:
            t = time.time(); rz = exact_distance(code, "Z", lower=args.dlow, upper=rZ.weight, verbose=True, code_id="finalistZ")
            dZ = rz.distance; print(f"exact d_Z = {dZ} ({time.time()-t:.1f}s)")
        rep["d_exact"] = {"dX": dX, "dZ": dZ, "d": min(dX, dZ)}
        rep["witness_X"] = [int(i) for i in np.flatnonzero(rx.witness)]
        d = min(dX, dZ)
        print(f"*** [[{n},{k},{d}]]  kd^2/n = {k*d*d/n:.2f}")
    # circuit
    if not args.no_circuit:
        from lp_lab.circuit import build_memory_circuit, circuit_distance_upper_bound
        LZ = f2.matrix_from_rows(logical_basis(code.H_Z, code.H_X, n), n)
        LX = f2.matrix_from_rows(logical_basis(code.H_X, code.H_Z, n), n)
        dref = rep.get("d_exact", {}).get("d", min(rX.weight, rZ.weight))
        cycles = args.cycles or dref
        rep["circuit"] = {}
        for basis, logs in (("Z", LZ), ("X", LX)):
            circ0 = build_memory_circuit(P, code.H_X, code.H_Z, logs, basis=basis, cycles=2, p=0.0)
            s = circ0.compile_detector_sampler().sample(shots=8, append_observables=True)
            ok = not s.any()
            circ = build_memory_circuit(P, code.H_X, code.H_Z, logs, basis=basis, cycles=cycles, p=1e-3)
            t = time.time()
            w, errs = circuit_distance_upper_bound(circ, max_set_size=args.circ_set_size)
            print(f"circuit basis {basis}: noiseless-deterministic={ok}; cycles={cycles}; d_circ upper bound = {w} ({time.time()-t:.1f}s)")
            rep["circuit"][basis] = {"deterministic": ok, "cycles": cycles, "dcirc_ub": w, "set_size": args.circ_set_size}
    tag = args.tag or f"cert_{'x'.join(map(str,G.orders))}_a{a}"
    Path("data/certs").mkdir(parents=True, exist_ok=True)
    Path(f"data/certs/{tag}.json").write_text(json.dumps(rep, indent=1, default=str))
    print("saved", f"data/certs/{tag}.json")


if __name__ == "__main__":
    main()
