"""Depth-7 syndrome circuit for an a=2 candidate: determinism check, DEM size, d_circ upper bound."""
import argparse, json, time
from pathlib import Path
import numpy as np
from lp_lab.lp import Group, lp_square
from lp_lab.structure import two_block_perms
from lp_lab.circuit import build_memory_circuit, circuit_distance_upper_bound
from lp_lab.distance import logical_basis
from lp_lab import f2

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--orders", required=True); ap.add_argument("--p", required=True); ap.add_argument("--q", required=True); ap.add_argument("--b", required=True)
    ap.add_argument("--cycles", type=int, default=4); ap.add_argument("--set-size", type=int, default=2); ap.add_argument("--basis", default="ZX")
    ap.add_argument("--tag", required=True); ap.add_argument("--p-noise", type=float, default=1e-3)
    args = ap.parse_args()
    orders = tuple(int(x) for x in args.orders.split(",")); G = Group(orders); one = G.mono(*([0] * len(orders)))
    p = frozenset(int(x) for x in args.p.split(",")); q = frozenset(int(x) for x in args.q.split(",")); b = frozenset(int(x) for x in args.b.split(","))
    A = {(0, 0): p, (0, 1): one, (1, 0): one, (1, 1): q}
    code = lp_square(G, A, 2, b); n = code.n; P = two_block_perms(G, A, 2, b, code)
    LZ = f2.matrix_from_rows(logical_basis(code.H_Z, code.H_X, n), n); LX = f2.matrix_from_rows(logical_basis(code.H_X, code.H_Z, n), n)
    rep = {"tag": args.tag, "n": n, "k": int(LZ.shape[0]), "cycles": args.cycles, "set_size": args.set_size}
    for basis, logs in (("Z", LZ), ("X", LX)):
        if basis not in args.basis: continue
        c0 = build_memory_circuit(P, code.H_X, code.H_Z, logs, basis=basis, cycles=2, p=0.0)
        s = c0.compile_detector_sampler().sample(shots=16, append_observables=True)
        det_ok = bool(not s.any())
        circ = build_memory_circuit(P, code.H_X, code.H_Z, logs, basis=basis, cycles=args.cycles, p=args.p_noise)
        dem = circ.detector_error_model()
        print(f"{basis}: deterministic={det_ok} qubits={circ.num_qubits} depth(ticks)={circ.num_ticks} detectors={circ.num_detectors} DEM errors={dem.num_errors}", flush=True)
        t = time.time(); w, errs = circuit_distance_upper_bound(circ, max_set_size=args.set_size)
        print(f"{basis}: d_circ ≤ {w}  (cycles={args.cycles}, set_size={args.set_size}, {time.time()-t:.0f}s)", flush=True)
        rep[basis] = {"deterministic": det_ok, "detectors": circ.num_detectors, "dem_errors": dem.num_errors, "dcirc_ub": w, "seconds": time.time() - t,
                      "fault_locations": [str(e) for e in errs][:40]}
    Path("data/certs").mkdir(exist_ok=True)
    Path(f"data/certs/{args.tag}_circuit_c{args.cycles}_s{args.set_size}.json").write_text(json.dumps(rep, indent=1))

if __name__ == "__main__":
    main()
