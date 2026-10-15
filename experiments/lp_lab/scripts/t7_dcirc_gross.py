"""Calibrate Stim circuit-distance search on gross [[144,12,12]] (paper: d_circ ≤ 10)."""
import time, sys
from lp_lab.lp import Group, lp_square
from lp_lab.structure import two_block_perms
from lp_lab.circuit import build_memory_circuit, circuit_distance_upper_bound
from lp_lab.distance import logical_basis
from lp_lab import f2
G = Group((12, 6)); A = {(0, 0): G.parse("x^3 + y + y^2")}; B = G.parse("y^3 + x + x^2")
code = lp_square(G, A, 1, G.conj(B)); n = code.n; P = two_block_perms(G, A, 1, G.conj(B), code)
LZ = f2.matrix_from_rows(logical_basis(code.H_Z, code.H_X, n), n); LX = f2.matrix_from_rows(logical_basis(code.H_X, code.H_Z, n), n)
set_size = int(sys.argv[1]) if len(sys.argv) > 1 else 3
for basis, logs in (("Z", LZ), ("X", LX)):
    circ = build_memory_circuit(P, code.H_X, code.H_Z, logs, basis=basis, cycles=12, p=1e-3)
    print(basis, "detectors", circ.num_detectors, "DEM errors", circ.detector_error_model().num_errors, flush=True)
    t = time.time(); w, errs = circuit_distance_upper_bound(circ, max_set_size=set_size)
    print(f"{basis}: d_circ upper bound = {w} (set_size={set_size}, {time.time()-t:.0f}s)", flush=True)
