"""Validate SAT circuit distance on BB [[72,12,6]] (Stim search gave d_circ ≤ 6 both bases)."""
import time, sys
from lp_lab.lp import Group, lp_square
from lp_lab.structure import two_block_perms
from lp_lab.circuit import build_memory_circuit
from lp_lab.circdist import circuit_distance_sat
from lp_lab.distance import logical_basis
from lp_lab import f2
G = Group((6, 6)); A = {(0, 0): G.parse("x^3 + y + y^2")}; B = G.parse("y^3 + x + x^2")
code = lp_square(G, A, 1, G.conj(B)); n = code.n; P = two_block_perms(G, A, 1, G.conj(B), code)
LZ = f2.matrix_from_rows(logical_basis(code.H_Z, code.H_X, n), n); LX = f2.matrix_from_rows(logical_basis(code.H_X, code.H_Z, n), n)
cycles = int(sys.argv[1]) if len(sys.argv) > 1 else 3
for basis, logs in (("Z", LZ), ("X", LX)):
    circ = build_memory_circuit(P, code.H_X, code.H_Z, logs, basis=basis, cycles=cycles, p=1e-3)
    t = time.time(); res = circuit_distance_sat(circ, lower=3, upper=8)
    print(f"BB72 basis {basis}, cycles={cycles}: d_circ = {res.distance} (unsat {res.unsat_weights}) {time.time()-t:.0f}s", flush=True)
