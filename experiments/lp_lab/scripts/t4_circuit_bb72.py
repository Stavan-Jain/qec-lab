"""Sanity: the depth-7 schedule on BB [[72,12,6]]: noiseless detectors all zero; circuit distance search."""
import time, numpy as np, stim
from lp_lab.lp import Group, lp_square
from lp_lab.structure import two_block_perms
from lp_lab.circuit import build_memory_circuit, circuit_distance_upper_bound
from lp_lab.distance import logical_basis
from lp_lab import f2

G = Group((6, 6))
A = {(0, 0): G.parse("x^3 + y + y^2")}; B = G.parse("y^3 + x + x^2")
code = lp_square(G, A, 1, G.conj(B)); n = code.n
P = two_block_perms(G, A, 1, G.conj(B), code)
LZ = f2.matrix_from_rows(logical_basis(code.H_Z, code.H_X, n), n)   # Z-logicals: ker H_X / rs H_Z
LX = f2.matrix_from_rows(logical_basis(code.H_X, code.H_Z, n), n)
print("k =", LZ.shape[0], LX.shape[0])
for basis, logs in (("Z", LZ), ("X", LX)):
    circ = build_memory_circuit(P, code.H_X, code.H_Z, logs, basis=basis, cycles=3, p=0.0)
    # noiseless: sample detectors -> must all be zero
    s = circ.compile_detector_sampler().sample(shots=20, append_observables=True)
    print(basis, "noiseless detectors all zero:", not s[:, :circ.num_detectors].any(), " observables zero:", not s[:, circ.num_detectors:].any())
    noisy = build_memory_circuit(P, code.H_X, code.H_Z, logs, basis=basis, cycles=6, p=1e-3)
    dem = noisy.detector_error_model(decompose_errors=False)
    print(basis, "noisy circuit: detectors", noisy.num_detectors, "errors in DEM", dem.num_errors)
    t = time.time()
    try:
        w, errs = circuit_distance_upper_bound(noisy, max_set_size=4)
        print(basis, "circuit distance upper bound:", w, f"({time.time()-t:.1f}s)")
    except Exception as e:
        print("search failed:", e)
