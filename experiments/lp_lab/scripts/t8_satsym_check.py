"""Validate symmetry-broken SAT ladder on the certified [[144,8,12]] (Z18×Z2 a=2) and gross."""
import time, numpy as np, sys
from lp_lab.lp import Group, lp_square
from lp_lab.distance import logical_basis
from lp_lab.satdist import x_distance_sym
from lp_lab import f2

def main():
    which = sys.argv[1] if len(sys.argv) > 1 else "144"
    if which == "144":
        G = Group((18, 2)); one = G.mono(0, 0)
        A = {(0, 0): frozenset([0, 5]), (0, 1): one, (1, 0): one, (1, 1): frozenset([7, 20])}; b = frozenset([0, 2, 4])
        code = lp_square(G, A, 2, b); N = 36; reps = [0, N, 2 * N, 3 * N]; lower = 10
    else:
        G = Group((12, 6)); A = {(0, 0): G.parse("x^3 + y + y^2")}; B = G.parse("y^3 + x + x^2")
        code = lp_square(G, A, 1, G.conj(B)); N = 72; reps = [0, N]; lower = 10
    n = code.n
    LZ = f2.matrix_from_rows(logical_basis(code.H_Z, code.H_X, n), n)
    t = time.time(); res = x_distance_sym(code.H_X, code.H_Z, LZ, reps, lower, 14, workers=4)
    print(which, "d_X =", res.distance, "unsat:", res.unsat_weights, f"total {time.time()-t:.0f}s")

if __name__ == "__main__":
    main()
