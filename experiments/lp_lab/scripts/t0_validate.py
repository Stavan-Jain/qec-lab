"""Validate lp_lab against the gross code and a 2x2 example."""
import time, numpy as np
from lp_lab.lp import Group, lp_square, monomial_decomposition
from lp_lab.distance import code_k, distance_upper_bounds

# gross code as LP(A, b) with a = 1 over Z12 x Z6: A = x^3 + y + y^2, B = y^3 + x + x^2
G = Group((12, 6))
A = {(0, 0): G.parse("x^3 + y + y^2")}
B = G.parse("y^3 + x + x^2")
code = lp_square(G, A, 1, G.conj(B))   # H_X = [A | b* I] so pass b = B* to get [A|B]
code.check_css()
print("gross:", code.n, "k =", code_k(code), code.weights())
t = time.time(); rX, rZ = distance_upper_bounds(code, trials=60, seed=1); 
print("gross ISD d_X<=", rX.weight, "d_Z<=", rZ.weight, f"({time.time()-t:.1f}s)")

# 2x2 example over Z_36: A = [[p,1],[1,q]]
G = Group((36,))
p = G.parse("1 + x^7"); q = G.parse("x^3 + x^20"); b = G.parse("1 + x^5 + x^12")
A = {(0,0): p, (0,1): G.mono(0), (1,0): G.mono(0), (1,1): q}
code = lp_square(G, A, 2, b)
code.check_css()
print("2x2:", code.n, "k =", code_k(code), code.weights())
print("decomp:", [{k: sorted(v) for k, v in P.items()} for P in monomial_decomposition(G, A, 2)])
