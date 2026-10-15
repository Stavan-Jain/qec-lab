import time, itertools, random
import numpy as np
from lp_lab.lp import Group, lp_square
from lp_lab.distance import code_k, distance_upper_bounds
from lp_lab.fastlp2 import hx_rows, hz_rows, k_of, to_numpy
from lp_lab import f2

# consistency: fast builder vs generic builder on random instances
rng = random.Random(0)
for trial in range(5):
    ell = rng.choice([9, 12, 15, 18])
    p = (0, rng.randrange(1, ell)); q = tuple(rng.sample(range(ell), 2)); b = (0,) + tuple(rng.sample(range(1, ell), 2))
    G = Group((ell,))
    A = {(0,0): frozenset(p), (0,1): G.mono(0), (1,0): G.mono(0), (1,1): frozenset(q)}
    code = lp_square(G, A, 2, frozenset(b))
    fast = to_numpy(ell, p, q, b)
    assert (code.H_X == fast.H_X).all(), "H_X mismatch"
    assert (code.H_Z == fast.H_Z).all(), "H_Z mismatch"
    fast.check_css()
    assert code_k(code) == k_of(ell, p, q, b)
print("fast builder agrees with generic LP construction")
t = time.time(); cnt = 0
for _ in range(300):
    ell = 36
    p = (0, rng.randrange(1, ell)); q = tuple(rng.sample(range(ell), 2)); b = (0,) + tuple(rng.sample(range(1, ell), 2))
    k_of(ell, p, q, b); cnt += 1
print(f"k throughput ell=36: {cnt/(time.time()-t):.0f} codes/s")
# ISD timing on gross-size
G = Group((12, 6)); A = {(0,0): G.parse("x^3 + y + y^2")}; B = G.parse("y^3 + x + x^2")
code = lp_square(G, A, 1, G.conj(B))
t = time.time(); rX, rZ = distance_upper_bounds(code, trials=60, seed=1); print("gross ISD", rX.weight, rZ.weight, f"{time.time()-t:.2f}s trials={rX.trials}")
