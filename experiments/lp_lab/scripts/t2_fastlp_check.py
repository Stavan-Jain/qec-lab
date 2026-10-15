import random
from lp_lab.lp import Group, lifted_product, lp_square
from lp_lab.distance import code_k
from lp_lab.fastlp import Ab, lp_ab_rows, lp_AB_rows, k_from_rows, to_code
rng = random.Random(3)
for trial in range(6):
    orders = rng.choice([(9,), (6, 3), (4, 5), (12,)])
    G = Group(orders); Gf = Ab(orders); N = Gf.N
    a = 2
    A = {}
    for i in range(a):
        for j in range(a):
            w = 2 if i == j else 1
            A[(i, j)] = frozenset(rng.sample(range(N), w))
    b = frozenset(rng.sample(range(N), 3))
    ref = lp_square(G, A, a, b); ref.check_css()
    rX, rZ = lp_ab_rows(Gf, {k: tuple(v) for k, v in A.items()}, a, tuple(b))
    fast = to_code(rX, rZ, ref.n)
    assert (ref.H_X == fast.H_X).all() and (ref.H_Z == fast.H_Z).all(), "lp_ab mismatch"
    assert code_k(ref) == k_from_rows(rX, rZ, ref.n)
    B = {}
    for i in range(a):
        for j in range(a):
            w = 2 if i == j else 1
            B[(i, j)] = frozenset(rng.sample(range(N), w))
    ref2 = lifted_product(G, A, (a, a), B, (a, a)); ref2.check_css()
    rX2, rZ2 = lp_AB_rows(Gf, {k: tuple(v) for k, v in A.items()}, {k: tuple(v) for k, v in B.items()}, a)
    fast2 = to_code(rX2, rZ2, ref2.n)
    assert (ref2.H_X == fast2.H_X).all() and (ref2.H_Z == fast2.H_Z).all(), "lp_AB mismatch"
    print(orders, "ok: n_ab =", ref.n, "k =", code_k(ref), "| n_AB =", ref2.n, "k =", code_k(ref2), fast2.weights())
print("fastlp agrees with generic construction")
