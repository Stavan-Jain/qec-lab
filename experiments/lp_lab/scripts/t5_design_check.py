"""Verify k = 2·dim R/(b, 1+pq) numerically (incl. even & bivariate) and that max_k_designs hits k = 2δ_b."""
import random, time
from lp_lab.fastlp import Ab, lp_ab_rows, k_from_rows
from lp_lab.design import Ideal, poly_int, mul, kernel_code, max_k_designs, b_table
from lp_lab import f2
rng = random.Random(5)
for orders in [(14,), (28,), (6, 6), (12, 3), (10,), (9, 4)]:
    G = Ab(orders); N = G.N
    for _ in range(25):
        p = (0, rng.randrange(1, N)); q = tuple(rng.sample(range(N), 2)); b = (0,) + tuple(rng.sample(range(1, N), 2))
        A = {(0, 0): p, (0, 1): (0,), (1, 0): (0,), (1, 1): q}
        rX, rZ = lp_ab_rows(G, A, 2, b)
        k = k_from_rows(rX, rZ, 4 * N, True)
        # dim R/(b, 1+pq) = N - rank(span of translates of b and of 1+pq)
        c = mul(G, poly_int(tuple(G.neg[g] for g in p)), poly_int(tuple(G.neg[g] for g in q))) ^ 1
        from lp_lab.design import translate
        rows = [translate(G, poly_int(b), g) for g in range(N)] + [translate(G, c, g) for g in range(N)]
        kk = 2 * (N - f2.rank(rows))
        assert k == kk, (orders, p, q, b, k, kk)
    print(orders, "k = 2·dim R/(b,1+p*q*) verified")
G = Ab((56,)); N = 56
t = time.time(); tab = b_table(G, delta_range=(3, 12), cap=14); print(f"b_table Z56: {len(tab)} classes ({time.time()-t:.1f}s)")
for bi in tab[:8]: print("  ", bi, G.poly_str(bi.b))
b = tab[0].b
t = time.time(); des = list(max_k_designs(G, b)); print(f"designs for b={b}: {len(des)} ({time.time()-t:.1f}s)")
for p, q in des[:3]:
    A = {(0, 0): p, (0, 1): (0,), (1, 0): (0,), (1, 1): q}
    rX, rZ = lp_ab_rows(G, A, 2, b); print("   ", p, q, "k =", k_from_rows(rX, rZ, 4 * N, True), "(expect", 2 * tab[0].delta, ")")
