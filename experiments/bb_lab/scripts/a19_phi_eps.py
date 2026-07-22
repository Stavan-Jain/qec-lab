"""A19 fibering session 2a: (A) Phi-sparsity spectrum, (B) exhaustive eps-layer.

(A) Phi = Bt * At^{-1} as a GF(16)[H] convolution element phi; joint site
    sparsity min(|supp U| + |supp Phi U|) for |supp U| <= 3 — the pure-delta
    sector floor input.
(B) All 2^18 f_eps: pair patterns (A_e f, B_e f), exact pair-weight
    distribution and minimum nonzero pair weight — the eps-layer census,
    exhaustive by construction.
"""
import sys
from pathlib import Path
import numpy as np
from itertools import combinations

# GF(16) tables (t^4 = t + 1)
EXP = [1]
for _ in range(14):
    e = EXP[-1] << 1
    if e & 0x10:
        e ^= 0x13
    EXP.append(e)
LOG = {e: i for i, e in enumerate(EXP)}


def gmul(a, b):
    return 0 if (a == 0 or b == 0) else EXP[(LOG[a] + LOG[b]) % 15]


ZETA = EXP[3]
H_SITES = [(q, b) for q in range(6) for b in range(3)]
SITE_IDX = {s: i for i, s in enumerate(H_SITES)}
N = 18


def conv_mat_gf16(supp):
    M = [[0] * N for _ in range(N)]
    for j, (qj, bj) in enumerate(H_SITES):
        for (dq, db), c in supp.items():
            M[SITE_IDX[((qj + dq) % 6, (bj + db) % 3)]][j] ^= c
    return M


def gf16_inv(Mat):
    n = N
    A = [row[:] + [1 if i == j else 0 for j in range(n)]
         for i, row in enumerate(Mat)]
    r = 0
    for c in range(n):
        piv = next(i for i in range(r, n) if A[i][c])
        A[r], A[piv] = A[piv], A[r]
        inv = EXP[(15 - LOG[A[r][c]]) % 15]
        A[r] = [gmul(x, inv) for x in A[r]]
        for i in range(n):
            if i != r and A[i][c]:
                cc = A[i][c]
                A[i] = [x ^ gmul(cc, y) for x, y in zip(A[i], A[r])]
        r += 1
    return [row[n:] for row in A]


def matvec(M, v):
    out = [0] * N
    for i in range(N):
        acc = 0
        Mi = M[i]
        for j in range(N):
            if v[j]:
                acc ^= gmul(Mi[j], v[j])
        out[i] = acc
    return out


A_d = conv_mat_gf16({(3, 0): EXP[12], (0, 1): 1, (0, 2): 1})
B_d = conv_mat_gf16({(0, 0): 1, (1, 0): 1, (2, 0): ZETA})
PHI = [[0] * N for _ in range(N)]
Ainv = gf16_inv(A_d)
for i in range(N):
    for j in range(N):
        acc = 0
        for k in range(N):
            acc ^= gmul(B_d[i][k], Ainv[k][j])
        PHI[i][j] = acc

phi_col = [PHI[i][0] for i in range(N)]
print(f"(A) |supp phi| = {sum(1 for x in phi_col if x)} / 18")

best = {1: 99, 2: 99, 3: 99}
for a in (1, 2, 3):
    for pos in combinations(range(N), a):
        # projective coefficients: first = 1, rest range over GF(16)*
        from itertools import product
        for coeffs in product(*([range(1, 16)] * (a - 1))):
            U = [0] * N
            U[pos[0]] = 1
            for p, c in zip(pos[1:], coeffs):
                U[p] = c
            V = matvec(PHI, U)
            joint = a + sum(1 for x in V if x)
            if joint < best[a]:
                best[a] = joint
    print(f"    min |supp U| + |supp Phi U| over |supp U| = {a}: {best[a]}")
print(f"    => pure-delta joint-sparsity floor (a <= 3): "
      f"{min(best.values())}; site cost >= 2 per delta-only site "
      f"=> pure-delta stabilizer weight >= {2 * min(best.values())}")

# ---------------- (B) exhaustive eps-layer ----------------
def conv_mat_f2(supp):
    M = np.zeros((N, N), dtype=np.uint8)
    for j, (qj, bj) in enumerate(H_SITES):
        for (dq, db) in supp:
            M[SITE_IDX[((qj + dq) % 6, (bj + db) % 3)], j] ^= 1
    return M


Ae = conv_mat_f2([(3, 0), (0, 1), (0, 2)])
Be = conv_mat_f2([(0, 0), (1, 0), (2, 0)])

bits = ((np.arange(1 << 18)[:, None] >> np.arange(18)[None, :]) & 1
        ).astype(np.uint8)
U_all = (bits @ Ae.T) % 2
W_all = (bits @ Be.T) % 2
wt = U_all.sum(1) + W_all.sum(1)
nz = wt > 0
print("(B) eps-layer exhaustive (2^18 f_eps):")
kerdim = int(np.log2((~nz).sum()))
print(f"    ker(pair map) dim = {kerdim} -> eps-pair code [36, {18 - kerdim}]")
hist = np.bincount(wt[nz])
print(f"    min nonzero eps-pair weight = {np.flatnonzero(hist)[0]}")
print("    eps-pair weight histogram (w: count) up to 23:")
print("   ", {int(w): int(hist[w]) for w in np.flatnonzero(hist) if w <= 23})
