"""A22 phase 0: the (eps,delta) CRT fibering of f2a6f17e.

Claims verified here (all asserted exactly; see notes/A22_analytic_classification.md):

  V1. With z = y^3, w = y^5:  F2[z]/(z^5-1) ~= F2 x GF16 via (p(1), p(zeta));
      the map is a bijection on the 32 fiber polynomials and induces the
      exact local weight table W(eps,delta) in {0,1,2,3,4,5}.
  V2. G = Z5(x) x Z15(y) ~= base Z15 (= Z5(x) x Z3(w)) x fiber Z5(z);
      every u in F2[G] decomposes into (u_eps in F2[Z15], u_delta in
      GF16[Z15]) with |u| = sum over the 15 base sites of W(u_eps, u_delta),
      and the boundary pair decomposes as
        u_eps = Abar * f_eps,  v_eps = xbar * Abar * f_eps   (Abar = 1+x+w^2)
        u_delta = Atil * f_del, v_delta = Btil * f_del       (GF16 coeffs)
      with (f_eps, f_del) INDEPENDENT free coordinates.
  V3. Abar is invertible in F2[Z15-base]  (no character zeros).
  V4. Atil-hat and Btil-hat have exactly ONE common zero eta0 among the 15
      GF16-valued characters of the base, and no other zeros; the transfer
      T = Btil-hat/Atil-hat has multiplicative orders in {1,3,5,15}.
  V5. tau (inverse DFT of T extended by 0 at eta0) satisfies tau*Atil = Btil.
  V6. All 113 enumerated classes decompose per the ansatz: v_eps = xbar u_eps,
      delta-parts satisfy the transfer relation, the site-weight formula
      reproduces b_weight; taxonomy (active sites, optimal cost, h-excess).
  V7. INDEPENDENT RE-DERIVATION: enumerate all light (h, alpha) via the
      <=7-active-site GF16 rank sweep + h-flip accounting; canonicalize by
      the 75 translations; the class set must equal the file's 113 exactly.

Usage: uv run --project experiments/bb_lab python experiments/bb_lab/scripts/a22_eps_delta_structure.py
"""

from __future__ import annotations

import itertools
import json
import sys
from collections import Counter
from pathlib import Path

import numpy as np

LAB_ROOT = Path(__file__).resolve().parent.parent
DATA = LAB_ROOT / "data" / "a17" / "f2a6_light_classes.jsonl"

ELL, M = 5, 15  # G = Z5(x) x Z15(y)

# ---------------------------------------------------------------- GF(16)
# F2[t]/(t^4+t+1), t primitive of order 15.  zeta := t^3 (order 5),
# omega := t^5 (order 3).
EXP = np.zeros(30, dtype=np.int64)
LOG = np.zeros(16, dtype=np.int64)
_v = 1
for _k in range(15):
    EXP[_k] = _v
    LOG[_v] = _k
    _v <<= 1
    if _v & 16:
        _v ^= 0b10011
EXP[15:30] = EXP[0:15]
ZETA = int(EXP[3])
OMEGA = int(EXP[5])


def gmul(a: int, b: int) -> int:
    if a == 0 or b == 0:
        return 0
    return int(EXP[(LOG[a] + LOG[b]) % 15])


def ginv(a: int) -> int:
    assert a != 0
    return int(EXP[(15 - LOG[a]) % 15])


def gpow(a: int, k: int) -> int:
    if a == 0:
        return 0
    return int(EXP[(LOG[a] * k) % 15])


ZP = [gpow(ZETA, k) for k in range(5)]   # zeta^k
WP = [gpow(OMEGA, k) for k in range(3)]  # omega^k
MU5 = set(ZP)                            # fifth roots of unity (incl 1)

# ------------------------------------------------- base sites & fibers
# j in Z15  <->  (a,b), j = 3a + 5b mod 15, a = 2j mod 5, b = 2j mod 3.
# base site = (i, b) in Z5 x Z3 (15 sites), fiber coordinate a in Z5.
SITES = [(i, b) for i in range(5) for b in range(3)]
SIDX = {s: k for k, s in enumerate(SITES)}


def jofab(a: int, b: int) -> int:
    return (3 * a + 5 * b) % 15


def site_sub(s1, s2):
    return ((s1[0] - s2[0]) % 5, (s1[1] - s2[1]) % 3)


def site_add(s1, s2):
    return ((s1[0] + s2[0]) % 5, (s1[1] + s2[1]) % 3)


def fiber_decompose(u: np.ndarray):
    """u: 5x15 F2 array (indexed [i, j]).  Returns (eps, delta):
    eps: len-15 int array (F2), delta: len-15 int array (GF16)."""
    eps = np.zeros(15, dtype=np.int64)
    dlt = np.zeros(15, dtype=np.int64)
    for k, (i, b) in enumerate(SITES):
        e, d = 0, 0
        for a in range(5):
            if u[i, jofab(a, b)]:
                e ^= 1
                d ^= ZP[a]
        eps[k], dlt[k] = e, d
    return eps, dlt


def fiber_glue(eps: np.ndarray, dlt: np.ndarray) -> np.ndarray:
    """Inverse of fiber_decompose (uses the unique-preimage property V1)."""
    u = np.zeros((5, 15), dtype=np.uint8)
    for k, (i, b) in enumerate(SITES):
        p = PREIMAGE[(int(eps[k]), int(dlt[k]))]
        for a in range(5):
            if (p >> a) & 1:
                u[i, jofab(a, b)] = 1
    return u


# ------------------------------------------------------- V1 weight table
PREIMAGE: dict[tuple[int, int], int] = {}
WTAB: dict[tuple[int, int], int] = {}
for p in range(32):
    e = bin(p).count("1") & 1
    d = 0
    for a in range(5):
        if (p >> a) & 1:
            d ^= ZP[a]
    key = (e, d)
    assert key not in PREIMAGE, "fiber (eps,delta) map not injective!"
    PREIMAGE[key] = p
    WTAB[key] = bin(p).count("1")
assert len(PREIMAGE) == 32
# the claimed table:
for d in range(16):
    if d == 0:
        assert WTAB[(0, 0)] == 0 and WTAB[(1, 0)] == 5
    elif d in MU5:
        assert WTAB[(1, d)] == 1 and WTAB[(0, d)] == 4, (d, WTAB[(1, d)], WTAB[(0, d)])
    else:
        assert WTAB[(1, d)] == 3 and WTAB[(0, d)] == 2
print("V1 PASS: (eps,delta) bijection + weight table (O:0/5, M:1/4, D:2/3)")


def site_type(d: int) -> str:
    return "O" if d == 0 else ("M" if d in MU5 else "D")


# per-(typeU,typeV) optimal-h cost and h*, flip penalty
def cost_and_h(alpha: int, beta: int) -> tuple[int, int, int]:
    """returns (cost, h*, flip_penalty)"""
    c0 = WTAB[(0, alpha)] + WTAB[(0, beta)]
    c1 = WTAB[(1, alpha)] + WTAB[(1, beta)]
    if c0 <= c1:
        return c0, 0, c1 - c0
    return c1, 1, c0 - c1


# ------------------------------------------------------------ upstairs conv
def conv_f2(P: np.ndarray, f: np.ndarray) -> np.ndarray:
    """(P*f)(g) = sum_h P(h) f(g-h) over Z5 x Z15, F2."""
    out = np.zeros((5, 15), dtype=np.uint8)
    for i in range(5):
        for j in range(15):
            if P[i, j]:
                out ^= np.roll(np.roll(f, i, axis=0), j, axis=1)
    return out


A_UP = np.zeros((5, 15), dtype=np.uint8)
for (i, j) in [(0, 0), (0, 1), (1, 0)]:          # 1 + y + x
    A_UP[i, j] = 1
B_UP = np.zeros((5, 15), dtype=np.uint8)
for (i, j) in [(1, 6), (1, 10), (2, 12)]:        # x y^6 + x y^10 + x^2 y^12
    B_UP[i, j] = 1

# ---------------------------------------------- downstairs polys (site level)
# Abar = 1 + xbar + wbar^2 over F2;  Atil = 1 + x + zeta^2 w^2,
# Btil = zeta^2 x + x w^2 + zeta^4 x^2 over GF16.
ABAR = {(0, 0): 1, (1, 0): 1, (0, 2): 1}
ATIL = {(0, 0): 1, (1, 0): 1, (0, 2): ZP[2]}
BTIL = {(1, 0): ZP[2], (1, 2): 1, (2, 0): ZP[4]}


def conv_site_gf16(P: dict, f: np.ndarray) -> np.ndarray:
    """(P*f)(s) = sum_{s'} P(s') f(s-s') over the 15 base sites, GF16."""
    out = np.zeros(15, dtype=np.int64)
    for sp, coef in P.items():
        for k, s in enumerate(SITES):
            out[k] ^= gmul(coef, int(f[SIDX[site_sub(s, sp)]]))
    return out


def conv_site_f2(P: dict, f: np.ndarray) -> np.ndarray:
    out = np.zeros(15, dtype=np.int64)
    for sp, coef in P.items():
        assert coef == 1
        for k, s in enumerate(SITES):
            out[k] ^= int(f[SIDX[site_sub(s, sp)]])
    return out


# ------------------------------------------------------------------ V2
rng = np.random.default_rng(2026)
for trial in range(50):
    f = rng.integers(0, 2, size=(5, 15)).astype(np.uint8)
    u = conv_f2(A_UP, f)
    v = conv_f2(B_UP, f)
    fe, fd = fiber_decompose(f)
    ue, ud = fiber_decompose(u)
    ve, vd = fiber_decompose(v)
    assert np.array_equal(ue, conv_site_f2(ABAR, fe)), "u_eps != Abar f_eps"
    assert np.array_equal(ud, conv_site_gf16(ATIL, fd)), "u_delta != Atil f_del"
    assert np.array_equal(vd, conv_site_gf16(BTIL, fd)), "v_delta != Btil f_del"
    # v_eps = xbar * u_eps
    sh = np.array([ue[SIDX[site_sub(s, (1, 0))]] for s in SITES])
    assert np.array_equal(ve, sh), "v_eps != xbar u_eps"
    # weight formula
    wu = sum(WTAB[(int(ue[k]), int(ud[k]))] for k in range(15))
    wv = sum(WTAB[(int(ve[k]), int(vd[k]))] for k in range(15))
    assert wu == int(u.sum()) and wv == int(v.sum()), "site weight formula"
    # glue round-trip
    assert np.array_equal(fiber_glue(ue, ud), u)
print("V2 PASS: decomposition, independence-form, exact site weight formula (50 rand)")

# ------------------------------------------------------------------ V3, V4
def dft_site(P: dict):
    """hat(P)(p,q) = sum_s P(s) zeta^{p i} omega^{q b}; 5x3 GF16 array."""
    out = np.zeros((5, 3), dtype=np.int64)
    for (i, b), coef in P.items():
        for p in range(5):
            for q in range(3):
                out[p, q] ^= gmul(coef, gmul(ZP[(p * i) % 5], WP[(q * b) % 3]))
    return out


AB_HAT = dft_site(ABAR)
assert np.all(AB_HAT != 0), "Abar has a character zero!"
print("V3 PASS: Abar invertible in F2[Z15-base] (15 nonzero character values)")

AT_HAT = dft_site(ATIL)
BT_HAT = dft_site(BTIL)
zerosA = [(p, q) for p in range(5) for q in range(3) if AT_HAT[p, q] == 0]
zerosB = [(p, q) for p in range(5) for q in range(3) if BT_HAT[p, q] == 0]
assert zerosA == zerosB and len(zerosA) == 1, (zerosA, zerosB)
ETA0 = zerosA[0]
orders = Counter()
T_HAT = np.zeros((5, 3), dtype=np.int64)
for p in range(5):
    for q in range(3):
        if (p, q) == ETA0:
            continue
        T_HAT[p, q] = gmul(int(BT_HAT[p, q]), ginv(int(AT_HAT[p, q])))
        orders[15 // np.gcd(LOG[T_HAT[p, q]], 15)] += 1
print(f"V4 PASS: unique common character zero eta0 = {ETA0}; T orders: {dict(orders)}")

# ------------------------------------------------------------------ V5 tau
TAU = np.zeros(15, dtype=np.int64)
for k, (i, b) in enumerate(SITES):
    acc = 0
    for p in range(5):
        for q in range(3):
            acc ^= gmul(int(T_HAT[p, q]),
                        gmul(ZP[(-p * i) % 5], WP[(-q * b) % 3]))
    TAU[k] = acc
TAU_D = {SITES[k]: int(TAU[k]) for k in range(15) if TAU[k]}
tauA = conv_site_gf16(TAU_D, np.array([ATIL.get(s, 0) for s in SITES]))
btil_vec = np.array([BTIL.get(s, 0) for s in SITES])
assert np.array_equal(tauA, btil_vec), "tau * Atil != Btil"
print(f"V5 PASS: tau*Atil = Btil; |supp tau| = {len(TAU_D)}")

# beta'_g = v_delta(g + (1,0)) = (tau * alpha)(g + (1,0))
THETA = np.zeros((15, 15), dtype=np.int64)  # beta' = THETA @ alpha (GF16)
for kg, g in enumerate(SITES):
    tgt = site_add(g, (1, 0))
    for kgp, gp in enumerate(SITES):
        THETA[kg, kgp] = int(TAU[SIDX[site_sub(tgt, gp)]])

# eta0 row: alpha-hat(eta0) = sum_s alpha(s) zeta^{p0 i} omega^{q0 b}
p0, q0 = ETA0
ETA0_ROW = np.array([gmul(ZP[(p0 * i) % 5], WP[(q0 * b) % 3])
                     for (i, b) in SITES], dtype=np.int64)


def apply_theta(alpha: np.ndarray) -> np.ndarray:
    out = np.zeros(15, dtype=np.int64)
    for kg in range(15):
        acc = 0
        for kp in range(15):
            acc ^= gmul(int(THETA[kg, kp]), int(alpha[kp]))
        out[kg] = acc
    return out


# ------------------------------------------------------------------ V6
classes = [json.loads(line) for line in DATA.read_text().splitlines() if line.strip()]
classes = [c for c in classes if "b_weight" in c]  # drop the completion marker
assert len(classes) == 113, len(classes)
wspec = Counter(c["b_weight"] for c in classes)
nk = [c for c in classes if c["coset_min"] >= 5]
print(f"file check: 113 classes, weights {dict(wspec)}, near-kernel {len(nk)} "
      f"(coset minima {sorted(set(c['coset_min'] for c in classes))})")

taxonomy = Counter()
class_canon = set()


def canon_uv(u: np.ndarray, v: np.ndarray) -> bytes:
    best = None
    for di in range(5):
        for dj in range(15):
            uu = np.roll(np.roll(u, di, 0), dj, 1)
            vv = np.roll(np.roll(v, di, 0), dj, 1)
            cand = uu.tobytes() + vv.tobytes()
            if best is None or cand < best:
                best = cand
    return best


for c in classes:
    u = np.zeros((5, 15), dtype=np.uint8)
    v = np.zeros((5, 15), dtype=np.uint8)
    for blk, i, j in c["b_support"]:
        (u if blk == 0 else v)[i, j] = 1
    assert int(u.sum()) == c["u_weight"] and int(v.sum()) == c["v_weight"]
    ue, ud = fiber_decompose(u)
    ve, vd = fiber_decompose(v)
    # v_eps = xbar u_eps
    sh = np.array([ue[SIDX[site_sub(s, (1, 0))]] for s in SITES])
    assert np.array_equal(ve, sh), "class violates v_eps = xbar u_eps!"
    # transfer relation on delta parts: vd = tau * ud
    assert np.array_equal(apply_theta(ud),
                          np.array([vd[SIDX[site_add(g, (1, 0))]] for g in SITES])), \
        "class violates delta transfer"
    # eta0 condition
    acc = 0
    for k in range(15):
        acc ^= gmul(int(ETA0_ROW[k]), int(ud[k]))
    assert acc == 0, "class violates alpha-hat(eta0) = 0"
    # site weight formula
    wu = sum(WTAB[(int(ue[k]), int(ud[k]))] for k in range(15))
    wv = sum(WTAB[(int(ve[k]), int(vd[k]))] for k in range(15))
    assert wu + wv == c["b_weight"]
    # taxonomy: active sites of (alpha, beta'), optimal cost, excess
    alpha = ud
    beta = np.array([vd[SIDX[site_add(g, (1, 0))]] for g in SITES])
    S = [k for k in range(15) if alpha[k] or beta[k]]
    cost = sum(cost_and_h(int(alpha[k]), int(beta[k]))[0] for k in S)
    excess = c["b_weight"] - cost
    prof = "".join(sorted(site_type(int(alpha[k])) + site_type(int(beta[k])) for k in S))
    taxonomy[(c["b_weight"], len(S), cost, excess, prof,
              c["coset_min"] >= 5)] += 1
    class_canon.add(canon_uv(u, v))

assert len(class_canon) == 113, "file classes not translation-distinct?!"
print("V6 PASS: all 113 classes decompose per ansatz; taxonomy "
      "(|b|, |S|, cost, excess, profile, near-kernel?):")
for key, cnt in sorted(taxonomy.items()):
    print(f"   {key}: {cnt}")

# ------------------------------------------------------------------ V7
print("V7: independent re-derivation via <=7-site GF16 rank sweep ...")


def gf16_kernel(mat: np.ndarray) -> list[np.ndarray]:
    """Kernel basis of mat (rows x cols) over GF16."""
    m = mat.copy()
    rows, cols = m.shape
    pivots = []
    r = 0
    for c in range(cols):
        piv = None
        for rr in range(r, rows):
            if m[rr, c]:
                piv = rr
                break
        if piv is None:
            continue
        m[[r, piv]] = m[[piv, r]]
        inv = ginv(int(m[r, c]))
        m[r] = [gmul(inv, int(x)) for x in m[r]]
        for rr in range(rows):
            if rr != r and m[rr, c]:
                coef = int(m[rr, c])
                m[rr] = [int(m[rr, k]) ^ gmul(coef, int(m[r, k])) for k in range(cols)]
        pivots.append(c)
        r += 1
        if r == rows:
            break
    free = [c for c in range(cols) if c not in pivots]
    basis = []
    for fc in free:
        vec = np.zeros(cols, dtype=np.int64)
        vec[fc] = 1
        for ri, pc in enumerate(pivots):
            vec[pc] = int(m[ri, fc])  # x_p = -sum m[r,free]*x_free; char 2
        basis.append(vec)
    return basis


# collect candidate alpha (delta-side configs), dedupe raw
alpha_set: set[tuple] = set()
n_special = 0
for size in range(1, 8):
    for S in itertools.combinations(range(15), size):
        Sset = set(S)
        comp = [k for k in range(15) if k not in Sset]
        # rows: beta' on comp must vanish; plus eta0 row; cols: alpha on S
        mat = np.zeros((len(comp) + 1, size), dtype=np.int64)
        for ri, k in enumerate(comp):
            for ci, kp in enumerate(S):
                mat[ri, ci] = THETA[k, kp]
        for ci, kp in enumerate(S):
            mat[len(comp), ci] = ETA0_ROW[kp]
        basis = gf16_kernel(mat)
        if not basis:
            continue
        n_special += 1
        dim = len(basis)
        for coeffs in itertools.product(range(16), repeat=dim):
            if all(cc == 0 for cc in coeffs):
                continue
            alpha = np.zeros(15, dtype=np.int64)
            for cc, vec in zip(coeffs, basis):
                if cc:
                    for ci, kp in enumerate(S):
                        alpha[kp] ^= gmul(cc, int(vec[ci]))
            if not alpha.any():
                continue
            alpha_set.add(tuple(int(x) for x in alpha))
print(f"   special site-sets (nonzero kernel): {n_special}; "
      f"distinct nonzero alpha with <=7 active sites: {len(alpha_set)}")

# filter by cost, translation-dedupe alpha level
def canon_alpha(alpha: tuple) -> tuple:
    best = None
    for di in range(5):
        for da in range(5):
            for db in range(3):
                sc = ZP[da]
                cand = tuple(
                    gmul(sc, alpha[SIDX[site_sub(s, (di, db))]]) for s in SITES)
                if best is None or cand < best:
                    best = cand
    return best


alpha_reps: dict[tuple, tuple] = {}
n_light_alpha = 0
for alpha in alpha_set:
    a = np.array(alpha, dtype=np.int64)
    beta = apply_theta(a)
    cost = sum(cost_and_h(int(a[k]), int(beta[k]))[0]
               for k in range(15) if a[k] or beta[k])
    if cost > 14:
        continue
    n_light_alpha += 1
    rep = canon_alpha(alpha)
    if rep not in alpha_reps:
        alpha_reps[rep] = alpha
print(f"   alpha configs with optimal cost <=14: {n_light_alpha}; "
      f"translation classes: {len(alpha_reps)}")

# h-flip enumeration around each alpha rep + the pure-h sector
found: set[bytes] = set()


def emit(alpha: np.ndarray, h: np.ndarray):
    beta = apply_theta(alpha)
    ue = h
    ud = alpha
    ve = np.array([h[SIDX[site_sub(s, (1, 0))]] for s in SITES])
    vd = np.zeros(15, dtype=np.int64)
    for kg, g in enumerate(SITES):
        vd[SIDX[site_add(g, (1, 0))]] = beta[kg]
    w = sum(WTAB[(int(ue[k]), int(ud[k]))] for k in range(15)) + \
        sum(WTAB[(int(ve[k]), int(vd[k]))] for k in range(15))
    if w == 0 or w > 14:
        return
    u = fiber_glue(ue, ud)
    v = fiber_glue(ve, vd)
    found.add(canon_uv(u, v))


# pure-h sector (alpha = 0): |h| = 1 only (10 per site pair)
h = np.zeros(15, dtype=np.int64)
h[0] = 1
emit(np.zeros(15, dtype=np.int64), h)

for rep, alpha in alpha_reps.items():
    a = np.array(alpha, dtype=np.int64)
    beta = apply_theta(a)
    hstar = np.zeros(15, dtype=np.int64)
    pens = []
    cost = 0
    for k in range(15):
        ck, hk, pk = cost_and_h(int(a[k]), int(beta[k]))
        cost += ck
        hstar[k] = hk
        pens.append(pk)
    budget = 14 - cost
    assert budget >= 0
    # enumerate flip sets F with sum of penalties <= budget
    flippable = [k for k in range(15) if pens[k] <= budget]
    # simple DFS
    def dfs(idx: int, rem: int, flips: list[int]):
        hh = hstar.copy()
        for k in flips:
            hh[k] ^= 1
        emit(a, hh)
        for pos in range(idx, len(flippable)):
            k = flippable[pos]
            if pens[k] <= rem:
                dfs(pos + 1, rem - pens[k], flips + [k])
    dfs(0, budget, [])

print(f"   re-derived class count: {len(found)}")
assert found == class_canon, (
    f"MISMATCH: derived {len(found)} vs file 113; "
    f"missing {len(class_canon - found)}, extra {len(found - class_canon)}")
print("V7 PASS: (h,alpha) enumeration reproduces the 113 classes EXACTLY")
print("ALL PASS")
