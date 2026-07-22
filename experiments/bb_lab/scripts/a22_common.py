"""A22 shared machinery: the (eps,delta) CRT fibering of f2a6f17e.

Verified end-to-end by a22_eps_delta_structure.py (V1-V7, ALL PASS).
This module re-exposes the delta-side objects for the structure scripts,
plus a cached recompute of the 94 alpha translation classes.

Import as:  from a22_common import *   (scripts add their dir to sys.path)
"""

from __future__ import annotations

import itertools
import json
from pathlib import Path

import numpy as np

LAB_ROOT = Path(__file__).resolve().parent.parent
A22_DATA = LAB_ROOT / "data" / "a22"

# ---------------------------------------------------------------- GF(16)
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


def gfrob(a: int, r: int = 1) -> int:
    """a^(2^r) — the Galois action."""
    if a == 0:
        return 0
    return int(EXP[(LOG[a] * pow(2, r, 15)) % 15])


ZP = [gpow(ZETA, k) for k in range(5)]
WP = [gpow(OMEGA, k) for k in range(3)]
MU5 = set(ZP)

# ------------------------------------------------- base sites & fibers
SITES = [(i, b) for i in range(5) for b in range(3)]
SIDX = {s: k for k, s in enumerate(SITES)}


def jofab(a: int, b: int) -> int:
    return (3 * a + 5 * b) % 15


def site_sub(s1, s2):
    return ((s1[0] - s2[0]) % 5, (s1[1] - s2[1]) % 3)


def site_add(s1, s2):
    return ((s1[0] + s2[0]) % 5, (s1[1] + s2[1]) % 3)


# weight table W(eps, delta) and fiber preimages
PREIMAGE: dict[tuple[int, int], int] = {}
WTAB: dict[tuple[int, int], int] = {}
for _p in range(32):
    _e = bin(_p).count("1") & 1
    _d = 0
    for _a in range(5):
        if (_p >> _a) & 1:
            _d ^= ZP[_a]
    PREIMAGE[(_e, _d)] = _p
    WTAB[(_e, _d)] = bin(_p).count("1")


def site_type(d: int) -> str:
    return "O" if d == 0 else ("M" if d in MU5 else "D")


def cost_and_h(alpha: int, beta: int) -> tuple[int, int, int]:
    """(optimal cost, h*, flip penalty) for one site pair."""
    c0 = WTAB[(0, alpha)] + WTAB[(0, beta)]
    c1 = WTAB[(1, alpha)] + WTAB[(1, beta)]
    if c0 <= c1:
        return c0, 0, c1 - c0
    return c1, 1, c0 - c1


# ------------------------------------------------- delta-side polynomials
ABAR = {(0, 0): 1, (1, 0): 1, (0, 2): 1}
ATIL = {(0, 0): 1, (1, 0): 1, (0, 2): ZP[2]}
BTIL = {(1, 0): ZP[2], (1, 2): 1, (2, 0): ZP[4]}
TRIANGLE = [(0, 0), (1, 0), (0, 2)]  # supp(Abar) = supp(Atil)


def conv_site_gf16(P: dict, f: np.ndarray) -> np.ndarray:
    out = np.zeros(15, dtype=np.int64)
    for sp, coef in P.items():
        for k, s in enumerate(SITES):
            out[k] ^= gmul(coef, int(f[SIDX[site_sub(s, sp)]]))
    return out


def dft_site(P: dict) -> np.ndarray:
    out = np.zeros((5, 3), dtype=np.int64)
    for (i, b), coef in P.items():
        for p in range(5):
            for q in range(3):
                out[p, q] ^= gmul(coef, gmul(ZP[(p * i) % 5], WP[(q * b) % 3]))
    return out


def dft_vec(f: np.ndarray) -> np.ndarray:
    """DFT of a site vector: hat(f)(p,q) = sum_s f(s) zeta^{pi} omega^{qb}."""
    out = np.zeros((5, 3), dtype=np.int64)
    for k, (i, b) in enumerate(SITES):
        if f[k]:
            for p in range(5):
                for q in range(3):
                    out[p, q] ^= gmul(int(f[k]), gmul(ZP[(p * i) % 5], WP[(q * b) % 3]))
    return out


def idft_vec(fh: np.ndarray) -> np.ndarray:
    out = np.zeros(15, dtype=np.int64)
    for k, (i, b) in enumerate(SITES):
        acc = 0
        for p in range(5):
            for q in range(3):
                acc ^= gmul(int(fh[p, q]), gmul(ZP[(-p * i) % 5], WP[(-q * b) % 3]))
        out[k] = acc
    return out


AT_HAT = dft_site(ATIL)
BT_HAT = dft_site(BTIL)
ETA0 = (4, 2)
assert AT_HAT[ETA0] == 0 and BT_HAT[ETA0] == 0

T_HAT = np.zeros((5, 3), dtype=np.int64)
for _pq in itertools.product(range(5), range(3)):
    if _pq != ETA0:
        T_HAT[_pq] = gmul(int(BT_HAT[_pq]), ginv(int(AT_HAT[_pq])))

TAU = idft_vec(T_HAT)

# beta'_g = v_delta(g+(1,0)) = (tau * alpha)(g+(1,0));  beta' = THETA @ alpha
THETA = np.zeros((15, 15), dtype=np.int64)
for _kg, _g in enumerate(SITES):
    _tgt = site_add(_g, (1, 0))
    for _kp, _gp in enumerate(SITES):
        THETA[_kg, _kp] = int(TAU[SIDX[site_sub(_tgt, _gp)]])

_p0, _q0 = ETA0
ETA0_ROW = np.array([gmul(ZP[(_p0 * i) % 5], WP[(_q0 * b) % 3])
                     for (i, b) in SITES], dtype=np.int64)


def apply_theta(alpha: np.ndarray) -> np.ndarray:
    out = np.zeros(15, dtype=np.int64)
    for kg in range(15):
        acc = 0
        for kp in range(15):
            acc ^= gmul(int(THETA[kg, kp]), int(alpha[kp]))
        out[kg] = acc
    return out


def config_cost(alpha: np.ndarray) -> int:
    beta = apply_theta(alpha)
    return sum(cost_and_h(int(alpha[k]), int(beta[k]))[0]
               for k in range(15) if alpha[k] or beta[k])


def active_set(alpha: np.ndarray) -> list[int]:
    beta = apply_theta(alpha)
    return [k for k in range(15) if alpha[k] or beta[k]]


def canon_alpha(alpha) -> tuple:
    """Canonical rep under the 75 translations (site shift x mu5 scalar)."""
    best = None
    for di in range(5):
        for da in range(5):
            for db in range(3):
                sc = ZP[da]
                cand = tuple(
                    gmul(sc, int(alpha[SIDX[site_sub(s, (di, db))]])) for s in SITES)
                if best is None or cand < best:
                    best = cand
    return best


def gf16_kernel(mat: np.ndarray) -> list[np.ndarray]:
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
            vec[pc] = int(m[ri, fc])
        basis.append(vec)
    return basis


def compute_alpha_classes(force: bool = False) -> list[tuple]:
    """The 94 canonical alpha translation classes with optimal cost <= 14,
    via the <=7-site sweep.  Cached in data/a22/alpha_reps.json."""
    cache = A22_DATA / "alpha_reps.json"
    if cache.exists() and not force:
        return [tuple(r) for r in json.loads(cache.read_text())]
    alpha_set: set[tuple] = set()
    for size in range(1, 8):
        for S in itertools.combinations(range(15), size):
            comp = [k for k in range(15) if k not in S]
            mat = np.zeros((len(comp) + 1, size), dtype=np.int64)
            for ri, k in enumerate(comp):
                for ci, kp in enumerate(S):
                    mat[ri, ci] = THETA[k, kp]
            for ci, kp in enumerate(S):
                mat[len(comp), ci] = ETA0_ROW[kp]
            basis = gf16_kernel(mat)
            if not basis:
                continue
            for coeffs in itertools.product(range(16), repeat=len(basis)):
                if all(cc == 0 for cc in coeffs):
                    continue
                alpha = np.zeros(15, dtype=np.int64)
                for cc, vec in zip(coeffs, basis):
                    if cc:
                        for ci, kp in enumerate(S):
                            alpha[kp] ^= gmul(cc, int(vec[ci]))
                if alpha.any():
                    alpha_set.add(tuple(int(x) for x in alpha))
    reps: dict[tuple, None] = {}
    for alpha in alpha_set:
        if config_cost(np.array(alpha, dtype=np.int64)) <= 14:
            reps[canon_alpha(alpha)] = None
    out = sorted(reps.keys())
    A22_DATA.mkdir(parents=True, exist_ok=True)
    cache.write_text(json.dumps([list(r) for r in out]))
    return out
