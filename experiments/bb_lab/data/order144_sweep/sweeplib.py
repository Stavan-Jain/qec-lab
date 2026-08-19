"""Order-144 breadth sweep — shared machinery.

Fast (numpy) canonical form for rank-2 abelian groups under the SAME
equivalence bb-lab enumerate uses (Aut x translation x block-swap,
lex-min bitset rep), exact for weight-3 x weight-3 pairs; verified
against bb_lab.canonical.canonical_bits on small groups by smoke.py.

Also: the canonical exponent-reduction quotient ladder
(order 144 -> 72 -> 36 -> 18, halve axis 1 if even else axis 0).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from bb_lab.group import AbelianGroup
from bb_lab.poly import Poly
from bb_lab.checks import bb_check_matrices
from bb_lab.codeparams import code_params
from bb_lab.store import canonical_hash


# --------------------------------------------------------------------------
# fast Aut x translation x swap canonicalization (weight-3 exact)


@dataclass
class GroupTables:
    orders: tuple[int, int]
    N: int
    PHI: np.ndarray   # (n_aut, N) int32: PHI[j, i] = index(phi_j(elem_i))
    ADD: np.ndarray   # (N, N) int32: ADD[i, h] = index(elem_i + elem_h)


_TABLES: dict[tuple[int, int], GroupTables] = {}


def _index_grid(L: int, M: int) -> tuple[np.ndarray, np.ndarray]:
    """(a, b) coordinate arrays for row-major index i = a*M + b."""
    i = np.arange(L * M)
    return i // M, i % M


def group_tables(L: int, M: int) -> GroupTables:
    key = (L, M)
    if key in _TABLES:
        return _TABLES[key]
    N = L * M
    a, b = _index_grid(L, M)  # (N,), (N,)

    # Candidate generator images. img0 = phi(e0) needs order | L:
    # L*(u,v) = 0  <=>  (L*v) % M == 0  (L*u % L == 0 always).
    # img1 = phi(e1) needs order | M: (M*u) % L == 0.
    all_u, all_v = _index_grid(L, M)
    cand0 = np.flatnonzero((L * all_v) % M == 0)   # indices usable as img0
    cand1 = np.flatnonzero((M * all_u) % L == 0)   # indices usable as img1
    u0, v0 = all_u[cand0], all_v[cand0]            # (C0,)
    u1, v1 = all_u[cand1], all_v[cand1]            # (C1,)

    # phi((a,b)) = ((a*u0 + b*u1) % L, (a*v0 + b*v1) % M), all pairs.
    # Shape juggling: pairs axis = (C0, C1) flattened.
    U = (a[None, None, :] * u0[:, None, None] + b[None, None, :] * u1[None, :, None]) % L
    V = (a[None, None, :] * v0[:, None, None] + b[None, None, :] * v1[None, :, None]) % M
    IDX = (U * M + V).reshape(-1, N)               # (C0*C1, N)
    # Bijective rows = automorphisms (order conditions make them homs).
    srt = np.sort(IDX, axis=1)
    bij = (srt == np.arange(N)[None, :]).all(axis=1)
    PHI = IDX[bij].astype(np.int32)

    AA = ((a[:, None] + a[None, :]) % L) * M + (b[:, None] + b[None, :]) % M
    ADD = AA.astype(np.int32)
    gt = GroupTables(orders=key, N=N, PHI=PHI, ADD=ADD)
    _TABLES[key] = gt
    return gt


def _desc_key(T: np.ndarray) -> np.ndarray:
    """T: (..., 3) index triples -> int64 key preserving bitset order
    (descending-sorted triple, base 2^9 digits; N <= 512)."""
    S = np.sort(T, axis=-1)[..., ::-1]
    return (S[..., 0].astype(np.int64) << 18) | (S[..., 1].astype(np.int64) << 9) | S[..., 2].astype(np.int64)


def canonical_pair_fast(
    L: int, M: int, A_idx: np.ndarray, B_idx: np.ndarray
) -> tuple[tuple[int, ...], tuple[int, ...], int]:
    """Lex-min (A_bits, B_bits) orbit rep under Aut x translation x swap,
    for weight-3 supports given as index triples. Returns
    (A_indices_desc, B_indices_desc, orbit_size) of the canonical rep.
    Exactly reproduces bb_lab.canonical.canonical_bits for weight-3 pairs
    (bitset int order == lex order on descending index triples)."""
    gt = group_tables(L, M)
    N = gt.N
    h = np.arange(N)[None, :, None]
    phiA = gt.PHI[:, A_idx]                       # (n_aut, 3)
    phiB = gt.PHI[:, B_idx]
    TA = gt.ADD[phiA[:, None, :], h]              # (n_aut, N, 3)
    TB = gt.ADD[phiB[:, None, :], h]
    kA = _desc_key(TA)                            # (n_aut, N)
    kB = _desc_key(TB)
    K1 = (kA << 27) | kB                          # (A,B) orientation
    K2 = (kB << 27) | kA                          # swapped
    orbit = np.unique(np.concatenate([K1.ravel(), K2.ravel()]))
    best = int(orbit[0])
    kb_best = best & ((1 << 27) - 1)
    ka_best = best >> 27
    dec = lambda k: tuple(int(x) for x in ((k >> 18) & 511, (k >> 9) & 511, k & 511))
    return dec(ka_best), dec(kb_best), int(orbit.size)


def idx_to_support(idx_desc: tuple[int, ...], L: int, M: int) -> tuple[tuple[int, int], ...]:
    return tuple((i // M, i % M) for i in idx_desc)


def canonicalize(L: int, M: int, A_supp, B_supp):
    """Full canonicalization: supports -> canonical Poly strings + id + orbit size."""
    G = AbelianGroup((L, M))
    A_idx = np.array([a * M + b for a, b in A_supp])
    B_idx = np.array([a * M + b for a, b in B_supp])
    cA, cB, orbit = canonical_pair_fast(L, M, A_idx, B_idx)
    A_str = Poly(support=frozenset(idx_to_support(cA, L, M)), group=G).canonical_string()
    B_str = Poly(support=frozenset(idx_to_support(cB, L, M)), group=G).canonical_string()
    label = G.label()
    iid = canonical_hash(label, A_str, B_str)
    return label, A_str, B_str, iid, orbit


# --------------------------------------------------------------------------
# quotient ladder (exponent reduction)


def halve_axis(orders: tuple[int, int]) -> int:
    """Canonical descent convention: halve axis 1 if its order is even,
    else axis 0 (which must then be even)."""
    if orders[1] % 2 == 0:
        return 1
    assert orders[0] % 2 == 0, f"no even axis in {orders}"
    return 0


def quotient_chain(orders: tuple[int, int], depth: int = 3) -> list[tuple[int, int]]:
    """The canonical exponent-reduction chain, `depth` halvings."""
    out = []
    cur = orders
    for _ in range(depth):
        ax = halve_axis(cur)
        nxt = list(cur)
        nxt[ax] = cur[ax] // 2
        cur = (nxt[0], nxt[1])
        out.append(cur)
    return out


def reduce_support(supp, axis: int, new_len: int):
    """Exponent reduction mod new_len on `axis`, with mod-2 cancellation."""
    from collections import Counter
    c = Counter()
    for g in supp:
        e = list(g)
        e[axis] = e[axis] % new_len
        c[tuple(e)] += 1
    return tuple(sorted(g for g, cnt in c.items() if cnt % 2 == 1))


def quotient_code(orders: tuple[int, int], A_supp, B_supp, q_orders: tuple[int, int]):
    """Push (A, B) down one halving step to q_orders (which must differ
    from `orders` in exactly one axis, by factor 2)."""
    ax = 0 if orders[0] != q_orders[0] else 1
    A_q = reduce_support(A_supp, ax, q_orders[ax])
    B_q = reduce_support(B_supp, ax, q_orders[ax])
    return A_q, B_q


def code_params_of(orders: tuple[int, int], A_supp, B_supp):
    """(n, k, checks) for a support pair on `orders`."""
    G = AbelianGroup(orders)
    A = Poly(support=frozenset(tuple(g) for g in A_supp), group=G)
    B = Poly(support=frozenset(tuple(g) for g in B_supp), group=G)
    checks = bb_check_matrices(A, B)
    p = code_params(checks)
    return p, checks, A, B


def poly_string(orders: tuple[int, int], supp) -> str:
    G = AbelianGroup(orders)
    return Poly(support=frozenset(tuple(g) for g in supp), group=G).canonical_string()


def support_from_string(orders: tuple[int, int], s: str):
    G = AbelianGroup(orders)
    return tuple(sorted(Poly.from_string(s, G).support))
