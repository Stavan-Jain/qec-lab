"""Fast bit-packed builder + k for the a=2 class LP(A, b) over Z_ℓ with
A = [[p, 1], [1, q]] (p, q weight-2, b weight-3).  Qubit/check indexing:
  L-qubits: block j∈{0,1}, element h  -> column j*ℓ + h
  R-qubits: block i∈{0,1}, element h  -> column 2ℓ + i*ℓ + h
  X-checks: (i, g);  Z-checks: (j, g).
H_X = [A | b* I],  H_Z = [b I | A†].  X-check (i,g) hits L(j, g−s) for s∈A[i,j]
and R(i, g+s) for s∈b  (since (b* I)[g,h]=1 iff g−h ∈ −b iff h = g+s).
"""
from __future__ import annotations

from . import f2


def rot(v: int, g: int, ell: int, mask: int) -> int:
    g %= ell
    if g == 0:
        return v
    return ((v << g) | (v >> (ell - g))) & mask


def pattern(poly: tuple[int, ...], ell: int, negate: bool) -> int:
    """Row pattern at g=0: bits h = -s (negate=True) or h = s."""
    v = 0
    for s in poly:
        h = (-s) % ell if negate else s % ell
        v |= 1 << h
    return v


def hx_rows(ell: int, p, q, b) -> list[int]:
    mask = (1 << ell) - 1
    # X-check (0,g): L0 via p (h = g - s), L1 via 1 (h = g), R0 via b (h = g + s)
    # X-check (1,g): L0 via 1, L1 via q, R1 via b
    P = pattern(p, ell, True); Q = pattern(q, ell, True); One = 1; Bp = pattern(b, ell, False)
    rows = []
    for g in range(ell):
        r0 = rot(P, g, ell, mask) | (rot(One, g, ell, mask) << ell) | (rot(Bp, g, ell, mask) << (2 * ell))
        rows.append(r0)
    for g in range(ell):
        r1 = rot(One, g, ell, mask) | (rot(Q, g, ell, mask) << ell) | (rot(Bp, g, ell, mask) << (3 * ell))
        rows.append(r1)
    return rows


def hz_rows(ell: int, p, q, b) -> list[int]:
    mask = (1 << ell) - 1
    # Z-check (j,g): L(j, g - s) for s in b ; R(i, h) with A†[j,i] = conj(A[i,j]): (g-h) ∈ -A[i,j] -> h = g + s, s∈A[i,j]
    Bm = pattern(b, ell, True); Pp = pattern(p, ell, False); Qp = pattern(q, ell, False); One = 1
    rows = []
    for g in range(ell):  # j = 0: L0 via b, R0 via conj(p) (i=0), R1 via conj(1) (i=1, A[1,0]=1)
        rows.append(rot(Bm, g, ell, mask) | (rot(Pp, g, ell, mask) << (2 * ell)) | (rot(One, g, ell, mask) << (3 * ell)))
    for g in range(ell):  # j = 1: L1 via b, R0 via conj(A[0,1]=1), R1 via conj(q)
        rows.append((rot(Bm, g, ell, mask) << ell) | (rot(One, g, ell, mask) << (2 * ell)) | (rot(Qp, g, ell, mask) << (3 * ell)))
    return rows


def k_of(ell: int, p, q, b) -> int:
    n = 4 * ell
    rX = f2.rank(hx_rows(ell, p, q, b))
    return n - 2 * rX   # A symmetric ⇒ rank H_X = rank H_Z


def to_numpy(ell: int, p, q, b):
    import numpy as np
    from .lp import CSSCode
    n = 4 * ell
    HX = f2.matrix_from_rows(hx_rows(ell, p, q, b), n)
    HZ = f2.matrix_from_rows(hz_rows(ell, p, q, b), n)
    return CSSCode(HX, HZ, meta={"ell": ell, "p": list(p), "q": list(q), "b": list(b), "a": 2})
