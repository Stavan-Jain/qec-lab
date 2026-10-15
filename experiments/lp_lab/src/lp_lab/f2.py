"""Bit-packed F₂ linear algebra on Python ints (rows = ints, bit i = column i).

Sized for n ≤ ~1000. Pure Python big-int XORs are fast enough for
rank/kernel computations of a few hundred rows and for the randomized
information-set distance heuristic in `lp_lab.distance`.
"""
from __future__ import annotations

import random
from typing import Iterable

import numpy as np


def rows_from_matrix(M: np.ndarray) -> list[int]:
    """Convert a 0/1 numpy matrix (rows × cols) to a list of row ints."""
    M = (np.asarray(M) & 1).astype(np.uint8)
    out = []
    for r in M:
        # little-endian: bit j <-> column j
        v = 0
        for j in np.flatnonzero(r):
            v |= 1 << int(j)
        out.append(v)
    return out


def matrix_from_rows(rows: Iterable[int], ncols: int) -> np.ndarray:
    rows = list(rows)
    M = np.zeros((len(rows), ncols), dtype=np.uint8)
    for i, r in enumerate(rows):
        j = 0
        while r:
            if r & 1:
                M[i, j] = 1
            r >>= 1
            j += 1
    return M


def popcount(x: int) -> int:
    return x.bit_count()


def rank(rows: Iterable[int]) -> int:
    """Rank over F₂ via elimination on the lowest set bit."""
    basis: list[int] = []  # kept in echelon form keyed by lowest bit
    piv: dict[int, int] = {}
    for r in rows:
        while r:
            lb = r & -r
            b = piv.get(lb)
            if b is None:
                piv[lb] = r
                break
            r ^= b
    return len(piv)


def echelon(rows: Iterable[int]) -> dict[int, int]:
    """Return {pivot_bit: row} in (non-reduced) echelon form."""
    piv: dict[int, int] = {}
    for r in rows:
        while r:
            lb = r & -r
            b = piv.get(lb)
            if b is None:
                piv[lb] = r
                break
            r ^= b
    return piv


def reduce(r: int, piv: dict[int, int]) -> int:
    """Reduce r against an echelon dict (lowest-bit pivots)."""
    while r:
        lb = r & -r
        b = piv.get(lb)
        if b is None:
            # lowest bit not a pivot: cannot be eliminated; move to next bit
            # by clearing it temporarily — instead do full reduction over all pivots
            break
        r ^= b
    # full reduction: eliminate any pivot bits present
    changed = True
    while changed:
        changed = False
        for lb, b in piv.items():
            if r & lb:
                r ^= b
                changed = True
    return r


def in_span(r: int, piv: dict[int, int]) -> bool:
    return reduce(r, piv) == 0


def nullspace(rows: list[int], ncols: int) -> list[int]:
    """Basis of {v : <row, v> = 0 for all rows} as ints over `ncols` bits.

    Standard: RREF the matrix, free columns give basis vectors.
    """
    # full RREF with pivot -> row, tracking pivot columns
    piv: dict[int, int] = {}
    for r in rows:
        # reduce r by existing pivots (any bit)
        for lb, b in piv.items():
            if r & lb:
                r ^= b
        if r:
            lb = r & -r
            # eliminate lb from existing pivot rows
            for k in list(piv.keys()):
                if piv[k] & lb:
                    piv[k] ^= r
            piv[lb] = r
    pivot_cols = {lb.bit_length() - 1 for lb in piv}
    free_cols = [c for c in range(ncols) if c not in pivot_cols]
    basis = []
    for f in free_cols:
        v = 1 << f
        for lb, b in piv.items():
            if b & (1 << f):
                v |= lb  # pivot column of this row gets a 1
        basis.append(v)
    return basis


def dot(a: int, b: int) -> int:
    return (a & b).bit_count() & 1
