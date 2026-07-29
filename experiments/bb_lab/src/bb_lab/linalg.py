"""F₂ linear algebra shared by codeparams.py and sat_distance.py.

Dense uint8 implementations sized for v0 (n ≤ 288). v1 will switch
to bitpacked uint64 columns when the L1 enumeration phase comes online.
"""

from __future__ import annotations

import numpy as np


def rref_f2(M: np.ndarray) -> tuple[np.ndarray, list[int]]:
    """Return (reduced-row-echelon form, list of pivot columns)."""
    A = (M & 1).astype(np.uint8, copy=True)
    rows, cols = A.shape
    pivots: list[int] = []
    pivot_row = 0
    for col in range(cols):
        if pivot_row >= rows:
            break
        sub = A[pivot_row:, col]
        nonzero = np.flatnonzero(sub)
        if nonzero.size == 0:
            continue
        r = pivot_row + int(nonzero[0])
        if r != pivot_row:
            A[[pivot_row, r]] = A[[r, pivot_row]]
        mask = A[:, col] == 1
        mask[pivot_row] = False
        if mask.any():
            A[mask] ^= A[pivot_row]
        pivots.append(col)
        pivot_row += 1
    return A, pivots


def rank_f2(M: np.ndarray) -> int:
    _, pivots = rref_f2(M)
    return len(pivots)


def nullspace_f2(M: np.ndarray) -> np.ndarray:
    """Return a basis for ker(M) over F₂ as rows.

    Shape `(n − rank(M), n)`. Uses the standard transpose-augment-reduce
    trick: reduce ``[Mᵀ | I_n]`` in column-by-column row-echelon, then
    the rows whose left half is zero have their right half as a kernel
    basis vector.
    """
    n = M.shape[1]
    A = np.concatenate(
        [(M & 1).T.astype(np.uint8), np.eye(n, dtype=np.uint8)], axis=1
    )
    rows = A.shape[0]
    left_cols = M.shape[0]
    pivot_row = 0
    for col in range(left_cols):
        if pivot_row >= rows:
            break
        sub = A[pivot_row:, col]
        nonzero = np.flatnonzero(sub)
        if nonzero.size == 0:
            continue
        r = pivot_row + int(nonzero[0])
        if r != pivot_row:
            A[[pivot_row, r]] = A[[r, pivot_row]]
        mask = A[:, col] == 1
        mask[pivot_row] = False
        if mask.any():
            A[mask] ^= A[pivot_row]
        pivot_row += 1
    null_rows = [
        A[r, left_cols:].copy()
        for r in range(pivot_row, rows)
        if not A[r, :left_cols].any()
    ]
    if not null_rows:
        return np.zeros((0, n), dtype=np.uint8)
    return np.stack(null_rows).astype(np.uint8)


def quotient_complement_basis(
    base: np.ndarray, extension: np.ndarray
) -> np.ndarray:
    """Return rows of `extension` that are linearly independent modulo
    `rowspan(base)`.

    Used for CSS logical representatives: ``base = H_Z``,
    ``extension = nullspace(H_X)`` returns a basis of
    ``ker(H_X) / rowspan(H_Z)``  — the k logical-Z operators.

    Returned rows are the *original* extension rows (not their reduced
    forms), so they still anticommute correctly with X-logicals as
    `<·, v>` witnesses.

    Implementation note (2026-07-29): the original joint-reduction with
    origin tracking could *misattribute* a pivot to an extension row
    when a swap displaced a base row below the extension block (first
    triggered by a Z7×Z12 corpus instance, returning k+1 rows). This
    version is order-robust: reduce each extension row against a fixed
    rref of `base`, then against the previously accepted (reduced)
    rows; accept iff nonzero. Count is exactly
    rank(base ∪ extension) − rank(base).
    """
    base = (base & 1).astype(np.uint8)
    extension = (extension & 1).astype(np.uint8)
    if base.size == 0:
        return extension.copy()

    base_rref, _ = rref_f2(base)
    base_rows = [r for r in base_rref if r.any()]
    base_pivots = [int(np.flatnonzero(r)[0]) for r in base_rows]

    accepted: list[np.ndarray] = []       # reduced forms
    accepted_pivots: list[int] = []
    chosen_ext: list[int] = []
    for i in range(extension.shape[0]):
        r = extension[i].copy()
        for b, p in zip(base_rows, base_pivots):
            if r[p]:
                r ^= b
        for b, p in zip(accepted, accepted_pivots):
            if r[p]:
                r ^= b
        if r.any():
            chosen_ext.append(i)
            accepted.append(r)
            accepted_pivots.append(int(np.flatnonzero(r)[0]))
    if not chosen_ext:
        return np.zeros((0, extension.shape[1]), dtype=np.uint8)
    return extension[chosen_ext].copy()
