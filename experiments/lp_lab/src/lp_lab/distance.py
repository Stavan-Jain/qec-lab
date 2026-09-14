"""Distance tools for CSS codes: exact k, randomized information-set upper
bounds (QDistRnd-style), and a wrapper around bb_lab's exact SAT ladder.
"""
from __future__ import annotations

import random
import time
from dataclasses import dataclass

import numpy as np

from . import f2
from .lp import CSSCode


def code_k(code: CSSCode) -> int:
    rX = f2.rank(f2.rows_from_matrix(code.H_X))
    rZ = f2.rank(f2.rows_from_matrix(code.H_Z))
    return code.num_qubits - rX - rZ


def logical_basis(H_stab_same: np.ndarray, H_stab_other: np.ndarray, n: int) -> list[int]:
    """Representatives of ker(H_other) / rowspan(H_same) as ints.

    For X-logicals: pass H_same = H_X (X stabilizers), H_other = H_Z.
    """
    ker = f2.nullspace(f2.rows_from_matrix(H_stab_other), n)
    piv = f2.echelon(f2.rows_from_matrix(H_stab_same))
    piv = dict(piv)
    reps = []
    acc = dict(piv)  # grows with accepted logicals
    for v in ker:
        r = f2.reduce(v, acc)
        if r:
            reps.append(v)
            # insert reduced r into acc echelon
            lb = r & -r
            acc[lb] = r
    return reps


@dataclass
class ISDResult:
    weight: int
    witness: int          # bit-packed vector
    trials: int
    seconds: float


def isd_min_weight(H_other: np.ndarray, L_dual: list[int], n: int,
                   trials: int = 200, seed: int | None = None,
                   pairs: bool = True, stop_at: int | None = None) -> ISDResult:
    """Randomized information-set search for a low-weight vector in
    ker(H_other) that is *nontrivial*, i.e. has odd overlap with some
    element of `L_dual` (a basis of the dual logicals).

    Each trial: random column order, RREF of a kernel basis in that order
    (pivots = information set), inspect rows and (optionally) pairwise sums.
    Returns the best weight seen (an upper bound on the distance).
    """
    rng = random.Random(seed)
    t0 = time.time()
    K = f2.nullspace(f2.rows_from_matrix(H_other), n)
    if not K:
        raise ValueError("trivial kernel")
    Ld = list(L_dual)

    def nontrivial(v: int) -> bool:
        for l in Ld:
            if (v & l).bit_count() & 1:
                return True
        return False

    best_w = n + 1
    best_v = 0
    cols = list(range(n))
    for t in range(trials):
        rng.shuffle(cols)
        # RREF with pivot preference following the permuted column order
        rows = list(K)
        basis: list[int] = []
        pivots: list[int] = []
        for c in cols:
            bit = 1 << c
            # find a row (not yet used) with this bit
            idx = -1
            for i, r in enumerate(rows):
                if r & bit:
                    idx = i
                    break
            if idx < 0:
                continue
            prow = rows.pop(idx)
            rows = [r ^ prow if (r & bit) else r for r in rows]
            basis = [b ^ prow if (b & bit) else b for b in basis]
            basis.append(prow)
            pivots.append(bit)
            if not rows:
                break
        # inspect rows
        for v in basis:
            w = v.bit_count()
            if w < best_w and nontrivial(v):
                best_w, best_v = w, v
        if pairs:
            m = len(basis)
            for i in range(m):
                bi = basis[i]
                wi = bi.bit_count()
                if wi - 1 >= best_w + 20:  # cheap prune: pair sums rarely help if both heavy
                    pass
                for j in range(i + 1, m):
                    v = bi ^ basis[j]
                    w = v.bit_count()
                    if w < best_w and nontrivial(v):
                        best_w, best_v = w, v
        if stop_at is not None and best_w <= stop_at:
            return ISDResult(best_w, best_v, t + 1, time.time() - t0)
    return ISDResult(best_w, best_v, trials, time.time() - t0)


def distance_upper_bounds(code: CSSCode, trials: int = 200, seed: int | None = None,
                          stop_at: int | None = None) -> tuple[ISDResult, ISDResult]:
    """(d_X upper bound, d_Z upper bound) via ISD."""
    n = code.num_qubits
    # X-logicals live in ker(H_Z) \ rs(H_X); nontrivial iff odd overlap with a Z-logical rep
    LZ = logical_basis(code.H_Z, code.H_X, n)   # Z-logical reps: ker(H_X) / rs(H_Z)
    LX = logical_basis(code.H_X, code.H_Z, n)   # X-logical reps: ker(H_Z) / rs(H_X)
    rX = isd_min_weight(code.H_Z, LZ, n, trials=trials, seed=seed, stop_at=stop_at)
    rZ = isd_min_weight(code.H_X, LX, n, trials=trials, seed=None if seed is None else seed + 1,
                        stop_at=stop_at)
    return rX, rZ


def exact_distance(code: CSSCode, which: str = "X", lower: int = 1, upper: int | None = None,
                   verbose: bool = False, code_id: str = "code"):
    """Exact d_X (or d_Z) via bb_lab's CryptoMiniSat ladder."""
    from bb_lab.sat_distance import x_distance
    if which == "X":
        c = code
    else:
        c = CSSCode(code.H_Z, code.H_X)
    return x_distance(c, weight_lower_bound=lower, weight_upper_bound=upper,
                      verbose=verbose, code_id=code_id)
