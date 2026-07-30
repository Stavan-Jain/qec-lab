"""The (w, b) normal form and twist-syndrome floors.

For a free axis deck σ (|⟨σ⟩| = 2), every cover vector decomposes
UNIQUELY and linearly as

    v = p*(b) + L₀(w),   w = p₊v (odd-fiber pattern),
                          b = v|sheet-1 (in base coordinates),

with cost |v| = |w| + 2·|b ∖ w|. The cycle condition splits:

    H_Z v = 0   ⟺   H̄_Z w = 0  ∧  H̄_Z b = T(w)

where T (the *twist*) is the sheet-0-lift composed with the cover
checks, read off either check-row lift — the two row choices agree on
ker H̄_Z (asserted). Consequences, each asserted numerically here:

- moving-sector cost F2(λ) obeys the *window inequality*

      |v| ≥ max( |w| , 2·K(T(w)) − |w| )

  where K(s) = min{|b| : H̄_Z b = s} is a BASE-code syndrome-coset
  minimum — light pushforwards force heavy twist corrections exactly
  when their twist syndrome is expensive (the A15 window/dangerous
  mechanism in exact linear-algebra form);

- stratified floors: if every base cycle w in class λ with
  |w| ≤ W₀ has max(|w|, 2K(T(w)) − |w|) ≥ m, then
  F2(λ) ≥ min(m, W₀ + 2) (weights are even under the verified parity
  premise).

For bb_288 along (0, 6) the base is literally the gross code, so K is
a gross syndrome-decoding minimum.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from pysat.card import CardEnc, EncType
from pysat.formula import IDPool

from .checks import CheckMatrices
from .descent_sat import DescentData
from .linalg import nullspace_f2

try:
    import pycryptosat
except ImportError:  # pragma: no cover
    pycryptosat = None


@dataclass(frozen=True, slots=True)
class TwistData:
    """The verified (w, b) normal form for one deck."""

    L0: np.ndarray        # (n, n̄) canonical sheet-0 section
    E1: np.ndarray        # (n̄, n) extract b = v|sheet-1
    T: np.ndarray         # (Nb, n̄) twist: cycles satisfy H̄_Z b = T w
    Hb_Z: np.ndarray      # base checks (Nb, n̄)


def compute_twist(
    checks: CheckMatrices, dd: DescentData
) -> TwistData:
    """Build and verify the normal form. σ must be an axis deck with
    the second coordinate halved or the first (matching dd.sigma)."""
    G = checks.group
    N = G.cardinality
    n = checks.num_qubits
    nb = dd.S.shape[0]
    Nb = nb // 2
    sigma = dd.sigma
    ell, m = G.orders
    elems = list(G)

    # Sheet of a cover position: exponent along the halved axis below
    # half ⟹ sheet 0.
    if sigma[0]:
        half, axis = ell // 2, 0
    else:
        half, axis = m // 2, 1

    Gb = dd.base_checks.group

    def red(g):
        return (
            (g[0] % half, g[1]) if axis == 0 else (g[0], g[1] % half)
        )

    L0 = np.zeros((n, nb), dtype=np.uint8)
    E1 = np.zeros((nb, n), dtype=np.uint8)
    for i, g in enumerate(elems):
        j = Gb.index(red(g))
        for beta in (0, 1):
            if g[axis] < half:
                L0[i + beta * N, j + beta * Nb] = 1
            else:
                E1[j + beta * Nb, i + beta * N] = 1

    # Verification: bijection + cost formula on random vectors.
    rng = np.random.default_rng(11)
    for _ in range(20):
        v = rng.integers(0, 2, size=n).astype(np.uint8)
        w = (dd.S @ v) % 2
        b = (E1 @ v) % 2
        v2 = ((dd.P_lift @ b) + (L0 @ w)) % 2
        assert np.array_equal(v, v2), "(w,b) is not a bijection"
        cost = int(w.sum()) + 2 * int((b & (1 - w)).sum())
        assert cost == int(v.sum()), "cost formula broken"

    # Twist rows: pick, per base check row, the lift with axis-exponent
    # below half; the other lift must agree on ker H̄_Z (asserted).
    HL0 = (checks.H_Z @ L0) % 2
    row_of = {}
    for r, g in enumerate(elems):
        row_of[g] = r
    T0 = np.zeros((Nb, nb), dtype=np.uint8)
    T1 = np.zeros((Nb, nb), dtype=np.uint8)
    belems = list(Gb)
    for rb, gb in enumerate(belems):
        g0 = list(gb)
        g1 = list(gb)
        g1[axis] += half
        T0[rb] = HL0[row_of[tuple(g0)]]
        T1[rb] = HL0[row_of[tuple(g1)]]
    ker_b = nullspace_f2(dd.base_checks.H_Z)
    assert not (((T0 + T1) % 2) @ ker_b.T % 2).any(), (
        "twist rows disagree on base cycles; normal form invalid"
    )

    # Cycle-condition split, asserted on random cover cycles.
    ker_c = nullspace_f2(checks.H_Z)
    Hb_Z = dd.base_checks.H_Z
    for _ in range(20):
        coeff = rng.integers(0, 2, size=ker_c.shape[0]).astype(np.uint8)
        v = (coeff @ ker_c) % 2
        w = (dd.S @ v) % 2
        b = (E1 @ v) % 2
        assert not ((Hb_Z @ w) % 2).any()
        assert np.array_equal((Hb_Z @ b) % 2, (T0 @ w) % 2), (
            "twist equation fails on a cover cycle"
        )
    # Falsify-first: a random non-cycle must violate one of the halves.
    for _ in range(10):
        v = rng.integers(0, 2, size=n).astype(np.uint8)
        if not ((checks.H_Z @ v) % 2).any():
            continue
        w = (dd.S @ v) % 2
        b = (E1 @ v) % 2
        ok = (not ((Hb_Z @ w) % 2).any()) and np.array_equal(
            (Hb_Z @ b) % 2, (T0 @ w) % 2
        )
        assert not ok, "non-cycle passed the split condition"

    return TwistData(L0=L0, E1=E1, T=T0, Hb_Z=Hb_Z)


# ---------------------------------------------------------------------------
# K(s): base syndrome-coset minima (budgeted, certified)


def syndrome_min_budgeted(
    Hb_Z: np.ndarray,
    syndrome: np.ndarray,
    cap: int,
    confl_budget: int = 200_000,
) -> tuple[int, int | None]:
    """(floor, exact): floor = certified K(s) ≥ floor; exact = K(s)
    when a model was found (then floor == exact). Ascending budgeted
    CMS on the base code: solve H̄_Z b = s with |b| ≤ t."""
    nb = Hb_Z.shape[1]
    syndrome = syndrome.astype(np.uint8)
    if not syndrome.any():
        return 0, 0

    def _try(t: int):
        pool = IDPool()
        bv = [pool.id() for _ in range(nb)]
        solver = pycryptosat.Solver(confl_limit=confl_budget)
        for r, row in enumerate(Hb_Z):
            idx = np.flatnonzero(row)
            if idx.size:
                solver.add_xor_clause(
                    [bv[i] for i in idx], bool(syndrome[r])
                )
            elif syndrome[r]:
                return False           # 0 = 1: no solution at all
        card = CardEnc.atmost(
            lits=bv, bound=t, vpool=pool, encoding=EncType.seqcounter
        )
        for cl in card.clauses:
            solver.add_clause(cl)
        sat, model = solver.solve()
        return sat

    floor = 1
    for t in range(1, cap + 1):
        res = _try(t)
        if res is True:
            return t, t
        if res is not False:
            break                      # budget
        floor = t + 1
    return floor, None


def window_floor(w_weight: int, K_floor: int) -> int:
    """The certified per-w cost bound max(|w|, 2K − |w|)."""
    return max(w_weight, 2 * K_floor - w_weight)


# ---------------------------------------------------------------------------
# Banded enumeration + masked-K certificates (the stratified floor sweep)


def enum_band(
    td: TwistData,
    Lb_Z: np.ndarray,
    kb: int,
    lam: int,
    wmin: int,
    wmax: int,
    limit: int = 3000,
    confl_budget: int = 2_000_000,
    wall_budget: float | None = None,
) -> tuple[list[np.ndarray], bool]:
    """AllSAT: base cycles w with μ̄(w) = λ and wmin ≤ |w| ≤ wmax,
    enumerated with blocking clauses. Returns (members, complete).
    `complete=True` REQUIRES the final call to return UNSAT within the
    conflict budget — a budget-out or limit/wall stop returns False and
    certifies nothing about the band."""
    import time as _time

    nb = td.Hb_Z.shape[1]
    lam_bits = np.array(
        [(lam >> j) & 1 for j in range(kb)], dtype=np.uint8
    )
    pool = IDPool()
    wv = [pool.id() for _ in range(nb)]
    s = pycryptosat.Solver(confl_limit=confl_budget)
    for row in td.Hb_Z:
        idx = np.flatnonzero(row)
        if idx.size:
            s.add_xor_clause([wv[i] for i in idx], False)
    for j in range(kb):
        idx = np.flatnonzero(Lb_Z[j])
        s.add_xor_clause([wv[i] for i in idx], bool(lam_bits[j]))
    if lam == 0:
        s.add_clause(wv)
    for cl in CardEnc.atmost(
        lits=wv, bound=wmax, vpool=pool, encoding=EncType.seqcounter
    ).clauses:
        s.add_clause(cl)
    if wmin > 1:
        for cl in CardEnc.atleast(
            lits=wv, bound=wmin, vpool=pool, encoding=EncType.seqcounter
        ).clauses:
            s.add_clause(cl)
    out: list[np.ndarray] = []
    t0 = _time.perf_counter()
    while True:
        if len(out) >= limit:
            return out, False
        if wall_budget and _time.perf_counter() - t0 > wall_budget:
            return out, False
        sat, model = s.solve()
        if sat is False:
            return out, True
        if sat is not True:               # conflict budget exhausted
            return out, False
        w = np.array([1 if model[v] else 0 for v in wv], dtype=np.uint8)
        out.append(w)
        s.add_clause([-wv[i] if w[i] else wv[i] for i in range(nb)])


def _consistent(Asub: np.ndarray, s: np.ndarray) -> bool:
    """Is A_sub·x = s solvable over F₂?"""
    aug = np.concatenate(
        [Asub % 2, s.reshape(-1, 1) % 2], axis=1
    ).astype(np.uint8)
    M = aug.copy()
    r = 0
    for c in range(Asub.shape[1]):
        nz = np.flatnonzero(M[r:, c])
        if nz.size == 0:
            continue
        p = r + int(nz[0])
        M[[r, p]] = M[[p, r]]
        for q in range(M.shape[0]):
            if q != r and M[q, c]:
                M[q] ^= M[r]
        r += 1
        if r == M.shape[0]:
            break
    return not any(
        M[i, :-1].max(initial=0) == 0 and M[i, -1]
        for i in range(M.shape[0])
    )


def masked_K_at_least(
    td: TwistData, w: np.ndarray, k: int
) -> bool:
    """Certified: every b with H̄_Z b = T(w) has |b ∖ supp(w)| ≥ k,
    by exhausting all outside-support patterns of size < k with
    F₂ linear solves (no SAT). k ≤ 3 is the practical range
    (sizes 0, 1, 2 → 1 + n̄ + n̄²/2 solves)."""
    supp = np.flatnonzero(w)
    outside = np.setdiff1d(np.arange(td.Hb_Z.shape[1]), supp)
    s0 = (td.T @ w) % 2
    import itertools as _it

    for size in range(k):
        for extra in _it.combinations(outside, size):
            cols = (
                np.concatenate([supp, np.array(extra, dtype=int)])
                if extra else supp
            )
            if _consistent(td.Hb_Z[:, cols], s0):
                return False
    return True


def masked_K_value(
    td: TwistData, w: np.ndarray, kmax: int = 3
) -> int:
    """min |b ∖ supp(w)| over the twist coset {b : H̄_Z b = T(w)},
    capped at kmax (returns kmax when ≥ kmax). Exact linear algebra:
    eliminate the supp(w) columns once; solvability with an
    outside-support pattern E ⟺ the residual of T(w) equals the sum
    of the residuals of E's columns — checked for |E| = 0, 1 by direct
    residual comparison and |E| = 2 by a hash of column residuals.
    Cross-checked against `masked_K_at_least` in tests."""
    H = td.Hb_Z
    supp = np.flatnonzero(w)
    outside = np.setdiff1d(np.arange(H.shape[1]), supp)
    aug = np.concatenate(
        [H[:, supp], H[:, outside], ((td.T @ w) % 2).reshape(-1, 1)],
        axis=1,
    ).astype(np.uint8)
    ns = len(supp)
    r = 0
    for c in range(ns):
        nz = np.flatnonzero(aug[r:, c])
        if nz.size == 0:
            continue
        p = r + int(nz[0])
        aug[[r, p]] = aug[[p, r]]
        mask = aug[:, c] == 1
        mask[r] = False
        if mask.any():
            aug[mask] ^= aug[r]
        r += 1
        if r == aug.shape[0]:
            break
    res = aug[r:]                      # rows with supp-part eliminated
    if res.shape[0] == 0:
        return 0                       # always solvable with b ⊆ supp
    s_res = res[:, -1]
    cols = res[:, ns:-1]               # residuals of outside columns
    if not s_res.any():
        return 0
    if kmax <= 1:
        return 1
    col_bytes = {cols[:, i].tobytes(): i for i in range(cols.shape[1])}
    if s_res.tobytes() in col_bytes:
        return 1
    if kmax <= 2:
        return 2
    for i in range(cols.shape[1]):
        if ((s_res ^ cols[:, i])).tobytes() in col_bytes:
            return 2
    return kmax if kmax <= 3 else 3


def masked_K_sat_budgeted(
    td: TwistData,
    w: np.ndarray,
    need: int,
    confl_budget: int = 150_000,
) -> tuple[int, np.ndarray | None]:
    """Certified floor on min |b ∖ supp(w)| over the twist coset, by
    ascending budgeted CMS with the cardinality restricted to the
    outside-support bits. Returns (certified_j, witness_b): UNSAT at
    j certifies ≥ j+1 (ascent continues to `need`); a SAT stops with
    the exact value and its witness; budget-out stops with the floor
    proven so far. Complements the exact hash path
    (`masked_K_value`, need ≤ 3) for larger targets."""
    nb = td.Hb_Z.shape[1]
    supp = set(np.flatnonzero(w).tolist())
    outside = [i for i in range(nb) if i not in supp]
    syn = (td.T @ w) % 2
    certified = 0
    for j in range(0, need):
        pool = IDPool()
        bv = [pool.id() for _ in range(nb)]
        s = pycryptosat.Solver(confl_limit=confl_budget)
        consistent = True
        for r, row in enumerate(td.Hb_Z):
            idx = np.flatnonzero(row)
            if idx.size:
                s.add_xor_clause([bv[i] for i in idx], bool(syn[r]))
            elif syn[r]:
                consistent = False
                break
        if not consistent:
            return need, None          # coset empty: vacuous
        card = CardEnc.atmost(
            lits=[bv[i] for i in outside], bound=j, vpool=pool,
            encoding=EncType.seqcounter,
        )
        for cl in card.clauses:
            s.add_clause(cl)
        sat, model = s.solve()
        if sat is True:
            b = np.array(
                [1 if model[v] else 0 for v in bv], dtype=np.uint8
            )
            return j, b
        if sat is not False:
            return certified, None     # budget
        certified = j + 1
    return certified, None
