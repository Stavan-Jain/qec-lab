"""Descent-strengthened shard SAT: the quotient/pushforward layer.

`shard_distance.py` strengthens the monolithic instance with facts from
the *translation action* (classes, transport, anchors, coset parity).
This module adds the next analytic layer from the repo's cover/descent
theory (A10–A15): the **quotient structure** of the code along a free
deck translation σ of order 2.

Let σ ∈ G with 2σ = 0, σ ≠ 0, and let p : G → Ḡ = G/⟨σ⟩ be the
quotient. p acts on qubits blockwise and induces

    p₊ : F₂^{2|G|} → F₂^{2|Ḡ|},   (p₊ v)(q̄, β) = Σ_{q ∈ p⁻¹(q̄)} v(q, β),

the chain-level pushforward. Three elementary facts — each *asserted
numerically at construction*, none trusted from algebra — power the
decomposition:

1.  **Intertwining.** R·H_Z = H̄_Z·S where S is the matrix of p₊ and R
    the row (check) pushforward. Hence H_Z v = 0 ⟹ H̄_Z (p₊ v) = 0:
    the pushforward of a cover cycle is a base cycle, so the *base
    code's parity system on auxiliary pushforward bits is an implied
    constraint* — pure strengthening, propagation for free.

2.  **Sector dichotomy.** For |⟨σ⟩| = 2, p₊ v = 0 ⟺ v is σ-invariant
    ⟺ v = p* v̄ for a base-sized v̄, and then |v| = 2|v̄| (σ acts
    freely). So each shard splits into
      (a) the invariant sector — an instance over n/2 variables at
          weight ⌊w/2⌋ with class pins μ(p* v̄) = c, and
      (b) the moving sector — the cover instance plus p₊ v ≠ 0 and the
          implied base rows on the pushforward bits.
    Refuting both refutes the shard; a model of either is a genuine
    member (models of (a) are expanded through p* and re-verified).

3.  **Class transport (the A14 face).** μ̄(p₊ v) = Λ·μ(v) for a
    computable Λ (well-defined because p₊ maps X-stabilizers into base
    X-stabilizers — asserted). Sector (b) therefore gets k̄ *pinned*
    parity constraints on the pushforward bits: the solver knows which
    base homology class the pushforward must land in. Sector (a) gets
    the complementary analytic filter: c is reachable only if
    c ∈ Pμ·ker H̄_Z (and only with even parity, since |v| = 2|v̄|) —
    reps failing the rank membership are refuted *with no SAT call*.
    By the A14 LES (`ker p₁ = im τ₁`, dim k/2) roughly half of each
    class space dies this way.

The decomposition composes with everything shard_distance already
does: the same parity-homogeneous coarse groups, the same Stab(c)
anchor transversals in sector (b) (translations commute with σ, so
they preserve both sectors and the pin values — asserted via the
Λ·ρ(t) = ρ̄(t̄)·Λ commutation), and the same coset-parity tightening
of refutation bounds.

Scope (v1): axis decks σ ∈ {(ℓ/2, 0), (0, m/2)} on Z_ℓ × Z_m (the
literal-lift direction the descent theory targets); CMS backend only;
X direction (BB transpose symmetry gives d_Z = d_X). Diagonal decks
and odd-order quotients (v = (1+σ)u substitution) are documented
upgrades, not implemented.
"""

from __future__ import annotations

import time
from dataclasses import dataclass, field

import numpy as np

from pysat.card import CardEnc, EncType
from pysat.formula import IDPool

from .checks import CheckMatrices, bb_check_matrices
from .group import ZmZn
from .linalg import nullspace_f2
from .poly import Poly
from .sat_distance import find_logical_z
from .shard_distance import (
    ClassAction,
    _bits_to_int,
    _greedy_upper_bound,
    _int_to_bits,
    _monolith_probe,
    _verify_witness,
    compute_class_action,
)

try:
    import pycryptosat
    _HAVE_CMS = True
except ImportError:  # pragma: no cover
    _HAVE_CMS = False


# ---------------------------------------------------------------------------
# F₂ helpers


def _rank_f2(M: np.ndarray) -> int:
    A = (M % 2).astype(np.uint8).copy()
    rank = 0
    rows, cols = A.shape
    for c in range(cols):
        nz = np.flatnonzero(A[rank:, c])
        if nz.size == 0:
            continue
        r = rank + int(nz[0])
        if r != rank:
            A[[rank, r]] = A[[r, rank]]
        mask = A[:, c] == 1
        mask[rank] = False
        if mask.any():
            A[mask] ^= A[rank]
        rank += 1
        if rank == rows:
            break
    return rank


def _in_span(M: np.ndarray, v: np.ndarray) -> bool:
    """Is v in the row span of M (over F₂)?"""
    if M.shape[0] == 0:
        return not v.any()
    return _rank_f2(np.vstack([M, v])) == _rank_f2(M)


# ---------------------------------------------------------------------------
# Descent data: the numerically-verified quotient package


@dataclass(frozen=True, slots=True)
class DescentData:
    """Quotient structure along one free order-2 deck σ, verified."""

    sigma: tuple[int, ...]
    base_checks: CheckMatrices          # the descended BB code on G/⟨σ⟩
    S: np.ndarray                       # (n̄, n) pushforward p₊ on qubit vectors
    P_lift: np.ndarray                  # (n, n̄) pullback p*
    Lb_Z: np.ndarray                    # (k̄, n̄) base pairing basis
    kb: int
    Pmu: np.ndarray                     # (k, n̄): μ(p* v̄) = Pmu·v̄
    Lam: np.ndarray                     # (k, k̄): μ̄(p₊ v) = μ(v)·Lam
    W_a: np.ndarray                     # basis of Pμ(ker H̄_Z): (a)-reachable classes
    rho_b_rows: tuple[np.ndarray, ...] = field(repr=False)  # base ρ̄ per t̄


def axis_decks(checks: CheckMatrices) -> list[tuple[int, ...]]:
    """The axis involutions of Z_ℓ × Z_m (order-2 deck candidates)."""
    ell, m = checks.group.orders
    out = []
    if ell % 2 == 0:
        out.append((ell // 2, 0))
    if m % 2 == 0:
        out.append((0, m // 2))
    return out


def compute_descent(
    checks: CheckMatrices,
    action: ClassAction,
    sigma: tuple[int, ...],
    A: Poly | None = None,
    B: Poly | None = None,
) -> DescentData:
    """Build and *verify* the quotient package for one axis deck σ.

    Every fact the SAT layer uses is asserted here on the actual
    matrices; an AssertionError means the input does not descend along
    σ and no decomposition is returned.
    """
    G = checks.group
    N = G.cardinality
    n = checks.num_qubits
    ell, m = G.orders
    assert G.add(sigma, sigma) == (0, 0) and sigma != (0, 0), (
        "σ must be a free involution of G"
    )
    if sigma == (ell // 2, 0) and ell % 2 == 0:
        Gb = ZmZn(ell // 2, m)
        red = lambda g: (g[0] % (ell // 2), g[1])          # noqa: E731
    elif sigma == (0, m // 2) and m % 2 == 0:
        Gb = ZmZn(ell, m // 2)
        red = lambda g: (g[0], g[1] % (m // 2))            # noqa: E731
    else:
        raise ValueError(f"σ = {sigma} is not an axis deck of Z{ell}×Z{m}")
    Nb = Gb.cardinality
    nb = 2 * Nb
    assert 2 * Nb == N

    # Position maps: fiber_of[i] = base index of cover position i.
    elems = list(G)
    fiber_of = np.array([Gb.index(red(g)) for g in elems], dtype=np.int64)
    # Each base position has exactly two preimages (σ acts freely).
    lift_pairs = np.full((Nb, 2), -1, dtype=np.int64)
    for i, b in enumerate(fiber_of):
        j = 0 if lift_pairs[b, 0] == -1 else 1
        lift_pairs[b, j] = i
    assert (lift_pairs >= 0).all(), "quotient fibers are not all size 2"

    # Pushforward S (both blocks) and pullback P_lift as explicit matrices.
    S = np.zeros((nb, n), dtype=np.uint8)
    P_lift = np.zeros((n, nb), dtype=np.uint8)
    for b in range(Nb):
        for i in lift_pairs[b]:
            for beta in (0, 1):
                S[b + beta * Nb, i + beta * N] = 1
                P_lift[i + beta * N, b + beta * Nb] = 1

    # The descended code, built from the image polynomials (collisions
    # XOR out mod 2), through the same trusted constructor as the cover.
    if A is None or B is None:
        raise ValueError("compute_descent needs the A, B polynomials")
    Ab = Poly(frozenset(
        h for h in {red(g) for g in A.support}
        if sum(1 for g in A.support if red(g) == h) % 2
    ), Gb)
    Bb = Poly(frozenset(
        h for h in {red(g) for g in B.support}
        if sum(1 for g in B.support if red(g) == h) % 2
    ), Gb)
    base = bb_check_matrices(Ab, Bb)

    # (1) Intertwining: R·H = H̄·S for both H_X and H_Z, where R is the
    # row pushforward (rows of the cover check at g and g+σ both map to
    # the base row at p(g); R has a single 1 per column).
    R = np.zeros((Nb, N), dtype=np.uint8)
    for i in range(N):
        R[fiber_of[i], i] = 1
    for Hc, Hb in ((checks.H_Z, base.H_Z), (checks.H_X, base.H_X)):
        assert np.array_equal((R @ Hc) % 2, (Hb @ S) % 2), (
            "pushforward does not intertwine the check matrices; "
            f"code does not descend along σ = {sigma}"
        )

    # (2) Sector dichotomy: ker S = im P_lift (σ-invariant vectors).
    assert not ((S @ P_lift) % 2).any()          # p₊ p* = 2·id = 0
    assert _rank_f2(S) == nb and _rank_f2(P_lift) == nb, (
        "pushforward/pullback ranks wrong; fibers are not free"
    )
    # dim ker S = n − nb = nb = rank P_lift and im P_lift ⊆ ker S ⟹ equal.

    # (3) Class transport. Base pairing basis; Λ well-defined needs
    # p₊(rowspan H_X) ⊆ rowspan(H̄_X): rows of S·H_Xᵀ are p₊ of cover
    # X-stabilizer generators — pair them against L̄_Z.
    Lb_Z = find_logical_z(base)
    kb = Lb_Z.shape[0]
    push_stab = (S @ checks.H_X.T) % 2           # columns → pushed rows
    assert not ((Lb_Z @ push_stab) % 2).any(), (
        "p₊ of a cover X-stabilizer pairs nontrivially with a base "
        "Z-logical; Λ would be ill-defined"
    )
    k = action.k
    Lam = np.zeros((k, kb), dtype=np.uint8)
    for i in range(k):
        Lam[i] = (Lb_Z @ ((S @ action.V[i]) % 2)) % 2

    # Pμ (sector-(a) pins) and the (a)-reachable class subspace
    # W_a = Pμ(ker H̄_Z).
    Pmu = (action.L_Z @ P_lift) % 2              # (k, n̄)
    ker_b = nullspace_f2(base.H_Z)
    Wa_vecs = (ker_b @ Pmu.T) % 2                # rows: μ(p* v̄), v̄ ∈ ker basis
    # reduce to a row basis
    W_rows = []
    for row in Wa_vecs:
        if row.any() and not _in_span(
            np.array(W_rows, dtype=np.uint8).reshape(-1, k), row
        ):
            W_rows.append(row.copy())
    W_a = (
        np.array(W_rows, dtype=np.uint8).reshape(-1, k)
        if W_rows else np.zeros((0, k), dtype=np.uint8)
    )

    # Base translation action ρ̄ (for the anchor/pin compatibility check).
    from .shard_distance import compute_class_action as _cca  # local reuse
    try:
        base_action = _cca(base)
        rho_b = base_action.rho_rows
        # Transport commutation: Λ·ρ̄(t̄) = ρ(t)·Λ for every cover t.
        for t_idx, t in enumerate(elems):
            tb_idx = Gb.index(red(t))
            lhs = (action.rho_rows[t_idx] @ Lam) % 2
            rhs = (Lam @ rho_b[tb_idx]) % 2
            assert np.array_equal(lhs, rhs), (
                "Λ does not commute with transport; anchors would be "
                "unsound in sector (b)"
            )
    except ValueError:
        # base has k̄ = 0: Λ is the zero map, commutation is trivial.
        rho_b = tuple()

    # Spot-check Λ on random cycle elements (falsify-first).
    rng = np.random.default_rng(7)
    ker_c = nullspace_f2(checks.H_Z)
    for _ in range(20):
        coeff = rng.integers(0, 2, size=ker_c.shape[0]).astype(np.uint8)
        v = (coeff @ ker_c) % 2
        c = (action.L_Z @ v) % 2
        want = (c @ Lam) % 2
        got = (Lb_Z @ ((S @ v) % 2)) % 2
        assert np.array_equal(want, got), "Λ fails on a random cycle"

    return DescentData(
        sigma=sigma, base_checks=base, S=S, P_lift=P_lift,
        Lb_Z=Lb_Z, kb=kb, Pmu=Pmu, Lam=Lam, W_a=W_a,
        rho_b_rows=tuple(rho_b),
    )


# ---------------------------------------------------------------------------
# Sector SAT instances (CMS, native XOR)


def _solve_sector_a(
    dd: DescentData, reps: tuple[int, ...], k: int, weight_half: int,
) -> np.ndarray | None:
    """Invariant sector: base-sized instance. Variables v̄ over the n̄
    base qubits; constraints H̄_Z v̄ = 0 and (selector-conditional)
    class pins Pμ v̄ = rep; cardinality |v̄| ≤ ⌊w/2⌋. A model expands
    to the genuine cover member p* v̄ (re-verified by the driver)."""
    Hb_Z = dd.base_checks.H_Z
    nb = Hb_Z.shape[1]
    pool = IDPool()
    qv = [pool.id() for _ in range(nb)]
    solver = pycryptosat.Solver()

    for row in Hb_Z:
        idx = np.flatnonzero(row)
        if idx.size:
            solver.add_xor_clause([qv[i] for i in idx], False)

    if len(reps) == 1:
        cbits = _int_to_bits(reps[0], k)
        for j in range(k):
            idx = np.flatnonzero(dd.Pmu[j])
            if idx.size == 0:
                assert not cbits[j], "unreachable rep passed the W_a filter"
                continue
            solver.add_xor_clause([qv[i] for i in idx], bool(cbits[j]))
    else:
        outs = []
        for j in range(k):
            idx = np.flatnonzero(dd.Pmu[j])
            if idx.size == 0:
                outs.append(None)  # pin value forced to 0 for every rep
                continue
            oj = pool.id()
            solver.add_xor_clause([qv[i] for i in idx] + [oj], False)
            outs.append(oj)
        sel = [pool.id() for _ in reps]
        solver.add_clause(sel)
        for s, rep in zip(sel, reps):
            cbits = _int_to_bits(rep, k)
            for j in range(k):
                if outs[j] is None:
                    assert not cbits[j], "unreachable rep in sector (a)"
                    continue
                solver.add_clause([-s, outs[j] if cbits[j] else -outs[j]])

    if weight_half < nb:
        card = CardEnc.atmost(
            lits=qv, bound=weight_half, vpool=pool,
            encoding=EncType.seqcounter,
        )
        for cl in card.clauses:
            solver.add_clause(cl)

    sat, model = solver.solve()
    if not sat:
        return None
    return np.array([1 if model[v] else 0 for v in qv], dtype=np.uint8)


def _add_scaffold(
    solver, pool: IDPool, qv: list[int], dd: DescentData, weight: int,
) -> tuple[list[int], list[int | None]]:
    """Install one deck's pushforward scaffolding (no sector clause):
    aux bits w = p₊ v, implied base parity rows on them, μ̄ output
    literals for pins, and |w| ≤ weight. Returns (wv, μ̄-outs)."""
    nb = dd.S.shape[0]
    wv = [pool.id() for _ in range(nb)]
    for b in range(nb):
        idx = np.flatnonzero(dd.S[b])
        solver.add_xor_clause([qv[i] for i in idx] + [wv[b]], False)
    for row in dd.base_checks.H_Z:
        idx = np.flatnonzero(row)
        if idx.size:
            solver.add_xor_clause([wv[i] for i in idx], False)
    b_outs: list[int | None] = []
    for Lb in dd.Lb_Z:
        idx = np.flatnonzero(Lb)
        assert idx.size, "zero row in the base logical basis"
        ob = pool.id()
        solver.add_xor_clause([wv[i] for i in idx] + [ob], False)
        b_outs.append(ob)
    if weight < nb:
        card = CardEnc.atmost(
            lits=wv, bound=weight, vpool=pool, encoding=EncType.seqcounter
        )
        for cl in card.clauses:
            solver.add_clause(cl)
    return wv, b_outs


def _pin_lam(
    solver, dd: DescentData, b_outs: list[int | None],
    cbits: np.ndarray, guard: list[int],
) -> None:
    """Transport pins μ̄(w) = Λ·c (guarded by the rep selector)."""
    lam = (cbits @ dd.Lam) % 2
    for jb in range(dd.kb):
        solver.add_clause(
            guard + [b_outs[jb] if lam[jb] else -b_outs[jb]]
        )


def _solve_sector_b(
    H_Z: np.ndarray, L_Z: np.ndarray, dd: DescentData,
    reps: tuple[int, ...], weight: int, case: str,
    transversals: tuple[tuple[int, ...], ...],
    scaffolds: tuple[DescentData, ...] = (),
    group_floor: int | None = None,
) -> np.ndarray | None:
    """Moving sector: the cover shard instance of `_solve_shard_cms`
    plus the pushforward scaffolding of the split deck — aux bits
    w = p₊ v (XOR definitions), implied base parity rows H̄_Z w = 0,
    the sector clause ⋁ w, class-transport pins μ̄(w) = Λ·rep, and
    |w| ≤ weight (implied by |w| ≤ |v|). `group_floor` adds the
    certified cardinality LOWER bound |w| ≥ group_floor (sound when it
    is ≤ every selected rep's base-coset floor — the driver's
    obligation, computed by `base_coset_floors`). `scaffolds` installs
    the same implied structure (minus sector clause/floor) for
    additional verified decks."""
    n = H_Z.shape[1]
    N = n // 2
    k = L_Z.shape[0]
    pool = IDPool()
    qv = [pool.id() for _ in range(n)]
    solver = pycryptosat.Solver()

    for row in H_Z:
        idx = np.flatnonzero(row)
        if idx.size:
            solver.add_xor_clause([qv[i] for i in idx], False)

    wv, b_outs = _add_scaffold(solver, pool, qv, dd, weight)
    # Sector (b): the pushforward is nonzero.
    solver.add_clause(list(wv))
    if group_floor is not None and group_floor > 1:
        card = CardEnc.atleast(
            lits=wv, bound=group_floor, vpool=pool,
            encoding=EncType.seqcounter,
        )
        for cl in card.clauses:
            solver.add_clause(cl)
    extra = [
        _add_scaffold(solver, pool, qv, s_dd, weight)
        for s_dd in scaffolds
    ]

    if len(reps) == 1:
        cbits = _int_to_bits(reps[0], k)
        for j, L in enumerate(L_Z):
            idx = np.flatnonzero(L)
            solver.add_xor_clause([qv[i] for i in idx], bool(cbits[j]))
        _pin_lam(solver, dd, b_outs, cbits, [])
        for s_dd, (s_wv, s_outs) in zip(scaffolds, extra):
            _pin_lam(solver, s_dd, s_outs, cbits, [])
        T = transversals[0]
        if case == "L":
            solver.add_clause([qv[p] for p in T])
        else:
            for p in range(N):
                solver.add_clause([-qv[p]])
            solver.add_clause([qv[N + p] for p in T])
    else:
        a = []
        for L in L_Z:
            idx = np.flatnonzero(L)
            aj = pool.id()
            solver.add_xor_clause([qv[i] for i in idx] + [aj], False)
            a.append(aj)
        sel = [pool.id() for _ in reps]
        solver.add_clause(sel)
        if case == "R":
            for p in range(N):
                solver.add_clause([-qv[p]])
        for s, rep, T in zip(sel, reps, transversals):
            cbits = _int_to_bits(rep, k)
            for j in range(k):
                solver.add_clause([-s, a[j] if cbits[j] else -a[j]])
            _pin_lam(solver, dd, b_outs, cbits, [-s])
            for s_dd, (s_wv, s_outs) in zip(scaffolds, extra):
                _pin_lam(solver, s_dd, s_outs, cbits, [-s])
            if case == "L":
                solver.add_clause([-s] + [qv[p] for p in T])
            else:
                solver.add_clause([-s] + [qv[N + p] for p in T])

    if weight < n:
        card = CardEnc.atmost(
            lits=qv, bound=weight, vpool=pool, encoding=EncType.seqcounter
        )
        for cl in card.clauses:
            solver.add_clause(cl)

    sat, model = solver.solve()
    if not sat:
        return None
    return np.array([1 if model[v] else 0 for v in qv], dtype=np.uint8)


# ---------------------------------------------------------------------------
# Certified base-coset floors (the b3 increment)


def base_coset_floors(
    dd: DescentData, lam_classes: set[int], cap: int,
) -> dict[int, int]:
    """For each base class c̄ in `lam_classes`, the exact minimum
    weight of the *nonzero* members of its coset, capped: floors[c̄]=f
    with f ≤ cap means the min is exactly f; f = cap+1 means the coset
    has no nonzero member of weight ≤ cap. Certified by witness-jumping
    SAT calls on the *base* code (the same trusted encoding, base
    size): solve at ≤ b, jump b to witness-weight − 1, repeat until
    UNSAT. c̄ = 0 (the stabilizer coset) gets an explicit ⋁ w clause
    to exclude the zero vector."""
    Hb_Z = dd.base_checks.H_Z
    nb = Hb_Z.shape[1]

    def _try(cbits: np.ndarray, cb: int, w: int) -> np.ndarray | None:
        pool = IDPool()
        qv = [pool.id() for _ in range(nb)]
        solver = pycryptosat.Solver()
        for row in Hb_Z:
            idx = np.flatnonzero(row)
            if idx.size:
                solver.add_xor_clause([qv[i] for i in idx], False)
        for j, Lb in enumerate(dd.Lb_Z):
            idx = np.flatnonzero(Lb)
            solver.add_xor_clause([qv[i] for i in idx], bool(cbits[j]))
        if cb == 0:
            solver.add_clause(qv)
        card = CardEnc.atmost(
            lits=qv, bound=w, vpool=pool, encoding=EncType.seqcounter
        )
        for cl in card.clauses:
            solver.add_clause(cl)
        sat, model = solver.solve()
        if not sat:
            return None
        return np.array([1 if model[v] else 0 for v in qv], dtype=np.uint8)

    out: dict[int, int] = {}
    for cb in lam_classes:
        cbits = _int_to_bits(cb, dd.kb) if dd.kb else np.zeros(0, np.uint8)
        bound, best = cap, cap + 1
        while bound >= 1:
            vb = _try(cbits, cb, bound)
            if vb is None:
                break
            best = int(vb.sum())
            assert 0 < best <= bound
            bound = best - 1
        out[cb] = best
    return out


# ---------------------------------------------------------------------------
# Driver


@dataclass(frozen=True, slots=True)
class SectorStat:
    sector: str          # 'A' | 'B-L' | 'B-R' | 'A-analytic'
    group: int
    n_reps: int
    weight: int
    status: str          # 'SAT' | 'UNSAT' | 'SKIP'
    seconds: float


@dataclass(frozen=True, slots=True)
class DescentDistanceResult:
    distance: int
    witness: np.ndarray
    sigma: tuple[int, ...]
    k: int
    kb: int
    num_orbits: int
    rounds: int
    initial_upper_bound: int
    a_reachable_reps: int
    a_killed_reps: int
    sector_stats: tuple[SectorStat, ...] = field(repr=False)
    wall_seconds: float = 0.0
    floor_seconds: float = 0.0     # certified base-floor precompute

    @property
    def solver_cpu(self) -> float:
        return sum(s.seconds for s in self.sector_stats)


def descent_shard_distance(
    checks: CheckMatrices,
    A: Poly,
    B: Poly,
    *,
    sigma: tuple[int, ...] | None = None,
    scaffold_rest: bool = False,
    use_floors: bool = False,
    verbose: bool = False,
    probe_confl_limit: int = 100_000,
) -> DescentDistanceResult:
    """Exact d_X via the descent-strengthened shard decomposition.

    `sigma=None` picks the first axis deck (prefer the ℓ-axis). The
    layer is sound for any verified axis deck; which deck is *fastest*
    is an empirical question — benchmark both when both exist.
    `scaffold_rest=True` additionally installs every *other* axis
    deck's implied pushforward structure (aux bits, base rows, Λ pins)
    inside sector (b).
    `use_floors=True` certifies per-base-class coset floors by tiny
    base-code SAT (`base_coset_floors`) and uses them to (i) drop reps
    whose sector (b) is provably empty at the round weight and (ii)
    add the group-min |p₊ v| ≥ floor cardinality.
    """
    if not _HAVE_CMS:
        raise RuntimeError("descent_sat requires pycryptosat")
    t_start = time.perf_counter()
    action = compute_class_action(checks)
    L_Z, V, k = action.L_Z, action.V, action.k

    if sigma is None:
        cands = axis_decks(checks)
        if not cands:
            raise ValueError("no axis deck: |G| has no even axis")
        sigma = cands[0]
    dd = compute_descent(checks, action, sigma, A=A, B=B)
    scaffolds: tuple[DescentData, ...] = ()
    if scaffold_rest:
        scaffolds = tuple(
            compute_descent(checks, action, s2, A=A, B=B)
            for s2 in axis_decks(checks) if s2 != sigma
        )

    witness = _greedy_upper_bound(V, checks.H_X)
    d_ub = int(witness.sum())
    initial_ub = d_ub
    if verbose:
        print(
            f"[descent] σ={sigma}  base n̄={dd.S.shape[0]} k̄={dd.kb}  "
            f"k={k} orbits={len(action.orbit_reps)}  initial d_ub={d_ub}",
            flush=True,
        )

    while probe_confl_limit > 0 and d_ub > 1:
        v = _monolith_probe(checks.H_Z, L_Z, d_ub - 1, probe_confl_limit)
        if v is None:
            break
        _verify_witness(checks, L_Z, v, int(v.sum()))
        witness, d_ub = v, int(v.sum())
        if verbose:
            print(f"[descent] probe improved d_ub → {d_ub}", flush=True)

    # Parity machinery (identical to shard_distance).
    hx_even = not any(int(row.sum()) % 2 for row in checks.H_X)
    v_parity = V.sum(axis=1) % 2
    rep_parity = {
        r: int(_int_to_bits(r, k) @ v_parity) % 2
        for r in action.orbit_reps
    }

    def _tight(w: int, parity: int) -> int:
        return w - ((w - parity) % 2) if hx_even else w

    groups: list[tuple[int, ...]] = []
    for p in (0, 1):
        block = tuple(r for r in action.orbit_reps if rep_parity[r] == p)
        if block:
            groups.append(block)
    group_transversals = [
        tuple(action.transversal[r] for r in grp) for grp in groups
    ]

    # Sector-(a) analytic filter: reachable ⟺ even parity ∧ c ∈ W_a.
    a_reps: list[tuple[int, ...]] = []
    killed = 0
    for grp in groups:
        keep = []
        for r in grp:
            cbits = _int_to_bits(r, k)
            if rep_parity[r] == 0 and _in_span(dd.W_a, cbits):
                keep.append(r)
            else:
                killed += 1
        a_reps.append(tuple(keep))
    if verbose:
        tot = sum(len(g) for g in groups)
        reach = sum(len(g) for g in a_reps)
        print(
            f"[descent] sector (a): {reach}/{tot} orbit reps reachable "
            f"({killed} killed analytically: parity + W_a rank)",
            flush=True,
        )

    # Certified base-coset floors for the Λ-image classes (b3).
    floors: dict[int, int] = {}
    floor_secs = 0.0
    if use_floors:
        t0 = time.perf_counter()
        lam_classes = {
            _bits_to_int((_int_to_bits(r, k) @ dd.Lam) % 2)
            for r in action.orbit_reps
        }
        floors = base_coset_floors(dd, lam_classes, cap=d_ub - 1)
        floor_secs = time.perf_counter() - t0
        if verbose:
            vals = sorted(set(floors.values()))
            print(
                f"[descent] base-coset floors: {len(floors)} classes, "
                f"values {vals} ({floor_secs:.2f}s)",
                flush=True,
            )

    def _lam_floor(r: int) -> int:
        if not floors:
            return 0
        return floors[_bits_to_int((_int_to_bits(r, k) @ dd.Lam) % 2)]

    refuted: dict[tuple[str, int], int] = {}
    stats: list[SectorStat] = []
    rounds = 0

    while True:
        w = d_ub - 1
        rounds += 1
        improved = False

        tasks: list[tuple] = []   # (kind, gi, payload)
        for gi, grp in enumerate(groups):
            wt = _tight(w, rep_parity[grp[0]])
            if wt < 1:
                continue
            if a_reps[gi] and refuted.get(("A", gi), -1) < w:
                tasks.append(("A", gi, a_reps[gi], wt // 2))
            # Floor filter: reps whose base coset has no nonzero member
            # of weight ≤ wt have a certifiably empty sector (b).
            b_reps = tuple(r for r in grp if _lam_floor(r) <= wt)
            b_trans = tuple(
                action.transversal[r] for r in b_reps
            )
            if len(b_reps) < len(grp) and verbose and rounds == 1:
                print(
                    f"[descent] group {gi}: {len(grp) - len(b_reps)} reps "
                    f"floor-killed in sector (b) at w ≤ {wt}",
                    flush=True,
                )
            gf = (
                min((_lam_floor(r) for r in b_reps), default=0)
                if floors else None
            )
            for case in ("L", "R"):
                if b_reps and refuted.get((f"B-{case}", gi), -1) < w:
                    tasks.append(
                        ("B" + case, gi, b_reps, wt, b_trans, gf)
                    )
        if not tasks:
            break

        for task in tasks:
            kind, gi = task[0], task[1]
            t0 = time.perf_counter()
            if kind == "A":
                _, _, reps, wh = task
                vb = _solve_sector_a(dd, reps, k, wh)
                dt = time.perf_counter() - t0
                if vb is None:
                    refuted[("A", gi)] = w
                    stats.append(
                        SectorStat("A", gi, len(reps), wh, "UNSAT", dt)
                    )
                    continue
                v = (dd.P_lift @ vb) % 2
            else:
                case = kind[1]
                _, _, reps, wt, trans, gf = task
                v = _solve_sector_b(
                    checks.H_Z, L_Z, dd, reps, wt, case, trans,
                    scaffolds=scaffolds, group_floor=gf,
                )
                dt = time.perf_counter() - t0
                if v is None:
                    refuted[(f"B-{case}", gi)] = w
                    stats.append(
                        SectorStat(
                            f"B-{case}", gi, len(reps), wt, "UNSAT", dt
                        )
                    )
                    continue
            stats.append(
                SectorStat(
                    kind if kind == "A" else f"B-{kind[1]}",
                    gi, len(task[2]), task[3], "SAT", dt,
                )
            )
            v = v.astype(np.uint8)
            _verify_witness(checks, L_Z, v, int(v.sum()))
            witness, d_ub = v, int(v.sum())
            improved = True
            break

        if verbose:
            done = [s for s in stats if s.status == "UNSAT"]
            print(
                f"[descent] round {rounds}: w ≤ {w}, "
                f"{len(done)} sector tasks refuted so far"
                + (f", improved d_ub → {d_ub}" if improved else ""),
                flush=True,
            )
        if not improved:
            break

    _verify_witness(checks, L_Z, witness, d_ub)
    return DescentDistanceResult(
        distance=d_ub,
        witness=witness,
        sigma=sigma,
        k=k,
        kb=dd.kb,
        num_orbits=len(action.orbit_reps),
        rounds=rounds,
        initial_upper_bound=initial_ub,
        a_reachable_reps=sum(len(g) for g in a_reps),
        a_killed_reps=killed,
        sector_stats=tuple(stats),
        wall_seconds=time.perf_counter() - t_start,
        floor_seconds=floor_secs,
    )


def invariant_floor(dd: DescentData, cap: int | None = None) -> int:
    """Exact minimum COVER weight of a σ-invariant nontrivial logical
    (= 2·min |v̄| over base v̄ with H̄_Z v̄ = 0, μ(p* v̄) ≠ 0), by
    witness-jumping CMS on the sector-(a) instance. Returns
    2·(cap+1) if the sector is empty at |v̄| ≤ cap (cap defaults to
    n̄: exhaustive, still tiny)."""
    Hb_Z = dd.base_checks.H_Z
    nb = Hb_Z.shape[1]
    if cap is None:
        cap = nb

    def _try(w: int) -> np.ndarray | None:
        pool = IDPool()
        qv = [pool.id() for _ in range(nb)]
        solver = pycryptosat.Solver()
        for row in Hb_Z:
            idx = np.flatnonzero(row)
            if idx.size:
                solver.add_xor_clause([qv[i] for i in idx], False)
        outs = []
        for j in range(dd.Pmu.shape[0]):
            idx = np.flatnonzero(dd.Pmu[j])
            if idx.size == 0:
                continue
            oj = pool.id()
            solver.add_xor_clause([qv[i] for i in idx] + [oj], False)
            outs.append(oj)
        if not outs:
            return None          # no invariant logicals at all
        solver.add_clause(outs)
        if w < nb:
            card = CardEnc.atmost(
                lits=qv, bound=w, vpool=pool, encoding=EncType.seqcounter
            )
            for cl in card.clauses:
                solver.add_clause(cl)
        sat, model = solver.solve()
        if not sat:
            return None
        return np.array([1 if model[v] else 0 for v in qv], dtype=np.uint8)

    bound, best = cap, cap + 1
    while bound >= 1:
        vb = _try(bound)
        if vb is None:
            break
        best = int(vb.sum())
        assert 0 < best <= bound
        bound = best - 1
    return 2 * best


def base_coset_floors_budgeted(
    dd: DescentData,
    lam_classes: set[int],
    cap: int,
    confl_budget: int = 100_000,
) -> dict[int, int]:
    """Sound floors under a per-call conflict budget: ascend w = 1..cap;
    UNSAT at w extends the floor to w+1, the first SAT stops with the
    exact floor, and a budget-exhausted call stops with the floor
    proven so far (= last UNSAT bound + 1). Never overclaims: every
    returned value f asserts only 'no nonzero coset member of weight
    < f', which was explicitly refuted. The exact engine
    (`base_coset_floors`) remains for small bases."""
    Hb_Z = dd.base_checks.H_Z
    nb = Hb_Z.shape[1]
    out: dict[int, int] = {}
    for cb in lam_classes:
        cbits = _int_to_bits(cb, dd.kb) if dd.kb else np.zeros(0, np.uint8)
        floor = 1
        for w in range(1, cap + 1):
            pool = IDPool()
            qv = [pool.id() for _ in range(nb)]
            solver = pycryptosat.Solver(confl_limit=confl_budget)
            for row in Hb_Z:
                idx = np.flatnonzero(row)
                if idx.size:
                    solver.add_xor_clause([qv[i] for i in idx], False)
            for j, Lb in enumerate(dd.Lb_Z):
                idx = np.flatnonzero(Lb)
                solver.add_xor_clause(
                    [qv[i] for i in idx], bool(cbits[j])
                )
            if cb == 0:
                solver.add_clause(qv)
            card = CardEnc.atmost(
                lits=qv, bound=w, vpool=pool, encoding=EncType.seqcounter
            )
            for cl in card.clauses:
                solver.add_clause(cl)
            sat, _ = solver.solve()
            if sat is True:
                break               # exact: min weight is exactly w
            if sat is None:
                break               # budget: keep the floor proven so far
            floor = w + 1           # UNSAT at ≤ w
        out[cb] = floor
    return out
