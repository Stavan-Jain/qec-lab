"""Structure-exploiting exact distance for BB codes.

The monolithic SAT instance in `sat_distance.py` asks for a nontrivial
logical of weight ≤ w in one query whose logical-class freedom (the
`⋁ a_j` clause ranging over all 2^k − 1 nonzero classes) gives the
solver almost no propagation; the UNSAT round at w = d − 1 is where it
dies. This module strengthens/decomposes that instance using the
translation symmetry, via four sound ingredients:

1.  **Classes.** For v ∈ ker(H_Z) the pairing μ(v) = L_Z · v ∈ F₂^k is
    constant on cosets of rowspan(H_X) and classifies them; v is a
    nontrivial logical iff μ(v) ≠ 0. Hence

        d_X = min over c ≠ 0 of  min{ |v| : H_Z v = 0, μ(v) = c }.

    A shard is one such c, encoded by pinning each ⟨L_Z[j], v⟩ to c_j —
    k *unit* XOR constraints instead of one wide disjunction.

2.  **Transport.** The translation group G acts on qubits by the
    block-diagonal permutation π_t (both blocks translated by t). π_t
    preserves ker(H_Z), rowspan(H_X), and weight, so it induces
    ρ(t) ∈ GL_k(F₂) on classes with μ(π_t v) = ρ(t) μ(v). Shard minima
    are constant on ρ(G)-orbits: only one shard per orbit is solved.

3.  **Anchors.** For a rep c, translations by Stab(c) = {t : ρ(t)c = c}
    map the shard to itself, so WLOG the support meets a transversal T
    of the Stab(c)-translation action on positions:
      case 'L' — some left-block position in T is set;
      case 'R' — left block ≡ 0 and some right-block position in T set.
    Both cases only *restrict* the shard, and every shard member has a
    Stab(c)-translate in one of them, so refuting both refutes the
    shard. When Stab(c) = G the case-L anchor is a unit literal.

4.  **Descent loop.** Start from a cheap upper bound (stabilizer-
    reduced logical basis rows, then conflict-budgeted monolithic SAT
    probes). Each round targets w = d_ub − 1: a SAT shard yields a
    strictly better witness (new d_ub, next round); a fully-UNSAT
    round proves d = d_ub. "UNSAT at weight ≤ w" is monotone downward
    in w, so refuted shards are never re-solved. Coset weight parity
    (see `_tight` in the driver) additionally lowers each group's
    refutation bound to the largest weight of the group's parity.

**Granularity (measured, 2026-07-29).** The original fine-grained
design — one task per orbit rep — LOSES to the monolith by 5–9× at
n = 72–108: per-task cost is dominated by shared structure, so task
count is the enemy. The winner at every size tested is the *coarse*
end of the dial (`orbits_per_shard=None`, the default): one selector-
encoded task per parity class covering all orbit reps. That instance
is the monolith *strengthened* by three sound extra constraint
families (classes restricted to orbit reps; per-rep support anchors;
parity-tightened cardinality), and it beats the monolithic baseline
by ~4–6× solver-CPU on bb_90/bb_108. Fine granularities remain
available for parallel refutation and for per-shard proof emission
(the future certificate mode), where small independent instances are
the point.

Soundness does not depend on the orbit/anchor computation being
*optimal*, only on facts asserted numerically at construction: μ
classifies; ρ comes from weight-preserving code automorphisms
(translation invariance of H_Z/H_X is checked entrywise); the orbits
partition the 2^k − 1 classes; the anchor cases cover each shard up
to Stab(c)-translation; parity tightening is applied only when every
H_X row has even weight (which makes weight parity a coset
invariant). Every SAT model is independently re-verified as a
nontrivial logical before it is accepted as a witness.

Scope (v1): BB inputs (`CheckMatrices` from `bb_check_matrices`),
translation symmetry only (Aut extensions beyond G, per the A9
machinery, are a later upgrade), X direction only (BB transpose
symmetry gives d_Z = d_X, as in `sat_distance.py`).
"""

from __future__ import annotations

import multiprocessing as mp
import time
from dataclasses import dataclass, field

import numpy as np

from pysat.card import CardEnc, EncType
from pysat.formula import CNF, IDPool
from pysat.solvers import Cadical195

from .checks import CheckMatrices
from .linalg import nullspace_f2, quotient_complement_basis
from .sat_distance import _xor_chain, find_logical_z

try:
    import pycryptosat
    _HAVE_CMS = True
except ImportError:
    _HAVE_CMS = False


MAX_K_FOR_CLASS_ENUM = 20  # 2^k class enumeration guard


# ---------------------------------------------------------------------------
# F₂ helpers local to the class-action computation


def _inv_f2(M: np.ndarray) -> np.ndarray:
    """Invert a k×k matrix over F₂ (raises if singular)."""
    k = M.shape[0]
    A = np.concatenate(
        [(M & 1).astype(np.uint8), np.eye(k, dtype=np.uint8)], axis=1
    )
    pivot_row = 0
    for col in range(k):
        nz = np.flatnonzero(A[pivot_row:, col])
        if nz.size == 0:
            raise ValueError("matrix is singular over F₂")
        r = pivot_row + int(nz[0])
        if r != pivot_row:
            A[[pivot_row, r]] = A[[r, pivot_row]]
        mask = A[:, col] == 1
        mask[pivot_row] = False
        if mask.any():
            A[mask] ^= A[pivot_row]
        pivot_row += 1
    return A[:, k:].copy()


def _int_to_bits(c: int, k: int) -> np.ndarray:
    return np.array([(c >> j) & 1 for j in range(k)], dtype=np.uint8)


def _bits_to_int(bits: np.ndarray) -> int:
    return int(sum(1 << j for j in range(bits.shape[0]) if bits[j]))


# ---------------------------------------------------------------------------
# Class action: ρ : G → GL_k(F₂), orbits, stabilizer transversals


@dataclass(frozen=True, slots=True)
class ClassAction:
    """The induced action of the translation group on logical classes.

    `rho_rows[t_idx]` is the k×k matrix R_t with *rows* μ(π_t v_i), so
    the action on a class bit-vector is `c_bits @ R_t (mod 2)`.
    `transversal[rep]` is the anchor transversal of Stab(rep) acting on
    the |G| positions (shared by both blocks).
    """

    k: int
    n: int                                  # 2|G| qubits
    L_Z: np.ndarray                         # (k, n) pairing basis
    V: np.ndarray                           # (k, n) X-logicals, μ(V[i]) = e_i
    orbit_reps: tuple[int, ...]             # one class int per orbit
    orbit_sizes: tuple[int, ...]
    orbit_rep_of: np.ndarray                # (2^k,) class -> its orbit rep
    transversal: dict[int, tuple[int, ...]]  # rep -> anchor positions
    stab_sizes: dict[int, int]              # rep -> |Stab(rep)|
    rho_rows: tuple[np.ndarray, ...]        # per t: R_t with c ↦ c·R_t


def _translation_perm(group, t: tuple[int, ...]) -> np.ndarray:
    """Qubit permutation of π_t: qubit (b, h) ↦ (b, h + t).

    Returned as `perm` with `perm[i] = image of qubit i`; a translated
    vector is `w` with `w[perm] = v`.
    """
    N = group.cardinality
    pos = np.empty(N, dtype=np.int64)
    for i in range(N):
        pos[i] = group.index(group.add(group.from_index(i), t))
    return np.concatenate([pos, pos + N])


def compute_class_action(checks: CheckMatrices) -> ClassAction:
    """Compute μ-normalized logical bases, ρ(G), orbits, and anchors.

    Asserts (numerically) every soundness prerequisite listed in the
    module docstring; raises AssertionError on any failure rather than
    returning a wrong decomposition.
    """
    G = checks.group
    N = G.cardinality
    n = checks.num_qubits
    elems = list(G)

    L_Z = find_logical_z(checks)
    k = L_Z.shape[0]
    if k == 0:
        raise ValueError("code has k = 0; distance is undefined")
    if k > MAX_K_FOR_CLASS_ENUM:
        raise ValueError(
            f"k = {k} exceeds the 2^k class-enumeration guard "
            f"({MAX_K_FOR_CLASS_ENUM}); shard decomposition needs a "
            "sparser class representation first"
        )

    # X-logical basis, then normalize so that μ(V[i]) = e_i.
    ker_Z = nullspace_f2(checks.H_Z)
    V0 = quotient_complement_basis(checks.H_X, ker_Z)
    assert V0.shape[0] == k, (
        f"X-logical count {V0.shape[0]} != Z-logical count {k}"
    )
    M = (V0 @ L_Z.T) % 2
    V = (_inv_f2(M) @ V0) % 2
    V = V.astype(np.uint8)
    assert np.array_equal((V @ L_Z.T) % 2, np.eye(k, dtype=np.uint8)), (
        "pairing normalization failed"
    )

    # Translation invariance of the check matrices (soundness of ρ):
    # H[t + r, π_t(c)] == H[r, c] for every translation t.
    rho_rows: list[np.ndarray] = []
    pos_perms: list[np.ndarray] = []
    for t in elems:
        perm = _translation_perm(G, t)
        pos_perms.append(perm[:N].copy())
        row_perm = perm[:N]  # checks are indexed by G in the same order
        for H in (checks.H_Z, checks.H_X):
            assert np.array_equal(H[np.ix_(row_perm, perm)], H), (
                f"H not translation-invariant at t = {t}; "
                "shard transport would be unsound for this input"
            )
        Vt = np.zeros_like(V)
        Vt[:, perm] = V  # rows are π_t v_i
        R = (Vt @ L_Z.T) % 2  # row i = μ(π_t v_i)
        rho_rows.append(R.astype(np.uint8))

    # ρ(0) = identity; spot-check the homomorphism property.
    id_idx = G.index(tuple(0 for _ in G.orders))
    assert np.array_equal(rho_rows[id_idx], np.eye(k, dtype=np.uint8))
    rng = np.random.default_rng(0)
    for _ in range(min(50, N * N)):
        i, j = int(rng.integers(N)), int(rng.integers(N))
        gij = G.index(G.add(elems[i], elems[j]))
        assert np.array_equal(
            (rho_rows[i] @ rho_rows[j]) % 2, rho_rows[gij]
        ), "ρ is not a homomorphism (bug in the perm/index bookkeeping)"

    # Orbits of the nonzero classes under ρ(G). Ascending scan makes
    # each orbit's first-seen member (the minimum) its rep.
    size = 1 << k
    orbit_of = np.full(size, -1, dtype=np.int64)
    orbit_reps: list[int] = []
    orbit_sizes: list[int] = []
    transversal: dict[int, tuple[int, ...]] = {}
    stab_sizes: dict[int, int] = {}
    for c in range(1, size):
        if orbit_of[c] != -1:
            continue
        cbits = _int_to_bits(c, k)
        stab: list[int] = []
        members: set[int] = set()
        for t_idx in range(N):
            img = _bits_to_int((cbits @ rho_rows[t_idx]) % 2)
            members.add(img)
            if img == c:
                stab.append(t_idx)
        for mmb in members:
            assert orbit_of[mmb] in (-1, c), "orbits are not disjoint"
            orbit_of[mmb] = c
        orbit_reps.append(c)
        orbit_sizes.append(len(members))
        assert len(members) * len(stab) == N, (
            "orbit-stabilizer mismatch on classes"
        )
        # Anchor transversal: orbit reps of Stab(c) translating the N
        # positions (the action is the same on either block).
        seen = np.zeros(N, dtype=bool)
        T: list[int] = []
        for p in range(N):
            if seen[p]:
                continue
            T.append(p)
            for s_idx in stab:
                seen[pos_perms[s_idx][p]] = True
        transversal[c] = tuple(T)
        stab_sizes[c] = len(stab)

    assert int(sum(orbit_sizes)) == size - 1, (
        "class orbits do not partition the nonzero classes"
    )

    return ClassAction(
        k=k, n=n, L_Z=L_Z, V=V,
        orbit_reps=tuple(orbit_reps),
        orbit_sizes=tuple(orbit_sizes),
        orbit_rep_of=orbit_of,
        transversal=transversal,
        stab_sizes=stab_sizes,
        rho_rows=tuple(rho_rows),
    )


# ---------------------------------------------------------------------------
# Single-shard SAT instances


def _solve_shard_cms(
    H_Z: np.ndarray, L_Z: np.ndarray, reps: tuple[int, ...],
    weight: int, case: str, transversals: tuple[tuple[int, ...], ...],
) -> np.ndarray | None:
    """CryptoMiniSat backend, parity rows as native XOR.

    Single-rep shard: class pins are unit XOR constraints (strongest).
    Multi-rep "super-shard": one selector s_r per orbit rep; s_r forces
    the pairing outputs to r and the anchor to r's transversal, and
    `⋁ s_r` requires some rep to be selected. Sound because any
    weight-≤w logical whose class-orbit rep lies in `reps` has a
    transport into one of the anchored cases (module docstring, items
    2–3); pins keep every model a genuine nontrivial logical.
    """
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

    if len(reps) == 1:
        cbits = _int_to_bits(reps[0], k)
        for j, L in enumerate(L_Z):
            idx = np.flatnonzero(L)
            assert idx.size, "zero row in the logical-Z basis"
            solver.add_xor_clause([qv[i] for i in idx], bool(cbits[j]))
        T = transversals[0]
        if case == "L":
            solver.add_clause([qv[p] for p in T])
        else:  # case R: left block ≡ 0, right-block support meets T
            for p in range(N):
                solver.add_clause([-qv[p]])
            solver.add_clause([qv[N + p] for p in T])
    else:
        # pairing indicators a_j ≡ ⟨L_Z[j], v⟩
        a = []
        for L in L_Z:
            idx = np.flatnonzero(L)
            assert idx.size, "zero row in the logical-Z basis"
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


def _solve_shard_cadical(
    H_Z: np.ndarray, L_Z: np.ndarray, reps: tuple[int, ...],
    weight: int, case: str, transversals: tuple[tuple[int, ...], ...],
) -> np.ndarray | None:
    """pysat CaDiCaL backend: Tseitin XOR chains; same shard semantics
    as `_solve_shard_cms` (unit pins when single-rep, selectors when
    multi-rep)."""
    n = H_Z.shape[1]
    N = n // 2
    k = L_Z.shape[0]
    pool = IDPool()
    qv = [pool.id() for _ in range(n)]
    cnf = CNF()

    for row in H_Z:
        idx = np.flatnonzero(row)
        if idx.size == 0:
            continue
        out = _xor_chain((qv[i] for i in idx), pool, cnf)
        if out is not None:
            cnf.append([-out])

    a = []
    for L in L_Z:
        idx = np.flatnonzero(L)
        assert idx.size, "zero row in the logical-Z basis"
        a.append(_xor_chain((qv[i] for i in idx), pool, cnf))

    if len(reps) == 1:
        cbits = _int_to_bits(reps[0], k)
        for j in range(k):
            cnf.append([a[j]] if cbits[j] else [-a[j]])
        T = transversals[0]
        if case == "L":
            cnf.append([qv[p] for p in T])
        else:
            for p in range(N):
                cnf.append([-qv[p]])
            cnf.append([qv[N + p] for p in T])
    else:
        sel = [pool.id() for _ in reps]
        cnf.append(sel)
        if case == "R":
            for p in range(N):
                cnf.append([-qv[p]])
        for s, rep, T in zip(sel, reps, transversals):
            cbits = _int_to_bits(rep, k)
            for j in range(k):
                cnf.append([-s, a[j] if cbits[j] else -a[j]])
            if case == "L":
                cnf.append([-s] + [qv[p] for p in T])
            else:
                cnf.append([-s] + [qv[N + p] for p in T])

    if weight < n:
        card = CardEnc.atmost(
            lits=qv, bound=weight, vpool=pool, encoding=EncType.seqcounter
        )
        cnf.extend(card.clauses)

    solver = Cadical195(bootstrap_with=cnf.clauses)
    try:
        if not solver.solve():
            return None
        truth = {abs(lit): lit > 0 for lit in solver.get_model()}
        return np.array(
            [1 if truth.get(v, False) else 0 for v in qv], dtype=np.uint8
        )
    finally:
        solver.delete()


def _row_basis(M: np.ndarray) -> np.ndarray:
    """Rows of `M` forming a basis of rowspan(M) (greedy rank filter)."""
    basis: list[np.ndarray] = []
    acc: list[np.ndarray] = []  # rref of the kept rows, maintained inline
    for row in M:
        r = row.copy()
        for b in acc:
            pivot = int(np.flatnonzero(b)[0])
            if r[pivot]:
                r ^= b
        if r.any():
            basis.append(row.copy())
            # insert reduced row, keep acc reduced
            for b in acc:
                pivot = int(np.flatnonzero(r)[0])
                if b[pivot]:
                    b ^= r
            acc.append(r)
    return np.stack(basis) if basis else np.zeros((0, M.shape[1]), np.uint8)


def _solve_shard_coset_cms(
    H_X_basis: np.ndarray, offsets: tuple[np.ndarray, ...],
    weight: int, case: str, transversals: tuple[tuple[int, ...], ...],
) -> np.ndarray | None:
    """Coset-leader encoding: shard(c) = v_c + rowspan(H_X), so the
    instance has NO feasibility constraints at all — free coefficients
    u over the H_X row basis, qubit bits *defined* by short XORs
    (column weight of H_X ≈ |A| + |B| = 6 for BB codes), plus anchors
    and the cardinality bound:

        w_i  ≡  v_c[i] ⊕ ⟨(H_X_basis)ᵀ[i], u⟩.

    Every assignment of u is a genuine member of the coset (CSS gives
    rowspan(H_X) ⊆ ker H_Z and μ is coset-invariant), so the solver's
    whole fight is weight-vs-cardinality — the coupling QDistSAT
    identified as the actual bottleneck — through 6-literal XORs
    instead of |G|/2-wide parity rows.

    Multi-rep groups: selector s_r mixes rep r's offset into each
    defining XOR (valid because selectors are exactly-one: at-most-one
    is enforced pairwise, at-least-one by the ⋁ s_r clause; without
    AMO two active offsets could silently combine into a *trivial*
    coset).
    """
    n = H_X_basis.shape[1]
    N = n // 2
    pool = IDPool()
    uv = [pool.id() for _ in range(H_X_basis.shape[0])]
    wv = [pool.id() for _ in range(n)]
    solver = pycryptosat.Solver()

    single = len(offsets) == 1
    sel: list[int] = []
    if not single:
        sel = [pool.id() for _ in offsets]
        solver.add_clause(sel)
        for i in range(len(sel)):          # pairwise at-most-one
            for j in range(i + 1, len(sel)):
                solver.add_clause([-sel[i], -sel[j]])

    col_support = [np.flatnonzero(H_X_basis[:, i]) for i in range(n)]
    for i in range(n):
        lits = [wv[i]] + [uv[r] for r in col_support[i]]
        if single:
            solver.add_xor_clause(lits, bool(offsets[0][i]))
        else:
            # w_i ⊕ ⟨col_i, u⟩ ⊕ Σ_r s_r·v_r[i] = 0
            lits += [s for s, off in zip(sel, offsets) if off[i]]
            solver.add_xor_clause(lits, False)

    if case == "R":
        for p in range(N):
            solver.add_clause([-wv[p]])
    for gi, T in enumerate(transversals):
        guard = [] if single else [-sel[gi]]
        if case == "L":
            solver.add_clause(guard + [wv[p] for p in T])
        else:
            solver.add_clause(guard + [wv[N + p] for p in T])

    if weight < n:
        card = CardEnc.atmost(
            lits=wv, bound=weight, vpool=pool, encoding=EncType.seqcounter
        )
        for cl in card.clauses:
            solver.add_clause(cl)

    sat, model = solver.solve()
    if not sat:
        return None
    return np.array([1 if model[v] else 0 for v in wv], dtype=np.uint8)


# Worker-process state (set once per worker by the pool initializer so
# the check matrices are not re-pickled for every task).
_WORKER: dict = {}


def _init_worker(
    H_Z: np.ndarray, L_Z: np.ndarray, backend: str,
    encoding: str = "pin",
    H_X_basis: np.ndarray | None = None,
    V: np.ndarray | None = None,
) -> None:
    _WORKER["H_Z"] = H_Z
    _WORKER["L_Z"] = L_Z
    _WORKER["backend"] = backend
    _WORKER["encoding"] = encoding
    _WORKER["H_X_basis"] = H_X_basis
    _WORKER["V"] = V


def _run_task(task: tuple) -> tuple:
    """One (super-)shard solve.
    `task = (group_idx, reps, case, weight, transversals)`; returns
    `(group_idx, reps, case, weight, witness_support | None, seconds)`.
    """
    gidx, reps, case, weight, transversals = task
    H_Z, L_Z = _WORKER["H_Z"], _WORKER["L_Z"]
    backend = _WORKER["backend"]
    t0 = time.perf_counter()
    if _WORKER["encoding"] == "coset":
        V = _WORKER["V"]
        k = L_Z.shape[0]
        offsets = tuple(
            (_int_to_bits(rep, k) @ V) % 2 for rep in reps
        )
        v = _solve_shard_coset_cms(
            _WORKER["H_X_basis"], offsets, weight, case, transversals
        )
    elif backend == "cms":
        v = _solve_shard_cms(H_Z, L_Z, reps, weight, case, transversals)
    else:
        v = _solve_shard_cadical(H_Z, L_Z, reps, weight, case, transversals)
    dt = time.perf_counter() - t0
    support = None if v is None else tuple(int(i) for i in np.flatnonzero(v))
    return gidx, reps, case, weight, support, dt


# ---------------------------------------------------------------------------
# Budgeted monolithic SAT probes (descent accelerator)
#
# During descent the expensive failure mode is refuting shards at a
# weight ABOVE the true distance (bench: bb_108 burned ~140 s at
# w ≤ 11 when d = 10). A conflict-budgeted monolithic CMS probe finds
# improving witnesses cheaply while they exist; once a probe fails to
# improve within budget we hand over to the shard rounds. Probes carry
# *improvement* authority only (every witness is re-verified);
# refutation authority stays exclusively with the shard rounds, so
# soundness is unaffected by the budget heuristic.


def _monolith_probe(
    H_Z: np.ndarray, L_Z: np.ndarray, weight: int, confl_limit: int
) -> np.ndarray | None:
    """One conflict-budgeted CMS call on the *monolithic* instance
    (any nontrivial class) at ≤ `weight`. Returns a witness or None
    (None covers both UNSAT and budget exhaustion)."""
    if not _HAVE_CMS:
        return None
    n = H_Z.shape[1]
    pool = IDPool()
    qv = [pool.id() for _ in range(n)]
    solver = pycryptosat.Solver(confl_limit=confl_limit)
    for row in H_Z:
        idx = np.flatnonzero(row)
        if idx.size:
            solver.add_xor_clause([qv[i] for i in idx], False)
    a_outs = []
    for L in L_Z:
        idx = np.flatnonzero(L)
        a = pool.id()
        solver.add_xor_clause([qv[i] for i in idx] + [a], False)
        a_outs.append(a)
    solver.add_clause(a_outs)
    if weight < n:
        card = CardEnc.atmost(
            lits=qv, bound=weight, vpool=pool, encoding=EncType.seqcounter
        )
        for cl in card.clauses:
            solver.add_clause(cl)
    sat, model = solver.solve()
    if sat is not True:  # False = UNSAT, None = budget exhausted
        return None
    return np.array([1 if model[v] else 0 for v in qv], dtype=np.uint8)


# ---------------------------------------------------------------------------
# Initial upper bound: stabilizer hill-climb on the X-logical basis


def _greedy_upper_bound(V: np.ndarray, H_X: np.ndarray) -> np.ndarray:
    """Reduce each X-logical basis row by greedily XORing X-stabilizer
    rows while the weight drops (class-preserving, so the result is
    always a genuine nontrivial logical). Returns the lightest vector.
    """
    best = None
    for row in V:
        v = row.copy()
        improved = True
        while improved:
            improved = False
            for s in H_X:
                cand = v ^ s
                if int(cand.sum()) < int(v.sum()):
                    v = cand
                    improved = True
        if best is None or int(v.sum()) < int(best.sum()):
            best = v
    return best


# ---------------------------------------------------------------------------
# Driver


@dataclass(frozen=True, slots=True)
class ShardStat:
    rep: int             # smallest class rep in the task's group
    n_reps: int          # orbit reps covered by this task
    case: str
    weight: int
    status: str          # 'SAT' | 'UNSAT'
    seconds: float


@dataclass(frozen=True, slots=True)
class ShardDistanceResult:
    distance: int
    witness: np.ndarray
    direction: str                     # 'X' (BB symmetry: d_Z = d_X)
    k: int
    num_classes: int                   # 2^k − 1
    num_orbits: int
    orbit_sizes: tuple[int, ...]
    rounds: int
    initial_upper_bound: int
    shard_stats: tuple[ShardStat, ...] = field(repr=False)
    wall_seconds: float = 0.0


def _verify_witness(
    checks: CheckMatrices, L_Z: np.ndarray, v: np.ndarray, d: int
) -> None:
    assert int(v.sum()) == d, "witness weight != claimed distance"
    assert not ((checks.H_Z @ v) % 2).any(), "witness not in ker(H_Z)"
    assert ((L_Z @ v) % 2).any(), "witness is a stabilizer, not a logical"


def shard_distance(
    checks: CheckMatrices,
    *,
    backend: str = "auto",
    jobs: int = 1,
    verbose: bool = False,
    initial_witness: np.ndarray | None = None,
    probe_confl_limit: int = 100_000,
    orbits_per_shard: int | None = None,
    encoding: str = "pin",
) -> ShardDistanceResult:
    """Exact d_X via the shard decomposition (see module docstring).

    `backend`: 'cms' (CryptoMiniSat, native XOR), 'cadical' (pysat
    in-process), or 'auto' (CMS when importable, matching
    `sat_distance`'s no-proof default).
    `jobs > 1` runs each round's shards in a spawn-based process pool,
    terminating the round early when a shard finds a better witness.
    `initial_witness` (optional) seeds the descent, e.g. from QDistRnd
    or a stored table; it must be a nontrivial X-logical.
    `probe_confl_limit`: conflict budget for the monolithic descent
    probes (0 disables probing).
    `orbits_per_shard`: granularity dial. None (default) = coarsest —
    one task per parity class covering all its orbit reps, the
    measured-best "strengthened monolith" configuration. 1 = one task
    per orbit rep (unit pins, most tasks — use for parallel refutation
    or per-shard proof emission). Intermediate values group that many
    reps per selector-encoded task. Exhaustiveness is preserved at any
    value: the groups partition the orbit reps and each task covers
    its reps up to transport.
    `encoding`: 'pin' (parity rows + class pins over qubit variables)
    or 'coset' (coset-leader form over free H_X-row coefficients — no
    feasibility constraints; see `_solve_shard_coset_cms`). 'coset'
    requires CMS.
    """
    t_start = time.perf_counter()
    if backend == "auto":
        backend = "cms" if _HAVE_CMS else "cadical"
    if backend == "cms" and not _HAVE_CMS:
        raise RuntimeError("backend='cms' requested but pycryptosat missing")
    if encoding not in ("pin", "coset"):
        raise ValueError(f"unknown encoding {encoding!r}")
    if encoding == "coset" and not _HAVE_CMS:
        raise RuntimeError("encoding='coset' requires pycryptosat")

    action = compute_class_action(checks)
    L_Z, V, k = action.L_Z, action.V, action.k
    H_X_basis = _row_basis(checks.H_X) if encoding == "coset" else None

    if initial_witness is not None:
        w0 = initial_witness.astype(np.uint8) % 2
        _verify_witness(checks, L_Z, w0, int(w0.sum()))
        witness = w0
    else:
        witness = _greedy_upper_bound(V, checks.H_X)
    d_ub = int(witness.sum())
    initial_ub = d_ub
    if verbose:
        print(
            f"[shard] k={k}  classes={2**k - 1}  orbits="
            f"{len(action.orbit_reps)}  initial d_ub={d_ub}  "
            f"backend={backend}  jobs={jobs}",
            flush=True,
        )

    # Descent probes: cheap budgeted monolithic SAT calls that improve
    # d_ub while improvements are easy to find (see _monolith_probe).
    # The last probe is always a miss (it runs at d − 1 or exhausts its
    # budget), so its cost is pure overhead — keep the budget modest.
    while probe_confl_limit > 0 and _HAVE_CMS and d_ub > 1:
        t0 = time.perf_counter()
        v = _monolith_probe(checks.H_Z, L_Z, d_ub - 1, probe_confl_limit)
        dt = time.perf_counter() - t0
        if v is None:
            if verbose:
                print(
                    f"[shard] probe at w ≤ {d_ub - 1}: no improvement "
                    f"within budget ({dt:.2f}s) → shard rounds",
                    flush=True,
                )
            break
        _verify_witness(checks, L_Z, v, int(v.sum()))
        witness, d_ub = v, int(v.sum())
        if verbose:
            print(
                f"[shard] probe improved d_ub → {d_ub} ({dt:.2f}s)",
                flush=True,
            )

    # refuted[(rep, case)] = highest weight at which this anchored shard
    # was proved empty. UNSAT at weight w implies UNSAT at all w' ≤ w.
    refuted: dict[tuple[int, str], int] = {}
    stats: list[ShardStat] = []
    rounds = 0
    ctx = mp.get_context("spawn")

    last_time: dict[tuple[int, str], float] = {}

    # Coset weight parity. Weight parity is linear over F₂
    # (|a⊕b| ≡ |a|+|b| mod 2), so when every H_X row has even weight
    # (|A|+|B| even — all BB instances here) the parity is constant on
    # each coset and equals parity(v_c) — the solver-side face of the
    # repo's `chainWeight_coset_even`. A rep of parity p cannot hold a
    # vector of weight ≢ p (mod 2), so a refutation bound w tightens to
    # the largest w' ≤ w with w' ≡ p — a strictly easier instance at
    # zero soundness cost.
    hx_even = not any(int(row.sum()) % 2 for row in checks.H_X)
    v_parity = V.sum(axis=1) % 2  # parity of each basis logical
    rep_parity = {
        r: int(_int_to_bits(r, k) @ v_parity) % 2
        for r in action.orbit_reps
    }

    def _tight(w: int, parity: int) -> int:
        return w - ((w - parity) % 2) if hx_even else w

    # Static partition of the orbit reps into task groups —
    # parity-homogeneous so each group gets the tightest sound bound.
    if orbits_per_shard is None:
        orbits_per_shard = len(action.orbit_reps)
    g = max(1, orbits_per_shard)
    groups: list[tuple[int, ...]] = []
    for p in (0, 1):
        block = [r for r in action.orbit_reps if rep_parity[r] == p]
        groups += [
            tuple(block[i:i + g]) for i in range(0, len(block), g)
        ]
    group_transversals = [
        tuple(action.transversal[r] for r in grp) for grp in groups
    ]
    group_of_rep = {
        r: gi for gi, grp in enumerate(groups) for r in grp
    }

    while True:
        w = d_ub - 1
        rounds += 1
        witness_rep = int(
            action.orbit_rep_of[_bits_to_int((L_Z @ witness) % 2)]
        )
        witness_group = group_of_rep[witness_rep]
        tasks = [
            (
                gi, groups[gi], case,
                _tight(w, rep_parity[groups[gi][0]]),
                group_transversals[gi],
            )
            for gi in range(len(groups))
            for case in ("L", "R")
            if refuted.get((gi, case), -1) < w
        ]
        # Order: the current witness's group first (likeliest home of a
        # better vector), then cheap-before-expensive so improvement
        # rounds terminate before sinking time into hard refutations.
        tasks.sort(
            key=lambda t: (
                t[0] != witness_group,
                last_time.get((t[0], t[2]), 0.0),
            )
        )
        if not tasks:  # everything already refuted at ≥ w
            break

        improved = False
        round_stats_start = len(stats)

        def _handle(res: tuple) -> bool:
            """Record one task result; True iff it improved d_ub."""
            nonlocal witness, d_ub
            gidx, reps, case, wt, support, dt = res
            last_time[(gidx, case)] = dt
            if support is None:
                # UNSAT at the parity-tightened bound `wt` refutes all
                # of (wt, w] too: those weights have the wrong parity
                # for this group's cosets. Record against `w`.
                refuted[(gidx, case)] = w
                stats.append(
                    ShardStat(min(reps), len(reps), case, wt, "UNSAT", dt)
                )
                return False
            stats.append(
                ShardStat(min(reps), len(reps), case, wt, "SAT", dt)
            )
            v = np.zeros(action.n, dtype=np.uint8)
            v[list(support)] = 1
            assert int(v.sum()) <= wt
            _verify_witness(checks, L_Z, v, int(v.sum()))
            witness, d_ub = v, int(v.sum())
            return True

        initargs = (checks.H_Z, L_Z, backend, encoding, H_X_basis, V)
        if jobs > 1:
            with ctx.Pool(
                processes=jobs,
                initializer=_init_worker,
                initargs=initargs,
            ) as pool:
                for res in pool.imap_unordered(_run_task, tasks):
                    if _handle(res):
                        improved = True
                        pool.terminate()  # abandon the rest of the round
                        break
        else:
            _init_worker(*initargs)
            for task in tasks:
                if _handle(_run_task(task)):
                    improved = True
                    break

        if verbose:
            n_unsat = sum(
                1 for s in stats[round_stats_start:] if s.status == "UNSAT"
            )
            print(
                f"[shard] round {rounds}: w ≤ {w}, "
                f"{n_unsat}/{len(tasks)} refuted"
                + (f", improved d_ub → {d_ub}" if improved else ""),
                flush=True,
            )
        if not improved:
            break  # full round refuted at w = d_ub − 1  ⟹  d = d_ub

    _verify_witness(checks, L_Z, witness, d_ub)
    return ShardDistanceResult(
        distance=d_ub,
        witness=witness,
        direction="X",
        k=k,
        num_classes=2**k - 1,
        num_orbits=len(action.orbit_reps),
        orbit_sizes=action.orbit_sizes,
        rounds=rounds,
        initial_upper_bound=initial_ub,
        shard_stats=tuple(stats),
        wall_seconds=time.perf_counter() - t_start,
    )
