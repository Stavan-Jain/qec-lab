"""Exact BB distance via branch-and-bound MaxSAT (MaxCDCL et al.).

Two WCNF emitters over the same qubit variables (soft: ¬v_i, weight 1,
so the optimum cost IS the distance):

- `naive`: the monolithic encoding QDistSAT benchmarks — H_Z parity
  rows as Tseitin XOR chains, pairing indicators a_j, one clause
  `⋁ a_j`. This reproduces (up to encoding details) the configuration
  under which MaxCDCL solves gross in 48.2 s in arXiv:2606.12445 —
  the same-machine baseline.

- `strengthened`: the shard_distance structure as *hard constraints*:
  classes restricted to the G-orbit reps (selector s_r per rep,
  exactly-covering via `⋁ s_r`; pins `s_r → a_j = rep_j`), per-rep
  support anchors (`s_r → support meets T_r` — the union-over-blocks
  clause, single-instance form of the L/R split: any shard member has
  a Stab(r)-translate whose support meets T_r in whichever block is
  nonempty), and, when every class has even weight parity (all H_X
  rows even ⟹ parity is a coset invariant), a global even-weight XOR
  row. The MaxSAT optimum over this space is still exactly d: reps
  cover all classes up to weight-preserving transport, anchors
  preserve per-shard minima, and parity removes nothing of minimal
  weight.

The point (2026-07-29): the 48 s MaxCDCL result is solver-level
progress on the *naive* encoding; shard_distance's ~3× is
encoding-level progress under a plain SAT backend. This module stacks
them. Every reported witness is independently re-verified.

`maxsat_distance_descent` (2026-07-29, experiment 1 of the descent→
Tandem transfer): the *orchestration-level* port of the descent-layer
sector dichotomy (`bb_lab.descent_sat`). For a verified axis deck σ,

    d = min( 2 · opt_a , opt_b )

where (a) is a **base-sized** MaxSAT over the quotient code (the
σ-invariant sector: v = p* v̄, |v| = 2|v̄|, nontriviality = ⋁ of the
Pμ pairing outputs) and (b) is the cover naive WCNF plus the sector
restriction p₊ v ≠ 0 (fiber-pair XOR definitions + one wide clause),
optionally strengthened by the implied base parity rows and the
Λ-transport link μ̄(w) = Λ·μ(v) — expressible without selectors
because the naive encoding's a_j literals *are* μ(v). Hard-constraint
scaffolding measurably taxes BnB (the `strengthened` lesson above),
so which (b) variant wins — and whether the decomposition beats plain
naive at all — is exactly what the A/B measures. Soundness of the
dichotomy is inherited from `compute_descent`'s numeric battery; both
sector optima are the solver's word, both witnesses re-verified here.
"""

from __future__ import annotations

import subprocess
import time
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from pysat.formula import CNF, IDPool

from .checks import CheckMatrices
from .sat_distance import _xor_chain, find_logical_z
from .shard_distance import (
    ClassAction,
    _bits_to_int,
    _int_to_bits,
    compute_class_action,
)


@dataclass(frozen=True, slots=True)
class MaxSatDistanceResult:
    distance: int
    witness: np.ndarray | None   # None iff seeded run proved the seed
    mode: str                 # 'naive' | 'strengthened'
    solver_seconds: float
    optimum_found: bool


def _build_hard(
    checks: CheckMatrices,
    mode: str,
    action: ClassAction | None,
) -> tuple[CNF, list[int], list[int]]:
    """Hard-clause CNF + qubit variable ids (1..n, allocated first) +
    the pairing-indicator literals a_j (μ(v) bit j ⟺ a_j true)."""
    H_Z = checks.H_Z
    n = H_Z.shape[1]
    N = n // 2
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

    L_Z = find_logical_z(checks) if action is None else action.L_Z
    a = []
    for L in L_Z:
        idx = np.flatnonzero(L)
        assert idx.size, "zero row in the logical-Z basis"
        a.append(_xor_chain((qv[i] for i in idx), pool, cnf))

    if mode == "naive":
        cnf.append(list(a))
        return cnf, qv, a

    assert action is not None
    k = action.k
    sel = [pool.id() for _ in action.orbit_reps]
    cnf.append(sel)
    for s, rep in zip(sel, action.orbit_reps):
        cbits = _int_to_bits(rep, k)
        for j in range(k):
            cnf.append([-s, a[j] if cbits[j] else -a[j]])
        T = action.transversal[rep]
        cnf.append(
            [-s] + [qv[p] for p in T] + [qv[N + p] for p in T]
        )

    # Global even-weight row when parity is a coset invariant and every
    # class is even (the corpus norm; toric-like codes skip this).
    hx_even = not any(int(row.sum()) % 2 for row in checks.H_X)
    v_parity = action.V.sum(axis=1) % 2
    all_even = hx_even and not any(
        int(_int_to_bits(r, k) @ v_parity) % 2 for r in action.orbit_reps
    )
    if all_even:
        out = _xor_chain(iter(qv), pool, cnf)
        cnf.append([-out])

    return cnf, qv, a


def write_wcnf(
    checks: CheckMatrices,
    path: Path,
    *,
    mode: str = "strengthened",
    action: ClassAction | None = None,
) -> tuple[list[int], list[int]]:
    """Emit new-format WCNF (h-prefixed hard lines, weight-1 softs on
    ¬v_i). Returns (qubit variable ids, pairing-indicator ids)."""
    if mode not in ("naive", "strengthened"):
        raise ValueError(f"unknown mode {mode!r}")
    if mode == "strengthened" and action is None:
        action = compute_class_action(checks)
    cnf, qv, a = _build_hard(checks, mode, action)
    with Path(path).open("w") as f:
        for clause in cnf.clauses:
            f.write("h " + " ".join(str(l) for l in clause) + " 0\n")
        for v in qv:
            f.write(f"1 -{v} 0\n")
    return qv, a


def emit_fiber_certificate(
    checks: CheckMatrices,
    A,
    B,
    sigma: tuple[int, ...],
    qv: list[int],
    a_lits: list[int],
    path: Path,
    floor_cap: int = 20,
    floor_budget: int | None = 100_000,
) -> Path:
    """Write the `-fiber-lb` certificate for the naive WCNF whose
    variable ids are (qv, a_lits): the σ fiber pairing, the invariant-
    sector floor (exact, witness-jumping base SAT), and the per-class
    moving-sector table fiberFloorB[c] = base-coset floor of Λ·c.
    Every structural premise is asserted inside `compute_descent`;
    floors are certified by the base-code SAT engine —
    conflict-budgeted ascending by default (sound at any stopping
    point; heavy classes get the floor proven within budget), exact
    witness-jumping with `floor_budget=None`. The file is the
    caller-verified half of the Tandem `-fiber-lb` contract."""
    from .descent_sat import (
        base_coset_floors,
        base_coset_floors_budgeted,
        compute_descent,
        invariant_floor,
    )

    action = compute_class_action(checks)
    dd = compute_descent(checks, action, sigma, A=A, B=B)
    k = action.k
    lam_of = {
        c: _bits_to_int((_int_to_bits(c, k) @ dd.Lam) % 2)
        for c in range(1, 1 << k)
    }
    if floor_budget is None:
        floors = base_coset_floors(
            dd, set(lam_of.values()), cap=floor_cap
        )
    else:
        floors = base_coset_floors_budgeted(
            dd, set(lam_of.values()), cap=floor_cap,
            confl_budget=floor_budget,
        )
    inv = invariant_floor(dd)
    pairs = [
        tuple(int(qv[i]) for i in np.flatnonzero(dd.S[b]))
        for b in range(dd.S.shape[0])
    ]
    assert all(len(p) == 2 for p in pairs)
    with Path(path).open("w") as f:
        f.write(f"p fiberlb {k} {len(pairs)} {inv}\n")
        f.write("a " + " ".join(str(int(x)) for x in a_lits) + "\n")
        tbl = [0] * (1 << k)
        for c in range(1, 1 << k):
            tbl[c] = floors[lam_of[c]]
        f.write("f " + " ".join(map(str, tbl)) + "\n")
        for x, y in pairs:
            f.write(f"{x} {y}\n")
    return Path(path)


def _parse_v_bits(stdout: str) -> str:
    bits = []
    for line in stdout.splitlines():
        if line.startswith("v "):
            bits.append(line[2:].strip().replace(" ", ""))
    return "".join(bits)


def maxsat_distance(
    checks: CheckMatrices,
    solver_binary: Path | str,
    *,
    mode: str = "strengthened",
    work_dir: Path | str = ".",
    timeout: float | None = None,
    extra_args: tuple[str, ...] = (),
    seed_ub: int | None = None,
    init_lb: int | None = None,
    phase_bits: np.ndarray | None = None,
    fiber_sigma: tuple[int, ...] | None = None,
    fiber_polys: tuple | None = None,
    fiber_reuse: bool = True,
) -> MaxSatDistanceResult:
    """Run a WCNF MaxSAT solver on the chosen encoding and return the
    verified distance. Trust model: the *witness* (weight = optimum
    cost) is re-verified here; the optimality claim is the solver's —
    cross-check against `shard_distance` where independence matters.

    `fiber_sigma` (naive mode only) emits the `-fiber-lb` certificate
    for that verified axis deck and passes it to the fork;
    `fiber_polys` = (A, B) is then required. `fiber_reuse=True` reuses
    an already-emitted certificate file (its ids are deterministic for
    a given instance).
    """
    work_dir = Path(work_dir)
    work_dir.mkdir(parents=True, exist_ok=True)
    wcnf = work_dir / f"{mode}_{checks.group.label()}.wcnf"
    action = (
        compute_class_action(checks) if mode == "strengthened" else None
    )
    qv, a_lits = write_wcnf(checks, wcnf, mode=mode, action=action)

    argv = [str(solver_binary), *extra_args]
    if fiber_sigma is not None:
        assert mode == "naive", "fiber-lb targets the naive encoding"
        assert fiber_polys is not None, "fiber_sigma needs fiber_polys=(A,B)"
        sig = "".join(str(int(x)) for x in fiber_sigma)
        flb = work_dir / f"fiber_{checks.group.label()}_{sig}.flb"
        if not (fiber_reuse and flb.exists()):
            emit_fiber_certificate(
                checks, fiber_polys[0], fiber_polys[1], fiber_sigma,
                qv, a_lits, flb,
            )
        argv.append(f"-fiber-lb={flb}")
    if init_lb is not None:
        # CALLER-CERTIFIED floor (analytic/kernel-checked). The solver
        # stops when an incumbent reaches it — a wrong floor yields a
        # wrong answer, so only pass certified values.
        argv.append(f"-init-lb={init_lb}")
    if phase_bits is not None:
        pf = work_dir / f"{mode}_{checks.group.label()}.phases"
        pf.write_text("".join("1" if b else "0" for b in phase_bits))
        argv.append(f"-phase-file={pf}")
    argv.append(str(wcnf))
    if seed_ub is not None:
        argv.append(str(seed_ub))  # positional initial-UB (stock feature)
    t0 = time.perf_counter()
    proc = subprocess.run(
        argv, capture_output=True, text=True, timeout=timeout,
    )
    dt = time.perf_counter() - t0

    out = proc.stdout
    optimum = "s OPTIMUM FOUND" in out
    if not optimum:
        raise RuntimeError(
            f"solver did not report OPTIMUM (exit {proc.returncode});"
            f" tail:\n{out[-1500:]}\n{proc.stderr[-500:]}"
        )
    cost = None
    for line in out.splitlines():
        if line.startswith("o "):
            cost = int(line[2:].strip())
    bits = _parse_v_bits(out)

    if len(bits) < max(qv):
        # No (usable) model printed. Legitimate exactly when the run
        # was seeded and the solver proved nothing better than the
        # seed exists: optimum = seed, and the CALLER holds the
        # weight-`seed_ub` witness that justified the seed.
        if seed_ub is not None and cost == seed_ub:
            return MaxSatDistanceResult(
                distance=int(cost), witness=None, mode=mode,
                solver_seconds=dt, optimum_found=True,
            )
        raise RuntimeError(
            f"solver printed OPTIMUM but no usable v-line "
            f"(cost={cost}, seed={seed_ub}); tail:\n{out[-1500:]}"
        )

    v = np.array(
        [1 if bits[q - 1] == "1" else 0 for q in qv], dtype=np.uint8
    )

    # Independent verification of the witness.
    assert int(v.sum()) == cost, "v-line weight != reported optimum"
    assert not ((checks.H_Z @ v) % 2).any(), "witness not in ker(H_Z)"
    L_Z = find_logical_z(checks)
    assert ((L_Z @ v) % 2).any(), "witness is a stabilizer"

    return MaxSatDistanceResult(
        distance=int(cost),
        witness=v,
        mode=mode,
        solver_seconds=dt,
        optimum_found=True,
    )


# ---------------------------------------------------------------------------
# Experiment 1: descent sector decomposition at the orchestration level


@dataclass(frozen=True, slots=True)
class DescentMaxSatResult:
    distance: int
    witness: np.ndarray
    sector: str                  # 'a' (invariant) | 'b' (moving)
    opt_a: int | None            # cover weight 2·|v̄| (None = sector empty)
    opt_b: int | None
    seconds_a: float
    seconds_b: float
    variant: str                 # 'min' | 'rows' | 'full'

    @property
    def solver_seconds(self) -> float:
        return self.seconds_a + self.seconds_b


def _run_wcnf(
    solver_binary: Path | str, wcnf: Path,
    extra_args: tuple[str, ...], timeout: float | None,
) -> tuple[int | None, str, float]:
    """Run the solver; return (optimum_cost | None-if-hard-UNSAT,
    stdout, seconds)."""
    t0 = time.perf_counter()
    proc = subprocess.run(
        [str(solver_binary), *extra_args, str(wcnf)],
        capture_output=True, text=True, timeout=timeout,
    )
    dt = time.perf_counter() - t0
    out = proc.stdout
    if "s UNSATISFIABLE" in out:
        return None, out, dt
    if "s OPTIMUM FOUND" not in out:
        raise RuntimeError(
            f"solver reported neither OPTIMUM nor UNSATISFIABLE "
            f"(exit {proc.returncode}); tail:\n{out[-1200:]}"
        )
    cost = None
    for line in out.splitlines():
        if line.startswith("o "):
            cost = int(line[2:].strip())
    assert cost is not None, "OPTIMUM without an o-line"
    return cost, out, dt


def _model_bits(out: str, qv: list[int]) -> np.ndarray:
    bits = _parse_v_bits(out)
    assert len(bits) >= max(qv), "no usable v-line in solver output"
    return np.array(
        [1 if bits[q - 1] == "1" else 0 for q in qv], dtype=np.uint8
    )


def _write_wcnf_lines(path: Path, cnf: CNF, soft: list[int]) -> None:
    with Path(path).open("w") as f:
        for clause in cnf.clauses:
            f.write("h " + " ".join(str(l) for l in clause) + " 0\n")
        for v in soft:
            f.write(f"1 -{v} 0\n")


def maxsat_distance_descent(
    checks: CheckMatrices,
    A,
    B,
    solver_binary: Path | str,
    *,
    sigma: tuple[int, ...] | None = None,
    variant: str = "min",
    work_dir: Path | str = ".",
    timeout: float | None = None,
    cost_step: bool = True,
    timing_only_naive_check=None,
) -> DescentMaxSatResult:
    """d_X = min(2·opt_a, opt_b) via the sector dichotomy (module
    docstring). `variant` controls sector (b)'s scaffolding: 'min' =
    fiber definitions + sector clause only; 'rows' = + implied base
    parity rows; 'full' = + the Λ-transport link. `cost_step=True`
    passes -cost-step=2 to each solve whose parity premise holds
    (verified here per instance, per side)."""
    from .descent_sat import axis_decks, compute_descent
    from .linalg import nullspace_f2, quotient_complement_basis

    if variant not in ("min", "rows", "full"):
        raise ValueError(f"unknown variant {variant!r}")
    work_dir = Path(work_dir)
    work_dir.mkdir(parents=True, exist_ok=True)
    action = compute_class_action(checks)
    if sigma is None:
        sigma = axis_decks(checks)[0]
    dd = compute_descent(checks, action, sigma, A=A, B=B)
    n = checks.num_qubits
    nb = dd.S.shape[0]
    label = checks.group.label()

    # --- sector (a): base-sized instance -------------------------------
    pool_a = IDPool()
    qa = [pool_a.id() for _ in range(nb)]
    cnf_a = CNF()
    for row in dd.base_checks.H_Z:
        idx = np.flatnonzero(row)
        if idx.size:
            out = _xor_chain((qa[i] for i in idx), pool_a, cnf_a)
            if out is not None:
                cnf_a.append([-out])
    pouts = []
    for j in range(action.k):
        idx = np.flatnonzero(dd.Pmu[j])
        if idx.size == 0:
            continue
        pouts.append(_xor_chain((qa[i] for i in idx), pool_a, cnf_a))
    if not pouts:
        opt_a, out_a, dt_a = None, "", 0.0   # no invariant logicals at all
    else:
        cnf_a.append(pouts)
        wcnf_a = work_dir / f"descent_a_{label}.wcnf"
        _write_wcnf_lines(wcnf_a, cnf_a, qa)
        # parity gate for the base side: every base cycle even ⟹ all
        # feasible costs share parity.
        extra_a: tuple[str, ...] = ()
        if cost_step and not any(
            int(r.sum()) % 2 for r in nullspace_f2(dd.base_checks.H_Z)
        ):
            extra_a = ("-cost-step=2",)
        opt_a, out_a, dt_a = _run_wcnf(
            solver_binary, wcnf_a, extra_a, timeout
        )

    wit_a = None
    if opt_a is not None:
        vb = _model_bits(out_a, qa)
        assert int(vb.sum()) == opt_a
        wit_a = (dd.P_lift @ vb) % 2
        L_Z = find_logical_z(checks)
        assert int(wit_a.sum()) == 2 * opt_a
        assert not ((checks.H_Z @ wit_a) % 2).any()
        assert ((L_Z @ wit_a) % 2).any(), "sector-(a) model is a stabilizer"

    # --- sector (b): cover naive + sector restriction ------------------
    pool_b = IDPool()
    qv = [pool_b.id() for _ in range(n)]
    cnf_b = CNF()
    for row in checks.H_Z:
        idx = np.flatnonzero(row)
        if idx.size:
            out = _xor_chain((qv[i] for i in idx), pool_b, cnf_b)
            if out is not None:
                cnf_b.append([-out])
    L_Z = find_logical_z(checks)
    a_lits = []
    for L in L_Z:
        idx = np.flatnonzero(L)
        a_lits.append(_xor_chain((qv[i] for i in idx), pool_b, cnf_b))
    cnf_b.append(list(a_lits))
    # fiber-pair definitions w_b := ⊕(fiber), sector clause ⋁ w.
    wv = []
    for b in range(nb):
        idx = np.flatnonzero(dd.S[b])
        wv.append(_xor_chain((qv[i] for i in idx), pool_b, cnf_b))
    cnf_b.append(list(wv))
    if variant in ("rows", "full"):
        for row in dd.base_checks.H_Z:
            idx = np.flatnonzero(row)
            if idx.size:
                out = _xor_chain((wv[i] for i in idx), pool_b, cnf_b)
                if out is not None:
                    cnf_b.append([-out])
    if variant == "full":
        # Λ-link: ⟨L̄_Z[j̄], w⟩ ⊕ Σ_i Λ[i,j̄]·a_i = 0 (μ̄(p₊v) = Λ·μ(v);
        # premise verified inside compute_descent).
        for jb in range(dd.kb):
            lits = [wv[i] for i in np.flatnonzero(dd.Lb_Z[jb])]
            lits += [
                a_lits[i] for i in range(action.k) if dd.Lam[i, jb]
            ]
            out = _xor_chain(iter(lits), pool_b, cnf_b)
            if out is not None:
                cnf_b.append([-out])
    wcnf_b = work_dir / f"descent_b{variant}_{label}.wcnf"
    _write_wcnf_lines(wcnf_b, cnf_b, qv)
    extra_b: tuple[str, ...] = ()
    if cost_step:
        hx_even = not any(int(r.sum()) % 2 for r in checks.H_X)
        V = quotient_complement_basis(checks.H_X, nullspace_f2(checks.H_Z))
        if hx_even and not any(int(v.sum()) % 2 for v in V):
            extra_b = ("-cost-step=2",)
    opt_b, out_b, dt_b = _run_wcnf(solver_binary, wcnf_b, extra_b, timeout)

    wit_b = None
    if opt_b is not None:
        wit_b = _model_bits(out_b, qv)
        assert int(wit_b.sum()) == opt_b
        assert not ((checks.H_Z @ wit_b) % 2).any()
        assert ((L_Z @ wit_b) % 2).any(), "sector-(b) model is a stabilizer"
        assert ((dd.S @ wit_b) % 2).any(), "sector-(b) model is invariant"

    cands = []
    if opt_a is not None:
        cands.append((2 * opt_a, "a", wit_a))
    if opt_b is not None:
        cands.append((opt_b, "b", wit_b))
    assert cands, "both sectors hard-UNSAT: the code has no logicals?"
    d, sector, witness = min(cands, key=lambda t: t[0])
    return DescentMaxSatResult(
        distance=d, witness=witness, sector=sector,
        opt_a=None if opt_a is None else 2 * opt_a,
        opt_b=opt_b, seconds_a=dt_a, seconds_b=dt_b, variant=variant,
    )
