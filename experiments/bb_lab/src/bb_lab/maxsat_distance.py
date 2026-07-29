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
) -> tuple[CNF, list[int]]:
    """Hard-clause CNF + qubit variable ids (1..n, allocated first)."""
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
        return cnf, qv

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

    return cnf, qv


def write_wcnf(
    checks: CheckMatrices,
    path: Path,
    *,
    mode: str = "strengthened",
    action: ClassAction | None = None,
) -> list[int]:
    """Emit new-format WCNF (h-prefixed hard lines, weight-1 softs on
    ¬v_i). Returns the qubit variable ids."""
    if mode not in ("naive", "strengthened"):
        raise ValueError(f"unknown mode {mode!r}")
    if mode == "strengthened" and action is None:
        action = compute_class_action(checks)
    cnf, qv = _build_hard(checks, mode, action)
    with Path(path).open("w") as f:
        for clause in cnf.clauses:
            f.write("h " + " ".join(str(l) for l in clause) + " 0\n")
        for v in qv:
            f.write(f"1 -{v} 0\n")
    return qv


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
) -> MaxSatDistanceResult:
    """Run a WCNF MaxSAT solver on the chosen encoding and return the
    verified distance. Trust model: the *witness* (weight = optimum
    cost) is re-verified here; the optimality claim is the solver's —
    cross-check against `shard_distance` where independence matters.
    """
    work_dir = Path(work_dir)
    work_dir.mkdir(parents=True, exist_ok=True)
    wcnf = work_dir / f"{mode}_{checks.group.label()}.wcnf"
    action = (
        compute_class_action(checks) if mode == "strengthened" else None
    )
    qv = write_wcnf(checks, wcnf, mode=mode, action=action)

    argv = [str(solver_binary), *extra_args, str(wcnf)]
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
