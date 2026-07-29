# third_party — patched external solvers

## MaxCDCL (qec-lab fork)

Upstream: MaxCDCL, Li–Coll–Li (MIT license), MaxSAT Evaluation 2023
source artifact —
<https://maxsat-evaluations.github.io/2023/mse23-solver-src/exact/MaxCDCL.zip>.
We do not vendor the source; `build_maxcdcl.sh` fetches the public zip,
builds the pristine baseline (`maxcdcl_stock`), applies
`maxcdcl-qeclab.patch` (~40 changed lines across `core/Solver.h`,
`core/Solver.cc`, `simp/Main.cc`) and builds the fork
(`maxcdcl_release`). Builds clean on macOS arm64 and Linux.

The patch adds two *caller-verified, instance-specific* hints; with
neither flag the fork is behaviourally identical to stock:

- **`-cost-step=N`** — declares that all feasible costs are congruent
  mod N. After each improving model (cost printed first), the exclusive
  search bound is tightened by N−1; the seeded-UB path (stock's
  positional `<wcnf> <UB>` argument) is tightened the same way.
  Soundness is the caller's obligation: for BB distance WCNFs
  (soft ¬v_i), N=2 is exactly the coset weight-parity theorem
  (`chainWeight_coset_even`) — every H_X row even ⟹ weight parity is a
  coset invariant; all-classes-even is verified per instance by
  `bb_lab.scripts.fork_ab._all_even` before the flag is passed.
- **`-prime-vars=a,b,…`** — 1-based DIMACS vars to branch on first
  with positive phase (bumped in both VSIDS and CHB heaps). Used for
  the translation-symmetry anchor qubits (L0, R0). A search *bias*,
  never a constraint — cannot affect correctness.

Measured (2026-07-29, Apple Silicon, naive BB distance WCNFs, under
background load): `-cost-step=2` gives a consistent ~2–3× speedup
(gross 5.4→1.9 s; [[168,·,14]] instances ~39→13–18 s); priming adds a
little at n≥168; incumbent seeding is neutral at n≤168 (the descent
phase warms the clause database, so cold-started proofs can lose what
they save). See `scripts/fork_ab.py` for the harness that produced
these numbers and asserts every distance against known values.
