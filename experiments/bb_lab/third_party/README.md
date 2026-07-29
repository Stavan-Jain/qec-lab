# third_party — patched external solvers

## Tandem (the qec-lab MaxCDCL fork)

**Tandem** — a tandem bicycle: two riders, the code's algebra and
CDCL branch-and-bound, on one (bivariate) bicycle. The analytic layer
verifies instance theorems (coset weight parity, anchor transversals)
and hands them to the solver as search hints; neither rider gets
anywhere alone at this speed.

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
- **`-init-lb=F`** — caller-certified lower bound on the optimum
  (an analytic or kernel-checked distance floor). Routed through the
  solver's existing `initLB → infeasibleUB` machinery: the search
  stops the moment an incumbent reaches the floor, deleting the proof
  phase entirely when the floor is tight. Composes with `-cost-step`
  (the parity shift encodes exactly the excluded gap). SOUNDNESS IS
  THE CALLER'S OBLIGATION — only pass certified floors.
- **`-phase-file=path`** — 0/1 string over DIMACS vars: initial
  branching phases, typically a known witness. Warm-starts the descent
  *without* cold-starting the clause database (the measured defect of
  plain UB seeding). Pure bias; cannot affect correctness.

Measured (2026-07-29, Apple Silicon, naive BB distance WCNFs, under
background load; full table in `scripts/fork_ab.py` output):

| code            | stock | step | step+prime | step+seed |
|-----------------|------:|-----:|-----------:|----------:|
| gross [[144,12,12]] | 5.48 | 1.80 | 1.99 | 3.51 |
| [[150,8,12]]    |  7.80 | 3.22 | 3.42 | 3.36 |
| [[168,6,12]]    |  7.67 | 3.69 | 3.44 | 4.16 |
| [[168,·,14]] ×4 | 34–49 | 16–21 | 13–22 | 14–21 |

`-cost-step=2` alone is a consistent **~2–3×**; priming is ±noise to
mildly positive at n=168; incumbent seeding is neutral at these sizes
(the descent phase warms the clause database, so cold-started proofs
can lose what they save — its value thesis is n≥288). Every distance
is asserted against known values by the harness.
