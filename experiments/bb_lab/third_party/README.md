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

The build script ends by asserting that `tandem` advertises all five
fork flags and that `maxcdcl_stock` advertises none of them. Two ways
it used to produce a mislabelled pair, both now fixed and both silent:

- patch-path arithmetic. Apple's `/usr/bin/patch` skips every hunk
  whose stripped path doesn't exist relative to its cwd (GNU patch's
  fuzzier prefix search hides the mismatch), and in the pre-guard era
  the build then continued and copied the *unpatched* binary to
  `tandem`. It bit twice: the original patch carried `a/code/core/...`
  paths and was applied with `-p1 -d code` (→ `code/code/core/...`,
  matches nothing); the regenerated v5 patch carries `a/core/...`
  paths and was applied with plain `-p1` from the MaxCDCL root
  (→ `MaxCDCL/core/...` — the sources live under `code/`). It is
  `-p1 -d code` today, and a costStep grep hard-fails the build if
  application is ever skipped again.
- on a *rerun*, `unzip -o` restores pristine sources carrying their 2023
  archive timestamps, which are older than the release objects (`*.or`)
  left by the previous patched build. make judged those stale objects up
  to date and relinked them into `maxcdcl_stock`, so the "pristine
  baseline" was itself the fork. Both builds now `make clean` first.

Neither corrupts the table below: with no flags passed the fork is
behaviourally identical to stock, so a `stock` column accidentally
measured on a fork binary measures the same thing. They matter for
reproducibility, and would matter for correctness the moment the patch
stops being inert by default.

The patch adds five *caller-verified, instance-specific* hints; with
none of the flags passed the fork is behaviourally identical to stock:

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
- **`-fiber-lb=path`** — the per-fiber lower-bound rounding transplant
  (descent theory → BnB bound arithmetic; no clauses, no aux vars).
  The certificate file carries a free-deck fiber pairing of the soft
  vars, the class-indicator (a-var) ids, a per-class table of
  certified base-coset floors, and the invariant-sector floor. At
  each search node the solver computes, from fiber states alone,
  min over the two completion types (σ-invariant / moving) of a
  certified completion cost; reaching the incumbent fires a soft
  conflict through the `UBconflictFlag`/`involvedLits` premise hook,
  so the learned clause is valid *given the caller's theorem*.
  Emitted by `bb_lab.maxsat_distance.emit_fiber_certificate` (all
  structural premises asserted numerically per instance; floors
  certified by base-code SAT). SOUNDNESS IS THE CALLER'S OBLIGATION —
  a wrong table yields wrong answers (demonstrated in the toy
  falsification test).

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
