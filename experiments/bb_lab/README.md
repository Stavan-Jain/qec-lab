# bb-lab — Bivariate-Bicycle code laboratory

Python substrate for the QEC moonshot pipeline. Computes properties of
bivariate-bicycle (BB) quantum LDPC codes outside Lean, so that conjecture
search can happen before formalization.

## Status

See the [research assessment](../../docs/research-status.md) for evidence
tiers, exact sources, researcher significance, and current limitations.

- **Gross `[[144,12,12]]`:** unconditional distance theorem on QECLean
  `main`, using only `propext`, `Classical.choice`, and `Quot.sound`.
  Its proof uses the reusable doubling layer and kernel-checked finite
  arguments; it does not await an LRAT importer.
- **Descent and tower certificates:** `[[360,12,24]]`, `[[288,8,20]]`,
  and `[[432,12,24]]` have end-to-end computational distance certificates.
  The A39 study records 97 further exact distances. These Python/C
  certificate results are not automatically Lean theorems.
- **Research branches:** the weight-six fractal family has a formal
  lower bound and unbounded actual `k d²/n`; mitten `[[150,30,10]]` has
  a formal endpoint using native decision procedures. Some older BB
  instances are parked off `main` pending removal of those procedures.
- **Analytic work:** original IBM homology/dimension results and exact
  modified families are proved in the notes. The unchanged IBM all-size
  distance formula and the fractal closed distance formula remain open.

The solver lane remains useful. CaDiCaL emits DRAT proofs, `drat-trim`
checks them independently, and the certificate JSON hashes each proof and
its DIMACS input. The compact `bb_72_12_6` bundle is retained; large proof
outputs are regenerable with `bb-lab bravyi-check --emit-proofs`.
The original v0 substrate and proposed LRAT route are historical context,
not the current gross-proof status. Storage policy and cleanup details are
in [repository maintenance](../../docs/repository-maintenance.md).

## Quick start

```bash
cd experiments/bb_lab
uv sync --extra dev
uv run pytest                                # default test selection
uv run pytest -m slow                        # optional long-running solver tests
uv run pytest -o addopts=''                  # everything including slow
uv run bb-lab bravyi-check --quick           # 3 small SAT distances
uv run bb-lab bravyi-check --quick --emit-proofs   # +DRAT proofs
uv run bb-lab bravyi-check --full --emit-proofs    # includes expensive larger SAT instances
uv run bb-lab verify-cert certificates/bb_72_12_6.cert.json   # full re-check

# Build drat-trim once (used by the verify-cert test):
( cd /tmp && git clone --depth 1 https://github.com/marijnheule/drat-trim.git \
  && cd drat-trim && cc -O2 -o drat-trim drat-trim.c )
```

## The browser UI

```bash
./third_party/build_maxcdcl.sh      # once: fetch + build maxcdcl_stock and tandem
uv run bb-lab ui                    # http://127.0.0.1:8765
```

Type a group and two polynomials. The default view is just the answer:
`[[n, k, d]]`, the check weight, and one line of secondary facts (group,
|G|, check counts, qubit degree, rate, CSS guard). **Compute distance**
fills in `d` — incumbents stream live, so a long run shows its current
upper bound rather than a spinner, and Cancel kills the solver process.

Everything else is one click away. **Solver options** holds the backend,
the encoding and the solver flags, with a summary line
(`Tandem · naive · -cost-step=2`) so you can see the selection without
opening it. **Full report** holds the ranks, weight profiles, canonical
polynomial forms and the solver premises.

`--binary <path>` (or `$BB_LAB_TANDEM`) points at a different solver
build; without one the UI falls back to the in-process CryptoMiniSat
ladder and says so. The Tandem run mirrors the sweep scripts exactly —
`mode="naive"`, `-cost-step=2` only when the coset weight-parity premise
is verified for that code — and reports the `d_method` string the run
would earn in the corpus.

The solver's options are read from `tandem --help` at startup rather than
hardcoded, so a flag added to the fork shows up in the UI as soon as you
rebuild. Only flags carrying a *soundness obligation* need a line of
Python: an entry in `webui/solver.py: FLAG_NOTES` naming the premise from
`webui/analysis.py: premises()` that must hold before the flag may be
passed. That check runs server-side on every request.

## Module layout

| Path | Purpose |
| --- | --- |
| `src/bb_lab/group.py` | Finite abelian group `ZMod ℓ × ZMod m` arithmetic |
| `src/bb_lab/poly.py` | F₂-polynomial parsing (`'x^3 + y + y^2'`) and canonicalization |
| `src/bb_lab/checks.py` | `H_X`, `H_Z` construction. **Single source of truth vs the Lean `conv` definition.** |
| `src/bb_lab/linalg.py` | Dense F₂ rank, nullspace, quotient-complement |
| `src/bb_lab/codeparams.py` | `n`, `k` from check matrices |
| `src/bb_lab/sat_distance.py` | Exact distance via SAT; dual backend (pysat fast / cadical CLI for proofs) |
| `src/bb_lab/maxsat_distance.py` | WCNF emitters (`naive` / `strengthened`) + the Tandem driver |
| `src/bb_lab/doubling_certify.py` | Doubling certification front end and hypothesis assembly |
| `src/bb_lab/tower.py` | Quotient-tower distance machinery |
| `src/bb_lab/cosetbz.py` | Complete coset censuses from information windows |
| `src/bb_lab/fibering.py` | Odd-fiber coordinates and safe-floor computations |
| `src/bb_lab/webui/` | `bb-lab ui`: analysis, flag discovery + premise gating, stdlib HTTP/SSE, static page |
| `src/bb_lab/certificate.py` | JSON witness + DRAT-reference certificate format `bb-cert/v2` |
| `src/bb_lab/lean_bridge.py` | `state.yaml` ↔ JSON ↔ `BBChainComplex.lean` skeleton |
| `src/bb_lab/store.py` | DuckDB schema (for v1 corpus) |
| `src/bb_lab/cli.py` | `bb-lab` entry points |
| `instances/bravyi_table.yaml` | Regression contract for the five published BB codes |
| `certificates/<code_id>.cert.json` | Per-instance witness + DRAT-hash certificate |
| `certificates/<code_id>/` | DRAT + CNF proof bundle (committed for bb_72_12_6 only) |

## The Lean handshake

The Lab agrees with `QEC/Stabilizer/Framework/Homological/BBChainComplex.lean`
in the sibling QECLean checkout (`QECLEAN_ROOT`, default `../QECLean`)
on the convolution convention:

> `conv A f g = ∑_h A(h) · f(g − h)` over F₂[G]

which means the check matrix is `H_X[g, h] = A(g − h)` for the X-block and
`B(g − h)` for the Z-block. **If a future change makes `bb_lab` and the
Lean side disagree on this, `tests/test_gross_agreement.py` is the
canary.**

The full bridge: `pipeline/attempts/<id>/state.yaml` (polynomial strings) →
`bb_lab.lean_bridge.descriptor_from_state_yaml` → JSON
(`schema_version: bb-instance/v1`) → `emit_skeleton` → compilable
Lean file importing `BBChainComplex`. See `tests/test_lean_roundtrip.py`
for the round-trip property test (compiles the emitted file under
`lake env lean` in the QECLean checkout).

## CI gate

```
pytest experiments/bb_lab/tests/
```

The foundational checks include the suites below; the current test tree
also covers the later certification and solver machinery. Some tests need
external solvers or the sibling library, so inspect skips as well as passes.

1. `test_gross_agreement` — Lab `H_X` for gross matches the Lean truth
   bitwise.
2. `test_bravyi_quick` — `(n, k)` matches the published table for all
   five instances.
3. `test_bravyi_sat` — SAT distance matches `{6, 10, 10}` for the three
   smallest instances.
4. `test_lean_roundtrip` — `gross` state.yaml → Lab → emitted `.lean`
   compiles cleanly under `lake env lean` (requires the sibling QECLean
   checkout).
5. `test_drat_emission` — `cadical` CLI subprocess emits well-formed
   DRAT + CNF pairs; certificate JSON records SHA256-stable references;
   **`drat-trim` independently verifies every proof** (skipped if the
   `drat-trim` binary is missing).
6. `test_verify_cert` — the committed `bb_72_12_6.cert.json` re-verifies
   end-to-end: matrix hashes match the rederived (G, A, B), witness
   passes the three-condition logical check, and every committed DRAT
   proof passes `drat-trim`. This is the v1-track (b) soundness gate.

The `bb-lab bravyi-check --full` larger SAT runs are **not in CI**. Run
them deliberately when regenerating that solver lane's artifacts; the
existence of a later analytic or descent proof does not make this SAT
workload cheap.

## Historical DRAT-to-LRAT proposal

The original v1 solver-verification roadmap proposed a Lean-side LRAT
consumer. The following describes that optional route, not unfinished
work required by the current gross theorem:

- Read `certificates/<code_id>.cert.json` — get `n`, `direction`,
  `witness_support`, `h_check_sha256`, `logical_space_sha256`,
  `unsat_proofs: [{weight_bound, drat_path, cnf_path, sha256s}]`.
- Reconstruct `H_check` and `L_logical` from the same `(G, A, B)` Lean
  already has (`bbChainComplex`). Hashes match → the certificate is
  about *our* code, not someone else's. Note `logical_space_sha256`
  covers the *space* `span(H_check ∪ L_logical) = ker(H_X)` in reduced
  row-echelon form, so Lean is free to pick its own logical basis —
  only the code has to agree.
- For each `UnsatProofRef`: convert DRAT→LRAT with `drat-trim -L`,
  then run an in-kernel LRAT verifier (Carneiro-style), conclude no
  logical of weight ≤ `weight_bound` exists.
- Combine with the witness model (in the certificate) for the SAT
  direction at weight `distance`. Together: a Lean-checked proof of
  exact distance.

Mario Carneiro's [LRAT-check project](https://github.com/digama0/mm0)
was a reference for the historical proposal. No LRAT integration or global
priority claim is asserted here. The current recommended work is a reusable
descent certificate format with an independent checker, as described in
the [research assessment](../../docs/research-status.md#recommended-next-research-direction).
