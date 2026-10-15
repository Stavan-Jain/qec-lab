# qec-lab

[![Dashboard](https://img.shields.io/badge/pipeline%20dashboard-live-success)](https://stavan-jain.github.io/qec-lab/)
[![Library](https://img.shields.io/badge/Lean%20library-QECLean-blueviolet)](https://github.com/Stavan-Jain/QECLean)
[![License: Apache 2.0](https://img.shields.io/badge/License-Apache_2.0-blue.svg)](LICENSE)

The research workbench behind [**QECLean**](https://github.com/Stavan-Jain/QECLean),
a Lean 4 / mathlib formalization of quantum error correction. QECLean holds the
polished, sorry-free library; this repo holds everything it took to get there —
and everything in flight:

- **`experiments/bb_lab/`** — the bivariate-bicycle (BB) code distance program:
  Python lab (`src/`, `scripts/`), SAT/linear-algebra certificates, phase
  write-ups, and the approach notes (`notes/A<N>_*.md`). The flagship outcome is
  the fully analytic, Lean-verified `d = 12` proof for IBM's `[[144,12,12]]`
  gross code.
- **`pipeline/`** — the catalog-driven formalization pipeline: per-code attempt
  state (`attempts/<code_id>/`), the priority queue (`queue.md`), and the
  moonshot research log (`research_log.md`, where failure write-ups are
  first-class outputs).
- **`catalog/`** — 267 quantum codes ingested from the
  [Error Correction Zoo](https://errorcorrectionzoo.org/) (`zoo.yaml`) with
  formalization-priority scores (`scoring.yaml`).
- **`docs/`** — the human-readable proofs and guides: the
  [gross-code distance proof](docs/gross-distance-proof.md) (+
  [LaTeX](docs/gross-distance-proof.tex)), the
  [doubling-template extensibility program](docs/gross-distance-extensibility.md),
  the [toric distance proof](docs/distance_proof.md), the
  [Steane distance-3 explainer](docs/steane-distance-proof.md), Lean pattern references,
  and the [pipeline architecture](docs/pipeline.md) /
  [operator's manual](docs/pipeline-usage.md).
- **`dashboard/`** — the static generator behind the
  [live dashboard](https://stavan-jain.github.io/qec-lab/) tracking all 267
  catalog codes (done / in-flight / queued / deferred) and moonshot attempts.
- **`.claude/agents/`** — operating specs for the four pipeline agents
  (prioritizer, skeleton-drafter, formalization-runner, moonshot).

## Research results

Start with the [research assessment](docs/research-status.md) for the strongest
results, their significance, proof assumptions, and remaining gaps. The main
outcomes are:

- **Gross distance and cover theory:** `[[144,12,12]]` has an unconditional
  Lean proof on QECLean `main`, using only its three standard axioms. The
  reusable doubling theorem has explicit sufficient conditions.
- **Descent certificates:** exact distances for `[[360,12,24]]`,
  `[[288,8,20]]`, and `[[432,12,24]]`, plus 97 new exact distances in the
  recorded A39 study. These computational certificates are distinct from
  Lean-checked endpoints.
- **A weight-six fractal family:** a formal lower bound and unbounded actual
  `k d²/n`, on a QECLean research branch. Its closed distance formula remains
  conjectural; this is not a hardware-performance claim.
- **Further results:** a non-abelian mitten distance proof on a branch that
  uses native decision procedures, an all-size IBM logical-dimension
  correction, and analytic exact-distance theorems for modified BB families.

The original IBM all-size distance conjecture remains open. Some older
handoffs report formal instances now parked off `main`; the assessment
distinguishes those from the current library. The preferred next direction
is a reusable descent engine with independently checkable certificates.

For local storage and reproduction, see
[repository maintenance](docs/repository-maintenance.md).

## Relationship to QECLean

The two repos share git ancestry (qec-lab was seeded from the pre-split QECLean
history, so all research provenance is in `git log` here) and work as a pair:

- **Sibling checkouts.** Tooling here expects the library at `../QECLean`
  (override with the `QECLEAN_ROOT` env var). The Lean generators under
  `experiments/bb_lab/` emit files into that checkout; the dashboard reads its
  git history for the recent-activity panel.
- **Upstreaming.** Work-in-progress formalization branches live in the QECLean
  checkout and can be backed up to this repo's remote (shared ancestry makes
  `git push lab <branch>` work). When a formalization is sorry-free and
  polished, the branch is PR'd to QECLean `main` — the mathlib/cslib model.
- **Generated files.** Lean files in QECLean carrying a
  `GENERATED FILE — DO NOT HAND-EDIT` banner are emitted by generators here;
  see [`experiments/bb_lab/GENERATORS.md`](experiments/bb_lab/GENERATORS.md).
  Change the generator and regenerate — landing both repos' changes together.

## Reading order

1. [`docs/research-status.md`](docs/research-status.md) — current assessment,
   evidence tiers, and recommended next work.
2. [`pipeline/research_log.md`](pipeline/research_log.md) — index of all
   moonshot approach lines (A3…), one rich paragraph each.
3. [`experiments/bb_lab/notes/A_HANDOFF.md`](experiments/bb_lab/notes/A_HANDOFF.md)
   — the canonical Phase-A handoff; §0 is the resume-here paragraph.
4. Per-approach notes `experiments/bb_lab/notes/A<N>_*.md` — approach numbers
   are claimed in [`experiments/bb_lab/notes/README.md`](experiments/bb_lab/notes/README.md).
5. [`docs/gross-distance-extensibility.md`](docs/gross-distance-extensibility.md)
   — the generalization program beyond the gross code.

## License

Released under the [Apache License 2.0](LICENSE).

## Maintainer

Maintained by Stavan Jain — a project from the
**University of Wisconsin–Madison**.
