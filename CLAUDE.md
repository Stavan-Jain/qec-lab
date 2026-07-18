# CLAUDE.md — agent orientation (qec-lab)

This is the **research workbench** for
[QECLean](https://github.com/Stavan-Jain/QECLean), the Lean 4 / mathlib QEC
library. This repo contains **no Lean library code**: the BB research program
(`experiments/bb_lab/`), the formalization pipeline (`pipeline/`, `catalog/`),
the narrative docs (`docs/`), and the dashboard (`dashboard/`).

## Two-repo layout (read this first)

- **QECLean** — the library: `QEC/`, lakefile, toolchain. Its `main` is
  sorry-free and mathlib-grade. All Lean editing, `lake build`, the lean-lsp
  MCP workflow, and formalization worktrees happen **in a QECLean checkout**,
  not here. Read **QECLean's `CLAUDE.md`** before any Lean work — naming/style
  conventions, layering rules, linter policy, build discipline, and the
  require-permission lake operations all live there.
- **qec-lab** (this repo) — everything else. Tooling here locates the library
  via the **sibling-checkout convention**: env `QECLEAN_ROOT`, default
  `../QECLean` relative to this repo's root.
- **Upstreaming flow**: WIP formalization branches are created in the QECLean
  checkout. Back them up to this repo's remote with `git push lab <branch>`
  (the repos share pre-split git ancestry, so branches move freely). When
  sorry-free and polished, PR the branch to QECLean `main`.
- The full pre-split history (research + library) lives in this repo's
  `git log` — provenance queries about retired approaches, deleted branches,
  and archive tags belong here.

## Where to look for what

- **[`docs/pipeline-usage.md`](docs/pipeline-usage.md)** — task-oriented
  recipes: weekly triage, start a new code, review a skeleton, run Stage 4,
  open a PR, refresh the catalog, initialize a moonshot. Start here to *do*
  something.
- **[`docs/pipeline.md`](docs/pipeline.md)** — pipeline architecture: stages,
  scoring rubric, artifact contents.
- **[`docs/lean-patterns.md`](docs/lean-patterns.md)**,
  **[`docs/mathlib-version-quirks.md`](docs/mathlib-version-quirks.md)**,
  **[`docs/lean-conversion-recipes.md`](docs/lean-conversion-recipes.md)** —
  Lean tactical references consulted during formalization work (the work
  itself happens in the QECLean checkout).
- **[`pipeline/research_log.md`](pipeline/research_log.md)** — index of
  moonshot attempts (failures are first-class outputs).
- **[`experiments/bb_lab/notes/A_HANDOFF.md`](experiments/bb_lab/notes/A_HANDOFF.md)**
  — canonical Phase-A handoff; §0 is the resume-here paragraph. Approach
  numbers are claimed in `experiments/bb_lab/notes/README.md` (the registry).
- **[`experiments/bb_lab/GENERATORS.md`](experiments/bb_lab/GENERATORS.md)** —
  the Python generators that emit Lean into the QECLean checkout. Generated
  files there carry `GENERATED FILE — DO NOT HAND-EDIT` banners: change the
  generator here and regenerate, landing both repos' changes together.

## Formalization pipeline

Key artifacts at a glance:

- `catalog/zoo.yaml` — 267 quantum codes from the Error Correction Zoo
- `catalog/scoring.yaml` — per-code formalization-priority scores
- `pipeline/queue.md` — top-of-queue + tracks (engineering / moonshot /
  defer / skip)
- `pipeline/attempts/<code_id>/` — per-code formalization state
  (`state.yaml`, `informal_spec.md`, `plan.md`, `result.md`, …)
- `.claude/agents/qec-{prioritizer,skeleton-drafter,formalization-runner,moonshot}.md`
  — operating specs for the four pipeline agents. Read when modifying agent
  behavior, not for day-to-day use.

Stage-2/4 agents draft and prove Lean **in the sibling QECLean checkout**
while tracking attempt state here — spawn them with both checkouts available.

## Python environment

`experiments/bb_lab/` is a uv project (`pyproject.toml`, `uv.lock`). Run lab
scripts as `uv run --project experiments/bb_lab python ...` from the repo
root, or `cd experiments/bb_lab && uv run python scripts/...`. The generators
refuse to run if the QECLean sibling checkout is missing.

`pipeline/cache/eczoo_data/` is a gitignored snapshot of the EC Zoo data,
reproducible from the SHA pinned in `pipeline/cache/PIN.md`
(`scripts/ingest_zoo.py` rebuilds `catalog/zoo.yaml` from it).

## Dashboard

`dashboard/build.py` renders the static site to `dashboard/dist/` from
`catalog/` + `pipeline/`; the recent-activity panel additionally reads git
history from the QECLean sibling when present (degrades gracefully without
it). Deployed to GitHub Pages by `.github/workflows/dashboard.yml` on pushes
touching `dashboard/`, `catalog/`, or `pipeline/`. See `dashboard/README.md`.

## Conventions carried over from the pre-split repo

- Parallel agent sessions run in worktrees — now **QECLean worktrees** (under
  the QECLean checkout's `.claude/worktrees/<name>/`; its CLAUDE.md documents
  the mathlib-sharing symlink trick). A worktree's checked-out branch name
  usually differs from its directory name — check `git branch --show-current`.
- Use `git cherry main <branch>` (patch-aware), not raw ahead/behind counts,
  to judge merge state.
- Sorry markers in draft Lean: `sorry  -- TODO(<tag>): <one-line note>` so the
  next session can grep for them.
