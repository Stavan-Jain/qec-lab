# FOLD_IN — order-144 sweep as cohort batch 2 (depth-4 stratum)

A background agent is (as of 2026-08-17) running an order-144 sweep
writing to `experiments/bb_lab/data/order144_sweep/`. That directory was
NOT read or written by the prediction session (contamination + shared-
machine rules). When the sweep finishes, its results become the missing
**2-adic depth-4 stratum** of this test via the mechanical procedure
below.

## Procedure

1. Wait until the sweep has fully stopped (no process writing
   `results.jsonl`).
2. From the worktree root:

       cd experiments/bb_lab && \
       uv run python data/descent_theory_test/fold_in.py

   (`--results PATH` if the sweep's output file is named differently;
   `--limit N` for a smoke run.)
3. The script emits `cohort_batch2.jsonl`, `predictions_batch2.jsonl`,
   `aggregates_batch2.json`, and `MANIFEST_batch2.sha256` (its own
   hashes + ISO timestamp) into this directory. Batch 1 files are not
   touched; rows whose instance_id already appears in the batch-1
   cohort are skipped and counted.
4. Check `aggregates_batch2.json` → `depth_regime_occupancy` for the
   depth-4 R1–R4 cells, and record in the Phase-2 notes whether any
   remain empty.

`fold_in.py` was written blind to the sweep's exact schema: it maps the
plausible key spellings (`FIELD_MAP`) and exits with the observed keys
if a required field cannot be mapped — in that case edit `FIELD_MAP`
(a pure renaming) and re-run; do not change any prediction semantics.
The prediction assembly mirrors `tools/make_predictions.py` /
PROTOCOL.md §4 exactly (same frozen route rule, same thresholds, same
wall model, same contamination guards).

## CRITICAL honesty rule — retrodiction vs prediction

The sweep itself computes outcomes for some of its codes. For any row
that **already carries an exact / certificate-tier outcome** (non-null
`d_exact`, or a floor + a method naming a certificate/solver lane), the
outcome PRE-DATES the prediction. Those rows are labeled

    "row_class": "retrodiction"

and they test the screen's **calibration only** — they are NOT
foresight, must NOT be counted toward the pre-registered criteria
(iii)–(iv) closure/agreement rates, and are reported separately.
The screen still runs BLIND on them (the row's own `d_exact` is never
fed into the window or cost machinery); the known outcome is attached
separately as `calibration_truth`, and an immediate `retro_check`
records whether the blind G5 window/ceiling contains it. A
`retro_check.ceiling_contains_d = false` with an all-exact chain is a
G5 falsification datum and must be surfaced, not averaged away.

Only rows with **unknown / bounded-only** outcomes (`"row_class":
"prediction"`) are genuine pre-registered predictions: they join the
batch-1 frontier under criteria (i)–(v) for Phase 2+.

## Scope notes

- Order-144 groups have v₂(|G|) = 4, so towers run 4 rungs / 3 pairs —
  each in-scope row contributes depth-4 regime cells (this cannot fill
  the open depth-2 R4 cell noted in PROTOCOL.md §6.5).
- Odd-axis subtleties (e.g. Z16xZ9: the ladder passes through
  degenerate x-axis orders 2 and 1) are handled by the same ladder
  truncation rules as batch 1; truncations are recorded per row.
- Compute per row is the batch-1 structural screen (~0.1–1 s) plus
  quotient lookups; thousands of rows are fine, but run it `nice`d and
  after the sweep has released the machine.

## Smoke test (done at freeze time)

`fold_in.py` was verified against a SYNTHETIC 4-row results file (in
the session scratchpad, never in the sweep dir): two corpus Z12xZ12
rows with d_exact → labeled retrodiction, blind windows contained the
truth (retro_check pass); one batch-1 cohort dupe → skipped; one
fabricated degenerate Z16xZ9 row (k = 0 at the target) → NO-ROUTE with
reason. The four batch-2 output files from that test were deleted
afterward; the batch-2 files present after the real sweep lands are
the real ones.
