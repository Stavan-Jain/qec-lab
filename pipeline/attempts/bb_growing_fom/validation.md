# Validation record

2026-09-13. All required checks passed in the isolated QECLean worktree
on branch `codex/bb-growing-fom`, based on `9e7d121`.

The new modules have zero compiler/linter diagnostics under the repository
settings. Existing modules emitted their pre-existing style warnings.
The full axiom-policy check reports:

```text
axiom policy OK — 7024 declarations across 190 modules depend on only [propext, Classical.choice, Quot.sound]
```

| Check | Result | Seconds |
|---|---|---:|
| `lake build QEC.Stabilizer.Codes.BivariateBicycle.Gross.SafeFloor.MImFloorY0` | Passed before the remaining suite | — |
| `lake build QEC.Stabilizer.Codes.BivariateBicycle.Gross.SafeFloor.MImFloorY1` | Passed | 94.44 |
| `lake build QEC.Stabilizer.Codes.BivariateBicycle.Gross.SafeFloor.MImFloorY4` | Passed | 103.05 |
| `lake build` | Passed | 494.22 |
| `lake build QECLight` | Passed | 6.54 |
| `lake build QECBlueprint` | Passed | 15.62 |
| `lake build QECWidgets` | Passed | 24.71 |
| `lake env lean Playground.lean` | Passed | 4.46 |
| `lake env lean scripts/AxiomCheck.lean` | Passed | 4.82 |
| `bash scripts/check-umbrellas.sh` | Passed | 1.7 |

The three memory-heavy existing gross leaves were built sequentially.
No dependency versions changed, and no mathlib cache download or rebuild
was requested. Logs remain at `/private/tmp/bb-validation/`.

The finite seed checker was also run successfully: its default checks
exhaustively determine distances through `s=4`, and two independent native
runs exhaustively checked `s=5`. These experiments supply no premise to
the universal Lean proofs.
