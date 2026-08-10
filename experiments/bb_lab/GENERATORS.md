# bb_lab → QEC generators

Python scripts that emit (parts of) tracked Lean files under
`QEC/Stabilizer/Codes/BivariateBicycle/` in the sibling QECLean checkout
(env var `QECLEAN_ROOT`, default `../QECLean` relative to the qec-lab repo
root). The Lean-side inventory (which
files are generated, class G/F/H taxonomy) lives in
`QECLean:QEC/Stabilizer/Codes/BivariateBicycle/README.md`; this file is the
operational side.

**Environment**: `cd experiments/bb_lab && uv sync` (Python 3.11, deps in
`pyproject.toml`; dev group for pytest). Run generators as
`uv run python <script> [args]`.

**Clobber policy**: every generator that writes into the QECLean checkout's
`QEC/` must (a) emit the
`GENERATED FILE — DO NOT HAND-EDIT` banner naming itself (as a
`qec-lab:experiments/...` path), its data source and
the regen command, and (b) refuse to overwrite an existing target without
`--force` (pattern: `gen_pair72_packaging_data.py`). A hand-edit found in a
Class-G file is a bug: move the change into the generator and land the
generator change (qec-lab) and the regenerated output (QECLean) together,
as cross-referencing PRs.

## Live generators

| Script | Emits | Notes |
|---|---|---|
| `phase5/gen_file.py` | `Gross/StabilizerCodeData.lean` | rewritten 2026-07-18: emits ONLY the §1 data module (7 defs, line-wrapped ≤100 chars) with banner + `--force` guard; the proof file `Gross/StabilizerCode.lean` is hand-maintained and never touched |
| `scripts/gen_floor_lean.py` | `Gross/SafeFloor/MImFloorData.lean` | cost tables `D3V`/`RCELL` + Γ data |
| `scripts/gen_yrep_module.py <i>` | `Gross/SafeFloor/MImFloorY<i>.lean` | i ∈ {0,1,4} only — Y11/Y12 diverged (hand-evolved analytic Tier-3 form, QECLean PR #58); the script refuses 11/12 without `--force` |
| `scripts/gen_assembly_2d.py` | fragments for `Gross/SafeFloor/MImAssembly.lean` | Class F: paste between the `BEGIN/END GENERATED` markers only |
| `scripts/gen_pair72_packaging_data.py` | `Z3Z6/StabilizerCodeData.lean` | 15-check ALL-PASS validation gate + `--force` guard; NOTE: the Lean-side split is done but the generator still emits the pre-split full file — retarget to data-only before the next regen |
| `scripts/gen_base_floor_lean.py` | `BaseFloors/<Name>.lean` | class-member small-cycle bundles (A15/T2) |
| `scripts/gen_f2a6_z5z30_data.py` | data feed for `Z5Z15F2A6/` | A17 line |
| `scripts/m150_gen_lean_data.py` | A32 mitten [[150,30,10]]: `rehearsal` mode = budget probe (scratchpad); `instance` mode = **`Codes/Mitten/M150/Data.lean`** (Class G: GAP dictionary, sets, pivot/packed-inverse tables both sides, 30+30 symplectic logical supports, weight-10 witness pair) | falsify-first: dictionary hom + closed-form H vs `a26.mitten_code` + pivot inverses + symplectic conditions + witness all validated in numpy before emission; carrier `Multiplicative (ZMod 5) × DihedralGroup 3`; `Defs.lean` is hand-maintained |

## Retired (`scripts/attic/`)

| Script | Why retired |
|---|---|
| `attic/gen_orbit_module.py` | emitted `MImFloorO*` modules that no longer exist (superseded by the Y-representative transport + `gen_yrep_module.py`) |
| `attic/gen_assembly.py` | emitted the 13-orbit `MImAssembly` (superseded by the 5-orbit 2-D dispatch, `gen_assembly_2d.py`) |
