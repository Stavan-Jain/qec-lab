# Post-freeze amendments

Amendments to files covered by `MANIFEST.sha256` (frozen
2026-08-17T17:27:38Z). Batch-1 *content* files (cohort, predictions,
aggregates, protocol, logs) are untouched — `shasum -c` passes on every
manifest entry except the one listed here.

## #1 — 2026-08-17 ~13:30 ET: degenerate-case fix in `tools/dtt_lib.py`

**Trigger.** First real run of `fold_in.py` (order-144 sweep, 58 rows)
crashed in `TS._preimage` (scripts/a32_tower_slice.py:624) via
`dtt_lib.screen_structure`:
`ValueError: matmul ... (size 4 is different from 0)`.

**Cause.** When the descended seam image `W = p_bot*(S)` spans the FULL
bottom H1 space, the annihilator of `span(Wb)` is empty;
`np.array(kernel_basis(Wmat))` on an empty list is 1-D of size 0, and
`F @ M` fails. Mathematically the preimage of the full space is the
whole domain. The case never arises in batch 1 (183 candidates
screened) nor in the A32/A33/A35 production runs — their W-spaces are
proper subspaces; several depth-4 sweep towers hit it.

**Fix.** Added `_preimage_safe()` in `tools/dtt_lib.py` and swapped the
single call site (`screen_structure`, pair loop). The helper is
byte-for-byte the `TS._preimage` computation (same
`TS.i2v`/`TS.v2i`/`TS.kernel_basis` conventions), plus the degenerate
branch: empty annihilator ⟹ return the standard basis of the domain.
No non-degenerate computation changes; no prediction semantics change;
batch-1 outputs were not regenerated.

**Hashes.**
- `tools/dtt_lib.py` frozen:
  `8372a4235ae37d1a3f6092fba6db782bc427aa5e7da387eb4cbbabd577f24106`
- `tools/dtt_lib.py` amended:
  `5c45da60947e4cfcc526592dce83791f19f7046a5b5e395c398e0fd4070f1fb4`

**Upstream.** The latent bug remains in
`scripts/a32_tower_slice.py::_preimage` (committed program code, not
edited by this data errand); flagged separately for a proper fix +
regression test in the main repo.
