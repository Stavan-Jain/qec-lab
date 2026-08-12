import QEC.Stabilizer.Codes.BivariateBicycle.Bb288.Defs
import QEC.Stabilizer.Codes.BivariateBicycle.Bb288.Witness
import QEC.Stabilizer.Codes.BivariateBicycle.Bb288.Distance

/-!
# The `[[288,12,18]]` record code — instance umbrella

`bb288` over `Z₁₂ × Z₁₂` (`A = x³ + y² + y⁷`, `B = y³ + x + x²`;
arXiv:2308.07915 Table 3) as a **twisted** free ℤ₂ cover of the gross
code (`y⁷ ≡ y mod 6`), through the target-floor (deficit-rung) assembly
of `Framework/Homological/BBDoubling.lean` — the first instance with
`d(cover) < 2·d(base)` (`18 < 24`) and the first whose base floor is a
kernel-checked library theorem (`gross_chain_distance_eq_12`).

- `Defs.lean`     — cover polynomials, `bb288Complex`, the
                    `XDoubleCoverData G288 GrossGroup` bundle (base pair
                    = the library's `grossA/grossB`)
- `Witness.lean`  — the weight-18 nontrivial cycle + dual witness
                    (`d ≤ 18`, kernel-checked)
- `Distance.lean` — `LogicalFloor 12` from the library +
                    `bb288_{chain,pauli}_distance_eq_18`, conditional on
                    the two A36 sector-floor certificates (the
                    `Z5Z15F2A6` two-tier shape; provenance in the module
                    docstring)
-/
