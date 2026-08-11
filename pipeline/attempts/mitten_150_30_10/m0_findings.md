# M0 findings — mitten_150_30_10 (A32)

**Date:** 2026-08-10 · **Tool:** `experiments/bb_lab/scripts/a32_m150_lean_scoping.py`
· **Data:** `experiments/bb_lab/data-a32/m150_scoping.json` · **Matrices:**
sha-pinned shipped bytes (census) / canonical Table XIII rebuild
(symmetry; A26 §5 established the Tanner isomorphism between the two).

## F1 — Symmetry: the code's automorphism group is EXACTLY C₅ (nauty-definitive)

- Exhaustive blockwise-L/R translation search (all 30 group elements ×
  2⁹ block assignments): exactly the 4 non-identity **central** elements
  act (uniform translation, all-L = all-R on the center). Center of
  C₅×S₃ = C₅.
- Aut(G) has order 24 (= C₄×S₃, matches theory); **only the identity**
  preserves the four base sets setwise. No Aut(G)-induced symmetries.
- **pynauty on the colored 270-vertex Tanner graph: |Aut| = 5.0e0 = 5,
  one generator.** There are no exotic symmetries of any kind. The
  orbit-reduction factor available to the M4 floor sweeps is exactly
  **5** (gross had 72; per-|G| this code is 6× poorer in transport).

## F2 — No X↔Z duality (nauty + census ⟹ airtight)

- Natural family (grid-slot permutation × antipode on/off per sector,
  96 candidates): none maps rows(H_X) onto rows(H_Z).
- pynauty certificates: colored-Tanner-graph iso (H_X,H_Z) ≅ (H_Z,H_X)
  is **False** — no qubit permutation bijects check supports across
  sides.
- Once F3 (census) confirms the weight-≤9 rowspace words are exactly the
  rows, no-duality is airtight at rowspace level too: a permutation
  carrying rowspace(H_X) to rowspace(H_Z) must biject their unique
  weight-9 words — the rows — hence be a Tanner iso, excluded above.
- **Consequence: M4 proves two independent floors.** (Consistent with
  the asymmetric canonical logical weights 18 vs 10.)

## F3 — Exact light-kernel census (the floor statement's ground truth)

Statement under test: *the only nonzero weight-≤9 vectors of ker H_X are
the 60 rows of H_Z, and mirror* (then "light cycle ⟹ generator row" is
the exact M4 target, A21's shape one notch up).

Pre-census evidence: ISD over 6,000 information sets per side found no
non-row word of weight ≤ 9; pair sums ≥ 14; sampled triples ≥ 15.

CMS-XOR blocking-clause enumeration at wt ≤ 9 (incremental solver,
XOR checks + seqcounter atmost-9 + per-solution blocking, final call
UNSAT ⟹ exhaustive):

- ker H_X: **60 solutions in 374.2 s (exhausted), all weight 9, = the
  60 H_Z rows exactly** (0 non-row extras, 0 rows missed)
- ker H_Z: **60 solutions in 179.6 s (exhausted), all weight 9, = the
  60 H_X rows exactly**

**VERDICT: `rows_only = True` — both sides.** There is no weight-≤8
kernel vector at all, and nothing at weight 9 beyond the generator
rows. The M4 floor statements are, verbatim:

```
∀ v ∈ ker H_X, 0 < |v| ≤ 9 → v is a row of H_Z     (⟹ v ∈ rowspace H_Z)
∀ v ∈ ker H_Z, 0 < |v| ≤ 9 → v is a row of H_X
```

## F4 — Elimination structure (from the scoping session, re-checked)

- rank L(a1) = rank R(b1) = rank R(b0) = 30 (invertible); L(a0)/R(a0)
  rank 28. On ker H_X the two check equations solve v₃, v₄ from
  (v₁, v₂, v₅); mirror side analogous. Kernel free rank 90 ✓ (= 150−60).
- ε-parity (all entries odd weight): |v₁|+|v₃|+|v₅| ≡ 0,
  |v₂|+|v₄|+|v₅| ≡ 0 (mod 2) ⟹ |v| ≡ |v₅| (mod 2).

## F5 — GAP → Lean group dictionary validated (`leandict` mode)

Explicit C₅×S₃ parameterization of `group_30_1.txt`: z = index 2
(central, order 5), r = index 3 (order 3), s = index 1 (order 2,
s·r·s = r⁻¹); every index = z^i·r^j·s^k uniquely, and the map is a
**homomorphism validated on all 30×30 products** against the GAP table.
Sets in (i, j, k) coordinates (for the M2 emitter; identity ∈ a1, b1 ✓):

```
a0 = {(0,0,0), (1,2,0), (2,1,1)}    a1 = {(0,0,0), (1,0,0), (0,1,1)}
b0 = {(1,1,0), (2,2,0), (4,1,0)}    b1 = {(0,0,0), (1,0,0), (4,1,1)}
```

The mathlib orientation (`r j` vs `sr j`, `DihedralGroup 3` mul
convention) is deliberately NOT baked in here — the emitter re-validates
it with an in-Lean `decide` at generation time.

## Consequences for the plan

1. M4 statement shape confirmed — classification onto a *named finite
   list* (60 rows/side), certificate-friendly.
2. Transport budget: normalize only mod C₅ (factor 5, not 30). Split-map
   cells ≈ 6× a naive |G|-normalized estimate; reinforces the
   certificate-first budget rules (raw sweeps stay offline).
3. Both sides proven independently; keep leaf files side-disjoint so
   lake parallelism recovers most of the 2×.
4. The census wall time (374 s + 180 s for exhaustive ≤9 reasoning by a
   SAT solver) is a scale datum: brute weight-≤9 reasoning at n = 150
   costs minutes even for CMS — reinforcing the budget rule that
   enumeration lives offline and the Lean build only checks
   certificates.

## F6 — Budget rehearsal (M0.4): DONE, all obligations ≈ 9.2 s

`scripts/m150_gen_lean_data.py rehearsal` emitted a self-contained probe
(`A32Rehearsal.lean`, imports Dihedral + ZMod + Fintype.Prod only) with
four representative obligations — dictionary hom (900 products,
native_decide), H_Z rows ∈ ker H_X + weight 9 via **closed-form H
entries** (540k-term sums, native_decide), Gaussian left-inverse pivot
certificate W·H_X[:,P] = I₆₀ (packed-Nat testBit, native_decide), and a
900-pair kernel `decide`. Compiled clean on the FIRST attempt via
`lake env lean` in the QECLean checkout: **12.4 s wall warm vs 3.2 s
import-only baseline ⟹ ≈ 9.2 s for everything.** Numbers + per-shape
table in `build_budget.md`.

Carrier note: the multiplicative carrier is
`Multiplicative (ZMod 5) × DihedralGroup 3` (plain `ZMod 5` is
additive); mathlib orientation resolved and in-Lean-validated:
z^i·r^j ↦ `(ofAdd i, .r j)`, z^i·r^j·s ↦ `(ofAdd i, .sr (−j))`.
The closed-form entry formulas (no stored H tables needed):

```
HX (β,h) (m,x) = [m=4] · [x⁻¹h ∈ b_β]  +  [m<4 ∧ m%2=β] · [hx⁻¹ ∈ a_{m/2}]
HZ (α,h) (m,x) = [m=4] · [xh⁻¹ ∈ a_α]  +  [m<4 ∧ m/2=α] · [h⁻¹x ∈ b_{m%2}]
```

validated entry-for-entry against `a26.mitten_code` before emission.

**M0 is COMPLETE.** Next: M1 (`Framework/Homological/LiftedProduct.lean`
in a QECLean worktree) ∥ growing the emitter into the M2 data modules.
