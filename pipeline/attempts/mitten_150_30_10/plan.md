# Formalization plan: \([[150,30,10]]\) mitten code

## Strategy summary

Six layers M0–M5. The framework and packaging layers (M1–M3, M5) are
ports of proven patterns (generic `HomologicalCode`, gross Phase-5
packaging minus drop sets, cover300 assembly); the research core is M4 —
the two chain-level floors d_X ≥ 10 and d_Z ≥ 10. M0 measurements put
the floors in the best shape the program knows: each floor is an
A21-style *light-cycle ⟹ generator-row* classification (cap 9, five
blocks, C₅-only transport, both sides separately).

**Two-tier discipline** (cover300 pattern): land M1–M3 + a *conditional*
assembly first — `mitten150_hasCodeDistance_10 (hX : FloorX) (hZ :
FloorZ)` with the floors as named hypotheses whose provenance is the
existing bb-cert/v2 SAT certificates — then discharge the floors. The
bundled `StabilizerCodeWithDistance 150 30 10` inhabitant exists only
once M4 closes.

**Gating constraint: total added build time ≤ 5 min (hard cap 10).**
Every layer below carries a budget line; `build_budget.md` is the ledger
and the PR gate. The budget effectively *mandates* the
certificate-first architecture: enumeration happens offline in Python;
Lean checks rank/pivot certificates and small batched decidable facts.
No raw 10⁸-mask sweeps in the final artifact (A15 history: 53 min as
filter-form masks → 3.9 s as Gaussian-pivot certificates).

## M0 — lab-side scoping (no Lean, no build cost) — IN PROGRESS

1. **Exact weight-≤9 kernel census, both sides** (SAT enumeration with
   blocking clauses over the CMS XOR encoding). Confirms/refutes the
   floor statement "≤ 9 ⟹ one of the 60 rows". *If extra words exist
   they join the classification target list.*
2. **True Tanner Aut group** (nauty, colored 270-vertex graph) + X↔Z
   iso test. Decides the real orbit-reduction factor (≥ 5) and whether
   a duality halves M4.
3. **Emitter skeleton** `m150_gen_lean_data.py` per GENERATORS.md:
   group dictionary (GAP index → `ZMod 5 × DihedralGroup 3`), sets,
   H/L/Lx/Lz tables (packed-Nat), pivot certificates, falsify-first
   validation of every emitted fact.
4. **Budget rehearsal**: prototype one representative certificate file
   and one batched `decide`/`native_decide` leaf in a QECLean worktree;
   measure wall cost; extrapolate before committing to the M4 design.

Exit criteria: census result recorded in `m0_findings.md`; measured
per-leaf costs recorded in `build_budget.md`.

## M1 — framework: `Framework/Homological/LiftedProduct.lean`

Non-abelian analog of `BBChainComplex.lean`:

- `lconv a f = fun g => ∑ h, a h * f (h⁻¹ * g)` and
  `rconv f b = fun g => ∑ h, f h * b (h⁻¹ * g)`-style right version
  (exact convention fixed against the J1 matrices by the emitter);
  bilinearity, `lconv_rconv_comm` (associativity — the chain law).
- `mittenChainComplex (A : Fin 2 → G → ZMod 2) (B : Fin 2 → G → ZMod 2) :
  HomologicalCode` with C2 = Fin 2 × G, C1 = Fin 5 × G, C0 = Fin 2 × G,
  boundary maps per Eq. (J1).
- Layering: file sits beside `BBChainComplex.lean`, imports only
  `Homological/Distance.lean` upstream; no BB imports.

Budget: pure-proof file, ≤ 10 s. Est. 1 session.

## M2 — instance + `StabilizerCode 150 30`

`QEC/Stabilizer/Codes/Mitten/M150/{Defs,Data,StabilizerCode}.lean`
(+ umbrella `Mitten.lean`, README per repo conventions):

- `Defs.lean`: `MittenGroup := ZMod 5 × DihedralGroup 3`; the four sets
  as explicit element lists (emitter-generated); `m150Complex :=
  mittenChainComplex ..`; small sanity `decide`s (set weights,
  identity ∈ a1, b1).
- `Data.lean` (generated): packed-Nat tables for H rows, Lx/Lz, the
  full-rank pivot certificates, decoder left-inverses.
- `StabilizerCode.lean`: gross Phase-5 pattern, *simplified*: both
  matrices full rank ⟹ no drop sets, no kernel corrections. Obligations:
  (1) closure equality (trivial: kept = all), (2) generator independence
  via decoder identities (φ_X·∂₂ = id on faces; 60×150 scale — gross's
  72×144 analog ran ~5 s per identity as one `native_decide`), (3) the
  30 logical qubits from Lx/Lz with `Lx·Lzᵀ = I₃₀`, (4) assembly →
  `m150StabilizerCode : StabilizerCode 150 30`.

Budget: ≤ 60 s (two decoder identities ~5–10 s each + table elaboration
+ symplectic-basis checks batched). Est. 1–2 sessions.

## M3 — witness (d ≤ 10)

`Witness.lean`: `uStar := Lz row 0` (weight 10, cycle), `zStar := Lx
row 0`; `⟨uStar, zStar⟩ = 1` + `zStar ∈ ker H_Z` ⟹ `uStar ∉
rowspace(H_Z)`. All kernel `decide` on packed vectors.

Budget: ≤ 10 s. Bundled with M2.

## M4 — the floors (research core)

Target statements (chain level, per side; shown for the ker-H_X side):

```
∀ v ∈ ker H_X, 0 < |v| ≤ 9 → v ∈ {the 60 rows of H_Z}
```

(M0 census permitting; extra census members get named and classified
the same way.) Architecture, in kill order:

1. **ε-parity layer** (pure proof, zero build cost): the two block
   parities + |v| ≡ |v₅| (mod 2) kill a slab of the (w₁..w₅) partition
   table outright.
2. **Elimination layer**: invertible L(a1) (resp. R(b0), R(b1) on the
   mirror side) solves v₃, v₄ from (v₁, v₂, v₅) — Lean-side as
   *certified* left-inverse tables (rank facts, not computation).
   Weight-transfer lemmas bound |v₃|, |v₄| below by rank-style facts
   about short products (the non-abelian small-weight theory: left-coset
   Sidon/difference arguments; `SidonConvBound` groundwork is
   abelian-stated but the two-sided difference-set half is
   frame-agnostic per the A11/A16 line).
3. **Split map**: case on the block-support pattern (2⁵, minus
   ε-killed, minus wlog under C₅ + any M0-discovered Aut) then on
   partitions of ≤ 9. Per case: either an analytic kill (weight
   arithmetic), or
4. **Certificate leaves**: per-case Gaussian-pivot / window-rank
   certificates (KernelCert pattern) emitted and falsify-first-validated
   offline; residual genuinely-finite checks as *batched* packed-Nat
   `native_decide` leaves, budget-capped (see ledger).

The 60 target rows appear in the (3,3,0,0,3)/(0,0,3,3,3)-type cells as
the classified survivors — same shape as A21's `sweepA3` "fires at
exactly the generator column" leaves.

Fallback ladder if a case resists certificate form within budget:
(a) more offline math (turn it into a rank fact); (b) restructure the
window decomposition; (c) spend the 5→10 min fallback *only here*;
(d) if still over: the case's sweep stays a named hypothesis backed by
its SAT certificate and the unconditional bundle waits — the hard cap
is never exceeded.

Budget: ≤ 3 min for both sides combined (target); this is the layer the
whole budget exists for. Est. 3–6 sessions, real risk of new theory
being needed for 1–3 stubborn cases. Parallel side bet (moonshot track,
separate note): F₂[C₅×S₃] ≅ F₂[S₃] × F₁₆[S₃]-component fibering along
the central C₅ (Maschke on the odd part) — an A22-mechanism port that
would replace layer 3–4 sweeps with two small component classifications.

## M5 — assembly + bundle

Transport the two chain floors + witness through the CSS bridge to
`HasCodeDistance m150StabilizerCode 10` (gross `Distance.lean` +
cover300 assembly pattern), then

```
def mitten150StabilizerCodeWithDistance : StabilizerCodeWithDistance 150 30 10 :=
  { m150StabilizerCode with hasDistance := ... }
```

Axiom audit: `#print axioms` = 3 std + `native_decide` only (allowed
per QECLean CLAUDE.md); no `ofNat`/instance drift; linter-clean.
Budget: ≤ 30 s. Est. ≤ 1 session.

## Build-budget governance

- Ledger: `build_budget.md` (this dir) — per-file measured wall cost on
  the dev machine, updated every session that touches Lean.
- Measurement protocol: `touch` the new module files, `time lake build`
  the instance umbrella, record wall + CPU; full-library delta measured
  once per phase.
- Design rules (binding): certificate-first; batch decidable facts
  (few `native_decide` invocations, table-driven ∀-forms); packed-Nat
  tables (`docs/lean-patterns.md` recipes); one heavy leaf per file
  only if it parallelizes; no `Finset`-per-mask in sweeps; prefer
  kernel `decide` for sub-second facts (avoids per-call compiler
  overhead).

## Risks

| Risk | Exposure | Mitigation |
|---|---|---|
| M4 case machine needs new non-abelian small-weight theory | schedule (sessions), not budget | A16/A21 analytic playbook; moonshot fibering side bet |
| A leaf only closes as a big sweep (> budget) | build budget | fallback ladder above; hypothesis + SAT cert keeps cap intact |
| Census finds non-row light words | statement shape | they're finitely many; classify alongside rows |
| `DihedralGroup 3` decide performance vs `ZMod` products | M2/M4 leaf costs | budget rehearsal in M0; fallback = table-driven mul on `Fin 30` |
| Elaboration cost of 150-col tables | M2 | packed-Nat literals (proven: 24 min → 0.4 s pattern) |

## Milestone order

M0 → (M1 ∥ emitter) → M2+M3 → conditional assembly (tier 1 shippable)
→ M4-Z, M4-X → M5 (tier 2 = the bundle). Lean work in a dedicated
QECLean worktree per repo conventions; attempt state tracked here.
