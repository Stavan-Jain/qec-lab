# Informal spec: \([[150,30,10]]\) mitten code

## Summary

The smallest **mitten code** of Bhardwaj–Ma–Meister–King–Bluvstein–
Preskill–Cain–Xu–Huang, *High-rate qLDPC processors* (arXiv:2607.28795):
a **non-abelian lifted-product** CSS code LP(A,B) with 1×2 base matrices
over 𝔽₂[G], G = C₅×S₃ (GAP `SmallGroup(30,1)`, |G| = 30). Five data
blocks of 30 qubits, rate 1/5, check weight 9. The paper reports d = 10
exact from a Gurobi-IP workflow with no checkable artifacts; the lab has
independently SAT-certified d_X = d_Z = 10 in two solver families
(qec-lab PR #10, `certificates/mitten_150_30_10_{X,Z}.cert.json`).

Formalization target:

```
mitten150StabilizerCodeWithDistance : StabilizerCodeWithDistance 150 30 10
```

(`QEC/Stabilizer/Framework/Core/Logical/CodeDistance.lean`), sorry-free,
axiom profile ≤ {3 std + native_decide}. This would be the library's
largest bundled `[[n,k,d]]` object, its first qLDPC one, and the first
formally verified non-abelian qLDPC code distance anywhere.

**Gating constraint (user, 2026-08-10): the whole certification may add
at most 5 minutes to a full QECLean `lake build` (hard cap 10). See
`build_budget.md`.**

## Parameters

- **Physical qubits**: `n = 150` (5 blocks of |G| = 30)
- **Logical qubits**: `k = 30` (= |G|; rank H_X = rank H_Z = 60, both
  full — verified on the vendored matrices)
- **Distance**: `d = 10`, with d_X = d_Z = 10 individually
- **Family**: non-abelian lifted product ("mitten" shape: 1×2 base
  matrices, every entry weight 3)
- **Check weight**: 9 (all 120 generators); qubit degrees 3 or 6 per side

## Construction (paper Eq. (J1) / Definition 4)

Base data: G = C₅×S₃ and four weight-3 sets (GAP `Elements(G)` 0-based
indices, paper Table XIII):

```
a0 = (0, 14, 23)   a1 = (0, 2, 11)   b0 = (7, 20, 24)   b1 = (0, 2, 29)
```

With L(s)/R(s) the left/right regular-representation sums over a set s
(`L(a)[y,x] = 1 iff y = g·x, g ∈ a`; `R(b)[y,x] = 1 iff y = x·g⁻¹, g ∈ b`)
and `s* = {g⁻¹ : g ∈ s}`, the check matrices over the qubit blocks
D1..D5 are

```
H_X = [ L(a0)   0     L(a1)   0     R(b0*) ]     (X-checks, 2 blocks)
      [  0     L(a0)   0     L(a1)  R(b1*) ]
H_Z = [ R(b0)  R(b1)   0      0     L(a0*) ]     (Z-checks, 2 blocks)
      [  0      0     R(b0)  R(b1)  L(a1*) ]
```

CSS commutation H_X·H_Zᵀ = 0 is blockwise **L/R commutation** — left and
right multiplications commute by associativity. This (not abelian
convolution commutativity) is the chain-complex law for the Lean layer.

Grid indexing: D_(α,β) for α,β ∈ {0,1} with D1=(0,0), D2=(0,1),
D3=(1,0), D4=(1,1); X-check block β acts on column (·,β) by L(a_α),
Z-check block α acts on (α,·) by R(b_β); D5 is the shared block.

## Data provenance

`experiments/bb_lab/instances/mitten_150_30_10/{Hx,Hz,Lx,Lz}.npy` —
bit-identical to the authors' release (github.com/a7b/yarn @ `82fb695a`,
MIT), sha-pinned in the two bb-cert/v2 certificates. Validated facts
(script `validate` mode + this attempt's M0):

- shapes 60×150 / 60×150 / 30×150 / 30×150; H_X·H_Zᵀ = 0
- all check rows weight 9; ranks 60/60 ⟹ k = 30
- `Lx` rows weight **18**, ∈ ker H_Z; `Lz` rows weight **10**, ∈ ker H_X;
  `Lx·Lzᵀ = I₃₀` — row i of each is the same logical qubit
- shipped matrices = Table XIII rebuild up to relabeling (A26 §5); the
  Lean instance will be **generated from the Table XIII sets + group
  table** (canonical labeling), not from the shipped byte order
- group table: `instances/mitten_groups/group_30_1.txt` (GAP export,
  `Elements(G)` order; identity = index 0)

## Distance, precisely

Chain-level, both CSS directions:

- **d_Z-side floor**: every v ∈ ker H_X \ rowspace(H_Z) has |v| ≥ 10
- **d_X-side floor**: every v ∈ ker H_Z \ rowspace(H_X) has |v| ≥ 10
- **witness**: any `Lz` row is a weight-10 element of ker H_X, and its
  unit pairing against the corresponding `Lx` row certifies it is not in
  rowspace(H_Z) (over 𝔽₂: rowspace(H_Z)ᗮ = ker H_Z ∋ Lx row)

Pauli-level `HasCodeDistance C 10` follows from the two chain floors +
the witness through the generic bridge
(`Framework/Homological/Distance.lean` xChain/zChain decomposition), as
in the gross and cover300 assemblies.

**Equivalent floor form (the one to prove).** Since d = 10 is the true
distance, weight-≤9 kernel vectors are exactly weight-≤9 rowspace
elements, and M0 measurement (ISD, 6k information sets/side; exact SAT
census in this attempt) supports:

> the only nonzero weight-≤9 vectors of ker H_X are the **60 rows of
> H_Z themselves** (all weight exactly 9), and mirror on the other side.

So each floor is an A21-shaped statement: *light cycle ⟹ generator row*
(pair sums ≥ 14, sampled triples ≥ 15 — comfortable slack above 9).

## Structural facts for the proof (M0 measurements, 2026-08-10)

1. **Symmetry = center C₅ only.** Exhaustive search over all 30 group
   elements × all 2⁹ blockwise L/R translation assignments: exactly the
   5 central elements (uniform translation, L = R there) are code
   automorphisms. Aut(G) ≅ C₄×S₃ (order 24) preserves the four sets only
   at the identity. Orbit reduction for sweeps: **5×** (vs |G| = 72 for
   the gross code). Full Tanner-graph Aut census: M0 item (nauty).
2. **No X↔Z duality found** in the natural family (all grid
   permutations × antipode on/off per sector): the two floors must be
   proven independently. (Consistent with canonical weights 18 vs 10.)
3. **Elimination structure**: L(a1), R(b0), R(b1) all invertible over 𝔽₂
   (rank 30); L(a0)/R(a0) have rank 28. On the ker-H_X side the two
   check equations solve the grid blocks v₃, v₄ in terms of (v₁, v₂, v₅)
   — the kernel is freely parameterized by 90 coordinates; mirror on the
   other side. (The paper only *requires* L(a1), R(b1) full rank; b0
   invertible is a bonus of this instance.)
4. **ε-parity survives non-abelianity**: all eight entries have odd
   weight (3), so augmentation gives, on ker H_X: |v₁|+|v₃|+|v₅| ≡ 0,
   |v₂|+|v₄|+|v₅| ≡ 0 (mod 2), hence |v| ≡ |v₅| (mod 2) — a partition
   pruner for the split map.

## Lean-side shape (see plan.md)

- Group: `ZMod 5 × DihedralGroup 3` (mathlib's S₃ with decidable
  eq/fintype); a generated dictionary maps GAP indices → elements.
- New framework file `Framework/Homological/LiftedProduct.lean`:
  left/right convolutions on a finite (possibly non-abelian) group,
  `lpChainComplex` for the mitten shape, chain law by L/R commutation.
  Instantiates the *generic* `HomologicalCode` — everything downstream
  (CSS, logical correspondence, Pauli distance bridge) is reused as-is.
- Instance under `QEC/Stabilizer/Codes/Mitten/M150/`, generated data via
  a bb_lab emitter (falsify-first, per GENERATORS.md).

## References

- arXiv:2607.28795 (paper; Table I row 1, Table XIII, Eq. (J1), Thm 4)
- `experiments/bb_lab/notes/mitten150_tandem_verification.md` (SAT
  certification of d = 10, both directions, two solver families)
- `experiments/bb_lab/notes/A26_mitten_descent.md` (family structure,
  no deck for [[150]], Table XIII conventions, GAP tooling)
- Memory/program context: A21 base floor (`[[150,8,8]]`, same n, cap 7),
  A15 window engine + KernelCert pivot certificates (build-cost model)
