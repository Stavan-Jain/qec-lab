# A22 — analytic light-boundary completeness for f2a6f17e (moonshot)

**Target.** Replace the 9.6 h SAT enumeration-completeness verdict behind
`LightClassification` (QECLean `…/Z5Z15F2A6/Classification.lean`, PR #61)
with analytic structure: prove that the [[150,8,8]] code
`bb_neigh_z5z15_f2a6f17e` (A = 1 + y + x, B = x·y⁶ + x·y¹⁰ + x²·y¹² over
Z₅×Z₁₅) has EXACTLY 113 translation classes of nonzero boundaries of
weight ≤ 14.

**Status: in progress (session 1, 2026-07-22).**

## 0. The (ε,δ) fibering ansatz (session-1 working hypothesis)

Coordinates: z := y³ (order 5), w := y⁵ (order 3); G = Z₅(x)×Z₁₅(y)
≅ Z₁₅(base) × Z₅(fiber) with base Z₁₅ = ⟨x̄⟩×⟨w̄⟩ ≅ Z₅×Z₃ and fiber ⟨z⟩.
Then y = z²w², and

- A = 1 + x + z²w², B = xz² + xw² + x²z⁴,
- σ: (x,z,w) ↦ (xz², z⁻¹, w) — reflection about the diagonal of the
  (x,z)-torus, fixing the base pointwise (σ ≡ id mod ⟨z⟩).

CRT: R := F₂[z]/(z⁵−1) ≅ F₂ (ε: z→1) × GF(16) (δ: z→ζ, ζ a primitive
5th root; q(z) = 1+z+z²+z³+z⁴ irreducible since ord₂ mod 5 = 4). Hence

    F₂[G] = R[Z₁₅] ≅ F₂[Z₁₅] × GF(16)[Z₁₅],   u ↔ (u_ε, u_δ).

**Exact local weight formula.** (ε,δ): F₂[z]/(z⁵−1) → F₂×GF(16) is an
iso, so each fiber value (ε,δ) has a *unique* fiber polynomial, of known
weight W(ε,δ):

| δ \ ε | 0 | 1 |
|---|---|---|
| 0 | 0 | 5 (fiber = N = 1+z+z²+z³+z⁴) |
| ∈ μ₅ | 4 (fiber = N + monomial) | 1 (monomial) |
| ∉ μ₅∪{0} | 2 (two monomials) | 3 (three monomials) |

and |u| = Σ_{sites g ∈ Z₁₅} W(u_ε(g), u_δ(g)).

**The boundary pair in (ε,δ) coordinates.** f ↔ (f_ε, f_δ) ranges over
ALL of F₂[Z₁₅] × GF(16)[Z₁₅] (independent coordinates). Then
∂₂f = (u,v) with

- ε-parts: u_ε = Ā f_ε, v_ε = x̄·Ā f_ε where Ā = 1 + x̄ + w̄²;
  **Ā is invertible** in F₂[Z₁₅] (checked: no character zero), so
  substitute h := Ā f_ε — a FREE F₂[Z₁₅] coordinate; u_ε = h, v_ε = x̄h.
- δ-parts: u_δ = Ã f_δ, v_δ = B̃ f_δ over GF(16)[Z₁₅], where
  Ã = 1 + x + ζ²w², B̃ = ζ²x + w²x + ζ⁴x². All 15 characters η of Z₁₅
  are GF(16)-valued; expect a single common zero η₀ of (Ã,B̃)
  (= the kernel component: dim ker ∂₂ = 4 over F₂ = 1 GF(16)-line).

Writing α := u_δ (site-indexed), β′_g := v_δ(x̄g) (so that v's site x̄g
pairs with u's site g through the shared h_g):

    |u| + |v| = Σ_{g∈Z₁₅} [ W(h_g, α_g) + W(h_g, β′_g) ]   (EXACT)

with moduli (h, α): h free, α̂(η₀) = 0 (one GF(16)-linear condition),
β′ the explicit linear image of α (transfer T = B̃/Ã off η₀). Dimension
audit: 15 + 4·14 = 71 = dim im ∂₂ ✓. Translations act as Z₁₅ site-shift
(x,w) × μ₅ scalar on (α,β′) (z), order 75 ✓.

**Optimal-h cost taxonomy.** Site types by δ-value: O (δ=0), M (δ∈μ₅),
D (else). Per-site cost c(α_g,β′_g) = min_h [W(h,α)+W(h,β′)]:

| pair | c | argmin h | flip penalty |
|---|---|---|---|
| (O,O) | 0 | 0 | +10 |
| (O,M) | 4 | 0 | +2 |
| (O,D) | 2 | 0 | +4 |
| (M,M) | 2 | 1 | +6 |
| (M,D) | 4 | 1 | +2 |
| (D,D) | 4 | 0 | +2 |

Budget 14 ⟹ ≤ 7 active sites (each ≥ 2), and the δ-side constraint
"supp(α) ∪ supp(β′) ⊆ S, β′ = Tα, α̂(η₀)=0" is an *overdetermined*
GF(16) linear system for |S| ≤ 7 (|S| unknowns, 16−|S| conditions) —
nonzero solutions only at special S. So the classification should reduce
to a 2¹⁴-subset rank sweep (Gaussian pivot certificates, the
KernelCert.lean pattern) + finite h-flip accounting.

Hand-verified on 3 classes from `f2a6_light_classes.jsonl` before any
code: the |b|=10 near-kernel class is exactly (h = single site, α = 0)
(two full N-fibers at sites (4,1) → x̄·(4,1) = (0,1) ✓); the |b|=6 class
is 3 (M,M) sites at cost 2 each; a |b|=12 near-kernel class shows a
4/5-fiber = W(0,μ₅) = 4 site (h-flip case). Also explains d = 6 = 3·2
prospectively.

## 1. Phase 0 — the fibering VERIFIED end-to-end; SAT re-derived in ~2 min
## (2026-07-22)

`scripts/a22_eps_delta_structure.py`, ALL PASS (V1–V7):

- **V1–V2**: the (ε,δ) bijection, the weight table, the decomposition
  identities, and the EXACT site-weight formula all hold (50 random f +
  algebra asserts).
- **V3**: Ā = 1+x̄+w̄² invertible (15 nonzero character values).
- **V4**: Ã, B̃ have exactly ONE common character zero η₀ = (4,2), no
  other zeros; transfer T̂ orders {1:1, 3:3, 5:5, 15:5}.
- **V5**: τ (inverse DFT of T̂, extended by 0 at η₀) satisfies
  τ⋆Ã = B̃; |supp τ| = 10.
- **V6**: all 113 classes decompose per the ansatz (v_ε = x̄u_ε,
  δ-transfer, η₀-condition, weight formula — 113/113). Taxonomy
  (|b|, #active sites, optimal cost, h-excess, profile):
  - the |b|=6 class = 3 (M,M) sites at cost 2 each — and its S is a
    translate of supp(Ā) (the "A-triangle"); the class = ∂₂(monomial).
  - ALL cost-tight (excess 0) classes = 94 small classes exactly.
  - ALL 19 near-kernel classes = small δ-configs + ONE ε-flip:
    |b|=10 = pure-h (two N-fibers, α=0); 3× |b|=12 = the |b|=6 α with
    one (M,M) site flipped (+6); 6× |b|=14 = (D,D)³ configs flipped
    (+2); 9× |b|=14 = (D,D)(M,M)⁴ flipped (+2). The "mysterious"
    near-kernel stratum is now fully demystified: it is the ε-sector
    (N-fiber insertions) over small δ-configs, and the 31–33 coset
    minima are the cost of realizing an ε-flip through Ā⁻¹.
- **V7 (headline)**: independent complete re-derivation. Rigorous
  active-site bound ≤ 7 (each active site costs ≥ 2 for ANY h) ⟹ sweep
  all C(15,≤7) = 16,384 site-subsets, solve the GF(16) linear system
  (supp α ⊆ S, supp Θα ⊆ S, α̂(η₀) = 0): 9,235 special subsets,
  1,377,525 nonzero α, 7,050 with optimal cost ≤ 14, **94 translation
  classes of α** (+ 1 pure-h) → h-flip enumeration → 113 canonical
  (u,v) classes **exactly equal** to the SAT list. Wall time ≈ 2 min
  vs 9.6 h SAT; and every ingredient is a finite algebraic fact
  (CRT iso, 32-entry table, arithmetic bound, rank certificates,
  flip accounting) — i.e. kernel-checkable in principle, no solver.

**Consequence for `LightClassification`:** the 2⁷⁵-quantified Prop
reduces to: (a) the fibering lemmas (algebra), (b) the weight table
(finite), (c) "active ≤ 7" (arithmetic), (d) a 16,384-subset GF(16)
rank sweep (pivot certificates, the exact KernelCert.lean pattern
Stage 4 shipped for the window sweeps), (e) finite flip/translation
accounting. This is the certificate-compression target of the moonshot
achieved at the numeric level; Lean feasibility now looks *plausible*
rather than aspirational.

## 2. Phase 1 — delta-side structure (2026-07-22)

`scripts/a22_delta_structure.py`, `scripts/a22_emit_dataset.py`;
dataset `data/a22/alpha_classes_full.json` (the canonical per-class
record: α, active sites, per-site types/costs, m*, flip-children,
file-class match — all 113 matched).

**P1 (m-smallness).** Every light α = Ãm for an m with ≤ 3 sites
(min over the 16 η₀-line preimages): histogram {1: 3, 2: 41, 3: 50}.
So the small stratum is literally "boundaries of site-support-≤3
δ-data" — the SAT-observed coset minima {1,2,3} explained.

**P2 (geometry).** 94/94 active sets ⊆ supp(m*) + Ā-triangle. m*
values are mostly μ₅ (profiles MMM:50, MM:35, DM:6, D:2, M:1).

**P3 (symmetry).** No Galois site-symmetry exists (correct in
hindsight: the F₂-classes are Galois-stable; α-coordinates are merely
covariant). The Φ-involution acts on δ-data as
**α ↦ ζ^{2i(g)}·Frob²(β′)** (fiber conjugation + x-twist + block
swap; 75 translated forms all verified to permute the 94 reps).
Orbits: 94 → 55 families (16 fixed + 39 pairs).

**P4 (|S| structure).** Light α classes: |S| ∈ {3: 3, 5: 21, 6: 32,
7: 38} — none at 1, 2 (⟹ d = 6 analytic, see below) and none at 4
(cost phenomenon, certified by the sweep).

**P5 (sweep bookkeeping).** Special subsets (nonzero kernel): 9,235
of 16,384, by size {3: 15, 4: 180, 5: 945, 6: 2830, 7: 5265}; kernel
dims {1: 6680, 2: 2250, 3: 290, 4: 15}. The 15 special |S|=3 subsets
are EXACTLY the translates of the Ā-triangle; its kernel line's three
μ₃-cosets give the |b|=6 class (c ∈ μ₅) and the two (D,D)³ |b|=12
classes (c in the ω/ω²-cosets) — the minimum-weight geometry is
completely forced.

**Line-profile collapse (verified).** Site types are governed by the
μ₃-character δ ↦ δ⁵: along a kernel line c·v the whole type vector
shifts uniformly by c⁵. Hence each kernel line needs only 3 cost
evaluations (μ₅-scalars are z-translations). Total cost evaluations
across all special kernels: **568,905**.

## 3. The classification tree (complete, all counts close)

| stratum | #α | children (|b|: count each) | classes |
|---|---|---|---|
| pure-h (α=0, |h|=1) | — | 10: 1 | 1 |
| |S|=3, cost 6, |m*|=1 (triangle, c∈μ₅) | 1 | 6: 1, 12: 3 | 4 |
| |S|=3, cost 12, |m*|=1 (triangle, ω-cosets) | 2 | 12: 1, 14: 3 | 8 |
| |S|=5, cost 10, |m*|=2 | 6 | 10: 1 | 6 |
| |S|=5, cost 12, |m*|=2 | 9 | 12: 1, 14: 1 | 18 |
| |S|=5, cost 14, |m*|=2 | 6 | 14: 1 | 6 |
| |S|=6, cost 12, |m*|∈{2,3} | 22 | 12: 1 | 22 |
| |S|=6, cost 14, |m*|=3 | 10 | 14: 1 | 10 |
| |S|=7, cost 14, |m*|=3 | 38 | 14: 1 | 38 |
| **total** | 94 | | **113** |

(|b| totals: 6 → 1, 10 → 7, 12 → 36, 14 → 69 ✓ = the file. The 19
near-kernel classes = exactly the flip-children with excess > 0.)

**Analytic bonus: d(im ∂₂) = 6 with forced geometry.** |S| ≤ 2 has no
special subsets (120 tiny full-rank certificates), so nonzero light
means pure-h (weight 10) or |S| ≥ 3 (cost ≥ 6); the |S| = 3 minimum
is attained only on Ā-triangle translates, kernel = the Ãm line
(m = δ-monomial), i.e. the minimum-weight codewords are exactly the
∂₂(fiber-monomial) family with ε-optimal dressing.

## 4. Lean proof architecture for `LightClassification` (assessment)

The 2⁷⁵ Prop reduces to, in order:

- **L1 (algebra):** the fibering iso F₂[G] ≅ F₂[Z₁₅] × GF16[Z₁₅]
  and the three decomposition identities (u_ε = Āf_ε, v_ε = x̄u_ε,
  δ-transfer via τ⋆Ã = B̃ — the last is a single finite convolution
  identity). GF16 can be avoided entirely: δ-parts are F₂⁴-vectors
  (reduction mod the quartic q(z)), all maps F₂-linear; the "GF16
  linear systems" become F₂ matrices of 4× the dimensions.
- **L2 (finite):** the 32-entry weight table; the per-site cost/flip
  table (6 cases).
- **L3 (arithmetic):** weight ≤ 14 ⟹ ≤ 7 active sites (each ≥ 2).
- **L4 (rank certificates):** the 16,384 subset systems — Gaussian
  pivot certificates exactly as in the shipped `KernelCert.lean`
  pattern (each ≤ 64×28 over F₂ after the 4× blowup); optionally 15×
  reduced by one translation-equivariance lemma to 1,096 orbits.
- **L5 (kernel dispositions):** per special subset, kernel basis cert
  + 568,905 line-coset cost evaluations (each: 15 site-costs of an
  explicit vector pair, sum, compare 14) — packed-table
  decide/native_decide chunks; same order as the Stage-4/5 sweep
  engineering that already shipped.
- **L6 (glue):** the 94+1 survivors → flip-children → identification
  with the 113 tabulated `repChain` reps (emitted decide tables), and
  the translation-orbit walk (`weight_floor_translate1_reduce`
  machinery already exists on the branch).

Verdict: **full kernel-checked LightClassification is
engineering-feasible** — no SAT, no 2⁷⁵ enumeration, total obligation
mass ≈ 6·10⁵ tiny checks + 1.6·10⁴ pivot certs, comparable to the
(M)-route sweep obligations already built in Stage 4/5. This answers
the moonshot's reachability question positively at the architecture
level; the remaining work is an emitter + a Lean statement layer
matching `Defs.lean` conventions (session 2).

## 5. Session log

- 2026-07-22: fibering ansatz drafted; Phase-0 script ALL PASS on
  first complete run (one loader fix). Phase 1: m-smallness, Φ-action
  on α, triangle geometry, line-profile collapse, dataset emitted
  (`data/a22/alpha_classes_full.json`), obligation counts computed.
  Next session: Lean statement layer + emitter prototype (own
  worktree, NOT a15-m-kernel-route), starting from L1/L2.
