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

## 5. Session 2 (2026-07-22): `lightClassification` PROVEN in Lean —
## sorry-free, axiom-clean, ~14 min build

**Headline: `LightClassification` is now a THEOREM** —
`lightClassification : LightClassification` compiles on branch
`claude/a22-light-classification` (worktree off `f9159e9`), zero
sorries, `#print axioms` = {propext, Classical.choice, Quot.sound} + 9
named `native_decide` obligations. The 9.6 h SAT-enumeration
certificate hypothesis of PR #61 is replaced by an analytic,
kernel-checkable proof.

### The architecture as SHIPPED (diverges from §4's plan in 3 ways)

1. **No GF(16), no CRT iso, no H-matrix anywhere.** δ-extraction is the
   per-fiber F₂-map `dₙ = tₙ + t₄` (reduction mod `q(z)`), ε is the bit
   sum, and the inverse is the *ring identity* `σbit` (`t₄ = ε + Σd`,
   `tₙ = dₙ + t₄`) — the entire (ε,δ)↔fiber bijection is one 160-case
   `decide` (`σbit_eps_delta`). The per-site weight/cost tables emerge
   from σ + popcount (32-entry `W5TAB`, decide-checked); μ₅/M/D
   taxonomy never appears.
2. **Two reductions cut 16,384 subsets → 429 certificates**: (a)
   supp-monotonicity ⟹ only *maximal* (size-7) subsets need
   certificates; (b) Z₁₅ translation-orbit reduction (429 = 6435/15,
   free action), realized in-proof by `gOfTau τ = (τ_i, 5·τ_b)` which
   shifts sites with **zero fiber rotation** (`2·5τ_b ≡ 0 mod 5`), so
   δ-data transports as a pure site permutation (`dU_translate`).
   `SHIFT_ANS` (16,384-entry table) gives each ≤7-site mask its (τ, j).
3. **Orthogonality via a transposed syndrome fold, not per-W basis
   checks**: pivot functionals W (120-bit, RREF row-combinations of the
   left-null-space of the δ-data map) are certified by
   `xorFoldCols W = 0` — XOR of `COLPACK` columns (the *transposed*
   basis-image matrix) at W's set bits — one 75-bit accumulator instead
   of 75 parities per W. One semantic `native_decide` bridges `COLPACK`
   (and `DGEN`, its transpose, used for the generator-preimage
   identities `rowFold (pre) = gen`) to `deltaData ∘ ∂₂` on the
   75-point basis; `funLiftF2` lifts to all boundaries.

### Proof chain (9 new Lean files + 2 generated)

`LightSite` (site geometry: `siteEquiv` regroup, `recon_eq` σ-bijection
reconstruction, `active_card_le` ≤7-sites, `flip_card_le` ≤1-outside-
flip — the two card bounds are 1024-case fiber-pair decides + a generic
`card·m ≤ Σ` count) → `LightCertData`/`LightTTData` (generated; 1,033
KB + 422 KB, **1,024-entry chunked arrays** — flat `#[...]` literals at
16K+ entries blow `maxRecDepth`/`isDefEq`; `NAME_CHUNKS : Array (Array
Nat)` dispatchers fix it) → `LightLinear` (maskFun unpackers, `xorSelTab`
generic fold + testBit-sum lemma, the 3 semantic natives, `packChain`
round trip) → `LightPeel` (`peel_delta` + `classify`: y = span element
of its own free-coordinate code `e < 2^{4k}`, via peel of `y + Σ eᵢ·genᵢ`
whose ∂₂-preimage correction `f + Σ eᵢ·preᵢ` keeps everything a
boundary's δ-data) → `LightChecks` (`subCert_all` over all 429 reps,
`shift_all`, `tt_decode`, fuel-bounded `memberTT` binary search with
one-sided soundness) → `LightSweep` (`checkA_all`: mincost prefilter
completeness vs the emitted survivor lists; `checkB_all`: all
(survivor, 2⁷ in-rep ε-patterns, ≤1 outside flip) light reconstructions
∈ translate table; `expandMask`/`extractS` round trip) → `LightBridge`
(`weight_eq_wtOf` semantic-weight = packed-table-weight, `wtOf_or_single`
+10 flip split, `nib_arith` mod-16 bit decomposition, translation
transport) → `LightAssembly` (`rep_classify` per-rep core +
`lightClassification`).

### Generator (`a22_gen_light_classification.py`)

Pure-F₂ pipeline in the exact Lean conventions (no GF(16)); emits
COLPACK/DGEN/REP7/SHIFT_ANS/certificates (21,996 pivots, 2,028
generators with ∂₂-preimages)/13,780 survivors/TT (8,475 sorted
translate rows). Hard-asserts: rank(Φ)=56, dim V'⊥=64, per-rep RREF
triangularity + kernel dims (histogram {0: 78, 4: 215, 8: 117, 12: 18,
16: 1}), σ/weight-formula on random chains, **sweep soundness** (2,200
distinct recon masks all ∈ TT) and **route completeness** (every one of
the 8,415 distinct TT masks re-derived through the exact
SHIFT_ANS→translate→sweep path the Lean proof takes).

### Build costs (wall, M-series laptop, shared with A23)

CertData 110 s (parse), TTData ~45 s, Site 26 s, Linear 43 s (3
natives), Peel 20 s, Checks 58 s, **Sweep 498 s** (checkA+checkB
natives — the dominant cost, ~44.5K reconstruction packs at 150
function-evals each), Bridge 58 s, Assembly 9 s ≈ **14.5 min total**
for the A22 layer.

### Traps hit (for lean-patterns.md)

- Array literals: >2K entries risk `maxRecDepth`; 16K entries hit
  `isDefEq` timeouts. Chunk at 1,024 + `Array (Array Nat)` dispatcher.
- `rcases h : boolExpr with _ | _` *substitutes* the scrutinee in the
  goal — a following `rw [h]` fails with "pattern not found".
- `Finset.filter`-instance whnf pits: mixing `set`-abstracted and raw
  forms of the same filter in `simp only [...] at`-style Bool
  extraction hit 2M-heartbeat `whnf` timeouts; explicit `by_cases` +
  targeted `rw` of individual Bool atoms (plus `rw [← hMdef]` to
  re-align `set`-atoms before `omega`) avoids it.
- `ite` at function type: ascribe the else-branch
  (`else (0 : G150 → ZMod 2)`) or `Pi.single` loses its family.
- Renames on this toolchain: `Nat.testBit_eq_decide_div_mod_eq`,
  `Finset.notMem_empty`, `Finset.card_insert_of_notMem`, `push Not`
  (not `push_neg`).

### Consequence

The `_of_classification` theorems in `Distance.lean` can now consume
`lightClassification` directly (NOT wired in this session per scope);
the [[300,8,16]] hypothesis set drops to {SeamCosetFloor 16} once
integrated. The moonshot's Lean-feasibility question is closed
**positively and completely** for LightClassification.

## 6. Session log

- 2026-07-22: fibering ansatz drafted; Phase-0 script ALL PASS on
  first complete run (one loader fix). Phase 1: m-smallness, Φ-action
  on α, triangle geometry, line-profile collapse, dataset emitted
  (`data/a22/alpha_classes_full.json`), obligation counts computed.
  Next session: Lean statement layer + emitter prototype (own
  worktree, NOT a15-m-kernel-route), starting from L1/L2.
- 2026-07-22 (session 2): full Lean chain designed, generated, and
  proven in one session — see §5. Branch
  `claude/a22-light-classification`; generator
  `scripts/a22_gen_light_classification.py`.
