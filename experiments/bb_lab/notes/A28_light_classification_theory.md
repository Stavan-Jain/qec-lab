# A28 — Theory for the light-stabilizer classification (LSC) problem: a certified BZ census lane, the shift-bound verdict, and the ε-trisection theorem

> **STATUS (2026-08-06, session 1).** The user-directed probe of the two
> classical residues from the 2026-08-06 teaching-session Q&A: van
> Lint–Wilson shifting as a (CLASS)-systematizer, and Brouwer–Zimmermann as
> an independent census certifier — plus new theory where the classical
> material runs out. Headlines: **(1)** a BZ-style certified census lane
> (solver-free, ~100-line C kernel + linear algebra) that re-derives the
> f2a6 113-class census in **4.5 s** (SAT: 9.6 h + 34,554 s final UNSAT)
> and settles the **A17-docket [[180,4,10]] censuses the SAT lane could not
> touch** (docket37: 2,203 classes; docket5e: 2,371 classes; W = 18,
> ≈ 11–15 min each, exact-node certificates); **(2)** the **ε-trisection theorem** — for G = Z₂ × H the
> LSC census provably decomposes into a doubled H-census, finitely many
> coset censuses, and bounded lifts — validated clause-by-clause on all
> 2,203 docket37 classes; **(3)** a **measured negative** on shifting as
> the kill engine: beam = exhaustive on test oracles, yet a
> BZ-certified-empty cell exists where the best possible game value is 8
> against the needed 15 — the gap is intrinsic to the technique, not to
> search. Failures are first-class; §7 has the verdict table.

Companions: `A22_analytic_classification.md` (the fibering census for
f2a6), `A27_safe_floor_generality.md` (the docket UNKNOWNs this feeds; its
§2 cost anatomy is assumed), `A4_writeup.md` §6.3 (the gross-base (CLASS)
hand proof, re-derived here in 0.2 s), `A17_d7plus_doubling_hunt_plan.md`
(docket provenance).

---

## 1. The LSC problem, in one place

For an abelian group G = Z_ℓ(x) × Z_m(y) and A, B ∈ F₂[G], the X-side
stabilizer space of the BB code is the **boundary code**

    C(A,B) = { ∂f := (f⋆A, f⋆B) : f ∈ F₂[G] } ⊆ F₂[G]² ,

a 1-generator quasi-abelian code of index 2, dim = |G| − k/2. The **LSC
problem** at threshold W: list all translation classes of nonzero ∂f with
|∂f| ≤ W. It is I2 of A27 §2 (min-distance-type, NP-hard), the (CLASS)
input of the doubling template, and the 2^{75}-quantified hypothesis that
A22 spent a moonshot deleting for one instance.

| instance | G | A; B | (n,k,d) | W | prior census |
|---|---|---|---|---|---|
| grossbase | Z₆×Z₆ | x³+y+y²; y³+x+x² | [[72,12,6]] | 11 | hand proof (A4 §6.3): hexagons + D-pairs |
| f2a6 | Z₅(x)×Z₁₅(y) | 1+y+x; xy⁶+xy¹⁰+x²y¹² | [[150,8,8]] | 14 | SAT 9.6 h: 113 classes; A22 fibering 2 min |
| docket37 | Z₁₅(x)×Z₆(y) | 1+y+x; y⁴+x+x¹¹y² | [[180,4,10]] | 18 | **OPEN** (SAT extrapolated days-plus) |
| docket5e | Z₁₅(x)×Z₆(y) | 1+y+x; y⁴+x⁸y²+x¹³ | [[180,4,10]] | 18 | **OPEN** |

Conventions locked against the a17 census file (`a28_v_foundations.py`,
ALL PASS): b = (f⋆A, f⋆B), block 0 = A-block, supports [block, x, y].

## 2. Lane 1 — the Brouwer–Zimmermann census certifier (shipped, validated)

**Method** (`a28_bz_census.py` + emitted C kernel). Take two DISJOINT
information sets I₁, I₂ for C(A,B) (full-rank column sets, found by
interleaved/randomized RREF orders — a greedy first set can starve the
complement: on docket37 taking all pivots in the A-block leaves a rank-80
complement because dim Ann(B) = 10), with systematic bases G₁, G₂. Every
codeword is reconstructed from its restriction to either set; by
disjointness, |c| ≤ W forces min(|c|_{I₁}, |c|_{I₂}) ≤ ⌊W/2⌋ =: r. So DFS
over all row-combinations of size ≤ r in both systematic matrices visits
every census word. The **completeness certificate** is finite and
auditable: pivot sets + disjointness + systematicity + the exact node
count Σ_{s≤r} C(κ,s) per set (asserted; any mismatch = enumeration hole).
No SAT, no DRAT, no fibering; the trusted base is a 100-line enumerator
plus F₂ linear algebra.

**Results.**

| run | κ | r | nodes/set | wall | classes (by weight) | cross-check |
|---|---|---|---|---|---|---|
| grossbase | 30 | 5 | 174,436 | 0.2 s | 7 = {6:1, 10:6} | **= (CLASS)** (A4 §6.3) — third independent derivation |
| f2a6 | 71 | 7 | 1,487,160,959 | 4.5 s | 113 = {6:1, 10:7, 12:36, 14:69} | **exactly** the SAT census; raw distinct words 8,415 = A22's TT mask count |
| docket37 | 88 | 9 | 642,559,165,633 | 661 s | **2,203** = {6:1, 10:6, 12:42, 14:54, 16:478, 18:1622} | NEW; trisection-validated (§4) |
| docket5e | 88 | 9 | 642,559,165,633 | 870 s | **2,371** = {6:1, 10:6, 12:43, 14:55, 16:503, 18:1763} | NEW; trisection-validated (§4) |

The f2a6 line is the calibration: 4.5 s vs 9.6 h + 34,554 s = a ~10⁴×
wall reduction over the SAT lane at identical output, with a certificate
a referee can re-run. The docket lines are **new mathematical data**: the
(180, 18) censuses A27 §2 wrote off as days-plus, each with exactly one
weight-6 class (the stamp) and six weight-10 classes (D-pair analogues) —
the gross-base small-pattern geometry recurs at d = 10.

Solver-cost anatomy note: this does NOT contradict the A27 inversion
story (absence-proving UNSAT is the SAT lane's wall); BZ sidesteps it
because completeness is a counting invariant, not a refutation.

## 3. Lane 2 — shifting: the theory, and what it is actually worth here

### 3.1 Theorem S (abelian shift bound; van Lint–Wilson transported)

G odd, Ĝ = characters into F̄₂ (values in GF(2^{ord₂(exp G)}); GF(16)
suffices for every instance here). For 0 ≠ c ∈ F̄₂[G] and χ ∈ Ĝ let
v_χ := χ|_{supp c} ∈ F̄₂^{|c|}. Define **independent sets w.r.t. an
oracle O ⊆ Z(c)** by the game: ∅; if I ⊆ O and b ∉ Z(c) then I ∪ {b};
any I ↦ μI (μ ∈ Ĝ). Then every reachable I has rank{v_χ : χ ∈ I} = |I|,
hence **|c| ≥ |I|**.

*Proof.* Close: I ⊆ O ⊆ Z(c) puts span(v_I) ⊆ c^⊥ while ⟨v_b, c⟩ =
ĉ(b) ≠ 0, so the rank grows by 1. Shift: v_{μχ} = v_μ ∘ v_χ entrywise and
v_μ has unit entries, so Hadamard-multiplication by it is invertible. ∎

For c = fA the zero set is exactly Z(f) ∪ Z(A) (semisimple pointwise
product), so the game runs on hypothesized spectral patterns of f.

### 3.2 Lemma B (block additivity)

For b = (u, v), u-game and v-game histories I_u, I_v (each on its own
block, closers valid for its own word) have jointly independent vectors:
rank = |I_u| + |I_v|, so |b| ≥ ℓ_u + ℓ_v. *Proof.* A pure-u vector
vanishing on the v-block cannot use v-block vectors in a dependency:
project any dependency to the u-block. ∎

### 3.3 Theorem J (the syzygy game) and its measured value

Enrich states to pairs (χ, λ), λ ∈ P¹(F̄₂), vectors V_{χ,λ} =
(λ₁χ|_{supp u}; λ₂χ|_{supp v}); pairing ⟨V_{χ,λ}, b⟩ = f̂(χ)·(λ₁Â +
λ₂B̂)(χ). Three orthogonality sources: χ ∈ Z_A ∩ Z_B (kernel, always);
χ ∈ Z(f) (hypothesis); and λ = **ann(χ) := [B̂(χ) : Â(χ)]** — the
syzygy line, free of any hypothesis on f. Its slope function
ρ = B̂/Â : Ĝ ∖ kernel → P¹ is the abstract home of the C-ratio
bookkeeping in gross's hand proof (A4 §§6.2–6.3): direction forcing and
co-point rigidity are statements about level sets of ρ on radical layers.
Shifts move (χ,λ) ↦ (μχ, λ) (Hadamard, unit); per-χ at most 2 lines
(λ ↦ V_{χ,λ} is linear); close under whole-state orthogonality. Rank =
|I| as before; |b| ≥ |I|. (`a28_joint_game.py`, replay verifier included.)

**Measured** (f2a6, 9 sample cells): the naive whole-state rule makes the
joint game *weaker* than Lemma-B additivity everywhere tested (Δ = −2 to
−8): mixed vectors obstruct the independent per-block walks. A refined
rule (per-block orthogonality with two-sided rank tracking) that would
subsume both is OPEN — not attempted beyond design.

### 3.4 The dichotomy tree (the systematizer shape)

Galois-orbit-level branching is exact (f̂(χ) ≠ 0 iff nonzero on the whole
2-cyclotomic orbit), so cells (P_in, P_out) — f̂ ≡ 0 on P_in, ≠ 0 on
P_out — partition the f-space. At each cell: certified bounds (Lemma B +
Theorem J with free closers P_out); bound > W kills the cell with a
replayable move-list certificate; otherwise branch on an undecided orbit.
Survivor cells = where census classes live; depth-capped cells = honest
giveups. Implemented in `a28_dichotomy.py`.

### 3.5 The verdict measurements

f2a6 numbers (`a28_m_f2a6.py`, `a28_diag_gap.py`, `a28_diag_joint.py`):

- **Spectral terrain**: |Z_A| = |Z_B| = |Z_A ∩ Z_B| = 4 = k/2 — the
  polynomials have NO zeros beyond the forced kernel; 20 Galois orbits
  (1+2+18×4); the 113 classes occupy only **39 distinct visible spectral
  patterns**, 22 classes at the empty pattern.
- **Soundness**: game bound ≤ true weight on all 226 census block words;
  0 violations anywhere; all kill certificates replay.
- **Root exactness**: ℓ(Z_A) + ℓ(Z_B) = 3 + 3 = 6 = the true census
  minimum — a ~6-move analytic certificate of d(im ∂₂) = 6 (matches A22
  §3's "analytic bonus" by a different route).
- **Search is not the bottleneck**: exhaustive BFS = beam on all test
  oracles (3/4/5/5).
- **The gap is intrinsic**: the cell "f̂ ≡ 0 on orbits {0,1,5} (+kernel)"
  is **BZ-certified empty** below 15 (ideal dim 62, both info sets,
  560,339,291 nodes/set, 0 hits) — yet the best game value there is
  4 + 4 = 8 against the needed 15. No search improvement can close this;
  the character-evaluation rank family is simply too weak in this regime.
- **Why the regime is hostile** (the structural insight): C(A,B) has
  dim 71 of 75 — the defining sets classical spectral bounds grip on have
  size 4, scattered, run-free. vL–W-class machinery excels in the
  opposite corner (low-rate, structured zero sets).
- **Kills exist but deep**: one 4-orbit-in cell reaches union bound
  8 + 8 = 16 ≥ 15; the depth-6 tree (125 nodes, 21 s) achieves 4 kills /
  59 giveups with all 113 classes safely in giveups; the depth-10 tree
  kills more but the partition cost grows exponentially while
  certified-empty-but-game-blind cells (above) can never be closed.

**Verdict on the user's scouted residue**: van Lint–Wilson shifting is
sound, cheap, exactly computable here, and delivers (i) tiny root-level
certificates (d = 6), (ii) the conceptual unification of gross's C-ratio
endgames via ρ, (iii) a partition scaffold with replayable kill
certificates. It is **not** a (CLASS) kill engine at threshold 2d − 2:
measured ceilings sit 3–7 below the needed values, with a certified
blind-spot cell. "Moderate odds" resolved: the systematizer role survives
only as scaffold + endgame-language; the census workhorse is Lane 1 (BZ)
or A22-style spatial fibering.

## 4. The ε-trisection theorem (new; the general Z₂-part chapter)

**Theorem T.** Let G = ⟨s⟩ × H, s² = 1, |H| odd; π : F₂[G] → F₂[H] the
s ↦ 1 quotient; sheets f = f₀ + s f₁, f̄ := π(f) = f₀ + f₁; ε = 1 + s;
K_H = ker ∂_H for the quotient pair (Ā, B̄) = (π A, π B), and sheet
polynomials A = A₀ + sA₁, B = B₀ + sB₁. Every nonzero b = ∂_G f with
|b| ≤ W falls in exactly one sector:

- **I (pure-ε).** f̄ ∈ π(K_G)-orbit of 0: a representative has f₀ = f₁,
  and b = ε ⊗ ∂_H f₁ with |b| = 2|∂_H f₁|. Classes biject with
  census_H(⌊W/2⌋) (s acts trivially on ε-words, so G-translation classes
  = H-classes).
- **II (kernel-coset).** f̄ ∈ K_H ∖ (π(K_G)·0): both quotient blocks
  vanish, b = ε ⊗ (u₀, v₀) with (u₀, v₀) ∈ (f̄A₀, f̄B₀) + ∂_H F₂[H] — an
  AFFINE-coset census at threshold W/2, one per f̄-class (finitely many:
  2^{dim K_H} − … quotient-kernel classes).
- **III (lift).** b̄ := ∂_H f̄ ≠ 0: then |b| ≥ |b̄| (sheet inequality
  |u| = |u₀| + |u₁| ≥ |u₀ + u₁|), b̄ is a census_H(W) member, and
  |b| = |b̄| + 2·excess with excess = |u₀ off supp ū| + |v₀ off supp v̄|
  minimized/enumerated over the same affine coset — the m(b)-species of
  A3 Entry 5.

*Proof sketch.* u = fA has sheets u₀ = f₀A₀ + f₁A₁, u₁ = f₀A₁ + f₁A₀;
u₀ + u₁ = f̄Ā; |u| = |u₀| + |u₁|; the boolean identity |x| + |x + b̄| =
|b̄| + 2|x off supp b̄| gives III; in I/II, u₀ = f̄A₀ + f₁Ā parametrizes
the affine coset as f₁ runs over F₂[H]; sector membership is
translation-invariant and (I vs II) well-defined modulo π(K_G). ∎

So **census_G(W) = doubled census_H(W/2) ⊔ coset censuses ⊔ bounded
lifts of census_H(W)** — the exponential core halves, and the pieces are
exactly the species the program already computes (censuses, coset
minima, off-support flip accounting = A22's ε-flips). For 2-parts Z₄ or
Z₂² the theorem iterates along a filtration (gross's 4-sheet layer frame
is the iterated case; not re-derived here).

**Validation on certified data** (`a28_trisection.py`, both docket codes,
all asserts pass — **4,574 classes total, zero violations**):

- docket37: π(K_G) = 0 (sectors sharp); sheet identities and ū = f̄Ā,
  v̄ = f̄B̄ hold for all 2,203 classes; sectors (I, II, III) =
  (2, 0, 2201); sector I = doubled census_H(9) = {6:1, 8:1}
  **bijectively** (independent H-side BZ run; the weight-12 and weight-16
  pure-ε classes); all 2,201 sector-III b̄'s land inside census_H(18)
  (4,883 classes; only 749 are hit — the lift constraint prunes 85% of
  the H-census); excess even always. H-instance: Z₁₅×Z₃, Ā = 1 + t + x,
  B̄ = t + x + x¹¹t².
- docket5e: sectors (1, 0, 2370); its H-quotient (B̄ = t + x⁸t² + x¹³)
  has census_H(9) = {6:1} only, matching the single pure-ε class
  (weight 12) bijectively; 2,370 sector-III classes over 1,104 of
  census_H(18)'s 4,417 H-classes. Same π(K_G) = 0, same identities.
  Note the two docket codes have DIFFERENT y-axis quotients here, while
  A27 §3's x-axis fibering gave them a shared ε-quotient — the two
  decompositions are complementary frames, not the same fact.

Relation to A27 §3: complementary decompositions — A27's probe fibers
the x-axis 5-part (sites Z₃×Z₆, shared ε-quotient [[36,4,2]]); Theorem T
quotients the y-axis 2-part (H = Z₁₅×Z₃, per-code quotients differ).
Both are instances of "frame first, floors second" (A27 §2); Theorem T is
the one with a clean general statement and a validated census reduction.

## 5. Lemma F (Frobenius) — a definitive small negative

For odd G, f² = σ₂(f) where σ₂ : g ↦ g² is a group automorphism — the
"Frobenius symmetry" of the census IS the monomial symmetry σ₂ and adds
nothing beyond the standard orbit accounting. On f2a6 the induced map
(with τ = σ₂⁻¹ = (a,b) ↦ (3a, 8b)) is the identity: 113 fixed classes,
zero compression (`a28_frobenius.py`). Matches A22's P3 ("no Galois
site-symmetry"). Closed.

## 6. What this changes for the program

1. **The docket census wall is deleted.** The [[180,4,10]] UNKNOWN-cell
   censuses are BZ-reachable in ~11 min each: docket37 shipped
   (`data/a28/census_docket37.json` + `bzcert_docket37.json`), docket5e
   in flight (§9 records its landing). A27 §5's gate 1 (V6-analog
   sampling) can run against real lists; the ε-chapter (gate 2) has its
   general theorem (§4) with sector II EMPTY and sector I trivial for
   docket37 — the classification tree is all sector III over 749
   H-classes.
2. **Any future (CLASS) target gets a same-day census.** The BZ lane is
   generic over the registry; cost ~C(κ, ⌊W/2⌋). Practical ceiling ≈
   10¹³ nodes/set on a laptop (κ ≈ 100 at W = 20); beyond that,
   Theorem T's halving is the designed escape.
3. **Certification story.** The BZ certificate (pivots + node counts) is
   a deterministic counting invariant — smaller trusted base than
   SAT+DRAT at these sizes, and the natural Lean shape is a verified
   enumerator or per-info-set weight-distribution recomputation
   (KernelCert-style pivot certificates cover the linear-algebra half
   already). Not attempted this session.
4. **Shifting takes the endgame-language role**, not the kill role:
   ρ-function level sets as the uniform vocabulary for A4-style survivor
   endgames; root-level d-certificates for free.

## 7. Verdict table (falsify-first ledger)

| claim probed | verdict | evidence |
|---|---|---|
| BZ as independent census certifier | **CONFIRMED, shipped** | §2 table; f2a6 = SAT exactly; grossbase = (CLASS) exactly |
| docket censuses reachable | **CONFIRMED** | 2,203 / 2,371 classes, 11–15 min each, exact node counts |
| vL–W shifting as (CLASS) kill engine | **REFUTED (measured)** | certified-empty cell with game ceiling 8 < 15; beam = exhaustive |
| shifting as scaffold/endgame language | confirmed-in-part | root d = 6 tight; 4-orbit kills; ρ ↔ C-ratios |
| syzygy game (naive rule) beats union | **REFUTED** | Δ ∈ [−8, −2] on all sampled cells |
| ε-trisection theorem | **VALIDATED** (2,203/2,203 clauses) | §4; sector I bijection via independent H-BZ |
| Frobenius census compression (odd G) | **REFUTED (trivial)** | §5 |

## 8. Open problems / next gates

1. **Refined mixed-rule joint game** (two-sided rank tracking) — would
   subsume Lemma B + Theorem J; only then re-test the kill ceiling.
2. **Survivor endgames**: per-pattern classification inside the 39 f2a6
   patterns via ρ-level-set arguments — the honest remaining distance
   between "census as data" and "census as theorem" on the analytic lane.
3. **Trisection converse execution**: enumerate sector-III lifts per
   H-class (the m(b)-species coset problems) for docket37 → a fully
   H-side derivation of the 2,203 list; then Lean-shape it (A22 L1–L6
   pattern).
4. **Iterate Theorem T** for Z₄ / Z₂² 2-parts (gross's frame as the
   worked instance); docket5e trisection re-run (same script, different
   quotient pair).
5. **Lean-certify the BZ lane** (verified enumerator or weight-enumerator
   recomputation) — would make every census in §2 kernel-checked.
6. Literature dedup before any external claim: the census specialization
   of BZ and Theorem T are believed new in this combination; the shift
   bound transport is classical; QC spectral bounds (Semenov–Trifonov,
   Lally, Jensen) are related but target d_min, not censuses (assessed
   from working knowledge this session, not re-surveyed).

## 9. Session log

- 2026-08-06: foundations locked (V-series ALL PASS on first full run
  after one wrong expectation: gross base is [[72,12,6]], k = 12 — k is
  (R)-preserved, of course). BZ lane built and validated (grossbase
  0.2 s = (CLASS); f2a6 4.5 s = SAT census exactly). docket37 census
  landed (2,203 classes, 661 s). Shift engine + measurements + gap
  certification same day; ε-trisection derived and validated on the
  docket37 data end-to-end. Lemma F closed trivial.
- 2026-08-06 (later): docket5e census landed (2,371 classes, 870 s) and
  its trisection validated (sectors 1/0/2370, H-census bijection ✓) —
  Theorem T now at 4,574/4,574 clauses across two codes. §6 gates 1–2
  are therefore fully data-backed for BOTH docket cells. Deep dichotomy
  tree (depth 10) left running for the kill-rate curve; its result
  refines §3.5 but changes no verdict.

## Appendix. Verification map

| claim | check |
|---|---|
| conventions, params, census membership | `scripts/a28_v_foundations.py` (ALL PASS) |
| BZ censuses + certificates | `scripts/a28_bz_census.py <code>`; outputs `data/a28/census_*.json`, `bzcert_*.json` |
| grossbase = (CLASS), f2a6 = SAT census | assertions inside the BZ driver (run exits nonzero on mismatch) |
| spectral layer, patterns, game soundness | `scripts/a28_m_f2a6.py` (M0–M4) |
| beam = exhaustive; survivor terrain | `scripts/a28_diag_gap.py` (D1, D2) |
| technique-gap cell certified empty | `scripts/a28_diag_joint.py` part (b) |
| joint-vs-union deltas | `scripts/a28_diag_joint.py` part (a) |
| dichotomy tree + landings + replay | `scripts/a28_dichotomy.py <depth>`; `data/a28/dichotomy_f2a6.json` |
| Theorem T validation (all clauses) | `scripts/a28_trisection.py docket37` |
| Lemma F | `scripts/a28_frobenius.py` |

Engine library: `scripts/a28_lsc_lib.py` (group algebra, registry, F₂
linear algebra), `scripts/a28_spectral.py` (GF(16) characters, orbits),
`scripts/a28_shift_engine.py` (Theorem S game + exact minima),
`scripts/a28_joint_game.py` (Theorem J + replay verifier).
