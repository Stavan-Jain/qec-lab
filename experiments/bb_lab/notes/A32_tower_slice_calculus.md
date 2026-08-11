# A32 — The tower slice calculus: two-level shadow/overflow reduction, and d([[360,12,≤24]]) = 24 at certificate tier

**Session 2026-08-10** (worktree `worktree-bravyi360-fibering-fit`; defensive
registry claim, note lives on this branch). Mission: analytic reductions that
move the open d = 24 obligations of C = Bravyi `[[360,12,<=24]]` (Z₃₀×Z₆,
A = x⁹+y+y², B = y³+x²⁵+x²⁶) off the n = 180 / W = 22 enumeration scale.
The intended product was the **machinery** — a two-level slice/carry calculus
down the tower C → BY → GB with a GB-sector trisection of every remaining
obligation. The calculus turned out to be so much cheaper than the baseline
(~9.4 min total vs the ~40–60 h enumeration estimate) that the entire
remaining program **closed live in this session**:

    d([[360,12,≤24]]) = 24  — certificate tier, end-to-end, no SAT on the
    critical path (X-side; Z-side by the BB transpose duality).

Bravyi et al.'s Table 3 entry "≤ 24" (arXiv:2308.07915) resolves to "= 24";
this exceeds the published solver-exact BB record ([[288,12,18]]) by 25% in
n and 33% in d (A19 §4's record context), at deterministic-certificate
grade. Everything below is falsify-first: each claim was validated against
banked data (A19 (M)@24 census, A24 band censuses, the ISD witness bank)
before entering this note; the verification map is §8. NOT Lean-checked
(§10 item 1).

Scripts `a32_*.py`; data `data/a32/`. Read-only inputs: main-checkout
`data/a19/`, `data/a24/`; worktree `data/a19_scope/`.

## §0 Starting state and the two open obligations (A19/A24 recap)

Tower (A19 §1): GB = (15,3) `[[90,8,8]]` →(x-deck)→ BY = (30,3)
`[[180,8,12]]` →(y-deck)→ C, both lifts twisted, k jumps 8 → 12 at the top
((R) fails on all top decks — A19 §2; the *lower* rung is k-preserving).
Certified before this session: 12 ≤ d(C) ≤ 24. Open for d = 24 after A19
§9–10 + A24 (banked: y-deck (M)@24 floors over stabilizer-shadow census
bands ≤ 20 — 8,310 SAT UNSATs, reproduced solver-free in this worktree as
8,461 rung passes in 735 s, `data/a19_scope/`; A24 reachable band ≤ 16
closures):

1. **SF24-y live strata** (b = p_y(v) a reachable-class BY *logical*):
   (|b|, m₁) ∈ {(18,2), (20,1), (22,0)} — per-b floors are cheap, but the
   *census* of reachable-class logicals at bands 18–22 was the
   n = 180 / W = 22 wall (A24 §3.2 phase-1 engine's band-22 scaling edge;
   coset-BZ ≈ 3.3e13 nodes, ~40–60 h).
2. **Flat band-22 dangerous residue** ([b] = 0 side): |b| = 22 BY
   *stabilizers* with flat top carry (m₁ = 0); the BY band-22 stabilizer
   census had stalled at 151 of ~20k projected classes.

## §1 The scouted map (all four candidate directions)

| direction | verdict | grounds |
|---|---|---|
| (a) tower slice calculus (compose the two deck slices) | **DEVELOP** — became this note | the A24 §2.6 W-shadow theorem is its rank-0 case ("reachability is decided one rung below"); A24's band-16 family split (2 τ_x-lifts + 3 W-lifts + 1 stab-lift) is exactly a GB-sector trisection read off the second slice — the organizing theorem was already visible in banked data |
| (b) ε-trisection port to coset floors (A28 Thm T → logical side) | **MERGED INTO (a)** | for this tower they are the *same* mathematics: BY = Z₂×(Z₁₅×Z₃) and the x-deck fold BY → GB *is* A28's s ↦ 1 quotient; the logical-side port = the cycle/carry version of Thm T = Lemma 1 + Theorem 4 below. Census-side Thm T applies verbatim to obligation 2, where its sectors I and II die by parity (Thm 5) |
| (c) Φ-transfer sparsity theory (GF(16)[H] m-smallness, A19 §9 inputs) | **KILL (superseded for this instance)** | (a) makes the n = 180 censuses it was designed to tame unnecessary; its target regime (A24 §2.4's \|S\| ∈ 9–16 floor minimum) is exactly what the tower bypasses by moving the census to n = 90. Residual value only for codes with *no* second deck (odd-\|G\| bases, e.g. bb_90's Z₁₅×Z₃) |
| (d) parity/deficit-wall structure at the top band | **ABSORBED AS LEMMAS** | the productive fragment is the augmentation parity theorem (Lemma 2 — in scope, \|A\| = \|B\| = 3 odd), which kills the flat-22 β = 0 branch outright and halves every stratum table; the A17-P3 pushforward mechanism (im p₁ = im δ₂) needs (R) and dies at the top rung — checked, not used |

## §2 The calculus

Throughout: lab convention (H_X = [M_A|M_B] rows = X-stabilizers,
H_Z = [M_B^T|M_A^T]; X-logicals = ker H_Z ∖ rowspace H_X); folds
p_y: C₁(C) → C₁(BY) (y mod 3, per block), p_x: C₁(BY) → C₁(GB) (x mod 15);
sheet embeddings e₀, e₁; transfer τ = e₀ + e₁ (full preimage).

**Lemma 1 (deck transport, twist-generic).** For any free Z₂-deck pair with
arbitrarily twisted polynomial lift (some terms carrying the deck element
σ): (i) the fold is a chain map (H_Z(base)∘p = p₀∘H_Z(cover)) and
transports stabilizer row spaces into stabilizer row spaces; (ii) τ
satisfies H_Z(cover)∘τ = τ₀∘H_Z(base), maps base stabilizers to cover
stabilizers, p∘τ = 0, τ∘p = 1+σ; (iii) every cover chain v has unique
sheets (v₀, v₁), shadow b := p(v) = v₀+v₁, and the **slice identity**
|v| = |b| + 2·|v₀ ∧ v₁|; (iv) v is a cycle ⟺ b is a cycle and the **carry
system** E v₀ = R b holds (E = H_Z(cover)∘(e₀+e₁), R = H_Z(cover)∘e₁).
*Proof.* The fold is the ring quotient σ ↦ 1, so lifted polynomials
descend term-by-term (twists die); for (ii), (1+σ) is central and
(1+σ)f̃ũ = τ(fu) for any lifts f̃, ũ — the antipode (transposed circulants)
commutes because σ⁻¹ = σ. (iii) is a per-cell count; (iv) unpacks
(i)+(ii). ∎

**Lemma 2 (cycle parity; scope: |A|, |B| odd).** Every element of ker H_Z
has even weight. *Proof.* The augmentation ε: F₂[G] → F₂ is a ring map and
the cycle condition pairs the blocks through A, B with ε(A) = ε(B) = 1, so
ε(u_L) = ε(u_R) and |u| = |u_L| + |u_R| is even. ∎ (Verified *exhaustively*
per instance: every kernel-basis vector of C/BY/GB/BX is even, and even
weight is closed under sums — `a32_tower_slice.py` Part 1.)

**Theorem 3 (two-level slice).** For any C-chain v, with b = p_y(v),
β = p_x(b), overflows m₁ (y-level) and m₂ (x-level):

    |v| = |b| + 2m₁ = |β| + 2(m₁ + m₂),

and if v is a cycle, so are b and β, with both carry systems in force.
The descent square commutes: descending C → BX → GB (x first) yields the
*same* β and the **overflow square** m₁ + m₂ = m_x + m_y′ (total overflow
is path-independent). [200 random-cycle checks incl. the square; 13
converse lift reconstructions.]

**Rank inputs (measured; third independent derivation).** rank p_y* = 6,
rank p_x* = rank τ_x* = 4; **exactness both ways**: im τ_x* = ker p_x* =:
K_x (BY side) and ker τ_x* = im p_x* (GB side); K_x ⊆ R_y := im p_y*;
W := p_x*(R_y) has dim 2 — its 3 nonzero classes form one translation
orbit; p_x*⁻¹(W) = R_y (the A24 §2.6 reachability theorem). R_y ∖ 0 = 63
classes in 11 translation orbits (6×3 + 5×9, matching A24 §2.2);
K_x ∖ 0 = 15 classes in 3 orbits (3+3+9). All reproduced from scratch in
this worktree's frame conventions (`a32_tower_slice.py` Part 2).

**Theorem 4 (GB-sector trisection of the SF24-y obligations).** Let v be a
nontrivial C X-logical in the y-safe sector ([b] ≠ 0, b = p_y(v)) with
|v| ≤ 22 (even, Lemma 2) — a d = 24 violation candidate. Then
[β] = p_x*[b] ∈ W, and exactly one of:

- **(A) [β] ∈ W ∖ 0.** β is a weight ≥ 14 member of one of the three
  W-cosets (coset minima = 14, §3), and m₁ + m₂ ≤ (22−|β|)/2 ≤ 4.
- **(B) β = 0.** b = τ′(γ) with γ a GB X-logical, [γ] ∉ im p_x*
  (⟺ [b] = τ_x*[γ] ≠ 0, by exactness), |γ| even, |v| = 2|γ| + 2m₁; hence
  (|γ|, m₁) ∈ {(8, ≤3), (10, ≤1)} only.
- **(C) β ∈ Stab(GB) ∖ 0.** |β| even, [b] ∈ K_x ∖ 0 (reachability
  automatic, K_x ⊆ R_y), b is an x-dangerous BY-logical over shadow β
  with m₂ ≥ (12−|β|)/2 (the d(BY) floor), and m₁ + m₂ ≤ (22−|β|)/2.

In all three sectors every cycle-lift v of a logical b is automatically a
nontrivial logical ([b] ≠ 0 ⟹ [v] ≠ 0 by stabilizer transport), so the
per-object certificates are **pure feasibility statements** — no class
dispatch, no triviality windows inside the rungs. ∎ (case split on [β];
the sector-B parity constraint is Lemma 2 on γ.)

**Theorem 5 (flat-22-dangerous trisection).** For obligation 2 (|b| = 22
BY-stabilizer, m₁ = 0): β = p_x(b) ∈ Stab(GB) and

- **β = 0 is impossible**: b = τ′(γ) forces γ ∈ ker H_Z(GB) with
  |γ| = 11 — odd, contradicting Lemma 2. (In A28-Theorem-T terms: census
  sectors I *and* II are parity-dead for this residue.)
- **β ≠ 0**: |β| ∈ {6, 10, 12, …, 22} (no weight-8 GB stabilizers — §3),
  m₂ = (22−|β|)/2 *exactly*, so b is an x-lift of a censused GB stabilizer
  with pinned overflow — per-β bounded fiber + flat-top rung (M = 1) per
  lift. ∎

The stalled n = 180 band-22 stabilizer census is thereby **replaced** by
the GB stabilizer census at ≤ 22 (n = 90, seconds) + per-β fibers — and
more generally (§5, census-completeness bonus) the *entire* BY stabilizer
census at ≤ 22 is re-derived solver-free as the union of these fibers plus
the analytically-determined β = 0 family τ′(Stab(GB) ≤ 10).

## §3 GB-level censuses (n = 90; the workhorse layer)

`a32_gb_census.py` (2.6 s) + the ≤ 22 extensions in `a32_sectorAC_full.py`,
all via the a30 coset-BZ C kernel over two disjoint 41-column information
windows (complete to W = r₁+r₂+1; node counts exact binomial sums,
kernel-asserted):

| species | depth | result (translation orbits) | nodes |
|---|---|---|---|
| logicals, all 255 classes | ≤ 10 | w8: 2, w10: 36 — **all 38 outside im p_x***; **d(GB) = 8 census-complete** (formerly SAT ladders) | 2.5e8 |
| W-cosets (3 classes) | ≤ 16 | w14: 6, w16: 68; **minima = 14 REPRODUCED** (A24 §2.6's SAT → BZ) | 4.5e8 |
| W-cosets | ≤ 22 | {14: 6, 16: 68, 18: 1,627, 20: 19,873, 22: 175,057}; 8,848,365 vectors (count corrected +2 per the A33 §8 empty-window erratum — orbit table unchanged, both extra elements' cap-0 fibers empty; `census()` fixed and this row re-run post-fix) | 1.9e10 |
| stabilizers | ≤ 16 | {6: 1, 10: 6, 12: 21, 14: 64, 16: 333} — **no weight-8 stabilizers**: 4th instance of the Prop-10 weight-gap pattern (gross, A20-Y4, BY, now GB) | 1.5e8 |
| stabilizers | ≤ 22 | {…, 18: 1,733, 20: 10,602, 22: 64,619}; 3,481,595 vectors | 6.3e9 |

Notable structure: **every** weight ≤ 10 GB logical lies outside im p_x* —
the "light logicals concentrate outside the pushforward image" phenomenon
A24 measured one rung up repeats at the bottom rung.

## §4 Live closures, round 1 (`a32_subclosures.py`, 1.1 s)

Per-object engine: the x-deck **lift fiber enumerator** (restricted MITM
lane over off-support columns of the level-2 carry system, complete to
off-support size 4 by the exact-subset-sum argument; every solution
re-verified; deck-translate v₀ ↦ v₀+β quotiented) + the **top rung** =
`scope_bravyi_rung.BravyiRungCell.rung` (a30 architecture: restricted
meet-in-the-middle + coset-BZ lanes; every violation candidate re-verified
end-to-end). For logical b, PASS is exactly m₁(b) ≥ M.

- **6a sector B closed.** 38 rungs (36 |γ| = 10 reps at M = 2; 2 |γ| = 8
  reps at M = 4, re-deriving A24's band-16 diagonal closure): **ALL PASS,
  0 violations.** Cross-check: the two |γ| = 8 diagonals coincide
  orbit-exactly with A24's two banked band-16 B-classes.
- **6b sector A, |β| = 14 closed.** 6 fibers, 5–13 lifts each (49 rungs):
  **ALL PASS.** Both falsify-first anchors hit: m₂ = 0 fibers EMPTY (a
  flat W-14 lift would be a reachable weight-14 BY logical, contradicting
  A24's band-14 census), and the m₂ = 1 lifts reproduce **exactly** the 3
  banked band-16 A-classes (orbit-level match).
- **6c sector A, |β| = 16 closed.** 68 fibers, 246 lifts (m₂-histogram
  {1: 22, 2: 46, 3: 178}): **ALL PASS.**
- **Covariance spot-check**: a translated fiber has the identical
  m₂-profile (transport soundness of per-orbit dispatch).

## §5 Live closures, round 2: everything else (546 s + 13.5 s + 0.04 s)

**`a32_sectorAC_full.py`** (shallow lanes; batched translation-orbit
canonicalization for the large censuses; 545.8 s wall):

- **Sector A at |β| ∈ {18, 20, 22} closed**: 196,557 fibers, caps
  (2, 1, 0); only 13,109 lifts exist at all (m₂-histogram
  {0: 6,929, 1: 4,844, 2: 1,336} — ~93% of fibers are carry-infeasible at
  small overflow, the structural reason the tower wins), **all 13,109
  rungs PASS, 0 violations**. With §4: **SECTOR A CLOSED at every shadow
  weight** (|β| ≥ 24 is free: |v| ≥ |β|).
- **Sector C at |β| ∈ {14, 16} closed**: 397 fibers (caps 4, 3), 4,132
  lifts, 2,770 rungs (logical lifts + 765 flat-22 rungs) **ALL PASS**;
  1,362 stabilizer lifts at |b| ≤ 20 checked present in the banked A19
  census.
- **Sector C at |β| ∈ {18, 20, 22} closed**: 76,954 fibers (caps 2, 1,
  0), 38,223 lifts, 31,613 rungs **ALL PASS** — of which **25,492
  flat-22-dangerous rungs**; 6,610 banked-census membership checks.

**`a32_deep_fibers.py`** (the predicted hard tail — measured cheap):
extends the lift enumerator to off-support size 8 by ordered-split MITM
(validation gate: reproduces the size-4 lane EXACTLY on shallow fibers
before any deep fiber runs). The 28 deep fibers (|β| ∈ {6, 10, 12}, caps
8/6/5) ran in **13.5 s total** (hexagon: 3.5 s — a cell the A24-era
planning priced at solver-hours):

- 2,371 lifts; 2,030 rungs **ALL PASS, 0 violations** — of which **895
  flat-22-dangerous rungs** and 1,135 logical-lift rungs at
  |b| ∈ {18, 20, 22}.
- 341 stabilizer lifts at |b| ≤ 20: all present in the banked census.
- **0 light logical hits**: no deep fiber produces a reachable-class
  logical at |b| ≤ 14, and all |b| = 16 lifts are stabilizers — so the
  reachable band-12/14 EMPTINESS and the band-16 reachable census are
  **independently re-derived by the tower** (A24's engine censuses demote
  from dependencies to cross-checks).

**`a32_dby_floor.py`** (the calculus one rung down; 0.04 s): a nontrivial
BY logical of weight ≤ 10 would be a lift in one of 45 GB fibers (38
logical-shadow at caps 1/0 by d(GB) = 8; hexagon at cap 2 + six w10
stab-shadows at cap 0; β′ = 0 dead: 2|γ| ≥ 16; parity kills odd weights).
Measured: 7 lifts, all stabilizers (exactly the predicted τ′-family
profile), **0 logicals** ⟹ **d(BY) ≥ 12 solver-free** (+ the explicit
banked weight-12 logical re-verified ⟹ = 12). This was the last SAT-tier
input on the d(C) = 24 critical path.

**Census-completeness bonus.** The fiber sweep enumerates, up to
translation and sheet flip, *every* BY stabilizer of weight ≤ 22 (its
GB-shadow is a censused β ≠ 0, or it is τ′(Stab(GB) ≤ 10) — the β = 0
family, = the banked '12,0' + six '20,0' records, verified present). All
8,313 enumerated stabilizer lifts at |b| ≤ 20 were found in the banked
A19 class list ⟹ **the banked census's completeness at bands ≤ 20 is now
solver-free** (previously: the lex-leader SAT engine's terminal UNSATs),
and the band-22 stabilizer census — the stalled object — is *generated*
by the sweep (27,152 flat-22 rungs = its members, all floor-passed).

### The assembly

For any nontrivial X-logical v of C with |v| < 24 (so |v| ≤ 22, Lemma 2),
with b = p_y(v):

1. [b] ≠ 0 ⟹ Theorem 4 sector A, B, or C ⟹ excluded by §4–§5
   (49,855 rungs, all PASS; census + fiber completeness per §2–§3).
2. [b] = 0, b ≠ 0 ⟹ b is a BY stabilizer: |b| ≤ 20 ⟹ excluded by the
   banked (M)@24 floors (solver-free scope rung reruns; census
   completeness per the bonus above); |b| = 22 ⟹ m₁ ≥ 1 gives |v| ≥ 24,
   and m₁ = 0 is excluded by Theorem 5 + the 27,152 flat-22 rungs.
3. b = 0 ⟹ v = τ_y(u) with [u] ≠ 0 (stabilizer transport), |v| = 2|u|
   ≥ 2·d(BY) = 24 (d(BY) = 12 now solver-free).

Hence d(C) ≥ 24; the verified weight-24 τ-lift witnesses (A19 §3; 100
weight-24 ISD vectors re-verified as nontrivial logicals in Part 4c) give
d(C) ≤ 24. **d(C) = 24.** Z-side by the BB transpose duality (antipode
relabeling + block swap). ∎

**Claim tier, stated exactly**: deterministic certified computation —
BZ censuses with exact node-count invariants, MITM fiber enumerations
complete by the exact-subset-sum argument with every solution re-verified,
a30-architecture rungs with end-to-end re-verification of every candidate
— the same tier as A30's d = 20 results and strictly more auditable than
SAT-without-DRAT. SAT survives only in historical cross-checks (A19's
original UNSATs, A24's engine censuses), not in the critical path.
NOT kernel-checked; Lean packaging is follow-on 1 in §10.

## §6 Cost model: before → after (all measured)

| obligation | before (A19/A24 architecture) | after (tower calculus) | outcome |
|---|---|---|---|
| reachable-logical census, bands 18–22 (the SF24-y wall) | n = 180, W = 22 coset-BZ ≈ 3.3e13 nodes (~40–60 h), or the A24 cell engine's band-22 C(35,≤9)-per-cell edge | n = 90 censuses (2.6e10 nodes, ~seconds) + 274k fibers at ms-scale; 93–97% of fibers carry-infeasible | **closed**, ~10 min |
| per-shadow m₁-floors, bands 18–22 | λ/rank certificates per censused b (A24 §3.2 phase 2; de-risked but census-blocked) | 49,855 top rungs at M ≤ 5, restricted lanes, ~2–10 ms each | **closed** (0 violations) |
| flat band-22 dangerous residue | BY band-22 stabilizer census ~20k classes (stalled at 151) + per-class flat queries near the free-UNSAT wall (A19 §10's "one query" died with its session) | β = 0 branch: **0 compute** (parity); β ≠ 0: pinned-overflow fibers, 27,152 flat-top rungs at M = 1 | **closed** |
| deep fibers |β| ∈ {6, 10, 12} (the predicted hard tail) | (priced at solver-hours by the A24-era per-b epistemics) | ordered-split MITM: **13.5 s** for all 28 | **closed** |
| d(BY) = 12 floor | floor@11 SAT ladders, 3.1 h (A19 §1) | 45 fibers, 0.04 s | **solver-free** |
| d(GB) = 8 floor | coset-SAT ladder, 187 s | BZ census ≤ 10, < 1 s | **solver-free** |
| BY stabilizer census completeness ≤ 20 | lex-leader XOR SAT engine (A19 §9, hours) | fiber-sweep membership (8,313 checks) + β = 0 analytics | **solver-free** |
| (M)@24 floors, bands ≤ 20 | 8,310 SAT UNSATs (153 s) / 8,461 solver-free rungs (735 s, pre-session scope rerun) | consumed as-is | banked |

## §7 What is C-specific and what generalizes (the mission deliverable)

**Fully portable (any BB code, any free Z₂-deck tower, twisted or not, no
(R) anywhere):**

1. Lemma 1 + Theorem 3 — transport, slice, carry, overflow square. The
   twist-genericity matters: *both* rungs of this tower are twisted and
   the calculus never sees it.
2. The trisection *shape* of Theorems 4/5: the case split on the second
   shadow ([β] nontrivial / zero / stabilizer) is exhaustive for any
   two-deck tower; only *which classes occur* is instance data.
3. Lemma 2 wherever |A|, |B| are odd (scope provably necessary — A11 E3).
4. The iterated-rung engine: censuses at the bottom rung (smallest n),
   bounded-overflow lift fibers per censused shadow, top rungs per lift —
   per-level completeness composing. Plus the two bonus patterns: distance
   floors of the *middle* code re-derive by the same fibers at a smaller
   budget (§5 dby), and the middle code's stabilizer census re-derives as
   the fiber union (census-completeness bonus) — i.e., **one machine
   serves the distance floor, the census, and the sector floors at every
   rung**.

**Instance-measured (cheap rank data per tower):** the subspace lattice
(R_y, K_x, W) and the single containment K_x ⊆ im p_y* that makes
reachability "decided one rung below" (p_x*⁻¹(W) = R_y follows by
dimension count). When it fails, the trisection survives but sectors B/C
lose their automatic-reachability pruning.

**C-specific:** the numerical values (minima 14, band caps, fiber
multiplicities) and the parity kill of the *entire* β = 0 branch of the
band-22 residue (uses 22 ≡ 2 mod 4; a band ≡ 0 mod 4 analogue would need
the [γ] ∉ im p_x* class argument instead).

**Next instances** (value order): (1) `[[288,8,20]]` (A20's cheapest
d ≥ 20 target — an all-(R) tower with bases at n = 144/72; the same
trisection would organize its (M)@20 census program); (2) the
`[[360,4,20]]` doubled codes of A30 (rung-2 re-doubling = three-level
towers); (3) `[[756,16,≤34]]` y-rung (n = 378 base — a tower census is
the only plausible route); (4) retrospective: gross over bb72 over (3,6)
(would reproduce the slot-frame dangerous census from n = 36 data — a
teaching-doc chapter).

## §8 Verification map

| claim | check |
|---|---|
| Lemma 1 (all transports, both decks, twisted) | `a32_tower_slice.py` Part 0: matrix identities + per-row stab transport + twist-invariance asserts |
| Lemma 2 exhaustive | Part 1: kernel bases of C/BY/GB/BX all even-weight |
| Theorem 3 + overflow square | Part 1: 200 random cycles, both descent orders; 13 converse lifts |
| rank inputs / exactness / W-theorem | Part 2: independent re-derivation (asserts); orbits 11 = 6×3+5×9; K_x 3 orbits; W one 3-orbit |
| banked (M)@24 census decomposes | Part 3: 8,461/8,461 transport + slice asserts; shadow-orbit compression 8,461 → 3,925; per-band table in `tower_validation.json` |
| Theorem 4 vs banked light logicals | Part 4a: band-16 sector split (3A + 2B + 1C) reproduced with exact (|β|, m₂); Part 4b: all 25 light unreachable logicals have nontrivial non-W GB shadows |
| ISD witness bank in tower coordinates | Part 4c: all 100 weight-24 witnesses classified; tight (m₁ ≥ 1) witnesses sit at total overflow M ∈ {2,3} — one unit above the excluded strata (certificates tight with slack 1) |
| GB censuses complete | BZ two-window r-pairs, exact node counts asserted; W-minima = 14 re-derived; ≤ 22 runs re-match the ≤ 16 orbit rows |
| sector B closure | `a32_subclosures.py` 6a: 38/38 PASS + A24 B-class orbit match |
| sector A closure | 6b/6c + `a32_sectorAC_full.py`: 6 + 68 + 196,557 fibers, 13,404 rungs PASS; m₂ = 0 empties; banked A-class match |
| sector C closure | `a32_sectorAC_full.py` + `a32_deep_fibers.py`: 397 + 76,954 + 28 fibers, 36,413 rungs PASS; 8,313 banked membership asserts; deep-enumerator validation gate (== size-4 lane) |
| flat-22 residue | Theorem 5 parity (β = 0) + 27,152 flat-top rungs PASS |
| d(BY) = 12 solver-free | `a32_dby_floor.py`: 45 fibers, 7 stabilizer lifts, 0 logicals; witness re-verified |
| translation covariance of fibers | spot-check (identical m₂-profiles under translation) |
| top-rung soundness | inherited: scope run's m12/A19 §8 reproduction + a30's f2a6 113/113 + planted-control validation |

## §9 Falsified-claims ledger (session-internal)

- "im τ_x* ≠ ker p_x* on this tower" — briefly believed during Part 2
  debugging; FALSE: the subspaces are equal, the comparison was
  ordering-sensitive (`rref_ints` returns an unsorted reduced basis).
  Lesson: span comparisons must sort.
- "The banked stabilizer census should compress ~an order of magnitude
  under GB-shadow quotienting" — measured: only 2.2× (8,461 → 3,925
  orbits). The tower's win is NOT compression of existing censuses; it is
  that census *generation* moves to n = 90 and ≥ 93% of fibers are
  carry-infeasible. Do not sell the calculus as a compression statement.
- "The deep sector-C fibers (hexagon at cap 8) are the expensive tail,
  plausibly solver-hours" — session-internal pricing, REFUTED by
  measurement: 13.5 s for all 28 (the carry system's overdetermination
  dominates the combinatorics there too).
- "Sector B might need censuses at |γ| = 9, 11" — dead on arrival by
  Lemma 2 (kept as the example of parity halving strata tables).
- "W-coset orbit counts grow ~3.5×/band like stabilizer censuses" —
  measured growth is much faster at the top (68 → 1,627 → 19,873 →
  175,057; ~9–24×/band): weight 22 at n = 90 approaches the bulk. The
  fiber layer absorbed it (emptiness rate rises in parallel), but census
  *storage*, not enumeration, becomes the scaling frontier for deeper
  towers.
- (Inherited, respected:) A23-style site sweeps at the top level and
  packing certificates were not re-attempted (A24 §2.4 measurement / A23
  refutation stand); no SAT witness weights were reported as floors.

## §10 Residue and next steps

1. **Lean packaging** (the natural follow-on; nothing here is
   kernel-checked). All certificate species have shipped analogues:
   BZ counting certificates (A28 §6.3's enumerator shape), pivot/rank
   certificates (KernelCert), rung dispatches (A15/A30 pattern), plus the
   new species — the lift-fiber enumeration (a `native_decide`-friendly
   bounded search per fiber; 274k fibers is large but the per-fiber
   statements are tiny). Strategy question for the S-track: certify the
   *assembly* with the fiber layer as data-carried certificates, vs
   re-run-in-Lean.
2. **Witness-side polish**: bank one explicit weight-24 logical + its
   verification data as the canonical d ≤ 24 certificate alongside the
   floor (currently: the ISD bank's 100, re-verified in Part 4c).
3. **Port to `[[288,8,20]]`** (first non-C instance; tests every
   portability claim in §7 on an (R)-holding tower).
4. **Deficit-wall postscript**: with d(C) = 24 = 2·d(BY), the top rung is
   a *perfect* doubling rung despite (R) failing on every deck — the
   first such instance; the A17 wall/deficit taxonomy (built on (R)
   rungs) now has a deck-nontrivial data point. Worth a paragraph in the
   teaching doc's outlook and a check against the A19 §6 open question
   (deficit rungs + deck-nontrivial top as a design principle).
5. **A19/A24 bookkeeping**: A19 §10's "one query" is superseded (the
   near-flat gate's stratum decomposition is exactly what §5 executed,
   solver-free); A24's SF24-y is closed by this note — its phase-1 engine
   and λ-certificate designs remain valuable for odd-|G| instances where
   no second deck exists.
