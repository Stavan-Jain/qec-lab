# A23 — Analytic seam-coset floor for f2a6f17e:y (`SeamCosetFloor 16`)

**Target.** Replace the CryptoMiniSat XOR-native UNSAT@14 verdict behind
`coverData.SeamCosetFloor 16` (QECLean
`QEC/Stabilizer/Framework/Homological/BBDoubling.lean`, consumed by
`QEC/Stabilizer/Codes/BivariateBicycle/Z5Z15F2A6/Distance.lean`, branch
`claude/a15-m-kernel-route`) with an analytic / certificate-style proof.
Instance: base `[[150,8,8]]` over `Z₅×Z₁₅`, `A = 1 + y + x`,
`B = xy⁶ + xy¹⁰ + x²y¹²`; cover `[[300,8,16]]` over `Z₅×Z₃₀`, deck `(0,15)`.

Statement: for every base 2-cycle `ζ` (`∂₂ζ = 0`) and every `f`, if
`seamC ζ + ∂₂f ∉ boundaries` then `|seamC ζ + ∂₂f| ≥ 16`.
Equivalently: each nonzero class in `im δ₂ ⊂ H₁(base)` has coset minimum
≥ 16 (δ₂ = Smith connecting map; `seamC ζ` is its chain-level
representative).

**Session-1 outcome (2026-07-22), headline:** three certificate-style
mechanisms (disjoint dual-detector packing, fractional packing,
region-confined parity floors) are all REFUTED with sharp numeric
certificates of failure (§§2–5) — the seam clutter has a genuine
covering-LP gap (LP ≈ 13.8 < 15).  But a structural reduction (§6)
collapses the whole Prop, quantifiers and all, to ONE self-contained
two-trinomial inequality

    ∀ f ∈ F₂[Z₅×Z₁₅]:  |A⋆f + e₀| + |B⋆f| ≥ 16,

`e₀` = the explicit weight-40 idempotent of the F₁₆ kernel block, and
the small strata of that inequality are verified exactly with slack
(§7).  Full statement of the reduction chain in §6; Lean bridge design
and next-session plan in §8.

Scripts: `a23_seam_calibration.py`, `a23_detector_packing.py`,
`a23_targeted_packing.py`, `a23_fractional_packing.py`,
`a23_lp_estimate.py`, `a23_lp_exact.py`, `a23_column_generation.py`,
`a23_w8_cycle_census.py`, `a23_census_packing.py`, `a23_region_probe.py`,
`a23_band_probe.py`, `a23_transfer_structure.py`, `a23_final_form.py`,
`a23_strata_check.py`.
Data: `data/a23/`.

All numerics in the REPO convention (`conv A f h = Σ_x A(x)f(h−x)`,
`∂₂f = (A⋆f | B⋆f)`, `∂₁c = B⋆c_L + A⋆c_R`, `sec` = the `y∈[0,15)`
window, deck `(0,15)`); cross-validated against the shipped weight-8
`ustar` witness (cycle ✓, non-boundary ✓) from
`data/a17/f2a6_z5z30_lean_data.json`.

## §1 Calibration (2026-07-22) — the seam-coset geometry

`a23_seam_calibration.py`, all checks green:

* `ker ∂₂(base)` has **dim 4** (15 nonzero elements, every one of weight
  40), and forms a **single G-orbit of size 15** with **stabilizer
  `S = ⟨(1,3)⟩` of order 5** (the `x·y³` diagonal — S fixes each kernel
  element *as a chain*).  So one certificate + G-translation covers all
  15 seam classes.  (The docket's "5 orbit reps" figure was the
  hit3/4/6 count; f2a6:y is even better placed: **one** orbit.)
* `δ₂` injective on `ker ∂₂ \ 0` ✓ (15 distinct nonzero classes; matches
  A14 Prop A14.1(2)).  `k = 8`, cycles dim 79, boundaries dim 71 ✓.
* **Raw seam weights `|seamC ζ|`: 18 (×10) and 20 (×5), all even.**
  Coset minima are 16 (SAT-tight), so the raw seam is *not* the coset
  minimum here — unlike gross, where S0 was tight on 18/63 classes.
* Dual side: dual cycles dim 79, dual homology dim 8; the seam class
  pairs odd with exactly **128/255** dual classes (the pairing functional
  restricted to `im δ₂ ≅ F₂⁴` is linear: any dual class pairs odd with
  either 8 or 0 of the 15 seam classes).

**Finding 1 (attack line 1 of the brief is dead as stated).** The
proposed reduction "|w| ≤ 14 ⟹ ∂₂f is light ⟹ apply the 113-class
`LightClassification`" needs the boundary part to have weight
≤ 14 + |seamC ζ|.  But `seamC ζ` is itself a coset element (`f = 0`), so
|seamC ζ| ≥ 16 whenever the floor *holds* — the premise "|seamC ζ| likely
small" is impossible.  Concretely |seamC ζ| ∈ {18, 20}: the reduction
would need a classification of boundaries up to weight **32–34**, far
beyond `LightClassification`'s cutoff 14.  No choice of coset
representative helps (any representative has weight ≥ 16 by the floor
itself; even using the weight-16 tightness witness as offset leaves a
cutoff of 30).  A conditional route through the 113 classes **cannot
close by triangle inequality alone**.  Composing with A22 therefore has
to happen through a different mechanism, not through "the boundary part
is light".

## §2 The disjoint dual-detector mechanism (new attack)

For a *fixed* class, a family `z₁, …, z_N` of **dual cycles** (`z ⊥ im ∂₂`,
i.e. `Ā⋆z_L + B̄⋆z_R = 0`) each pairing **odd** with the offset
(`⟨seamC ζ, zᵢ⟩ = 1`) and with **pairwise disjoint supports** forces every
coset element to meet every `zᵢ` (odd ⟹ nonempty), hence
`|seamC ζ + ∂₂f| ≥ N` for *all* `f` — no case analysis, no boundary
classification.  With N = 15 the parity lemma (all coset weights even;
`boundary_sum_zero` + the four basis seam parities) lifts to the tight 16.

Assets making this plausible here:

* the order-5 stabilizer `S` fixes `ζ` and hence the class, so an
  *S-clean* shape (support differences avoiding `S∖0`) contributes **5
  mutually disjoint odd-pairing translates for free** — 3 shapes with
  disjoint S-saturations would give 15;
* reflect-swap (`z_L(g) = w_R(−g)`, `z_R(g) = w_L(−g)`) maps primal
  cycles/boundaries to dual ones weight-preservingly, so the shipped
  weight-8 base witness yields a weight-8 logical-carrying detector ✓;
* G-covariance: translating a certificate by `g` certifies class `g·c`,
  and the per-class pairing checks can be re-verified directly by
  `decide` — **no transport lemma needed in Lean**.

Lean shape (generic lemma, ~30 lines): `N` disjoint odd-pairing dual
cycles ⟹ coset floor `N`; instance side = per-class finite `decide`
facts (dual-cycle check, pairing, disjointness) + kernel-enumeration
certificate (16 elements of `ker ∂₂`, KernelCert-style) + parity.

**Status of the integer packing hunt:** two rounds
(`a23_detector_packing.py`, `a23_targeted_packing.py`):

* light dual cycles are *scarce*: per-dual-class annealing over all 128
  odd-pairing classes found weight ≤ 10 representatives in only ~9
  shapes (2× w8, 7× w10); most odd-pairing dual classes appear to have
  min weight ≥ 12 (mirrors the bimodal coset-min spectrum of the
  113-class data on the primal side);
* best disjoint family so far: **12** (total support 100/150) — three
  short of the target 15.

## §3 The weight-8 cycle census (standalone structural fact)

`a23_w8_cycle_census.py`, MITM over splits (4,4)/(5,3)/(3,5), exact:

* the base code has **exactly 75 weight-8 cycles, all of split (4,4),
  all logicals, forming a single free G-orbit** (the `ustar` translates);
  no weight-8 boundaries exist.  (Caveat: splits (6,2)/(7,1)/(8,0) not
  yet closed — 201M+ streams; no structural reason to expect them.)
* they occupy 15 homology classes, 5 translates each;
* 40 of the 75 reflect-swap to detectors pairing odd with `c1`, and all
  40 are S-clean, forming **8 S-orbits of 5**.

## §4 Packing refutations (the covering-LP gap)

The complete verdict chain, all on the fixed class rep `c1 = seamC K[0]`
(one G-orbit ⟹ WLOG):

| mechanism | script | verdict |
|---|---|---|
| disjoint packing, w8 universe (complete) | `a23_census_packing.py` | **exact max = 10** (B&B, exact) |
| disjoint S-orbit triples (would give 15) | same | **none** (disjointness graph on the 8 S-orbits has only 4 edges — pairs (0,6),(0,7),(1,5),(3,7) — no triangle) |
| disjoint packing, enriched pool w≤12 (1688 cands) | same | ≥ 12, stalls (time-capped B&B) |
| fractional packing LP, pool w≤12 (1360 cands) | `a23_lp_exact.py` (HiGHS) | **LP = 13.7037** |
| column generation over all dual classes (μ-annealed columns, 12 rounds) | `a23_column_generation.py` | LP creeps 13.70 → **13.84**, violated columns → mass 1: converging ≈ 13.9 |
| PB-SAT multiset certificate M=2, Σm ≥ 29 | `a23_fractional_packing.py` | (superseded: LP < 14.5 makes it UNSAT; run killed) |

**Finding 2 (the packing mechanism is dead, quantitatively).** The seam
clutter (blocker pair: minimal seam-coset supports at weight 16 vs
minimal odd-pairing dual-cycle supports at weight 8) has a genuine
fractional packing/covering gap: max fractional packing ≈ 13.9 < 15 <
16 = min blocker weight.  No detector-counting certificate — integer,
multiset, or LP-rounded — can reach the floor.  The optimal dual
fractional covers (value 13.70, support 92, values in 27ths) are
interesting objects: they block every odd-pairing dual cycle at cost
strictly below the code's own coset minimum.  This is the
Lehman/Seymour binary-clutter non-idealness phenomenon materializing in
a BB seam coset.

## §5 Region-confined parity floors also cap (12)

`a23_region_probe.py`, `a23_band_probe.py`: partition the 150 qubits,
lower-bound `|w ∩ R|` per region by the parity system of dual cycles
*confined* to R (this sees odd intersections and even-pairing
constraints — strictly more than the LP per region):

* S-orbit saturations (40 pts): confined dim exactly 5 = the five
  translates; local floor exactly 5.  Disjoint pairs exist (4 of 28),
  triples don't.  Best partition sum: 5 + 5 + 2 = **12**.
* y-bands: bands shorter than 8 rows carry **zero** confined dual
  cycles (B's y-spread {6,10,12} + A's 1 need ≥ 8 rows); one 8-row band
  fits in 15 rows, dim 9, floor ≤ 2.  A 13-row band has dim 59 and
  annealed floor ≤ 21 > 16 — but a 130-point region is uncertifiable by
  subset enumeration, and it is not a decomposition (single region).
* seam geometry (used later): `sC_L` = 2 (K[3]: 4) points, all in row
  y=0; `sC_R` = 16 points, rows y ∈ {6..11}.

**Finding 3:** confinement destroys dual cycles; local certificates
lose ≥ 4 of the floor.  Every "sum of local reasons" formalism tested
caps at 12–13.9 against the true 16.  The floor's mechanism is
genuinely global — consistent with the deficit-wall picture, and
explains a-posteriori why the docket reached for XOR-native SAT here.

## §6 THE REDUCTION: seam cosets = graph lifts of the kernel block

`a23_transfer_structure.py`, `a23_final_form.py` — every step verified
numerically end-to-end; all of it is elementary linear algebra over
`R = F₂[Z₅×Z₁₅]`.

**(1) Kernel coincidence.** `ker(A⋆·) = ker(B⋆·) = ker ∂₂` (all dim 4 —
the single F₁₆ character block: `Â` vanishes exactly at the Frobenius
orbit of `(χ(x), χ(y)) = (ξ, 1+ξ)`, `ξ` primitive 5th root; `B̂` vanishes
there too, and 4 = 4 forces equality).  Character count: `Â(χ) = 0` has
no solution with `χ(x) = 1`, and for each of the 4 primitive `ξ`,
`1 + ξ ∈ F₁₆* = μ₁₅` always — one orbit of size 4.

**(2) Block idempotent.** The unique `e₀ ∈ ker A` acting as identity on
the block: `|e₀| = 40`, `e₀⋆e₀ = e₀`.  Support: rows `y ≡ 0 (mod 3)`
carry 4 points (all x except the S-diagonal point `x = y/3`), the other
10 rows carry 2 points; exactly 8 of the 15 S-cosets, 5 points each.
`ker ∂₂ ∖ 0` = the 15 G-translates of `e₀` (single orbit, stab
`S = ⟨(1,3)⟩`).

**(3) Quasi-inverse and transfer.** `P` with `P⋆A = 1 + e₀` (|P| = 19),
`Q = B⋆P` (|Q| = 19, support in rows y ∈ {5..11}).  On the seam:
**`sC_R = Q⋆sC_L` exactly (t' = 0)** — the seam offset itself satisfies
the homogeneous transfer relation.  Hence

    w ∈ seam-coset(ζ)  ⟺  w_R = Q⋆w_L  ∧  e₀⋆w_L = e₀⋆sC_L(ζ) ≠ 0.

**(4) Graph form.** `C_graph := {(v, Q⋆v)} = boundaries ⊕ {(β,0) : β ∈
ker ∂₂}` — note `(β,0)` is a cycle since `β ∈ ker B`.  The union of all
15 seam cosets is `C_graph ∖ boundaries`, i.e. every coset element is

    (β + A⋆f | B⋆f),   β ∈ ker ∂₂ ∖ 0.

The Lean-side dictionary needs **no** P, Q, e₀ theory: for each kernel
basis element `K[i]`, `seamC K[i] = (bᵢ + A⋆fᵢ | B⋆fᵢ)` with explicit
`bᵢ ∈ ker ∂₂ ∖ 0`, `fᵢ` (|fᵢ| = 18–22) — four `decide` facts
(`data/a23/final_form.json`) — plus additivity of `seamC` (already in
`BBCover.lean`) and the 16-element kernel enumeration
(KernelCert-style rank certificate).

**(5) Orbit collapse.** `β = g·e₀` and translating `f` by `g` shifts
both blocks: WLOG `β = e₀`.  **Final form:**

    SeamCosetFloor 16  ⟺  ∀ f: |A⋆f + e₀| + |B⋆f| ≥ 16.

**(6) σ-symmetrization.** `σ: x ↦ xy⁶, y ↦ y⁴` is an order-2
automorphism with `B = xy⁶·σ(A)` and `σ(e₀) = e₀` (σ permutes the
kernel characters within their Frobenius orbit).  Substituting `h = σf`:

    ∀ h: |σ(A)⋆h + e₀| + |A⋆h| ≥ 16,   σ(A) = 1 + y⁴ + xy⁶

— an inequality between TWO weight-3 trinomial images of the same `h`,
with a single explicit weight-40 offset.  Tight: annealing independently
finds min = 16 at split `(a,b) = (10,6)` (matches SAT).  Annealed Pareto
frontier (a = |A⋆f+e₀|, b = |B⋆f|): (2,16),(4,16),(5,13),(6,12),(7,11),
(8,10),(9,9),(10,6)*,(11,7),(12,8),… — all totals even (parity), the
tight point isolated at (10,6); the seams themselves sit at (2,16)/(4,16).

## §7 Strata verification of the final form (exact, small strata)

Proof skeleton: if `a ≥ 8` and `b ≥ 8`, done (16, no parity lemma
needed).  Otherwise WLOG (σ) `b ≤ 7`: then `u := B⋆f`... more precisely
`u := A⋆h` lies in the ≤7-layer of `V = im(A⋆·)`, `h = P⋆u + κ` with
block part `κ` invisible (`σ(A)⋆κ = 0`), so the check depends on `u`
alone:

    ∀ u ∈ V, |u| = b ≤ 7:  |Ψ⋆u + e₀| ≥ 16 − b,   Ψ := σ(A)⋆P  (|Ψ| = 19).

Layer sizes (exact, MacWilliams from `V⊥` = reflected block, enumerator
`1 + 15z⁴⁰`): N₂ = 150, N₃ = 4 375, N₄ = 76 200, N₅ = 1 076 250,
N₆ = 12 584 250, N₇ = 124 077 000.  (V has minimum distance **2**: the
150 S-diagonal pairs `{u, u+s}`, `s ∈ S∖0`, have syndrome 0.)

Exact checks (`a23_strata_check.py`, GF(16)-syndrome enumeration):

| stratum b | elements | exact min a | min total | slack |
|---|---|---|---|---|
| 2 | 150 | 24 | 26 | 10 |
| 3 | 4 375 | 21 | 24 | 8 |
| 4 | 76 200 | 20 | 24 | 8 |
| 5 | 1 076 250 | 15 | 20 | 4 |
| 6 | 12 584 250 | **10** | **16 — tight, 0 violations** | 0 |

Zero violations anywhere.  The b=6 stratum carries **exactly 15 tight
pairs (a,b) = (10,6)** — precisely one G-orbit (75 translates / stab 5),
i.e. the SAT witness orbit and nothing else.  b=7 (124M; needs a
numpy-batched enumerator, ~10 min) and the σ-side strata (a ≤ 7,
syndrome-coset layers of the same sizes) remain for the next session —
together they would complete a NO-SAT numerical verification of the
inequality, replacing the CMS UNSAT@14 with structured enumeration.

## §8 Assessment + next session

**What died (with certificates):** the brief's attack line 1 as stated
(|seamC| ≥ 16 always — Finding 1); all detector-counting certificates
(Findings 2–3: exact packing 10, LP 13.84, regions 12).

**What lives:** the final-form inequality is a clean, self-contained,
tight combinatorial statement about two trinomials on `Z₅×Z₁₅`, with
(i) a verified quantifier-collapse reduction from `SeamCosetFloor 16`
that is Lean-cheap (4 decide facts + seamC additivity + kernel
enumeration + translation), and (ii) an exact-enumeration path whose
cost profile matches the a15 branch's proven sweep technology
(C(75,≤7) syndrome-filtered ≈ 2×10⁹ raw cells vs the branch's 560M-mask
sweeps; Gaussian-certificate or packed-Nat native_decide batching are
both plausible).  The elegant target remains a structure theorem for
the ≤7 strata (A-thin lines à la `a17_f2a6_athin_lines.py` — the w≤5
layers are already censused there) that collapses the per-u checks to
per-family lemmas.

**Next session plan, in order:**
1. finish b=6/b=7 exact checks + the σ-side (a ≤ 7) strata (numeric
   completion of the inequality);
2. structure theorem attempt for the ≤7 layers of V (S-diagonal pairs +
   A-lines composition; the layers look like unions of few "cells");
3. Lean bridge: worktree `a23-seam-transfer`; prove the reduction chain
   (§6(4)–(5)) against `coverData`; state the final-form inequality as
   the single remaining certificate Prop;
4. decide brute vs structured discharge of the inequality in Lean based
   on (2)'s outcome.

**Reachability:** good.  This no longer needs new mechanisms — the
remaining risk is Lean-scale engineering of a 10⁸–10⁹-cell sweep (or
the strata structure theorem removing it).
