# A42 — the spectral/module lane on the tour-de-gross column

**Charter.** Close (part of) the lower gap ⌈3r/2⌉ ≤ d(C_{r,1}) ≤ 12r
(A40 §7.1/§8.1) ALGEBRAICALLY: the common-zero variety V(A,B), its
torsion points and char-2 scheme structure (multiplicities), CRT +
chain-ring decompositions of the cylinder/torus group algebras, and
per-spectral-point weight floors (BCH/trace-count species bounds).
Sibling lane to A40 (same branch, same family, disjoint files); the
combinatorial transfer/wall lane stays with A40.  House rules: claim
tiers exact, witness weights are upper bounds only, falsified claims
logged first-class, every consumed vector re-verified end-to-end.

Family (A40 §1, verbatim from arXiv:2506.03094): fixed Laurent pair
A = 1 + y + x³y⁻¹, B = 1 + x + x⁻¹y⁻³ on Z_ℓ × Z_m; members
(ℓ, m) = (6(r+b), 6r), b ∈ {0,1}; conjectured d = 6(2r+b−1).

Unit-cleared polynomial forms used throughout (same ideals, same
torus zeros): Ã = y·A = y + y² + x³ (monic in y, deg 2),
B̂ = x·y³·B = 1 + (x+x²)y³ (deg_y 3, lc = x+x²).

Tooling: `scripts/a42_lib.py` (exact F₂[x] arithmetic, F_{2^d} towers,
Artin–Schreier quadratic extensions, the variety enumeration, the
local-multiplicity engine, the frame k-formula), consumed by staged
scripts `a42_s*.py`; data in `data/a42/`.

---

## §1 S0 — the variety anchors (certificate tier, all gates green)

`a42_s0_variety.py` (3.7 s, `s0_variety.json`, log `s0_run.log`);
anchors `a42_s0b_anchor_checks.py` (`s0b_anchors.json`).

### §1.1 The variety table (exact, over F̄₂)

Res_y(Ã, B̂) recomputed from scratch (5×5 Sylvester over F₂[x]) equals
the banked A40 S3 certificate 1 + x + x² + x⁴ + x⁵ + x¹¹ + x¹³ (G1),
and factors as

> **Res = (x² + x + 1)³ · (x⁷ + x⁶ + x⁴ + x² + 1)**

(deg 2·3 + 7 = 13 ✓).  Since Ã is monic in y, every resultant root
carries ≥ 1 torus point (product formula), and the full variety is
**three Frobenius orbits** (G2, each verified against both defining
equations end-to-end in its residue field):

| orbit | rep (α, β) | residue field | D = [F₂(α,β):F₂] | (ord α, ord β) | plane mult |
|---|---|---|---|---|---|
| **T** (diagonal) | (ω, ω) | F₄ | 2 | (3, 3) | 1 (transverse) |
| **G** (antidiagonal) | (ω, ω²) | F₄ | 2 | (3, 3) | **2 (tangency)** |
| **W** | (ξ, ξ⁹¹) | F₁₂₈ | 7 | (127, 127) | 1 (transverse) |

ξ = root of x⁷+x⁶+x⁴+x²+1 (primitive, ord 127); β_W's min poly is
x⁷+x⁶+x⁵+x⁴+x²+x+1; the discrete log is β = α⁹¹ (s0b).  Scheme
consistency (G3): per x-factor, Res-exponent = Σ (D/deg)·mult exactly
((x²+x+1)³ = 1 + 2 from T + G; the deg-7 factor: 1 from W).

**The tangency, by hand (verified by the engine):** at (ω, β) with
x = ω(1+u), y = β(1+v): Ã_loc = u + βv + O(2), B̂_loc = ωu + v + O(2);
det [[1, β], [ω, 1]] = 1 + βω, which vanishes iff β = ω² — so the
antidiagonal point is a simple tangency (mult 2), tangent line
u = ω²v, and the diagonal point is transverse.  This is the exact
scheme-theoretic content behind "the ω-eigenclass" of A40 §10.9.

### §1.2 The spectral k-formula (theorem + machine validation)

For any frame G = Z²/L, L = ⟨(ℓ,0),(d,p)⟩:

> **k = 2 · Σ_{orbits O: α^ℓ = 1 ∧ α^d β^p = 1} D_O ·
> dim_{k(P)} k(P)[G₂]/(A_loc, B_loc)**

where G₂ is the 2-Sylow part of G and A_loc, B_loc are the images
under x ↦ α·x₂, y ↦ β·y₂ (odd/2-part CRT of F₂[G]; k = 2·dim R/(A,B)
because group algebras are self-injective, so rank H_X = rank H_Z =
|G| − dim R/(A,B)).  Exactness of truncation caps in the engine:
colength-c ideals contain m^c (Nakayama chain), and every orbit here
has c ≤ 2, so any live-axis truncation ≥ 2 is exact; ramified
(one-axis) frames use cap 2⁵ with a stabilization assert.

**Validation (G4–G6, all PASS):**
- 240 frames of the banked triage grid (ℓ ∈ {12,18}, p ≤ 8, all
  shears): formula == TowerCode on honestly-cancelled transported
  supports, every frame.
- 32 spot frames at ℓ ∈ {24,30}, p ∈ {3,5,6,7}: all match.
- Rectangular members: k = 12 at ALL (both columns) r ≤ 8,
  direct-rank-verified r ≤ 6 (n up to 3024).
- The 127-line certified at code level: frame (127, 1, 36)
  (36 = −91 mod 127) has k = 14 = 2·7 by formula AND TowerCode
  (n = 254); shear controls d = 37, 99 give k = 0 (s0b).

### §1.3 The k-laws of the family (∀r, theorem tier)

The shape table (G5 + `f4_shape_table`) gives the local dims at every
2-part shape: **T contributes 1 always** (any shape, including odd
axes); **G contributes 2 iff both axes are even (a₁, a₂ ≥ 1), else
1**; W contributes 1 whenever its character condition holds.  Since
members have 6 | ℓ, m (both axes even) and the truncation-redundancy
argument (§1.2) makes the dims shape-independent beyond a_i ≥ 1:

> **k-law.** For every member: k = 2·(2·1 + 2·2) = 12 + 2·7·[W on].
> W is on iff ord α = 127 | ℓ and ord β = 127 | m, i.e. (CRT with
> s₁ = t₁ = 127 coprime to 6): **b = 1: never (r ≡ 0 and ≡ −1
> mod 127 are incompatible) — k = 12 for the WHOLE b = 1 column,
> ∀r.  b = 0: k = 12 except at 127 | r, where k = 26** (first at
> (ℓ,m) = (762,762), r = 127, n = 1,161,288).

Claim tier: theorem (the formula is proven; the finite inputs are the
variety table and the two shape values, machine-certified; the
127-line's reality is TowerCode-certified at (127,1,36)).  The
(762,762) value itself is formula-tier (n too large for direct rank
here); flagged for an independent check.

**Consequence for the conjecture's shape:** any d-law argument for
the b = 0 column must break at 127 | r (the code there has 14 extra
logical dimensions living on the W-line, with their own — a priori
smaller — weights); the b = 1 column is spectrally uniform ∀r.  The
paper's k = 12 family row, if intended ∀r at b = 0, is FALSE at
r ≡ 0 (mod 127); at b = 1 it is a theorem.  (A40 §2 verified k = 12
at r ≤ 6 only; no contradiction with anything banked.)

**Literature position (agent sweep, 2026-08-28; A31 discipline).**
Liang–Liu–Song–Chen (arXiv:2503.04699, PRL 2025) analyzed THIS pair
(their BB(1̄,1̄,3,3̄)) by Gröbner aggregates: k_max = 26 and anyon
periods (ℓ₀,m₀) = (762,762) — exactly consistent with (and
anticipating, in aggregate form) the k = 26 switch-on above; they
also prove the gcd-periodicity law k(ℓ,m) = k(gcd(ℓ,ℓ₀), gcd(m,m₀)).
What §1 adds over 2503.04699 / Postema–Kokkelmans 2502.17052 (factor
sums + general Gröbner staircase) / Wang–Pryadko 2606.08771
(resultant multiplicity counting, surface setting) /
Eberhardt–Steffan 2407.03973 (semisimple CRT, odd ℓ,m only): the
per-point decomposition with 2-adic local multiplicities at
char-dividing frames, the tangency identification, the per-orbit
switch-on arithmetic (b = 1 blocked ∀r), and the species dictionary.
Frame the k-formula as a refinement unifying those four, NOT as the
first variety connection.  **Sourcing flag for A40:** the sweep's
full-text probes of arXiv:2506.03094 found only gross + two-gross
defined there ((12,6) and (12,12)) — no general (r,b) family, no
d = 6(2r+b−1) — in tension with A40 §1's "sourced verbatim"; A40
owns §1, flagged for re-verification.  ((12,12) = two-gross with
certified d = 18 also calibrates §2's gates.)

### §1.4 The species dictionary (the §10.9 "ω-eigenclass" answered at k-level)

- **W7 chirality is tangency chirality** (s0b): the frames (ℓ,7,ℓ−2)
  that carry the weight-8 winding species select exactly the TANGENT
  orbit G (character condition α^{ℓ−2}β⁷ = 1 forces β = α²); the
  empty +2 mirror selects the transverse orbit T.  Both have k = 4;
  the chirality lives in the minimum WEIGHT of the eigenclass, not in
  k — the sharpened form of A40's open chirality question, now
  pinned to: *the tangent F₄ class carries weight 8; the transverse
  one floors ≥ 14 (censused); explain via the local structure.*
- The banked "ω-point z²+z+1 at every k = 4 shear" (A40 §10.9 item
  3) is the F₄ pair-cluster seen univariately; the bivariate table
  above refines it into T vs G and adds the W-line nobody had seen
  (invisible at p ≤ 8: it needs 127 | p or 127 | ℓ-side conditions).

### §1.5 Finding for A40: 10 phantom triage rows (artifact, flagged)

G5's diff pass: exactly 10 rows of `data/a40/s4_phase_triage.json`
disagree with the true k, all at support-collision shears, all
banked k = 4 vs true k = 0: (12,1,0), (12,1,9), (12,2,9), (12,3,1),
(12,3,2), (18,1,0), (18,1,15), (18,2,15), (18,3,1), (18,3,2).
Cause: `a40_s4_phase_triage.py`'s `tr()` transports supports into a
frozenset, merging colliding terms instead of cancelling mod 2 (at
(12,1,0), y ↦ 1 makes A ↦ x³, a unit — true k = 0).  Direction is
conservative everywhere (banked ≥ true): A40's censuses at those
frames were unnecessary, not unsound; no true-k > 0 frame was
missed.  Fix chip spawned for the A40 lane owner; not touching A40
files from this lane.

---

## §2 The Stage-1 program: the ∀p compact-phase floor via the ω-tower

Target (L-P of A40 §9.7): no nontrivial x-compact phase of the
period-p straight cylinder (either lane) has weight < 2p, ∀p ≥ 2.
The spectral decomposition collapses this to a two-parameter family.

### §2.1 Theorem A (barren periods) — proof sketch, to be written up

y^p − 1 = Π_{d | odd(p)} f_d^{2^{v₂(p)}} over F₂; the cylinder ring
CRT-splits into Λ-Laurent factors Λ_f[x^±], Λ_f = F₂[y]/(f^{2^a}).
For a factor whose root η has A(·,η), B(·,η) coprime in F_Q[x^±]
(⟺ no variety point with that β = η — only ord η ∈ {1, 3, 127}
fail, and ord η = 1 gives the unit x³), the residue Bezout identity
lifts (nilpotent correction) to uÃ + vB̂ = unit over Λ_f[x^±], and a
unimodular row kills homology outright:

> **Unimodular lemma.** Over any commutative ring, if uA + vB = 1
> then ker(f,g ↦ Af+Bg) = im(h ↦ (Bh, Ah)): given Af = Bg, take
> h = vf + ug; then Bh = f + u(Af+Bg) = f and Ah = g + v(Af+Bg) = g.

> **Theorem A.** For every p with 3 ∤ p and 127 ∤ p, the period-p
> compact cylinder homology vanishes identically — no nontrivial
> compact phases exist at ANY weight, both lanes, all 2-parts.

This explains the atlas's emptiness at p ∈ {2,4,5,7,8} structurally
and settles all non-multiples of 3 (and of 127) forever; the old
"unbounded chain depth" wall is dodged because nilpotency of the
2-part radical is all that is used.  For 3 | p (and 127 | p), the
homology concentrates on the (y²+y+1)^{2^a}-factor (resp. the
127-factor): **the ∀p problem reduces to the ω-chain family Λ_a,
indexed only by a = v₂(p), plus classical block-weight arithmetic.**

### §2.2 S1 executed state (2026-08-28): instruments, falsifications,
### realized 2p objects

Machinery (`a42_s1_cylfloor.py`, `a42_s1_syzygy.py`,
`a42_s1_omegafloor.py`; data `s1_*.json`):

- **The windowed cylinder engine** (X-sector convention): exact
  cycle space + exact cylinder-boundary space on any window (the
  no-telescoping argument: B̄'s and Ā's extreme-x coefficients are
  unit multiples, so window-supported boundaries have
  window-bounded trivializers), hence EXACT cylinder-triviality
  functionals; CMS-SAT minimum-weight search over nontrivial
  classes, upward (each UNSAT banks "floor > w"); every find
  re-verified (cycle, pairing, embedding torus).  Calibration: p=3
  min = 6 ✓ (= d((21,3)) banked), p=6 min = 12 ✓ (atlas), p ∈
  {5,7}: dim Z = dim B on windows — **Theorem A machine-checked**.
  A control subtlety worth banking: the p=3 embed-check at Lx = 19
  first returned "torus-trivial" — correct, since 3 ∤ 19 blocks the
  ω-line's switch-on; at Lx = 21 nontrivial ✓.  Spectral switch-on
  arithmetic must pick the control frames.
- **The ω-window engine** over Λ_a = F₂[y]/((y²+y+1)^{2^a}):
  syzygies/trivials/classes of (Ã, B̂) with exact class counts:
  4 (a=0), 6 (a=1), 6 (a=2) — matching the full-ring window classes
  at p = 3, 6, 12 exactly (the ω-factor carries ALL classes;
  Theorem A's CRT concentration, seen in data).
- **The (12,12) gap probe**: 12/12 CMS-enumerated weight-18 minimal
  logicals of two-gross have cyclic x-gaps ≤ 2 AND y-gaps ≤ 2 — no
  Lemma-K windowing exists; two-gross minima are genuinely toroidal
  and give NO cheap compact phase at p = 12.

**Falsified in S1 (first-class):**
1. **The S-law S_ℓ(a) = 3·2^{a−ℓ} is FALSE**: exhaustive
   pattern search gives min slots = 3 at EVERY level 0..2^a−1 for
   a ≤ 2 (always the residue pattern f@{0}, g@{0,1}).  No doubling.
2. **The naive product assembly min_ℓ S_ℓ·c_ℓ is FALSE**: with
   valuation-bucketed costs it predicts 6 at p = 6 (true floor 12).
   Costs must be charged per exact ω-content λ, not per valuation
   (the specific λ's of a syzygy can be forced expensive).
3. (Inherited-scope note: c_ℓ = 2m·2^ℓ also fails as a law — at
   p = 12 the level-2 block cost Frobenius-collapses to 4, exact
   tables in `s1_syzygy.json`.)

**The per-content assembly instrument** (`s1_omegafloor.json`),
exact and SAT-free, for p = 3·2^a ≤ 12: over all class-nontrivial
ω-syzygies with ≤ smax slots (gap-≤3 pruning by the splitting
argument; free-costs ≥ 1/slot bound the pattern size):
- LB(p) = min Σ_slots free(λ) — a certified lower bound (free(λ) =
  min weight of ANY y-content with ω-component λ);
- UB(p) = min Σ_slots pure(λ) — realized by the pure lift
  (complementary content zero IS a cycle), each constructed and
  re-verified.

| p | LB (free) | UB (pure, realized) | probe status |
|---|---|---|---|
| 3 | 3 | **6** ✓ = banked | exact 6 (calibration) |
| 6 | 6 | **12** ✓ = atlas | exact 12 (calibration) |
| 9 | 3 | **18 — NEW object** | SAT climbing (floor > 11 at kill-check) |
| 12 | 12 | **24 — NEW object** | SAT climbing (floor > 9) |

The UBs equal 2p at every tested period — the pure lift of the
3-slot ω-syzygy IS the 2p species, ∀a ≤ 2, and its construction
scales (the m-dilation multiplies pure costs by the trace-count
law).  The LB side is weak because the free relaxation ignores the
**barren coupling**, and here is the mechanism (hand-verified at
p = 6, the purity lemma candidate):

> **Purity mechanism.** On the 3-slot pattern the barren factors
> admit NO nonzero syzygy (their syzygy modules are free on
> (B̂_η h, Ã_η h), and B̂_η has no monomial multiples — its roots
> are nonzero — so any nonzero barren content needs ≥ 2+2 slots
> with a (B̂,Ã)-shape).  Hence on small patterns every column is
> forced PURE, costs revert to pure(λ), and Σ pure = 2p.  Mixed
> patterns large enough to host barren syzygies pay ≥ #slots and
> the specific joint contents (the p = 6 SAT floor 12 shows they
> never win there).  The quantitative "mixed patterns never beat
> 2p" statement is the named open lemma of this stage.

### §2.2.1 The σ-conditioned h-DP: exact m = 1 floors, SAT-free
### (`a42_s1_jointdp.py`, `s1_jointdp.json`)

The relaxation gap closed.  For p = 3·2^a the CRT
R_p ≅ Λ' × Λ is a bijection on contents, so for ANY cycle v:

> wt(v) = Σ_cols tab[(z'_col, λ_col)]   (identity, not a bound)

where tab is the exact 2^p-element weight table, the Λ'-contents are
(B̂'h, Ã'h) for a unique h (freeness ⟸ unimodularity of the barren
row — Theorem A's lemma), and the λ-contents form a class-nontrivial
ω-syzygy σ.  Fixing σ, min over h is a shortest path (state =
(h_{c−3}, h_{c−2}, h_{c−1}) ∈ Λ'³; structure constants y³ and y+y²),
with EXACT λ=0 tail costs on both sides (Dijkstra to/from the
all-zero state — free-boundary relaxation demonstrably leaks:
it returned 10 < 12 at p = 6; the exact tails return 12).

Results (exhaustive over class-nontrivial σ with ≤ smax slots, gaps
≤ 3 by splitting, up to x-translation):
- p = 3: floor 6 over 4,170 σ (smax 6) — matches banked exactly.
- p = 6: **floor 12 over 55,287 σ (smax 7)** — the atlas value,
  re-derived with no SAT and no BFS.
- p = 12 (a = 2): running; the pure 3-slot σ prices at 24.

Claim tier: certificate (exact linear algebra + shortest paths;
independently re-runnable) with the STATED SCOPE ω-support ≤ smax
slots; the σ-size gap (smax < s < 2p) is closed today only
solver-tier (the CylWindow SAT minimizes over ALL supports at once:
p=6 UNSAT ≤ 11 window-complete).  Closing it certificate-tier needs
a per-slot cost ≥ 2-type lemma for large mixed patterns — open.

**Stage-2 consequence (B12 sketch, contingent on floor(12) = 24):**
by Lemma K + the L1 spanning branch, d((ℓ,12)) ≥ min(24, ⌈ℓ/4⌉) on
the windowed/spanning dichotomy, so d((ℓ,12)) = 24 for all
6 | ℓ ≥ 93 (UB: two stacked L12's, x-local); (12,12) = 18 is the
TOROIDAL exception — its minima are gap-free in both axes (measured,
12/12 witnesses), i.e. two-gross beats the cylinder floor only by
wrapping.  The intermediate ℓ need B6-style per-ℓ descents.  This is
the second row of the B6 ladder; the b = 1 column consumes rows
p = 6r, so the ∀a floor theorem (deformation cascade) is the
gateway to the column.

### §2.2.2 The deformation calculus (the ∀a mechanism, in progress)

Λ_a ≅ F₄[π]/(π^{2^a}) (equal characteristic; ζ = y + π + π² + ...
Artin–Schreier–Teichmüller-style), and the pair becomes

> Ã_ω = (x³+1) + π,   B̂_ω = (1+x+x²) + π(x+x²)(y+1).

Deforming the residue syzygy σ₀ = (1, 1+x) order-by-order in π: the
order-1 obstruction is h₀·[ζ²(x+x³) + 1] ≡ 0 mod (1+x+x²), and the
bracket vanishes at x = ζ² (the TANGENT point) but not at x = ζ (the
transverse point) — so lifting forces (x+ζ) | h₀.  The tangency of
V(A,B) is exactly what lets deformations through on one branch; the
transverse branch blocks them and forces support growth
(supp((x+ζ)^k) = 2^{popcount(k)}).  The support-doubling intuition
lives at EXACT class level (not the level-≥ℓ filtration the S-table
measured — that is falsified-item 1's resolution).  Working out the
full cascade (all orders, both brackets) is the named route to
floor(3·2^a) = 6·2^a ∀a, hence to the b = 1 member rows.

### §2.3 S1 remaining plan

1. Theorem A write-up (done in substance: §2.1 + the machine checks;
   polish owed).
2. The purity lemma (mixed-pattern assembly bound) — the single
   named gap between the instruments and "floor(p) = 2p for all
   p = 3·2^a" at certificate tier; then the m-dilation and the
   127-line margin for full ∀p.
3. p = 9 / p = 12 probes to completion (running; checkpointed,
   RSS-capped after the 2026-08-28 OOM sweep killed the first runs).
4. p = 15 (m > 1 discriminator) once 1-3 settle the m = 1 spine.

---

## §3 Falsified / corrected (running ledger)

1. S-law S_ℓ(a) = 3·2^{a−ℓ}: FALSE (S1, exhaustive; min slots = 3
   at all levels a ≤ 2).
2. Naive assembly min_ℓ S_ℓ·c_ℓ: FALSE (predicts 6 at p=6 vs true
   12); replaced by per-content assembly.
3. c-law c_ℓ = 2m·2^ℓ: FALSE at p = 12 level 2 (Frobenius
   collapse; exact tables banked).
4. §1.3's k = 26 at (762,762) "discovery" framing: RETRACTED as a
   novelty — aggregate-anticipated by 2503.04699 (k_max = 26,
   periods (762,762)); our addition is the per-point/2-adic
   refinement and the b = 1 blocking arithmetic.
5. (cross-lane, conservative direction) 10 banked A40 triage rows
   are frozenset-collision phantoms (§1.5).

## §4 Residue / next

1. Independent verification of k = 26 at (762,762) by a second
   method; engage 2503.04699's gcd-law as the cross-check
   (k(762,762) = k(gcd, gcd) consistency).
2. §2.3 items 2-4 (purity lemma; probes; p = 15).
3. Stage 2 (members): the 2-D ω-bi-block structure on tori; the
   tangency direction vs the member lattice (the b-bit); gates
   d((18,12)) = 24, d((12,12)) = 18 (= two-gross), d((18,6)) = 12.
4. The W7 chirality weight question, now spectrally posed: why does
   the TANGENT F₄ class carry weight 8 at (ℓ,7,ℓ−2) while the
   transverse class floors ≥ 14 — the local-multiplicity-2
   structure as weight-softener (Stage 2 mechanism candidate).
5. Lean targets: the unimodular lemma + Theorem A (decide-shaped);
   the variety table as kernel certificates.
