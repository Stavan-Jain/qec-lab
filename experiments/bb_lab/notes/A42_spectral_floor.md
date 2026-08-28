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

**State of the lane (2026-08-28, session 1):**
- V(A,B) exact: T (transverse) ⊕ G (tangency, mult 2) F₄ orbits +
  W (F₁₂₈, ord 127); spectral k-formula validated on 272 frames +
  members r ≤ 8; k-law: b = 1 column k = 12 ∀r (theorem), b = 0
  jumps to 26 at 127 | r (formula tier; ingredients code-grounded).
- Compact-phase floors: floor(3) = 6, floor(6) = 12,
  **floor(12) = 24** — certificate tier over ω-support ≤ 7 slots
  (h-DP, no SAT), solver-corroborated; barren periods (3∤p, 127∤p)
  EMPTY at all weights (Theorem A, machine-checked p = 5, 7);
  weight-2p objects realized at p = 9, 12 (new; atlas wall was
  p ≤ 8).
- **Theorem H (Tor-purity)**: H(period-p ω-homology) = M[π^{2^a}],
  M = F₄ ⊕ F₄[t]/(t²) = transverse ⊕ tangency scheme structure —
  the class inventory STABILIZES at depth 1 (dim 6 ∀a ≥ 1).
- Cost side: pure(λ) = unique-lift weight; m-scaling exact
  (pure_{3m·2^a} = m·pure_{3·2^a}); UB = 2p ∀ 3 | p.
- Named open: L-pure (mixed patterns), L-band (8..2p−1 slot band) —
  the two quantitative lemmas between here and floor(p) = 2p
  ∀ 3 | p, hence the b = 1 member windowed branch at 12r.

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
missed.  RESOLVED same-day: the A40 lane landed the parity-aware
transport fix + regenerated table (commit 66e0e7b); this lane's G5
diff pass is the detection record.

---

## §2 The Stage-1 program: the ∀p compact-phase floor via the ω-tower

Target (L-P of A40 §9.7): no nontrivial x-compact phase of the
period-p straight cylinder (either lane) has weight < 2p, ∀p ≥ 2.
The spectral decomposition collapses this to a two-parameter family.

### §2.1 Theorem A (barren periods) — the proof

> **Unimodular lemma.** Over any commutative ring R, if
> uA + vB = 1 then ker((f,g) ↦ Af+Bg) = im(h ↦ (Bh, Ah)), i.e. the
> 3-term complex R → R² → R is exact in the middle.  *Proof.* A
> boundary is a cycle (char-free: A(Bh) − B(Ah) = 0 by
> commutativity; in char 2, +).  Conversely, given Af = Bg, set
> h = vf + ug: then Bh = Bvf + Bug = (1−uA)f + uBg = f − u(Af−Bg)
> = f, and Ah = Avf + Aug = vAf + (1−vB)g = g + v(Af−Bg) = g.  ∎
> (It follows that the syzygy module of a unimodular row is FREE on
> (B, A) — consumed by the h-DP and the cost assembly.)

> **Theorem A.** Let p ≥ 1 with 3 ∤ p and 127 ∤ p.  Then the
> period-p compact cylinder homology of the tour-de-gross pair
> vanishes identically — no nontrivial compact phases at any
> weight, both lanes, every 2-part.  *Proof.*  Write p = q·2^a, q
> odd.  y^p − 1 = Π_{d | q} f_d^{2^a} over F₂ with f_d the
> irreducible factors of order-d roots; CRT splits the cylinder
> ring R_p[x^±] = Π_d Λ_d[x^±], Λ_d = F₂[y]/(f_d^{2^a}), and
> homology is the direct sum over factors.  Fix d and a root η of
> f_d, F_Q = F₂(η).  If d = 1: Ã(x, 1) = x³ (a unit times 1), so
> the row is unimodular over the residue and, since the correction
> terms lie in the nilpotent ideal (f_d) of Λ_d, any residue Bezout
> identity uÃ + vB̂ = 1 + n lifts (n nilpotent because f_d^{2^a} =
> 0, so 1 + n is a unit): the row is unimodular over Λ_d[x^±].
> If d > 1: coprimality of Ã(x,η) = x³ + (η+η²) and B̂(x,η) =
> 1 + (η³)(x+x²)·η-units in F_Q[x^±] fails only if they share a
> root ξ ∈ F̄₂; any common root is nonzero (constant terms/leading
> coefficients: B̂'s constant term is 1; Ã's roots are nonzero
> unless η+η² = 0, i.e. d = 1), hence (ξ, η) ∈ V(Ã, B̂) — but the
> variety table (§1.1) has ord β ∈ {3, 127} only, and 3 ∤ p,
> 127 ∤ p exclude ord(η) ∈ {3, 127}.  So the residues are coprime,
> Bezout lifts as above, and the unimodular lemma kills the
> factor's homology.  All factors barren ⟹ H = 0.  The θ′-lane is
> the antipode-automorphic image (its variety's β-orders are the
> α-orders {3, 127} — same barren set).  ∎

Machine checks: dim Z_W = dim B_W on windows at p = 5, 7 (banked,
`s1_cylfloor.json`), consistent with the atlas's trivial-only
verdicts at p ∈ {2,4,5,7,8}.  The old "unbounded chain depth" wall
is dodged: only nilpotency of the 2-part radical is used.  For
3 | p (and 127 | p) the homology concentrates on the
(y²+y+1)^{2^a}-factor (resp. the 127-factor): **the ∀p problem
reduces to the ω-chain family Λ_a plus classical block-weight
arithmetic** — and Theorem H (§2.2.3b) then collapses the depth
ladder too.

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
- p = 12 (a = 2): **floor 24 over 99,399 σ (smax 7, 364 s)** — the
  FIRST exact compact floor beyond the atlas's memory wall (p ≤ 8),
  = 2p on the nose; independent SAT corroboration climbing
  (floor > 12 at last checkpoint) and the realized weight-24 object
  is the matching UB.

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

### §2.2.3 THE STABILIZATION (S1e, `a42_s1_hfiltration.py`)

> **H(a) ≅ F₄[π]/(π²) ⊕ F₄ for a = 1, 2, 3** (dim_F₂ 6; im π = 2;
> π² = 0 on H; window-stable at two widths per depth; a = 0: F₄²).

The cylinder ω-homology STABILIZES at depth 1: the class inventory
is the SAME six dimensions at every period p = 3·2^a, a ≥ 1 (and,
since the ω-factor ring does not see the odd dilation, at every
p = 3m·2^a).  The ∀p floor problem is thereby reshaped: not an
unbounded ladder of new classes, but ONE fixed inventory (a height-2
π-tower plus a height-1 summand) whose per-class minimal realization
costs must scale like 2p.  Claim tier: exhaustive linear algebra at
a ≤ 3, two windows each; the ∀a statement is the natural
quotient-functoriality theorem to prove next.

### §2.2.3b THEOREM H (Tor-purity: the stabilization proven)

Λ_a ≅ F₄[π]/(π^{2^a}) (equal characteristic; Teichmüller ζ = y + π
+ π² + ... solves ζ²+ζ+1 = 0).  Work over the Noetherian 2-dim ring
S = F₄[π][x^±] with the (per-a truncated) pair; Ã = (x³+1) + π
exactly.  Since V(Ã, B̂) ⊂ Spec S is finite and S is CM, (Ã, B̂) is
a regular sequence, so the length-2 Koszul-shaped complex
S → S² → S resolves M := S/(Ã, B̂), and base change to
S/(π^{2^a}) gives

> **H(a) = Tor₁^S(M, S/π^{2^a}) = M[π^{2^a}]** (π-power torsion of
> the intersection module).

Localize M at the π-nilpotent locus (π = x³+1 forces x ∈ {1, ζ,
ζ²}; truncation-artifact roots x₀ of the eliminated ρ(x) =
B̂|_{π=x³+1} have π = g₁(x₀) ≠ 0 acting invertibly — zero torsion):
- x = 1: ρ(1) = ĝ(1) = 1, unit — dies.
- x = ζ (transverse): ρ = t·unit (linear coefficient
  1 + ζ²·ζ² = ζ² ≠ 0) — M_ζ ≅ F₄, π acts as 0.
- x = ζ² (the TANGENT point): the linear terms cancel exactly
  (1 + ζ²·ζ = 0) and the t²-coefficient is 1 + ζ² + ζ² = 1 —
  M_{ζ²} ≅ F₄[t]/(t²), π acts as t·unit.

Hence **H(0) = M[π] = F₄ ⊕ tF₄ ≅ F₄² (dim 4) and H(a) = M ≅
F₄ ⊕ F₄[π]/(π²) (dim 6) for every a ≥ 1** — exactly the measured
filtration (§2.2.3), now a theorem.  The compact-phase class
inventory IS the scheme structure of the ω-slice of V(A,B): the
height-2 π-tower = the tangency; the height-1 summand = the
transverse point.  The same argument at the 127-factor (transverse,
residue field F₁₂₈) gives its H = skyscraper ∀ depth; barren factors
give 0 (Theorem A).  Uniformly:

> **H(period-p cylinder) = ⊕_{variety points P with ord(β_P) | p}
> (local ring of the ω-slice at P)[π_y^{2^{v₂(p)}}]**.

Claim tier: hand proof over the measured a ≤ 3 verification; the
write-up owes the CM/regular-sequence lemma, the localization
bookkeeping, and the truncated-Teichmüller details — no gaps
expected.  Consequence: the ∀p floor needs only the COST half (how
the fixed classes' minimal realizations grow with p) — conjecturally
floor(p) = 2p exactly for 3 | p (measured/DP: 6, 12, 24 at p = 3, 6,
12; SAT: 18 at p = 9 pending), +∞ barren, ≥ 2p + 2/127-margin on
the 127-line.

### §2.2.4 The member synthesis: the conjecture's b-bit is the
### K-dichotomy

2p at p = m = 6r is EXACTLY the conjectured b = 1 distance 12r.
Lemma K splits any member logical into:
- **windowized branch** (some x-gap ≥ 4): unrolls to a compact
  period-m cylinder cycle, same weight, cylinder-nontrivial whenever
  torus-nontrivial (trivializers reduce); so the cylinder floor
  floor(m) = 2m = 12r bounds this branch — the whole b = 1
  conjecture value, from this lane alone if floor(m) = 2m holds ∀m
  of the form 6r;
- **toroidal branch** (gap-free in both axes): where two-gross
  actually lives — its d = 18 = 12r − 6 (r = 2, b = 0) minima are
  measured gap-free (12/12), i.e. b = 0's conjectured 6(2r−1) is the
  TOROIDAL branch's floor, one wrap-discount (−6) below the
  cylinder.

So the family law d = 6(2r + b − 1) = [cylinder floor at b = 1] /
[toroidal floor at b = 0], and the b-bit's entire job is to forbid
the wrap discount — the spectral face of A40 §11's
boundary-coupling wall.  The division of labor is now clean: this
lane owns the cylinder half (stabilized H + cost scaling); the
toroidal half is the wall/momentum program (A40), with the
(12,12)-style gap-free minima as its witnesses.

### §2.2.5 The cost side (S1f): m-scaling exact; the ∀p architecture

The ω-block ≅ Λ_a via the content map (injective — verified), so
pure(λ) is the weight of the UNIQUE complementary-zero lift, and
(`a42_s1_mscaling.py`, exact at (m,a) = (3,0), (5,0), (7,0), (3,1)):

> **pure_{3m·2^a}(λ) = m · pure_{3·2^a}(λ) for every λ** — with the
> base valuation buckets a=0: {2}, a=1: {2,4}, a=2: {2,4,4,8}
> (pure(π^ν·unit) = 2m·2^{popcount(ν)}: the CMSS digit shadow).

Combined with Theorem H (σ-side m-independent): the pure-lift UB is
2p for EVERY p ≡ 0 mod 3, all m, all a ≤ 2 realized (and ∀a once a
pure level-top σ is exhibited per a — it is: π^{2^a−1}σ₀).  The ∀p
floor architecture now reads:

  floor(p) = 2p for 3 | p, 127 ∤ p, provided
  (L-pure)  mixed patterns never beat pure ones [open: the purity
            lemma; certificate-true at m = 1, a ≤ 2 via the h-DP
            with exhausted σ ≤ 7 slots; solver-corroborated at
            p = 9, 12],
  (L-band)  ω-supports in the band (7, 2p) slots never beat 2p
            [open; counting gives only ≥ #slots].

Both open lemmas are quantitative-combinatorial, not structural —
the structural half (classes, costs, scaling) is done.

### §2.3 S1 remaining plan

1. Theorem A + Theorem H write-ups (substance done; polish owed —
   CM/regular-sequence lemma, localization bookkeeping, truncated
   Teichmüller details).
2. The purity lemma (L-pure) and the band lemma (L-band) — the two
   named quantitative gaps to floor(p) = 2p ∀ 3 | p at certificate
   tier; the 127-line margin (≥ 4·64 per 127-period vs 2·127) rides
   the same assembly.
3. p = 9 / p = 12 probes to completion (running; checkpointed,
   RSS-capped after the 2026-08-28 OOM sweep killed the first runs).
4. The B12 row (d((ℓ,12)) = 24 ∀ 6 | ℓ ≥ 93, toroidal exception at
   (12,12)) — write up once the p = 12 corroboration lands; then
   the (ℓ, 6r)-row ladder toward the b = 1 column.

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
