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

**State of the lane (2026-08-28, session 2 — read §2.4–§2.8 first;
two session-1 items corrected there):**
- The b = 1 WINDOWED-BRANCH THEOREM is assembled (§2.8): windowed
  member logicals floor at 12r via floor_cyl(6r); r = 1
  unconditional at 12; **r = 2 at >= 18 certificate (floor(12) >=
  18 unrestricted all-classes, the jet runs)**, = 24 modulo the
  L-band corner; general r modulo (L-pure, L-band).  The b = 1
  problem reduces to the toroidal sector (handed to A40).
- S1h's residue register PROVEN blind at a >= 1 (ker = ker(pi);
  identically zero at a = 2): the S1h p = 12 "floor > 11 racer
  certificate" RETRACTED to solver tier; floor(9) (a = 0)
  re-certified.  Repair SHIPPED: the two-branch full-depth jet
  engine, jointly class-complete, ground-truth-validated on all
  315 atlas cycles (§2.5, §2.5b).
- The named 11+11 sandwich FALSIFIED as a certificate (block
  double-pay with witnesses; the sound variant's frontier shell is
  2^{2p+1} at cost 0) — §2.4.1.
- New machine checks: Theorem A at p = 10, 11, 13, 14; Theorem H
  inventory at p = 15, 21 (dim 4) and p = 18, 24 (dim 6 — the r = 3
  and r = 4 member periods; p = 24 is the first full-window a = 3
  contact); parity columns = 3 at all eight (§2.6).
- L-band pure half at s >= p CLOSED (cost-1 rigidity + pure >= 2);
  the beta-lemma (barren-only columns cost >= 3 at m = 1) pinned
  (§2.7).

**State of the lane (2026-08-28, session 1):**
- V(A,B) exact: T (transverse) ⊕ G (tangency, mult 2) F₄ orbits +
  W (F₁₂₈, ord 127); spectral k-formula validated on 272 frames +
  members r ≤ 8; k-law: b = 1 column k = 12 ∀r (theorem), b = 0
  jumps to 26 at 127 | r (formula tier; ingredients code-grounded).
- Compact-phase floors: floor(3) = 6, floor(6) = 12,
  **floor(9) = 18 (FULL certificate — racer + parity + UB)**,
  **floor(12) = 24** (certificate over ω-support ≤ 7 slots via the
  h-DP; unrestricted > 11 racer-certified, ≥ 14 solver; UB 24
  realized); barren periods (3∤p, 127∤p) EMPTY at all weights
  (Theorem A, machine-checked p = 5, 7); parity lemma (all compact
  cycles even) halves every ladder; the atlas's p ≤ 8 wall is gone
  (racer: p = 9 in 9 s).
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

### §2.1b Parity lemma (probe-halver; coordinator prompt + mechanical
### verification)

> **Lemma (even weights).** On the period-p cylinder (any p), every
> compact cycle of either lane has EVEN weight.  *Proof.*  Each
> qubit variable appears in exactly 3 = |Ā| = |B̄| checks (the three
> support offsets are distinct), so the sum of ALL check rows is the
> all-ones functional; a cycle zeroes every check, hence its weight
> is 0 mod 2.  ∎  Verified mechanically on the exact window systems
> at p ∈ {3, 6, 9, 12} (column weights identically 3); consistent
> with every banked spectrum (atlas {12, 14}, species {8, 10, 12}).

Consequences: weight-stepped probes run even weights only (halves
every ladder), and UNSAT at 2p − 2 yields floor ≥ 2p outright: at
p = 9, floor > 15 (banked) + UNSAT at 16 ⟹ floor = 18 = 2p with
the realized UB — no w = 17 step needed; at p = 12 the ladder is
14, 16, 18, 20, 22.

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
- p = 12 (a = 2): **floor 24 over 99,399 σ (smax 7, 364 s)**
  [S2: scope is smax 7 slots within a 10-column span — ledger
  §3.8] — the
  FIRST exact compact floor beyond the atlas's memory wall (p ≤ 8),
  = 2p on the nose; independent SAT corroboration climbing
  (floor > 12 at last checkpoint) and the realized weight-24 object
  is the matching UB.

Claim tier: certificate (exact linear algebra + shortest paths;
independently re-runnable) with the STATED SCOPE ω-support ≤ smax
slots [S2: within a 10-column span — ledger §3.8]; the σ-size gap (smax < s < 2p) is closed today only
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

### §2.2.5b THE RACER (S1h): floor(9) = 18 certified in 9 seconds
### (`a42_s1_atlas9.py`, `s1_atlas9.json`)

The S4 atlas automaton, re-engineered (coordinator direction): uint64
bit-packed states, Dial-bucket min-cost BFS with zero-cost closure,
np.unique/np.union1d dedup, y-rotation quotient (verified to
commute), NO parent tracking — the state carries a 4-bit ω-class
register pair (F₄ evaluations at the two x-points of the ω-gcd,
per-step x₀⁻¹-rescaled so march-position drops out; v1's cofactor is
a UNIT — B̄_ω ~ ĝ — so the input-block registers detect the syzygy
parameter h(x₀) exactly, both points).  Registers are valid exactly
when 3 | p (the switch-on condition; the p = 5 control caught the
scope error — barren p classify by Theorem A instead, first-class
falsification of the naive "always-on" functional).  Controls: p = 3
min nontrivial 6 ✓; p = 5 returns {6, 10} all trivial-by-Theorem-A ✓;
p = 6 boundary returns at {6, 10} classified TRIVIAL and first
nontrivial exactly 12 ✓ (the decisive classifier control).

> **floor(9) = 18 = 2p, CERTIFICATE tier.**  The p = 9 run completed
> level 16 in 9 s: every compact cycle of weight ≤ 16 enumerated
> with its class value (returns at {6, 10, 12, 14, 16}, ALL
> registers zero), so no nontrivial compact cycle of weight ≤ 16
> exists; the parity lemma excludes 17; the realized weight-18
> object (§2.2) is the matching UB.  Lane AB; the θ′-lane by the
> banked duality (assumption recorded in the JSON).  Ops note: peak
> RSS 3.2 GB for seconds during the final level (the 2 GB check
> fires at level granularity; the abort flag raised AFTER level 16
> completed, so the certificate is unaffected).  This supersedes the
> SAT probe (killed at floor > 15, checkpoint banked) at a stronger
> tier and ~600× the speed; the S5 "memory wall at p = 9" is gone.

### §2.2.6 The class-weight profile (S1g): the complete cost law
### (`a42_s1_classprofile.py`, `s1_classprofile.json`)

Min realization weight PER nonzero H-class (63 classes), by the
h-DP, same σ-scope:

> p = 6:  {12: 45 classes, 14: 12, 18: 6}
> p = 12: {24: 45 classes, 28: 12, 36: 6}
> [S2: per-class minima over the smax-7-slot, 10-column-span scope
> — ledger §3.8]

The strata counts (45, 12, 6) are IDENTICAL (Theorem H's fixed
inventory) and the weights are exactly **{2p, 7p/3, 3p}** at both
periods — the closed-form class-weight law candidate for all
p = 3·2^a, a ≥ 1.  Cross-check: the atlas's p = 6 nontrivial
spectrum {12: 66, 14: 444} confirms the 7p/3-stratum's objects at
14.  The heavy 6-class stratum at 3p and the 12-class stratum at
7p/3 are the levers Stage 2's member assembly can pull: a member
whose logical classes are forced into a heavy stratum floors ABOVE
2p — a candidate spectral mechanism for the b-bit.

**Strata decode (p = 6, full span closure, in-session):** the heavy
3p stratum lies entirely in ker(π) and contains ALL of im(π)
(3 classes) plus 3 further kernel classes (the transverse-summand
candidates); the middle 7p/3 stratum is entirely OUTSIDE ker(π);
the floor stratum = 36 non-kernel + 9 kernel classes.  Economy:
generic classes are cheap; the scheme-distinguished rays —
π-multiples of honest classes and the pure transverse skyscraper —
are the expensive ones.  (Note [πσ₀] is cheap and in ker π but NOT
in im π: σ₀ itself is not a syzygy — the lifting obstruction is
what separates the fresh socle from im(π).)  Identifying the
member-forced classes inside this stratification is the Stage-2
lever.

### §2.3 S1 close-out state and remaining plan

Session-1 final instrument states (all processes terminated clean,
checkpoints banked):
- p = 9: **floor(9) = 18 = 2p, CERTIFICATE tier** (racer level 16
  complete + parity + realized UB; the SAT probe, killed superseded,
  had independently reached floor > 15).
- p = 12: floor = 24 exactly at certificate tier over the ≤ 7-slot
  ω-scope (h-DP, 99,399 σ); unrestricted: racer level 11 complete
  (all returns trivial) ⟹ floor > 11 certificate [S2 CORRECTION:
  the racer's classifier is provably vacuous at a = 2 — this
  sub-claim is retracted to solver tier; see §2.5, ledger §3.7; the
  jet engine's runs supersede it]; SAT ⟹ floor ≥ 14
  solver-tier; UB 24 realized.  The racer's p = 12 wall is TIME
  (×2.7/level ⟹ cap 22 ≈ days) plus an identified memory bug-let
  (bucket insertions un-deduped across inputs — fix known); the
  right tool is a **meet-in-the-middle sandwich** (forward and
  backward level-11 frontiers meet at 22 — both ends already
  computed cheaply) — the named next move, which would close
  floor(12) = 24 unrestricted.

Remaining:
1. Theorem A + Theorem H write-up polish (CM/regular-sequence lemma,
   localization bookkeeping, truncated Teichmüller details).
2. L-pure and L-band — now needed only where the sandwich doesn't
   reach; the 127-line margin rides the same assembly.
3. The sandwich racer (closes p = 12 and plausibly p = 15, 18, 24 —
   the whole small-period ladder at unrestricted certificate tier).
   [S2: FALSIFIED-as-named and retired — §2.4.1, ledger §3.6; the
   jet engine (§2.5b) is the successor instrument.]
4. The B12 row (d((ℓ,12)) = 24 ∀ 6 | ℓ ≥ 93, toroidal exception at
   (12,12)); then the (ℓ, 6r)-row ladder toward the b = 1 column.
5. Decode the heavy strata's H-coordinates fully (transverse vs
   im(π) identification at p = 12; Stage-2 lever).

---

## §2.4 SESSION 2 (2026-08-28): the sandwich autopsy, the register
## audit, the jet repair, and the windowed-branch assembly

Scope: the named p = 12 closer (§2.3 item 3) is dissected and
FALSIFIED as a complete certificate; a deeper audit then finds a
soundness gap in the S1h register itself at 2-adic depth a >= 1,
retracting one banked p = 12 sub-claim; the repaired (class-complete)
jet engine is built, validated on ground truth, and run; new-period
fills machine-check Theorems A and H; L-pure/L-band are sharpened
with one new closed piece; the b = 1 windowed-branch theorem is
assembled with exact conditionality.  All data in `data/a42/s2_*`;
scripts `a42_s2_{sandwich,racer_deep,fills,registers,jet_racer}.py`.

### §2.4.1 The sandwich, falsified-as-named (s2_sandwich_analysis)

Ground work, all asserted mechanically: the REVERSED automaton is the
same generic Automaton class on the x -> x^-1 pair (derived shape
(nf, no, forced_blk, adv_f, adv_o, top_j) = (2, 3, 0, 1, 0, 3));
every atlas cycle at p in {3, 6} (354 cycles) replays through BOTH
directions — fwd consumes v1-columns left-to-right, bwd consumes
v2-columns right-to-left, each march emitting the other block's
columns exactly and returning to zero at accumulated cost = weight —
with the per-cut state contents and the cost formulas
P_F(g) = |v1[<=g]| + |v2[<=g+2]|, Q_B(g) = |v1[>=g-1]| + |v2[>=g]|
asserted at every cut of every cycle.

The join trilemma (the reason "11 + 11 = 22" cannot certify):
1. BOTH-STANDARD marches double-pay the shared 5-column block:
   P_F(g) + Q_B(g) = w + s(g) at EVERY cut (identity, verified).  A
   (11, 11)-join therefore reaches only cycles with a cut where
   s(g) <= 22 - w.  Witnesses: ALL 66 nontrivial weight-12 cycles at
   p = 6 have min_g max(P_F, Q_B) in [10, 12] — excess over w/2 up
   to +6; a (6, 6)-join finds NONE of them.  Minimal nontrivial
   cycles are block-dense (min interior s >= 2, most cuts heavy).
2. EXACT-PARTITION pairing (bwd charged on exit): P_F + Bx = w with
   Bx = Q_B - s; the max-cut argument is airtight (the F-jump columns
   lie in the shared block, so Bx(g*) <= Wcap - CF - 1: caps
   (11, 10) cover all w <= 22), and measured balance on real cycles
   is fine (min-max <= w/2 + 2) — but the exit-charged frontier
   contains the FREE-SUFFIX SHELL: states reachable at exclusive
   cost 0 = 2^{2p+1} exactly (measured complete: 128 / 512 / 2048 at
   p = 3 / 4 / 5; the recent two v2-inputs ride unpaid and the
   third is confined to ker(1 + Y^{-1}), dim 1).  At p = 12 that is
   2^25 = 33.5M states at cost ZERO, before any level growth — the
   deep side is unenumerable.
3. Every mixed charging scheme interpolates between 1 and 2: the
   heavy-jump mass sits in the 5-column seam and must be paid by a
   side (coverage hole) or by nobody (frontier explosion).  The seam
   is the automaton's memory depth — irreducible for exactness.

Verdict: §2.3 item 3 RETIRED (ledger §3.6).  What survives of the
two-sided idea: nothing needed — the deep standard racer plus the
class-register repair below supersedes it.

### §2.5 The register audit: S1h's classifier is blind at a >= 1
### (s2_registers.json; the decisive session-2 finding)

The S1h racer's nontriviality functional R_A = the pair of RESIDUE
evaluations of the input-block omega-content at the two x-points of
the residue gcd (y -> omega kills pi).  Exact computation of R_A on
the OmegaWindow class space (invariance under trivial syzygies
asserted mechanically; class reps for all 63 labels):

- a = 0 (p = 9 scale): ker R_A = 0, both blocks.  The residue
  register is COMPLETE at a = 0 — **floor(9) = 18 STANDS as
  banked** (re-certified end-to-end by the new engine below).
- a = 1 (p = 6): **ker R_A = ker(pi) — 15 of 63 classes — with
  h-DP min weight 12 = the floor itself.**  The racer misses
  floor-weight classes; the S1h p = 6 control passed only because
  R_A-visible classes also realize 12 (the control conclusion stands
  on the atlas census, not on the register).
- a = 2 (p = 12): **R_A = 0 on ALL 63 classes** — and so is any
  functional factoring through content mod pi^2 (including the
  first-draft jet): the Tor-representative of a pi-torsion class of
  order pi^nu is pi^{2^a - nu}-divisible, so at a = 2 every class
  rep is pi^2-divisible and every mod-pi^2 evaluation dies.
  Consequence: S1h's "p = 12 racer level 11, all returns trivial =>
  floor > 11 (certificate)" claimed a classifier that was provably
  vacuous — RETRACTED to solver tier (the independent SAT probe's
  floor >= 14 is unaffected; the h-DP <= 7-slot floor 24 is
  unaffected).  Ledger §3.7.

### §2.5b The jet repair: full-depth branch registers, class-complete
### in two runs (s2_registers.json; engine a42_s2_jet_racer.py)

Theory, then machine verification:

- Over Lambda_a the v1-block boundary multiplier Pbar_omega =
  1 + ybar^{-1} + x^{-3} ybar has exactly 3 unit roots (Teichmueller
  cube roots x 3rd-root rigidity); the two with residue != 1 are the
  usable branch points x*_T, x*_V (residues omega, omega^2).  The
  FULL-DEPTH register R = sum_c x*^{+c} ev(v1_c) in Lambda_a kills
  boundaries exactly (x* is an exact Lambda_a-root) — the +c
  convention is forced: the -c convention evaluates at x*^{-1},
  which is a root on the OTHER branch, and loses one block's
  invariance (found mechanically; a first-draft R_B failed its
  assert exactly there).  Explicit constants at a = 2: Teichmueller
  zeta = y + pi + pi^2, (1+pi)^{1/3} = 1 + pi + pi^2 + pi^3, both
  verified to full depth in the engine's init asserts.
- The tangency is visible as arithmetic: the two curves share their
  branch point mod pi^2 exactly at the tangent branch (computed:
  both deformed roots have t = 1 there) and separate at the
  transverse branch — this is why residue registers exist at all,
  and why full depth is needed beyond it.
- Kernels on H (exact, both a): ker R(x*_tangent) = the 3
  transverse classes (min weight 3p: 18 / 36 at p = 6 / 12);
  ker R(x*_transverse) = 15 classes (min weight 14 / 28);
  **joint kernel = 0 at a = 1 AND a = 2: two runs — one register
  per branch — detect every nonzero class.**  Any two invariant
  x*-twisted evaluations of the same block at the same branch
  differ by a unit, so the engine's registers have these kernels.

Engine: the S2b deep core (insertion-time dedup — the S1h "memory
bug-let" fixed: 4 GB of duplicate bucket flow at p = 12 level 11
collapses to ~0.5 GB; level 16 in 343 s vs level 11 in 151 s for
S1h), states uint64 (5p bits), register a parallel uint8
(dim Lambda_a <= 8), 9-byte void keys for dedup/membership,
y-rotation canon with register co-rotation by ybar^s (commutation
asserted), RSS-capped clean aborts, one branch per process.

Validation stack (all green, `jet_racer.py controls`):
- p=3 union-min nontrivial 6; p=6 union-min 12 (branch minima each).
- p=9 both branches to level 16: no nonzero-register return —
  floor(9) = 18 re-certified by the new engine.
- **GROUND TRUTH: all 315 atlas cycles at p = 6 (66 nontrivial)
  satisfy [embed-nontrivial <=> some branch register nonzero],
  every single cycle** (engine-exact scalar replay).
- 60 random explicit p = 12 boundaries (Pbar t, Qbar t) replay to
  register 0 in both branches (a = 2 engine-exact invariance).
- The p = 12 runs' per-level state counts and return spectra match
  the register-free deep run (s2_racer12_deep) level by level.

Production (p = 12, cap 18, budget 1.3 h/branch, RSS-capped;
after a lexsort rework of the pair-key path — the first flight spent
80% of its time in numpy's VOID_compare, sampled and replaced by
native u64/u8 lexsort + sort-merge membership, ~2.1x faster/level):
- branch T: **level 18 COMPLETE** (7.57M novel, 12.9M seen, 4013 s,
  peak 1.75 GB), returns {6,10,12,14,16,18} all even, NO nonzero
  register;
- branch V: level 17 complete (1519 s), then the time budget stopped
  level 18 (banked per the coordinator's repricing directive —
  per-level time grows ~x3.2, level 19+ is out of day-budget);
  NO nonzero register.

> **floor(12) >= 18 unrestricted, ALL classes — CERTIFICATE tier.**
> Both branches enumerate every compact cycle of weight <= 17 with
> their branch registers; no nonzero register appears; the joint
> kernel is empty (s2_registers), so no nontrivial compact cycle of
> weight <= 17 exists; parity excludes odd weights.  Lane AB; the
> theta'-lane by the banked duality.  Sharpening from T's level 18:
> a weight-18 nontrivial cycle, if any, must lie in the 3 T-blind
> transverse classes (whose h-DP-envelope minimum is 36).  With the
> realized UB, floor(12) in {18, 20, 22, 24}; = 24 modulo the
> corner [{18 (transverse-only), 20, 22}] x [outside the h-DP
> envelope].  Supersedes the retracted S1h sub-claim (§3.7) AND the
> solver-tier >= 14.

### §2.6 The compact-floor table after session 2

| p | floor / H | tier | mechanism |
|---|---|---|---|
| 3 | 6 = 2p | certificate | atlas + DP + jet controls |
| 5, 7 | H = 0 | theorem + machine | Theorem A; window LA |
| 6 | 12 = 2p | certificate | atlas census + DP (s<=7) + jet |
| 9 | 18 = 2p | **certificate (re-certified)** | jet engine both branches, a=0-complete register + parity + UB |
| 10, 11, 13, 14 | H = 0 | theorem + machine (NEW) | Theorem A; dim Z_W = dim B_W at W in {6,9,12}; parity columns = 3 |
| 12 | **>= 18 unrestricted, ALL classes**; = 24 over the (s <= 7, span <= 10, gaps <= 4) envelope (§3.8); UB 24 | certificate | jet runs (T level 18, V level 17) + parity; h-DP; pure lift |
| 15 | dim H = 4 (NEW check) | Theorem H corroborated (a=0, m=5) | window LA; UB 30 = 2p via m-scaling |
| 18 | dim H = 6 (NEW check) | Theorem H corroborated (a=1, m=3) | window LA; UB 36 = 2p; the r = 3 member period |
| 21 | dim H = 4 (NEW check) | Theorem H corroborated (a=0, m=7) | window LA |
| 24 | dim H = 6 (NEW check) | Theorem H corroborated (a=3, m=1 — first full-window a=3 contact) | window LA; the r = 4 member period |

(The register-free deep run independently enumerated returns
{6, 10, 12, 14, 16} through level 16, and the jet runs' per-level
novel counts match it exactly at levels 0..15, +6 states at 16 —
the first fiber separation of the finer 8-bit register.)

[S3 pure-half column: pure_floor(p) = 2p EXACTLY at every
p in {3, 6, 9, 12, 15, 18, 21} — §2.10; the p = 15, 18, 21, 24
rows' "UB 2p" entries now carry certified PURE-sector floors at
15, 18, 21 as well.]

### §2.7 L-pure / L-band after session 2: statements, one new closed
### piece, and the honest obstruction

Fix p = 3m 2^a (3 | p, 127 coprime for now), the CRT frame of §2.2.
Every compact cycle is (sigma, h) with sigma a class-nontrivial
omega-syzygy (slots = columns with nonzero omega-content) and h the
free barren parameter; wt = sum_cols tab[(z'(h), lambda(sigma))].

Machine-pinned cost facts (s1_jointdp tables, re-verified):
- **beta-lemma (m = 1)**: every column with zero omega-content and
  nonzero barren content has weight >= 3 (a weight-2 such column
  needs (y^2+y+1)^{2^a} | y^d + 1, i.e. 3*2^a | d < p — impossible
  at m = 1; at m > 1 the bound drops to 2 and is realized).
- **cost-1 rigidity**: free(lambda) = 1 exactly on the p monomial
  contents lambda = y^j; every weight-1 column carries a NONZERO
  barren monomial content as well.  Hence: an h = 0 (pure) cycle
  has NO cost-1 slots — every slot costs pure(lambda) >= 2.
- pure(lambda) >= 2 for every lambda (nonzero multiples of the
  barren cofactor are never monomials).

L-pure (unchanged statement): over sigma with <= S0 slots, no h
beats the pure lift; certificate-true at (m = 1, a <= 2, S0 = 7)
by DP exhaustion (55,287 / 99,399 sigma at p = 6 / 12).  ForAll-a:
open; the deformation cascade (§2.2.2) remains the named route.

L-band (slots s in [8, 2p-1]) — session-2 split:
- **[NEW, closed] pure half at s >= p**: an h = 0 cycle with s
  slots costs >= 2s >= 2p.  QED by cost-1 rigidity + pure >= 2.
- pure half at s in [8, p-1]: open; DP-enumerable per p in
  principle (needs smax = p - 1; the sigma-pattern count at
  smax >= 8 with gap <= 4 spans is the wall).
  [S3: CLOSED at every p in {3,...,21}, 3 | p, ALL s — the
  omega-quotient pure racer dodges the pattern wall entirely
  (§2.10); the sigma-pattern enumeration is retired as the
  pure-half instrument.]
- mixed half (h != 0): open.  Session-2 structure: cost-1 slots
  force h != 0 with MATCHED monomial barren content at the slot;
  the barren pair (B'h, A'h) of any h != 0 has forced nonzero
  columns at its run edges (f at c0 and c*+2 always, g at c*+3
  always, g at c0 unless h_{c0} is socle) — each costing >= 3
  (beta) when outside the slot pattern.  The taxes are O(1) per
  h-run; the missing piece is a mechanism converting slot-count +
  h-support into the 2p scaling.  The obstruction, precisely: h may
  hide its entire barren support inside a wide slot pattern (the
  linear system (B'h)|_complement = 0 has solutions for s >= ~8),
  so per-slot accounting alone cannot exceed s + O(1).

Corollary of the split: for the p = 12 target, the open corner
after the jet runs is exactly [weight in {JET+2..22} even] x
[s in [8, 23]] x [any h]; the h-DP covers its (s <= 7, span <= 10, gaps <= 4) envelope at
24 (ledger §3.8), the pure band covers h = 0, s >= 12 at >= 24.

### §2.8 THE WINDOWED-BRANCH THEOREM (b = 1 column), with exact
### conditionality — the session-2 assembly

**Theorem W (b = 1 windowed branch).**  Let r >= 1, C_r the member
on Z_{6r+6} x Z_{6r} (b = 1), m = 6r.  Let v be a nontrivial
X-logical of C_r whose support has some cyclic x-gap >= 4 (the
windowed branch of Lemma K).  Then wt(v) >= floor_cyl(m), the
period-m compact-cylinder floor.  Consequently:
- **r = 1: wt(v) >= 12 = 12r, UNCONDITIONAL** (floor_cyl(6) = 12,
  certificate).
- **r = 2: wt(v) >= 18 unconditional at certificate tier** (the
  two-branch jet runs, joint level 17 + parity); **wt(v) >= 24 =
  12r under the L-band corner** ([{18 (transverse classes only, by
  T's level 18), 20, 22}] x [outside the (s <= 7, span <= 10,
  gaps <= 4) h-DP envelope] only — ledger §3.8; UB realized).
- **general r: wt(v) >= 12r = 2m under (L-pure + L-band at p = m)**
  — plus, when 127 | r, the same two lemmas at the 127-factor (the
  W-line skyscraper joins H(m) by Theorem H's uniform form; its
  conjectured margin is >= 2m + the 2/127-margin of §2.2.3b).
Proof: the Lemma-K unroll (an x-gap >= 4 logical lifts x-compactly
to the period-m cylinder at the same weight, cylinder-nontrivial
since trivializers reduce — banked, standard) + floor_cyl(m).  QED

**Reduction handed to the A40 lane:** with Lemma K, the b = 1
lower-bound problem is now EXACTLY: (i) the cylinder floor
floor_cyl(6r) = 12r (this lane; conditional as above, r <= 2
effectively closed), plus (ii) the doubly-spanning (toroidal,
gap-free-both-axes) sector — where the conjectured b-bit mechanism
is that the b = 1 twist forbids the two-gross wrap discount (A40
§11's boundary-coupling wall).  Any b = 1 member counterexample to
d = 12r must be toroidal.  The measured two-gross minima (12/12
gap-free, §2.2) are the b = 0 witnesses that the toroidal sector is
where the -6 discount lives.

## §2.9 SESSION 3 (2026-09-01): the corridor jet racer — sound
## register-reachability pruning (the p = 12 unrestricted closer)

Charter item: the "register early-exit / earlier-terminating
argument" of §4.3, implemented as the conservative form of
register-reachability pruning: a state that cannot CLOSE at all
within the remaining budget a fortiori cannot reach a
nonzero-register closure — and closure-feasibility is exactly
computable from the falsified sandwich's sound remnant.  Script
`a42_s3_corridor_racer.py`; data `s3_*`.

### §2.9.1 The corridor lemma (soundness, exact)

Ingredients, all mechanically asserted in S2a (354/354 atlas cycles,
every cut) and re-asserted in this session's battery:
- the forward state at cut g is the 5-column window
  (v2[g], v2[g+1], v2[g+2], v1[g-1], v1[g]) at forward cost
  P_F(g) = |v1[<=g]| + |v2[<=g+2]|;
- the backward (x -> x^-1) automaton visits, at backward cost
  Q_B(g) = |v1[>=g-1]| + |v2[>=g]|, the SAME five columns in
  reversed field order — the state correspondence phi = packed
  field reversal (asserted per cut);
- the cut identity P_F(g) + Q_B(g) = w + s(g) with s(g) =
  POPCOUNT OF THE STATE (asserted per cut in popcount form).

> **Corridor lemma.**  Let G(t) = the backward racer's min arrival
> cost of backward state t (the true minimum over all closure data;
> rotation-quotiented, complete through level g_max).  Every compact
> cycle of weight w <= W through forward state s at forward cost c
> satisfies  w = c + Q_B - s(g) >= c + G(phi(s)) - popcount(s).
> Hence a state may be DROPPED whenever
>     c + G(phi(s)) - popcount(s) > W,
> and, when phi(s) is absent from a table complete to g_max
> (i.e. G > g_max), whenever W - c + popcount(s) <= g_max.
> No state lying on ANY compact cycle of weight <= W is ever
> dropped; the splice induction (kept representatives lie on spliced
> cycles of no larger weight, which are themselves corridor-immune)
> shows the pruned jet run's returns <= W, with their branch
> register values, are IDENTICAL to the unpruned run's.  The
> certificate semantics of §2.5b are preserved verbatim.  QED

This is a PRUNE, not a join: the S2 trilemma (§2.4.1 — double-pay,
zero-cost shell) does not apply.  The double-pay of the shared block
is exactly compensated by the +popcount(s) slack, and the backward
racer uses standard charging (no exit-charged shell).

### §2.9.2 The battery (ALL GREEN, `s3_battery.log`)

1. validate_banked PASS.
2. BackRacer step/canon vs the generic Automaton on the x^-1 pair
   at p in {3, 6, 9, 12} (derived shape (2, 3, 0 | adv 1, 0 |
   top_j 3); 500 random state x input steps each + canon
   commutation).
3. phi correspondence + popcount cut identity + G <= Q_B asserted at
   630 cuts (all 39 p = 3 atlas cycles) and 5,829 cuts (all 315
   p = 6 cycles) against full-depth backward tables; absent states
   verified absent only where Q_B > g_max (the contrapositive
   completeness check).
4. 60 random explicit p = 12 boundaries replay through the backward
   automaton (forced columns match, cost = weight, return to zero).
5. Pruned-vs-unpruned RETURN EQUALITY, weights AND register values:
   p = 6 cap 13 both branches (the nontrivial register sets at
   w = 12 reproduced exactly: T {0,7,9,14}, V {0,1,4,6,9,10,12,14});
   p = 9 cap 16 both branches at W in {16, 18} (162,743 states
   pruned at W = 18, zero return drift) — floor(9)'s certificate
   reproduced THROUGH the pruned engine end-to-end.

Engine deltas vs the S2 jet: membership by searchsorted against the
lex-sorted seen store (the stock member() re-merge-sorted the full
12.9M-pair store once per 1K chunk — the dominant cost of the S2
level-18 flight); per-level prune stats; per-level counts of
distinct omega- and barren-window projections of the novel frontier
(the quotient-growth instrument for the S3 hypothesis).

## §2.10 The omega-quotient pure racer: THE PURE HALF FALLS AT
## EVERY PERIOD <= 21 (session 3, `a42_s3b_omega_racer.py`)

The register-quotient racer in exact-cost form on the h = 0 sector.
At p = 3m·2^a the pure-lift cycles (all non-omega content zero) are
in bijection with compact omega-syzygies over Lambda_a, with true
weight sum_cols pure(lambda_col); so the S4 automaton rebuilt with
columns IN Lambda_a (5·dim <= 40 bits of state), per-column cost
pure(lambda) (CRT-idempotent table; pure >= 2 asserted; rotation-
invariance asserted), and the S2e branch registers acting on the
columns themselves, enumerates ALL pure cycles by weight over ALL
supports — no slot/span/gap scope at all.  Controls: p = 3 pure min
nontrivial 6; p = 6 all 21 PURE atlas cycles replay column-by-column
(cost = weight, omega registers == the jet replay registers, class
verdicts match — the entire p = 6 pure sector sits at weight 12);
racer returns == pure-atlas spectrum.  m-scaling checked TABLE-WIDE
at (m, a) = (3, 1): pure_18 = 3·pure_6 on all 16 contents.

> **Pure-half floor theorem (machine certificate, both branches,
> joint kernel empty — s2_registers).  For every p in {3, 6, 9, 12,
> 15, 18, 21}: no nontrivial PURE compact cycle of weight < 2p
> exists, and the racer itself returns nonzero registers at exactly
> 2p — pure_floor(p) = 2p EXACTLY.**  Caps run: 8/13/20/26/32/42/44.
> Consequences: the L-band pure half at s in [8, p-1] (§2.7, open)
> is CLOSED at these periods, s >= p re-derived, and p = 18 — the
> r = 3 member period, where the full racer's 90-bit state is
> infeasible — gets its first quantitative floor piece (the omega
> state is 20 bits; 220 states at level 42).  p = 24 (a = 3) is
> blocked only by u64 packing (5·16 = 80 bits) — a 2-word state
> engine reaches it.

Quotient growth (the S3 hypothesis measured): the p = 12 pure
frontier at level 16 is 1,300 states vs the full racer's 2.17M
(~1,700x compression); growth ~2.8x per even level; pure levels are
EVEN-ONLY (every pure column cost is even — the parity lemma's pure
face).  The p = 12 pure return spectrum: trivial returns at
{12, 20, 22}, nontrivial from 24 (and 26); the p = 18 spectrum is
{36, 42} nontrivial only — the pure sector rigidifies as m grows
(at a = 0, m > 1: every pure column costs exactly 2m, quantizing
the ladder to multiples of 2m; 3 slots force 6m = 2p on the nose).

## §2.11 The halving lemma + the universal sigma-inequality: the
## mixed half analytically, up to a factor 2 (session 3)

Write e = e_p for the CRT idempotent of the omega-factor at
p = 3m·2^a (e ≡ 1 mod (y²+y+1)^{2^a}, ≡ 0 on the complement), and
for a column with barren content z' and omega-content λ let
tab(z', λ) = wt(CRT(z', λ)) as in §2.2.1.

> **Halving lemma (m = 1, any a).**  tab(z', λ) >= pure(λ)/2 for
> every z'.  *Proof.*  e·CRT(z', λ) = CRT(0, λ) (e kills the barren
> part and fixes the omega part), and multiplication by e spreads
> each term of a vector to at most wt(e) terms (convolution bound),
> so pure(λ) = wt(CRT(0, λ)) <= wt(e)·tab(z', λ).  At m = 1,
> wt(e) = 2 for every a: e is the weight-2 idempotent
> y^{p/3} + y^{2p/3} (idempotency is one line; machine-checked at
> p = 3, 6, 12, 24, 48, i.e. a <= 4).  ∎
> Table-checked at p = 12: min over z' of tab equals pure(λ)/2
> EXACTLY in all four buckets ({2,4,6,8} -> {1,2,3,4}) — the bound
> is tight, and finer than the reverse triangle
> tab >= |pure(λ) − pure'(z')| in the 4- and 8-buckets.
> (General m: wt(e) = 2m by m-scaling, so the same argument gives
> tab >= pure/2m — degrading with m; the odd part is semisimple and
> should be handled by its own CRT instead.)

> **Universal sigma-inequality (the pure racer's theorem, restated).
> For every class-nontrivial compact omega-syzygy sigma at
> p in {3, 6, 9, 12, 15, 18, 21}: sum_slots pure(λ_c) >= 2p** —
> no slot-count, span, or gap scope (§2.10).

> **Corollary (the mixed half, up to a factor 2; m = 1).**  Every
> class-nontrivial compact cycle (sigma, h) at p = 3·2^a
> (a <= 2 certified) has
>     wt = sum_c tab(z'_c, λ_c) >= sum_slots pure(λ_c)/2 >= p,
> ANY h, ANY support.  With the parity lemma: >= p + (p even? 0 : 1),
> i.e. >= 12 at p = 12 — analytic, enumeration-free given §2.10.

The remaining factor 2 is exactly the hiding-mass fight, now in a
sharper form than §2.7's obstruction paragraph.  Reverse-triangle
per column gives  wt >= sum_S pure(λ) − D(h, S)  with the hiding
mass  D = sum_{S∩H} pure'(z') − sum_{H\S} pure'(z')  (H = barren-
active columns; pure' buckets at p = 12: {3, 6, 9, 12} — all
multiples of the beta-bound 3).  So L-band mixed reduces to: *the
barren syzygy (B'h, A'h) cannot concentrate more than the sigma-side
pure-slack of its pure'-mass on the slots.*  The halving lemma
bounds the per-column concentration at exactly half; columns where
halving is TIGHT have their barren content FORCED by the omega
content (the halving partner is determined), so a near-p cycle
overdetermines h — the rigidity direction to push next.  At p = 12
the corner is being closed by exhaustion regardless (§2.9); the
analytic route matters for forall-a: the deformation cascade's
target (the pure sigma-floor at all a) now yields the MIXED half to
within the factor 2 for free.

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
6. (S2) The named meet-in-the-middle sandwich ("fwd 11 + bwd 11,
   join on states") is NOT a complete certificate: both-standard
   joins double-pay the 5-column block (all 66 nontrivial w = 12
   cycles at p = 6 have min-max balance >= 10 > 6 — witnesses), and
   the coverage-complete exact-partition variant has a 2^{2p+1}
   zero-cost frontier shell (measured exactly at p <= 5).  §2.4.1;
   `s2_sandwich_analysis.json`.
7. (S2, retraction) The S1h racer's residue register is blind to
   ker(pi) on H at depth a >= 1 (15/63 classes at a = 1, including
   floor-weight classes; identically zero on ALL classes at a = 2 —
   Tor-reps are pi^2-divisible).  S1h's "p = 12 racer floor > 11
   (certificate)" is RETRACTED to solver tier (SAT >= 14 stands;
   h-DP scoped 24 stands; floor(9) at a = 0 is UNAFFECTED — kernel
   empty, machine-checked, and re-certified by the jet engine).
   Repair: the two-branch full-depth jet registers, jointly
   class-complete (§2.5, §2.5b; `s2_registers.json`).
8b. (S3, hypothesis verdict) The session-3 charter hypothesis "the
   quotient's growth is subexponential in the omega-direction" is
   FALSE AS STATED for the mixed frontier: the corridor runs'
   instrumentation shows the omega-projection nearly INJECTIVE on
   the frontier (n_omega/novel = 0.84..1.0 at p = 12 levels 8..19 —
   the exponential growth lives in the omega-direction, 2.9M
   distinct omega-windows at level 19), while the BARREN projection
   saturates at ~6.3K windows and then falls.  The compact direction
   is barren, not omega.  What IS true: the PURE sector's
   omega-frontier is tiny (1,300 states at level 16 vs the full
   racer's 2.17M), and the mixed frontier's omega-territory at
   level c matches the pure racer's at level ~2c - 10 — the halving
   lemma's factor-2 discount, measured (§2.11).  Engine moral: the
   profitable quotient is over the barren direction; the omega
   direction must be enumerated (or theory-collapsed via Theorem H
   classes, which windows do NOT factor through).
8. (S2, scope correction) The S1 sigma-enumerations (s1_jointdp,
   s1_omegafloor, s1_classprofile) confined patterns to a 10-column
   window (Wx = 10), but the splitting argument bounds only
   consecutive-slot gaps (<= 3), not the span: a minimal 7-slot
   sigma can span up to 1 + 6*3 = 19 columns.  Moreover the
   enumerations also SKIP patterns with a consecutive-column gap
   > 4 even inside the window, and the split justification for that
   prune is incomplete for the (sigma, h) floor: a split sigma's
   halves are separate syzygies, but a mixed h can bridge the gap,
   and tab(z', lambda) can undercut tab(z', 0) by up to 2 per
   column, so W(sigma, h) >= DP(sigma_L) does NOT follow columnwise.
   Every banked "omega-support <= smax slots" certificate scope is
   therefore really "<= smax slots, span <= 10, consecutive-column
   gaps <= 4" — §2.2.1's floor-24, §2.2.6's class profile, and
   §2.2.5's L-pure instance all carry the qualifier (annotated in
   place).  The racer/jet certificates are scope-clean (they
   enumerate CYCLES by weight over all supports) — the repair lane
   for the corner is jet depth + L-band theory, not a wider
   enumeration.

## §4 Residue / next

1. Independent verification of k = 26 at (762,762) by a second
   method; engage 2503.04699's gcd-law as the cross-check
   (k(762,762) = k(gcd, gcd) consistency).
2. L-pure forall-a (deformation cascade) and the L-band mixed half.
   [S3 state: the pure-half window is CLOSED at every 3 | p <= 21 by
   the omega racer (§2.10) — the smarter-enumeration item is
   retired; the mixed half is reduced to the hiding-mass factor 2
   (§2.11), and the cascade's forall-a payoff is doubled: pure
   sigma-floor(3·2^a) = 2p would now yield mixed >= p for free.
   The concrete next rungs: (i) the p = 24 pure racer (2-word
   states; register soundness at a = 3 BANKED — s3d joint kernel
   empty), which with halving gives floor_cyl(24) >= 24 toward
   r = 4; (ii) the hiding-mass rigidity argument (halving-tight
   columns determine their barren content).]
3. The p = 12 remaining corner: weights {18 (transverse-only),
   20, 22} x outside the h-DP envelope.  Per-level jet time grows
   ~x3.2 (level 18 = 44 min): levels 19-22 cost ~2.3/7.5/24/77 h —
   beyond day-budget (coordinator repricing, 2026-08-28); the
   certificate path forward is an earlier-terminating argument
   (register early-exit / omega-structure routing / L-band).
4. Stage 2 (members): the 2-D omega-bi-block structure on tori; the
   tangency direction vs the member lattice (the b-bit); gates
   d((18,12)) = 24, d((12,12)) = 18 (= two-gross), d((18,6)) = 12;
   the toroidal-sector reduction of §2.8 is the interface to A40.
5. The W7 chirality weight question, now spectrally posed: why does
   the TANGENT F4 class carry weight 8 at (l,7,l-2) while the
   transverse class floors >= 14 — the local-multiplicity-2
   structure as weight-softener (Stage 2 mechanism candidate).
6. Lean targets: the unimodular lemma + Theorem A (decide-shaped);
   the variety table as kernel certificates; NEW: the register
   soundness pair (invariance + joint-kernel emptiness) is
   finite-linear-algebra shaped — a natural KernelCert candidate.
7. The twisted-lattice extension of Theorem H (TC63 at <(3,6)>):
   does Tor-purity cover twisted cylinders?  (Untouched in S2.)
8. A40 §1 family-sourcing re-verification (flagged in §1.3):
   still owed by the A40 lane.
