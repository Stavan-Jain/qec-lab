# A21 — analytic base floor for f2a6f17e: `LogicalFloor 8`

> **Target.** Analytic (no exhaustive enumeration, no SAT) proof that
> `bb_neigh_z5z15_f2a6f17e` — G = Z₅×Z₁₅, A = 1 + y + x,
> B = xy⁶ + xy¹⁰ + x²y¹², [[150,8,8]] — has no non-boundary 1-cycle of
> weight ≤ 7.  Exact Lean Prop: `coverData.LogicalFloor 8` in
> `QEC/Stabilizer/Codes/BivariateBicycle/Z5Z15F2A6/` (branch
> `claude/a15-m-kernel-route`, PR #61).  Ground truth d = 8 exact
> (witness kernel-checked; ≥ 8 half currently CaDiCaL UNSAT@7).
> Discharges one of the three hypotheses of
> `cover300_pauli_distance_eq_16_of_classification`.

## 0. Conventions (pinned against the Lean target, 2026-07-22)

Repo (`BBChainComplex.lean`): C₂ = C₀ = G, C₁ = G × {L,R};
∂₂ f = (A⋆f | B⋆f) (LEFT block gets A); ∂₁(u_L,u_R) = B⋆u_L + A⋆u_R.
**Cycle condition: B⋆u_L = A⋆u_R.**  Boundaries: u_L = A⋆f, u_R = B⋆f.
Generator column ∂₂δ_g = (A+g, B+g), weight 6.
`chainWeight` = |supp u_L| + |supp u_R|; split (a,b) := (|u_L|,|u_R|).
Lab H_Z-kernel vectors ↦ repo cycles by global reflection g ↦ −g,
blocks unchanged (`Defs.lean` header).  All math below is in the repo
convention unless tagged lab.

## 1. Status of the target vs. existing machinery (session 1 triage)

* **PAR**: ε(A) = ε(B) = 1, so cycles have |u_L| ≡ |u_R| (mod 2) —
  every odd total weight (1,3,5,7) is dead for free.  The target
  reduces to: (α) no nonzero cycle of weight ∈ {2,4}; (β) no
  *non-boundary* cycle of weight 6.
* **A16 class theorem does NOT apply**: certifier run 2026-07-22 —
  W3 ✓, D1 ✓, D2 ✓ (coord-disjoint (False, True)), FRM ✓ (semisimple),
  **(iii) ✗**: A = 1+y+x collapses to a monomial under BOTH parity
  projections (π_x A = x, π_y A = y; maA = [0,1], B mono-in-x only).
  So even the d ≥ 6 layer needs a per-instance argument here.  The
  A16 *mechanisms* (Sidon/difference-multiset counting, projections,
  triangle-image taxonomy §5, σ-match analysis §6) remain the toolbox.
* **T5 of the A15 plan** (past-2w depth: the seven weight-6 splits)
  was never executed for any code; A21 is the first depth rung.
* Weight-6 splits after PAR: (0,6),(6,0),(1,5),(5,1),(2,4),(4,2),(3,3).
  (3,3) contains the 75 generator columns — the statement at weight 6
  is a **boundary classification**, not a kill.
* Weight-≤5 splits after PAR: (0,2),(2,0),(0,4),(4,0),(1,1),(1,3),
  (3,1),(2,2).
* Semisimple frame (|G| = 75 odd): Ann(A) = I(V_A) exactly; every
  convolution equation P⋆u = c is solvable iff ĉ vanishes on V_P, and
  then the solution set is one coset of Ann(P) — |Ann(P)| = 2^{dim}.
  This makes every split's solution set a *small affine space*, which
  is the instance's structural gift (bb_90 Entry-4 style).

## 2. Approach log

### 2026-07-22 session 1 — setup, numeric census (falsify-first)

Plan: (1) pin conventions [done, §0]; (2) certifier triage [done, §1];
(3) numeric session: k/V_A/V_B/component table, μ(Ann A), μ(Ann B),
exhaustive weight-≤6 cycle census by split with boundary tests,
family classification; (4) derive the lemma list; (5) attack.

**Census results (`a21_numeric_session1.py`, `a21_numeric_probes2.py`):**

1. k = 8 ✓; V_A = V_B = **one shared Frobenius orbit** (rep (1,14),
   size 4, F₁₆ component); Â, B̂ vanish nowhere else.  dim Ann(A) =
   dim Ann(B) = 4.
2. **μ(Ann A) = μ(Ann B) = 40**, and ALL 15 nonzero annihilator
   elements have weight exactly 40.  One-sided splits dead to
   weight 39.
3. **Exhaustive weight-≤6 cycle census: exactly the 75 generator
   columns ∂₂δ_g = (A+g, B+g); nothing else.**  0 non-boundary
   cycles of weight ≤ 6.  All 67,525 weight-3 u_L are solvable for
   u_R (single-orbit vanishing ⟹ everything solvable) but only the
   75 A-translates admit weight-3 u_R.
4. Coset A⋆w = B: min weight 19 (profile 19,31,35,…).

**Structural discoveries (all machine-verified in probes2):**

* **(ψ) The mirror symmetry.**  ψ(gx,gy) := (gx, 6gx+4gy) ∈ Aut(G),
  and B = xy⁶·ψ(A), ψ⁻¹(B) = x·A, ψ(dA) = dB.  The induced map
  Φ(u_L,u_R) = (ψ⁻¹u_R, x²y⁶·ψ⁻¹u_L) sends cycles to cycles with the
  split swapped (a,b)→(b,a), and Φ(∂₂f) = ∂₂(x·ψ⁻¹f) — boundaries to
  boundaries, weights preserved.  **Every mirror split is free**; the
  code is a "twisted self-pairing" of the single trinomial 1+x+y.
* **(Lemma W) Weight formula.**  For Sidon P, |P| = 3:
  wt(P⋆T) = 3|T| − 2·E_P(T) + 4·N_P(T), where E_P = #unordered pairs
  of T with difference in dP, N_P = #cells c with c−P ⊆ T.
  (0/8000 random violations; proof = incidence count, Σ_c C(mult_c,2)
  = E_P by Sidon, mult ≤ 3.)
* **dA-cluster tallies** (dA-graph = triangular lattice on the 5×15
  torus): 4-sets max E = 5; 5-sets max E = 7; and over ALL 5-sets
  **E − 2N ≤ 5** (never 6; (6,0) does not occur — 6-edge 5-clusters
  always contain a stencil triangle).  Hence wt(A⋆T) ≥ 5 for every
  |T| = 5 ⟹ **(1,5) dead** (image would need weight 3).
* **(2,4)**: 404 anchored 4-sets have the right image weight
  (E−2N ∈ {3,4}), but **0/404 images lie in the σ_B family**
  (all B⊕(B+δ) translates).  Analytic route: D2-pigeonhole for the
  disconnected (3,0) case (A-translate can't fit in σ_B: ≥2 cells in
  one B-copy ⟹ dA ∩ dB ≠ ∅); connected cases by the spread bound
  (image cells within diam(T)+2 ≤ 5 in the dA-word metric, but a full
  B-translate inside the image forces a length-7 pair (1,6)-class);
  size-4 images: residual finite table (δ ∈ ±(1,6)-class only).
* **(2,2)**: size-6 σ-match dies by D2-pigeonhole instantly (A-copy's
  3 cells into 2 B-copies).  Size-4: 0/36 translate matches
  (δ_L ∈ dB, δ_R ∈ dA); analytic = finite σ-pair-value table.
* **(3,3) classes of u_L by (E_B, N_B)**: (0,0) 2083 / (1,0) 585 /
  (2,0) 27 / (3,0) 3 / (3,1) 3 anchored; **only (0,0) admits
  weight-3 u_R, and only at the generators**.
  - (3,0): u_L = ±A-translate…  3-set classification: dA-triangles =
    {A+t} ∪ {−A+t} exactly (no O3: dA has no order-3 elements; no
    PROG: dA ≠ ±{c,2c,3c}); same for dB.
  - (3,1) [u_L = −B+t]: image = {0}∪dB+t (wt 7); A-side (3,1) gives
    {0}∪dA+t′; equality of symmetric-set translates forces 2s = 0 ⟹
    s = 0 ⟹ dA = dB, D2-dead.  Cross class (1,0) needs a small table.
  - (1,0)/(2,0): both-sides weight-7/5 images; pigeonhole +
    finite σ/cherry tables (to be written).
  - **(0,0) crux = sumset rigidity**: B⊕u_L = A⊕u_R (both direct,
    |c| = 9).  KEY: B+t ⊆ ⊔_r(A+r) with ≤1 cell per A-copy (D2) ⟹
    exactly one ⟹ **mutual transversality**: ∀(t,r): t−r ∈ A⊖B
    (|A⊖B| = 9, distinct by D2).  Localization: fixing r₀ = 0 ∈ u_R:
    **u_L ⊆ A⊖B** (and u_R ⊆ t+B⊖A ∀t) — the whole case reduces to a
    surveyable table over 3-subsets of the fixed 9-set A⊖B.
* **(1,3)**: A⋆T = B+g needs (E,N) = (3,0) ⟹ T = A+t ⟹ image =
  A²+t; 2dA ∩ dB = ∅ ⟹ never a B-translate.  DEAD.  ((3,1) gives
  weight 7.)  Mirror (3,1) free by ψ.
* **One-sided μ = 40, analytic**: Ann(A) ≅ ker ε ⊂ F₂[Z₅] via the row
  recurrence z_{v−1} = (1+x)z_v (rows of A⋆z = 0), closure
  (1+x)¹⁵ = 1 + ring; ker ε ≅ F₁₆ (Φ₅ irreducible), 1+x primitive ⟹
  single multiplicative orbit ⟹ all nonzero z are translates/orbit
  images with equal weight; weight = Σ_{j=0}^{14}|(1+x)^j p| = 40.
  Equivalent global view: Ann = minimal ideal of the shared orbit;
  χ surjective onto μ₁₅ = F₁₆ˣ ⟹ all 15 nonzero elements are
  monomial translates of the idempotent e_O; wt(e_O) =
  5·#{u ∈ F₁₆ˣ : Tr(u) = 1} = 5·8 = 40.

**Proof skeleton (full split map at ≤ 7).**  PAR kills odd totals;
one-sided by μ = 40; (1,1) by D2; (1,3)+(3,1), (1,5)+(5,1) as above;
(2,2), (2,4)+(4,2) as above; (3,3) classification ⟹ weight-6 cycles
= generator columns = boundaries.  Remaining substantive work:
(a) the 5-set inequality E−2N ≤ 5 at analytic grade (or surveyable
reduction; currently 926-cluster machine census),
(b) the (0,0) sumset-rigidity finite table + write-up,
(c) the small residual tables for (2,2)size-4 / (2,4)size-4 /
(3,3)(1,0),(2,0) classes.

**Session-1 late developments (probes 3–4)** — all three residues
above DISSOLVED into clean arguments:
(a) replaced by the **y-span lemma** (row-block analysis; kills
(1,1),(1,3),(1,5) uniformly, no cluster census at all);
(b) replaced by the **overlap lemma** (d(A⊖B) mult-≥3 values = dA∪dB
exactly) + triangle classification + three row-structure certificates;
(c) replaced by T1/T2/T3 pigeonholes + the word-metric spread bounds.
CAUTION filed: hand word-lengths were wrong (|(1,6)| = 6 not 7 —
BFS-corrected in `a21_numeric_probes4.py`; spread kills re-verified
with the correct values, margins hold).
The full proof is written out in §3 below; the math for
`LogicalFloor 8` is **complete at analytic grade** as of session 1.

---

## 3. The proof (analytic, complete) — d(f2a6f17e) ≥ 8 at chain level

Throughout: G = Z₅×Z₁₅ written additively, e_x = (1,0), e_y = (0,1);
A = {0, e_y, e_x}, B = {(1,6), (1,10), (2,12)}; repo conventions (§0):
a 1-chain is (u_L, u_R), the cycle condition is B⋆u_L = A⋆u_R,
boundaries are ∂₂f = (A⋆f, B⋆f).  Write |·| for Hamming weight /
support size, d(S) for the difference (multi)set of a support S,
σ_P(δ) := P ⊕ (P+δ).

### 3.0 Instance constants (all finitely checkable, `probes2–4`)

C0.1  dA = {±(0,1), ±(1,0), ±(1,14)}, dB = {±(0,4), ±(1,6), ±(1,2)};
      both Sidon (D1); dA ∩ dB = ∅ (D2); no 2-torsion, no periods
      (|G| odd).
C0.2  Element orders: dA-orders {15,5,15}, dB-orders {15,5,15}: no
      order-3 elements ⟹ **no O3 triangle class** for either set;
      neither dA nor dB is ±{c,2c,3c} ⟹ **no PROG class**; for
      δ ∈ dA: 2δ ∉ dA, for β ∈ dB: 2β ∉ dB ⟹ **no 3-AP** has its
      difference pair inside either set.  Hence (A16 Lemma 5.3.1
      taxonomy, instance form): the 3-sets T with d(T) ⊆ dA are
      exactly {A+t} ∪ {−A+t}, and with d(T) ⊆ dB exactly
      {B+t} ∪ {−B+t}; N_A(A+t) = 0, N_A(−A+t) = 1 (mirror for B).
C0.3  2dA = {±(0,2),±(2,0),±(2,13)} is disjoint from dB;
      2dB = {±(0,8),±(2,12),±(2,4)} is disjoint from dA.
C0.4  **ψ-mirror**: ψ(gx,gy) := (gx, 6gx+4gy) is an involutive
      automorphism of G with ψ(dA) = dB, ψ(A)·xy⁶ = B, ψ(B) = x·A.
      Φ(u_L,u_R) := (ψu_R, x²y⁶·ψu_L) maps cycles to cycles with
      blocks (and the two stencil geometries) swapped, and
      Φ(∂₂f) = ∂₂(x·ψf).  All weights preserved.
C0.5  **Word metric**: on the Cayley graph Cay(G, dA) (a 5×15
      triangular-lattice torus), the dB elements have lengths
      {(0,4): 4, (1,2): 3, (1,6): 6} (BFS-verified); every full
      B-translate contains a cell pair at distance 6; every σ_B(β),
      β ∈ dB, contains a cell pair at distance 7 (6-row table, C4).
C0.6  **σ pair tables** (T1): for every β ∈ dB, the 6 pairwise
      difference values of the 4-set σ_B(β) contain NO dA element
      and all six dB elements (3 dB-pairs); mirror: σ_A(α), α ∈ dA,
      contains NO dB-pair.
C0.7  **S-set table** (T2): S_A := {0}∪dA has no dB-pair;
      S_B := {0}∪dB has no dA-pair.
C0.8  **Cherry table** (T3): for every A-cherry T = {0, α₁, α₁+α₂}
      (α_i ∈ dA, α₁+α₂ ∉ dA ∪ {0}), the 5-cell image A⋆T contains
      ≤ 1 dB-pair (18 shapes: twelve 0s, six 1s); mirror for
      B-cherries and dA-pairs.
C0.9  **Overlap lemma** (T5): A⊖B has 9 distinct cells (⟸ D2), and
      the multiplicity-≥3 values of the 72-element multiset d(A⊖B)
      are EXACTLY dA ∪ dB (profile: 28 ones, 4 twos, 12 threes).
C0.10 **Row structure**: B's nonzero rows (y-levels) are {6,10,12}:
      cyclic gap word (4,2,9), min-arc span 6, min gap 2.  The three
      9-sets A⊕B, B⊖A, A⊖B each have three 2-cell rows and three
      1-cell rows; their (doubled-row gap word, singleton offset)
      invariants are ((4,2,9), +1), ((4,2,9), −1), ((2,4,9), +1)
      respectively — pairwise non-translate (translation preserves
      the gap word up to rotation and the singleton offset;
      {D+1} = {D−1} would force a period 2 of D = {6,10,12} in Z₁₅,
      impossible for a 3-set since ord(2) = 4).
C0.11 A⊕B is aperiodic (an order-3 period would shift its row set
      {6,7,10,11,12,13} by 5; it doesn't match).
C0.12 ker(1+x) in F₂[Z₅] = {0, ring}, ring := 1+x+x²+x³+x⁴.

### 3.1 Parity (PAR)

ε: F₂[G] → F₂ is a ring hom; ε(A) = ε(B) = 1.  A cycle satisfies
ε(u_L) = ε(u_R), so |u_L| ≡ |u_R| (mod 2): **no cycle has odd total
weight**.  This kills weights 1, 3, 5, 7.  ∎

### 3.2 One-sided floor: μ(Ann A) = μ(Ann B) = 40

Â and B̂ vanish exactly on the single Frobenius orbit O of the
character χ = χ_{(1,14)} (F₁₆ component; k = 2·4 = 8 ✓).  Hence
Ann(A) = Ann(B) = M_O, the minimal ideal of spectral support O,
dim 4.  χ: G → F₁₆ˣ is onto (its restriction to y already covers
μ₁₅ = F₁₆ˣ), so translations act transitively on M_O ∖ {0}: all 15
nonzero elements share one weight.  That weight is
wt(g ↦ Tr_{F₁₆/F₂}(χ(g))) = |ker χ| · #{u ∈ F₁₆ˣ : Tr u = 1}
= 5·8 = 40.  (Equivalent row view: A⋆z = 0 ⟺ z_{v−1} = (1+x)z_v with
ε_x(z_v) = 0; 1+x is a primitive root of the F₁₆ factor of F₂[Z₅];
ladder weight Σ_{j=0}^{14}|(1+x)^j p| = 40 for every even p ≠ 0.)
**Every one-sided split (a,0)/(0,a) with 1 ≤ a ≤ 39 is dead**; in
particular all of weight ≤ 7.  ∎

### 3.3 The y-span lemma: (1,1), (1,3), (1,5) are dead

**Lemma (span).**  If A⋆T = B+g then T (as a set) cannot have
|T| ≤ 5.

*Proof.*  Decompose chains into rows T_v ∈ F₂[Z₅] (v ∈ Z₁₅).  Since
A = 1 + x + y: (A⋆T)_v = (1+x)T_v + T_{v−1}.  The image B+g has
exactly three nonzero rows, with min gap 2 and min-arc span 6
(C0.10).  Let the maximal cyclic blocks of consecutive nonzero rows
of T be [a₁,b₁], … .  For each block: the image row b+1 equals
T_b ≠ 0; the image row a equals (1+x)T_a, nonzero unless T_a = ring
(C0.12) — which forces |T| ≥ 5 with T = the x-ring, whose image
y·ring has one nonzero row of weight 5 ≠ B+g's three rows: excluded.
Now:
* a length-1 block gives adjacent nonzero image rows a, a+1 —
  impossible (min gap 2);
* hence all blocks have length ≥ 2, so |T| ≤ 5 admits at most two
  blocks;
* two blocks produce ≥ 4 pairwise-distinct nonzero image rows
  (a₁, b₁+1, a₂, b₂+1; separation by zero rows keeps them distinct)
  — impossible (exactly 3);
* one block [a,b] of length ℓ ≤ |T| ≤ 5: all nonzero image rows lie
  in the arc [a, b+1] of span ℓ ≤ 5 < 6 — impossible (span 6).  ∎

(1,1) and (1,3) are the |T| = 1, 3 sub-cases ((1,1) also dies by D2
directly); (1,5) is |T| = 5.  **Dead.**  Mirrors (3,1), (5,1) follow
by Φ (C0.4), which carries them to (1,3), (1,5) cycles of the same
code.  ∎

### 3.4 Weight formula and 4-set classes

**Lemma W.**  For P Sidon with |P| = 3 and any finite T:
wt(P⋆T) = 3|T| − 2E_P(T) + 4N_P(T), where E_P(T) = #unordered pairs
of T with difference in dP and N_P(T) = #cells c with c−P ⊆ T.
*Proof.*  Cell multiplicity m_c = |T ∩ (c−P)| ≤ 3; Σ_c C(m_c,2) =
E_P(T) (Sidon: each T-pair with difference in dP determines exactly
one coincidence cell); N = #{m_c = 3}; combine
Σm_c = 3|T|, wt = #{m_c odd}.  ∎

4-sets have E_A ≤ 5 (max = rhombus+diagonal), and E_A − 2N_A = 4
only for class (4,0), = 3 only for classes (3,0) and (5,1)
(exhaustive class tally, probes2; (4,1) gives 2, etc.).  Components:
(3,0) is a connected tree (diam ≤ 3) or an A-oriented-triangle ⊔
singleton; (4,0), (5,1) are connected with diam ≤ 2.

### 3.5 (2,2) is dead

The equation is σ_B(δ_L)+t = σ_A(δ_R)+t′ with |σ| ∈ {4,6} on each
side (never 0 or 2: no periods, no 2-torsion), so sizes match.
* **6–6** (δ_L ∉ dB, δ_R ∉ dA): the left side is B+t ⊔ B+t+δ_L; the
  right side contains the full A-translate A+t′; its 3 cells lie in
  the two B-copies, so two share a copy, giving a difference in
  dA ∩ dB = ∅.  Dead.
* **4–4** (δ_L ∈ dB, δ_R ∈ dA): σ_B(δ_L) has 3 dB-pairs, σ_A(δ_R)
  has 0 (C0.6); translate-equal sets have equal pair-difference
  multisets.  Dead.  ∎

### 3.6 (2,4) and (4,2) are dead

Equation: A⋆T = σ_B(δ)+t with |T| = 4, so E_A − 2N_A ∈ {3,4} (§3.4).
* **Size 6** (δ ∉ dB; classes (3,0), (5,1)):
  - (3,0) = A-triangle ⊔ singleton: image = (A²+s) ⊔ (A+f); both the
    A²-translate (pair differences 2dA, C0.3) and the A-translate
    (pair differences dA) have 3 cells with pairwise differences
    outside dB, so each fits ≤ 1 cell per B-copy: 3 cells in 2
    copies — pigeonhole.  Dead.
  - (3,0) connected (diam ≤ 3) and (5,1) (diam ≤ 2): every image
    cell is t+a with a ∈ A of word length ≤ 1, so the image spread
    is ≤ diam+2 ≤ 5; but the image ⊇ B+t contains a pair at
    distance 6 (C0.5).  Dead.
* **Size 4** (δ ∈ dB; class (4,0), diam ≤ 2): image spread ≤ 4, but
  every σ_B(β) contains a pair at distance 7 (C0.5).  Dead.
(4,2) follows by Φ.  ∎

### 3.7 The (3,3) classification

Let B⋆u_L = A⋆u_R =: c with |u_L| = |u_R| = 3.  3-set classes
(E, N) ∈ {(0,0), (1,0), (2,0), (3,0), (3,1)} (E ≤ 3; N = 1 forces
the reflected translate, hence E = 3), with |c| = 9 − 2(E−2N) equal
on both sides.
* **|c| = 1**: needs E−2N = 4 > 3.  Impossible.
* **|c| = 3** ((3,0)×(3,0)): u_L = B+t (C0.2; the −B class has
  N = 1), image B²+t; u_R = A+t′, image A²+t′.  B²+t = A²+t′ means
  δ₂(B) is a translate of δ₂(A) (δ₂: doubling automorphism, |G| odd)
  ⟹ B a translate of A ⟹ dA = dB.  Dead by D2.
* **|c| = 5** ((2,0)×(2,0), cherries): the B-side image retains one
  full dB-pair inside each of the two outer B-copies: ≥ 2 dB-pairs.
  The A-side image has ≤ 1 dB-pair (C0.8).  Dead.
* **|c| = 7** ({(1,0),(3,1)} × {(1,0),(3,1)}):
  - (3,1)×(3,1): u_L = −B+t, u_R = −A+t′; c = S_B+t = S_A+t′ with
    both S-sets symmetric: S_B = S_A + s and (negating)
    S_B = S_A − s force S_A = S_A + 2s, a period; |S_A| = 7 and G
    has no order-7 elements, so s = 0 and dA = dB.  Dead.
  - (1,0)×(1,0): c = σ_B(β)+t ⊔ (B+t+s) contains the A-side's
    untouched full A-translate; pigeonhole into the two parts: the
    σ_B part has no dA-pair (C0.6), the B-copy has none (D2).  Dead.
  - (1,0)×(3,1): the untouched B-copy of the left side sits inside
    S_A+t′, which has no dB-pair (C0.7).  Dead.  (3,1)×(1,0):
    mirror via C0.7 for S_B.  Dead.
* **|c| = 9** ((0,0)×(0,0), both sumsets direct — the crux):
  1. *Transversality.*  For t ∈ u_L: B+t ⊆ c = ⊔_r (A+r), and by D2
     each A-copy holds ≤ 1 cell of B+t, hence exactly one:
     ∀(t,r) ∈ u_L×u_R: t−r ∈ A⊖B.
  2. *Localization.*  u_R ⊆ (B⊖A)+t for every t ∈ u_L, so for
     t ≠ t′: 3 ≤ |(B⊖A)+t ∩ (B⊖A)+t′| = mult_{d(A⊖B)}(t−t′).
     By C0.9: d(u_L) ⊆ dA ∪ dB; the class gives d(u_L) ∩ dB = ∅,
     so d(u_L) ⊆ dA.  Mirror: d(u_R) ⊆ dB.
  3. *Rigidity.*  No AP/coset degeneracies (C0.2) force u_L Sidon
     with d(u_L) = dA, so u_L ∈ {A+t₀, −A+t₀}; u_R ∈ {B+s₀, −B+s₀}.
  4. *Combos.*  (A+t₀, B+s₀): c = A⊕B+t₀ = A⊕B+s₀ and A⊕B is
     aperiodic (C0.11) ⟹ t₀ = s₀: **the generator column ∂₂δ_{t₀}**.
     The other three combos equate two of A⊕B, B⊖A, A⊖B up to
     translation — impossible by the row invariants (C0.10).  Dead.
  (Equivalently: survivors of the finite localization table T0 = the
  3 generator rows; the analytic derivation above replaces the
  table.)
**Conclusion**: the only weight-6 (3,3) cycles are the generator
columns ∂₂δ_g = (A+g, B+g), which are boundaries.  ∎

### 3.8 Assembly: LogicalFloor 8

Let u = (u_L, u_R) be a 1-cycle with chainWeight ≤ 7.  By §3.1 the
total weight is even.  Weight 0: u = 0 ∈ boundaries.  Weights 2, 4:
splits (1,1), (2,2), (1,3), (3,1) are dead (§§3.3, 3.5), one-sided
are dead (§3.2) — no such cycle exists.  Weight 6: one-sided dead
(§3.2); (1,5), (5,1) dead (§3.3); (2,4), (4,2) dead (§3.6); (3,3)
forces u = ∂₂δ_g ∈ boundaries (§3.7).  Hence every cycle of weight
≤ 7 lies in the boundary space, i.e. every non-boundary cycle has
weight ≥ 8: **`coverData.LogicalFloor 8`**.  ∎

Combined with the kernel-checked weight-8 witness (`Witness.lean`),
d(base) = 8 exactly.

### 3.9 Verification map (confirmation, not ingredients)

| claim | script |
|---|---|
| census: all ≤6-weight cycles = 75 generators; coset profiles | `a21_numeric_session1.py` |
| ψ/Φ, Lemma W, cluster tallies, (2,4)/(2,2) zero-match, (3,3) classes | `a21_numeric_probes2.py` |
| T0 localization table, T1/T2/T3, overlap profile, μ-ladders | `a21_numeric_probes3.py` |
| C1 set comparisons, BFS word lengths, span-lemma constants, spreads | `a21_numeric_probes4.py` |

Independent cross-check: the census (exhaustive at weight ≤ 6, plus
PAR at 7) re-derives the CaDiCaL UNSAT@7 verdict with no solver.

---

## 4. Lean formalization plan (next sessions)

Target Prop (exact): `coverData.LogicalFloor 8` in
`QEC/Stabilizer/Codes/BivariateBicycle/Z5Z15F2A6/` — own worktree
(`a21-*`), NEVER `a15-m-kernel-route`.

**Session-1 Lean status (2026-07-22):** worktree
`QECLean/.claude/worktrees/a21-logical-floor`, branch
`claude/a21-logical-floor` (based on `claude/a15-m-kernel-route`
`de4b547`, which owns `coverData`; mathlib shared by symlink).  New file
`Z5Z15F2A6/BaseFloor.lean` (+ umbrella import):

* `floorData : SmallCycleData G150` — the weights-2/4 layer via the
  parametric T2 bundle, BB108-pattern `native_decide` obligations
  (`check_four` is 2·150³ ≈ 6.75M tuples — biggest such sweep yet;
  watch build time).
* `strong_floor` — nonzero cycles weigh ≥ 6 (sharp).
* `weight6_cycle_is_boundary` — STATED, `sorry -- TODO(a21-w6)`.
  This is the single remaining sorry; it carries the whole §3.7
  split map.
* `logicalFloor_8 : coverData.LogicalFloor 8` — assembly PROVEN modulo
  the weight-6 sorry (parity + strong floor + omega squeeze weight to
  exactly 6, then the classification).

**BUILD GREEN (commit `c239d1f` on `claude/a21-logical-floor`)**:
`lake build …Z5Z15F2A6.BaseFloor` succeeds, 3361 jobs; the only sorry
is `weight6_cycle_is_boundary` (line 114, `TODO(a21-w6)`).  All four
bundle obligations verified — `check_four` (6.75M tuples) runs ~10 min
in `native_decide`; per-file recompiles re-pay that cost, so session 2
should consider splitting the obligations into a separate leaf file
before iterating.  Build traps hit & solved: `rw [coverData_baseComplex]`
under `∈ boundaries` breaks the rewrite motive (membership type depends
on the complex) — use `show … from` defeq instead.

Session-2 charter: discharge `weight6_cycle_is_boundary` per the §3
split map.  Suggested order: (i) split bookkeeping
(`card_filter_split` from BBDoubling gives the L/R partition);
(ii) one-sided via a `KernelCert`-pattern pivot certificate for
`rank(conv a150) = 71` + the 16-element weight table; (iii) the y-span
lemma (row recurrence — the only genuinely new framework piece);
(iv) pigeonhole splits with the C0.5–C0.8 constants as `decide`s;
(v) the (3,3) localization.  Fallback for stubborn splits:
translation-normalized finite sweeps (the (3,3) census is 900
translation classes — kernel-decidable with the coset trick).

**Two-stage discharge option (recommended).**  Stage A (engineering
leaf, ~1 session): translation-normalize (`translate1` +
`bbBoundary1Fn_translate1` already in BBSmallCycle), then per split use
the coset structure: the solution set of `A⋆u_R = c` is
`u₀(c) + Ann(A)` with `|Ann(A)| = 16`; provide the particular-solution
matrix and the 4-element kernel basis as data tables, certify kernel
completeness via the `KernelCert` Gaussian-pivot pattern (already in
this instance directory!), and `native_decide` the ≈2,701-class (3,3)
sweep + the smaller splits.  This makes `logicalFloor_8` sorry-free
quickly — discharging one of the three hypotheses of
`cover300_pauli_distance_eq_16_of_classification` — and the analytic
§3 arguments then replace the sweep leaf at leisure (exactly the
gross-proof floorOK staging).  Stage B (analytic flip): §4 items
(i)–(v).

Proposed lemma decomposition (bottom-up):
1. `parity_even`: ε as a `ZMod 2`-algebra hom; cycles have
   |u_L| ≡ |u_R|.  (BBSmallCycle may already have the pattern.)
2. `ann_floor_40`: kernel of conv-by-A has the 16-element explicit
   basis (KernelCert-style Gaussian pivot certificate) and all 15
   nonzero members have weight 40 (kernel `decide` over 16 vectors).
   Consumes: nothing analytic — certificate + finite check.
3. `span_lemma`: row decomposition (G ≅ ZMod 15 → rows in
   F₂[ZMod 5]); the block analysis of §3.3.  Medium difficulty; the
   payoff is (1,1),(1,3),(1,5) in one lemma.
4. `lemW`: incidence-count weight formula.  The heaviest single
   formalization; alternatively per-split finite decides can replace
   its uses (translation-reduced class tallies are small).
5. `phi_mirror`: the chain map Φ; kills (3,1),(5,1),(4,2) given
   3/6/7.
6. `split_22`, `split_24`: pigeonhole + spread arguments; the
   instance constants C0.5–C0.8 are tiny `decide`s.
7. `split_33`: transversality + localization + C0.9–C0.11 decides.
8. `logicalFloor_8`: assembly (case on chainWeight ≤ 7, split
   bookkeeping).
Fallback engineering route where analytic formalization stalls:
translation-reduced finite sweeps per split (900-class (3,3) etc.)
with the equivariance lemma — legitimate per the two-grade doctrine,
but the analytic route above is the target.

## 5. Honest status & risks

* Math: complete at analytic grade; every finite residue is a fixed
  instance-constant table (≤ 27 rows), machine-confirmed.
* The proof is *instance-specific* but the mechanisms (span lemma,
  overlap localization, σ pair tables) look portable to other
  neighbor-family BB codes — candidate generalization target for the
  A15 T5 depth program.
* Lean: not started; est. 2–4 sessions for the analytic route.
* Known correction during session 1: word-length constants (C0.5)
  were initially mis-computed by hand; all downstream margins
  re-verified by BFS.  No other reversals.
