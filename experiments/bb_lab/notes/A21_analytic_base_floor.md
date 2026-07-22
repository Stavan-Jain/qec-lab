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
