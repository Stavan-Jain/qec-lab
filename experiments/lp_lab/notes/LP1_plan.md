# LP1 — Non-BB lifted-product codes meeting the Bravyi et al. (arXiv:2308.07915) desiderata

**Goal.** Find a new stabilizer code (ideally a family) that is *not* a bivariate-bicycle
(BB) code — nor, preferably, any two-block group-algebra (2BGA) code — yet satisfies the
BB paper's selection criteria (its §2, verbatim list):

1. large distance `d` and high encoding rate `r ≫ 1/d²`;
2. short-depth syndrome-measurement (SM) circuit (BB: 7 CNOT layers, 8 rounds incl. init/measure),
   preferably distance-preserving (`d_circ = d`; BB only *conjectures* this for [[72,12,6]]);
3. pseudo-threshold near 1% under circuit noise (BB: 0.5–0.7%);
4. Tanner graph of thickness ≤ 2 (degree 6 = weight-6 checks, each qubit in 3 X + 3 Z checks);
5. fault-tolerant load/store, per-qubit readout/initialisation (BB: via an ancilla system
   attached to a surface code, thickness-2 preserved).

Benchmarks (weight 6, degree 6): BB [[72,12,6]], [[90,8,10]], [[108,8,10]], [[144,12,12]]
(kd²/n = 12), [[288,12,18]] (13.5), [[360,12,≤24]]; Aydin–Tamo–Barg coset codes
(arXiv:2606.17268) [[48,8,6]], [[96,8,10]], [[224,12,16]] (13.7); 2BGA covers
[[112,12,8]], [[168,16,10]], [[280,12,16]]. Lin–Pryadko (arXiv:2306.16400) Table I: no
weight-6 2BGA code with n < 102 has kd ≥ n and d > 6.

## The class

Lifted products over R = F₂[G], G finite abelian (PK convention, conjugate transpose †):

    LP(A,B):  H_X = [A ⊗ I | I ⊗ B†],   H_Z = [I ⊗ B | A† ⊗ I].

* `LP(A, b)`: A an a×a R-matrix with exactly three monomials in every row and column,
  b a weight-3 polynomial (1×1).  H_X = [A | b* I_a], H_Z = [b I_a | A†], n = 2a|G|.
  This is a two-block code [A|B],[Bᵀ|Aᵀ] with B = b⊗I scalar; BB codes are a = 1.
* `LP(A, B)`: both a×a (3,3)-regular; n = 2a²|G|.

Every X-check has weight 6 (3 A-monomials + 3 B-monomials), every qubit degree 6.

**Equivalences used to normalise the search (all are qubit/check relabelings):**
* A ↦ D₁ A D₂ for diagonal monomial D₁, D₂ (shift X-checks & R-qubits of copy i by s_i,
  Z-checks & L-qubits of copy j by t_j; b⊗I is invariant).  For a=2 with pattern
  [[2,1],[1,2]] this gives A = [[p,1],[1,q]] with one residual shift (p,q) ~ (x^s p, x^{-s} q),
  so p = 1 + x^{p₂}.  For a=3 all-monomial A, first row and column = 1.
* b ↦ x^s b; (A, b) ↦ (A(φ), b(φ)) for φ ∈ Aut(G); (A, b) ↦ (A*, b*) (X↔Z swap).
* For LP(A,B) the two factors normalise independently.

**Symmetric A ⇒ d_X = d_Z and rank H_X = rank H_Z** (relabel L↔R, conjugate).
For A = [[p,1],[1,q]] this holds, halving the distance work.

**k in the semisimple case (|G| odd, cyclic):** F₂[G] ≅ ⊕ fields; k = 2·deg gcd(b, pq+1, x^N+1)
for LP(A,b), a=2; k = 2·deg gcd(pq+1, rs+1, x^N+1) for LP(A,B) 2×2 with det A = pq+1,
det B = rs+1.  (Even orders: non-semisimple, compute numerically.)

## Why the hardware package transfers (this is the point of the class)

Lemma 2 (thickness 2 via wheel graphs) and the Table-5 depth-7 circuit correctness proof in
arXiv:2308.07915 use exactly one property beyond "A = A₁+A₂+A₃, B = B₁+B₂+B₃ permutation
matrices": every B_q commutes with every A_p.  In LP(A, b⊗I) each B_q = b_q ⊗ I is a scalar
translation, and in LP(A,B) the A_p⊗I and I⊗B_q act on different tensor factors — so both
proofs go through verbatim.  `lp_lab.structure.thickness2_certificate` checks planarity of
the two wheel layers with networkx as an independent certificate; `lp_lab.circuit` builds the
8-round Stim circuit and verifies (noiselessly) that all detectors are deterministic.

(Aydin–Tamo–Barg Lemma III.2 proves θ = 2 for weight-6 coset codes the same way; LP(A, b⊗I)
is formally in their framework with G' = Z_ℓ ≀ S_a acting on [a]×Z_ℓ and b central, but their
search covered small GAP groups, not these lifts; LP(A,B) with two non-scalar factors is
outside it — see §Novelty.)

**Toric layouts for a = 2 (new).** With A₃ = antidiag(1,1), A₁ = diag(p₁,q₁):
σ = A₃A₁⁻¹ swaps copies and σ² = translation by −(p₁+q₁).  If a monomial ratio p_i q_j
generates a cyclic factor Z_ℓ and a b-ratio generates the complementary Z_m, then
⟨σ, t⟩ ≅ Z_{2ℓ} × Z_m acts regularly and Lemma 4's construction gives a spanning
Cayley graph of Z_{4ℓ} × Z_{2m}: a toric layout in which the two long-range edges per
vertex are *copy-dependent* shifts (not a single cyclic shift as in BB).

## Non-BB / non-2BGA certificate

A 2BGA code's Tanner graph carries a regular group action, so Aut is transitive on X-checks.
`lp_lab.structure.tanner_automorphisms` (nauty, vertex classes X/Z/qubit) reports the number
of X-check orbits; ≥ 2 orbits ⇒ not a 2BGA code for any group, hence not BB.

## Pipeline

`scripts/s1_sweep_a2.py` (cyclic a=2), `scripts/s2_sweep.py` (modes ab2 / AB2 / ab3, any G):
k by bit-packed rank (~60 µs per n=144 code) → connectivity → randomized information-set
distance upper bound (QDistRnd-style, `lp_lab.distance.isd_min_weight`, early-stop) →
JSONL.  Finalists: exact d via bb_lab's CryptoMiniSat ladder (`exact_distance`), nauty orbits,
wheel-layer planarity, Stim circuit + `search_for_undetectable_logical_errors`.

Tool validation (2026-09-13): gross code reproduced ([[144,12]], ISD hits 12 both sides in
40 ms); Stim schedule on [[72,12,6]] gives all-zero noiseless detectors and d_circ ≤ 6 both
bases (paper: ≤ 6); wheel layers planar for gross (9/9 choices).

## Log

### 2026-09-13 — cyclic a = 2, LP([[p,1],[1,q]], b) over Z_ℓ, exhaustive (normalised)
| ℓ | n | best (k, d_ub) | fom | benchmark |
|---|---|---|---|---|
| 12 | 48 | (4,8) (8,4) | 5.3 | ATB [[48,8,6]] |
| 15 | 60 | (8,6) | 4.8 | — |
| 18 | 72 | (4,10) (8,6) | 5.6 | BB [[72,12,6]] |
| 21 | 84 | (6,10) | 7.1 | — |
| 24 | 96 | (4,12) (8,8) | 6.0 | ATB [[96,8,10]] |
| 36 | 144 | (8,12) | 8.0 | BB [[144,12,12]] |
Verdict: not competitive.  Structural reason: B = b⊗I makes ker(Bᵀ) two copies of the
cyclic code ker(b), and the X-logical (u,0), u ∈ ker(b), caps d by that cyclic code's
distance; over Z_ℓ that classical code is too weak at these lengths.

### 2026-09-13 — n = 144 wider classes (bivariate a=2, both-2×2 over Z₁₈, all-monomial 3×3 over Z₂₄)
All cap at (k, d_ub) = (8, 12), k=12 codes have d ≤ 8.  Same for every abelian group of order 36.

**Structural lemma (why).** For a two-block code H_X=[A|B], H_Z=[Bᵀ|Aᵀ] with A = ΣA_p, B = ΣB_q
permutation matrices, let Γ_A = ⟨A_p⟩, Γ_B = ⟨B_q⟩ ≤ Sym(m).  Since Γ_B ⊆ C(Γ_A): if Γ_A is
transitive its centraliser is semiregular, so Γ_B transitive ⇒ both regular ⇒ the code is a
2BGA code over Γ_A (A left-, B right-multiplication).  Hence *every non-2BGA weight-6 two-block
code has an intransitive Γ_A or Γ_B*, i.e. one constituent classical code ker A or ker B is a
direct sum over orbits, and d ≤ min over orbit codes.  In LP(A, b⊗I) the B-orbits have length
N = n/(2a); in LP(A,B) both factors are intransitive with orbit length aN = n/(2a).  This is the
price of leaving 2BGA, and it is why n = 144 with a = 2 (orbit length 36) cannot reach d = 12
with k = 12.  Aydin–Tamo–Barg's [[224,12,16]] has exactly the shape "Γ_A transitive, Γ_B with
two orbits of length 56", so the class is not hopeless at n ≥ 224: the sweep moves there.

### 2026-09-13 — exact k formula and the design principle (a = 2)
The left kernel of H_X = [A | b* I] with A = [[p,1],[1,q]] is z = (z₀, p*z₀) with
z₀ ∈ Ann(b) ∩ Ann(1 + p*q*); F₂[G] is Frobenius so dim Ann(J) = |G| − dim J.  Hence

    k = 2·dim_F₂ F₂[G] / (b*, 1 + pq)      (verified numerically on Z₁₄, Z₂₈, Z₆×Z₆, Z₁₂×Z₃, Z₉×Z₄),

maximal, k = 2δ_b with δ_b = dim ker b, iff pq ≡ 1 (mod (b*)).  (An earlier unconjugated
version passed only because random odd-N instances have k = 0; the s3 odd-N run used it and
is discarded.)  Distance is capped by the constituent code: (u,0), u ∈ ker b, is an X-logical,
so d ≤ D_b := d(ker b) and **fom ≤ δ_b·D_b²/(2|G|)**.  D_b is exact by 2^δ enumeration.
Z₅₆ table: b = 1+x+x³ → [56,3,32] (k=6); b = 1+x²+x⁶ = (1+x+x³)² → [56,6,16], cap 13.71 —
exactly the Aydin–Tamo–Barg [[224,12,16]] figure, so that code is (equivalent to) a member
of this class.  `scripts/s4_design_a2.py` enumerates only the max-k designs (pq ≡ 1 mod (b*))
per b, so the search is complete for maximal k at each (G, b).

### 2026-09-13 — design search results (a = 2, maximal k, all abelian groups of order 56 and 72)
`scripts/s4_design_a2.py`; d_ub from ISD (60–80 trials, early stop).  n = 4|G|.

| G | b (δ_b, D_b) | designs | best (k, d_ub) | fom≤ | X-orbits |
|---|---|---|---|---|---|
| Z₅₆ | 1+x²+x⁶ (6,16) | — | (12, ≤16) expected | 13.7 | — |
| Z₂₈×Z₂ | 1+x+x³y (6,24) | 1152 | (12, 16) | 13.71 | 2 (non-2BGA) |
| Z₁₄×Z₂² | all δ=12 (12,12) | 336 | (24, 6) | 3.9 | 1 |
| Z₇₂ | 1+x²+x⁴ (4,24) | 3456 | (8, 18) | 9.0 | — |
| Z₇₂ | 1+x²+x¹⁰ (4,24) | — | (8, 22) | 13.4 | 2 |
| Z₁₈×Z₂² | 1+x+x²z (8,18) | 216/b | (16, 16) | **14.22** | 2 (non-2BGA), \|Aut\| = 72 |
| Z₆×Z₆×Z₂ | 1+y+x (8,18) | 216/b | (16, 16) | 14.22 | 1 (transitive, \|Aut\|=144) |
| Z₃₆×Z₂, Z₁₂×Z₆, Z₂₄×Z₃ | δ=4, D=36 b's (cap 36) | running | | | |

Known weight-6 landscape at n = 288 (arXiv:2606.02418 exhaustive-ish BB search + Liang et al.):
[[288,12,18]] (fom 13.5), [[288,16,12]] (8.0), [[288,24,12]] = gross ⊕ gross.  So a certified
[[288,16,16]] would be a new best-known weight-6 parameter point at n = 288 and the first k = 16
indecomposable weight-6 code at this length with d > 12.

~~Headline candidate C288~~ (G = Z₁₈×Z₂×Z₂, k = 16): **RETRACTED** — the screening ISD used
`stop_at = 17`, so every reported value ≤ 17 means only "some logical of weight ≤ 17 was found";
a 300-trial ISD finds weight 12.  Lesson recorded: never read an early-stopped ISD value as an
estimate; at n = 288 only values > dstop are estimates.  The k = 16 designs at n = 288 therefore
all have d ≤ 17 (most ≤ 12) — consistent with the known BB [[288,16,12]].
Genuine n = 288 leads are k = 8: Z₇₂ with b = 1+x²+x¹⁰ (D_b = 24) and Z₂₄×Z₃ with b = 1+y+x²
(D_b = 24), both with ISD d_ub = 22 (> dstop, full trials) ⇒ [[288,8,≤22]], fom ≤ 13.4;
d = 22 would beat the known weight-6 BB [[288,8,20]] at equal (n, k).
Second candidate **C224**: G = Z₂₈×Z₂, b = 1+x+x³y, p = 1+xy, q = x¹³ + x²²y; k = 12, d_ub = 16;
2 X-orbits.  (Same parameters as Aydin–Tamo–Barg's coset code; ours is a monomial LP.)

Exact distance: `lp_lab.satdist.x_distance_sym` — CMS ladder with translation-orbit symmetry
breaking (v_r = 1 for one representative qubit r per G-orbit; four instances per rung, sound
because translations are code automorphisms).  Validated on the plain-SAT-certified
[[144,8,12]]: same d, rungs 14 s / 31 s vs 98 s / 348 s.

### 2026-09-13 — n = 288, k = 8 leads after deep re-screen (3000 ISD trials)
Screened "24"s all fall to 22 or 20; stable 22s exist in Z₇₂ (b = 1+x²+x¹⁰, 1+x²+x²²), Z₃₆×Z₂
(b = 1+x+x⁵y), Z₁₂×Z₆ (b = 1+y+xy²).  Headline **C288**: G = Z₇₂, b = 1 + x² + x¹⁰,
p = 1 + x⁵, q = x²¹ + x⁵⁸ (indices p=(0,5), q=(21,58), b=(0,2,10)): [[288, 8, ≤22]],
fom ≤ 13.4; 2 X-orbits (non-2BGA); **toric layout** Z₃₂ × Z₁₈ (4 layouts); wheels planar.
If d = 22 holds it beats the known weight-6 BB [[288,8,20]] (class Y, repo A33) at equal (n,k)
and is within 1% of [[288,12,18]]'s fom with two-thirds the logical qubits.

### Benchmarks refined (2026-09-13, evening)
Liang et al. (arXiv:2503.03827) optimal weight-6 generalized toric codes on twisted tori (still
2BGA/BB-type), n ≤ 400: [[120,8,12]], [[186,10,14]], [[210,10,16]], [[248,10,18]] (13.1),
[[254,14,16]] (14.1), [[294,10,20]] (13.6), [[310,10,≤22]] (15.6), [[340,16,18]] (15.2).
No [[288,8,22]] weight-6 code appears in the literature searched; the known BB at (288, 8) is
d = 20 (repo A33 certificate).  So C288 = [[288,8,≤22]] would be a new parameter point but not
a fom record; C224 ties Aydin–Tamo–Barg.  Fom records at n ≈ 250–340 remain with twisted-torus
BB codes; the LP class's added value is the non-2BGA structure with the full hardware package.
