# LP1 results — monomial lifted-product codes with the BB hardware package (draft, 2026-09-13)

**One-line claim.** The lifted products LP(A, b⊗I) with A = [[p,1],[1,q]] a 2×2 matrix over
F₂[G] (two monomials on the diagonal entries, single monomials off-diagonal) and b a weight-3
polynomial form a family of weight-6, degree-6 CSS codes that are **not** two-block
group-algebra codes (hence not bivariate bicycle codes), yet inherit every hardware property
of arXiv:2308.07915: thickness-2 Tanner graphs (explicit wheel layers), the 8-round /
depth-7-CNOT syndrome cycle, and — for suitable (p,q,b) — a toric layout on Z_{2μ}×Z_{2λ}.
Two concrete codes: **[[224, 12, ≤16]]** (Z₂₈×Z₂; torus 16×28) and **[[288, 8, ≤22]]**
(Z₇₂; torus 32×18).  Exact distances are being certified by SAT (see status table).

## 1. The family and its normal form
R = F₂[G], G abelian.  H_X = [A | b* I₂], H_Z = [b I₂ | A†] (BB codes are the a = 1 case
H_X = [A|B]).  Relabelings give the normal form A = [[p, 1], [1, q]], p = 1 + x^{p₂}, q any
weight-2 polynomial; b normalised under shifts and Aut(G).  A is symmetric ⇒ d_X = d_Z.

## 2. Theorems used (all elementary; proofs in LP1_plan.md)
* **k = 2·dim_F₂ R/(b*, 1 + pq)**; maximal k = 2·dim ker b iff pq ≡ 1 mod (b*).
  (Frobenius annihilator duality; verified numerically on 150 random codes incl. even orders.)
* **d ≤ D_b := d(ker b)** (the X-logical (u,0), u ∈ ker b), hence fom ≤ δ_b D_b²/(2|G|).
* **Thickness ≤ 2**: layers {A₂,A₃,B₃} and {A₁,B₁,B₂} are disjoint unions of planar wheel
  graphs because every B_q = b_q⊗I commutes with every A_p (the only property the BB proof uses).
  Checked with networkx planarity for every code reported (9/9 layer choices planar).
* **Depth-7 syndrome cycle**: Table 5 of the BB paper, with the same correctness proof
  (uses only B_q A_p = A_p B_q); Stim confirms all detectors deterministic in the noiseless
  circuit for every code reported.
* **Toric layout** (generalised Lemma 4): if σ = A_j A_i⁻¹ and τ = B_g B_h⁻¹ generate a group
  acting regularly on the L-qubits as an internal direct product Z_μ × Z_λ, the Cayley graph of
  Z_{2μ} × Z_{2λ} is a spanning subgraph (explicit coordinates; verified edge-by-edge).  Only
  commutation of B's with A's is used, so the lemma applies to all LP(A, b⊗I).
* **Non-2BGA certificate**: Aut(Tanner graph) has 2 orbits on X-checks (nauty); a 2BGA code is
  check-transitive.  (Transitivity lemma: commuting permutation sums with both permutation
  groups transitive ⇒ 2BGA.)

## 3. Codes
| code | G | b | p | q | [[n,k,d]] | fom | toric | X-orbits | status |
|---|---|---|---|---|---|---|---|---|---|
| C224 | Z₂₈×Z₂ | 1+x+x¹²y | 1+x¹⁷ | x³y+x⁴ | [[224,12,≤16]] | ≤13.7 | 16×28 | 2 | SAT ladder running |
| C288 | Z₇₂ | 1+x²+x¹⁰ | 1+x⁵ | x²¹+x⁵⁸ | [[288,8,≤22]] | ≤13.4 | 32×18 | 2 | SAT ladder running |
(index form: C224 = --orders 28,2 --p 0,34 --q 7,8 --b 0,2,25;  C288 = --orders 72 --p 0,5 --q 21,58 --b 0,2,10)

Context: BB [[144,12,12]] 12.0, [[288,12,18]] 13.5, [[288,8,20]] 11.1 (repo A33), [[288,16,12]] 8.0;
Aydin–Tamo–Barg coset [[224,12,16]] 13.7 (same parameters as C224; ATB also prove θ = 2 and
give a depth-8 cycle; C224 adds an explicit polynomial description and a toric layout);
Liang et al. twisted-torus BB [[248,10,18]] 13.1, [[294,10,20]] 13.6, [[254,14,16]] 14.1.

## 4. What did not work (recorded so it is not redone)
* Every LP class at n = 144 (a=2 over all groups of order 36, both-2×2 over Z₁₈, 3×3 monomial
  over Z₂₄) caps at [[144,8,12]]; k = 12 forces d ≤ 8.  Cause: constituent codes have length
  n/4 (orbit lemma).
* Two-factor LP(A,B) at n = 224: best [[224,6,≤18]].
* ISD early-stop values are not estimates (see LP1_plan.md, retracted k=16 "d=16" at n=288).

## 5. Status of certificates (fill in)

## 6. Artifacts
* Code cards (definition, permutation parts for the circuit, wheel layers, toric coordinates,
  certificates): `data/codes/cand_224_toric.json`, `data/codes/cand_288_toric.json`.
* Toric layout figures (Fig. 1B style): `data/figs/C224_toric_layout.svg` (Z₁₆×Z₂₈),
  `data/figs/C288_toric_layout.svg` (Z₃₂×Z₁₈) — `scripts/f1_toric_svg.py`.
* Theory with proofs: `notes/LP1_theory.md`.  Search logs and design tables: `data/logs/`,
  `data/design_a2_*.jsonl`, `data/sweep_*.jsonl`.
* Reproduce a code: `uv run python scripts/c1_certify.py --orders 28,2 --p 0,34 --q 7,8 --b 0,2,25 --skip-exact --no-circuit`
  (structure), `scripts/c2_satsym.py` (exact d, symmetry-broken SAT), `scripts/c4_dcirc_bposd.py`
  (d_circ upper bound, BP-OSD), `scripts/c5_code_card.py` (export).
* Both codes have minimum stabilizer weight 6 (no weight-4 stabilizers; relevant for BP decoding).

## 5. Status of certificates (2026-09-13 20:40)
* **C224**: symmetry-broken SAT ladder UNSAT at w ≤ 10, 11, 12, 13 (46 s, 237 s, 816 s, 2011 s
  under heavy load) ⇒ **d ≥ 14**; ISD (3000 trials) gives d ≤ 16.  Rungs 14, 15 running.
* **C288**: ladder restarted from w = 12 on an idle machine; ISD (3000 trials) d ≤ 22.
* **Circuit-level distance**: BP-OSD random-η method (paper §6) implemented and validated on
  [[72,12,6]] (finds 6 = Stim exact search).  Calibration on gross with 150 trials / 4 cycles
  gives only d_circ ≤ 12 (Z) / 13 (X), i.e. it does *not* recover the paper's ≤ 10 at this
  budget; for C224 (4 cycles) ≤ 20, C288 ≤ 33 — above the trivial bound d_circ ≤ d, hence
  uninformative.  A 2-cycle, 600-trial pass is running; unless it finds sub-d fault sets the
  honest statement is d_circ ≤ d with the depth-7 cycle verified correct (deterministic
  detectors), and distance preservation is open — exactly the status the BB paper has for
  its own codes other than [[72,12,6]].

## 7. The n = 288, k = 8 family (deep re-screen, 3000 ISD trials, first 488 of 1087 screened-22 designs)
Distribution: 59 designs at d_ub = 22, 402 at 20, 32 at 18.  The 22s split 43 (Z₃₆×Z₂, b ∈
{1+x+x²y, 1+x+x⁵y}) / 16 (Z₇₂, b = 1+x²+x¹⁰).  Twelve of them carry toric layouts:
Z₇₂: (p,q) ∈ {(1+x⁵, x²¹+x⁵⁸), (1+x³⁷, x²¹+x²⁶), (1+x⁴¹, x⁵⁷+x⁵⁸), (1+x⁷¹, x²⁷+x⁵⁸)} torus 32×18;
Z₃₆×Z₂, b = 1+x+x²y: (1+x¹⁴y, x⁹+x²), (1+x¹⁷y, x⁶+x²), (1+x¹⁹y, x⁹+x²), (1+x²²y, x⁶+x²) torus 16×36,
and two check-transitive ones (|Aut| = 144, torus 8×72; possibly 2BGA); b = 1+x+x⁵y: (1+x⁷, y+x²⁰y),
(1+x¹¹, y+x¹⁹) torus 8×72.  So C288 is one member of a small family, not an isolated hit.

### Certificate update (2026-09-13 21:10)
* C224: SAT ladder UNSAT through w ≤ 13 ⇒ d ≥ 14 (rung 14 running); ISD 3000 trials: d ≤ 16.
* C288: SAT ladder UNSAT at w ≤ 12 (599 s) ⇒ d ≥ 13 (rung 13 running); ISD 3000 trials: d ≤ 22.
* BP-OSD circuit-distance search (paper §6 method), 4 cycles / 150 trials: gross ≤ 12 (Z), ≤ 13 (X);
  C224 ≤ 20 / 20; C288 ≤ 33.  2 cycles / 600 trials: gross ≤ 12 (Z), **≤ 11 (X)** — below d = 12, so the
  method does detect circuit-level degradation at this budget, though not the paper's ≤ 10;
  C224 ≤ 19 (Z), X pending.  Verdict for C224/C288: no circuit-level fault set below the code
  distance found; d_circ ≤ d trivially; distance preservation open (as for most BB codes in the paper).
* 21:25 — C288 ladder stopped at its certified floor: w ≤ 12, 13 UNSAT (599 s, 1800 s) ⇒ **d(C288) ≥ 14**,
  ISD d ≤ 22.  C224 2-cycle BP-OSD X basis: ≤ 20 (inconclusive, as Z).  All cores now on C224 rung 14.
* 21:40 — C224: SAT witness at w = 16 (121 s, verified nontrivial) ⇒ d ≤ 16 by SAT as well as ISD;
  certificate now 14 ≤ d ≤ 16 (data/certs/cand_224_12_16_toric_w16_sat.json).  Rungs 14 (main ladder)
  and 15 (separate process) running in parallel.

### 2026-09-14 08:00 — solvers paused
The two C224 SAT processes (rung w = 14, 11.4 h; pre-computed rung w = 15, 10.4 h) were paused with
SIGSTOP at the user's request without finishing; resume with `kill -CONT` on the lp_lab python
processes or rerun `scripts/c2_satsym.py --lower 14 --upper 16`.  Certified state: C224
14 ≤ d ≤ 16 (SAT both sides), C288 14 ≤ d ≤ 22 (SAT floor, ISD ceiling).
