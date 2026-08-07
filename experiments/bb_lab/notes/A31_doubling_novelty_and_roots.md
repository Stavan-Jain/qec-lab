# A31 — Is the doubling argument novel? A literature verdict and a genealogy

Date: 2026-08-07. Question (verbatim from the session brief): *is the doubling
argument novel; what does it resemble; is it a special case of a more general
problem; does the approach resemble an approach to a different, better-understood
problem?*

**Method.** Four parallel web sweeps (quantum covers/lifts; BB/2BGA frontier
2004–2026; classical coding-theory roots; topology/geometry/algebra roots),
~270 live searches/fetches total, plus targeted checks run directly by the
synthesizing session: complete forward-citation graphs of the two prior-art
anchors (9 citing papers, the two closest body-checked), and the **v3 of the
closest competitor paper, posted 2026-08-06 — the day before this sweep** —
re-read at section level. Everything dated 2025–2026 was verified live, not
from model memory. Grade: researcher-grade synthesis; frontier-critical claims
(§2) independently re-fetched by the synthesizing session; classical/topology
theorem statements (§5) extracted from primary PDFs by the lane agents
(downloaded copies listed in §9). Internal cross-references: `A1_*` (the June
2026 four-lane sweep — this note re-verifies its frontier two months later and
extends the genealogy), the teaching doc `docs/teaching/bb-doubling-theorem.tex`
(template C1–C4, failure gallery, instance table), `A12`, `A17`, `A30`.

---

## 0. Headline verdict

**The doubling *theorem* is novel as of 2026-08-07. The doubling *setting* is
not, and has not been since 2025-11-17.** Precisely:

- The literal-lift cover construction (same polynomials, one cyclic factor
  multiplied), the observation that **the gross code is a double cover of
  [[72,12,6]]**, the inequality `k_h ≥ k` (any h), and odd-h distance transfer
  are **published**: Symons–Rajput–Browne (SRB), *Sequences of Bivariate
  Bicycle Codes from Covering Graphs*, arXiv:2511.13560 (v1 2025-11-17), with
  the odd-degree transfer proven concurrently at chain-complex generality by
  Guémard–Zémor (GZ), arXiv:2502.20297v2 (posted 38 minutes apart; documented
  contact). SRB v3 additionally credits a QIP 2025 talk by V. Guémard for the
  gross-is-a-double-cover observation. The A1 sweep found all of this in June
  2026; the program's doubling-theorem work was a deliberate attack on the
  gap they left, not an independent rediscovery of the setting.
- What SRB can prove stops at odd h. For even h they write (v3, verified
  2026-08-07): *"For even h, we empirically observe the same distance bounds
  and therefore conjecture that any h-cover BB code obeys the distance bounds
  d ≤ d_h ≤ hd."* **No proven distance bound of any kind exists in the
  literature for even-degree covers — in either direction.** Both the lower
  half (d₂ ≥ d) and upper half (d₂ ≤ 2d) of their conjecture are open there;
  the program's template proves the exact top value `d₂ = 2d` under checkable
  conditions, i.e., strictly more than the h = 2 case of their conjecture.
- The window was **re-verified open through yesterday**: SRB v1→v3 (v2
  2026-06-19, v3 2026-08-06) adds no even-h result; all 9 papers citing SRB
  and all 10 citing GZ were enumerated and the two title-level threats
  body-checked and cleared (§2.3); an arXiv full-text query for
  `"bivariate bicycle" AND "double cover"` returns **zero** papers.
- **But the territory is crowding.** Three groups now build on BB/2BGA covers
  without distance theorems: SRB (UCL; actively revising), Aydin–Tamo–Barg
  (arXiv:2606.17268, cover sequences of coset-2BGA codes), Hirasaki–Lee
  (arXiv:2607.28621, lift transport of logicals and surgery gadgets, posted
  2026-07-30). Any of them is one lemma away from asking the distance
  question. Publication urgency is real.

Ownership split, in one table:

| Claim | Status | Owner / citation |
|---|---|---|
| Literal-lift BB covers as a construction (voltage graphs, h-fold) | published | SRB Thm 3.1; extended by 2606.17268 §VII |
| "Gross = double cover of [[72,12,6]]" observation | published | SRB Ex. 5 (v1); Guémard QIP 2025 talk (per SRB v3) |
| k_h ≥ k for all h | published | SRB |
| Odd-h transfer: d_h ≤ hd; and d ≤ d_h when k_h = k | published ×2 | SRB Thm 4.5; GZ Prop 3.5 (free-Γ-module generality) |
| **Any proven even-h distance bound** | **open** (conjecture only) | — |
| **Exact doubling d₂ = 2d + checkable sufficient conditions (template C1–C4)** | **ours only** | teaching doc §10; A4/A8/A17/A30 |
| **k-preservation criterion: (R) ⟺ k₂ = k ⟺ 1+x^ℓ ∈ (A,B)** | **ours only** | A12 (SRB treat k_h = k as an empirical accident) |
| **Failure cartography: presentation sensitivity, k-jumps, safe-floor failures, tower freeze, 2d−2 deficit wall** | **ours only** | A11/A14/A17 (IBM's own gross→two-gross step, d 12→18 < 24, is an unremarked non-doubling cover in their tables) |
| **Solver-free analytic distance proofs for BB-proper codes at d ≥ 12–16** | **ours only** | gross d=12, cover300 d=16 (Lean, no solver); nearest prior: Lean-QEC, verified-SAT tier (§2.4) |

---

## 1. The object under assessment (one paragraph)

Base = BB code over ℤ_ℓ×ℤ_m with polynomials A, B; cover = same polynomials
over ℤ_{2ℓ}×ℤ_m; fold p (impose x^ℓ = 1), swap σ (multiply by x^ℓ), ε = 1+σ.
Template theorem (teaching doc §10): C1 (literal lift) ∧ C2 ((R): 1+x^ℓ ∈
(A,B)) ∧ C3 (both base-side floors at 2d: dangerous |b|+2m(b) ≥ 2d for every
stabilizer b; safe: every cycle with class in imΔ∖0 has weight ≥ 2d) ∧ C4
(tight witness outside imΔ) ⟹ d(cover) = 2·d(base). Verified instances:
[[72,12,6]]→[[144,12,12]] (analytic + Lean, unconditional);
[[36,4,4]]→[[72,4,8]] (Lean end-to-end); [[168,12,6]]→[[336,12,12]] both axes
(solver-exact, partial analytic); [[150,8,8]]→[[300,8,16]] (Lean,
unconditional, **no solver**); two [[180,4,10]]→[[360,4,20]] (A30 certificate
tier). Failure gallery: toric (floors fail), k-jump 8→16 over ℤ₃×ℤ₃ ((R)
fails), [[288,12,18]] (safe floor fails at scale), presentation sensitivity,
tower freeze, deficit wall at 2d−2.

---

## 2. The frontier, verified item by item

### 2.1 SRB (arXiv:2511.13560) — the mandatory comparator

v1 abstract (agent-verified verbatim): *"we prove that k_h ≥ k for any h, and
d_h ≤ hd when h is odd. Furthermore if h is odd and k_h = k, we prove the
distance lower bound d ≤ d_h. We conjecture it is always true that an h-cover
BB code of a base [[n,k,d]] BB code has parameters [[n_h = hn, k_h ≥ k,
d ≤ d_h ≤ hd]]."* Their machinery: covering map → chain map, projection and
lifting on homology, weight-preserving lifts; voltage-graph formalization
(free Γ = ℤ_u×ℤ_t deck action — the hypothesis our Smith-theory layer needs,
as A1-L2 already noted). Three independent full-text passes (v1 and v3)
confirm what is absent: no exact doubling statement anywhere (their own
Examples 6–7, [[18,8,2]]→[[36,8,4]], display exactly-doubled values without
remark, as does the gross pair); no even-h bound; no k_h = k criterion; no
account of what obstructs even h; no ε/transfer/Smith/Bockstein/parity
content; no sector split; no failure theory. The even-h obstruction they hit
is stated in their Lemma 4.4 (p∘τ = h·I = 0 mod 2 for even h; their Remark 9,
per A1-L2: even h "may require different tools").

Version history (fetched directly): v1 2025-11-17, v2 2026-06-19, **v3
2026-08-06**. v3 section-level re-read: even h still conjecture. *Watch item:
line-by-line v1↔v3 diff before any submission.*

### 2.2 GZ (arXiv:2502.20297v2) and Guémard's lift papers

Prop 3.5: Galois t-lift, t odd: ñ = tn, k̃ ≥ k, d̃ ≤ td; if k̃ = k then
d̃ ≥ d — proved via the cellular transfer homomorphism, generalized (Rem 3.6)
to any chain complex of free Γ-modules (includes BB covers; they never
mention BB). Identical even-t obstruction (π#∘τ# = t·I). Guémard's solo paper
(arXiv:2404.16736, *Lifts of quantum CSS codes*, IEEE-TIT) builds arbitrary
covers of the Tanner cone-complex of any CSS code, proves lifts of HGP = LP
codes, reports "improved relative parameters" for selected covers — **no
distance law**.

### 2.3 The 2026 citing network — cleared

All 9 citers of SRB and 10 of GZ enumerated (Semantic Scholar, 2026-08-07).
Body-checked the two threats: **Hirasaki–Lee 2607.28621** (*Lifting Lifted
Product Codes*, 2026-07-30): defers all parameter bounds to SRB ("in several
settings the parameters of the lift are lower bounded by those of the base
code [76]"), no odd/even distinction, no distance theorem — transports
surgery gadgets, not distances. **Aydin–Tamo–Barg 2606.17268** (2026-06-15):
Thm VII.1 = group-theoretic conditions for cover sequences of coset-2BGA
codes ("recovering and extending" SRB's construction); distances via
QDistRnd + integer programming; "no explicit formula for how distance scales
under covering." Remaining citers (fracton compactification 2605.19298,
copy-cup 2602.23307, multivariate multicycle 2601.18879, BB-surface algebra
2606.08771, check-weight bounds 2601.15446, thresholds 2607.21160, small
Tanner 2512.20532): none touches cover distance.

### 2.4 Adjacent frontier facts that constrain positioning

- **Every published BB-proper distance is solver-derived**: BCGMRY = ILP
  (Landahl method, their §"Code distance was computed by the mixed integer
  programming approach"); tour de gross (arXiv:2506.03094) **conjectures**
  its family law d = 6(2r+b−1) (verbatim: "We find n = 72r(r+b), k = 12 and
  conjecture d = 6(2r+b−1)"); IBM's discovery pipeline (2606.02418) calls
  exact MILP certification "prohibitively expensive" beyond a few hundred
  survivors; Eberhardt–Steffan, Lin–Pryadko, mitten codes: numeric/estimated.
- **Lean-QEC (arXiv:2605.16523, 2026-05-15)** must be cited and
  distinguished: Lean 4 **verified-SAT** distance certificates (bv_decide →
  CaDiCaL → LRAT replayed in-kernel) for BB codes up to [[90,8,10]]
  in-kernel; [[144,12,12]] handled only outside the kernel. Their own framing
  of the gap is quotable: distances "come either from non-scaling hand
  proofs, or from unverified solvers." Our tier is different (solver-free
  proof content, n = 300/d = 16 in-kernel) but their trusted base is smaller
  (LRAT checker vs `native_decide`); a paper must state both deltas honestly.
- **Tour de gross's own tables contain a doubling failure**: gross
  [[144,12,12]] → two-gross [[288,12,18]] is a literal-lift y-cover (A =
  x³+y²+y⁷ reduces mod y⁶−1 to x³+y+y²), d: 12 → 18 < 24 — unremarked by
  them, consistent with our freeze/deficit findings. Their family's r = 1
  b-step (6,6)→(12,6) *is* the gross doubling pair; our theorem proves that
  step's law; their r-steps are not covers, so the family conjecture stays
  open beyond it.
- **Prior distance-vs-extension theorems go the wrong direction**:
  Kovalev–Pryadko 2013 (hyperbicycle) Thm 5: D ≥ ⌊d/c⌋ (division, not
  multiplication); their exact c = 2 theorem (Thm 8) yields d =
  min(d₁,d₂,d̃₁,d̃₂) — preservation of a classical min. Lin–Pryadko 2BGA
  bounds likewise divide by the extension parameter (d_Z ≥ ⌈d₀/c⌉). The only
  prior *exact multiplication* of a quantum distance is **distance balancing**
  (Hastings; EKZ; Wills–Lin–Liu): D_X = d_X·d(C), D_Z = d_Z under tensoring
  with a classical code — a different operation (one-sided, weight-changing,
  block-length ×ℓ), the correct contrast class.

---

## 3. Why the literature stalled exactly at even h (the root, part 1)

The uniform stall point of SRB and GZ — transfer∘projection = h·I = 0 mod 2 —
is not an accident of their methods; it is **Maschke's theorem**. For odd h,
𝔽₂[ℤ_h] is semisimple: the cover splits by characters (CRT), the transfer is
a section up to the unit h, and distance/dimension transport is
character-theoretic bookkeeping. At h = 2, 𝔽₂[ℤ₂] = 𝔽₂[ε]/(ε²) is local
non-semisimple: the cover does **not** split, eigen-decomposition is replaced
by the radical filtration 0 → εC̃ → C̃ → εC̃ → 0, and the classical replacement
for character theory is **Smith theory** (P.A. Smith 1938; Bredon III.3.3;
the sequences are stated verbatim for free involutions in Hambleton–Savin
arXiv:1007.0495 §11, and at chain level in Hausmann, *Mod Two Homology and
Cohomology*, Prop 4.3.9). The program's machinery is exactly the
non-semisimple replacement kit:

- safe/dangerous split = the transfer LES of the double cover
  (im tr = ker p_* = εC̃ — "p_* = 0 ⟹ fiber-balanced" is its exactness);
- the Bezout criterion = vanishing of the connecting map, which for double
  covers is cup product with the characteristic class w(p) ∈ H¹ (Gysin =
  transfer sequence, Hausmann 4.7.36) — "1+x^ℓ ∈ (A,B)" is that Euler-class
  computation in ring coordinates (A12's (R) ⟺ Bezout);
- δ₁δ₂ = 0 via ℤ/4 lifts = a page-2 identity of the Bockstein spectral
  sequence (McCleary Ch. 10; Browder 1961); the deck towers ℤ₂ ⊂ ℤ₄ ⊂ ℤ₈ are
  its pages;
- the ε-filtration/Jordan-cell bookkeeping = modular representation theory of
  cyclic 2-groups (𝔽₂[ℤ_{2^r}] local; Alperin/Benson standard).

**No QEC paper has ever imported any of this for distance** (dedicated
searches, §9): Smith-type sequences appear in 2024–26 physics only for
SPT/anomaly classification (Debray et al.); the Bockstein appears in QEC only
on the logical-gate side (Cups-and-Gates arXiv:2410.16250; arXiv:2511.15224);
"doubled color codes" (Bravyi–Cross) is a name collision (d → d+2 recursion
for magic states). One sentence for a paper introduction: *the doubling
theorem works in the modular (non-semisimple) regime that both the quantum
cover literature and the classical QC structure theory explicitly route
around — even-degree covers are where character theory dies and Smith theory
begins.*

---

## 4. The classical shadow (the root, part 2): this argument in one variable

The exact statement pattern — split by sheet-sum, floor the nonzero-pushforward
branch, double the balanced branch — has a complete classical life in one
variable, where it is *easy* because weight is genuine Hamming weight, not
weight-modulo-stabilizers:

- **(u|u+v)/Plotkin (MacWilliams–Sloane §2.9):** d = min(2d₁, d₂).
- **Forney's squaring construction (IEEE-IT 1988, Lemma 1):** d(|S/T|²) =
  min[d(T), 2d(S)]; iterated: d = min[d(S_m), 2d(S_{m−1}), …, 2^m d(S₀)] —
  the RM/Barnes–Wall distance-doubling towers. His *twisted* squaring gives
  only ≥ and "can indeed improve minimum distance" — the classical
  presentation-sensitivity: the twist matters, exactly as in A11.
- **Repeated-root cyclic codes (Castagnoli–Massey–Schoeller–von Seemann,
  IEEE-IT 1991, Thm 1 + Lemma 2):** for length p^δ·n̄, d = min_t P_t·d(C̄_t)
  with P_t = ∏(digit_j+1); the proof is *literally* the doubling argument's
  skeleton: factor c = (x^n̄−1)^t v at maximal filtration level t, push
  forward by coefficient-class sums (sheet-sum), land in the descendant code,
  multiply by the level weight. **Specialized to p = 2, δ = 1: d =
  min{d(C̄₀), 2·d(C̄₁)} — an exact 1D doubling criterion.** Cleanest modern
  form: the "generalized van Lint theorem" (Chen–Ding arXiv:2402.02853,
  Thm 2.1): the char-2 double-length cyclic code *is* the Plotkin sum, d =
  min{2d(C₁), d(C₂)}. (van Lint 1991; priority note: Chen–Peterson–Weldon
  1969 already treated even-length binary cyclic codes.)
- **QC-LDPC graph covers (Smarandache–Vontobel, IEEE-IT 58(2):585–607, 2012;
  arXiv:0901.4129v2 Appendix I, Lemma 31, Eq. (26))** — *the* classical prior
  art for the sandwich, found by this sweep and not by A1. *Version-pinning
  caution (2026-08-07, verified against both arXiv PDFs): in the 2009 v1 the
  same result is **Lemma 30, Eq. (23), Appendix J** — v1 has no "Lemma 31",
  and Semantic Scholar/mirror links often serve v1; the IEEE-published
  version is typeset from the v2/accepted files.* Statement: for any
  splitting H = H⁽¹⁾+H⁽²⁾, the double-cover code H̃ = [[H⁽¹⁾,H⁽²⁾],
  [H⁽²⁾,H⁽¹⁾]] satisfies **d ≤ d̃ ≤ 2d** ("Lemma 31. The minimum Hamming
  distances of C and C̃ satisfy d_min(C) ⩽ d_min(C̃) ⩽ 2·d_min(C)"), proven
  in four lines by the pushforward split (c := c⁽¹⁾+c⁽²⁾ ∈ C; if c ≠ 0, w(c̃) ≥ w(c) ≥ d; if c = 0
  then c̃ is a repetition, w = 2w(c⁽¹⁾) ≥ 2d). Their §VII even documents a
  cover **freezing exactly at the base's permanent cap** (32 → 32) because
  the lifted blocks commute, with escape via non-commutativity, and
  pre-lifting (Mitchell–Smarandache–Costello 2014) as the systematic fix.
- **Tail-biting ↔ convolutional saturation:** d_min(TB at length M) ≤ d_free
  always (Bocharova et al. 2012 Thm 1; Tanner 1987), with equality once M
  exceeds an active-distance-slope bound (Bocharova–Handlery–Johannesson–
  Kudryashov 2002); Lally 2006 proves d_free ≥ d(QC) by the identical
  pushforward-or-recurse argument. This is the classical **tower freeze**:
  the circulant tower saturates at a polynomial-parent ceiling
  (Smarandache–Vontobel's permanent caps are r-independent).
- **Involutions on binary codes (Bouyuklieva IEEE-IT 2000):** φ(v) =
  orbit-sums, Ker φ = fixed subcode, im and fixed code dual to each other
  (φ(B)^⊥ = π(C_σ)), reconstruction bound d ≤ min{d(D*), 2d(B′)} — the
  ε-machinery in self-dual-classification clothing; Aksu et al. 2022's
  dim F_β ≥ ⌈k/2⌉ is im ε ⊆ ker ε (ε² = 0) used verbatim.
- **Graph avatar:** girth(bipartite double cover) = min(even girth, 2·odd
  girth) (Waller 1976; Gross–Tucker voltage-graph lifted-cycle theorem;
  Fossorier's mod-r cycle criterion is its QC-LDPC form). Bilu–Linial 2-lifts
  give the exact *spectral* dichotomy (old ⊔ new eigenvalues).

**Why the quantum case is not a corollary.** In the CSS setting both branches
of the SV-Lemma-31 argument break: the fold p(c̃) of a logical can be a
*stabilizer* (nonzero but homologically trivial — so the safe branch needs
minima over **cosets** of the base code, not the base distance), and c̃ can be
fiber-balanced without being a repetition of a logical (the ε-sector — so the
dangerous branch needs a genuine parity/robustness argument, A17's
(M)-machinery, whose odd-weight hypothesis is provably necessary by A11
Entry 3). This is exactly the content of the program's floors, and it is why
even d₂ ≥ d — trivial classically — is unproven in the quantum literature for
even h.

---

## 5. Genealogy summary: five bloodlines

| Bloodline | Ancestor statements | What the program adds |
|---|---|---|
| Covering-space topology | transfer/Gysin LES of double covers (Hausmann 4.3.9/4.7.36); Smith sequences (Bredon III.3.3; Hambleton–Savin §11); Bockstein SS (McCleary Ch. 10); Milnor/Wang sequence; modular rep theory of ℤ_{2^r} | quantitative weight content on the LES: floors per sector, exact value; first QEC import of the kit for distance |
| Classical 1D exact doubling | CMSS 1991 (exact repeated-root formula, δ=1: min{d̄₀, 2d̄₁}); van Lint/Chen–Ding; Plotkin; Forney squaring (+ twist-improves) | the 2-variable, homological (coset-valued) upgrade — the quadrant classical structure theory abandons (gcd(m,q)=1 assumed everywhere; repeated-root QC has no classical distance formula) |
| QC-LDPC lift engineering | SV Lemma 31 (classical d ≤ d̃ ≤ 2d sandwich); permanent caps; tail-biting saturation; pre-lifting; girth-of-cover monotonicity | the sandwich made *exact at the top* under certificates; freeze given a mechanism (pushforward d_safe ≤ d̃_safe) and a wall *value* (2d−2), which has no classical counterpart (classical freezes sit at the cap itself) |
| Systolic geometry / QEC cover-engine | sys(X̃) = min(sys_H, 2·sys_{∉H}) (folklore; explicit for orientation covers / bipartite double covers); Freedman ℤ₂-systolic freedom (cyclic covers as engine); Guth–Lubotzky, EKZ, golden codes (congruence towers); HHO fiber bundles (vertical/horizontal split = safe/dangerous, asymptotic) | the ℤ₂-*homological* systole-of-cover formula, which geometry **does not have** in exact form (only the three-term exactness); per-instance exactness instead of asymptotics |
| Certification lineage | Brouwer–Zimmermann floors (Lisoněk–Trummer Eq. (5)); Prange ISD; BMvT coset-weight NP-completeness; Landahl ILP → SAT ecosystem → Lean-QEC verified-SAT | **certified coset floors** (coset-BZ, A30) — the BZ floor formula transplanted onto the classically-hard coset object; solver-free Lean endgame at n=300/d=16 |

Direct answers to the two framing questions:

**Is it a special case of a more general problem? Yes, of four nested ones.**
(i) The h = 2 case of SRB/GZ's open conjecture d ≤ d_h ≤ hd (our theorem
proves more than the conjecture asks, on characterized instances, and our
failure data show the conjecture's interior — freeze at d, wall at 2d−2 — is
where non-doublers actually live). (ii) Distance of group-algebra codes under
base change along a quotient G̃ ↠ G — semisimple part classically solved
(CRT/Jensen/Lally floors), modular part open; the doubling theorem is the
smallest modular case (ℤ₂), and A22/A29's CRT-fibering + ε-machinery split is
exactly the semisimple/modular split. (iii) The ℤ₂-homological systole of a
free involution — whose Riemannian version has no exact theory (the reason
Freedman's freedom exists at all). (iv) The pre-saturation regime of the
QC → convolutional story: one doubling step is provable precisely because
δ = 1 has two filtration levels (P₁ = 2); iterating enters the regime CMSS
Thms 3–4 prove asymptotically bad — **the tower freeze is the classical
badness theorem seen from below**, and the ℤ-direction convolutional parent's
ceiling is the "toric bottleneck."

**Does the approach resemble a better-understood problem's approach? Yes,
two.** The proof is (a) the CMSS repeated-root argument upgraded from
subspaces to homology cosets — same factor-at-maximal-ε-level, same
sheet-sum pushforward, same level-weight multiplication, with the min-formula
replaced by two certified floors because quantum weight is coset-valued; and
(b) the Smith/transfer analysis of a free involution with, for the first
time, minimum-weight (systolic) conclusions extracted from the ε-filtration
rather than rank inequalities.

---

## 6. Pieces with no found ancestor (the genuinely new list)

1. **Coset-valued weight layer**: no classical theorem performs the
   (u|u+v)/CMSS split modulo a subcode; both floors are new objects.
2. **Certified coset floors** (A28/A30 coset-BZ): classical literature has BZ
   floors for code minima and ISD/NP-completeness for coset minima, but no
   certified-lower-bound theory for coset minima. New hybrid.
3. **The Bezout k-criterion as an iff** ((R) ⟺ k preserved ⟺ 1+x^ℓ ∈ (A,B)):
   classically k never "collides" under doubling (dimension is arithmetic);
   SRB have only k_h ≥ k. Its topological reading (Gysin connecting map = ∪
   w(p) vanishes) appears in no coding paper.
4. **The deficit wall as a value** (failures cluster at exactly 2d−2, parity
   theorem + pushforward mechanism, A17): classical freezes land *at* the
   inherited cap, not two below twice the base; nothing like it found.
5. **Exact doubling in ≥ 2 variables**: the classical exact formula is
   strictly univariate; the 2-variable char-dividing quadrant has no
   classical distance formula at all (structure theory assumes gcd(m,q)=1).
6. **Dangerous-sector parity ≥ 2d surviving the quotient** — the step with no
   analog even in Riemannian geometry (where no exact ℤ₂-systole-of-cover
   formula exists); its odd-row-weight hypothesis is provably necessary
   (A11 Entry 3), so this is a theorem with sharp classical-shadow-free
   content.

---

## 7. Recommended positioning (for the paper)

- **Frame**: "We resolve, in strengthened exact form, the even-degree case of
  the cover-distance problem posed by [SRB] and [GZ]" — never "we introduce
  BB covers." Concede explicitly: construction (SRB Thm 3.1 / GZ / Guémard
  2404.16736), the gross observation (SRB Ex. 5; Guémard QIP 2025 talk), k_h
  ≥ k, odd-h transfer (×2). The strongest single sentence available: *for
  even-degree covers, no distance bound in either direction was previously
  proven; we give checkable conditions under which the top of the conjectured
  band is attained exactly, and a structure theory for when it is not.*
- **Must-cite, must-distinguish list**: SRB 2511.13560 (+ diff v3 at
  submission); GZ 2502.20297; Guémard 2404.16736; Aydin–Tamo–Barg 2606.17268;
  Hirasaki–Lee 2607.28621; Lean-QEC 2605.16523 (verified-solver vs
  solver-free; trusted-base tradeoff stated both ways); distance balancing
  (Hastings / EKZ / 2305.00689) as the exact-multiplication contrast; KP-2013
  Thms 5/8 and LP 2306.16400 (division-direction bounds); HHO 2009.03921
  (sector-split precedent, asymptotic); tour de gross 2506.03094 (family law
  is conjectural; our theorem proves its r=1 b-step). Classical roots
  paragraph: CMSS 1991 + van Lint/Chen–Ding + Forney squaring + SV 2012
  Lemma 31/§VII + tail-biting saturation. Topology paragraph: Smith/Bredon/
  Hausmann + Bockstein (McCleary) + the Maschke stall-point observation.
- **Honesty guards** (from this sweep): don't claim "first exact distance
  statement in QEC" (lattice/toric exact distances exist — Kovalev–Pryadko
  1202.0928; univariate GB families 2508.09082; hyperbicycle Thm 8); don't
  claim "first machine-checked BB distance" (Lean-QEC, May 2026); the correct
  superlatives are "first proven distance law for an even-degree cover,"
  "first exact cover-to-base distance equality in the qLDPC/cover
  literature," and "largest solver-free machine-checked BB distances."
- **Watch items**: SRB v1↔v3 diff; the three active cover groups; locate the
  QIP 2025 Guémard talk abstract for the observation's exact public wording.

---

## 8. Bottom line

The doubling argument is **novel where it matters and derivative where that
is a strength**. Its setting became public in Nov 2025 and its skeleton is
textbook 1938–1991 mathematics — Smith theory on one side, repeated-root/
Plotkin distance analysis on the other — which is precisely what makes the
result robust and explainable. What has no precedent anywhere found, quantum
or classical, is the quantitative content on that skeleton: exact doubling
with checkable certificates in the non-semisimple two-variable regime, the
iff k-criterion, certified coset floors, the parity floor that survives the
stabilizer quotient, and the failure structure (freeze, wall, presentation
sensitivity). As of 2026-08-07 nobody else has proven *any* distance bound
for an even-degree cover of *any* quantum code family. The window is open;
three groups are adjacent to it; the classical-roots story (§4) is itself a
publishable expository asset.

---

## 9. Sweep provenance

Four lanes (background agents, 2026-08-07), each with a preserved
negative-query log in its report: **Lane Q** (quantum covers/lifts/products,
~30 searches): fiber-bundle/balanced/lifted-product distance state of the
art, distance balancing exact statements, Smith/Bockstein/transfer QEC
absence checks, Lean-QEC. **Lane B** (BB/2BGA 2004–2026, ~45 fetches incl. 4
local PDF extractions): BCGMRY Table 3 verbatim, tour-de-gross conjecture
passage verbatim, KP-2013/LP-2023 theorem texts, SRB v1 full-text absence
checks, the 40 newest arXiv "bivariate bicycle" abstracts individually
cleared, arXiv API zero-hit for "bivariate bicycle"+"double cover". **Lane C**
(classical roots, primary PDFs downloaded): CMSS 1991 from Massey's ETH
archive, SV 2012, Forney 1988, Lisoněk–Trummer, Bouyuklieva 2000, Güneri–
Özkaya–Solé chapter (scratchpad copies; session-local, re-download from the
cited URLs if needed). **Lane T** (topology/geometry/algebra): Hausmann/
Hambleton–Savin/Bredon chapter-verse, systole-cover folklore trace
(orientation covers, Waller, Gross–Tucker), Freedman 1999 construction
extraction, Fox–Weber/Milnor cyclic-cover algebra, negative results on exact
ℤ₂-homological systole formulas in geometry. Synthesizer's own checks:
Semantic Scholar citation graphs (both anchors), SRB abs page (version
history), SRB v3 HTML, Hirasaki–Lee abs+HTML, Aydin–Tamo–Barg abs+HTML,
Guémard title verification. Registry: this note claims **A31**; the
underlying question was posed as a follow-on to A30's docket close-out.
