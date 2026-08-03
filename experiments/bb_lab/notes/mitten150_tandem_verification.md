# [[150,30,10]] mitten code — independent distance certification

**Result: d = 10 confirmed** (d_X = 10 and d_Z = 10, each to solver
optimality in two independent lanes), matching Table I of
arXiv:2607.28795. First non-BB, non-abelian code through the lab SAT
stack. Artifacts: `certificates/mitten_150_30_10_{X,Z}.cert.json`
(bb-cert/v2), inputs vendored in `instances/mitten_150_30_10/`, driver
`scripts/mitten150_tandem.py`.

## Context

The mitten codes (Bhardwaj et al., *High-rate qLDPC processors*,
arXiv:2607.28795) are non-abelian lifted-product codes: 1×2 base
matrices over 𝔽₂[G], five data blocks, rate 1/5, check weight 9.
[[150,30,10]] is the smallest, over G = C₅×S₃ = `SmallGroup(30,1)`.
The paper reports exact distances for its five smallest codes from a
Gurobi integer-programming workflow — solver trust, no checkable
artifacts published. This note is the "certify their exacts" quick win
from the 2026-08-03 paper review: reproduce the distance claim with a
*different* solver family (two, in fact) and leave verifiable
witnesses behind. It doubles as the first test of the lab stack off
the bivariate-abelian reservation.

## Inputs

Check matrices and canonical logical bases are bit-identical copies of
the authors' data release (github.com/a7b/yarn @ `82fb695a`, MIT); see
`instances/mitten_150_30_10/README.md`. Re-deriving them from the
paper's Table XIII element indices needs GAP (not installed here);
instead `validate` mode checks the matrices against the paper's
structure — all 15 checks pass:

- shapes, CSS orthogonality, uniform check weight 9, k = 30
  (both check matrices full rank 60);
- released bases: `Lx` rows ∈ ker Hz, `Lz` rows ∈ ker Hx,
  `Lx·Lzᵀ = I₃₀`, canonical weights 18/10 exactly as Table I
  (so each `Lz` row already witnesses d_Z ≤ 10);
- Eq. (2) block anatomy: 2×5 pattern `[[1,0,1,0,1],[0,1,0,1,1]]` /
  `[[1,1,0,0,1],[0,0,1,1,1]]`, every nonzero 30×30 block a sum of 3
  permutation matrices, all same-column (Hx-block, Hz-block) pairs
  commute (the L/R regular-representation mechanism), and some
  same-side pair does NOT commute (the non-abelian signature — an
  abelian group algebra would commute everywhere).

## Method

Naive WCNF encoding from `bb_lab.maxsat_distance` (soft ¬v_i, H
parity rows as XOR chains, pairing indicators against a logical basis),
one run per CSS direction with the roles of (H_X, H_Z) swapped for the
Z side. Two independent lanes per side:

1. **Tandem** (MaxCDCL fork) to `s OPTIMUM FOUND`; witness re-verified
   in-process (weight = optimum, zero syndrome, anticommutes with a
   recomputed basis AND with the authors' released basis).
2. **CMS ladder** (`x_distance`, pycryptosat native XOR): UNSAT at
   w = 1..9, SAT at w = 10 — different engine, different encoding of
   the cardinality side.

Two soundness notes worth keeping:

- **`-cost-step=2` is UNSOUND here and was not passed.** The coset
  weight-parity premise requires every opposite-side check row even;
  mitten check weight is 9. With no flags the fork is behaviourally
  stock. (The parity gate in `fork_ab._all_even` would refuse it
  anyway; the point is that the mitten family is a live example of the
  premise failing.)
- The `strengthened` encoding was not used — `compute_class_action`
  is BB-translation-specific. Porting the orbit/anchor layer to
  group-algebra codes (here the G-action is left/right multiplication,
  single orbit on the canonical basis) is the natural next lever if
  bigger mitten instances get slow.

## Measurements (Apple Silicon, this checkout)

| side | tandem (naive, no flags) | CMS ladder w = 1..9 UNSAT + w = 10 SAT |
|---|---|---|
| d_X = 10 | OPTIMUM in 6.8 s | 137 s (w=9 rung 80 s) |
| d_Z = 10 | OPTIMUM in 8.0 s | 214 s (w=9 rung 135 s) |

(The run that emitted the committed certificates, 2026-08-03; an
earlier scratchpad run reproduced the same optima within noise. Rung
times grow ~2.5× per weight — consistent with the corpus experience
that the last UNSAT rung dominates.)

For scale: the same tandem configuration proves the gross code
(n = 144, d = 12) in ~5.5 s stock. A rate-1/5 code with k = 30 costs
about the same at d = 10 — the 2³⁰-classes worry is irrelevant to this
encoding (the pairing clause is a disjunction over k indicator bits,
and BnB cost tracks n and d, not 2^k).

## Trust model

Same split as everywhere in the lab: the **witnesses are checkable
artifacts** (support lists in the cert JSONs; `verify_certificate`
re-checks syndrome/pairing/weight against the pinned matrix hashes).
The **lower bound d ≥ 10 is the word of two independent solvers**
(MaxCDCL BnB optimality + CMS UNSAT rungs). No DRAT/LRAT artifacts:
CMS's XOR reasoning has no DRAT emission, and the CaDiCaL-on-Tseitin
route is the A15-documented intractable corner. A proof-producing
lower-bound lane for parity-heavy instances remains an open tooling
item (same status as the BB corpus).

## Follow-ups

- The other four "exact" mitten codes: [[200,40,12]] (C₄×D₁₀),
  [[300,60,14]] (C₁₀×S₃), [[500,100,16]] (C₅⋊C₂₀), [[540,108,18]]
  (C₉⋊C₁₂). d = 18 at n = 540 is the real stress test; expect the
  ladder lane to die first (extrapolating the rung growth), tandem to
  need the orbit layer.
- Descent probe: C₁₀×S₃ has a central involution with
  (C₁₀×S₃)/⟨ι⟩ ≅ C₅×S₃ — check whether [[300,60,14]] descends toward
  this code (cover-detect-and-descend, off-abelian port).
- Lean floor for this code (first formally verified non-abelian qLDPC
  distance) — scoped in the 2026-08-03 paper-review session notes;
  needs a registry number before any note lands.
