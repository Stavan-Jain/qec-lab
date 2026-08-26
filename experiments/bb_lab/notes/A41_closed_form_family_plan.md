# A41 — A constructed BB family with closed-form proven distance (charter)

**Claimed 2026-08-25.** Charter only; sessions log against it (the A38
pattern). Motivation: external feedback (M. Beverland, via the user) that
a BB family with a closed-form distance formula would get the QEC
community's attention — and the observation that the repo already owns
every ingredient except the one that matters most (a *class-wide*, rather
than per-instance, certificate for the value-carrying rung conditions).

**Companion thread**: A40 (parallel worktree session, claimed the same
day) attacks the *external* closed-form family — IBM's tour-de-gross
conjecture d = 6(2r+b−1) (arXiv:2506.03094) — rung-wise. This charter is
the *internal* route: design our own family so the proofs are possible by
construction. The two share technology (tower engine, class floors,
Bezout mechanism) and should be read together at merge.

## §0 Goal and claim shape

Construct an explicit infinite family of BB codes, indexed by odd
polynomial weight w (w = 3, 5, 7, …), with proven closed-form distance.
Two tiers, separable — Tier F is publishable alone:

- **Tier F (base family).** Members C_w = (G_w, A_w, B_w) in the
  mirrored-Sidon class (A16 hypotheses at weight w), with
  **d(C_w) = 2w exact**: the floor from the w-parametric class theorem,
  the upper bound from a constructive weight-2w witness. Since D1∧D2
  forces |G| ≥ 2w(w−1)+1 (A15 T4's structural bound), n = 2|G| = Θ(w²)
  at Sidon-efficient constructions, i.e. **d = Θ(√n) with closed-form
  constants** — toric scaling, BB-proper, k > 0 stated per member.
- **Tier D (doubled family).** One literal-lift Z₂ rung per member on a
  designed axis, certified class-wide: **d = 4w**. This is where the
  program's named frontier lives — the d-analog of A13's level-free
  Bezout witness: a *family-wide* certificate for SeamCosetFloor +
  DangerousFloor, not a per-instance census.

**Prioritization (learned, not assumed): horizontal first — vertical is
NOT excluded.** What the refutations actually establish is that no
*universal* per-rung law in d alone exists: the freeze law is refuted at
certificate tier (A38 S2: d([[720,4]]) ≥ 24 > 20), the number-only
ε-recursion is refuted on banked m\* (F2b), and deficits vary per rung
within single towers in both orders (A17 §8b, A20 §2). Those are
∀-statements about all covers; a closed-form *vertical* family is an
∃-statement about one designed tower, which they do not touch — the
weight-2 toric tower (A = 1+x, B = 1+y, alternating axes,
d_j = 2^⌊j/2⌋·L, k = 2) is a classical existence witness, proved by
geometry outside this program's odd-weight scope, and the gross chain's
6→12→18 prefix is itself consistent with a family-specific "+6 per
alternating-axis cover" law (at d = 6, ×2 and +6 coincide — the A40
thread's territory). What the refutations DO constrain is the descent
*proof strategy*: any rung-by-rung proof must carry the level's census
(the F2b kernel-shift form), not the bare number d — so a vertical
family is descent-provable only if its carried invariant is
self-similar up the tower (Q6). That is a harder design problem than
certifying one rung, hence this charter targets **horizontal** (a
parametric base class plus at most one designed rung) and scopes the
vertical route as the stretch front Q6 rather than excluding it.

## §1 Inventory (what exists, with provenance)

| ingredient | status | where |
|---|---|---|
| w = 3 class theorem: mirrored-Sidon class ⟹ d ≥ 6, μ ≥ 6, all free-Z₂ covers inherit | UNCONDITIONAL | `A16_class_theorem_writeup.md` |
| w = 5 conjecture (d ≥ 10) + port map + attack sequence | drafted; empirical gate PASSED: **450 members exactly at d = 2w = 10, zero falsifiers** (T4.2, commit 56e2cf6) | `A17_w5_class_plan.md` |
| D1∧D2 alone ⟹ floor ≥ 2w is **FALSE** (char-2 Frobenius squares; [[98,12,4]] wt-4 logical); the (iii) mirrored condition + `is_frobenius_related` gate (1887dd0) are load-bearing | falsified/guarded | two-sided/Frobenius line (memory + A16 §5.2) |
| odd-weight scope is provably necessary (parity is the engine; A11 E3: C-safe ⟹ doubling fails weight-agnostically — (M)-robustness must consume parity) | settled | A11 E3, A15-P3 |
| k-side class-wide mechanism: (R) ⟺ ε ∈ (A,B) ⟺ k preserved; level-free Bezout witnesses exist ((1+x²)B² = 1+x⁶ served every gross-tower level) | THEOREM + precedent | A12, A13 T3 |
| per-rung floor machinery: deficit-wall theorem (Lean), safe sector = ker τ₁ (L0), pushforward T2 | proved | `A17_deficit_wall.md` |
| (M)-robustness line (dangerous sector never binds alone ⟹ SF ⟺ doubles): A37.1–A37.4 stack; any-weight edge-bound falsification pass 704 D1&D2 pairs / 8,448 checks / 0 violations | active, parallel, unmerged | A37 (registry row) |
| analytic census/floor technology, Lean-complete per-instance (LightClassification, SeamCosetFloor 16, base floor 8) | done at f2a6/cover300 | A21, A22, A23 |
| `SidonConvBound` (Lean) — Sidon convolution weight bound, the natural class-wide census lever | exists, commit fa62959, **reachability unverified** (possibly unpushed) | goal-2 line (A5/A16) |
| the model Tier-D instance: [[150,8,8]] → [[300,8,16]] fully Lean-proven, unconditional | done | A15 line, `cover300_*` in QECLean |
| certification engine: `bb_lab.tower` (validate_banked, RungCell, deep fibers), `bb_lab.doubling_certify`, a5/a14 cover-cascade battery, A29 fibering (odd part) | promoted, gated | A38 S1/S2 |
| discovery-side sibling: A36(b) constructed doubling pairs (q > 13.5 by design), T5 family 5/5 exactly closed | unmerged branches `claude/a36-constructed-doubling`, `claude/high-fom-bb-codes-2356f1` | registry row |
| corpus support: 97 new certificate-tier exacts; corpus-merge apply owed | A39 done | `A39_descent_theory_validation.md` |

## §2 Fronts

**Q1 — the w = 5 class theorem.** Execute `A17_w5_class_plan.md` §2 from
wherever it actually stands (S1 audits the a17_w5 artifacts first). The
named big item is the coincidence-table rebuild (A16 §5.3's analog as a
cancellation-pattern classification). Known hazard, do not walk into it:
the pentagonal/torsion "Ann = 0" kills were w=3-specific — the T4.2 zd
data shows GF(16)/GF(64) *do* admit vanishing 5-term sums.

**Q2 — the ∀-odd-w statement.** Uniformize: which A16 pieces are
w-generic (PAR, (1,1)-by-D2, no-period, no-2-torsion — already marked
generic in the A17 port map) vs w-specific (coincidence tables, (2,2)
size stratification, the new (2,4)/(3,3)/… splits). Target:
*"G = Z_ℓ×Z_m, 4∤ℓ, 4∤m, |A| = |B| = w odd, D1 ∧ D2 ∧ (iii) ∧ (a′_w:
μ(Ann A), μ(Ann B) ≥ 2w) ⟹ d ≥ 2w, covers inherit."* Plus the
**witness half**: a mechanism producing weight-2w logicals class-wide
(the canonical cycle (A,B) is ∂₂(1) — trivial — so the witness is
genuinely non-canonical; the 450/450-exact w=5 population says witnesses
are generic in the class; find the pattern, then close it or scope Tier F
to the witness-carrying subfamily).

**Q3 — the constructive member sequence.** An explicit (G_w, A_w, B_w)
for every odd w satisfying D1/D2/(iii)/(a′_w)/non-Frobenius **by
construction** (Sidon pairs with disjoint difference sets near the
|G| = 2w(w−1)+1 bound — Singer/Bose-style B₂ constructions are the
starting point; (iii) constrains how the sets sit across the two axes).
Deliverables per member: k_w (closed form, or at minimum k_w > 0 proven
with the measured value stated), the (a′_w) discharge (per-member
classical Ann minimum now, analytic later), and the design levers of §3
built in.

**Q4 — Tier D, the safe half class-wide.** Choose the doubling axis by
design: doubling an odd axis ℓ ↦ 2ℓ keeps 4∤2ℓ, so **the cover stays
inside the class frames** and inherits d ≥ 2w for free (Theorem-B); the
rung buys the second factor of 2. Sub-items: (a) a class-wide Bezout
identity for (R)/k (design A_w, B_w so ε = 1+x^ℓ ∈ (A_w, B_w) via a
uniform char-2 identity — Frobenius squaring B² = B(x²,y²) is the
telescoping lever the gross witness used); (b) characterize the safe
sector im δ₂ = ker τ₁ for the class on that axis and prove
d_safe ≥ 4w family-wide — **this is the open frontier**; the deficit
ledger says it cannot be unconditional, so it must be a checkable
family-wide hypothesis satisfied by the Q3 design. Fallback (weaker,
still publishable): family-of-certificates — every member to some N
certified by the engine + the family law as a conjecture with a
per-member decision procedure.

**Q5 — Tier D, the dangerous half class-wide.** Two routes, race them:
(i) consume A37's (M)-robustness if it lands (SF ⟺ doubles kills Q5
outright for k-preserving literal lifts); (ii) analytic dangerous
censuses class-wide via Sidon structure (SidonConvBound + the A21–A23
pattern generalized from per-instance to parametric).

**Q6 (stretch) — vertical via census self-similarity.** The surviving
recursion carries censuses, not numbers (F2b's kernel-shift form), so a
vertical closed-form family is descent-provable exactly when the design
makes that recursion inductive: the level-r census structure mapping to
level-(r−1)'s under a uniform weight-rescaling, with the F2c
punctured-code coset-distance object as the candidate fixed point (F2c
§4 named the self-similarity verbatim, as a *risk*; here it would be
the *mechanism*). Toric is the degenerate model (its geometry makes the
carried object trivially self-similar). Falsify-first probe, cheap and
banked-data-only: compare stab/coset census band histograms and orbit
keys across the levels of the banked towers (gross x-ladder, c37xx
L3/L2/L1, a36 bb288, a33 IBM-288) — if nothing renormalizes even
approximately on any known tower, park the front with the numbers. A
positive signal feeds the Q3 design (choose members whose bottom-level
census is the renormalization fixed point). Do not let Q6 displace
Tier F/D work — it is scoped as one probe in S1 and otherwise idle
until the horizontal tiers are decided.

## §3 Design-for-certification principles (the actual novelty of this charter)

The A13 lesson, promoted to a design rule: the k-row was provable for
*every* tower level because the witness identity was level-free. Build
the family so every certificate has a level-free analog:

1. **Bezout by construction** — pick B_w so a fixed-shape identity
   places 1+x^ℓ in (B_w) (Frobenius-square telescoping).
2. **Frames**: ℓ odd (2-part of G in the m-axis or trivial), so the one
   designed rung stays class-eligible and a second same-axis rung is
   structurally excluded (4 | 4ℓ leaves the class) — the family is
   horizontal by frame arithmetic, not by discipline.
3. **Witness heredity**: design the weight-2w base witness so its class
   survives the fold (∉ im δ₂) — then the τ-lift is the weight-4w cover
   witness for free (upper bound half of Tier D).
4. **Seam uniformity**: prefer constructions where ker ∂₂ and the seam
   representatives have a w-independent shape (the A23 two-trinomial
   collapse is the model), so the SF proof is one parametric argument.

## §4 Falsify-first gates (S1 is all-empirical; no proof effort before these)

- **G1**: re-run the w = 5 gate; extend a first sweep to w = 7 frames
  (|G| ≥ 85, 4∤ℓ, 4∤m): does d = 2w = 14 hold exactly on the member
  population? Any member with d < 2w (hypothesis audit first — likely a
  gate violation) or d > 2w (kills *exactness*, not the floor) reshapes
  Tier F immediately.
- **G2**: enumerate the first 3–5 Q3 constructive members; certify
  d = 2w exactly (tower engine / ladder; claim tiers stated).
- **G3**: Tier-D screen over the w = 3 and w = 5 member populations on
  the designed axis (a5/a14 battery: k-gate, S0–S4): measure the SF pass
  rate. If SF fails class-systematically on every axis, Tier D dies
  early — record with numbers; Tier F stands alone.
- **G4**: k_w census across the populations → formula candidate before
  any proof.
- **G5**: audit `SidonConvBound` reachability (fa62959) and the a17_w5
  artifact state; audit A37's branch state (what is consumable today).

## §5 Walls and honesty guards

- **W1**: coincidence-table blow-up at w ≥ 7 (the A16 §5.3 analog grows
  combinatorially). Mitigation: census-first (the tables are finite per
  w; the ∀w statement needs the *pattern*, not each table).
- **W2**: torsion hazard — vanishing weight-w sums exist in GF(16)/GF(64)
  (T4.2), so Ann-vanishing kills are not free at w ≥ 5; (a′_w) carries
  that load as a hypothesis.
- **W3**: SF presentation-sensitivity (A11) — class-wide SF must fix a
  canonical presentation per member or quantify over the orbit.
- **W4**: the deficit ledger (walls attained at 2d−2; deficits 4, 6
  measured) proves unconditional Tier D is impossible — the theorem must
  be conditional on checkable design hypotheses. This is a feature
  (checkable ⟹ Lean-able), not a bug.
- **Superlative guards (A31)**: exact-distance code families exist
  (toric/lattice Kovalev–Pryadko 1202.0928; univariate GB 2508.09082;
  hyperbicycle Thm 8) — the claim is *BB-proper bivariate weight-w
  family, closed-form growing d, proof-carrying (and Lean-checkable)*,
  with the Lean-QEC (2605.16523) verified-SAT distinction stated both
  ways. d = Θ(√n) is toric scaling — do not oversell parameters; the
  content is the proof, the closed form, and k > 2 where achieved.

## §6 Known-false ledger (inherited; do not re-propose)

Everything in `A38_descent_generalization_plan.md` §5, plus
thread-specific: D1∧D2 ⟹ 2w floor (Frobenius counterexample); C-safe ⟹
doubling weight-agnostic (A11 E3); freeze-as-law (A38 S2); number-only
recursion (F2b); "deficit wall exactly 2d−2" as premise (A17-P3);
SAT/ISD witness weights as floors (standing).

## §7 Session map

- **S1 — gates** (all-empirical, no new claims): G1–G5 + Q3 candidate
  construction drafted and enumerated + the Q6 self-similarity probe
  (banked censuses only, no new computation). Output: verdict table +
  the re-scoped Tier F/D targets + the Q6 park-or-pursue call.
- **S2 — Q1**: the w = 5 theorem (coincidence-table census →
  classification → kills), per the A17-w5 attack sequence.
- **S3 — Q2 + Q3**: ∀-odd-w uniformization; constructive sequence
  existence + witness mechanism.
- **S4 — Q4/Q5**: Tier D on the designed axis; decide theorem vs
  family-of-certificates; consume A37 if landed.
- **S5 — write-up + Lean staging**: the w-parametric class theorem in
  QECLean (A16's Lean line generalized; SidonConvBound wired); first
  doubled members packaged on the cover300 pattern; paper positioning
  per `project_bb_paper1_positioning` (this charter IS the "family
  route" that positioning wants).
- Checkpoints: re-rank after S1 and after S2 (if Q1 stalls at the
  coincidence table, the fallback arc is Tier F at w ∈ {3, 5} +
  family-of-certificates Tier D — still a complete paper).

## §8 Relations

- **A40** (parallel): external family (tour de gross). If its P3
  chain-ring lane matures, port; if this charter's class-wide rung
  certificate lands first, it applies there.
- **A37**: Q5 route (i) consumes it; coordinate, don't duplicate.
- **A36(b)/A35(b)** (unmerged): discovery-side siblings — they search
  for instances, this charter designs them; their T5 family is a
  candidate Q3 seed if its members happen to sit in the class.
- **A38 F3/F5**: engine + complexity ledger; this charter's G-gates run
  through `bb_lab.tower.validate_banked()` first, per §6.0 discipline.
- **A39**: corpus-merge (once applied) supplies certified-d ground truth
  for G1/G3 populations.

## §9 Discipline

A38 §6.0 verbatim: falsify-first hard asserts against banked numbers
before anything new; claim tiers stated exactly; RED/AMBER/GREEN are
cost verdicts, never distance claims; no SAT on any floor's critical
path; witness weights never reported as floors; every session logs a
falsified-claims section. Scripts `a41_*.py`, data `data/a41/`,
checkpoint anything near the ~1 h background-kill wall.
