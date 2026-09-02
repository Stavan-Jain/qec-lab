# A40 — the tour-de-gross d-column: prove/certify the family's distance law rung-wise

**Claimed 2026-08-25** (session 1; worktree branch
`worktree-agent-aafc19cdb2dd71190`). Thread goal: the distance column of the
IBM tour-de-gross 2D-local BB family (arXiv:2506.03094) — decide, member by
member, what the descent machinery (A30/A32/A33/A36 calculus, `bb_lab.tower`,
A38 kernel-shift lane) can certify of the conjectured law d = 6(2r+b−1), and
what the ∀r analytic lane (A15 T-plan) would still owe. Session 1 scope: P0
source the family from the paper; P1 the internal cover lattice, mechanical;
P2 price (and close if GREEN) the first unproven member; P3 feasibility
verdict on the analytic lane. Discipline: A38 charter §6.0 verbatim
(falsify-first `validate_banked()` before any new claim; claim tiers exact;
RED/AMBER/GREEN are cost verdicts, never distance claims; no SAT on any
certificate-tier critical path; §5 known-false ledger respected).

Scripts `scripts/a40_*.py`; data `data/a40/`.

## §1 P0 — the family as the paper actually defines it (primary source)

Fetched arXiv:2506.03094v1 ("Tour de gross: A modular quantum computer based
on bivariate bicycle codes", Yoder–Schoute–Rall–Pritchett–Gambetta–Cross–
Carroll–Beverland, v1 2025-06-03, 68 pp; arXiv HTML retrieved 2026-08-25 and
text-extracted locally). The family lives in ONE passage — Section "Future
directions", item 3 "Increasing code and circuit distances":

> "Luckily, the BB code family contains some examples [Bra+24]. Also, it is
> reasonable to expect multiple "gross code families" to exist based on
> descriptions of BB codes as coupled copies of the toric code [Lia+25]. As
> an example, if we select integer r ∈ ℤ and bit b ∈ {0,1}, we can define a
> BB code with ℓ = 6(r+b), m = 6r and the same polynomials that we used for
> the gross and two-gross codes: A = 1 + y + x³y⁻¹, B = 1 + x + x⁻¹y⁻³. We
> find n = 72r(r+b), k = 12 and conjecture d = 6(2r+b−1). With fixed A, B
> polynomials independent of r, the code family is actually 2D-local, and
> therefore must satisfy the Bravyi-Poulin-Terhal bound [BPT10] on code
> parameters kd² = O(n)."

They add: assuming the conjectured distance, kd²/n → 24 as r grows (12 and
13.5 at gross/two-gross; surface code = 1).

**Facts the repo notes did not previously record, all load-bearing:**

1. **b is a BIT, b ∈ {0,1}** — the family is a single zigzag
   (1,0),(1,1),(2,0),(2,1),(3,0),(3,1),… — NOT a 2-parameter grid. In
   particular there is **no (2,2) member and no [[576,12]] member of any
   kind**: n = 72r(r+b) = 576 forces r(r+b) = 8 with b ∈ {0,1}, which has
   no solution. (The A14 §16 anti-instance — literal-lift doubling covers
   of two-gross die at the safe-floor level — therefore needs NO
   reconciliation with the family: those [[576,12]] covers are not family
   members. Consistency is automatic, not delicate. See §2.4.)
2. **The polynomials are FIXED Laurent polynomials** shared by every member
   (that is exactly what makes the family 2D-local): support of A =
   {(0,0),(0,1),(3,−1)}, support of B = {(0,0),(1,0),(−1,−3)} as exponent
   vectors of (x,y). The repo's standard presentations are unit shifts of
   these: y·A = x³ + y + y² (m-independent) and x·B = x + x² + y^{m−3}
   (m-DEPENDENT through the shift). At m = 6 this is exactly the stored
   bb72/gross pair (A = x³+y+y², B = x+x²+y³); at m = 12 it is
   (x³+y+y², x+x²+y⁹), which the group automorphism y ↦ y⁷ carries to the
   stored BCGMRY two-gross presentation (x³+y²+y⁷, x+x²+y³) — verified
   mechanically in §2.1.
3. **The membership grid**: (r,b) ↦ Z_{6(r+b)} × Z_{6r}:

   | (r,b) | (ℓ,m) | n | conj. d | status |
   |---|---|---|---|---|
   | (1,0) | (6,6) | 72 | 6 | = bb72; d = 6 proven (repo, Lean) |
   | (1,1) | (12,6) | 144 | 12 | = gross; d = 12 proven (repo, Lean, unconditional) |
   | (2,0) | (12,12) | 288 | 18 | = two-gross; d = 18 certificate tier (A36) |
   | (2,1) | (18,12) | 432 | 24 | open — this session's P2 target |
   | (3,0) | (18,18) | 648 | 30 | open |
   | (3,1) | (24,18) | 864 | 36 | open |
   | (4,0) | (24,24) | 1152 | 42 | open |
4. **Distance evidence in the paper** (their Appendix A.1, Fig. 12 context):
   gross and two-gross are "two codes identified in Ref. [Bra+24]" (BCGMRY);
   distances quoted as [[144,12,12]] and [[288,12,18]]. Their own numerics:
   BP+OSD logical-operator search per BCGMRY's method — "Usually the decoder
   will output a minimum weight logical, weight d = 12 and d = 18 for the
   144 qubit gross and 288 qubit two-gross codes, respectively" — plus a
   count of **336 weight-18 X-logicals** in two-gross (decoder-sampled,
   prior-randomized; solver/heuristic grade). **No member beyond (2,0) has
   any published distance, table, or verification anywhere in the paper**;
   the d-column beyond 18 is conjecture only, and even k = 12 is stated as
   "we find" (numerical), not proven, for the general member.
5. The toric-layout definition (their §2.1): unit cell with L/R qubits and
   X/Z checks on an ℓ×m torus; X check touches both cell qubits + L to the
   north + R to the east (the toric-code part) + L at offset
   (a→,a↑) = (3,−1) + R at (b→,b↑) = (−1,−3). Both gross and two-gross use
   these offsets; the family passage fixes them for all members. This
   matches the Laurent supports in item 2 (A ↔ L-connections of the X
   check, B ↔ R-connections).

**Cross-check against repo quotes**: A31 §2.4's verbatim quote ("We find
n = 72r(r+b), k = 12 and conjecture d = 6(2r+b−1)") reproduced exactly;
A31's reading "their family's r = 1 b-step (6,6)→(12,6) is the gross
doubling pair" confirmed ((1,0)→(1,1)); its "their r-steps are not covers"
is confirmed and sharpened in §2 (no consecutive step beyond (2,0) is a
cover of ANY index).

**Repo-fact correction found while sourcing (A13 note)**: A13
`deck_tower_plan.md` §4's consistency remark calls x-ladder level 2
(Z₂₄×Z₆, n = 288) "the known [[288,12,18]]". It is not: the two-gross is
Z₁₂×Z₁₂ = family (2,0), while the x-ladder (Z_{6·2^j}×Z₆, literal lift)
leaves the family at j = 2 — (24,6) solves 6(r+b) = 24, 6r = 6 only with
b = 3 ∉ {0,1} — and A14 §13's same-axis freeze battery bounds the literal
(24,6) re-double at d ≤ 12 < 18, so it is provably a DIFFERENT code. A13's
T3 mathematics (k ≡ 12 + full deck-triviality along the x-ladder, one
level-free Bezout witness) is untouched; only that parenthetical
identification is wrong. The x-ladder is not the tour-de-gross family;
the family is the (r,b) zigzag. (Mechanical check: §2.1.)

## §2 P1 — the cover lattice (mechanical)

`scripts/a40_family_lattice.py` (gate `validate_banked()` GREEN first, 1.9 s;
4.2 s total) + `scripts/a40_swap_equivalence.py` (0.1 s);
`data/a40/{family_lattice,swap_equivalence}.json`.

### §2.1 Members and presentations

All 9 members (r ≤ 4, both b, plus (5,0); n = 72 … 1800) built from the
fixed Laurent supports; **k = 12 mechanically confirmed for every one**
(the paper's "we find k = 12", now machine-checked through n = 1800 —
beyond A13 T3, which covers only the x-ladder). Presentation identities
verified mechanically: bb72/gross = the paper supports times the unit
shifts (y·A, x·B); two-gross = paper supports under (y ↦ y⁷; y⁷·A, x·B) =
the stored BCGMRY form (x³+y²+y⁷, y³+x+x²).

### §2.2 The lattice

Method per ordered pair: (i) integer covering degree n'/n (a non-integer
kills covers of EVERY kind); (ii) existence of any group quotient
G' ↠ G (per-prime invariant criterion, brute-force cross-checked on all
|G'| ≤ 300 pairs); (iii) literal-lift test on the axis-aligned projection
(the fixed Laurent supports reduce with no collisions); (iv) for pairs
with a quotient but no axis-aligned lift: enumerate ALL isomorphic
divisor frames, transport through an explicit group isomorphism, and
search the full standard monomial move set (Aut(G) incl. shears ×
per-block monomial shifts × block swap × X↔Z duality), verifying any hit
end-to-end on the parity matrices (HX and HZ rowspaces after the explicit
qubit permutation).

**Cover edges found (deck as Z_dx × Z_dy on (x,y)):**

| cover → base | deck | kind |
|---|---|---|
| (1,1) → (1,0) | Z₂ × 1 | literal (the proven perfect doubling, d 6→12, Lean) |
| (2,0) → (1,1) | 1 × Z₂ | literal (the A36 edge, deficit 6, d 12→18 certificate) |
| (2,0) → (1,0) | Z₂ × Z₂ | literal (composite) |
| **(2,1) → (1,1)** | **Z₃ × 1** | **equivalence-mediated: the (6,12) literal quotient of (2,1) IS the gross code** — shear automorphism x ↦ x⁹y, y ↦ x¹⁰y³ + shifts, verified on HX/HZ rowspaces. **[[432,12]] is a free Z₃ cover of gross.** |
| (2,1) → (1,0) | Z₃ × Z₂ | literal |
| (3,0) → (1,0) | Z₃ × Z₃ | literal — ALL-ODD deck (2-part trivial) |
| (3,1) → (1,1) | Z₂ × Z₃ | literal |
| (3,1) → (1,0) | Z₄ × Z₃ | literal |
| (4,0) → (2,0) | Z₂ × Z₂ | literal — the pure-2 chain continues |
| (4,0) → (1,1) | Z₂ × Z₄ | literal; (4,0) → (1,0): Z₄ × Z₄ |
| **(4,1) → (1,1)** | **Z₅ × Z₂** | equivalence-mediated via the same (6,12)-frame identity |
| (4,1) → (1,0) | Z₅ × Z₄ | literal |
| (5,0) → (1,0) | Z₅ × Z₅ | literal — all-odd |

**Non-edges (decisive):** no consecutive zigzag step beyond (1,1)→(2,0)
is a cover of ANY degree — (2,0)→(2,1), (2,1)→(3,0), (3,0)→(3,1),
(3,1)→(4,0), (4,0)→(4,1), (4,1)→(5,0) all have non-integer degree
(3/2, 3/2, 4/3, 4/3, 5/4, 5/4). **(3,1) is NOT a cover of (2,1)** under
literal quotients + the standard monomial move set (degree 2; the sole
isomorphic divisor frame (12,18) fails all 864 automorphisms × 4
variants × shifts) — the twisted-descent (A10-style) question is left
open, and group theory does NOT exclude it (a y-axis Z₂ extension of
(18,12) lives on Z₁₈×Z₂₄ ≅ Z₂₄×Z₁₈). (3,1)→(2,0) and (4,1)→(2,0) are
impossible outright (no group quotient: 2-part exponent drops).

### §2.3 Answers to the P1 key questions

**(a) Does the family continue as a cover chain beyond (2,0)?** As a
CHAIN, no — consecutive steps stop being covers immediately after
two-gross. As a TREE rooted at gross, yes, in a precise sense: **every
member except (3,0) and (5,0) is a free cover of gross** ((2,0): Z₂;
(2,1): Z₃; (3,1): Z₂×Z₃; (4,0): Z₂×Z₄; (4,1): Z₅×Z₂). The exceptions are
exactly the (r odd, b = 0) members, whose 2-part (2,2) is too small to
cover gross's (4,2) — they hang off bb72 by all-odd decks (Z₃×Z₃, Z₅×Z₅).
The [[432,12]] member is NOT a cover of [[288,12,18]] (degree 3/2) — the
two d-conjecture steps 18→24 CANNOT be a single cover step; but
[[432,12]] IS a Z₃ cover of the kernel-checked [[144,12,12]], with
d-conjecture ratio 24/12 = 2 over a degree-3 cover (so "perfect doubling
arithmetic" does not even apply — the odd-deck coupling law, charter
F1-Q1, is the relevant unknown).

**(b) The [[576,12]] question dissolves.** b ∈ {0,1} means **no
(2,2) member and no [[576,12]] member exists** (§1 item 1). A14 §16's
safe-floor negatives on [[576,12]] literal-lift doubles of two-gross
concern non-members; no reconciliation is needed. (Had a (2,2) member
existed, SF-death would still only have killed the doubling-template
ROUTE, not bounded d — but the premise is empty.)

**(c) Uniform decomposition.** Every member decomposes as
bb72 → (odd part) → (2-part tower): the deck of (r,b) over bb72 is
Z_{r+b} × Z_r, whose 2-part is a free iterated-Z₂ tower (the calculus
lane) and whose odd part is Z_{odd(r+b)} × Z_{odd(r)} (the F1 lane).
Equivalently: (r,b) = a free Z₂-tower of depth v₂(r+b)+v₂(r) over the
"odd anchor" member/frame (6·odd(r+b), 6·odd(r)). The d-column therefore
splits into (i) certifying odd-anchor distances (F1 / direct census) and
(ii) climbing 2-towers (the working calculus). The zigzag ORDER of the
paper is irrelevant to the machinery.

### §2.4 Corrections logged

- The A13 §4 parenthetical ("x-ladder level 2 is the known
  [[288,12,18]]") is wrong — §1. (24,6) has k = 12 (checked; T3 intact)
  but d ≤ 12 (A14 §13 SAT witness, solver/witness grade), and it is not
  a family member.
- The task-level premise "(2,2) predicted d = 30" (and the memory note's
  implicit 2-parameter grid) is out of family: b is a bit.
- My first verifier compared `rref_ints` bases as ordered lists; the
  reduced basis is canonical only as a SET (found via a
  false-negative on a true equivalence; fixed, then the equivalence
  verified). Recorded here because any future consumer of
  `rref_ints` for span equality should use `span_eq`/sorted bases.

## §3 P2 — pricing, and the [[432,12]] closure

### §3.1 The pricing table (`scripts/a40_price_open_members.py`, 8.2 s;
`data/a40/pricing.json`; gate re-passed first)

Cost verdicts only (A35 gates + the S2 n-blind-cap correction as an
extra MITM column), per member, at the partial-floor W's and the
conjecture W = d̂ − 2:

| member | tower (n chain) | k chain | W = 22 | conjecture W | verdict at conjecture |
|---|---|---|---|---|---|
| (2,1) [[432,12]] d̂ 24 | 432/216/108/54 (yxy) | 12/12/8/8 | GREEN | 22 | **GREEN — executed, §3.2** |
| (3,0) [[648,12]] d̂ 30 | 648/324/162 (xy) | 12/8/8 | AMBER | 28 | RED (cap 11) |
| (3,1) [[864,12]] d̂ 36 | 864/432/216/108/54 | 12/12/12/8/8 | GREEN | 34 | RED (cap 14) |
| (4,0) [[1152,12]] d̂ 42 | 1152/576/288/…/18 (depth 6) | 12×5/8/8 | GREEN | 40 | RED (cap 17) |
| (4,1) [[1440,12]] d̂ 48 | 1440/720/360/180/90 | 12/12/12/8/8 | GREEN | 46 | RED (cap 20) |

Notable structure: the (4,0) tower passes literally through the
two-gross (12,12) and bb72 (6,6) frames — its L2 IS the two-gross, so
two more Z₂ rungs sit on top of the certified d = 18; the (4,1) tower
passes through a (30,6) frame. Every open member is GREEN at W = 22
(d ≥ 24 partial floors are within the demonstrated envelope on all of
them); the conjecture-W wall is W3 (fiber caps), where the S2
kernel-shift lane covers light shadows and window population is the
honest residual cost.

### §3.2 The [[432,12]] closure (`scripts/a40_tdg432_close.py`;
`data/a40/tdg432/`)

Two runs, `validate_banked()` green before each; every census a
node-exact BZ walk or a gated composed-fiber derivation; every rung
candidate re-verified in-line; no SAT anywhere.

**W = 16 shakedown (v1 descent-primary architecture, 48.8 s, banked)**:
all gates EXACT (L2 direct-vs-descent; the independent (18,3)
y-quotient re-derivation of the L1-stab census); 701 dangerous rungs
ALL PASS at target 18, seam census EMPTY ≤ 16 ⟹ **d ≥ 18 at
certificate tier**, plus the τ₀-witness (below) ⟹ 18 ≤ d ≤ 24.
By-products, all census-complete EXACT: **d((9,3)) = 6** ([[54,8,6]]),
**d((9,6)) = 10** ([[108,8,10]]), **d((18,6)) = 12** ([[216,12,12]]).

**The τ₀-witness**: the top rung is (R) with exact_base (ker τ₀* =
im p₀* = SEAM), so τ₀ of the weight-12 minimum L1-logical (whose class
is outside SEAM) is a VERIFIED nontrivial [[432,12]] X-logical of
weight 24 ⟹ d ≤ 24. (Checked end-to-end: cycle, non-stabilizer, slice
identity |τu| = 2|u|.)

**W = 22 production (v2 direct-L2-primary architecture; census phase
744 s + rung phase, two-phase checkpoint per the 1-h ops rule)**. The
v2 architecture exists because the v1 L3-fiber layer explodes at
W = 22 while the composite-rank fact **rank(p₁* ∘ p₀*) = 2** collapses
the seam-shadow classes to 3 of 255 — the direct BZ walks then price at
~7e11 nodes total. Census phase, all gates EXACT again (L2
descent-vs-direct at ≤ 14: stab 85 = 85, nontrivial 666 = 666; the
(18,3) gate; d(L2) = 10 and d(L3) = 6 and d(L1) = 12 re-derived
in-run): L2-stab ≤ 22 = 33,691 orbit reps; S1'-cosets ≤ 22 = 72,977;
im-p₁*-cosets ≤ 18 = 19,411; all-class ≤ 16 = 7,780; L1 obligations =
**44,093 dangerous + 68 seam (all weight-22) + 5,727 nontrivial ≤ 18**
orbit reps; d(L1) = 12 EXACT ⟹ the b = 0 branch at W = 22 is dead
(2·12 = 24 > 22).

**Rung phase (470 s; checkpoint reloaded with every vector re-verified:
weight, cycle, stab/seam-class membership)**: dangerous rungs
**44,093/44,093 PASS** at target 24 (453 s; lanes r≤0: 34,184,
r≤1: 7,117, r≤2: 2,091, r≤3: 589, r≤4: 54, **kernel-shift: 58** — every
deep-cap cell went through the S2 F2b lane with ALL-CYCLE windows
≤ 18, sound unconditionally because the window census is complete
there; 6 lane cross-validations equal); seam rungs **68/68 PASS**;
covariance 3/3; X↔Z duality spot-check.

> **RESULT. d([[432,12]]) = 24 — EXACT, certificate tier.** Floor:
> the complete W = 22 sweep (b = 0 branch dead by 2·d(L1) = 24 > 22;
> dangerous + seam branches closed by the batteries + G-transport;
> no SAT anywhere on the critical path; inputs d((18,6)) = 12,
> d((9,6)) = 10, d((9,3)) = 6 all census-complete in-run, gates
> exact). Upper: the verified weight-24 τ₀-witness. Z side by
> transpose duality. Wall: 744 s census + 470 s rungs ≈ 20 min
> (+ 49 s W16 shakedown).
>
> This is the **conjectured tour-de-gross value at (r,b) = (2,1)**:
> the first member beyond the BCGMRY-identified pair (gross,
> two-gross) with a proven distance anywhere (the paper has no
> numerical value for it at all), the first confirmation instance of
> the d = 6(2r+b−1) conjecture beyond the solver-known members, and
> the largest-n exact BB distance at certificate tier in this repo
> (n = 432 > 360). kd²/n = 16.
>
> Mechanism note: the closure decomposes the (2,1) value as
> 6 →(Z₃ x-cover, ×2.0)→ 12 →(Z₂ y-cover, PERFECT doubling, deficit
> 0)→ 24 — the top rung (18,6)→(18,12) is a new perfect-doubling
> instance (at a non-member base frame), while the family's own
> zigzag steps are not covers at all (§2.2).

By-products banked (all census-complete exact, certificate tier):
d((3,3)) = 2 ([[18,8,2]], exhaustive), d((9,3)) = 6, d((9,6)) = 10,
d((18,6)) = 12 — the last is the "b = 2 frame" refuting the extended
formula (§4.1.3).

## §4 P3 — the ∀r analytic lane: feasibility verdict

Time-boxed assessment (no proof attempts), against the A15 T3c route
("the tour-de-gross ∀r column runs entirely through Z₄+ frames") and the
A38 charter fronts. **Verdict: RED on "the full d-column is provable
rung-wise with current technology" — with the burden itemized below and
two AMBER/GREEN sub-lanes.** RED/AMBER/GREEN are cost/feasibility
verdicts, never distance claims.

### §4.1 The structural facts the verdict rests on (P1, mechanical)

1. **The column is not a rung-chain.** No consecutive step beyond
   (1,1)→(2,0) is a cover of any degree (§2.2), so no induction along
   the zigzag exists. "Rung-wise" can only mean: per-member 2-part
   towers + odd-deck rungs off gross/bb72 + a uniformity argument over
   the deck family Z_{(r+b)} × Z_r.
2. **Odd decks are unavoidable.** (r odd, b = 0) members have all-odd
   decks (Z₃×Z₃, Z₅×Z₅ …) and 2-part (2,2) too small to cover gross;
   (2,1)/(3,1)/(4,1) hang off gross with odd deck factors. Wall W1 is
   load-bearing for the family, not incidental.
3. **The b-bit is forced by the values, not just stated**: the "b = 2"
   frame (18,6) — same fixed polynomials, r = 1 — has **d = 12 exact**
   (measured census-complete this session, §3), while the extended
   formula 6(2r+b−1) would demand 18. Extrapolating the formula off the
   bit is FALSE. (Constraint on any future closed-form family theorem.)
4. **The +6 law forces specific floors on NON-member intermediates**:
   along the pure-Z₂ chain (2^j,0), d̂ = 2·d̂_prev + 6, so e.g. the
   (4,0) = [[1152,12,d̂=42]] tower needs its non-member mid (12,24) to
   carry d ≥ 21 (the per-rung τ-ceiling d ≤ 2·d(mid)); A14 §16's SF
   hunts on exactly that intermediate (as a two-gross double) were
   SF-negative (sampled) — so the value-carrying content on the pure-Z₂
   chain must come from the census/tower machinery (or its F2 analytic
   replacement), NOT the doubling template. The doubling theorem covers
   exactly one edge of the family ((1,0)→(1,1)) and cannot cover any
   r ≥ 2 step even in principle (the increments are +6, not ×2).

### §4.2 The named missing pieces (with owners)

- **F1-Q1, the odd-deck coupling law** (charter A38 F1; session S3 is
  next and unstarted). The family supplies best-case structure — every
  family cover is k-preserving (k ≡ 12 measured through r = 5, §2.1),
  so by Maschke the twisted sectors are homology-free (the odd analog
  of deck-triviality/(R)) — and THIS SESSION supplies exact calibration
  data the law must reproduce (all census-complete, §3):

  | Z₃ x-cover | d(base) → d(cover) | ratio |
  |---|---|---|
  | (3,3) [[18,8,2]] → (9,3) [[54,8,6]] | 2 → 6 | 3.0 |
  | (3,6) [[36,8,4]] → (9,6) [[108,8,10]] | 4 → 10 | 2.5 |
  | (6,6) [[72,12,6]] → (18,6) [[216,12,12]] | 6 → 12 | 2.0 |
  | gross → (2,1) [[432,12,·]] (via §2.2 equivalence) | 12 → 24 (conj.; ≥ 18 certified, → §3) | 2.0 |

  The ratio drifts 3.0 → 2.0 as d(base) grows: the coupling is
  element-dependent, not a deck-order multiple — consistent with A38
  S1's F2a falsification (any odd coupling law must be
  cancellation-shaped). The k-preserving norm map gives the easy half,
  d(cover) ≤ |deck|·d(base) (τ of a minimum-weight base logical);
  the FLOOR half is the open calculus.
- **The bivariate chain-ring atoms at unbounded depth** (A15 T3c's two
  named-open lemmas, inherited from A1 §4): minimum weights of
  rad^t(F₂[Z₄×Z₂]) for general (non-monomial-like) ideals, and a
  weight-distortion bound between Hamming weight and chain-ring DFT
  coordinates. The family needs them not just at Z₄-depth: the 2-parts
  of the frames are Z_{2^{1+v₂(r+b)}} × Z_{2^{1+v₂(r)}}, unbounded
  along r = 2^j. Özadam–Özbudak Thm 3.6 supplies the per-axis
  ⟨(x−1)^i⟩ ⊂ F₂[Z_{2^s}] ladders at every depth; the BIVARIATE
  general-ideal gap is the blocker, at every depth ≥ 4.
- **The odd-part growth**: odd parts 3·odd(r+b) × 3·odd(r) are
  unbounded along the column; the semisimple side needs per-character
  (BCH-style) floors per isotypic block — F1-Q1 again, now with S1's
  constraint that census weight does NOT factor along the odd CRT
  (additivity 0/33,588 at gross) — the floors must survive
  cross-sector cancellation.
- **Uniformity in r**: the k-row precedent (A13 T3's level-free Bezout
  witness) shows what a uniform statement looks like; no analog exists
  for any value-carrying floor. The charter's F2b resolution (S2) says
  the true recursion is CENSUS-CARRYING, which is exactly what does NOT
  scale to ∀r without new theory (the window-population residual).

### §4.3 The achievable sub-lanes

- **GREEN — the family k-row ∀(r,b) as a theorem.** k ≡ 12 is
  machine-checked through n = 1800 here (§2.1); the statement reduces
  to dim F₂[Z_{6(r+b)}×Z_{6r}]/(A,B) = 6 uniformly — a
  commutative-algebra statement about one fixed Laurent pair, CRT over
  the odd part × chain rings over the 2-part, with A13 T1/T2/T3
  supplying the Z₂-tower half and Maschke the odd half. Named as the
  natural first ∀r theorem of the thread (it is the paper's "we find
  k = 12" made rigorous).
- **AMBER — the bounded prefix, member by member, at certificate
  tier.** Pricing (§3): [[432]] GREEN at its conjecture-W (EXECUTED,
  §3); [[648]] AMBER at W = 22 (d ≥ 24 partial) / RED at its
  conjecture W = 28 — BOTH walls bind there (bottom (9,9) census
  1e15.1 nodes AND cap 11; the depth-2 tower has the least descent
  leverage in the family); [[864]] GREEN at W = 22 / RED at W = 34;
  [[1152]]/[[1440]] RED at conjecture W. The kernel-shift lane (S2)
  moves light-shadow cells past the cap wall, so the honest frontier
  cost is window population, as the charter already records.
- **RED — everything ∀r about d.** Two hard names (F1-Q1 floor half;
  bivariate chain-ring atoms at depth), plus the structural facts of
  §4.1. A closed-form ∀r distance law is a different research program
  (exactly the A38 F1+F2 arc), not a session.

## §5 Falsified claims (session 1)

- **The A13 §4 x-ladder parenthetical** ("level 2 is the known
  [[288,12,18]]") — WRONG; corrected in place in
  `A13_deck_tower_plan.md` (Z₂₄×Z₆ is not a family member; d ≤ 12 by
  A14 §13's witness vs two-gross's certified 18; T3's mathematics
  untouched).
- **The "(2,2) member [[576,12]] with d̂ = 30" premise** (this thread's
  own task framing + the memory note's implicit 2-parameter grid) —
  no such member exists: b ∈ {0,1} is a bit in the paper (§1).
- **The extended formula at b = 2** — the (18,6) frame (r = 1, "b = 2",
  same fixed polynomials) has d = 12 EXACT (census-complete, §3.2),
  not 6(2r+b−1) = 18: the b-bit restriction is forced by the values,
  not notation.
- **Generic swap-symmetry of the family frames** — FALSE as a generic
  claim: the swapped literal quotient is the member code at
  (6,12)↔(12,6) (verified equivalence, §2.2) but NOT at (12,18)↔(18,12)
  under the full standard monomial move set (864 auts × 4 variants ×
  shifts, exhaustive). Presentation/frame-sensitivity again (A11's
  lesson, new instance).
- **Session-internal tooling bugs caught by the gates** (recorded per
  the falsify-first discipline): (i) my first equivalence verifier
  compared `rref_ints` bases as ordered lists — the reduced basis is
  canonical only as a set (false-negative on a true equivalence;
  fixed, then verified); (ii) my first L2 gate compared a ≤ 14 direct
  census against a ≤ 12 descent collection (window mismatch — the
  666-vs-65 "failure" was the gate's own bug; fixed, and the W = 22 v2
  gate then passed EXACTLY at ≤ 14 with the same 666 nontrivial
  orbits derived independently by both routes). Neither reached any
  claim.
- (Respected, inherited: no SAT anywhere; witness weights reported as
  upper bounds only; RED/AMBER/GREEN are cost verdicts; the §5
  known-false ledger of A38 untouched.)

## §6 Residue / next steps (session 1)

1. **Lean packaging of the [[432,12,24]] certificate** — the same
   census+rung data-carriage species as A36 §10.1/c37xx (charter F5);
   this instance adds the composite-rank fact (seam shadows in
   im(p₁*∘p₀*)) as a species worth designing once.
2. **The (3,1)→(2,1) twisted-descent question** (A10-style: 2^6 sheet
   twists on the y-extension of (18,12), compared to (3,1) under
   monomial moves) — decides whether the b = 1 members chain by
   twisted Z₂ covers; group theory does not exclude it (§2.2).
3. **[[648,12]] = (3,0)**: the cheapest un-closed member. AMBER at
   W = 22 (a d ≥ 24 partial floor is in the envelope today); its
   conjecture W = 28 is cap-RED but kernel-shift-eligible on light
   shadows — price the window populations first (S2's rule). Its
   all-odd deck (Z₃×Z₃ over bb72) also makes it the natural F1
   flagship after the charter's odd controls.
4. **The family k-row ∀(r,b) theorem** (§4.3 GREEN target): reduce
   dim F₂[Z_{6(r+b)}×Z_{6r}]/(A,B) = 6 to CRT × chain-ring analysis of
   the one fixed Laurent pair; A13 T1–T3 supply the Z₂-tower half,
   Maschke the odd half. The natural home for the thread's first ∀r
   statement.
5. **F1-Q1 calibration bank** (§4.2): the Z₃-cover ratio table
   (3.0 → 2.5 → 2.0 → 2.0) + the k-preserving norm upper bound
   d(cover) ≤ |deck|·d(base); any odd coupling law must reproduce the
   drift and be cancellation-shaped (S1's F2a constraint).
6. **Corpus merge**: new certificate-tier exacts from this session —
   [[54,8,6]] (9,3), [[108,8,10]] (9,6), [[216,12,12]] (18,6),
   [[18,8,2]] (3,3), and the [[432,12]] result (§3.2) — plus the
   member k-row; feed A39's corpus-merge pipeline when its apply
   lands.
7. **Promote into `bb_lab.tower`**: the composite-rank screen output
   (rank p_{i+1}*∘p_i* per adjacent pair — it collapsed the seam-shadow
   fan-out 63 → 3 here), and the v2 direct-L2-primary census pattern as
   a library lane (with the W16 v1 run as its regression battery).
8. **Session 2 of this thread**: (a) close [[648]] partials; (b) the
   (4,0) tower over the certified two-gross (its L2) — the first
   member whose closure would consume an already-certified member
   mid-level, A36-style; (c) the k-row theorem; (d) coordinate with
   the A38 S3 odd-controls session (this thread's odd rungs are its
   best-case instances).


## §7 SESSION 2 (2026-08-26) — the strongest honest attack on the conjecture

Directive: attempt the proof of d = 6(2r+b−1) ∀(r,b). Session-1's P3
verdict (RED ∀r with current technology) was the standing prior; the
attack decomposed into the upper half (P4), floor extensions (P5), and
the lower-half skeleton (P6). Scripts `a40_s2_*.py`, data
`data/a40/s2_*.json`; `validate_banked()` green before every stage.

### §7.1 P4 — the upper half: THEOREM for the b = 1 column; open for b = 0, r ≥ 3

**Falsify-first archaeology first** (`a40_s2_staircase.py`,
`s2_staircase_archaeology.json`): the "staircase of translated weight-6
bb72 blocks" reading of 6(2r+b−1) is REFUTED as a description of the
proven minimum witnesses — gross's 1,884 weight-12 logicals include 258
single-cell forms and every ragged 2-cell split (1+11 … 6+6); the a36
two-gross w18 witness crosses ALL FOUR deck cells (weights 4,4,8,2, and
all 144 translates cross 4); and NO weight-18 sum of three w6-translates
in staircase cells is even a cycle at (2,0) (0 finds over the full
84×6⁴ same-block space + 20k mixed samples, `a40_s2_ub_hunt.py`). What
the witnesses actually show is a **y-band structure**: band weights
(12) at gross, (12,6) at two-gross, (12,12) at [[432]] — i.e.
6(2r+b−1) = one 6-band + (r−1+b) 12-bands.

**The building block**: L12, an explicit x-local weight-12 pattern
(support in an x-window of width 4; 12 exponent triples recorded in
`s2_ub_bands.json`), which is a nontrivial logical of EVERY
(ℓ ≥ 12, 6)-frame flat (kernel-checked at ℓ = 12, 18, 24, 30; the
ℓ-uniformity is structural — no x-wrap is used). L12 is precisely the
class of object behind A14 §13's freeze ("the undoubled-direction
logical, lifted"): the freeze carrier and the family's witness block
are the same thing.

> **Theorem UB(r,1) [upper half, b = 1 column].** For every r ≥ 1,
> d(C_{r,1}) ≤ 12r = 6(2r+1−1). Witness: **v_{r,1} = Σ_{j=0}^{r−1}
> y^{6j}·L12** (the y-transfer of L12 along the Z_r y-deck
> (6(r+1), 6r) → (6(r+1), 6)).
> Proof ingredients and their tiers:
> (i) weight 12r — bands y-disjoint by construction;
> (ii) v is an X-cycle — L12 is an (ℓ,6)-cycle (finite local identity,
> x-window; kernel-checked ℓ = 12..30) and transfers of cycles along
> literal-lift decks are cycles (chain-map identity);
> (iii) nontriviality — the **cylinder dual certificate** u: an
> x-PERIODIC (period 6), y-LOCAL (one band, no y-wrap) Z-side cycle
> (12 triples per period, `s2_ub_dual3.json`), which is a Z-cycle of
> every member (periodicity absorbs the x-wrap at every ℓ ≡ 0 (6);
> y-locality absorbs the y-wrap at every m ≥ 12; m = 6 = r=1 checked
> directly) and pairs ⟨u, v_{r,1}⟩ = 1 (u meets only band 0; one
> finite computation). u ∈ ker H_X and odd pairing ⟹ v ∉ rowspace
> H_X.
> **Claim tier**: closed-form witness + closed-form dual certificate;
> every ingredient kernel-verified at r = 1..6 (n ≤ 3024,
> `s2_ub_bands.json`, `s2_ub_dual3.json`); the ∀r extension rests on
> the locality/periodicity uniformity arguments stated above
> (hand-proof grade, this note; Lean packaging = residue). No SAT
> anywhere.

Corollaries banked NOW (witness grade for the codes, the witnesses
themselves exact): **d([[864,12]]) ≤ 36, d([[1440,12]]) ≤ 48,
d([[2160,12]]) ≤ 60, d([[3024,12]]) ≤ 72** — the first upper bounds of
any kind on the (3,1), (4,1), (5,1), (6,1) members.

**The b = 0 column stays OPEN for r ≥ 3, with the obstruction mapped**:
the natural generalizations of the (2,0) witness all fail —
U0 (its 12-band) x-winds and does not self-stack (the U0-stack + U1
forms are non-cycles at (3,0)..(6,0), kernel-checked); the pure
y-transfer and the diagonal/helical transfers of L12 both overshoot to
12r (the −6 is exactly one band's discount); a weight-6 x-winding seam
band would need x-spacing ℓ/6 = r, i.e. an r-DEPENDENT pattern, so no
fixed local identity supplies it. The b = 0 upper half needs a genuinely
new seam construction (named residue), or the conjecture's b = 0 values
are wrong above r = 2 — no evidence either way.

### §7.2 P6 — the lower half: the per-cell lemma is dead; what survives

The proposed reduction "(i) staircase cell-count [geometry] + (ii) any
nontrivial logical pays ≥ 6 = d(bb72) per crossed cell" **fails at
(ii), in every form, on certified data** (`a40_s2_band_audit.py`):

- per-cell-each ≥ 6: the a36 two-gross witness has a cell of load 2;
  certified (18,6)-logicals have cells of load 1;
- total ≥ 6·(cells crossed): the a36 witness crosses 4 cells at weight
  18 < 24; five certified (18,6)-logicals of weight 16 cross 3 cells
  (16 < 18; e.g. cell loads (2,6,8), (1,6,9), (2,3,11)).

What the audits DO support (shadow evidence, not a lemma): the
**locality dichotomy**. Min weight of certified (18,6)-logicals by
crossed-cell count t: {1: 12, 2: 12, 3: 16}; at (2,1) the session-1
sweep proves there is NO nontrivial logical ≤ 22 at all — in
particular none touching fewer y-bands than the witness. The surviving
lower-half skeleton is therefore NOT cell accounting but:
(i) y-band-sparse logicals are CYLINDER logicals (frame-independent
objects; their x-local minimum is 12 = L12, measured), and
(ii) y-spanning logicals pay per-band with band minima governed by
boundary-defect matching — the b = 0 witness's (12,6) profile shows
band minima below 12 exist when the bands interact. A lower-half proof
needs cylinder-code distance theory (min weights of the strip/cylinder
codes and their defect-coupled stacks) — connected to, but not
supplied by, F1/F2: the F2b census-carrying refutation does not block
this route (it concerned ε-recursion in d alone), but nothing in the
current toolbox proves the cylinder minima either. **Lower half ∀r:
open; the cleanest route is now named (cylinder minima + defect
coupling), with its first two required constants measured (x-local
min 12; band loads down to 1 exist).**

### §7.3 P5 — floors: corrected pricing (no new executions this session)

`a40_s2_reprice.py` (`s2_reprice.json`), v2-architecture arithmetic
(deepest BZ-able anchor + composite-rank fan-out): **(3,0)**: k chain
12/8/8 over (18,18)→(9,18)→(9,9), rank(p₁*∘p₀*) = 2 again ⟹ 3 S1'
classes; bottom-walk totals ~1e10.1 (W=14, d ≥ 16) / ~1e11.0 (W=16,
d ≥ 18) / ~1e11.9 (W=18, d ≥ 20) / ~1e13.6 (W=22, d ≥ 24: RED today).
A d((3,0)) ≥ 18 partial floor is GREEN — **EXECUTED in-session**
(`a40_s2_t30_close.py`, `s2_t30_W16.json`, 87.7 s): TWO-ROUTE design —
both mid orderings ((18,18)→(9,18)→(9,9) and (18,18)→(18,9)→(9,9)) run
as complete independent closures whose agreement replaces the
independent-quotient gate (neither mid has a second Z₂ deck). Per
route: (9,9) direct BZ censuses (node-exact; stab ≤ 16, S1'-3 ≤ 16,
im-p₁*-15 ≤ 12, all-class ≤ 8), one descent layer, complete ntrv ≤ 12
window EMPTY ⟹ d(mid) ≥ 14 both mids, b = 0 branch dead; 1,043
dangerous rungs ALL PASS at target 18 per route, seam census EMPTY
≤ 16. **d([[648,12]]) ≥ 18, certificate tier, no SAT** — the first
floor on the (3,0) member (conjecture value 30; no upper bound known). **(3,1)**: anchor (6,9)
walks are cheap (1e9.5–1e11.4) but TWO stacked descent layers (n = 216,
432) above the anchor are the real cost — architecture-bound, not
walk-bound. The b = 1 members now carry banked upper bounds (§7.1), so
partial floors there close two-sided gaps: e.g. (3,1) at W = 22 would
give 24 ≤ d ≤ 36.

### §7.4 The distance to "the conjecture is proven" (exact statement)

- **Proven / certified now**: d = 6(2r+b−1) EXACTLY at (1,0), (1,1),
  (2,0), (2,1) — four members, certificate tier or better (two are
  Lean-backed). Upper half at b = 1: theorem (§7.1) — d ≤ conjecture
  for the whole column. Upper half at b = 0: r ≤ 2 only. Lower half:
  nothing ∀r; per-member floors are certificate-executable to roughly
  d ≥ 18–24 at n ≤ 900 with today's engine.
- **Open**: b = 0 upper half r ≥ 3; ALL lower halves r ≥ 3; the ∀r
  lower half has no working route (the per-cell route died in §7.2;
  the cylinder route is named but unbuilt). The session-1 P3 verdict
  (RED ∀r with current technology) STANDS, now with the upper half's
  b = 1 column carved out as proven.

### §7.5 Falsified claims (session 2)

- **The staircase-of-w6-blocks witness mechanism** — refuted three
  ways (§7.1): witness archaeology (cell shapes ragged; 4-cell
  two-gross witness), the exhaustive same-block staircase search at
  (2,0) (0 cycles), and the translate-orbit scan (no 3-cell form).
- **The per-cell cost lemma, both forms** (the proposed lower-half
  ingredient (ii)) — refuted on certified logicals: per-cell-each ≥ 6
  dies at cell loads 1–2; total ≥ 6·(cells) dies on the a36 witness
  (18 < 24) and on five (18,6)-logicals (16 < 18) — `s2_band_audit`.
- **The band-stack dual certificate** (u = swapbar(v)) — pairs EVEN at
  every r (structural: r·c + even cross terms); replaced by the
  cylinder dual (x-periodic, y-local), which works.
- **The naive window-kernel dual** — the (18,12)-window kernel used
  order-12 y-wraps and is not frame-portable (caught when flat
  placement failed at r = 3; replaced by the true cylinder kernel).
- **U0-stacking for the b = 0 column** — the a36 witness's 12-band
  winds x and does not self-stack; (r−1)·U0 + U1 is a non-cycle at
  every (r,0), r = 3..6 (kernel-checked).
- **Session-internal tooling bug caught**: the first band audit read
  the (18,6)-frame checkpoint vectors in the (18,12) frame — its
  output was discarded and the audit re-run in the correct frame
  (§7.2); no claim consumed the bad numbers.
- (Respected: witness weights reported as upper bounds only; no SAT;
  RED/AMBER/GREEN cost verdicts; the P3 RED verdict not inflated.)

### §7.6 Residue / next steps (session 2)

1. **The b = 0 seam construction** (upper half r ≥ 3): a weight-6
   x-winding band with r-dependent spacing, or a different −6
   mechanism; the (3,0) partial floors (next item) would also bound
   how wrong the conjecture could be there.
2. **d((3,0)) ≥ 18 partial floor** — GREEN-priced (~1e11 walks + one
   descent layer, §7.3); the v2 architecture ports with (9,9) as the
   BZ anchor (composite rank 2 again). First execution slot next
   session.
3. **The cylinder-code distance theory** (the named lower-half route,
   §7.2): min weights of the y-cylinder codes at each ℓ and the
   defect-coupled band stacks; the freeze connection (L12 = the A14
   §13 carrier) suggests the A14/A17 safe-floor machinery is the
   right starting toolbox.
4. **Lean packaging of Theorem UB(r,1)**: v_{r,1}/u are closed-form;
   the kernel checks are small decides; the ∀r uniformity arguments
   are the interesting formalization content (window locality +
   periodicity) — a natural QECLean target next to BBDeckTower.
5. **The k-row ∀(r,b) theorem** (session-1 residue, unchanged) — now
   also an ingredient of UB(r,1)'s cleanest packaging.
6. Corpus-merge additions from S2: the four new b = 1 upper bounds
   (witness grade) + the L12/dual closed forms as reusable patterns.

## §8 SESSION 3 (2026-08-26) — the cylinder lane: a proven ∀(r,b) growing floor (ladder level L1)

Directive: build the cylinder distance theory (the S2 §7.6 lower-half
route), ladder L1 = any ∀r linear floor / L2 = matching constants / L3
= the full conjecture. Scripts `a40_s3_*.py`, data `data/a40/s3_*`;
`validate_banked()` green before all work. **Level reached: L1, as a
theorem — plus the two-sided Θ(√n) pin for the b = 1 column when
combined with §7.1.** L2 formulated with its gates run and its one
missing lemma named. L3 untouched (per scope).

### §8.1 Theorem L1 — the gap-dichotomy floor

Formulating the cylinder theory's sector split exposed an elementary
route that closes L1 outright before any transfer machinery: the
two sectors (y-local vs y-wrapping) are governed by a single
no-fully-local-logicals lemma, which for translation-invariant codes
is a commutative-algebra fact.

**Machine-verified constants (gate K1)**: from the fixed Laurent
supports, every Z-check reads qubit rows spanning ≤ 4 and columns
spanning ≤ 4, and every X-stabilizer row likewise (spans s_x = s_y
= 4).

**Gate K2 (the regularity certificate)**: Res_y(y·A, x·y³·B) =
1 + x + x² + x⁴ + x⁵ + x¹¹ + x¹³ ≠ 0 in F₂[x] (exact 5×5 Sylvester
determinant over F₂[x], `a40_s3_l1_gates.py`), and the y-coefficient
contents are 1 — so A and B are coprime up to units, V(A,B) is
finite, ht(A,B) = 2 = dim F₂[x^±,y^±], and the Laurent ring is
Cohen–Macaulay: **(A,B) — hence (B̄,Ā) — is a regular sequence**, so
the plane Koszul homology H₁ vanishes: every finitely-supported plane
X-cycle is a finitely-supported plane X-boundary.

> **Lemma K (no windowed logicals).** No nontrivial X-logical of any
> member has both a cyclic x-gap ≥ 4 and a cyclic y-gap ≥ 4 in its
> support. *Proof.* With both gaps ≥ 4 = span, no member check reads
> support across either wrap, so v's member syndrome equals the plane
> syndrome of its window reading ṽ; ṽ is a plane cycle, hence by
> regularity ṽ = ∂̃₂(w̃) with w̃ of finite support; reducing mod
> (x^ℓ−1, y^m−1) intertwines ∂̃₂ with the member stabilizer map, so
> v ∈ rowspace(H_X) — trivial. ∎

> **Theorem L1.** For any BB code on Z_ℓ×Z_m whose polynomial pair is
> plane-regular with per-axis check spans (s_x, s_y):
> d ≥ min(⌈ℓ/s_x⌉, ⌈m/s_y⌉). For the tour-de-gross pair (spans 4,4):
> **d ≥ min(⌈ℓ/4⌉, ⌈m/4⌉)**. *Proof.* Let v be a nontrivial X-logical.
> If its y-support has all cyclic gaps ≤ 3, consecutive occupied rows
> are ≤ 4 apart around Z_m, so ≥ ⌈m/4⌉ rows are occupied and
> |v| ≥ ⌈m/4⌉. Otherwise v has a y-gap ≥ 4; by Lemma K it has no
> x-gap ≥ 4, so ≥ ⌈ℓ/4⌉ columns are occupied. The Z side is the
> transpose code with the same spans. ∎

> **Corollary (the family floor).** d(C_{r,b}) ≥ ⌈6r/4⌉ = ⌈3r/2⌉ for
> every r ≥ 1, b ∈ {0,1} — a growing, Θ(√n) floor for the whole
> tour-de-gross family (≈ 0.18·√n at b = 0), both columns, no solver
> anywhere. **Combined with Theorem UB(r,1) (§7.1): the b = 1 column
> has PROVEN two-sided bounds ⌈3r/2⌉ ≤ d(C_{r,1}) ≤ 12r — the
> distance is pinned to Θ(√n) = Θ(r) by proof.**

**Claim tier**: hand proof; its two computational inputs are the
machine-verified span constants and the exact resultant; the
falsify-first battery (gate K3) audited **7,765 certified nontrivial
logicals** across five banked populations (bb72 w6, gross ≤ 12
recomputed census-complete, (18,6) ≤ 22, the a36 two-gross witness,
the [[432]] w24 witness) — **zero both-windowed logicals** and every
counting inequality holds. Lean packaging (resultant + counting) is a
natural, small target.

**Positioning (guarded, per A31 discipline)**: the mechanism —
no-local-logicals for translation-invariant codes via a module/variety
condition — is Haah-adjacent commutative algebra, and BB-as-
coupled-toric readings suggest such floors are folklore-adjacent; the
A31 sweep found every published BB-proper distance solver-derived and
no proven floor for this family. Claim only: *the first stated and
certified growing distance floor for the tour-de-gross family, with
the b = 1 column's Θ(√n) pinned two-sided*. A dedicated literature
pass (Haah 2013 module framework; [Lia+25]; quasi-cyclic classical
floors) is owed before any external superlative — residue.

### §8.2 The cylinder/defect theory (L2) — formulated, gated, one lemma short

The L2 target (d ≥ 12r − 6) needs the per-band cost 12, beyond
counting. The theory as formulated: minimal logicals split into
(i) the y-local sector = strip objects (y-window, no wrap), which by
Lemma K must x-wind — **gate G1-lite: the strip code has strip-k = 6
independent of ℓ ∈ {12,18} and window height h ∈ {1,2,3}** (the
x-winding half of k = 12; the y-local sector is a fixed 6-class
object at every size); their exact minima (the strip-distance
function μ(ℓ)) are the sector's constants — floor ⌈ℓ/4⌉ proven,
exact values unmeasured (the transfer computation needs the
compression lemma below or heavy pruning: named residue); and
(ii) the y-wrapping sector, whose objects touch ≥ ⌈m/4⌉ rows (proven)
and which the defect graph models: nodes = interface states (the
restriction to 4 consecutive rows), edges = per-band continuations,
weight = per-band cost; d(wrapping) ≥ r · (min mean cycle) − C.
**Gate G2 (the killers reproduced)**: the wt-16 3-cell (18,6)
logicals and the a36 witness are y-wrapping-sector objects at m = 6
and 12 (sector histogram: 6,491 y-wrap-only + 1,274 both-wrap of
7,765; zero x-wrap-only in these populations); per-row loads go down
to **1** (measured minimum), so no per-row counting can reach 12r —
the mean-cycle formulation is genuinely necessary, and the low-load
rows must be paid for by their neighborhoods (exactly what a transfer
argument captures). **Gate G3 (x-localization license)**: x-extents
of all 7,765 certified logicals concentrate at 3–10 (histogram
banked); no unboundedness observed — resolution (a) is licensed as an
observation, but the **compression lemma** ("minimal ⟹ x-extent
≤ w₀, or x-periodic-reducible") is UNPROVEN and is the one missing
ingredient between the defect graph and an L2 theorem. **Gate G4: not
run** — without the compression lemma the interface state space is
not truncatable honestly (4 rows × w₀ columns × 2 blocks ≈ 2^96 raw);
running it on an unlicensed truncation would measure nothing
citable. The F2b refutation does not block any of this, precisely
because the defect graph's nodes carry interface STATES (census-like
data), not bare numbers — the S2 kernel-shift lesson is the same
lesson.

### §8.3 Falsified claims (session 3)

- "Per-occupied-row cost ≥ 2" (an L2 shortcut candidate) — dead
  before claiming: measured minimum row load = 1 on the banked
  populations.
- "U1 witnesses small strip minima" (an early G1 reading) — U1's
  cycle-ness at (6,6) uses the y-wrap, so it is not a strip cycle;
  no strip-minimum claim survives it.
- The Koszul lemma's falsification attempt (gate K3: hunt a
  both-windowed certified logical) — 0/7,765; the lemma stands.
- (Respected: cost verdicts vs distance claims; no SAT; tier
  statements exact; the S2 per-cell refutations untouched.)

### §8.4 Residue / S4

1. **Exact strip minima μ(ℓ)** (the y-local sector constants): a
   pruned x-transfer/Viterbi at strip-k = 6, or an algebraic route
   through the 6-dim strip homology; feeds L2's x-sector.
2. **The compression lemma** (minimal ⟹ x-local/x-periodic-
   reducible): the gate to G4 and the wrapping sector's mean-cycle
   theorem; G3's histogram is its evidence base.
3. **The literature pass** (Haah's module framework, [Lia+25],
   quasi-cyclic floors) before any external claim about L1's novelty.
4. **The corpus sweep corollary**: stamp d ≥ min(⌈ℓ/s_x⌉, ⌈m/s_y⌉)
   on every corpus BB code (spans + resultant are cheap per pair) —
   free floors for every frame, including the odd-|G| rows the
   Z₂-calculus cannot touch (W1-independent!).
5. **Lean**: Theorem L1 end-to-end (span constants by `decide`, the
   resultant by exact arithmetic, the counting argument); would make
   the family floor the library's first ∀r distance statement.
6. b = 0 seam construction: untouched this session (per scope); the
   strip/defect machinery, once L2-grade, is expected to hand over
   the −6.

## §9 SESSION 4 (2026-08-26) — the compression lemma dissolved: combs,
## stacks, the phase atlas, and Theorem B6

Directive: prove the compression lemma (§8.2's one missing ingredient)
and close the b = 1 gap. What actually happened: falsify-first testing
of the §8.2 formulation itself REFUTED it before any compression proof
could matter — the defect-graph min-mean-cycle plan is vacuous as
stated — and the correct replacement architecture was then built and
gated. In it, the compression lemma's two halves DISSOLVE: the
finiteness half is free (light slabs are finite by weight, heavy slabs
pay pointwise), and the "x-periodic-reducible" half becomes a phase
classification, executed this session as an exhaustive certificate
(the compact-phase atlas, p ≤ 8, both transfer lanes). By-product: the
strip-minima residue (§8.4 item 1) is CLOSED at m = 6 for all ℓ
simultaneously, giving **Theorem B6 — d((ℓ,6)) = 12 for every
6 | ℓ ≥ 12 — the thread's first ∀ℓ exact-distance statement, all
certificates closed in-session**. Scripts
`a40_s4_*.py`, data `data/a40/s4_*.json`; `validate_banked()` green
before every stage.

### §9.1 Gate C — the comb refutation: §8.2's formulation was vacuous

`a40_s4_comb_gate.py` (`s4_comb_gate.json`). "Combs" — sums of single
X-stabilizer generators (H_X rows) tiled with y-period 6 at x-fixed
position — are TRIVIAL y-spanning cycles of weight m (rate 1/row,
max cyclic y-gap 1, max 4-row slab weight 5) at every member frame
tested ((12,12), (18,12), (18,18)); at frames with 8 | m the period-8
comb reaches rate 3/4 ((12,16), (12,24); C4 tooth sweep: 3/4 is the
cheapest tiled boundary over all 2×2-window teeth and periods 4..9).
Since the §8.2 defect graph contains these cycles, its min mean cycle
is ≤ 3/4 — far below the target rate 2 — so "d(wrapping) ≥ r·(min
mean cycle) − C" was VACUOUS as formulated, no matter what state
truncation (compression) licensed it. The compression lemma was not
the missing ingredient; triviality-blindness was. (G4 as specced in
§8.2 would have measured ≤ 1 and proven nothing.)

Also pinned mechanically (gate C2, 20 random cycles at (18,12)): the
row recurrence
E_j: (1+x^{−1})v₁[j] + x·v₁[j−3] + v₂[j] + v₂[j+1] + x^{−3}v₂[j−1] = 0
(convention `circulant(P)[g,h] = P(g−h)`, X-cycles B̄v₁ + Āv₂ = 0),
and the **block-2 determinism**: v₂[j+1] is forced by
(v₁[j], v₁[j−3], v₂[j], v₂[j−1]) — the member is a convolutional code
in y over F₂[x]/(x^ℓ−1) with block-1 rows as free inputs (and,
mirrored, a convolutional code in x with block-2 columns forced at
offset +3). Gate C3: the slab telescope Σ_j W_j = 4|v| (W_j = weight
of rows j−3..j, both blocks; cyclic, exact).

### §9.2 The reformed cost frame: slab amortization kills the
### finiteness half of compression

Reparametrize the y-walk cost as ŵ_j := W_j/4. On closed (torus)
walks Σ_j ŵ_j = |v| EXACTLY (C3), and the b = 1 target d = 12r = 2m
becomes "mean slab weight ≥ 8". Consequences:

- **Heavy steps pay pointwise.** Any row j with W_j ≥ 8 contributes
  ŵ_j ≥ 2 by itself. No state bookkeeping, no truncation license, no
  compression is needed for heavy slabs — the original lemma's
  "x-extent ≤ w₀" demand was aimed at exactly these states and is
  unnecessary.
- **The light core is finite for free.** States with W ≤ 7 carry ≤ 7
  qubits; two occupied columns interact through a check only at
  x-distance ≤ 4, and bridging an x-gap of size g inside 4 rows costs
  ≥ g/4 ≥ the whole light budget for g > 28 — so light states, taken
  up to x-translation with inter-cluster gaps saturated at "far", form
  a FINITE, r-independent set. This is the compression lemma's
  finiteness content, obtained free of charge below the heavy line.
- **L1 is the trivial case**: y-spanning ⟹ every slab nonempty ⟹
  ŵ ≥ 1/4 ⟹ |v| ≥ m/4 — §8.1 recovered in one line.

The remaining difficulty is exactly the light core's cheap cycles —
and gate C shows the cheapest of them are TRIVIAL (combs live at slab
weight ≤ 5). The floor argument must therefore be minimality- and
class-aware, never purely geometric:

- **Local reduction (proven, elementary).** A class-minimal logical v
  satisfies |v + z| ≥ |v| for EVERY stabilizer z, in particular every
  x-local one: no local generator's boundary may overlap v in more
  than half its support. This prunes the comb states (a comb overlaps
  its own tooth fully) and is sound for ANY catalog of local
  generators — an incomplete catalog only weakens, never falsifies,
  the resulting bound.

### §9.3 Gate S — quotient stacks: any transfer floor must be
### class-aware (`a40_s4_stack_gate.py`)

The y-deck transfer τ (stacking) of certified minima, all end-to-end
re-verified:

| stack | weight (rate) | nontrivial |
|---|---|---|
| (18,6) w12 → (18,12) ×2 | 24 (2/row) | YES — equals the certified d = 24; the b = 1 minimum IS a stack of the rate-2 base phase |
| (12,12) a36 w18 → (12,24) ×2 | 36 (1.5/row) | YES |
| (12,12) a36 w18 → (12,36) ×3 | 54 (1.5/row) | YES |
| (12,6) w12 → (12,12) ×2 | 24 (2/row) | YES — but the (12,12) witness at 18 < 24: at b = 0 the stack is not minimal |

Rate-1.5 NONTRIVIAL y-spanning logicals therefore exist at (12,12k):
**no frame-free rate-2 transfer bound holds on the family's frame
grid** — the floor is a property of (frame class + homology), not of
the transfer system. Cheap sub-rate-2 behavior must be classified and
killed per class: trivial recurrent behavior by local reduction
(§9.2), and periodic nontrivial behavior — "phases" — by the frame
analysis below.

### §9.4 The phase theory: sheared frames, winding, the b-bit

A y-periodic walk segment with period p and x-drift d per period is
exactly an X-cycle of the BB code on the quotient lattice
Z²/⟨(ℓ,0),(d,p)⟩. `a40_s4_phase_triage.py` (Smith normal form +
transported supports, cross-checked against the rectangular frames)
maps k over all shears for ℓ ∈ {12,18}, p ≤ 8: k ∈ {0,4,8,12}
depending on (p, d mod structure); k > 0 frames exist at EVERY p
(including p = 1, where weight < 2p = 2 is impossible anyway, and
p = 6, where d ∈ {0,6,12} gives the k = 12 rectangular-equivalent
frames). So rank alone kills nothing: phases must be excluded by
WEIGHT (censuses), which splits by winding character:

- **x-winding phases pay Θ(ℓ).** Gate W (`a40_s4_winding_gate.py`):
  the a36 witness pattern uses BOTH wraps (naive re-placement at
  (24,12) and (12,24) is not a cycle), so as an x-period-12 object it
  enters a frame only through x-stacking, at cost ∝ ℓ/12 per y-period
  — 2D-stacking arithmetic (18·(r/2)² = Θ(r²) ≫ 12r at (6r,6r))
  removes the witness species from every large member. Sharper, the
  **b-bit mechanism**: tiling a member (6(r+b), 6r) by the witness
  species needs 12 | 6r AND 12 | 6(r+b) — at b = 1 this forces r even
  and r odd simultaneously: IMPOSSIBLE at every b = 1 member (pure
  arithmetic, asserted to r < 200). The witness species — the −6
  discount carrier at (2,0) — never touches the b = 1 column.
- **x-compact phases are the real threat** — and they are now
  certified away below rate 2 for p ≤ 8 (§9.5).

### §9.5 The compact-phase atlas: exhaustive, both lanes, p ≤ 8
### (`a40_s4_phase_atlas.py`, `_ext.py`, `a40_s4_compact_triviality.py`)

Block-2 determinism makes every x-compact period-p object a
zero-to-zero walk of a finite x-automaton (state = 6 columns' worth of
F₂^p data; the free input is one column, the forced column solves the
unique most-advanced term — uniqueness asserted per pair). Exhaustive
BFS over (state, accumulated cost ≤ Wcap) with DAG path readout
enumerates ALL x-compact cycles of weight ≤ Wcap up to x-translation;
every readout is independently re-verified as a cycle on an embedding
torus. Lanes: pair (A,B) (the y-lane) and pair (B,Ā) (= the θ′-image,
θ′: (x,y) ↦ (y^{−1},x), θ′(A) = B, θ′(B) = Ā — the x-lane/strip
sector in rotated coordinates; the 90° rotation is NOT a code
equivalence, consistent with S1's exhaustive (12,18)↔(18,12) failure,
but it intertwines the two cycle systems exactly, boundaries to
boundaries).

> **Atlas verdict (certificate tier, exhaustive).** For BOTH lanes
> and every period p ≤ 8: NO nontrivial x-compact cycle of weight
> < 2p exists. Counts (identical between the lanes — the θ′ duality
> visible in data): 0/0/4/5/42/350/884 compact cycles at
> p = 2..8, Wcap = 2p−1, ALL trivial. Positive control: at p = 6,
> Wcap = 12, the engine finds the nontrivial minimum at exactly 12
> (66 objects — the L12 species); at Wcap = 14 the nontrivial
> spectrum is {12: 66, 14: 444} — the rate-2 floor is achieved
> exactly, is isolated (nothing at 13), and the next band sits at
> rate 7/3.
>
> **Compact triviality (closing ∀ℓ).** Every one of the 2×1285
> sub-rate-2 compact cycles has a COMPACT stabilizer generator s̃
> (unique, produced by the deterministic x-march on the generator
> equation, verified against both blocks, and re-verified as an
> H_X-row combination on an embedding torus; max s̃-extent 7). Hence
> each is trivial at EVERY torus (ℓ, p) that fits it — the atlas's
> verdict is ℓ-uniform, not an artifact of one embedding.

Consequence for the architecture: a class-minimal y-spanning logical
cheaper than rate 2 cannot spend its length in x-compact periodic
behavior of period ≤ 8 (no nontrivial such phase exists below rate 2,
and trivial ones are barred by local reduction), nor in x-winding
periodic behavior at large ℓ (Θ(ℓ) cost), nor in heavy slabs
(pointwise ≥ 2). What is NOT yet excluded: x-compact phases of period
> 8 (the pumping bound on periods is the light-core size, not 8), and
aperiodic light behavior — the wall/domain accounting. These are the
two named residual lemmas (§9.7).

### §9.6 THEOREM B6 — the (ℓ,6) row of the family, ∀ℓ

> **Theorem B6.** For every ℓ ≡ 0 (mod 6), ℓ ≥ 12:
> d((ℓ,6)) = 12, with the minimum achieved x-locally (the L12
> species). **Status: PROVEN — every finite certificate closed
> in-session** (`a40_s4_b6_close.py`, `s4_b6_close.json`, 23 s).
> - Upper: L12 places as a nontrivial weight-12 logical at every
>   (ℓ,6) — kernel-checked ℓ = 12..42 (gate W2), ∀ℓ by x-locality
>   (no x-wrap is used; the S2 §7.1 uniformity argument).
> - Lower, x-windowed branch (ALL ℓ at once): an (ℓ,6)-logical of
>   weight ≤ 11 with an x-gap ≥ 4 lifts, by the one-axis windowing
>   move of Lemma K, to an x-compact period-6 cylinder cycle of the
>   same weight; the atlas is exhaustive there and every ≤ 11
>   compact cycle is compactly trivial ⟹ the original was trivial.
> - Lower, x-spanning branch: weight ≥ ⌈ℓ/4⌉ ≥ 12 for ℓ ≥ 45
>   (Theorem L1); ℓ = 12 (gross, Lean-grade) and ℓ = 18
>   (census-complete, §3.2) certified.
> - Lower, ℓ ∈ {24, 30, 36, 42}: **full d ≥ 12 floors by
>   Z₂-descent, certificate tier, no SAT** (direct W = 11 windows
>   exceed the walk kernel's n ≤ 192 cap — `a40_s4_b6_frames.py`
>   records the abort; the descent route walks only n ≤ 192
>   frames). Towers: (24,6)→(12,6); (30,6)→(15,6) (n = 180
>   direct); (36,6)→(18,6)→(9,6); **(42,6)→(21,6)→(21,3), MIXED
>   axes** — the x-mid (21,6) is uncensusable at n = 252, but its
>   y-fold (21,3) has n = 126: the axis-generic machinery dissolves
>   what looked like the odd-rung/F1 case. Per frame: complete base
>   stab (μ ≥ 6, 7 orbit reps {6:1, 10:6} at every base) and
>   all-class ≤ 11 censuses (EMPTY ⟹ the seam lane is empty), the
>   b = 0 lane by τ-injectivity + 2·min(μ, d(L1)) ≥ 12, and 7
>   dangerous rungs at target 12 per frame (restricted lanes,
>   M ≤ 3, ALL PASS, covariance spot-checked). In-run
>   re-derivations agree with banked values: d((9,6)) = 10,
>   d((18,6)) ≥ 12, d((21,6)) ≥ 12 (new), d((21,3)) = 6 — the last
>   is 2p at p = 3: another exact rate-2 base row, whose ≤ 11
>   logicals all die on the way up the tower.

B6 closes §8.4 residue item 1 (the strip minima) at m = 6 in the
strongest possible form (∀ℓ, exact), and it is precisely the base
floor the stack analysis (§9.3) demands: the rate-2 phase row the
b = 1 minima stack from is now a theorem row, not a per-frame
observation.

### §9.7 The reformed lower-half architecture, and what remains

For a class-minimal nontrivial X-logical v of a member (Lemma K
dichotomy: y-spanning, else x-spanning and the θ′-mirror of the same
analysis applies):

1. Slab-amortize (§9.2): |v| = Σ ŵ_j; heavy steps pay ≥ 2.
2. Light steps live in the finite pruned light core (local reduction
   against a generator catalog — sound for any catalog).
3. Cheap light cycles = periodic phases: x-winding ones pay Θ(ℓ)
   (§9.4); x-compact ones of period ≤ 8 do not exist below rate 2
   (§9.5, atlas + compact triviality).
4. **[L-P, open]** x-compact phases of period 9..N₀ (N₀ = the
   light-core pumping bound): the strong-induction shape is "a cheap
   period-p phase is itself a y-spanning cheap object at a p-row
   frame", but the induction constants must be EXACT (a −C leak per
   level multiplies along stacks — the same exactness discipline the
   tower calculus already enforces); the atlas engine extends
   mechanically (p = 8 cost 85 s; cost grows ~2^p per period, so
   p ≤ 10−11 is the honest reach of the current Python engine).
5. **[L-W, open]** the wall/domain accounting: a cheap walk decomposes
   into phase domains and transition walls; the floor needs wall cost
   ≥ 0 at b = 1 (the b = 0 witness realizes a wall discount of
   exactly −6, so the accounting must be b-aware — the b-bit
   mechanism of §9.4 is expected to be the discriminator). The
   mechanical target is a Karp/potential certificate on the pruned
   light core (min mean = 2 with explicit potentials); its
   fixed-point structure (potential radius vs. the heavy threshold)
   is one scalar iteration, not a circularity.

With L-P and L-W, the assembly gives d(C_{r,1}) ≥ 12r − C with C from
the potential radius — Theorem UB(r,1) then pins d(C_{r,1}) ∈
[12r − C, 12r], and C = 0 on the pure-phase branch would close the
b = 1 conjecture exactly. This is the honest restatement of the L2/L3
ladder after this session: the compression lemma is no longer on it.

### §9.7.1 The periodic leg of L-W: decomposition, pilots, survivors

A closed walk's cheap recurrent behavior decomposes into periodic
phases (§9.4); for the wall lemma L-W the periodic leg is "no usable
y-spanning phase below rate 2". This session closed its structure:

> **Decomposition (exhaustive, elementary).** A y-spanning period-p
> phase at x-order ℓ (any shear d) either (i) has an x-gap ≥ 4 — then
> the one-axis windowing move lifts it to an x-COMPACT cylinder phase,
> the atlas's territory (§9.5: none below rate 2 for p ≤ 8, ∀ℓ); or
> (ii) has all x-gaps ≤ 3 — then it occupies ≥ ⌈ℓ/4⌉ columns, so its
> rate is ≥ ⌈ℓ/4⌉/p ≥ 2 whenever ℓ ≥ 8p − 3. **The whole periodic
> leg at p ≤ 8 therefore reduces to the FINITE frame list
> {(ℓ, p): 6 | ℓ, 12 ≤ ℓ ≤ 8p − 4}**: 11 frame families within the
> walk kernel's reach ((12, 2..8), (18, 3..5), (24, 4) — pilots
> below) and 22 needing the descent lane ((18, 6..8) through
> (60, 8); named residue).

**The pruning pilot** (`a40_s4_prune_pilot.py`): for each in-reach
frame family, complete node-exact censuses of ALL cycles ≤ 2p over
every class of every shear (stabilizers included), filtered to
y-spanning, then tested against the local-reduction catalog (single
H_X rows + row-pair sums ≤ 10: any cycle holding > half of a catalog
generator cannot be class-minimal). Survivors below rate 2 are
classified by triviality, x-winding character (re-placement at
doubled ℓ), and the b = 1 CLOSURE ARITHMETIC: a (p,d)-phase tiles a
member (6r+6, 6r) only if p | 6r and (6r/p)·d ≡ 0 mod (6r+6) — the
generalization of §9.4's b-bit mechanism to every shear.

**Pilot results** (`s4_prune_pilot_l{12,18,24}.json`; ℓ=12 all
p ≤ 7 all shears zero skips, 428 s; ℓ=18 p ≤ 5 (+ the counting-covered
p = 2, and (18, 6..7) recorded as descent-lane residue, one
no-info-set skip); ℓ=24 p ≤ 4 (three no-info-set skips at p = 4)):

- **ℓ=18: 14,400 y-spanning orbits, 3,132 below rate 2 — ALL
  PRUNED, zero survivors.  ℓ=24: 4,272 orbits, 2,160 below rate 2 —
  ALL PRUNED, zero survivors.** On the b = 1-relevant swept frames,
  class-minimality alone eliminates every sub-rate-2 periodic orbit.
- **ℓ=12: 936 survivors — every one NONTRIVIAL** — concentrated in
  exactly four shear families: (p,d) = (4,4) w6 (rate 1.5, 48),
  (5,7) w8 (rate 1.6, 120), (6,3) w10 (rate 5/3; 402 x-winding + 30
  twist-compact), (7,2) w12 (84), (7,10) w8 (rate 8/7, 84) + w12
  (168).  **NONE of the 936 closes around any b = 1 member
  (r ≤ 2000, closure arithmetic): the b = 1 walk system has no usable
  sub-rate-2 periodic phase anywhere in the swept range.**
- Correction logged (§9.8): the §9.4/§9.5 dichotomy "x-compact
  (atlas) vs x-winding (Θ(ℓ))" is NOT exhaustive on shear frames —
  30 of the (12,6,3) survivors are **twist-compact** (compact on the
  TWISTED cylinder Z²/⟨(d,p)⟩, equivalently straight-cylinder compact
  phases of a SHEARED polynomial pair): a new object class the
  compact atlas (twist 0 only) does not cover.  They too are
  b1-non-closers; their systematic classification (the atlas engine
  over the shear orbit of the pair) is named residue.
- Consistency controls: the tight rate-2 populations reproduce known
  censuses exactly — (12,6,0) and (12,6,6) each carry 1,884 tight
  survivors = the banked count of gross's weight-12 logicals; the
  (7,10) frame's 1,680 rate-2 objects sit beside its cheap family as
  the k-triage predicted.

The aperiodic half of L-W (domains-with-walls: finite stretches of
cheap phases glued by transition walls) remains the open core; the
periodic classification above fixes its input alphabet.

### §9.8 Falsified claims (session 4)

- **The §8.2 defect-graph formulation** ("nodes = interface states,
  d(wrapping) ≥ r·(min mean cycle) − C") — VACUOUS as stated:
  stabilizer combs put trivial y-spanning cycles at rate ≤ 3/4 in
  every such graph (gate C1). Compression could not have saved it.
- **The compression lemma as "the one missing ingredient" (§8.2,
  §8.4 item 2)** — re-scoped out of existence: finiteness is free
  below the heavy line, "x-periodic-reducible" is the phase
  classification, and the actual missing pieces (triviality- and
  class-awareness) were invisible to it.
- **"Small-p phase frames have k = 0"** (this session's working
  hypothesis before the triage) — FALSE: k > 0 shear frames exist at
  every p ≥ 1 (`s4_phase_triage.json`); phases die by weight, not
  rank.
- **Correction (2026-08-28, found by A42 S0 gate G5)**: the triage's
  original `tr()` transported supports into a frozenset, merging
  quotient-collided terms instead of cancelling them mod 2 (F2); the
  10 rows (ℓ,p,d) = (12,1,0),(12,1,9),(12,2,9),(12,3,1),(12,3,2),
  (18,1,0),(18,1,15),(18,2,15),(18,3,1),(18,3,2) were banked-k = 4
  ARTIFACTS with true BB-quotient k = 0.  All diffs conservative
  (banked 4 > true 0), so the phantom frames were merely censused
  unnecessarily — every census verdict stands.  Script made
  parity-aware and `s4_phase_triage.json` regenerated (survivor rows
  26 → 24); the "k > 0 at every p" bullet above survives on the
  corrected grid.
- **The 90° rotation θ as a code equivalence** — the one-bar
  obstruction (θ(A) = B̄, θ(B) = A; every axis-swapping monomial map
  leaves exactly one antipode unpaired) shows it is NOT one,
  explaining S1's exhaustive (12,18)↔(18,12) failure; it survives as
  an exact intertwiner of cycle systems (used for the x-lane).
- Session-internal: the first B6 frame script assumed the walk
  kernel takes n = 288..504 — it caps at 192 (3×64-bit words);
  caught by the kernel's own assertion before any census was
  claimed.
- **The compact/winding dichotomy as stated in §9.4–§9.5** — not
  exhaustive on shear frames: TWIST-COMPACT phases (compact on
  Z²/⟨(d,p)⟩) are a third class, found by the pilot at (12,6,3)
  (30 objects, w10).  The atlas's verdict stands for twist 0; the
  twisted classification is opened as residue.  (No claim consumed
  the incomplete dichotomy: the pilot that exposed it also verified
  the b1-closure kill for every survivor.)
- "The local catalog suffices on the periodic landscape" (the
  pilot's working hypothesis) — FALSE at ℓ = 12: four shear
  families survive it; they die instead by the closure arithmetic.
  At ℓ = 18/24 the catalog does suffice.
- POST-FIX ADDENDUM (S8-era re-audit, `a40_s4_p3audit.py`): after
  the tr() collision fix (66e0e7b), the pilot's p ≤ 3 rows were
  re-censused on the CORRECTED quotient codes at ℓ ∈ {12, 18, 24} —
  zero sub-rate-2 y-spanning orbits exist at ANY p ≤ 3 shear (3,456
  orbits; `s4_p3audit_l{12,18,24}.json`).  The collision-frame
  soundness gap is closed; no pilot-era conclusion changes.
- (Respected: witness weights as upper bounds only; RED/AMBER/GREEN
  as cost verdicts; no SAT anywhere; every gate script re-verifies
  its vectors end-to-end.)

### §9.9 Residue / S5

1. ~~The four B6 certificates~~ — CLOSED in-session
   (`a40_s4_b6_close.py`, §9.6): all four frames certified d ≥ 12 by
   Z₂-descent in 23 s, including (42,6) via the mixed-axis
   (21,6)/(21,3) tower. Theorem B6 carries no residue.
2. **L-W, the wall certificate**: build the pruned light-core graph
   (states = W ≤ 7 slabs up to x-translation, gap-saturated;
   transitions by the C2 recurrence; pruning by the generator
   catalog seeded from the atlas's 2×1285 compact trivial cycles)
   and compute the min-mean/potential certificate. The b = 0 −6 and
   the b = 1 closure arithmetic must fall out of the wall terms.
3. **L-P beyond p = 8**: push the atlas to p = 10−11 (engine as is)
   and formulate the exact-constant induction for general p; the
   (ℓ,12) row (= B12: d((ℓ,12)) = 24?) is the natural next base
   theorem, with (18,12) = 24 already certified as its anchor.  A
   module-theoretic route (compact-support H₁ of the Laurent-in-x /
   Artinian-in-y complex) could give ∀p in one shot.
3b. **The twisted-compact atlas** (the §9.7.1 correction): run the
   compact engine over the shear orbit of the pair (transformed
   supports per twist class); finite per period.  Also the periodic-
   leg descent list (22 frame families, (18,6) … (60,8)) and the
   p = 8 @ ℓ = 12 layer (W = 16 cost wall — the one swept-range gap,
   priced ~hours with today's kernel).
4. **The x-lane run of the same program** (pair (B,Ā)): atlas done;
   the mirrored comb/stack/triage gates are cheap and owed for
   symmetry of the record.
5. **Lean**: B6 is a natural target next to L1 (the atlas is a
   finite certificate; the windowing lift and L12 uniformity are the
   same species as UB(r,1)'s arguments); the compact-generator
   marches are `decide`-shaped.
6. Corpus notes: the (12,16)/(12,24) rate-3/4 combs and the p ≤ 8
   phase-frame k-map are reusable structure for any BB pair with
   spans (4,4); the atlas engine is pair-generic.

## §10 SESSION 5 (2026-08-27) — the periodic leg closed at ℓ = 18/24,
## two chiral species, and the drift-blind wall certificate refuted

Directive: advance the L-W wall certificate — build the ℓ = 18 pruned
light-core transfer graph, extract its bi-recurrent part, and compute
the min-mean/potential certificate.  What actually happened,
falsify-first: the mandated pricing REFUTED the buildability of the
literal Stage-1 graph (a structural obstruction, §10.1), and the
periodic content the graph was meant to certify was obtained instead
by a cheaper, stronger instrument — fiber-complete Z₂-descent censuses
of every sheared frame (18/24, p ≤ 8), which CLOSED the periodic leg
on both columns and discovered that the hoped-for certificate target
is FALSE as stated: the pruned light core at ℓ = 18 contains a
nontrivial y-spanning cycle at rate 8/7, and no constant-C rate-2
statement can hold uniformly on the (18,·) column (explicit sub-2m
logicals at (18,63)/(18,36)/(24,48), end-to-end).  What survives — and
is now sharply supported — is the closure-constrained (member) form:
both cheap species are CHIRAL, closure-dead on the entire b = 1 column
by pure arithmetic, and the wall lemma becomes a momentum-budget
statement (§10.7).  Scripts `a40_s5_*.py`, data `data/a40/s5_*`;
`validate_banked()` green before every stage.

### §10.1 Stage-1 pricing: the light-core graph is not materializable
### (`a40_s5_lightcore.py`, `s5_lightcore.json`)

The mandated RED/AMBER/GREEN gate, with the counts computed exactly
and the structure lemma verified mechanically:

- Universe: 4+4-row windows (the dynamics+cost state) number
  Σ_{w≤7} C(8ℓ, w) = 2.3e11 at ℓ = 18 (1.3e10 per x-translation);
  the 5+5-row window (the smallest that can SEE a tooth — a single
  H_X generator spans 5 consecutive rows: block-1 cy−1..cy+1,
  block-2 cy..cy+3) is bounded by 2.8e20 (~6e13 after translation
  and the internal constraint).  RED at every fidelity.
- **Constraint lemma (mechanical, 200/200).** A 4-row window
  carries NO internal E_j constraint: every random light 4+4 window
  extends to an admissible bi-infinite walk (forced-v₂ forward,
  forced-v₁ backward — the x·v₁[j−3] pivot is monomial, hence
  solvable backward).  A 5-row window carries exactly ONE (E_{j−1}).
  Consequence: reachability prunes (essentially) nothing; the only
  cut is BI-RECURRENCE, a greatest fixed point that requires
  materializing the universe.  There is no seed set: light cycles
  need not contain any enumerable anchor (dense cycles have no
  3-gap; the y-spanning states are the whole universe).
- The charter's cluster reduction is exact at ℓ = 18 (no gap can
  saturate: the light-bridge bound 28 exceeds ℓ) — but it only
  divides the count by ℓ.  The reduction that would make the space
  small is a PER-CLUSTER POTENTIAL DECOMPOSITION (Φ = Σ_clusters φ
  + gap terms, verified on the single-cluster transition system) —
  named as the S6+ route, not executed.
- Verdict: RED; and §10.3–§10.5 show the plain certificate the
  graph was to compute is FALSE anyway — the pricing gate did its
  job twice over.

The session's engine deliverable is instead the exact WINDOW
MACHINERY the future certificate needs, fully validated: the
drift-periodic unroller (a (p,d)-phase vector on the SNF-normalized
quotient torus → walk rows with row(y+p) = x^d row(y), re-verified
against E_j for every j), slab weights with the exact telescope
Σ W_j = 4|v|, and the H = 5 window-prune rule (≥ 4 of 6 cells of a
fully visible tooth ⟹ not class-minimal).  Controls: combs are
window-pruned (1 event); all 936 banked ℓ = 12 pilot survivors
re-verify and are window-UNpruned (consistent with their global
verdicts); the a36 witness re-verifies at rate 1.5.

### §10.2 The fiber-complete descent census machine
### (`a40_s5_price.py`, `a40_s5_dense_close.py`)

The periodic-leg frames (ℓ, p, d) with n = 36p resp. 48p > 192 are
censused COMPLETELY (weight ≤ 2p−1, every class, stabilizers
included) through an even-axis Z₂ fold: every cover cycle v with
|v| ≤ W has base shadow b = P(v) with |v| = |b| + 2·overflow, so the
union of (i) the τ-lane over all base cycles u with 2|u| ≤ W (b = 0:
v = τ(u), free-deck elementary, τ-injectivity asserted per deck) and
(ii) `enumerate_lifts_deep` fibers with cap (W−|b|)//2 over base
translation-orbit reps is the full census up to cover translation
(fold-equivariance: cover translations surject onto base
translations; covariance spot-checked per frame).  For (24,8) a
second fold (depth 2, bottom n = 96) repeats the pattern one level
down.  No SAT anywhere; every produced vector re-verified
(`is_cycle`, weight, class).

Controls, all PASS:
- **A (engine).** (18,5,1): the descent census, translation-expanded,
  equals the direct kernel census EXACTLY (90 elements).
- **B (banked).** (18,6,0) ≅ rectangular (18,6): nontrivial ≤ 11
  EMPTY, matching d((18,6)) = 12 (census-complete, §3.2); trivial
  reps {6:1, 10:6} = the B6 base pattern.  Likewise (24,6,0) EMPTY
  matches Theorem B6 at (24,6).
- **C (depth-2).** (18,6,3) run depth-2 ((3,36)→(3,18)→(3,9))
  reproduces the depth-1 census exactly (counts and survivor orbit).
- Ops: the naive per-class node bill overprices by ~150× —
  `census_pass` walks one shared tree per ≤ 51 offsets (the pricing
  script's calibration row); all 72 frames of §10.3 cost 11 min
  total where the naive bill said 7+ hours.

### §10.3 THEOREM P18/P24 — the periodic leg at ℓ ∈ {18, 24}, p ≤ 8

> **Theorem P (periodic-leg closure; certificate tier, no SAT).**
> Over ALL sheared frames Z²/⟨(ℓ,0),(d,p)⟩ with ℓ ∈ {18, 24},
> 2 ≤ p ≤ 8, d ∈ Z_ℓ, the complete list of NONTRIVIAL X-cycles of
> weight < 2p (= sub-rate-2 periodic phases of the y-walk system) is:
>
> | frame | species | weights (orbit reps) | catalog-unpruned |
> |---|---|---|---|
> | (18,6,3) | twist-compact TC63 | 10 ×1 | 1 |
> | (18,7,16) | winding W7 (d ≡ −2) | 8 ×1, 10 ×4, 12 ×40 | 3 |
> | (24,6,3) | TC63 | 10 ×1 | 1 |
> | (24,7,22) | W7 (d ≡ −2) | 8 ×1, 10 ×4, 12 ×40 | 3 |
>
> Every other (ℓ, p ≤ 8, d) frame is EMPTY (k = 0 frames by rank;
> k > 0 frames by census: p = 8 empty at BOTH ℓ, all 28 frames).
> Combined with the S4 pilot (p ≤ 5 at 18; p ≤ 4 at 24 — the three
> (24,4,·) skips are k = 0 frames, nontrivial-empty by rank) and the
> compact atlas, this closes §9.7.1's descent-lane residue for the
> (18,·) and (24,·) families through p = 8, including (24,5).
>
> **Trivial side (the pruning-catalog sufficiency, extended).** At
> every frame swept, EVERY trivial y-spanning sub-rate-2 orbit
> (~4,000 orbit reps per ℓ across p = 6..8, combs included) is pruned
> by the local-reduction catalog (single H_X rows + pair sums ≤ 10).
> The S4 claim "the catalog suffices on the ℓ = 18/24 periodic
> landscape" now holds through p = 8.

The exact spectrum replication between ℓ = 18 and ℓ = 24 — same
weights, same orbit-rep counts, same pruned/unpruned split — is a
measured fact about both species (§10.4).

**ℓ = 30 extension (same session, same machine).**  The p ∈ {5, 6, 8}
layers of the (30,·) column (p = 7 is the F1 odd-fold gap: (30,7)
folds only to n = 210 > 192 and (1,105) has no even axis): p = 5 all
20 k>0 frames EMPTY; p = 6: exactly TC63 at (30,6,3) (w10 ×1,
unpruned, closure-dead — the twisted-cylinder descent prediction
d ≡ 3 (mod ℓ) confirmed at a FOURTH column); p = 8: all 20 k>0
frames EMPTY (depth-2, bottom n = 120, 264 s).  Trivial side: all
pruned, every frame swept.

### §10.4 The two species (banked objects)

**TC63 — the twist-compact (6,3) family.**  Weight 10, rate 5/3,
y-period 6, x-drift +3; slab profile [5,9] (NOT light: it crosses the
heavy boundary); catalog- and window-unpruned; present at ℓ = 12
(the pilot's 30), 18, 24 — it lives on the lattice ⟨(3,6)⟩ and
descends exactly to the frames with d ≡ 3 (mod ℓ), which is where
(and only where) the censuses find it.  The TWISTED-COMPACT ATLAS
(`a40_s5_twisted_atlas.py`, §9.9 item 3b executed): with
g = gcd(t,p), the twisted cylinder Z²/⟨(t,p)⟩ ≅ Z × Z_g by a
unimodular change of basis, so each twist LATTICE reduces to the S4
straight automaton on a transformed pair at period g.  SCOPE: the
sweep covers the lattice window 0 ≤ t < p (NOT "t mod p" — the
y ↦ y·x^k change maps ⟨(t,p)⟩ to ⟨(t+kp,p)⟩ only by TRANSFORMING
the pair, so ⟨(t,p)⟩ and ⟨(t+p,p)⟩ are genuinely different systems
for the fixed pair, and a frame (ℓ,p,d) receives descent from every
⟨(d+ℓk, p)⟩; the frame-level truth is protected by Theorem P's
censuses regardless).  Verdict (exhaustive on the window, both
lanes, p ≤ 8): **⟨(3,6)⟩ is the ONLY twist lattice with |t| < p
carrying any sub-rate-2 nontrivial compact object** (min weight 10,
found in both lanes — the θ′ duality again — and re-verified
end-to-end on the (42,6,3) torus); (4,2), (6,2), (6,4), (8,2),
(8,4), (8,6) are empty ((6,4)-BAbar has a genuinely unsolvable
column march — tied top terms with singular sum — and is closed by
its θ′-image (4,2)-AB at the transported window Wcap 11: empty).
g = 1 classes are LINE systems: over F₂[u,u⁻¹] the cycle module is
(P̄_u, Q̄_u)/ĝ · R and the boundary module (P̄_u, Q̄_u) · R, so
H₁ ≅ R/(ĝ): nontrivial line-compacts exist iff ĝ = gcd(P̄_u, Q̄_u)
is non-monomial.  Measured: exactly 12 of 21 directions are
non-monomial — ALL through the same factor u² + u + 1 — and every
generator sits at rate ≥ 3/row: **no sub-rate-2 line species**
(cross-validated: a cheap line element of x-extent < ℓ would appear
in the §10.3 censuses; none does.  The unbounded-extent corner —
|h·gen| < |gen| cancellations — is thin formal residue).

**W7 — the winding (7,−2) family.**  Weights {8, 10, 12} (rates
8/7..12/7), y-period 7, x-drift ≡ −2 (mod ℓ); present at ℓ = 12
((12,7,10), the pilot), 18, 24 with IDENTICAL spectrum and weights —
**the "x-winding pays Θ(ℓ)" heuristic is false in general** (§9.4's
argument was and remains specific to the witness species' stacking
arithmetic).  The w8 member is ALL-LIGHT (slabs [2,7]) and unpruned:
**μ_light(ℓ=18) ≤ 8/7 and μ_light(ℓ=12) ≤ 8/7** — the pruned strict
light core is NOT cycle-free at either ℓ.  Not twist-compact (the
(7,2)-line class has rate-≥3 floors), so it genuinely uses the x^ℓ
relation at each ℓ while keeping constant weight; note its frames
(ℓ,7,ℓ−2) are CYCLIC groups Z_{7ℓ} — a univariate/cyclotomic-gcd
mechanism is the natural suspect (open; the (30,7,28) prediction —
same spectrum — is stated for a future odd-fold session: (30,7)
folds to n = 210 > 192 and (1,105) has no even axis, the F1 case).

**Chirality (measured, new).**  Both species are ONE-SIGNED in
drift: (18,7,2), (24,7,2) (the +2 mirror) and (18,6,15), (24,6,21)
(the −3 mirror) are all EMPTY.  Cheap drift exists in one direction
only — this is what gives the closure arithmetic its teeth (§10.7).

### §10.5 Stack certificates: the uniform column statement is FALSE;
### the member column is arithmetically protected

Each species closes exactly on the tori its closure arithmetic
permits, and the stacks are genuine (end-to-end torus verification,
weight/cycle/class):

- W7 ×9 at (18,63): **d((18,63)) ≤ 72 < 126 = 2m**;
- TC63 ×6 at (18,36): **d((18,36)) ≤ 60 < 72 = 2m**;
- TC63 ×8 at (24,48): **d((24,48)) ≤ 80 < 96 = 2m**.

(UB tier only — witness weights are upper bounds.)  Hence **no
constant C makes "d((18,m)) ≥ 2m − C for all y-spanning-sector m"
true** — at m = 63k the deficit grows linearly.  The charter's
Stage-1 target ("μ = 2 with finite Φ ⟹ d ≥ 2m − C for EVERY (18,m)")
is refuted as a statement, not merely unreached.

> **Member protection (arithmetic, ∀r).**  A (p,d)-phase tiles the
> b = 1 member (6(r+1), 6r) only if p | 6r and
> (6r/p)·d ≡ 0 mod 6(r+1).  For W7 (p = 7, d ≡ −2): 7 | r and
> (7s+1) | 2s — impossible for every s ≥ 1.  For TC63 (p = 6,
> d = 3): 2(r+1) | r — impossible for every r.  **Neither species
> touches any b = 1 member, at any r** — the b-bit mechanism of
> §9.4, now carried by the only two sub-rate-2 species in existence
> at p ≤ 8.

### §10.6 Light-core measurements (`s5_lightcore.json`)

- The a36 witness (rate 1.5 at (12,12)) has slabs [2,12]: it is NOT
  a light-core cycle — its cheapness lives partly in the heavy
  region.  The charter control "μ_pruned(ℓ=12) = 1.5 via the
  witness" is superseded: the true light-core bound is ≤ 8/7 (W7's
  w8), at both ℓ = 12 and 18.
- The L12 species at (18,6) (rate 2): slab range [5,11] across the
  66-object family — even the tight species oscillates through the
  heavy boundary; the future potential certificate must price
  boundary crossings, not avoid them.
- Combs: window-pruned within one tooth alignment, as designed.

### §10.7 The sharpened L-W target (what survives, exactly)

For the b = 1 member (ℓ, m) = (6r+6, 6r), the periodic input
alphabet at p ≤ 8 is now: rate-≥2 phases, plus two closure-dead
chiral species — proven at the censused columns ℓ ∈ {12, 18, 24}
(members r ≤ 3), and conjectured ∀ℓ on the strength of the exact
ℓ-replication (the ∀ℓ alphabet needs either per-ℓ censuses, which
the descent machine makes ~minutes per column while even-axis folds
reach n ≤ 192, or the L-P induction).  The correct wall statement
is a MOMENTUM BUDGET:

- Any drift-blind window-potential certificate (ŵ + ΔΦ ≥ 2 on the
  pruned light core) is impossible — W7's w8 cycle is light,
  unpruned, and has mean 8/7; it also closes on the ℓ-drift cover
  (9 periods, drift −18 ≡ 0), so grading by drift alone cannot
  save rate 2 either (the (18,63) logical is real).  The grading
  that works must see the y-length: the species needs 7 | m
  (and 63 | m with the drift), which no member m = 6r with
  ℓ = 6r+6 satisfies compatibly.
- Chirality makes this quantitative: at a member, a W7 domain of
  T rows banks deficit (6/7)·T and drift −2T/7, and NO cheap
  structure with opposite drift exists (the +2 frames are empty) —
  the return drift must be paid at non-species rates.  The wall
  lemma's member form: species domains cannot wrap (closure
  arithmetic), so each domain's deficit + return-drift cost is
  chargeable against its walls and complement.  This is the open
  aperiodic half, unchanged in kind but with its alphabet, deficits
  (6/7 and 1/3 per row), and the chirality lever now pinned by
  data.

### §10.8 Falsified claims (session 5)

- **The charter's Stage-1 certificate target** ("min mean 2 with
  potential radius C on the pruned ℓ=18 light core ⟹ d ≥ 2m − C
  for every (18,m)") — FALSE: μ_light(18) ≤ 8/7 (explicit light
  unpruned nontrivial cycle), and the (18,63k) stacks kill every
  constant C.  The salvageable statement is the member/closure-
  constrained form (§10.7).
- **"μ_pruned(ℓ=12) = 1.5, realized by the a36 witness"** (charter
  control) — the witness is not light (slab 12); the light-core
  minimum is ≤ 8/7 via W7's w8 member (which IS all-light).
- **"x-winding phases pay Θ(ℓ)"** as a general mechanism — false:
  W7 has weight 8 at ℓ = 12, 18, 24.  (§9.4's claim survives only
  as proven: for the witness species via stacking arithmetic.)
- **The S4 `x_winds` classifier is representative-dependent** (found
  by the depth-2 control: the same (18,6,3) orbit reports
  False/True across translates).  Its re-placement lifts one
  section per point; a False verdict (a closing section exhibited)
  is sound, a True verdict is not.  The S4 (12,6,3) split "402
  x-winding + 30 twist-compact" degrades to "≥ 30 twist-compact-
  representable; winding count an upper bound"; no banked claim
  consumed a True verdict alone (b1-closure is arithmetic; the
  witness-species cost is the b-bit argument).  The twisted atlas
  is the sound classifier.
- **The naive census pricing** (nodes × classes) — overprices ~150×;
  `census_pass` shares one walk per ≤ 51 offsets.  Pricing scripts
  must calibrate against a real chunked run.
- Session-internal: the first universe-pricing draft claimed "no
  internal constraint in ANY window" — corrected mechanically: 4-row
  windows none, 5-row windows exactly one (E_{j−1}).
- (Respected: witness weights as upper bounds; RED/AMBER/GREEN as
  cost verdicts; no SAT on certificate paths; every consumed vector
  re-verified end-to-end.)

### §10.9 Residue / S6

1. **The member wall lemma (momentum-budget form, §10.7)** — the
   aperiodic half of L-W: charge species-domain deficits + one-signed
   return drift against walls at m = 6r, ℓ = 6r+6.  The first
   mechanical step: the wall/interface cost table between the three
   phase families (rate-2 tight, TC63, W7) under the window engine.
2. **p ≥ 9 periodic leg**: dense censuses priced RED (κ ≥ 77,
   W ≥ 17: ~days/frame); the straight atlas at p = 9 (Wcap 17) was
   run and KILLED at ~8 GB RSS pre-cap — the pure-Python BFS+DAG
   engine's memory envelope ends at p = 8 (measured; §9.9's "honest
   reach ~p ≤ 10" was optimistic — a compiled or disk-backed engine
   is needed; `s5_atlas_p9.json` records the abort).  The L-P
   induction (exact constants) is still the ∀p route.
3. **The W7 mechanism — now half-pinned (measured in-session).**
   The (ℓ,7,d) k = 4 frames are cyclic (Z_{7ℓ}) and their univariate
   images share the factor z² + z + 1 (the F₄ ω-point; it divides
   z^{7ℓ}−1 since 3 | 7ℓ for every 6 | ℓ) at EVERY k = 4 shear, at
   ℓ = 18 and 24 — the ω-eigenspace IS the k = 4 (2 blocks × deg 2),
   and the SAME ĝ = u²+u+1 runs through all 12 non-monomial line
   directions: one ω-point of the pair unifies the line H₁, the
   (ℓ,7,·) ranks, and (plausibly) the two-sided thread's Frobenius
   motif.  What remains open is the CHIRALITY selection: only
   d ≡ −2 (mod ℓ) has weight-8 representatives in its ω-class — a
   cyclic-code minimum-distance question on Z_{7ℓ} per shear (the
   sharpest next question).  The (30,7,28) same-spectrum prediction
   stands (F1 odd-fold territory for the census; the ω-rank part is
   already checkable univariately).
4. **ℓ ∈ {30, 36, ...} columns**: the descent machinery covers any
   ℓ with even-axis folds reaching n ≤ 192; (30,7) does not fold —
   F1 or a bigger kernel.
5. **The per-cluster potential decomposition** (the only identified
   route to a materializable light-core certificate; §10.1).
6. **Small formal residues**: the (6,4)-BAbar direct march (tied
   pivots — needs a 2-column block-march; currently closed via θ′);
   the line-species unbounded-extent corner; window-realization of
   the trivial-side pruning for p = 6..8 (verified global-catalog-
   wise; the H = 5 window realization is verified for combs only);
   the twist lattices with |t| ≥ p (transformed-pair atlas runs —
   finite per (p, ℓ-window); the swept frames are already protected
   by Theorem P, so this matters only for un-censused ℓ).
7. **Lean**: Theorem P's census certificates are decide-shaped
   (finite lists + lift arithmetic); the member-protection
   arithmetic ((7s+1) ∤ 2s, 2(r+1) ∤ r) is a two-line Lean lemma —
   natural first S6 formalization targets alongside B6.

## §11 SESSION 6 (2026-08-27) — the momentum budget instrumented:
## drift defined, the light frontier measured, ghost passages
## discovered, and the boundary-coupling wall named

Directive: advance L-W in the §10.7 momentum-budget form — define
drift soundly, measure the deficit-vs-drift Pareto frontier of light
pruned walk segments by weight-capped BnB, assemble a member floor,
and cash it out at (ℓ,m) = (24,18) = [[864,12]].  What happened: the
drift definition and its species verification are DONE and banked;
the frontier engine is built and validated (never materializing the
light core); the measured frontier REFUTED two of the three charter
hypotheses as stated — not by species geometry but by a genuinely
new object this session verified end-to-end: GHOST PASSAGES, cheap
near-vacuum light interludes admitted by the run-boundary
relaxation.  The assembled floors are correspondingly modest and
honestly scoped; the binding obstruction is now a named, concrete
lemma (the boundary coupling), not a vague "aperiodic half".
Scripts `a40_s6_*.py`, data `data/a40/s6_*`; `validate_banked`
green before every stage.

### §11.1 The drift definition (banked) and its verification
### (`a40_s6_drift.py`, `s6_drift.json`)

**Definition.**  Work on the universal cover of the x-circle (rows =
subsets of Z × {blk}, no wrap; the recurrence E_j wrap-free).  For a
fragment on rows [t0, t1], the per-slab anchor A_j := min occupied
column of rows [j−3, j]; the fragment drift
δ := A_{t1} − A_{t0+3} ∈ Z.  Additivity is definitional: fragments
glued on a shared 4-row window telescope, δ(AB) = δ(A) + δ(B)
(verified mechanically).  On drift-periodic loops (row(y+p) =
row(y) + s on the cover) the per-period anchor difference equals s
for EVERY anchor convention — loop drift is gauge-free (verified:
min- and max-anchor agree per period).  A torus phase's INTEGER
drift is determined mechanically by cover-lift trial: exactly one
s ≡ d (mod ℓ) admits a wrap-free lift (no false accepts possible —
a straddling rep assignment breaks E), and a phase with NO
admissible lift for any s genuinely winds x.

**Species verification (the charter's demand, exact).**
| object | frame | cover drift | slabs (1 period) | light |
|---|---|---|---|---|
| W7 w8 | (18,7,16) | **−2 per 7 rows** | [7,6,5,4,2,3,5] min 2 | ALL |
| W7 w8 | (24,7,22) | **−2 per 7 rows** | same profile | ALL |
| TC63 w10 | (18,6,3) | **+3 per 6 rows** | [9,6,6,6,5,8] max 9 | no |
| TC63 w10 | (24,6,3) | **+3 per 6 rows** | same | no |

3-period extents: W7 10–12 (embeds at every ℓ ≥ 16); TC63 17–18
(at ℓ = 18 a ≥ 3-period TC63 run already wrap-interacts; at ℓ = 24
runs of ≤ 4 periods are cover-visible).  Window-pruned: none.
Deficits per period: 6 (W7), 2 (TC63).

**Winding detection (new, load-bearing for control (c)).**  The a36
witness at (12,12) has NO compact cover lift for any s ∈ {0, ±12} —
mechanical confirmation of gate W's "uses both wraps".  All 48
all-light (4,4)-w6 pilot survivors at ℓ = 12 are the same: the
witness species IS the (4,4) family (3 periods = the w18 witness),
and its −6 discount lives in the WRAPPED/winding corner, invisible
to cover fragments.  The ℓ = 12 W7-twin (7,10)-w8 lifts at −2.

### §11.2 The connectivity lemma (the teleport killer)

E_j evaluated at column c touches exactly v1[j]{c, c+1},
v1[j−3]{c−1}, v2[j]{c}, v2[j+1]{c}, v2[j−1]{c+3}: footprint spans
(Δx, Δy) ≤ (4, 4).  Supports separated beyond one footprint SPLIT:
E(v) = E(S1) ⊕ E(S2) on disjoint constraint sets, so each component
of a closed cycle is a cycle (verified on fragments and on an
embedding torus).  A class-minimal nontrivial logical therefore has
(4,4)-CONNECTED support: trivial components could be subtracted;
several nontrivial components would each be lighter logicals.
Consequence: "teleport" adversaries (far cheap flashers faking
drift at O(1) cost) do not exist in minimal logicals, and the
per-fragment drift bookkeeping is anchored to a single connected
object.

### §11.3 The frontier engine (`a40_s6_frontier.py`,
### `s6_frontier_u*.json`)

The mandated enumeration lane, built and validated:

- **Min-window stratification** (the pricing move that beats the
  seed wall): every fragment splits at its first minimum-weight
  4-row window u.  For u ≥ 5, every slab ≥ 5 gives the ANALYTIC
  bound D ≤ (2 − u/4)h ≤ 0.75h < (6/7)h — below the species rate,
  no enumeration needed.  For u ≤ 4 the split window is a connected
  ≤ 4-point full-content window (~4e5 seeds vs ~2e8 for u ≤ 7 —
  the 4-row shadow of §10.1's materialization wall, avoided).
- **Two marches** from the split window: forward (v2-forced, free
  v1 inputs) and backward (v1-forced via the monomial pivot, free
  v2 inputs), each a cost-ordered Dijkstra on the slab-sum g with
  dominance on the normalized 7/8-row dynamic state, per-slab
  W ∈ [u, 7], the H = 5 window rule checked on the forced row
  BEFORE input branching, (4,4)-local growth, y-spanning by
  construction.  Tables minG(h, δ) are COMPLETE below the cost cap
  (each piece of a split costs at most the whole).
- **Composition**: loose join (min over seeds independently per
  side) can only overstate deficits — sound for the assembly;
  MATCHED (per-seed) join at u = 1 (8 seeds) measures the
  looseness: 253 buckets differ, by up to 3.25 deficit.
- **Universal envelope** (unconditional, all h): slabs ≥ 1 give
  g ≥ h, so D(h, δ) ≤ 1.75h always.
- **Validation**: every species sub-fragment whose min slab lands
  in the stratum is dominated by the tables (57/57 at u = 2); the
  drift script's fragments re-verify through the engine's
  conventions.
- **Pricing verdicts** (measured): u = 1, 2 exact strata GREEN
  (18M nodes / 213 s at u = 2, g ≤ 24); u = 3, 4 exact strata RED
  in this engine — the min-slab floor forces multi-point rows,
  branching ~120/node (24M pushes at g ≤ 8) — their analytic
  bounds (1.25h, h) are dominated by the u ≤ 2 tables anyway, so
  the master table loses nothing.

### §11.4 GHOST PASSAGES — the discovery, end-to-end verified
### (`a40_s6_ghost_verify.py`, `s6_ghost_verify.json`)

The frontier is NOT species-dominated in the transient regime.  The
run-boundary relaxation (E_{t0}..E_{t0+2} unenforced — exactly
right for a maximal light run entered from a heavy slab) admits
**ghost passages**: near-vacuum stretches sustaining slab weights
1–3 for many slabs.  Mechanism: with zero v1, the forward recurrence
is v2[t+1] = v2[t] + x^{−3} v2[t−1] — a Sierpinski/(1+u)^t weight
law 2^{popcount(t)} that self-limits forward coasts to ~8 slabs;
the BACKWARD march's free v2 inputs run the recursion in reverse
(each step can cancel), coasting much longer.  Verified specimen
(reconstructed from the search, independently checked through
CoverFragment: E_j, slabs, window rule, anchors): rows [−9, 3],
weight 6, TEN slabs [2,2,4,2,2,3,1,1,1,1], drift +5, unpruned.
A run entered from a heavy boundary can genuinely fade to weight
~1/slab and revive.

Matched-join frontier at u = 1 (real fragments only, g ≤ 30):
D grows ~1.5/slab through h ≈ 11, then BENDS — marginal deficits
+1.0, +0.75, +0.75, +1.0, +0.75 over h = 12..16 — and the
certified transient D(h) − (6/7)h SATURATES at 7.25–7.5 across
h = 12..16 (T0 = 7.50, also realized by a u = 2 bucket at
(14, +6)): within the certified range the frontier reads exactly
as "species rate + a bounded boundary transient",
D(h) ≈ (6/7)h + 7.4.

**Hypothesis verdicts (charter Stage 1):**
- (i) "D > 0 only on one drift sign" — **REFUTED**: ghost deficits
  carry BOTH signs (the verified specimen drifts +5; best-D buckets
  appear at δ ∈ [−6, +6] and beyond).  Chirality is a property of
  the sustained species, not of the transient frontier.
- (ii) D ≤ a·|δ| + C0 — reshaped: the drift-extreme region carries
  ghost deficit too (D ≈ 21 at |δ| ≥ 16 within caps); a winding
  kill must come from the closed-branch buckets, not a per-|δ|
  slope.
- (iii) "the frontier is achieved by species runs" — **split
  verdict**: transient regime (h ≲ 12) ghost-dominated (species
  beaten by ~2×); sustained regime consistent with species
  extremality (the matched bend drops the marginal rate through
  6/7 by h = 14), but the caps end before an asymptotic claim.

### §11.5 Stability and scope of the tables

- smax (new points/row) 3 → 4 at g ≤ 16: IDENTICAL tables — the
  cap is not binding.
- Growth radius dil 4 → 6 at g ≤ 16: 218 buckets get cheaper, 65
  new (mostly drift-extremes, some near δ = 0 by 0.5 deficit) —
  the (4,4)-adjacent growth restriction IS binding; production
  tables are scoped "radius-4 prefix-connected growth" and every
  floor below carries that condition.  (Fully sound closure of the
  scope = chained same-row placements + connect-later pioneers —
  converging strands — both named residue; dil-6 production priced
  ~7× = RED today.)
- Wrapped corner: fragments of torus extent ≥ ℓ−3 that
  wrap-interact without winding are NOT enumerated (the ℓ = 12
  witness species is the mechanical proof such objects exist and
  carry deficit).  At ℓ = 24 this is a listed conditionality;
  winding-by-±24 IS covered (extent cap 34).

### §11.6 The assembly and the (24,18) floors
### (`a40_s6_assemble.py`, `s6_assembly.json`)

**The assembly (y-spanning sector).**  For a class-minimal
nontrivial y-spanning X-logical v of (24,18): v is connected
(§11.2); 4|v| = Σ_j W_j over the 18 slabs (exact telescope); heavy
slabs (W ≥ 8) pay ≥ 2 pointwise; maximal light runs are boundary-
relaxed light fragments, deficit-bounded by the frontier tables;
the walk partitions cyclically into r ≥ 0 runs and ≥ r heavy
blocks.  Drift: anchor slips across heavy slabs are NOT weight-
bounded by any lemma proven this session (a slip can ride branches
bridged elsewhere in the connected walk), so the MIXED branch is
assembled WITHOUT the closure lever — the drift/momentum constraint
survives only on the all-light branch (r = 0 heavies: the walk is
one 19-slab fragment with equal end windows and total cover drift
≡ 0 mod 24, winding = ±24 included).

**Floors (member (24,18) = [[864,12]], y-sector):**
- MIXED branch: max Σ D over run partitions (heavy slips free)
  = 26.25 at one 17-slab run ⟹ **w ≥ 10**.
- ALL-LIGHT branch: per min-window stratum: u = 1 (matched,
  g ≤ 30): (19, δ ≡ 0 mod 24) buckets empty ⟹ w ≥ 8; u = 2
  (loose, g ≤ 28): ⟹ w ≥ 7; u ≥ 5 analytic ⟹ w ≥ 23 ⟹ **w ≥ 7**.
- **T1 (certificate-shaped, scope-listed): d_Y(24,18) ≥ 7**, the
  binding branch being the all-light u = 2 cost cap.  Scope
  conditions: radius-4 prefix-connected growth (§11.5), no
  wrap-interacting non-winding fragments, loose join at u = 2,
  smax 3 (measured non-binding), boundary-relaxed run semantics
  (exact for runs, the source of ghosts).
- **T2 (conjectural: D(h, δ) ≤ min(1.75h, (6/7)h + T0) with
  T0 = 7.50 — the measured saturation, extrapolated beyond caps):
  mixed ≥ 8, all-light ≥ 21 (sustained species rate) ⟹
  d_Y(24,18) ≥ 8**; the T2 gap to 2m = 36 is dominated by the
  per-boundary ghost transient (~7.5 per run) which the T2 form
  grants freely.
- Sector combination: the x-spanning sector is L1's mirror branch
  (⌈24/4⌉ = 6, untouched today) ⟹ **d(24,18) ≥ 6 overall** —
  +1 over the banked ⌈18/4⌉ = 5, with the y-sector at +2 (T1).
  The flagship ≥ 20 was NOT reached; the reason is structural and
  now precisely named (§11.8 residue 1).

**Controls (all PASS, `s6_assembly.json`):**
- (a) the S5 stacks: the species D-rate is present in the tables
  through the caps and the mixed IP never invokes closure — the
  (18,63)/(18,36)/(24,48) sub-2m stacks are PERMITTED by the
  assembly (at (18,63) the same arithmetic gives floor ≤ 126 − 54
  = 72 = the stack weight: the inequality degrades exactly where
  it must).
- (b) (18,12) instantiation: mixed floor 7 ≤ 24 = d((18,12)) ✓.
- (c) the b = 0 (12,12) witness (w18 = 2m − 6): mechanically
  confirmed to live in the wrapped/winding corner (no compact
  cover lift), i.e. the −6 is the admitted wrapped-corner term the
  cover-scoped assembly must and does exclude from its own scope.

### §11.7 Falsified claims (session 6)

- **Charter hypothesis (i)** ("D > 0 only on one drift sign") —
  REFUTED: ghost deficits at both signs, near-symmetric
  (max-D at δ > 0 vs δ < 0 differs by ≤ 0.75 through h = 17).
- **Charter hypothesis (iii)** ("the frontier is achieved by
  species runs; aperiodic segments never beat them") — REFUTED in
  the transient regime (ghosts beat species ~2× below h ≈ 12);
  the sustained regime is species-consistent (the saturation).
- **"First-row boundary-debt charging trims the ghosts"** (this
  session's own working idea) — FALSE: ghosts coast to low-debt
  endpoints (penalized tables ≈ pure tables at u = 1); the debt
  lives deeper than one forced row.
- **"u ≤ 4 exact strata are uniformly feasible"** — u ≥ 3 exact
  enumeration is RED (the min-slab floor forces multi-point rows,
  branching ~120/node, 24M pushes at g ≤ 8); their analytic
  bounds are dominated by u ≤ 2 tables, so nothing is lost.
- **The loose min-window join as a near-exact proxy** — it
  over-grants materially (383 buckets at u = 1, g ≤ 30, up to
  ~5 deficit); matched joins are required wherever affordable.
- Session-internal: the first drift-script lift gate (per-period
  extent ≤ ℓ − 5) would have falsely rejected compact lifts —
  removed after noting wrap-free admissibility already excludes
  false accepts; the first assembly run treated absent u = 3, 4
  strata as absent-not-analytic (unsound direction) — fixed the
  same hour.
- (Respected: witness weights as upper bounds; RED/AMBER/GREEN
  verdicts before every big run; no SAT anywhere; every consumed
  vector re-verified end-to-end; the S4 x_winds classifier not
  used.)

### §11.8 Residue / S7

1. **The boundary-coupling lemma — the sharpest next question.**
   The ghosts are real fragments of the relaxed run semantics, and
   the relaxation is EXACT per run — but two runs sharing one
   heavy slab cannot both coast freely: the ghost's backward
   E-extension forces specific content on the boundary rows, which
   are the neighbouring run's rows.  The loose IP ignores this
   cross-boundary E-consistency entirely; charging it needs the
   interface-matched composition (run tables keyed by 3-row
   boundary states, heavy transitions as weight-≥8 interface
   moves) — the finite-interface form of §10.1's per-cluster
   decomposition.  This is where the T2 gap (~7.5 per boundary)
   lives, and with it the road from floor ~8 to 36 − C.
2. **Deeper caps + dil-6 production**: every +4 in g-cap = +1 on
   the all-light floor mechanically; dil-6 production (~7×) and
   chained/connect-later growth close the enumeration scope.
3. **The wrapped corner at ℓ = 24** (extent ≥ 21 non-winding
   wrap-interacting fragments): torus-mode enumeration (mod-24
   columns, lifted-slip drift) — the ℓ = 12 witness demonstrates
   the corner is real.
4. **The x-sector mirror pass** (pair (B, Ā), 24 rows of width
   18): same engine, mirrored; would lift the overall floor to
   min(d_Y, mirror-floor) beyond 6.
5. **The closed-branch per-seed search** at (24,18) (no relaxed
   E at all — ghosts die): even moderate caps could push the
   all-light branch well past the open-bucket bound; priced RED
   at deep caps today, unpriced at shallow.
6. **Stage 4 (W7 ω-mechanism) untouched** — §10.9 item 3 stands
   as stated.
7. **Lean**: the drift additivity (telescope), the footprint/
   splitting lemma, and the (2 − u/4)h analytic strata are
   two-line-to-decide-shaped; the frontier tables are finite
   certificates once the engine's soundness contract (complete-
   below-cap) is stated.

## §12 SESSION 7 (2026-08-28) — the boundary coupling charged:
## links, the pinch, closed marches, and d_Y(24,18) >= 11

Directive: prove the boundary-coupling lemma — adjacent light runs
cannot both extract full transient deficit through a heavy
interface — and reassemble.  What happened: the lemma got THREE
independent teeth this session, each mechanically verified: (i) the
LINK DECOMPOSITION, which makes the interface an enumerable object
(forward-only marches through the heavy block, E-consistency
enforced end to end); (ii) the PINCH LEMMA, a three-line
combinatorial bound killing double-vacuum interfaces at short
blocks; (iii) the CLOSED MARCH, which replaces the s6 open-bucket
all-light readout with exact cycle closure and EXHAUSTS the u = 1
stratum outright.  The measured interface tax is 2.0-2.75 per
boundary in the certified range (vs the ~7.5 the s6 accounting
granted); the reassembled y-sector floor is d_Y(24,18) >= 11 (T1',
scope-listed), +4 over s6.  En route: an anchor-aliasing gap in the
s6 delta-resolved readouts was found, fixed, and retroactively
discharged; and the session's first parallel launch OOM'd the
machine — the ops rules that now govern every march are recorded in
§12.9.  Scripts `a40_s7_*.py`, data `data/a40/s7_*`;
`validate_banked` green before every stage.

### §12.1 The link decomposition (the bookkeeping that makes the
### interface enumerable)

Split every maximal light run at its first minimum-weight window
(the SEED).  The cyclic walk = seeds + LINKS, where link_i runs
from seed_i (inclusive) upward through run i's forward piece, heavy
block i, run (i+1)'s below-seed piece, ENDING AT seed_{i+1}
(inclusive).  Mechanical consequences (all verified):
- every slab lies in exactly one link except seed slabs (exactly
  two), so sum g_link = 4|v| + sum u_i and the deficit telescope
  D_total = sum_i [D_link_i - (2 - u_{i+1}/4)] is exact;
- drift telescopes with NO uncounted anchor steps — the s6
  assembly's free heavy slips are now MEASURED inside links, so
  closure (sum delta == 0 mod 24, winding included) applies to
  every all-enumerated configuration;
- every E_j is enforced by exactly one link (a forward march
  reaching final row s enforces E_{s-1}; the next link, seeded on
  rows [s-3, s], enforces from E_s): the cross-boundary
  E-consistency the s6 loose join ignored is now structural;
- the forward recurrence is monic in v2[t+1], so a forward march
  from a full-content seed window enumerates EVERY admissible
  continuation: the long BACKWARD ghost coasts of s6 appear as
  post-phase content driven from the block below.  Forward-only
  marches suffice — and the forward tree is the cheap one.

### §12.2 The pinch lemma (proven; the combinatorial tooth)

**Lemma.**  For a heavy block of L <= 3 slabs between light slabs
b (exit) and b+L+1 (entry): rows(slab b+1) = {b-2..b+1} is
contained in rows(slab b) u rows(slab b+L+1) exactly when L <= 3,
hence W_b + W_{b+L+1} >= W_{b+1} >= 8.  Double-vacuum interfaces
are impossible at short heavy blocks — a ghost cannot fade to
vacuum against its own interface.  Scope: the covering fails at
L >= 4 (an interior fat row can feed every heavy slab), though
lone rows are E-dead (E_{t-1} = v2[t] forces them to 0), so the
L >= 4 escape is not free either — unquantified, listed residue.
Verified mechanically (row-covering identity for L = 1..6) and on
the only sub-rate-2 species that crosses the heavy line: TC63's
real blocks have exit + entry = 5 + 6 = 11 >= 8, and all 60 of its
crossing segments are dominated by the analytic link grant
(`s7_validate.json`).

### §12.3 The anchor-aliasing finding (s6 ledger item, discharged)

The s6 March's dominance key (7-row state, h, delta) OMITS the
current slab's anchor: v2[t-3] was dropped from the state, and the
anchor of rows [t-3, t] is not a function of the remaining rows.
Two merged paths can differ in anchor, so their futures differ by a
delta shift: **delta-RESOLVED bucket completeness fails** (at u=1
g<=24, de-aliasing adds 2 reachable buckets), while delta-BLIND
per-h minima are EXACT (the merged representative reaches every
(h', g') the discarded path could, at <= its cost).  Consequences:
s6's mixed floor (delta-blind) stands as banked; s6's all-light
(19, delta == 0) emptiness readouts carried a latent caveat — now
RETROACTIVELY DISCHARGED, because the s7 closed march (alias-free:
its state is the full 4-row window = the slab) proves strictly
stronger emptiness (u=2 to g<=55 vs s6's g<=28).  All s7 engines
carry the anchor in the dominance key; the s7 kmax=0 regression
DOMINATES the s6 fwd table with per-h minima equal.

### §12.4 The closed march and the r in {0, 1} branches

Per-seed forward march whose state is the full last-4-rows window
(8 masks — anchors alias-free by construction); readout after
exactly m steps demands normalized-state equality with the seed and
total drift == 0 (mod ell) (winding +-24 covered).  Emptiness below
a g-cap gives w >= ceil((gcap + 1 - u)/4) for the stratum; the seed
slab's double count is subtracted exactly.

**Controls.**  (i) Positive: seeded at the L12 rate-2 species' own
minimum window at (18, 6) (profile [11,10,7,7,6,7], one L=2 block,
smax/dil relaxed to 8), the march detects EXACTLY ONE closure, at
g = sum slabs + u = 54, delta = 0, phase POST — the certified
d((18,6)) = 12 object, found through its own heavy block at exact
cost.  (ii) The heap and layered engines produce identical tables
and closed sets at equal caps.  (iii) (18,12): the certified
minimum (L12 x2, every slab W = 8) sits in the ALL-HEAVY branch at
exactly 2m = 24 — the s7 assembly admits it by construction — and
the closed marches instantiated at (m, ell) = (12, 18) find
nothing cheaper in the light branches (all-light u=1 empty to
g <= 40 in 0.3 s, one-block u=1 empty to g <= 30; no aborts):
the machinery is coherent at the certified frame.

**Results at (24, 18), y-sector, scope = radius-4/smax-3 growth,
extent <= 34, |delta| <= 30 (every truncation counter ZERO in every
production run):**

| branch | stratum | result | floor |
|---|---|---|---|
| all-light (r=0) | u=1 | **EXHAUSTED at g<=76, trunc_g False** — no walk at ANY weight | infinity |
| | u=2 | empty g <= 55 | >= 14 |
| | u=3 | empty g <= 63 | >= 16 |
| | u=4 / u>=5 | analytic | 18 / 23 |
| r=1, block L<=2, W<=14 | u=1 | empty g <= 42 | >= 11 |
| | u=2 | empty g <= 47 | >= 12 |
| | u>=3 | analytic 17u+8 | >= 15 |
| r=1, L>=3 | any | Dbest(18-L) delta-blind | >= 14 |
| r=1, fat (W>=15) | any | Dbest credit -1.75 | >= 13 |

The u=1 all-light EXHAUSTION is the session's sharpest single fact:
the u=1 stratum's forward tree from every seed dies of
E-consistency + the window rule before 19 slabs at EVERY budget
(53,372 nodes total, 0.5 s) — the s6 "+4 g-cap = +1 floor"
staircase is simply over on that stratum.  STABILITY-CHECKED on
the growth scope: at smax 4, at dil 6, and at both jointly
(g <= 56-60), the tree still exhausts with zero closures
(`s7_closed_k0_stability.json`) — the exhaustion is not an
artifact of the production caps.  (The g46 refinement of
r=1 u=1 was half-run when the OOM incident killed its two heaviest
seeds; 6/8 seeds are banked empty at g<=46 and the two others stand
at g<=42 — the honest claim stays 11, and g>=45 would give 12.)

### §12.5 The link tables and the measured interface tax

Engine: the LinkMarch — phase-scheduled forward march (light
[u,7] -> heavy [8,whcap] x L<=2 -> light [1,7]) with the anchor in
the key, run as a layered BFS (atomic h-layers, two layers in
memory, tables complete through the last finished layer) under the
§12.9 ops guards.  Stratum u=1 tables are COMPLETE to g <= 26 at
heavy class W <= 14 and to g <= 24 at W <= 16 (byseed, no aborts,
all counters zero).  Stratum u=2 is RED at today's budget (the
heavy-entry fan: 33M states at layer 3 whole-tree, ~2 h byseed) —
u=2 links enter the assembly through loose grants only.

**The certified tax table (u=1, L=1).**  J = max certified link
deficit; s6grant = what the s6 accounting granted the same slabs
(fwd table + free heavy + bwd table, delta-blind):

| h | J(h) | s6grant | tax |
|---|---|---|---|
| 4 | 3.25 (g=19, delta=-6) | 5.25 | **2.00** |
| 6 | 6.00 (g=24, delta=-2) | 8.75 | **2.75** |
| 7 | 7.75 (g=25, delta=-2) | 10.50 | **2.75** |

- **There is NO h=5 crossing at g <= 26 at all** — a parity-like
  rigidity of minimal interfaces (cheapest h=4 and h=6 crossings
  exist, h=5 does not); and every certified crossing is L=1 (a
  two-slab block costs >= 27 with its tails).
- The cheapest crossing weighs g = 19 against the naive
  seed+block+post floor of 11: erecting a heavy slab out of a
  stratum-1 seed costs ~8 extra weight in E-consistency alone.
- T_link := max certified [D - (6/7)(h - L)] = **2.61** — the
  per-boundary transient the s6 T2 read as ~7.5 (T0 per run end,
  two ends per interface) collapses by ~2/3 when the interface E's
  are enforced.
- **The two-ghost verdict (charter Stage 1): two ghosts cannot
  share an interface at full deficit.**  Standalone, the s6 matched
  frontier pays ~1.5/slab through the transient; through one
  enforced heavy slab the joint object pays J(h)/h ~ 1.1/slab with
  a per-interface tax of 2-2.75, and the deficit-maximal specimens
  are NOT vacuum-vacuum compositions (§12.6): the vacuum side must
  re-inflate to make the block's weight, exactly the pinch
  mechanism, and the E-forced debt row does the rest.

### §12.6 Specimens (`a40_s7_tax.py specimens`,
### `s7_specimens_*.json`) — the lemma in one object

Deficit-maximal links replayed from march parents and verified
INDEPENDENTLY through CoverFragment (every E including the heavy
block, slab classes, window rule, weight, drift — exact matches
asserted), search-narrowed to |delta| <= 10 (specimen search only;
any found object is fully verified, so narrowing is sound).  The
certified h=6 optimum (D = 6.0, g = 24, delta = -2):

    slabs [1, 2, 2, 4, 8, 7], slip across the block = -1,
    ghost stretches (pre, post) = (3, 0).

Read it: a genuine 3-slab ghost approach ([1,2,2]) must RE-INFLATE
(slab 4) before the heavy slab is even erectable — the pinch
mechanism in vivo (exit 4 + entry 7 = 11 >= 8) — and the post side
is NOT a ghost at all: the E-forced debt row through the block
makes the second coast impossible.  The two-ghost composition of
the charter is refuted in the strongest concrete sense: the
deficit-maximal crossings have ONE ghost side; vacuum-vacuum
interfaces never appear.  The h=7 optimum (D = 7.75, g = 25,
delta = -2, from a different seed) repeats the pattern one slab
longer: slabs [1, 1, 2, 2, 4, 8, 7], slip -1, ghost (4, 0) — a
four-slab coast, the same forced re-inflation, the same non-ghost
entry side.  Both deficit-maximal certified crossings share the
shape; it is the boundary-coupling mechanism, photographed.

### §12.7 Stage-2: the heavy-slip status (measured, not proven)

Within every enumerated link the anchor slip across the block is
MEASURED by the march's drift telescope — the closure constraint on
the r>=2 all-enumerated branch and on the r in {0,1} closed
branches therefore needs no separate slip lemma.  Outside the
enumerated scope the general weight-bound lemma remains OPEN: a
scope-free "slip <= f(block weight)" is not provable by the pinch
alone (left spurs built over many rows can absorb large jumps at
light-granted cost), so long blocks (L >= 3) and fat blocks stay
closure-free in the assembly — the honest boundary of today's
closure lever.  Measured slip on the verified deficit-maximal
specimen: -1 column across its L=1 block — minimal interfaces
barely slip; adversarial slips, if any, live in the fat/long
blocks outside today's enumeration.

### §12.8 The reassembly: d_Y(24, 18) >= 11 (T1')

Branch floors (min = the floor):

- all-heavy: >= 36.  all-light: >= 14 (u=2 cap).  r=1: >= 11
  (u=1 short cap; fat 13, long 14).
- r>=2 (the link DP, stratum-coupled pieces (u -> u_next), seed
  double-count subtracted, per-piece grants = certified J /
  capbound-capped-by-loose / loose with fat credit): DP-closure
  (all links enumerated, sum delta == 0 mod 24) grants 23.5;
  DP-free (>= 1 non-enumerated link, no closure) grants 25.25 at
  the binding composition [(u=2 piece, loose-fat) + (u=1 piece at
  the loose cap)] => **floor 11**.  Evaluated at BOTH heavy-class
  semantics (W<=14 full-depth, W<=16 with fat credit 2.25): both
  give 11 — the fat-heavy escape through the LOOSE side of the
  free branch is the binding leak, 0.25 deficit short of 12.
- **T1' (certificate-shaped, scope-listed): d_Y(24,18) >= 11** —
  +4 over s6's 7; the overall floor stays d(24,18) >= 6 through
  the untouched x-sector mirror (unchanged residue).
- T2' (conjectural tier): with the MEASURED T_link = 2.61 replacing
  the s6 per-run transient 7.5, the extrapolated mixed floor is
  >= 11 as well — the certificate and conjectural tiers now agree,
  because the boundary transient is no longer the dominant unknown.

**Controls (all PASS, `s7_assembly.json`):** (a) the S5 stacks:
W7 x9 at (18,63) is all-light and every closed-march claim is
(24,18)-specific g-cap emptiness far below the stack's g = 288 —
admitted; TC63's stacks cross heavy blocks whose links are
dominated by the analytic grant (60/60) and whose closure sums are
0 mod ell at the stack frames ((24,48): 8 x (+3) = +24 == 0);
(b) (18,12): floor <= 24 with the certified minimum in the
all-heavy branch at exactly 2m; (c) the b=0 (12,12) witness stays
in the wrapped/winding corner (no compact cover lift,
`s6_drift.json`), the -6 an admitted scope term; at l=24 the
wrapped corner remains a listed condition.

### §12.9 Falsified claims and ops incidents (session 7)

- **The first parallel launch OOM'd the machine** (concurrent
  heap marches at ~10M entries each; the user force-quit them).  Standing rules for
  every s7+ march, enforced IN CODE: sequential runs only; hard
  2 GB RSS budget checked between AND inside layers; frontier
  > 3M states = RED abort; atomic layers so an abort still leaves
  tables complete through `complete_h`; reprice at reduced caps
  before any relaunch.  The heap engine was replaced by a layered
  BFS (identical tables, proven on-shell) for exactly this reason.
- **Whole-tree link marches are RED** at production caps (u=1 g24:
  6.2M states at layer 4; u=2 g26: 33M at layer 3, caught only at
  7 GB by the between-layer check — the mid-layer check exists
  because of it).  byseed (per-seed sequential, tables merged) is
  the production mode; dominance across seeds is only an
  optimization, so the merged tables are identical.
- **u=2 link enumeration RED at today's budget** (byseed ~2 h);
  killed twice; u=2 pieces enter via loose grants.  W<=16-class
  u=1 enumeration RED at g26 (4.3M fan on one seed), complete at
  g24.
- **s6 delta-resolved bucket completeness** — the anchor-aliasing
  item of §12.3: the METHOD is retired; no banked s6 number
  changes (mixed floors delta-blind; all-light superseded by
  stronger s7 results).
- **The uncoupled r>=2 free DP** (every post piece granted the
  stratum-1 bwd table) leaked to floor 9 before the (u, u_next)
  stratum coupling was added — never claimed, caught same hour.
- byseed `complete_h` initially conflated "tree died" (complete at
  every h) with "aborted at h" — fixed; a dead frontier has no
  descendants.
- (Respected: witness weights as upper bounds; RED/AMBER/GREEN
  before every big run; no SAT anywhere; every consumed vector
  re-verified end-to-end; nothing re-proposed from §9.8/§10.8/
  §11.7.)

### §12.10 Residue / S8

1. **r=1 u=1 at g >= 45** (2 seeds, ~1-3 h guarded) lifts r=1 to
   12; the r>=2 free branch then binds at 11 — its levers, in
   order: boundary-weight-RESOLVED fwd/bwd tables (charge the
   pinch inside the loose split — the 0.25 gap), kmax=3
   enumeration, u=2 links via a heavy-entry-restricted seed fan.
2. The h=5 interface gap (no crossing at g <= 26): understand the
   parity mechanism — it smells like the two-step debt row forcing
   an even/odd obstruction; a proof would turn the tax table into
   a tax THEOREM for minimal interfaces.
3. Deeper u=1 link caps certify J at h in [8, 13] where today the
   loose cap rules (needs g ~ 30-35: RED whole, AMBER byseed
   overnight).
4. The x-sector mirror pass (lifts the overall floor past 6); the
   wrapped corner at l=24; the W7 omega-mechanism (§10.9 item 3)
   — all unchanged.
5. Lean: the pinch lemma and the row-covering identity are
   decide-shaped; the link/closed tables are finite certificates
   once the layered-march soundness contract (complete-below-cap,
   anchor-keyed dominance) is stated; the L12 closure control is
   the natural first mechanized witness.

## §13 SESSION 8 (2026-08-28) — the x-sector mirror: d_X(24,18) >= 12
## and the first double-digit OVERALL member floor, d(24,18) >= 11

Directive: instantiate the S6/S7 transfer machinery on the
theta'-image pair (B, Abar) and lift the OVERALL (24,18) floor past
min(11, 6); then attack the heavy-slip weight bound.  What happened:
the mirror turned out to be STRUCTURALLY CLEANER than the y-lane —
the sector definition itself removes the wrapped corner, the winding
branch, and the drift-mod-ell subtlety, and the walk's 24-slab
telescope makes every stratum except (u = 1, small n_H) die by
COUNTING at the target — so the whole x-sector floor reduces to
per-seed closed marches that run in minutes-to-hours.  Result:
**d_X(24,18) >= 12 (T1-mirror, scope-listed), hence
d(24,18) >= min(d_Y, d_X) >= 11 — the first double-digit
unconditional floor for a tour-de-gross member at (r,b) = (3,1),
with the y-sector now the binding side**
([[864,12]]; the conjectured value is 36).  Engines
`a40_s8_xlane.py` (mirror system + envelope-keyed closed march +
replay verifier), `a40_s8_slip.py` (Stage 2); data `data/a40/s8_*`;
`validate_banked` green before every stage.

### §13.1 The mirror reduction (exact; the sector is the theorem)

theta': (x, y) -> (y^-1, x) is a ring automorphism of
F2[x^pm, y^pm] with (A o theta', B o theta') = (B, Abar): the AB
code on Z_l x Z_m is isomorphic to the BAbar code on Z_m x Z_l
(support map (a, b) -> (b, -a), blocks preserved, weights, cycles,
stabilizers, and classes preserved — control-verified end to end).
Lemma K transports (regularity of (A, B) maps through the
automorphism), so a nontrivial X-logical of AB-(24,18) with a
cyclic y-gap >= 4 (the x-sector = the complement of the S6/S7
y-sector) has no x-gap >= 4, and its theta' image is a y-walk
object of BAbar-(18,24) with:

- all 24 walk slabs NONEMPTY (the telescope 4|v| = sum of 24 slab
  weights, cyclic);
- content on an x-INTERVAL of extent <= 18 - 4 = 14 — the y-gap
  cuts the content circle, so the interval embedding is the
  canonical cover lift, there is NO wrapped corner, NO winding
  branch, and closure demands total anchor drift == 0 EXACTLY
  (|dlt| <= extent - 1 <= 13 < 18 makes mod-18 closure literally
  equal to exact closure).  The three scope conditions the y-lane
  carries for these phenomena are THEOREMS at the mirror.

The mirror recurrence (derived generically from the supports and
machine-checked against the bit kernel): X-cycles of (B, Abar)
satisfy A u1 + Bbar u2 = 0, i.e.

  E'_t: u1[t] + u1[t-1] + x^3 u1[t+1] + (1+x^-1) u2[t]
        + x u2[t-3] = 0

— monic in u1[t+1] up to the unit x^3 (block 1 forced forward at
offset -3, block-2 rows free: §9.1's mirrored convolutional
structure).  Footprint spans (4,4) (connectivity lemma verbatim);
the mirror tooth is (Bbar, A): blk1 (0,0),(-1,0),(1,3); blk2
(0,0),(0,1),(3,-1); the H = 5 window rule is the same local
reduction.  Forced-row reach: supp(u1[t+1]) sits within [-4, -2]
of the window's columns (y-lane: [-3, +1]) — machine-derived.

### §13.2 The counting collapse (why the mirror is cheap)

For the (24,18) member the walk is 24 slabs of content width <= 14
(vs the y-lane's 18 slabs of width <= 34 + wrapped corner).  With
every slab >= 1 and heavy slabs >= 8 (n_H = # heavy slabs):
4|v| >= 24 + 7 n_H, and min-slab u >= 2 gives 4|v| >= 48.  At
target floor 11: only u = 1, n_H <= 2, per-heavy W <= 19 survives
(u >= 2 -> 12; n_H >= 3 -> 12; any heavy W >= 20 -> 11); at target
12 only u = 1, n_H <= 3, W <= 23 survives.  The infinite tail of
strata the y-lane had to march is counting-dead at the mirror.

### §13.3 The engine and its verification chain

`MClosedMarch`: per-seed layered closed march (S7 semantics: full
last-4-rows window state — alias-free anchors; completion prune;
RSS <= 2 GB and frontier <= 3M ops guards in code; per-slab class
light [u,7] / heavy [8, whcap] with an n_H counter instead of a
block phase machine, so adjacent AND separated blocks are one
lane).  Two mirror-specific soundness points, both new:

- **the envelope (xlo, xhi) is part of the dominance key**: the
  extent cap is the SECTOR DEFINITION here, so merging paths with
  different envelopes under a binding extent test would
  under-enumerate (kept-wide representative dies at the cap where
  a discarded-narrow path would have closed).  The y-lane could
  afford envelope-blind dominance because extent-34 was a listed
  scope, not a claim boundary.
- **trunc_dcap is benign by a width argument**: dlt is the anchor
  offset from the seed in a common frame, so |dlt| <= envelope
  width - 1; any candidate with |dlt| > 20 has width >= 22 > 14
  and is extent-doomed regardless of check order (R0 counted 329
  such candidates; extent prunes are the sector theorem working,
  7.4M in R0).

Verification chain (all green, `s8_xlane_selftest.json`,
`s8_xlane_control.json`):

- selftest: the recurrence terms, unique pivot, and (4,4) footprint
  re-derived from the pair supports; the forced-row bit kernel ==
  the generic set solve on 2000 random states; **400 random forced
  walks re-verified admissible through the independent
  MirrorFragment (generic E')**; drift additivity telescope; pinch
  row-covering (pair-free, L <= 3).
- **L12' control**: the BAbar p = 6 straight atlas has its own
  66-member weight-12 family (the theta' twin of L12, banked S4
  data).  Lifted: profile [7,8,9,7,8,9], extent 3, drift 0.  The
  closed march seeded at its own minimum window detects EXACTLY
  its closure at exact cost g = 55 = sum(profile) + u, dlt = 0,
  n_H = 4; the parent-tracked replay reconstructs the object, the
  MirrorFragment re-verifies every E', and the torus embedding
  certifies it nontrivial on BAbar-(24,6) AND its theta' image
  nontrivial on AB-(6,24) — the reduction validated end to end.
- **TC63' control (the mirrored species, §10.4's demand)**: the
  twisted-atlas transform for the BAbar <(3,6)> class regenerates
  the w10 species; its mirror cover lift is admissible with
  gauge-free anchor drift **+3 per 6 rows**, slabs [9,8,5,6,6,6],
  deficit 2/period, re-verified nontrivial on the sheared (42,6,3)
  torus — the mirror drift instrument measures the known mirrored
  species correctly.
- **(18,12) coherence via its x-sector**: the mirror march at
  BAbar-(12,18) (content interval <= 8), u = 1, n_H <= 2,
  W <= 19, g <= 48: **zero closures** (4 of 8 seeds exhaust
  outright; 120 s) — no contradiction with the certified
  d((18,12)) = 24, and the marched branches are empty far below
  the certificate's weight.
- **the b = 0 witness**: the a36 (12,12) witness occupies 4/12
  y-rows with max cyclic y-gap 2 <= 3: it is a Y-SECTOR object —
  the mirror scope excludes nothing about it and needs no wrapped
  admitted term (its -6 stays charged to the y-side scope, §11.6
  control (c) unchanged).
- **member protection, mirrored**: TC63' needs 6 | 24 (yes) and
  (24/6)*3 = 12 == 0 mod 18 (no) to wrap the member's mirror walk
  — and in-sector winding is impossible, so the only sub-rate-2
  species of the mirror alphabet at p <= 8 (the twisted atlas is
  a both-lanes census) cannot appear.  No species input is
  CONSUMED by the floors (they are enumerative + counting), so
  the open W7-analog question for the BAbar lane is not a scope
  condition — listed residue.

### §13.4 The floors: d_X(24,18) >= 12, OVERALL d(24,18) >= 11

Two production passes, both ALL 8 seeds, ZERO closures, zero
aborts (every guard green; all numbers re-derived from the banked
JSONs at close):

- **R0** (`s8_prod2418_u1_nh2_g42_s0_8.json`): u = 1, n_H <= 2,
  W <= 19, g <= 42, extent <= 14, smax 3 / dil 4 — per-seed trees
  EXHAUST by h ~ 7 (largest layer 1.30M), 685 s, peak RSS 786 MB.
- **R1'** (`s8_prod2418v2_u1_nh3_g46_s0_8.json`, the streaming
  engine of §13.6): u = 1, n_H <= 3, W <= 23, g <= 46, same
  scope — 31.2M states, layer peak 180,800 (the flush threshold),
  5551 s, RSS ~0.25 GB throughout.  The run that had FAILED
  operationally before the streaming rework (5/8 seeds aborted;
  §13.6) completes clean.

Branch table for a class-minimal nontrivial x-sector logical
(4|v| = sum of the 24 slab weights; every slab >= 1):

| branch | source | floor |
|---|---|---|
| u >= 2 | counting (2/slab x 24 = 48) | >= 12 |
| u = 1, n_H >= 4 | counting (24 + 28)/4 | >= 13 |
| u = 1, some heavy W >= 24 | counting (>= 47)/4 | >= 12 |
| u = 1, n_H <= 3, W <= 23 | R1' march empty g <= 46 | >= 12 |

**T1-mirror (certificate-shaped, scope-listed):
d_X(24,18) >= 12.**  Scope: smax-3 / dil-4 prefix-connected input
placement (the growth scope, same discipline as S7) — and NOTHING
ELSE: no wrapped corner, no winding branch, no boundary
relaxation, no loose join, no conjectural tier.  Stability: the
smax-4 and dil-6 batteries (§13.7 item 1) ran at the R0
configuration (g <= 42, n_H <= 2), so the >= 11 floor is
stability-backed on both scope axes; the 12th unit rests on R1'
at production scope only (g46-level stability = listed residue).
Combined with S7:

> **d(24,18) >= min(d_Y, d_X) = min(11, 12) >= 11** —
> [[864,12, >= 11]] unconditional (X side; Z side by transpose
> duality as in §3.2), up from 6.  The first double-digit member
> floor at (3,1); the sector split is two-sided instrumented —
> 11 (y, S7 boundary-coupled) vs 12 (x, S8 counting + march) —
> and the Y-SECTOR is now the binding side: the next unit of
> overall floor must come from §12.10's y-levers, not from x.

### §13.5 Stage 2 — the heavy-slip bound (caps proven; census built)

- **Per-step anchor caps (LEMMA, machine-derived constants,
  scope dil-D growth)**: every new cell (forced or input) sits
  within [min(window) - max(D, F_left), max(window) + max(D,
  F_right)] where the forced-row reaches are y-lane [-3, +1] and
  mirror [-4, -2] (derived from the supports; verified on 6000
  random states).  With D = 4: **A_{j+1} >= A_j - 4 (min gauge)
  and Amax_{j+1} <= Amax_j + 4 (max gauge), per step, both
  lanes.**  The caps are ONE-SIDED PER GAUGE: a min-anchor RIGHT
  jump is a strand death, bounded by the slab span, and no local
  footprint constant bounds it by local weight — the two-gauge
  closure bookkeeping (both telescopes vanish on closed walks) is
  the assembly route that survives this asymmetry.
- **The exhaustive slip census** (`a40_s8_slip.py census`,
  `s8_slip_u1_g24/g26.json`): the S7 link march with the
  heavy-entry drift IN the dominance key from HEAVY onward, so
  tab_slip[(L, gH, slip)] is complete below cap (the S7
  instrumentation was representative-only and its banked tables
  carried empty slip dicts).  Regression: the banked byseed
  pre/link tables reproduced EXACTLY at g <= 26, and the node
  count matches the banked run to within the seed count
  (17,005,749 vs 17,005,741) — the slip keys add nothing at these
  caps.  **The complete census: every certified crossing is
  L = 1, gH = 8, and the realized slips are {-3, -2, -1, 0}
  (min_g 24, 26, 19, 24) — |slip| <= 3, and NO POSITIVE
  (rightward) slip exists below g <= 26 at all.**  §12.7's
  "measured slip -1 on the optima" was representative-thin: the
  anchor can move 3 columns left across a minimal interface at
  near-optimal cost.
- **Species verification** (`species`, `s8_slip_species.json`):
  per-step anchor increments of the banked lifts — W7 at 18/24
  [-2, +2], TC63 at 18/24 [0, +2], TC63' mirror [0, +2] — all
  within the left cap.  Zero violations anywhere in the enumerated
  data: the charter's verification demand is met for the LEFT
  direction.
- **THE VERDICT (the honest split).**  (i) LEFTWARD slip is a
  LEMMA: telescoping the per-step cap across an L-slab block
  (L + 1 steps) gives leftward slip >= -4(L + 1), i.e.
  |leftward slip| <= 4(L+1) <= W_block/2 + 4 using W >= 8L — the
  charter's affine shape with machine constants (c, c0) = (1/2, 4),
  scope dil-4.  (ii) RIGHTWARD slip has NO footprint cap (a
  rightward jump is a strand DEATH onto a pre-existing right
  neighbour) and the two-sided charter lemma is OPEN with two
  named candidate resolutions pulling opposite ways:
  - REFUTATION mechanism: the right neighbour is built earlier by
    a dil-4 staircase (>= 1 cell per 4 columns of reach) whose
    cost lands in LIGHT slabs far from the block — slip unbounded
    in the block's own weight;
  - PROOF mechanism (debris accounting): (1+x^-1) has no finite
    kernel, so every staircase row leaves >= 2 forced v2-debris
    cells in the adjacent slabs, and the debris wake costs
    ~1-2/row until killed — hand simulation prices a rightward
    relocation at ~4 weight per 4 columns CONCENTRATED at the
    transition rows, i.e. slip <= ~1.0 x W(transition) + c0: the
    relocation rows themselves go heavy and the jump is charged
    to a BLOCK after all.
  The census's empty positive side (g <= 26) is evidence for the
  debris mechanism's strength at low cost; the discriminating
  experiment is the census at g ~ 28-32 (queued): positive slips
  appearing at light-adjacent cost refute, their continued absence
  at growing caps supports the debris-accounting lemma.
  **EXECUTED at g <= 30** (h <= 12, |dlt| <= 16;
  `s8_slip_u1_g30.json`, streaming-chunked byseed, no aborts,
  regression PASS, 34M nodes / 203 s): 14 classes, all L = 1,
  gH in {8, 9}; realized slips fill the ENTIRE left range
  {-8..0} at gH = 8 and STILL contain NO positive slip.  Two
  sharp consequences: (a) **the left cap is TIGHT** — slip -8 =
  -4(L+1) = -(W/2 + 4) is realized at g = 29 (table-certified;
  specimen replay is residue): the lemma saturates, and the slip
  is controlled by the STEP COUNT, not the block weight; (b) any
  |slip| <= c (gH - 8) + c0 form needs c0 >= 8 — at the MINIMUM
  block weight the anchor already jumps the full step allowance,
  so the charter's small-constant weight form is refuted while
  the L-form (equivalently W/2 + 4 via W >= 8L) stands, tight.
  Rightward emptiness now extends through g <= 30.
  (iii) Assembly consequence: the S7 free branch's binding
  compositions bank NEGATIVE (leftward) drift on enumerated
  pieces and need the free link to RETURN rightward — exactly the
  unresolved direction.  If the debris route proves out, the
  closure lever extends to ALL branches (the charter's Stage-2
  goal) with constants ~(1, c0); if the staircase route wins, the
  remaining gap is the L-W return-cost wall of §10.7 in
  per-interface form.  Either way the missing piece is now ONE
  named inequality about rightward relocation cost.

### §13.6 Falsified claims and incidents (session 8)

- **The floor-12 mirror run as FIRST parametrized was RED** (seeds
  0-1 > 3M frontier at h = 4; seed 5 RSS; seeds 6-7 spuriously
  killed at h = 1 by the ru_maxrss flaw below) — banked as partial
  (seeds 2-4 complete-and-empty at g <= 46), no claim made from
  it; SUPERSEDED same session by the streaming R1' (all 8 seeds
  clean, §13.4).
- **ru_maxrss is a LIFETIME PEAK, not current RSS**: the S7-style
  guard `getrusage(...).ru_maxrss` trips every subsequent seed
  once one seed spikes — R1's seeds 6-7 were killed at h = 1
  without running.  All S8 engines now guard on ps-based CURRENT
  RSS (the S7 engines retain the flaw; their banked runs were
  whole-run aborts where it made no difference — flagged for any
  future byseed reuse).
- **§12.7's "measured slip -1" was representative-thin** — the
  exhaustive census realizes slip -3 at the same block class
  (L = 1, gH = 8) barely above optimal cost.  No banked number
  changes (the S7 assembly consumed dlt-resolved tables, not the
  slip observations).
- **Layer chunking took three tries to be sound-AND-in-budget**
  (all failures caught by the in-code guards, no OOM): (1) a
  split check at the top of the layer loop splits too LATE — the
  fat child layer is already materialized (RSS abort); (2)
  disk-spilling the chunk slices still holds every ancestor
  layer across the recursion (RSS abort); (3) the fix that works:
  STREAMING — flush the child dict to disk whenever it crosses
  the threshold DURING expansion, then recurse per spilled part
  (cross-part dominance lost = over-exploration only, sound;
  g26 regression byte-identical, node count equal).  This is the
  engine change that completed the g30 census and unlocks the
  floor-12 re-run.
- Session-internal: the first Control-B seed-window indexing used
  the subfragment slab offset off by 3 (caught by the exact-cost
  closure assertion before any bank); the module initially
  computed the tooth's blk2 support twice through a redundant
  double-bar expression (cleaned, no numeric effect); the Stage-2
  species/`main` dispatch was appended after the `__main__` guard
  (NameError on first invocation, moved, no data effect).
- (Respected: witness weights as upper bounds; RED/AMBER/GREEN
  cost verdicts before every big run — the coherence and R0 runs
  were priced by the (18,12) pilot; sequential runs only, 2 GB
  RSS + 3M frontier guards IN CODE, no /tmp, no SAT anywhere;
  every consumed vector re-verified end to end; nothing
  re-proposed from §9.8/§10.8/§11.7/§12.9.)

### §13.7 Residue / S9

1. **Mirror stability battery COMPLETE**: smax-4 at g <= 42 —
   zero closures, all 8 seeds, NO aborts (75 min,
   `s8_stabA_u1_nh2_g42_s0_8.json`); dil-6 at g <= 42 — zero
   closures, all 8 seeds, NO aborts (32 min via the streaming
   engine, `s8_stabB_u1_nh2_g42_s0_8.json`).  NEITHER growth
   scope cap is binding for the floor-11 claim — the same
   stability shape as §12.4's y-lane exhaustion.
2. **The floor-12 x-push: EXECUTED** — R1' clean (§13.4),
   d_X >= 12 banked.  The next x-rungs, priced: (a) g46-level
   stability re-runs (smax 4 / dil 6 at nH <= 3 — the streaming
   engine makes them memory-flat; hours-scale); (b) x >= 13 needs
   the u = 2 stratum marched (284 seeds, every slab >= 2, budget
   ~50: tiny trees) plus nH <= 4 at g <= 50 — but the y-sector
   binds the overall floor at 11, so x-rungs are generalization
   evidence, not floor progress.
2b. **The slip -8 specimen replay**: reconstruct the saturating
   (L=1, gH=8, slip=-8, g=29) link through CoverFragment
   (per-seed parent-logged rerun at gcap 29) — turns the
   table-certified saturation into a verified witness.
3. **The all-light exhaustion probe** (n_H = 0, gcap effectively
   unbounded): does the mirror all-light stratum exhaust CAP-FREE
   like the y-lane's (§12.4)?  Cheap, bankable as a mirrored
   structural fact.
4. **The rightward-slip witness** (settles Stage 2's open half):
   construct the two-strand fragment — a staircase-built right
   strand + the anchor strand dying into the block — explicitly at
   g ~ 30+ (census extension or hand construction through
   CoverFragment); a verified witness refutes the two-sided local
   lemma OUTRIGHT, and its measured cost-per-column IS the
   empirical return-cost constant the L-W wall needs.
5. **The W7-analog question for the BAbar lane** (does a winding
   constant-weight species exist in the mirror alphabet beyond
   the twist-compact census?) — irrelevant to the banked floors,
   relevant to the mirror's momentum-budget narrative.
6. **Stage 3 (the h=5 / L=2 interface rigidity, §12.10 item 2)
   NOT REACHED** this session — no proof work done, the
   parity-like puzzle stands as S7 left it.  One new data point
   from the census bears on the L-half: at g <= 30 (h <= 12)
   every certified crossing class is STILL L = 1 (`slip_rows` of
   `s8_slip_u1_g30.json`: no L = 2 entry) — the two-slab-block
   exclusion extends 4 cost units beyond where S7 measured it.
   (The h = 5 half is not re-derivable from the banked census
   fields — the run does not store the h-resolved link table —
   so its extension is unverified, listed as-is.)
   **Stage 4 ((18,12) coherence + ell = 30) PARTIALLY reached**:
   the (18,12) x-sector coherence ran as a CONTROL (§13.3, zero
   closures, 4/8 seeds exhaust); the ell = 30 column's sector
   floors were not attempted.
7. Lean: the theta' reduction is a monomial-substitution
   isomorphism (formalizes as a support bijection); the counting
   collapse is arithmetic on the slab telescope; the mirror
   closed tables are finite certificates under the same
   layered-march contract as §12.10 item 5.

## §14 SESSION 9 (2026-09-01) — the rightward-relocation inequality
## DECIDED (census to g <= 40 + the mu-echo lemma), and cashed

Directive: decide S8's named open inequality — does a dil-4
staircase's forced E-debris charge ~1 weight per column of
rightward anchor transport, or does a two-strand ghost-carrier
evade it above g ~ 30 — then cash the two-sided slip control into
the y-sector assembly.  What happened: the census route ran CLEAN
to g <= 40 at full production caps (rightward slip: still EMPTY —
now 21 cost units above the cheapest crossing), and the debris
hand-argument came back STRONGER than commissioned: its zero-input
core is a PROVEN three-line lemma (R0, the mu-echo) with a 3.1M-
pair exhaustive confirmation, so the two-strand carrier is refuted
outright in the zero-input sector and every rightward unit is
input-charged.  The same deep run doubles as the u = 1 link-table
deepening the S7 assembly's free branch needed: r >= 2 lifts 11 ->
12 (the probe-predicted 96-quarter bc-wall), r = 1 becomes the
unique binding branch, and the streaming ClosedMarch port (S8's
rework applied to the closed lane, 3 port controls) reruns the two
seeds the S7 frontier guard killed.  Engines `a40_s9_slipext.py`
(census/specimen driver over the UNCHANGED S8 SlipLinkMarch),
`a40_s9_debris.py` (V0-V5 + R0), `a40_s9_closed.py` (streaming
ClosedMarch), `a40_s9_probe.py` (DP sensitivity, no claim tier),
`a40_s9_assemble.py` (glob-shim assembly, output s9_assembly.json;
the banked s7_assembly.json untouched); data `data/a40/s9_*`;
`validate_banked` green before every stage.

### §14.1 The census verdict (the measured side, now 4 rungs deep)

The S8 slip census extended with the S8 engine UNCHANGED (my
driver adds only output naming, the g30 slip-table regression, and
a strictly additive tab_link_L[(h, L, dlt)] recording — dominance
key untouched; the s7-g26 pre/link/link_L regressions and the g30
slip regression pass on every run):

| gcap | hcap | dcap | nodes | wall | classes | positive slips |
|---|---|---|---|---|---|---|
| 30 (S8) | 12 | 16 (1272 trunc) | 34M | 203 s | 14 | none |
| 32 | 13 | 20 (0 trunc) | 42M | 240 s | 16 | none |
| 34 | 14 | 22 (0 trunc) | 50M | 381 s | 23 | none |
| **40** | **19** | **30 (0 trunc)** | **79M** | **561 s** | **42** | **NONE** |

- **Rightward block slip does not exist at g <= 40** (u = 1 links,
  L <= 2, W <= 14, smax 3 / dil 4, extent <= 34; |dlt| <= 30 with
  ZERO truncations, so the g30 run's dcap-16 scope caveat is
  DISCHARGED through g <= 40).  The cheapest crossing costs 19;
  the census now spans 21 units above it with the positive side
  empty.
- **The pre-phase wall is frozen at +3 = the borrow radius,
  absolutely**: at g <= 40, h <= 19 the pre table's rightward rows
  are exactly dlt +1 at h 2..7 (min-g 2, 3, 4, 6, 20, 28), +2 at
  h 3..5 (5, 6, 9), +3 at h 3..4 (6, 11) — NOTHING at h >= 8 or
  dlt >= +4 at ANY cost <= 40.  Holding a +1 relocation for >= 5
  rows exceeds the whole 40 budget: the holding rows are forced
  heavy-class (~8/row measured) — "the relocation rows themselves
  go heavy", S8's debris prediction, photographed in a table.
- **Whole-link net drift caps at +2** (h 6 g 34, h 7 g 35): the
  pre transient survives a crossing only in its first two columns.
- The left side keeps saturating: L = 1 slip -8 = -4(L+1) (S8,
  tight); L = 2 now realizes slip -8 (g 35) of its -12 allowance.

### §14.2 LEMMA R0 (the mu-echo) and the mechanized debris battery

`a40_s9_debris.py`, `s9_debris.json` — the checkable steps of the
hand argument, in the pinch style (V1-V4 random-state identities
on the production bit kernel; V5 an exhaustion):

- **V0 (borrow radius, structural)**: reading E at column c, the
  only sources strictly right of c are v1[t](c+1) and
  v2[t-1](c+3); the radius is 3 = the x^-3 term of B (mirror: 4).
- **V1 (kill set)**: every v2 vacation event is killed by an odd
  subset of {v1[t](c), v1[t](c+1), v1[t-3](c-1), v2[t-1](c+3)}
  (151,546 events).  **V2 (leftward spawn)**: a frontier v1 cell
  forces v2[t+1](c-1) = 1 XOR a borrow at v2[t-1](c+2) (47,683
  events) — leftward motion is the free direction.  **V3**: an
  anchor rise is an exiting-row event (10,647 rising steps).
- **V4 (photos)**: a single v2 cell persists verbatim (1 weight/
  row); the one-cell K1 input at c+1 kills it and re-seeds
  EXACTLY at c+1 (retreat <= 1 column/row, >= 1 input/column,
  killers pairwise distinct by column); the K2 borrow at
  (c+3, t-1) kills with no re-seed; a 2-row v2 strand sheds x^-3
  debris LEFTWARD (the B-term powers left motion); and A's x y^3
  term gives every input a +1-column echo four rows up — the
  structural germ of S8's "staircase".
- **V5 (the K2-only channel, EXHAUSTED)**: with v1 == 0 the
  evolution v2[t+1] = v2[t] + x^-3 v2[t-1] is deterministic, so
  the zero-input retreat question is FINITE: over all 3,145,728
  translation-reduced initial pairs (width 11, 12 rows) the max
  net rise of mu_t := min supp(v2[t] u v2[t+1]) is **0**.
- **LEMMA R0 (mu-echo; PROVEN, three lines).**  In zero-input
  dynamics mu never rises: if the min cell c = mu_t lies in
  v2[t+1] it survives into the next pair; otherwise it is v2[t]-
  only and E at (t+2, c-3) reads v2[t+2](c-3) = v2[t+1](c-3) +
  v2[t](c) = 0 + 1 = 1 UNCONDITIONALLY (nothing in v2[t+1] lies
  below mu_t) — mu_{t+1} <= c - 3, a strict >= 3 DROP.  QED.
  (53,484 random case-2 events machine-checked.)  Consequence:
  **every column of rightward anchor transport is input-charged —
  the two-strand ghost-carrier is REFUTED OUTRIGHT in the
  zero-input sector**, and an input-financed evader must pay the
  K1 diagonal (>= 1 input per column, column-distinct) while its
  borrowed cells recurse onto their own columns' deaths.

**The inequality's S9 status (the honest grade)**: leftward slip
>= -4(L+1) is S8's lemma (tight); rightward slip <= 0 holds
CENSUS-EXACT through g <= 40 in-scope; R0 makes the zero-input
half a THEOREM; the per-input constant of the financed half is
measured (+1 costs 3-8 g held, wall at +4 absolute <= 40) and
hand-derived (K1 re-seed diagonal), not yet a theorem — that
residual constant is the only daylight left in the charter's
two-sided local lemma, and no assembly claim below consumes it.

### §14.3 The cash-out machinery (probe -> deepening -> streaming
### seeds)

- **The DP sensitivity probe** (`a40_s9_probe.py`, offline on
  banked tables, NO claim tier) located the cheapest sufficient
  deepening before any big run: S7's binding DP-free composition
  [(u=2 h=3 fat-loose, 6 q) + (u=1 h=17 at the loose cap, 95 q)]
  = 101 quarters is DOUBLY fictional — the fat piece violates the
  proven pinch (exit u=2 seed + entry weight-1 slab = 3 < 8), and
  the h=17 loose grant sits ~20 q above the real link-cost trend —
  and u=1 cap-deepening ALONE moves it: at g38 free drops to 97 q
  (floor 12); at g40 the DP rebalances onto a STABLE two-piece
  bc-wall [(u1->u1, h=10, loose, 48 q) x 2] = 96 q that no deeper
  cap, u=2 short-h table, or modeled cert touches.  Floor 13 would
  need attacking the bc loose grants themselves (pinch-charged
  splits / delta-resolved closure) — not needed, r=1 binds first.
- **The deepening = the census run**: one engine pass (g <= 40,
  hcap 19, dcap 30) serves both Stage 1 (the census of §14.1) and
  Stage 2 (`s9_link_u1k2g40.json`, an s7_link-schema table:
  complete_h 19, ZERO dcap/extent truncs, certs 93 buckets).  The
  real certs stayed under the caps exactly as the probe's J-trend
  model predicted; the class-16 lane also gains the g40 certs but
  its own completeness stays at g24 (floor 11 there; the
  max-over-class-semantics rule takes 14's 12 — each class split
  is individually complete, so the max is sound, as in S7).
- **The streaming ClosedMarch port** (`a40_s9_closed.py`): the S8
  spill-and-recurse rework applied to the closed lane, horizon
  readout per spilled part (sound: streaming loses only cross-part
  dominance = over-exploration; zero closures on an over-explored
  tree is still zero; every part reaching h = m + 1 is read out).
  The S7 module's ru_maxrss lifetime-peak guard (§13.6) is patched
  to ps-based current RSS.  Port controls ALL PASS: C1 exact
  node-count identity vs the banked seed-3 shard (355,673 popped,
  0 closed, threshold off); C2 spill-stability (threshold 20k,
  same tree); C3 the L12 (18,6) positive control detected THROUGH
  forced spills (threshold 64) at exact g = 54, dlt = 0.
- **The two S7 frontier victims fall in minutes**: seed 1 = 12.0M
  popped / 115 s, seed 2 = 18.0M / 134 s, both ZERO closures,
  zero truncs, no aborts, RSS ~0.3 GB throughout
  (`s9_closed_pk2sh_g46_s1_2/s2_3.json`) — the g46 cover
  {0..7} is complete and r=1 short_u1 lifts to
  ceil((46+1-1)/4) = 12.

### §14.4 The floors: d_Y(24,18) >= 12 (T1''), OVERALL
### d(24,18) >= 12

`a40_s9_assemble.py` (the S7 assembly logic verbatim, inputs =
banked s7 files + the S9 deepenings via a glob shim; output
`s9_assembly.json`, the banked `s7_assembly.json` untouched):

| branch | floor | what moved |
|---|---|---|
| all-light (r=0) | 14 | unchanged (u=1 EXHAUSTED, u=2 g55) |
| r=1 | **12** | short_u1 11 -> 12 (the two streamed seeds) |
| r>=2 | **12** | free 101 q -> 96 q (the g40 deepening) |
| all-heavy | 36 | pointwise |

**T1'' (certificate-shaped, scope-listed): d_Y(24,18) >= 12** —
+1 over S7's T1'.  Scope: radius-4 / smax-3 prefix-connected
growth (the 12th unit rests on production scope; the S7 stability
batteries backed the 11-level), extent <= 34, wrapped corner at
l = 24 a listed condition, |delta| <= 30 per fragment with every
truncation counter ZERO in every consumed run, heavy classes
evaluated at both W<=14 / W<=16 semantics.  Controls re-run TODAY,
all PASS: (a) S5 stacks admitted (TC63's 60 segments dominated,
W7 all-light above every cap); (b) (18,12) floor <= 24 with the
certified minimum in the all-heavy branch at exactly 2m; (c) the
b=0 witness stays in the wrapped/winding corner (its -6 an
admitted scope term at l=12); (d) the g40 link table clean; plus
the full S7-link control (ClosedMarch-vs-LinkMarch cross-check +
L12 closure at exact cost) and the S8 xlane selftest + control
battery (L12' theta'-both-tori closure replay, TC63' drift +3/6,
b=0 Y-sector classification, mirrored member protection
4 x 3 = 12 != 0 mod 18) re-verified green in today's environment.

> **d(24,18) >= min(d_Y, d_X) = min(12, 12) = 12** —
> [[864, 12, >= 12]] unconditional (X side; Z by transpose
> duality as in §3.2), up from 11.  Both sectors now stand at 12
> with DIFFERENT engines (y: boundary-coupled links + closed
> marches + the g40 deepening; x: the mirror counting collapse +
> closed marches) — the first two-sided-instrumented member floor
> where neither sector is the straggler.  The conjectured value
> is 36; the next y-unit needs r0/r1 pushes (g >= 49/50) AND a
> bc-wall attack, the next x-unit needs the u=2 mirror stratum
> (§13.7 item 2b's pricing).

### §14.5 Specimens, falsified claims, and incidents (session 9)

- **The slip -8 saturation witness (§13.7 item 2b DISCHARGED)**:
  parented replay at gcap 29 (`specimen` lane of
  `a40_s9_slipext.py`, S7 ParentedLinkMarch + CoverFragment
  verbatim; class recomputed FROM THE FRAGMENT): seed 0, g = 29,
  slabs [1, 3, 5, 5, 8, 7], anchors [0, 0, -1, -1, -4, -9] —
  slip -8 across the gH = 8, L = 1 block is TWO consecutive steps
  each at the full -4 per-step cap; fragment E-admissible
  end-to-end, 11 cells, drift -9 (`s9_specimen_slipm8.json`).
  The left cap's tightness is now a verified object, not a table
  row.
- **FALSIFIED (same-session, mine): "K2 borrows finance small
  rightward rises for free."**  En route to R0 I reasoned that
  zero-input rises up to the borrow radius should exist.  V5
  refuted it exhaustively (max rise 0 over 3.1M pairs) and R0
  explains why (the mu-echo): the measured +1..+3 pre-transients
  are ALL input-financed; the borrow radius 3 caps the financed
  transient's reach, not a free allowance.
- **The h = 5 interface hole is STRUCTURAL and deepens** (§12.10
  item 2, §13.7 item 6): the g40 link table has crossings at
  every h in {4, 6, 7, 8, 9, 10} (min-g 19, 24, 25, 30, 35, 40)
  and NONE at h = 5 through g <= 40 — 21 cost units above the
  h = 4 optimum.  The dead shapes are exactly the two 5-slab
  patterns [1,a,b,8,c] (3-pre) and [1,a,8,b,c] (2-post), while
  2-, 4-, 5-pre all live at post 1: NOT a parity — a hole at
  exactly 5.  By contrast the L = 2 exclusion DISSOLVED at
  g = 33 (first two-slab crossings, slips -3/-4, inside the
  -4(L+1) cap): a cost phenomenon, not structure — the S7 puzzle
  splits cleanly into one theorem target (h = 5) and one closed
  measurement (L = 2).
- **Incidents**: (i) `load_closed` has a latent schema fragility —
  `s7_closed_k0_stability.json` (battery schema, no 'params')
  breaks any RERUN of the S7 assembly loader; the banked
  s7_assembly.json predates that file; the S9 shim filters
  run-shaped files (flagged for future s7 reruns).  (ii) the
  first replay-march draft (frame bookkeeping half-built) was
  scrapped unused in favor of the S7 parented machinery — dead
  code excised same session.  (iii) seed-1/seed-2 watcher +
  specimen runs all exited clean; no stray processes at close.
- (Respected: witness weights as upper bounds — the one witness
  claimed is fragment-verified; RED/AMBER/GREEN pricing before
  every run — g32 calibrated g34, g34 calibrated g40, the probe
  priced the deepening before it ran; sequential heavy runs only;
  2 GB RSS + frontier guards in code, streaming spills deleted on
  use; no SAT anywhere; every consumed vector re-verified;
  nothing re-proposed from the §9.8/§10.8/§11.7/§12.9/§13.6
  ledgers.)

### §14.6 Residue / S10

1. **The financed-half constant** (the last daylight in the
   charter's two-sided lemma): prove "net anchor rise <= 3 + #K1
   inputs" — R0's mu-argument extended to count input
   interventions (the K1 photo says each raises the frontier by
   exactly 1).  Everything below the floors is already
   input-charged by R0; this would make the wall's SLOPE a
   theorem too.
2. **The h = 5 hole**: prove the two dead 5-slab shapes from the
   E-anatomy of a minimal crossing (S7's debt-row hint).  The
   refined data (item above) makes it a two-case combinatorial
   target.
3. **The ladder-rate conjecture**: J(h+1) - J(h) = 3 quarters for
   h >= 6 (measured at h = 6..10 on real certs, +5 g per slab).
   A proof collapses the r >= 2 DP into closed form and yields an
   ANALYTIC d_Y >~ 5m/4 (~7.5r) for every b = 1 member — the L-W
   momentum budget in per-slab form, and the natural route to the
   for-all-r toroidal statement with every constant traced.
4. **The u = 2 census widening** (priced AMBER ~30-60 min: 284
   seeds byseed-streaming at hcap ~8): moves the rightward
   inequality's measured scope from u = 1 to the assembly's full
   enumerated strata.  Not consumed by any current floor.
5. **Next floor units**: y-13 needs r0 u2 g >= 51, r1 u1 g >= 50
   (both AMBER with the streaming engine) AND a bc-wall attack
   (pinch-charged loose splits or delta-resolved closure);
   x-13 needs the u = 2 mirror stratum (284 seeds, tiny trees —
   §13.7 item 2 pricing).  The W7-analog question for the BAbar
   lane stays open, floors-irrelevant.
6. **Lean**: R0 is a three-line induction over F2 Laurent
   supports (decide-shaped once finitized to a width window); the
   census emptiness tables are finite certificates under the
   §12.10 layered-march contract extended by one refinement
   lemma: streaming = a sound over-approximation of the layer
   relation (C1's node-identity is its concrete witness).

## §15 SESSION 10 (2026-09-01) — the financed half as a theorem
## (LEMMA F), the ladder's fifth rung and its asymptotic refutation,
## and the first analytic for-all-r toroidal floor (~3.6 r, both
## sectors)

Directive: prove S9's two named pieces (the financed-half constant
"net anchor rise <= 3 + #K1 inputs" and the ladder rate J' = 3q/slab)
and cash the for-all-r toroidal statement.  What happened: the
financed half is a THEOREM — but not in the charter's form.  The
per-cell count is FALSE as a local lemma (a verified zero-input
fragment raises the anchor 12 columns with no input cell in the
stretch); the true statement is LEMMA F, a one-column reading of E
that charges every non-dropping step of the pair-min to an input cell
LEFT of the frontier, with R0 as its zero-input case and a holding
cost of >= 1 input cell per 4 rows as its corollary — verified
exhaustively over the step alphabet and on every replayed census
crossing.  The ladder lands its FIFTH rung EXACTLY (g <= 46 census:
h = 11 at min_g 45 = +5 g / +3 q per slab, zero positive slips
through 46) and is then REFUTED asymptotically by a verified witness
family: the W7 species coasts at 32/7 g per slab, and a u = 2 link
[W7 seed -> light -> block -> post] exists at g = 28 (fragment-
verified), so J grows 24 quarters per 7 slabs = 3.43 q/slab beyond
the census — the "~7.5 r" extrapolation is dead.  The prize came
from a different direction: the u = 1 all-light forward tree is
FINITE cap-free in BOTH lanes (h <= 10, exact cost tables), the light-
window graphs of weight <= 2, <= 3, <= 4 are ACYCLIC in both lanes
(longest paths 5, 8, 13 slabs), and a drift-blind discharging
assembly turns these into THEOREM S10.3: d(C_{r,1}) >= f(r) with
f(r) = 11, 14, 18, 22, 25, 29, ... at r = 3, 4, 5, ..., slope 29/48
per row = 3.625 r — 2.4x the L1 floor, no closure, no ell, both
sectors, every constant traced.  Engines `a40_s10_financed.py`
(F: lemma / naive / verify lanes), `a40_s10_ladder.py` (census,
strat, bwdprobe, w7fwd), `a40_s10_w7link.py` (the backward W7
search — over-RSS, retired, §15.6), `a40_s10_forall.py` (preexact,
theorem, mirror, all2, allk, rate); data `data/a40/s10_*`;
`validate_banked` green before every lane.

### §15.1 LEMMA F (echo-charge) — the financed half, as a theorem

Setting: the y-lane solved recurrence v2[t+1] = (1+x^-1) v1[t]
+ x v1[t-3] + v2[t] + x^-3 v2[t-1]; m_t := min supp v2[t];
mu_t := min supp(v2[t] u v2[t+1]) (R0's pair-min).  Read E at the
single cell (t+2, c-3) with c = m_t:

    v2[t+2](c-3) = v1[t+1](c-3) + v1[t+1](c-2) + v1[t-2](c-4)
                   + v2[t+1](c-3) + v2[t](c),    and v2[t](c) = 1.

**LEMMA F (PROVEN, one line).**  For every t with v2[t] nonempty:
EITHER mu_{t+1} <= m_t - 3, OR an odd number of the three K1 cells
(t+1, c-3), (t+1, c-2), (t-2, c-4) is occupied.  Proof: with even
kill parity v2[t+2](c-3) = v2[t+1](c-3) + 1, so one of v2[t+1],
v2[t+2] holds c-3.  QED.  The three kill cells sit at columns
<= c-2, strictly LEFT of the v2 frontier.

Corollaries (all theorem tier):
- **(F1 = R0, with the rate)**: zero-input, mu_{t+1} <= m_t - 3
  whenever v2[t] != 0, hence two consecutive free steps drop mu by
  >= 3 — the free leftward drift is >= 1.5 columns/row.
- **(F2, holding)**: over T consecutive steps in which mu never
  drops below its start, no two consecutive steps are free, so
  >= floor(T/2) are charged; a K1 cell (tau, c') charges at most two
  steps (t = tau-1 via m_{tau-1} in {c'+2, c'+3}; t = tau+2 via
  m_{tau+2} = c'+4), so a T-row hold costs >= ceil(T/4) input cells,
  all left of the frontier.
- **(F3, window-anchor delay)**: a charging cell at (t+1, <= c-2)
  lies in the slab windows t+1..t+4, so the WINDOW anchor (v1 u v2)
  cannot exceed c-2 before slab t+5; a financed pair-min rise is
  invisible to the slab anchor for four slabs and becomes visible
  only after the killer's own descendants ((1+x^-1) debris at
  c-3/c-2, the x-echo at c-1 four rows up) are dealt with — S9's K1
  diagonal, now a consequence.

Verification (`s10_financed_lemma_naive.json`,
`s10_financed_verify.json`): the identity and F hold on ALL
16,515,072 width-6 four-row states (v1[t-2], v1[t+1], v2[t],
v2[t+1]; 227,328 even-kill states), on 832,971 random width-14
states (157,693 charged), and the zero-input two-step drop on all
4,032 width-6 pairs; row by row on the S9 slip -8 witness, on the
W7 and TC63 lifts at both ell, and on all 14 completed crossings
replayed from a parented u = 1 link march at g <= 25 (every fragment
re-verified through CoverFragment first): 28 steps, 21 drops, 14
charged, ZERO violations, longest no-drop run 2.
- **Mirror lane** (x-sector, u1[t+1] = x^-3(...)): the two-row
  pair-min drop is FALSE there — cancellation can hold the min for
  three rows (({m},{m,m+3}) -> ({m,m+3},{m}) -> ({m},{m}) -> ({m},0);
  168 of 4,032 width-6 pairs).  The mirror's R0 lives in the MAX
  gauge instead: every zero-input mirror row lies entirely in
  [min-3, max-3] of the pair (exhaustive, gauge-free).  The mirror
  financed half is listed residue (§15.7).

### §15.2 The charter's per-cell form is FALSE locally; the +3
### wall re-scoped

- **"Net anchor rise <= 3 + #K1 cells in the stretch" — REFUTED
  (verified witnesses).**  Zero-input continuation of the weight-2
  window {v2 (t-3, 0), v1 (t-2, 11)}: slabs [2,1,1,2,4,5,7], anchors
  [0, 11, 12, 12, 9, 9, 6] — the anchor rises 12 with NO input cell
  in the stretch (the seed strand exits, the right strand's x-echo
  lands at 12), E-admissible, tooth-clean.  Among 296,940 random
  windows marched zero-input, rises >= 4 occur at EVERY window
  weight 2..9 (max rises 11-12).  A hand window of weight 12 (one
  v1 kill cell at c-1 financing a pair-min jump from c+3 to c+10)
  rises 5 with zero inputs.  Mechanism: once a right strand exists,
  ONE echo-kill lets the pair-min jump onto it, and the jump is the
  strand's DISTANCE, not an input count.  The only local per-cell
  bound is the max-gauge one (dil-4 scope): from a weight-1 seed,
  A_t - A_seed <= 4 #K1 + 1.
- **What the census says instead (cost, not count)**: the pre wall
  is +3 and the whole-link wall +2 through g <= 46 (this session's
  census; 57 slip classes, ZERO positive), and on every replayed
  crossing the from-seed count A_h - A_seed <= 3 + #K1 holds with
  slack >= 2 — a COST fact (right strands and kills are paid for
  from the seed), consistent with F, implied by no per-cell lemma.
- **Re-scoping §14.6 item 1**: the financed half is F + F2 + F3;
  "the wall's slope" is not a per-cell theorem and the item is
  CLOSED in this corrected form.  (F2's ceil(T/4) is the provable
  holding cost; the measured ~8 g/row of holding rows is 4-8x above
  it — the recursion of the killers' own debris, unquantified.)

### §15.3 The ladder: fifth rung exact, then refuted asymptotically

**Census g <= 46** (`s10_slip_u1_g46.json`, `s10_link_u1k2g46.json`;
S8 engine UNCHANGED, hcap 19, dcap 30, 115.9M nodes, 1788 s,
extent/dcap truncs ZERO, regression vs the g40 tables EXACT on pre,
link, link_L and the collapse identity):

| h | min_g | J (q) | dg | dJ (q) | tight (L, dlt) |
|---|---|---|---|---|---|
| 4 | 19 | 13 | – | – | (1,-6),(1,-5) |
| 6 | 24 | 24 | – | – | (1,-2..0) |
| 7 | 25 | 31 | +1 | +7 | (1,-2..0) |
| 8 | 30 | 34 | +5 | +3 | (1,-5..-3) |
| 9 | 35 | 37 | +5 | +3 | (1,-6..-4) |
| 10 | 40 | 40 | +5 | +3 | (1,-3) |
| **11** | **45** | **43** | **+5** | **+3** | (1,-2),(1,-1) |

h = 12 is absent at g <= 46 (+5 predicts 50).  **J(h) = 3h + 10
quarters EXACTLY at h = 7..11** — five consecutive rungs, all at
L = 1, none at the cap.  Rightward slip: still EMPTY at g <= 46 (the
positive side now 27 units above the cheapest crossing); the L = 2
classes reach slip -10 (within -4(L+1) = -12); pre wall +3 and link
wall +2 unchanged.

**The asymptotic verdict: REFUTED by a verified witness family.**
(`a40_s10_ladder.py w7fwd`, `s10_w7fwd_g36.json`.)  The S8 march in
stratum u = 2, seeded ONLY at the W7 species' own weight-2 window
(the (24,7,22) lift, drift -2 per 7 rows, slabs [5,4,2,3,5,7,6] =
32 per period), finds crossings at h = 5 (g 28), 6 (g 29), 7 (g 36,
cap); the cheapest, slabs [2, 4, 7, 8, 7], g = 28, anchors
[0,-1,-1,-3,-3], is replayed (ParentedLinkMarch) and E-admissible
end to end (CoverFragment, weight 9).  Every E straddling the seed
window reads only rows at or below it, so k further W7 periods glue
BELOW the seed for every k: a u = 2 link with h = 5 + 7k, g = 28 +
32k, J = 12 + 24k quarters.  **J grows 24 q per 7 slabs = 3.43 q/slab
> 3 q/slab**: the per-slab ladder extrapolation fails in the u = 2
table for every h beyond the first crossing, and with it the
"analytic ~7.5 r" of §14.6 item 3 (the assembly's DP admits u = 2
links).  (The u = 1 table's own asymptotics is untouched: its post
coast would need a block-to-W7 transient, not constructed — the
backward search over-ran RSS, §15.6.)  Also measured: from the W7
window the cheapest 8-slab light path costs 34 = W7 itself (PRE
min_g 2, 5, 10, 14, 19, 25, 31, 34).
- **Corrected conjecture (C-W7)**: the drift-blind deficit rate of
  light coasting is 24/7 q/slab (weight 8/7 per row, the W7 rate) —
  the light-rate wall; any per-slab ladder must have J' >= 24/7.
  Under C-W7 the drift-blind member floor is ~(8/7) m = 6.86 r, and
  the conjectured 2m needs the momentum budget (drift) on top.
- **Bwd u = 1 probe** (`s10_bwdprobe_g28.json`, S6 March bwd, g <=
  28, hcap 16, 5.0M nodes): min g 1,2,3,4,7,9,11,15,17,19,25 at
  h = 1..11 — NOT exhausted (trunc_g); the +2/slab of h = 8..10 is a
  transient (+6 at h = 11).  Not consumed below (the theorem uses
  the first-minimum convention, §15.5, not the bwd table).

### §15.4 The u = 1 exhaustion, cap-free, both lanes; the light-
### window graphs are acyclic

- **Forward u = 1 light tree, cap-free** (`a40_s10_forall.py
  preexact`, `s10_forall_preexact_theorem.json`; S6 March fwd, all 8
  weight-1 seeds, gcap 600, hcap 40, extent/dcap caps 60, NO node
  truncation): 50,212 nodes, **longest light path 10 slabs**, max
  cover extent 16; at smax 4 / dil 6: 171,263 nodes, still 10 slabs,
  extent 22.  Exact table fwd_min(h) = 1, 2, 3, 4, 6, 10, 15, 22, 33,
  40 (h = 1..10); slack S := max_h [2h - fwd_min(h)] = 4 (h = 4, 5).
  (§12.4's "dies before 19 slabs" was the closed march's hcap; the
  tree in fact dies at 10.)
- **Mirror lane** (`mirror` lane, MClosedMarch nHmax 0, gcap 600, m
  40, both extent semantics 14 / 34): 121,421 nodes, **longest 9
  slabs**, fwd_min 1, 2, 3, 4, 6, 10, 15, 22, 33 (extent 34; 26 at
  h = 8 under extent 14), S = 4 — §13.7 item 3 ANSWERED: the mirror
  all-light u = 1 stratum exhausts cap-free like the y-lane's.
- **Light-window graphs** (`allk` lane, `s10_forall_allk*.json`):
  nodes = translation-normalized 4-row windows of weight <= k, edges
  = the forced row plus <= 3 input cells in dil-4 keeping the next
  window <= k, tooth rule on.  ACYCLIC in both lanes for k = 2, 3, 4
  with longest paths P = 5, 8, 13 slabs (y: 340 / 13,703 / 543,934
  nodes; mirror: 340 / 13,573 / 530,821).  For k = 2, 3 the start set
  was widened to EVERY window of weight <= k in an 8 x (4k+4) box
  (no internal connectivity assumed; 783 / 65,713 nodes): P
  unchanged.  Consequence (scope: growth dil-4/smax-3 + tooth rule):
  **any 6 consecutive slabs of any walk contain a slab >= 3, any 9 a
  slab >= 4, any 14 a slab >= 5** — in both lanes.  (k = 5 is ~13M
  windows: not attempted.)

### §15.5 THEOREM S10.3 — the for-all-r toroidal floor (drift-blind)

Sigma_j W_j = 4|v| over the m slabs of a class-minimal doubly-
spanning walk (each row lies in 4 slabs).  Discharging: write
4|v| = 2m + Sigma_j (W_j - 2).  Partition the slabs into heavy slabs
(W >= 8: excess >= 6) and maximal light runs; a run's SEED is its
lowest minimum-weight slab, of weight u:
- slabs of the run BELOW the seed have weight >= u+1 >= 2; slabs of a
  u >= 2 run have weight >= 2: excess >= 0, and on any such stretch of
  h slabs the window facts of §15.4 force excess >= pm(h) (the exact
  minimum over the automaton of "slabs since the last >= 3 / >= 4 /
  >= 5", window lengths 6 / 9 / 14: pm = 0,0,0,0,0,1,1,1,2,2,2,3,3,4,5,
  5,5,6,6 at h = 1..19; min-mean cycle 5/12 per slab, i.e. **every
  sustained light coast weighs >= 2 + 5/12 = 2.417 per slab**);
- the slabs from a u = 1 seed upward form a light path from a weight-
  1 window: at most 10 slabs, excess >= -S = -4 (exact table);
- every u = 1 run is followed by a heavy block, so #u1-runs <= B, and
  each block's excess >= 6 > 4: more blocks never help the adversary;
- all-light u = 1 is impossible for m + 1 > 10 (a closed light walk is
  a light path of m + 1 slabs from its seed); all-heavy >= 8m.
Hence 4|v| >= min( 2m + pm(m) [all-light, u >= 2, open-path bound],
min_{h_f <= 10} [2m + 6 - (2h_f - fwd_min(h_f)) + pm(m - 1 - h_f)] )
and, for the family (m = 6r), with identical constants in both
lanes (`rate` lane, `s10_forall_rate.json`):

| r | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 | 11 | 12 | 40 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| f(r) | 7 | **11** | 14 | 18 | 22 | 25 | 29 | 32 | 36 | 40 | 43 | 145 |
| L1 | 3 | 5 | 6 | 8 | 9 | 11 | 12 | 14 | 15 | 17 | 18 | 60 |

> **THEOREM S10.3 (scope-listed).  d(C_{r,1}) >= f(r) for every
> r >= 2, both sectors (d_Y and d_X separately), with f(r) as
> tabulated and f(r) = (29/48) m + O(1) = 3.625 r asymptotically.**
> Tiers: the arithmetic and the automaton are theorem; the
> exhaustion tables and the window-graph acyclicity are certificates
> under the growth scope (dil-4 / smax-3 inputs, stability-checked
> at 4 / 6 for the tree), the tooth rule (a class-minimality prune),
> and wrap-free cover semantics (the objects have extent <= 22, so
> literal for ell >= 24, i.e. r >= 3; r = 2 is the exact 24 anyway).
> No closure, no drift, no ell enters: it is a statement about every
> class-minimal y-spanning (resp. x-spanning) walk of m rows.

Controls (all PASS): the ell-free floor at m = 12 is 7 <= 18 (the
b = 0 (12,12) witness is admitted — its -6 is invisible to a drift-
blind bound, which is the honest reason it never binds here) and
7 <= 24 at (18,12); at m = 6 it is 2 <= 12 (r = 1 stays with the
exact values); at (24,18) it is **11**, below the S9 certified 12 as
it must be — and only one unit below it, with no closed march and no
closure: the discharging assembly recovers 11/12 of the S7-S9
machinery's value at m = 18 from tables that cost seconds.

The conditional forms, for the record: under C-W7 (§15.3) the same
partition gives 4|v| >= (32/7) m - O(1), i.e. ~6.86 r; the ladder's
7.5 r is unavailable even conditionally; 12 r needs the momentum
budget on top of a light-rate wall.  Stage 4 (the h = 5 hole) was
NOT reached.

### §15.6 Falsified claims and incidents (session 10)

- **FALSIFIED (the charter's, re-scoped)**: "net anchor rise <= 3 +
  #K1 inputs" as a local lemma — verified zero-input counterexamples
  at every window weight >= 2 (§15.2).  The theorem in its place is
  F / F2 / F3.
- **FALSIFIED (same-session, mine)**: "the mirror lane obeys the
  same two-row pair-min drop" — exhaustive refutation (168 holding
  pairs); the mirror R0 is a max-gauge statement (§15.1).
- **FALSIFIED (§14.6 item 3's extrapolation)**: "J' = 3 q/slab for
  all h" — the W7-coast link family gives 24/7 q/slab in the u = 2
  table (§15.3).  The census-range statement (h <= 11) is exact.
- **Incidents**: (i) the backward W7 search (`a40_s10_w7link.py`,
  Dijkstra with a parents dict) reached 3.5 GB RSS with < 200k pops —
  the heavy-entry fan-out from a W7 window (~10^4 children per node)
  outran a pop-count-gated RSS check; killed by hand, no data
  consumed, the forward S8-engine form (`w7fwd`, 16 s) replaced it.
  (ii) the first census launch failed on argparse (a top-level
  `--log` placed after the subcommand); shell redirection is
  unavailable in this environment, so every long lane now tees its
  own log.  (iii) `experiments/bb_lab/data` is git-ignored: S10 data
  files are force-added like their predecessors.  (iv) the exact
  cyclic automaton minimum (756 start states x m steps) was too slow
  at r = 40 and was replaced by the open-path bound (sound, at most
  a few quarters weaker).  (v) the parented replay lane peaked at
  2.27 GB RSS on one seed at g <= 25 (under the 3 GB rule; noted for
  any deeper replay).
- (Respected: witness weights as upper bounds — the W7 crossing and
  the naive-form fragments are verified objects; RED/AMBER/GREEN
  pricing — the g46 census was priced off the g40 growth (1.6x nodes,
  ran 3.2x wall under a concurrent sibling racer); sequential heavy
  runs; RSS <= 3 GB per process with one over-run killed; no /tmp; no
  SAT; every consumed vector re-verified; nothing re-proposed from
  the §9.8/§10.8/§11.7/§12/§13.6/§14.5 ledgers or A42's.)

### §15.7 Residue / S11

1. **The u = 1 post-coast asymptotics**: a block-to-W7 transient
   (a u = 1 link whose post phase is a W7 coast) — the S8 forward
   engine seeded at u = 1 with a W7-window target check in POST, or
   a LAYERED backward march from W7 (the Dijkstra form is out).
   Decides whether the u = 1 table's own ladder bends at 24/7 too.
2. **The light-rate wall (C-W7)**: the k = 5 window graph (~13M
   windows) is the next discharging rung (+4 per 14+ slabs would
   lift the slope); the full min-mean cycle over <= 7 windows is
   the wall itself (2e8 windows, §10.1's non-materializable core).
   Any proven slope above 29/48 lifts f(r) for every r at once.
3. **The momentum budget on top of S10.3**: f(r) uses no drift;
   the census's empty rightward side (g <= 46) and F2's holding cost
   are the ingredients for a per-column rightward cost c; with the
   W7 coast drifting -2 per 7 rows, c >= ~3 per column would push
   the sustained rate above 2 per row.  The quantity to prove is c.
4. **The mirror financed half** (max-gauge form of F with the u2
   kill set), and the mirror all-k check under the extent-14 sector
   semantics (done here at extent 34 and 14 for the tree; the window
   graphs used the cover).
5. **The h = 5 hole** (Stage 4, not reached): with the g46 table the
   two dead 5-slab shapes stand 27 units above the h = 4 optimum.
6. **Lean**: F is a single-cell E reading (decide over a 4 x 8 window
   after finitization); the exhaustion tables and the acyclicity
   certificates are finite DFS certificates; S10.3's arithmetic is a
   sum-rearrangement lemma over a slab partition.

## §16 SESSION 11 (2026-09-01) — THE COMPARISON THEOREM: can
## toroidal minima beat windowed minima at b = 1?

Directive: decide whether, at every b = 1 member, every doubly-
spanning nontrivial X-logical weighs at least the member's windowed-
sector minimum — the shortcut that would reduce the b = 1 lower half
to Theorem W alone.  Load-bearing negative control: the statement is
FALSE at b = 0 ((12,12): the a36 witness is doubly-spanning at 18
while floor_cyl(12) = 24; (6,6): d = 6 while floor_cyl(6) = 12), so
any proof must consume b = 1 arithmetic, and any sweep must reproduce
the b = 0 violation.  Namespace: `a40_s11_*.py`, `data/a40/s11_*`.

### §16.0 Stage 0 — the pinned statement (before any computation)

**Setting.**  Member (ℓ, m) = (6(r+b), 6r), Λ := ⟨(ℓ,0), (0,m)⟩ ⊂ Z²,
T := Z²/Λ.  For u ∈ Λ \ {0} write C_u := Z²/⟨u⟩ (a cylinder; by a
unimodular change of basis C_u ≅ Z × Z_{g(u)}, g(u) = gcd of the
coordinates of u, and the fixed pair transports to a pair on it —
the S5 twisted-atlas reduction, §10.4) and π_u : C_u → T the covering.

**Definition (u-compact).**  A torus X-cycle v is *u-compact* if
v = π_{u*}(w) for a compactly supported X-cycle w of C_u with
|w| = |v| (the support of w maps bijectively onto that of v: an
injective compact lift).  Sectors of the nontrivial X-logicals of T:
- **windowed** W := W_x ∪ W_y, W_x := (0,m)-compact, W_y :=
  (ℓ,0)-compact.  (Lemma K's gap-defined windowed set — cyclic
  x-gap ≥ 4, resp. y-gap ≥ 4 — is CONTAINED in W: the gap supplies
  the lift.  W is the right boundary because Theorem W's proof uses
  only the existence of the compact lift, so **Theorem W bounds all
  of W_x by floor_cyl(m)**, and its θ′-mirror bounds W_y by
  floor_cyl'(ℓ).)
- **toroidal** D := the complement of W among nontrivial logicals
  ("doubly-spanning": no injective compact lift to either straight
  cylinder).  Sub-sectors: **helical** H := u-compact for some
  u ∈ Λ outside the two axes; **2D** := u-compact for NO u ∈ Λ.
- min_W(ℓ,m) := min{|v| : v ∈ W}, min_D(ℓ,m) := min{|v| : v ∈ D}.

**COMPARISON STATEMENT (global form, the target).**  For every b = 1
member (6r+6, 6r), r ≥ 1:   min_D(6r+6, 6r) ≥ min_W(6r+6, 6r).

**Class-wise form (the stronger probe).**  H₁(T) = W_x^cl ⊕ W_y^cl
where W_x^cl := span of the classes of x-windowed cycles (support in
an x-window of width ℓ − 4) and W_y^cl likewise, and for every
nonzero class c, min over c ≥ min_W with the minimum of a class in
W_x^cl attained in W_x (resp. W_y^cl in W_y).  The decomposition is
a linear-algebra fact per member (tested below); the attainment
half is what the b = 0 witness must violate.

**The norm reformulation (Lemma U + the wrapping norm).**  For
u ∈ Z² \ {0} let N(u) := min weight of a nontrivial compact X-cycle
of C_u (∞ if H₁(C_u) = 0).  N((0,p)) = floor_cyl(p); N((p,0)) = the
mirror's floor_cyl'(p); N((3,6)) = 10 (TC63, §10.4).
> **LEMMA U (free unwinding; proven, one line).**  If v is a
> nontrivial X-logical of T and v is u-compact, then |v| ≥ N(u).
> Proof: the lift w is a compact cycle of C_u; a compact trivializer
> of w would push forward to a trivializer of v; so w is nontrivial
> on C_u and |v| = |w| ≥ N(u).  ∎
Hence the comparison statement at a b = 1 member is IMPLIED by
  (H1) N(u) ≥ min_W(ℓ,m) for every u ∈ Λ \ (Z(0,m) ∪ Z(ℓ,0)), and
  (H2) every 2D nontrivial logical weighs ≥ min_W(ℓ,m);
and it is EQUIVALENT to (H1′) ∧ (H2) with (H1′) the same inequality
restricted to the u realized by class-minimal helical logicals.  The
b = 0 control must appear as a FAILURE of (H1) or (H2): the question
"is the (12,12) witness u-compact for some u, and which" is the
first computation of the session (if it is (12,12)-compact the
discount is a diagonal-norm fact, N((12,12)) ≤ 18 < 24 = N((0,12)),
and the b-bit is the period of the diagonal cylinders: g((m,m)) = m
at b = 0 versus g((aℓ, bm)) = 6·gcd(a(r+1), br) at b = 1).

**What the session will NOT claim.**  Theorem W's conditionality
(floor_cyl(m) = 12r certified only at r ≤ 2, pure half to p = 21) is
inherited verbatim; every "min_W = 12r" below is conditional on it
unless r ≤ 2.  At r ≤ 2 the comparison is TRIVIALLY true (d = 12r is
certified, so min_D ≥ d = min_W), and the empirical content at
r ≤ 2 is the EQUALITY question (is min_D = min_W attained?) plus the
b = 0 violations; the theorem's content lives at r ≥ 3, where no
b = 1 population is banked and only the floors of §13–15 exist.

### §16.1 The structure theorems (hand proofs; the machine checks
### are in §16.2)

Fix a nontrivial X-logical v of T = Z²/Λ, connected in the footprint
sense (two cells adjacent iff one Z-check reads both; class-minimal
logicals are connected by the S6 connectivity lemma, §11.2).  Let
V ⊂ Z² × {blk} be its plane lift (Λ-periodic), C one footprint-
component of V, and Π(v) := {λ ∈ Λ : C + λ = C} the PERIOD LATTICE
of v (well defined up to conjugation; the components of V are the
Λ/Π-translates of C).

> **THEOREM T (trichotomy by period rank).**  (i) rank Π(v) ≥ 1.
> (ii) rank Π(v) = 1, Π = Zu, iff v is u-compact; then C/Π ⊂ C_u is
> the injective compact lift.  (iii) rank Π(v) = 2: Π has finite
> index in Λ, C/Π is a cycle of the cover torus T′ := Z²/Π (a
> covering of T of degree [Λ:Π]) projecting bijectively onto v and
> nontrivial there, so |v| ≥ d(T′); Π = Λ iff V is connected
> ("fully 2D").
> *Proof.*  Adjacency lifts through coverings (a check's footprint
> is a footprint upstairs), so π(C) is closed under adjacency in v
> and π(C) = v; and c₁, c₂ ∈ C with c₂ = c₁ + λ force C + λ = C, so
> π|_C is injective modulo Π: C/Π → v is a bijection.  (i) If Π = 0,
> C is a finite plane cycle (a component of the cycle V splits off
> by the connectivity lemma), hence a plane boundary by the K2
> regularity certificate (§8.1), and v = π_*(C) is a torus boundary
> — contradiction.  (ii) With Π = Zu, C/Π is a compact cycle of C_u
> bijective onto v (u-compact); conversely an injective compact
> lift S̃ has plane preimage whose components have stabilizers in
> Zu, nonzero by (i).  (iii) C/Π ⊂ T′ is bijective onto v and a
> cycle; a trivializer on T′ would push down to one of v.  ∎

So D = {rank 2} ∪ {rank 1 off the axes}, and the fully-2D objects
are those whose plane lift is one connected net.  For a twisted
u = (aℓ, bm) with gcd(a,b) = 1, g(u) = gcd(|a|ℓ, |b|m); a rank-1
object in direction u must connect c to c + u by footprint steps
(each ≤ 4 in both coordinates), so |v| ≥ max(|a|ℓ, |b|m)/4: **at a
member, helices of weight < 12r have |a|, |b| ≤ 7 — (H1) is a FINITE
computation per member** (N(u) over ≤ ~70 directions).

> **PROPOSITION O (the cut obstruction is the class).**  For any
> u ∈ Λ and any torus cycle v, the following are equivalent: (a) v
> = π_{u*}(w) for SOME compact cycle w of C_u (of any weight); (b)
> [v] ∈ W_u^cl := im(π_{u*} : H₁(C_u) → H₁(T)).  Cut form: for any
> fundamental-domain chain w_F of the deck action on C_u (|w_F| =
> |v|), syn(w_F) = (1 + τ)·a with τ the deck generator and a
> supported at one cut, and every compact w with π_* w = v is
> w_F + (1 + τ) z with syn(z) = a; the obstruction is the class of
> a modulo compact syndromes, which vanishes iff (b).
> *Proof.*  (b)⇒(a): [v] = [π_* w₀] gives v = π_* w₀ + ∂ζ =
> π_*(w₀ + ∂ζ̃) for any lift ζ̃ of ζ.  (a)⇒(b) trivial.  The kernel
> of π_* on compact chains is (1 + τ)·(compact chains) (reduction
> modulo the deck), which gives the cut form.  ∎

Define the **u-lift weight** λ_u(v) := min{|w| : w compact cycle of
C_u, π_{u*} w = v} (∞ unless [v] ∈ W_u^cl) and the **closure excess**
e_u(v) := (λ_u(v) − |v|)/2 ∈ Z_{≥0} (parity: |w_F + (1+τ)z| ≡ |v|
mod 2); e_u(v) = 0 iff v is u-compact.

> **PROPOSITION E (surgery inequality).**  If [v] ∈ W_u^cl then
> |v| ≥ N(u) − 2 e_u(v).  *Proof.*  The minimizing w is a compact
> cycle of C_u pushing to the nontrivial v, hence nontrivial on C_u
> (Lemma U's argument), so |v| + 2e_u = λ_u ≥ N(u).  ∎

**The class decomposition.**  H₁(T) = W_x^cl ⊕ W_y^cl at every
member tested (r ≤ 4, both columns, §16.2), with both summands
6-dimensional; so every class is c = c_x + c_y uniquely, and by
Proposition O: x-surgery is available exactly on the pure-x classes
(c_y = 0), y-surgery exactly on the pure-y classes, and NEITHER
straight surgery exists on a mixed class (c_x ≠ 0 ≠ c_y).  For a
mixed class only twisted directions u with c ∈ W_u^cl can serve.

**The b = 1 cushion (the first place the arithmetic enters).**  At
a b = 1 member the two straight floors DIFFER: floor_cyl(m) = 12r
versus floor_cyl′(ℓ) = 12r + 12 (θ′-mirror of floor_cyl at period
ℓ = 6r + 6; both conditional as in Theorem W).  Hence on a pure-y
class, Proposition E gives |v| ≥ 12r + 12 − 2e_y(v) ≥ min_W whenever
e_y(v) ≤ 6 — a 12-unit cushion that does not exist at b = 0, where
ℓ = m and the two floors coincide.  On a pure-x class the cushion is
0 and the comparison needs e_x(v) = 0 (i.e. a windowed minimizer) or
an independent argument.  This is the exact shape of the "surgery
accounting" the charter asked for: the discount term is 2e_u(v),
and what b = 1 buys is the mismatch of the two windowed floors plus
the thinness (g = 6·gcd(a(r+1), br)) of every twisted cylinder.

### §16.2 Stage 1 — the empirical verdict (`a40_s11_compare.py
### classify / surgery / hunt / norm`, `s11_*.json`; validate_banked
### green before every lane)

The classifier decides u-compactness EXACTLY: an injective lift is
searched among the preimages in an X-window of C_u wide enough to
hold every connected compact lift (|v| × the maximal footprint
X-step, per component), so a positive verdict is a verified compact
cycle of C_u projecting bijectively onto v and a negative verdict is
solver-UNSAT over an exhaustive window (descriptive labels; no floor
consumes a negative verdict).  Gap sectors (Lemma K) are computed
alongside; conventions (READ/BND offset tables) are asserted against
the code's own H_Z/H_X rows.  Directions tested: 16 (all
gcd(a,b) = 1 with |a|, |b| ≤ 3, up to sign).

| frame | population (all re-verified: cycle, non-stab, class) | sector minima (n) | H₁ = W_x ⊕ W_y |
|---|---|---|---|
| (6,6) = member (1,0), b = 0 | all 84 nontrivial w ≤ 6 (census-complete) | **2D: 6 (84)**; W_x, W_y, helical: none | (window def. vacuous at ℓ = 6; cylinder images 6+6 ✓) |
| (12,12) = member (2,0), b = 0 | a36 witness + 48 (12;4,4)-w6 shear pullbacks = 49 distinct w18 | **2D: 18 (49)**; nothing else | ✓ 6 ⊕ 6 |
| (12,6) = member (1,1), gross | all 1,884 nontrivial w ≤ 12 (census-complete, 246 classes) | **2D: 12 (1,092)**, W_x: 12 (792); W_y, helical: none | ✓ (W_y window-vacuous at m = 6; cylinder images 6+6) |
| (18,12) = member (2,1) | the L12×2 τ₀-witness (w24) | W_x: 24 (1); 2D at 24: **OPEN** (SAT hunts timed out, below) | ✓ 6 ⊕ 6 |
| (24,18), (18,18), (30,24) | no population | — | ✓ 6 ⊕ 6 each |

**The b = 0 control (REQUIRED, reproduced twice).**  At (12,12) the
witness and all 48 shear pullbacks are 2D (no injective compact lift
in any of 16 directions; gap-dense both axes) at weight 18 <
24 = floor_cyl(12) = min_W: the comparison statement FAILS at b = 0
as it must.  At (6,6) all 84 minima are 2D at 6 < 12 = min_W: a
second b = 0 violation at r = 1.  Both violations are in the 2D
sector — NOT helical: the diagonal-helix reading of §16.0 is
REFUTED (the witness is not (12,±12)-compact, and N((12,±12)) ≥ 12,
N((4,4)) ≥ 12, N((6,6)) ≥ 12 at Wcap 11 by the twisted atlas).

**b = 1, r = 1 (gross): equality, massively.**  min_D = 12 = min_W:
1,092 of the 1,884 minimum-weight logicals are 2D, and 201 of the
246 classes that contain a weight-12 logical have ONLY 2D minimizers
(45 have windowed ones).  The class-wise "attainment" form of §16.0
is therefore FALSE already at r = 1: the comparison holds globally
with equality, not class by class.  Consistency: the SAT hunt for
gap-dense nontrivial logicals of weight ≤ 12 at gross finds them
immediately (12 translation-orbits, all 2D on re-classification).

**Helices do not exist in any banked population**: zero objects in
W_y, zero helical, at every frame — every non-windowed minimizer is
fully 2D.  Twisted-cylinder floors (`norm` lane, Wcap 11, compact
triviality by the generator march, TC63 control N((3,6)) = 10
reproduced): N(u) ≥ 12 for u ∈ {(4,4), (6,6), (12,6), (12,−6),
(24,−6), (12,12), (12,−12)}; (24,6) exceeded the 4M-state cap.
Cost: (12,6)-lattice directions at Wcap 11 run 3–14 s; the g = 12
directions 200–290 s; Wcap 23 (the r = 2 helix question) is RED for
this engine.

**The (18,12) equality question — OPEN.**  SAT existence hunts for
gap-dense nontrivial logicals (pycryptosat, XOR-native, cardinality
by sequential counter; nontriviality as an OR over the 12 pairings;
witnesses verified end to end, UNSAT never consumed): (12,12) ≤ 18
finds the witness family in 0.2 s and (12,12) ≤ 17 / (18,12) ≤ 23
time out (300 / 600 s) — the solver cannot certify emptiness where
d is already certified, so its silence at (18,12) ≤ 24 (1,500 s)
and in the mixed-class-restricted rerun (§16.3; 1,800 s) is
uninformative.  Whether min_D(18,12) = 24 or > 24 is listed residue
(§16.6); the L12-stack minimum is x-windowed and its class is pure-x.

### §16.3 THEOREM P — the parity law of cylinder images, and the
### mixed sector (`a40_s11_structure.py`, `s11_parity_law.json`,
### `s11_pullback.json`)

Computing W_u^cl = im(π_{u*}) by exact linear algebra (kernel of the
cylinder syndrome on a K-period X-window, pushed forward; stable
between K = 4 and K = 6) for 22 directions u = (aℓ, bm) at (6,6),
(12,6), (12,12), (18,12), (24,18):

> **THEOREM P (parity law; certificate at the five frames, proof
> sketch below).**  W_u^cl depends on u only through (a mod 2,
> b mod 2).  The three subspaces W_x := W_{(0,1)}, W_y := W_{(1,0)},
> W_d := W_{(1,1)} are each 6-dimensional and PAIRWISE COMPLEMENTARY
> (every pair spans H₁(T)); W_d is the graph of an isomorphism
> W_x → W_y (every basis vector of W_d has both components nonzero).
> Hence H₁(T) ≅ F₄ ⊗_{F₂} U with U 6-dimensional and W_{(a,b)} =
> (a + bω) ⊗ U: the 3·63 = 189 nonzero classes with a compact
> cylinder representative in SOME direction are exactly the rank-one
> tensors, and the remaining 4095 − 189 = **3,906 classes (95.4%)
> are MIXED: no compact cylinder representative in ANY direction, so
> no surgery of any kind (straight or twisted) exists for them.**
> *Why (the structural reason, with one measured input).*  π_u :
> C_u → T factors through the index-2 cover T_u := Z²/⟨u, 2Λ⟩,
> which depends only on u mod 2Λ; so W_u^cl ⊂ im(τ_u), the transfer
> from T_u, whose rank is k/2 = 6 by the A35 sheet-SES (universal).
> Equality by dimension.  Measured (same script): for each of the
> three Z₂-covers T_x = (2ℓ, m), T_y = (ℓ, 2m), T_d =
> Z²/⟨(ℓ,m),(2ℓ,0)⟩ the pullback ρ^* : H₁(T) → H₁(T_u) has rank 6 and
> **ker ρ_u^* = W_u EXACTLY** at (12,6), (12,12), (18,12) — the
> sheet-SES with σ_* = id in the form ρ^*τ_* = 1 + σ_* = 0.

> **COROLLARY P1 (mixed classes are seen by every Z₂-cover).**  If
> [v] is mixed, ρ_u^* v (weight 2|v|) is a nontrivial logical of
> EACH of T_x, T_y, T_d; hence |v| ≥ ½ max(d(T_x), d(T_y), d(T_d)).
> If [v] is rank one of parity p, it dies on T_p and is nontrivial on
> the other two covers.  (The (12,12) witness: pullbacks of weight 36,
> nontrivial on all three covers — verified.)

**Where the b = 0 discount lives — decided.**  Every cheap 2D object
found sits in a MIXED class: all 84 minima of (6,6) (84 mixed
classes), all 49 witnesses of (12,12) (13 mixed classes), and the
full gross census splits (`s11_class_kinds.json`, all 1,884 minima):

| gross w12 minima | sector | classes |
|---|---|---|
| 1,008 | 2D, MIXED class | 153 |
| 84 | 2D, pure-y class | 48 |
| 792 | W_x (x-windowed), pure-x class | 45 |
| 0 | any object in a pure-d (diagonal-parity) class | 0 of 63 |

The pure-y 2D minima have λ_y(v) > 17 at the 120-s solver cap
(closure excess e_y ≥ 3; by Prop E with floor_cyl′(12) = 24 they
need e_y ≥ 6 — the cushion is fully consumed at r = 1), and the 63
diagonal-parity classes have NO logical at weight 12.  Their exact
minima (`a40_s11_structure.py diag`, `s11_diag_classes.json`; walk
kernel, census-complete per class): **48 of the 63 pure-d classes
have minimum EXACTLY 16, the other 15 have minimum ≥ 18** — at
r = 1 the (odd, odd) parity is the "long" one, 4 above the member
distance, the first quantitative sign that diagonal directions are
expensive at b = 1 (compare b = 0, where the mixed/diagonal witness
is the CHEAP one).  So
(H2) fails at b = 0 inside the mixed sector, where NO cylinder floor
applies, and the comparison theorem at b = 1 is EQUIVALENT to:
  (M) every class-minimal logical of a MIXED class of the b = 1
      member weighs ≥ min_W = 12r,   plus
  (R1) the rank-one sector: Prop E with (H1) [finite per member] and
      the excess control.
At r = 1, (M) holds with equality (1,092 objects).  At r = 2 it holds
by d = 24 (certified) and the equality question is open.  At r ≥ 3
it is the conjecture's lower half, restated: **the shortcut does not
bypass the closure-aware program — it relocates it exactly to the
mixed sector**, the 95.4% of classes invisible to every cylinder.

### §16.4 Stage 2 — the surgery attempt: outcome and the precise
### obstruction

What is PROVEN (hand, with the certificate inputs named):
1. **Theorem T** (trichotomy: windowed / helical / 2D by the rank of
   the period lattice), **Lemma U** (u-compact ⟹ |v| ≥ N(u)),
   **Prop O** (the cut obstruction is the class: surgery in direction
   u exists iff [v] ∈ W_u^cl), **Prop E** (|v| ≥ N(u) − 2e_u(v)).
2. **Theorem P** at the five frames (certificate) with the transfer
   factorization as its mechanism; **Corollary P1**.
3. **(H1) is finite** per member (helices below 12r have
   |a|, |b| ≤ 7), and holds at r = 1 for every direction tested
   (N(u) ≥ 12 = min_W); at r ≤ 2 the whole comparison is implied by
   the certified d.
4. **The b = 1 cushion**: on pure-y classes |v| ≥ 12r + 12 − 2e_y(v),
   conditional on floor_cyl′(6r+6) = 12r + 12.

What the surgery CANNOT do, precisely: a cut-and-close surgery in
direction u is available exactly on the 63 classes of W_u^cl
(Prop O), so at most 189 classes admit any surgery at all; the
b = 0 discount is realized in the 3,906 MIXED classes (witness at
(12,12); all of (6,6)); and on the rank-one classes the surgery
inequality carries the excess −2e_u(v) which is NOT small (e_y ≥ 3
at gross, needing ≥ 6 for the pure-y floor to be met), so even there
"e = 0 for class-minimal 2D objects" is false in general and the
comparison on rank-one classes needs the cushion (pure-y, diagonal:
N ≥ 12r + 12 conjecturally) rather than the excess.  Largest special
case that falls: **rank-one classes of parity ≠ x at b = 1 with
excess ≤ 6 (pure-y) — and every helical object at r = 1** (N(u) ≥ 12
for all tested u).  The obstruction, in one sentence: *the b = 1
comparison theorem is the statement that MIXED classes of the member
have no logical below 2m, and mixed classes are exactly those that
every Z₂-descent detects but no cylinder sees — the doubling-deficit
problem (A38's F2b lane) in its cleanest form: |v| ≥ ½ d(T_u) for
all three covers T_u, and d(T_u) itself is a torus distance.*  The
b-bit enters through the three covers of a b = 1 member: T_x =
(12r+12, 6r), T_y = (6r+6, 12r), T_d = Z²/⟨(6r+6, 6r), (12r+12, 0)⟩,
whose shortest lattice vectors in the wrapping norm are all "long"
(no diagonal (m,m) direction of period m exists, unlike b = 0 where
T_d contains (m,m) with N((m,m)) ≤ 2m − 6 realized by the witness
species).

### §16.5 Falsified claims and incidents (session 11)

- **FALSIFIED (mine, §16.0's reading)**: "the b = 0 witness is a
  diagonal helix, N((12,12)) ≤ 18" — the witness is fully 2D (no
  injective compact lift in 16 directions), its class is MIXED
  (outside all three cylinder images, pulls back nontrivially to
  all three Z₂-covers), and N((12,±12)) ≥ 12 at Wcap 11 with no
  weight-≤ 11 nontrivial compact cycle; the discount is a
  mixed-class phenomenon, not a twisted-cylinder one.
- **FALSIFIED (the class-wise attainment form of §16.0)**: at gross,
  201 of 246 populated classes have ONLY 2D minimizers at weight 12.
  The comparison holds globally with equality, never class by class.
- **NOT a floor: SAT emptiness.**  Every UNSAT/timeout here is an
  observation; the (12,12) ≤ 17 and (18,12) ≤ 23 controls (where d
  is certified) TIMED OUT, which calibrates the (18,12) ≤ 24 silence
  as uninformative.  Witnesses found (gross w12 gap-dense ×12; the
  (12,12) witness family) are verified objects.
- **Observation, not claim**: the closure excess is not universally
  small (e_y ≥ 3 at gross pure-y 2D minima, solver-observed at a
  120-s cap); Prop E's discount term cannot be bounded by the borrow
  radius.
- **Incidents**: (i) the (18,6) ≤ 22 checkpoint files that
  `a40_s3_l1_gates.py` reads (`tdg432/ckpt_W22_{ntrv1,seam1}.jsonl`)
  are ABSENT from `data/a40/tdg432/` (only sweeps/rungs are banked) —
  the (18,6) population could not be re-audited this session and the
  S3 gate script no longer runs as written (flagged; n = 216 exceeds
  the walk kernel, so it is not recomputable in-session).  (ii) the
  twisted-atlas direction (24,6) exceeded the 4M-state cap (RSS
  guard never tripped; peak < 1 GB); the g = 12 directions cost
  200–290 s at Wcap 11 — Wcap 23 is RED for this engine.  (iii)
  shell heredocs and `&` are blocked in this environment; all lanes
  tee their own logs (S10 rule) and run via the harness's
  background mode.  (iv) log monitors time out silently at 5 min
  regardless of the requested timeout; polling by `until` loops
  replaced them.
- (Respected: validate_banked green before every lane; every
  consumed vector re-verified end to end — cycle, non-stabilizer,
  class, weight; witness weights as upper bounds only; no SAT on any
  certificate path — SAT appears only as a witness hunter and as the
  exact decision procedure for the DESCRIPTIVE u-compactness labels,
  never in a floor; RED/AMBER/GREEN pricing (the norm lane priced
  Wcap 11 off the S5 twisted atlas; Wcap 23 declared RED without
  launching); sequential heavy runs, RSS guard 2.5 GB in code, no
  /tmp; nothing re-proposed from the §9.8/§10.8/§11.7/§12.9/§13.6/
  §14.5/§15.6 ledgers or A42's.)

### §16.6 Residue / S12

1. **The (18,12) equality question**: does a MIXED class of (18,12)
   contain a weight-24 logical?  A certificate route: the tdg432
   descent census restricted to the 3,906 mixed classes at W = 24
   (the mixed classes are those outside ker ρ^* for all three
   covers — computable from the rung data), or the S7 closed march at
   (m,ℓ) = (12,18) with g-cap 96 restricted to all-heavy-8 profiles
   (the L12 stack's profile) — both priced, neither run.
2. **Theorem P as a theorem**: the transfer factorization + rank
   τ* = 6 + σ_* = id gives W_u^cl = ker ρ_u^*; write it out and Lean
   it (finite linear algebra per member; the ∀r statement needs the
   A35 sheet-SES).
3. **(H1) at r = 2, 3**: N(u) for u = (18a, 12b), (24a, 18b),
   |a|, |b| ≤ 7, needs the corridor-racer engine (A42 §2.9) on the
   TRANSPORTED pair at period g(u) ∈ {6, 12, 18, …}; the S4 automaton
   is RED beyond Wcap ~13 at g = 6.
4. **The mixed-sector floor (M)** — the whole remaining gap, now
   with a tool: Corollary P1's three-cover descent, |v| ≥ ½ d(T_u),
   iterated (the pullback class on T_u is again rank-one or mixed
   there) — the A38 deficit machinery applies verbatim, and the
   b = 1 arithmetic is visible as the absence of any period-m
   diagonal in the covers' lattices.  The sharpest question: **at a
   b = 1 member, is every mixed class's minimum ≥ 2m?**  (True at
   r = 1 with equality, r = 2 by certificate; the conjecture at
   r ≥ 3.)
5. **Lean**: Theorem T / Lemma U / Prop O / Prop E are finite-free
   hand lemmas over the covering-space formalism (a natural
   companion to the BBCover layer); Theorem P's per-member instance
   is a rank computation.

**Session-11 state (for the program map).**  The comparison
statement min_D ≥ min_W at b = 1 is TRUE at r ≤ 2 (by the certified
d; with equality at r = 1 through 1,092 two-dimensional minimizers
in 201 classes) and at r ≥ 3 is EQUIVALENT to the mixed-sector floor
(M): every class outside the three transfer images W_x, W_y, W_d
(3,906 of 4,095 classes) has minimum ≥ 2m.  Rank-one classes are
covered by Lemma U / Prop E with the finite (H1) and the b = 1
cushion; the b = 0 discount is a mixed-class phenomenon (witness at
(12,12), all minima at (6,6)) invisible to every cylinder, so no
surgery — straight or twisted — can reach it, and the b = 1 lower
half is not shortened by this route: it is relocated to (M), which
Corollary P1 hands to the three-cover descent (A38's deficit
machinery).  New certified facts: Theorem P at five frames; ker ρ_u^*
= W_u at three; gross diagonal-parity minima 16/≥18; N(u) ≥ 12 for
the r = 1 helix directions.  Open: the (18,12) mixed-sector equality
question; (H1) at r ≥ 2 (engine-bound); Theorem P ∀r.
