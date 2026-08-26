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

## §5 Falsified claims (session)

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

## §6 Residue / next steps

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
