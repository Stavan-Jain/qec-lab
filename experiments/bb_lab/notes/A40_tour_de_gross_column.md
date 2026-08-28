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
