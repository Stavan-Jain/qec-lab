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

*(filled by `scripts/a40_family_lattice.py`; see data/a40/)*

## §3 P2 — pricing (and closure attempt) of the first unproven member

*(filled after §2)*

## §4 P3 — the ∀r analytic lane: feasibility verdict

*(time-boxed assessment; no proof attempts)*

## §5 Falsified claims (session)

*(kept per charter; filled at close)*

## §6 Residue / next steps

*(filled at close)*
