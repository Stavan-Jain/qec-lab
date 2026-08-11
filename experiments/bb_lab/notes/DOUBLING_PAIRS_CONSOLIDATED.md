# Consolidated report — every doubling pair in the program, and how each floor was established

**Compiled 2026-08-11**, revised same day against `main` @ `c0cbb7b` (A33 +
the A34 mitten merge; supersedes the first pass, which was cut at `525d4d7`).
This is an *index and comparison* document, not a new approach line — no
A-number claimed. Sources are cited per row; where a note and its data
disagree, the data file wins and the discrepancy is flagged.

**Revision delta:** A33 adds the `[[144,8,10]] → [[288,8,20]]` pair (Tier C-ii),
an eighth safe-floor species (the *lift-aware per-element* seam floor, §2.8 —
the first pair whose naive coset floor is **false** and which doubles anyway),
a sixth dangerous-floor method (§3.7), the (R)-vs-non-(R) portability scorecard
(§4), the first measured instance *attaining* the deficit wall (§5), and an
erratum against A32's census counts (§3.6). It also surfaced one cross-note
priority overclaim, **reconciled in the same commit** — Tier C and A33 §8.

Companions: `A_HANDOFF.md` (program handoff), `README.md` (approach registry),
`docs/teaching/bb-doubling-theorem.tex` (the template, exposition of record),
`docs/gross-distance-proof.md` (the gross analytic proof).

---

## 0. What a "doubling pair" is here

A base BB code over `G = Z_ℓ × Z_m` with polynomials `(A, B)`, and its free-ℤ₂
cover over `Z_{2ℓ} × Z_m` (x-axis) or `Z_ℓ × Z_{2m}` (y-axis) with the **same**
polynomials. The pair *doubles* when `d(cover) = 2·d(base)`.

The doubling template (teaching doc Thm `thm:template`) reduces that to four
base-side conditions:

| | condition | how it is discharged |
|---|---|---|
| **C1** | cover uses the same polynomials | automatic in this setup |
| **C2** | **(R)**: `1 + x^ℓ ∈ (A,B)` in the cover ring | A12 theorem: (R) ⟺ `k(cover) = k(base)` ⟺ Bezout. Rank check or one polynomial identity |
| **C3a** | **dangerous floor**: `\|b\| + 2m(b) ≥ 2d` for every base stabilizer `b` | §3 below — the varied half |
| **C3b** | **safe floor**: every base cycle with class in `im Δ ∖ 0` has weight `≥ 2d` | §2 below — the cost centre |
| **C4** | **witness**: a min-weight base logical outside `im Δ` | free once C3b holds at `2d > d` |

Two facts govern everything downstream:

- **Parity (A17-P3 L1 / A32 Lemma 2).** If `\|A\|, \|B\|` are odd, every cycle
  weight is even. So a floor at `2d` is decided by a query at `2d − 2`, never
  `2d − 1`, and `2d − 2` is the unique maximal *failing* value (the deficit
  wall). This halves every stratum table in the program. It is **provably
  scope-critical**, not decoration — see the even-weight counterexample in §5.
- **Prop A14.1 (A14 §2).** Under (R): `p₂ = 0`, `Δ = δ₂` is injective,
  `im p₁ = im δ₂`, safe classes are exactly `Δ(ker ∂₂ ∖ 0)`, `δ₂[ζ] = [seamC ζ]`
  by a closed x-carry formula, and **coset minima are constant on G-translation
  orbits** — so one representative per orbit suffices. Every safe-floor method
  below consumes this.

---

## 1. The pairs

Tiers are by *strength of the distance claim*, strongest first.

### Tier A — kernel-checked in Lean, unconditional

| pair | axis | groups | A ; B | floors |
|---|---|---|---|---|
| `[[72,12,6]] → [[144,12,12]]` **gross** | x | Z₆×Z₆ → Z₁₂×Z₆ | `y + y² + x³` ; `y³ + x + x²` | both analytic (§2.1, §3.1) |
| `[[36,4,4]] → [[72,4,8]]` **pair72** | x | Z₃×Z₆ → Z₆×Z₆ | `x² + y + y³` ; `1 + x + y²` | both by finite sweeps (§2.2, §3.2) |

Lean capstones on QECLean `main`:
`grossStabilizerCode_hasCodeDistance_12_uncond` (Gross/Distance.lean),
`pair72_pauli_distance_eq_8` (Z3Z6/Distance.lean) — both hypothesis-free.
Axiom audit: standard three + `native_decide`.

### Tier B — Lean, all three floor Props discharged on an open branch

| pair | axis | groups | A ; B | status |
|---|---|---|---|---|
| `[[150,8,8]] → [[300,8,16]]` **f2a6f17e** | y | Z₅×Z₁₅ → Z₅×Z₃₀ | `1 + y + x` ; `xy⁶ + xy¹⁰ + x²y¹²` | see below |

`cover300_pauli_distance_eq_16` on QECLean `main` still carries its three
hypotheses (`LogicalFloor 8`, `DangerousFloorNZ 16`, `SeamCosetFloor 16`).
All three were discharged as sorry-free Lean theorems (`logicalFloor_8` A21,
`lightClassification` A22, `seamCosetFloor_16` A23, plus
`dangerousFloorNZ_of_lightClassification` A17 §6.5) — **on PR #61, which is
still OPEN and not merged**. The teaching doc's "Lean, unconditional" row
describes the PR branch, not `main`. This is the program's first solver-free
`d = 16` and the deepest single instance.

### Tier C — certificate tier, end-to-end (both floors certified, no Lean)

Four `d = 10 → 20` doublings across two independent lines.

**C-i — literal lifts, A30 (2026-08-07).** The deficit wall cleared.

| pair | axis | groups | A ; B | wall |
|---|---|---|---|---|
| `[[180,4,10]] → [[360,4,20]]` **37a70e02** | x | Z₁₅×Z₆ → Z₃₀×Z₆ | `1 + y + x` ; `y⁴ + x + x¹¹y²` | 17.4 min |
| `[[180,4,10]] → [[360,4,20]]` **5e50a976** | x | Z₁₅×Z₆ → Z₃₀×Z₆ | `1 + y + x` ; `y⁴ + x⁸y² + x¹³` | 23.1 min |
| same code | y | Z₁₅×Z₆ → Z₁₅×Z₁₂ | (same) | (same pass) |

Full scoreboard per cell: (R) ✓, `LogicalFloor 10` both CSS sides ✓,
`d_base = 10` exact (refutation@11) ✓, `LightClassification@18` ✓ (2,203 resp.
2,371 classes), `SeamCosetFloor 20` ✓, `(M)@20` ✓ (all classes PASS),
lifted weight-20 witness ✓. X-side; Z by the BB transpose duality.

The same machine, driven from the **cover** as sole input, re-derives gross
(1.1 s) and pair72 (0.1 s) at this tier — `data/certify_runs/`.

**C-ii — twisted lift on an all-(R) tower, A33 (2026-08-11).**

| pair | axis | groups | A ; B | critical path |
|---|---|---|---|---|
| `[[144,8,10]] → [[288,8,20]]` **IBM class Y** | y | Z₁₈×Z₄ → Z₁₈×Z₈ | base `1 + x + x¹⁴y` ; `1 + xy² + x²y³`<br>cover `1 + xy⁴ + x¹⁴y` ; `1 + xy² + x²y⁷` | **~105 s** |

Rung 2 of the tower Y2 `[[72,8,6]]` (Z₁₈×Z₂) →y→ Y4 `[[144,8,10]]` →y→ Y8.
**Both rungs are twisted lifts** — the cover's polynomials are not the base's,
so template C1 does not apply literally — but **(R) holds on both y-rungs**
with banked Bezout witnesses (`1+y⁴ ∈ (A₈,B₈)`, `1+y² ∈ (A₄,B₄)`), and `k = 8`
throughout. Hypothesis ledger H1–H6 all discharged at certificate tier; the
one that was open, H5, is the *lift-aware* seam floor (§2.8). IBM's MILP-exact
`d = 20` is reproduced **with no SAT on the critical path**, against ~16.5 h of
banked SAT partials — a ~600× cut on that obligation.

> **Citation, disambiguated (2026-08-11).** The source is arXiv:2606.02418
> **supplemental** Table II — *"All verified codes at n = 288 (98 polynomial
> representations, 49 distinct codes), sorted by FOM = kd²/n descending"*,
> caption *"d and FOM use MILP exact distances"* — class Y row:
> `Y (18,8) | 1+xy⁴+x¹⁴y | 1+xy²+x²y⁷ | k = 8 | d = 20 | FOM = 11.1`.
> The supplement numbers its own tables from I (I–III = the CSS catalog by
> block length, n = 144/288/360; IV–X = the non-CSS PBB catalog, whose codes
> take **four** polynomials `(A,B,C,D)`) and refers to the main text as
> "main text Table N". **Main-text Table 2 does not list this code** — it is a
> much shorter table in the four-polynomial format. Every note in the repo
> carried the bare "Table II", which reads as the main-text table; corrected
> in place at A20 §1 (with a full citation note), A20 §7, and A33 §0.
> Independently re-derived here: `(18,8)` with those polynomials gives
> `n = 288`, `k = 8`, and `FOM = kd²/n = 8·400/288 = 11.1` — matching the row.

> **Priority — reconciled 2026-08-11 (this pass).** Three sites called this
> "the program's first certified `d = 20`" outright: the A20 §7 addendum, the
> `A_HANDOFF.md` A33 block, and the `research_log.md` entry. A30 §5.5, dated
> four days earlier, had already taken both `[[360,4,20]]` codes to end-to-end
> certificate tier — same claim tier. All three now carry the qualified form:
> **first certified `d = 20` for a previously published code** (an independent,
> solver-free reproduction of IBM's MILP value) and **first via the tower
> calculus**. Source of truth: the priority paragraph at the head of the A33
> note; erratum logged in A33 §8. No computational result was affected.

### Tier D — safe floor certified, dangerous floor **not run**

Eighteen `(code, axis)` cells over fifteen distinct `d = 8` bases; targets
`[[252,8,16]]` and `[[300,8,16]]`. All SF-CERTIFIED at 16 (A17 §6.1).

**Z₂₁×Z₃ `[[126,8,8]] → [[252,8,16]]`, all with `A = y + y² + x³`:**

| id | axes certified | B |
|---|---|---|
| `68d8f03d` | x | `x² + x⁹ + x¹³` |
| `16884e06` | x, y | `x⁴y² + x¹¹ + x¹⁵y` |
| `3935e2cb` | x | `x⁷ + x¹⁴ + x¹⁸` |
| `71220ff9` | x | `x⁴ + x¹¹ + x¹⁵` |
| `63603d9b` | x | `x⁵ + x¹² + x¹⁶` |
| `e21c6389` | x, y | `xy² + x⁸ + x¹²y` |
| `8494e61f` | x | `x⁶ + x¹³ + x¹⁷` |
| `af479cf3` | x | `1 + x¹⁰ + x¹⁷` |
| `b1e272d5` | x | `x + x⁸ + x¹²` |
| `e9224b61` | x | `x³ + x¹⁰ + x¹⁴` |
| `873f8daa` | x | `1 + x⁷ + x¹¹` |
| `f35534f0` | x | `x³ + x¹³ + x²⁰` |

**Z₅×Z₁₅ `[[150,8,8]] → [[300,8,16]]`, all with `A = 1 + y + x`:**

| id | axes certified | B | extra |
|---|---|---|---|
| `f2a6f17e` | y | `xy⁶ + xy¹⁰ + x²y¹²` | promoted to Tier B |
| `38d3c884` | x, y | `y⁹ + y¹² + x²y⁴` | `d_safe = 16` **exact** on x; `≥ 18` on y |
| `ac46bbea` | y | `1 + y¹¹ + xy²` | `d_safe = 16` **exact** |

What is owed on these seventeen non-f2a6 cells: the `(M)`/dangerous half and
the cover-side value. The A30 rung engine that closed `(M)@20` for Tier C was
never pointed at them — `data/a30/` contains no `16884e06` rung file. A30 §5's
"the `[[126,8,8]] → [[252,8,16]]` bonus cell completes its template the same
way" states the route, not a run.

### Tier E — solver-exact doubling, floors partial or refuted

| pair | axis | groups | A ; B | floors |
|---|---|---|---|---|
| `[[168,12,6]] → [[336,12,12]]` (A8) | x **and** y | Z₆×Z₁₄ → Z₁₂×Z₁₄ / Z₆×Z₂₈ | `1 + y + x³y³` ; `1 + x + x²y⁷` | safe: certified 12 (A29, x-axis, 10.2 s). dangerous: `≥ 12` binding only at `b = 0`, transferred analytically from gross; one rung open |
| `[[72,12,6]] → [[144,12,12]]` **hit3** | y | Z₆×Z₆ → Z₆×Z₁₂ | `y³ + x + x²` ; `y + xy² + x²` | safe: SF-CERTIFIED 12 (A14 S4, 26 s). dangerous: not run |
| **hit4** | y | Z₆×Z₆ → Z₆×Z₁₂ | `y³ + x + x²` ; `y² + xy³ + x²y` | safe: SF-CERTIFIED 12 (32 s) |
| **hit6** | y | Z₆×Z₆ → Z₆×Z₁₂ | `y³ + x + x²` ; `xy + x²y² + x³` | safe: SF-CERTIFIED 12 (21 s) |
| **hit2** | x (anchorable pres.) | Z₆×Z₆ → Z₁₂×Z₆ | `y³ + x + x²` ; `1 + xy⁵ + x²y` | d(cover) = 12 exact; rescued by presentation / descent twist `εB=001` |
| **hit5** | x (anchorable pres.) | Z₆×Z₆ → Z₁₂×Z₆ | `y³ + x + x²` ; `y⁵ + xy + x²` | d(cover) = 12 exact; rescued by the mixed extension class at zero twist |

All six Z₆×Z₆ anchorable classes are gross twins. **Caveat (A11 Entry 1):**
literal-lift doubling is *presentation-sensitive within the same (code, axis)*
— the stored-form ladders that read hit2/hit5 as non-doubling were artifacts.
Descent-space verdicts *are* code-level (A10 Lemma L1).

**The 152-pair T1 sweep (A9).** Direct-sweep frames, `d = 4` bases, all with
`d(cover) = 8` SAT-exact:

| base | axis | pair | count |
|---|---|---|---|
| Z₃×Z₄ | y | `[[24,4,4]] → [[48,4,8]]` | 10 |
| Z₃×Z₅ | x / y | `[[30,4,4]] → [[60,4,8]]` | 25 / 8 |
| Z₃×Z₆ | x / y | `[[36,4,4]] → [[72,4,8]]` | 20 / 31 |
| Z₄×Z₆ | y | `[[48,4,4]] → [[96,4,8]]` | 58 |

Polynomials per row: the table in `A9_lean_target_screen.md`. Every row is
(R) ✓ + linchpin ✓ + tight witness ✓. **The split that matters: 111 rows have
safe-class coset minima `(8,8,8) ≥ 2d`; 41 rows have minima `(6,6,6) < 2d` and
double anyway.** Those 41 are the *overlap rescues* — the sheet-overlap term
`2\|v₀ ∧ v₁\|` carries them. They are why C3b is **sufficient, not necessary**,
and they are the only pairs in the program whose distance rests on nothing but
a SAT ladder.

### Tier F — (R)-free tower doubling

| pair | groups | A ; B | claim |
|---|---|---|---|
| `[[180,8,12]] → [[360,12,24]]` Bravyi (A32) | Z₃₀×Z₃ → Z₃₀×Z₆ | cover: `x⁹ + y + y²` ; `y³ + x²⁵ + x²⁶` | `d = 24` = 2·d(BY), certificate tier |

The lift is **twisted**, not literal: `B` descends by exponent reduction, so
BY carries `B = 1 + x²⁵ + x²⁶` (and GB `1 + x¹⁰ + x¹¹`, BX `y³ + x¹⁰ + x¹¹`)
against the cover's `y³ + x²⁵ + x²⁶`. The whole calculus is twist-generic and
never sees it.

Not a template instance: (R) fails on **all three decks** and `k` jumps 8 → 12,
so C2 is violated outright. It doubles anyway. Resolved Bravyi et al. Table 3's
"≤ 24" to "= 24", exceeding the published solver-exact BB record
(`[[288,12,18]]`) by 25 % in n and 33 % in d. The tower is
GB `[[90,8,8]]` (Z₁₅×Z₃) →x→ BY `[[180,8,12]]` (Z₃₀×Z₃) →y→ C. Machinery in
§3.6. First *perfect* doubling rung with a deck-nontrivial top — a genuinely
new data point for the deficit-wall taxonomy, which was built entirely on (R)
rungs.

**Trust basis, leg by leg** (the question "certificate or solver?" answered by
tracing the assembly rather than quoting A32's summary line):

| assembly leg | what carries it | trust |
|---|---|---|
| `[b] ≠ 0` — sectors A/B/C | GB coset-BZ censuses (exact binomial node counts, kernel-asserted) + MITM fiber enumerations complete by the exact-subset-sum argument, every solution re-verified + 49,855 rungs all PASS | certificate |
| `[b] = 0, b ≠ 0`, `\|b\| ≤ 20` | banked `(M)@24` floors = **8,461 solver-free rung passes**, `non_pass: []` (`data/a19_scope/scope_results.json`) | certificate |
| census completeness ≤ 20 | fiber-sweep membership (8,313 checks) + the analytic `β = 0` family — *replaces* the lex-leader SAT engine's terminal UNSATs | certificate |
| `\|b\| = 22` flat residue | parity (Thm 5 kills `β = 0` at **zero compute**) + 27,152 flat-top rungs | certificate |
| `b = 0` ⟹ `\|v\| = 2\|u\| ≥ 2·d(BY)` | `d(BY) = 12`: 45 fibers, 7 lifts all stabilizers, 0 logicals, 0.04 s; weight-12 witness re-verified | certificate |
| upper side `d ≤ 24` | 100 weight-24 τ-lift vectors re-verified as nontrivial logicals | certificate |

So **Tier F is certificate trust end to end**, at the same level as Tier C and
strictly above Tier D (§2.4: CMS UNSATs are proof-less by construction — DRAT
is disabled under its Gauss-Jordan XOR reasoning). SAT survives only as
historical cross-check (A19's original UNSATs, A24's engine censuses), and A32
§5 explicitly retired the last SAT-tier critical-path input (`d(BY) = 12`)
during that session. Certificate tier is still **not** kernel tier: nothing here
is Lean-checked, same as Tier C.

Two clarifications this table is designed to pre-empt:

- **ISD is not a trust dependency.** The weight-24 witnesses were *found* by a
  heuristic search, but a witness only ever needs checking, and the check is
  deterministic F₂ linear algebra. Search method is irrelevant to the trust tier
  of an upper bound — the same asymmetry (cheap witnesses, expensive UNSAT) that
  runs through the whole program, and the reason it never contaminates a floor.
- **The genuine caveat is reproducibility, not trust.** A32 reads `data/a19/`
  and `data/a24/` from the main checkout, and `experiments/bb_lab/.gitignore`
  ignores `data/` wholesale — subdirectories are force-added case by case.
  `data/a32/`, `data/a33/`, `data/a30/` and `data/a19_scope/` are committed;
  **`data/a19/` and `data/a24/` are not.** The 8,461-rung *result* is therefore
  in the repo but the per-class detail behind it is local-only, so an
  independent re-runner must regenerate the A19 census first. A32 §5's
  census-completeness bonus softens this — the banked census is re-derived by
  the fiber union, demoting it from input to cross-check — but the floors leg
  still points at ungitted data. Worth force-adding a per-class digest.

---

## 2. How the **safe** floor was established — eight distinct methods

Every method certifies the same statement (`SeamCosetFloor 2d`: no element of a
safe-class coset `t₀ + C_AB` weighs `< 2d`), one G-orbit representative per
class by Prop A14.1(4). They differ in *what makes the search finite*.

### 2.1 Analytic slot-frame / confined floor — **gross only**

`docs/gross-distance-proof.md` §§10–13. The CRT frame splits
`F₂[Z₆×Z₆] = F₂[Z₂²] ⊗ F₂[Z₃×Z₃]` into `F₂[Z₂²] × F₄[Z₂²]⁴`; the ρ-links
`ρᵢ = B̂ᵢÂᵢ⁻¹` are radical, `ρ² = 0`, and the confined frame gets a *slot* cost
calculus from the layer weight dictionary `d₃` (flat 6 over Z₃×Z₃), closed by a
118-achiever ρ-link kill. Lean: the MIm engine — Γ-membership + 13 y-orbit
transports + a 64-case dispatch (`MImClassify/MImFloor*/MImMembership/
MImTransport/MImAssembly`), `LightStab.mimBound_holds`.

**Why only gross.** The dictionary is F₄-specific: co-point rigidity ("one zero
and three pairwise-distinct nonzero values") uses the fact that F₄ has exactly
three nonzero elements. Over Z₃×Z₇ the split is `F₂ × F₄ × F₈² × F₆₄²` and the
rigidity is simply false — A8 §4.3 declared re-deriving it with a heterogeneous
8/12/14 dictionary "genuinely undeveloped mathematics" and made it the note's
open core. **A29 §5.2 closed that core without it** (see 2.5): the Tier-3 wall
was never on the path, because the fibering needs no value rigidity at all.

### 2.2 Direct finite kernel sweep

Enumerate the safe-class coset over the base cells: `2^(leaf bits)` per class,
12–24 bits on the A9 frames, `2^18` for pair72. Kernel-grade in Lean by
`native_decide` / packed-Nat tables. Used for pair72 and all 152 T1 rows.
Dies immediately past `n ≈ 50` base cells.

### 2.3 S4 — budgeted per-orbit coset SAT with membership-row augmentation

A14 §11. Augment `H_X` with `D @ L_X` membership rows so a SAT model *is* a
safe-coset element; f-vars + XOR-chained output cells + `seqcounter`
cardinality; CaDiCaL; every SAT answer re-verified in numpy. Query per orbit
rep: "∃ coset element of weight `≤ floor − 1`". Certification = all-UNSAT.

Cost is wildly asymmetric — SAT witnesses are cheap, UNSAT is the expensive
side. Certifies gross-x in 24 s (independently reproducing the whole Lean MIm
engine at solver grade) and hit3/4/6-y in 21–32 s; refutes bb288-y in 60 s with
a weight-34 witness.

### 2.4 XOR-native CMS at `floor − 2`, plus the parity lemma

A17 §6.1. Same query shape as S4, but (i) parity makes the `floor − 1` query
vacuous so the decisive bound is `floor − 2`, and (ii) the XOR rows are handed
to **cryptominisat5 as native DIMACS `x`-lines** with no Tseitin encoding.
This is what cracked the docket: queries that ground for hours under
pysat-CaDiCaL at 10M conflicts decided in 3–20 min each. Fourteen cells on the
first pass, four more on a 7200 s retry — the hardest landed query
(`38d3c884:x`) took 7,489 s.

**Trust caveat, stated in the source and repeated here:** CMS UNSATs ride its
Gauss-Jordan XOR reasoning, which is sound but *proof-less* — DRAT is disabled
under Gauss. One cell (`f2a6f17e:y`) was re-proved by kissat on the Tseitin CNF
in 9,506 s, emitting a 6.85 GB DRAT proof; drat-trim/LRAT checking was never
run (build declined + disk). So Tier D's floors sit at solver-correctness
trust, the same tier as the corpus's CaDiCaL `d_exact` values.

### 2.5 Fibering / state-sweep certificates (`bb_lab.fibering`, A29)

The portable engine. Pick `z ∈ G` of odd order `q`, sites = `G/⟨z⟩`:

- **L1 (residue coordinates)** gives an exact closed-form weight
  `\|w\| = Σ W(rₛ, εₛ)`, `W(r,e) = \|r\|` if `\|r\| ≡ e (2)` else `q − \|r\|` —
  for *every* odd `q`, killing A27 §3.2's 2-primitive-prime restriction (the
  field splitting was scaffolding; the invariant content is the residue map).
- **L2**: the fiber augmentation `ε` is a ring map; if `B_ε = x̄^m·A_ε` (the
  **link**) the fiber parities pair up and the twist `τ` is f-independent — an
  invariant of the safe class.
- **L3/L4**: per-site cost bounds by state (OFF/ON, or U0/V0/ON when `τ = 1`),
  then a complete state sweep at budget `2d − 2` with an affine solve per
  assignment modulo the cost-invariant kernel `K₀`.
- **L5**: exact ε-minimisation is a linear functional over an affine space.

Decides `ac46bbea:y` in 1.4 s (link on the standard `y³` fiber, `τ ≡ 0` — a
straight f2a6-shape port) and `38d3c884` in ~20 s/axis (**no link on the
standard fiber**; the engine's enumeration surfaces two diagonal linked fibers,
`z = (1,3)` with τ-weight 8 and `z = (1,12)` with `τ ≡ 0` — hand-picking fibers
would have missed both). And it certifies the A8 base floor at 12 in 10.2 s on
the `q = 7` unpaired frame with relaxed δ-costs alone (minima 12/15/15, zero
ε-suspects), closing A8's open core.

**Where it refuses** — with named caps, never a false claim:
Z₂₁×Z₃ (no ε-link on any fiber, ε-image ~2²¹ over the enumeration cap),
Z₁₅×Z₆ floor-20 (no link; unpaired DFS ≈ C(36,≤18) ≈ 3·10¹⁰ nodes),
gross (S = 12 vs budget 22 — vacuous on its only fiber, consistent with its
slot-frame history). Negative control passed: it does **not** certify by90's
known-false floor 16.

### 2.6 Coset-BZ two-window enumeration (A30)

Don't fiber — **enumerate**, with a counting invariant as the certificate. Take
disjoint information sets `I₁, I₂` for `C_AB`; every coset element is the
unique one matching its own window restriction, so the A28 census walk applies
verbatim, merely **seeded at `c_∅ⱼ` instead of 0**. Completeness is the
asymmetric pair bound: `r₁ + r₂ + 2 > W` ⟹ any `\|c\| ≤ W` is light in some
window. At `W = 18` the schedule `(r₁,r₂) = (9,8)` runs 9× cheaper than `(9,9)`.

Certificate = windows disjoint (listed) + systematic identity blocks + **exact
node counts** `Σ_{s≤r} C(κ,s)` asserted in-run + coset parity + ∅-pattern base
words re-checked in Python. **No UNSAT anywhere** — a refutation would return a
verified witness, and at a complete pair the found minimum is exact.

This is what decided the three Z₁₅×Z₆ floor-20 cells (κ = 88, ~7·10¹¹ nodes,
0 hits, 6.3–8.3 min/code). Validation: 10/10 against every recorded coset
ground truth in the program, including both A29-new decisions — **an
independent certificate species**, no shared trust base beyond the seam-offset
construction. Honest engineering note: the first production attempt used the
naive `(9,9)` schedule and hit the 15-minute budget on both codes; the 1.8× cut
came from arithmetic alone.

### 2.7 Analytic seam-coset floor via CRT fibering (A23, Lean)

f2a6's `seamCosetFloor_16`, reduced to the single inequality
`\|A⋆f + e₀\| + \|B⋆f\| ≥ 16` and proven by an A22-fibering site sweep
(`z = y³`, `w = y⁵`). Covering proofs were *refuted* first; the site sweep is
what worked. Sorry-free.

### 2.8 The lift-aware per-element seam floor (A33 H5) — when the naive floor is FALSE

Every method above certifies `SeamCosetFloor 2d`: *no* safe-coset element weighs
`< 2d`. On IBM class Y that statement is **false** — the 12-orbit seam coset
contains a weight-18 element, and A20 §5 recorded the naive floor as refuted.
The pair doubles anyway, because the floor was the wrong obligation: what
matters is not the weight of the base-side coset element `b` but the weight of
the cover chains lying over it. A33 replaces the coset-minimum statement with a
**per-element** one:

1. **Census the seam coset completely** to `2d − 2` (parity), by a 3-offset
   coset-BZ pass — stabilizer offset plus one offset per seam orbit, all in one
   kernel run (κ = 68, `W = 18`, r-pair (9,8), 6.62·10¹⁰ nodes, exact
   node-count asserts, 61 s). Result: 12-orbit coset = **1,680 elements**
   `{14:6, 16:84, 18:1,590}`; **3-orbit coset EMPTY** — so the naive floor is
   actually *true* there, and the A20 refutation is orbit-local, not code-level.
2. **Run a seam rung on every censused element.** These are
   **feasibility-only** rungs: by stabilizer transport every cover cycle over a
   nonzero-class element is automatically a nontrivial logical, so there is no
   class dispatch and no triviality window inside the rung (A32's "pure
   feasibility" observation, held verbatim). **1,680/1,680 PASS in 1.8 s**,
   94.6 % of them flat-top at `M = 1`.
3. **G-transport** covers the other 14 classes: the fold `G(Y8) ↠ G(Y4)` is
   onto, so every class translates to an orbit rep.

Cross-derived a second way (§3, "by descent"): 105,328 Y2-shadow fibers, 90.6 %
carry-infeasible, yielding **exactly the same 1,680 elements**. Banked
agreement: all 278 SAT-census elements present (the dying SAT census was two
classes short of the true 280), plus the SAT@18 witness.

**Why this matters beyond the instance.** It separates two things the template
conflates. `SeamCosetFloor 2d` is *sufficient*; it is not the real obligation.
The real obligation is that no cover chain over a safe-class element is both
light and a logical — and when the coset floor fails by a little (here, one
weight-18 element against a target of 20), the per-element form still closes.
This is the same gap the 41 overlap-rescued A9 rows exhibit (§1 Tier E) — but
here it is *certified* rather than merely observed by SAT.

### 2.9 The comparison

| method | scaling knob | certificate species | reached |
|---|---|---|---|
| slot frame (2.1) | needs F₄-uniform CRT split | hand proof + Lean dispatch | gross |
| direct sweep (2.2) | `2^(base cells)` | `native_decide` | n ≤ ~50 |
| S4 SAT (2.3) | solver, UNSAT-side | solver trust | n ≤ ~150 |
| CMS XOR (2.4) | solver + parity halving | solver trust (Gauss ⟹ no DRAT) | n ≤ 180, hours |
| fibering (2.5) | link present? τ-weight? S vs `2d−2` | state-sweep + rank certs | seconds — when it fits |
| coset-BZ (2.6) | `C(κ, W/2)` node count | exact counting invariant | κ ≈ 88, `W = 18`, minutes |
| A23 fibering (2.7) | site count | Lean, sorry-free | f2a6 |
| lift-aware per-element (2.8) | coset size × per-element rung | BZ census + feasibility rungs | the only method that survives a false coset floor |

**The one structural bit that decides which applies:** whether some odd-order
cyclic fiber carries the ε-link `B_ε = x̄^m·A_ε`. With it, 2.5 runs in seconds.
Without it (Z₂₁×Z₃, Z₁₅×Z₆), you fall back to 2.4 (solver-hours, no proof
object) or 2.6 (minutes, with a counting certificate) — and 2.6 is strictly
better on both cost and auditability, which is why A30 superseded the SAT lane.

---

## 3. How the **dangerous** floor was established — seven distinct methods

Statement: for every base stabilizer `b`, `\|b\| + 2m(b) ≥ 2d`, where `m(b)` is
the sheet-overlap of a cover chain over `b` (slice identity
`\|v\| = \|b\| + 2m(b)`). It splits into `b = 0` (the zero rung — a base logical
floor, easy) and `b ≠ 0`, which needs a **classification of all light
stabilizers** (`0 < \|b\| ≤ 2d − 2`) plus a per-class rung bound.

The census depth is the whole story:

| instance | `2d` | census depth | classes | note |
|---|---|---|---|---|
| pair72 | 8 | ≤ 6 | 24 | base chosen for **zero** seam-hostile classes |
| gross | 12 | ≤ 11 (parity: ≤ 10) | hexagons + D-pairs (Prop 10); **7** translation classes `{6:1, 10:6}` by the A30 census | flat, hand-surveyable |
| A8 Z₆×Z₁₄ | 12 | ≤ 11 | 84 hex + 504 D-pair + **21 weight-8 (new)** | the new class is *global* |
| f2a6 | 16 | ≤ 14 | **113** (94 small-preimage + 19 near-kernel) | two-strata |
| Z₁₅×Z₆ ×3 | 20 | ≤ 18 | 2,203 / 2,371 / 2,371 | engine-generated |
| IBM Y8 | 20 | ≤ 18 | 1,655 | **four independent derivations** |
| Bravyi C | 24 | ≤ 22 | ~20k projected — **stalled at 151** | forced the tower |

### 3.1 gross — light-stabilizer classification + m-rungs (analytic, Lean)

Classification to weight 11 (hexagons `\|A\|+\|B\| = 6`, D-pairs), then coset-
counting rungs `m(hexagon) ≥ 3`, `m(D-pair) ≥ 1`, `m(0) ≥ 6`. In Lean both
named Props are kernel-discharged (`lightStabilizerClassification_holds`,
`mimBound_holds`). Note the `b ≠ 0` dangerous minimum is actually **14**, not
12 — the hexagon/D-pair rungs are not binding.

### 3.2 pair72 — finite dispatch on a hand-picked base

`pair72_dangerousFloorNZ`, unconditional in Lean. The base was *selected*
because all 24 light classes have seam-good translates, so the single-shape
rung closes everything with no window machinery. This is the reason pair72 is
cheap and the reason it does not generalise: A9's stage-3 gate found every
`q₁ = 0` candidate seam-hostile (6–18 of 24 classes failing), which is exactly
why the T1 target was switched to the doc-verified pair.

### 3.3 f2a6 / cover300 — the rung dispatch (the deepest instance)

`dangerousFloorNZ_of_lightClassification : LogicalFloor 8 → LightClassification
→ DangerousFloorNZ 16`. Three sessions of structure:

1. **§6.3**: 94 light classes with a support-≤4 preimage; dispatch over
   94 × 75 = 7,050 cells into single-shape (5,765), pair-shape (1,170) and a
   new window rung (115).
2. **The cutoff conjecture was FALSE.** A completeness certificate came back
   SAT in 0 s: there are `\|b\| = 10` boundaries whose minimal preimage weighs
   31–33 — *balanced simultaneous near-annihilation*, distance ≥ 31 from both
   kernels, not kernel-plus-perturbation. **No support-bounded census can
   complete the classification**; the light-boundary set of a `d ≥ 7` base is
   genuinely two-strata.
3. **§6.4 classified the stratum.** The pair is secretly σ-correlated:
   `B = xy⁶·σ(A)` with `σ:(x,y) ↦ (xy⁶, y⁴)` of order 2. Exhaustive enumeration
   (9.6 h SAT with translation-orbit blocking; later re-derived by the A29
   engine in 65 s) gives **113 classes** — the 94 reproduced exactly, plus
   19 near-kernel. The coset-min spectrum is strictly bimodal
   `{1,2,3} ∪ {31,32,33}`: the `[5,30]` rep gap is **empty**. Unified dispatch:
   8,475/8,475 light cells covered, 0 uncovered, with exact (not probe-grade)
   window conditions at `t = 1, 2, 3`.

**Build-cost arc** (worth remembering): ∀-Fin ball `native_decide` dies
(`decidableBallLT` overflows the C stack past ~2²³, ~9+ CPU-h below it) →
falsifier-filter cores (53 min) → **Gaussian pivot certificates**
(`KernelCert.lean`: sweeps are *rank facts*; all 37 window systems admit
no-row-op elimination orders) → the full 560M-mask obligation set builds in
**3.9 s**.

### 3.4 A30 rung pass — `(M)@20` as an affine carry system

The Tier-C method, and the one that generalises. Per class: the cover-cycle
condition is the affine system `E·v₀ = rhs(b)` with `E` **class-independent**
(one 180×180 matrix per cell); solutions decompose into `≤ 2^{k/2}` sectors of
`ker E / im S`; violations are hunted by **two complete lanes** —
restricted-support enumeration with meet-in-the-middle subset sums over `E`'s
reduced columns (the heavy strata, `M − 1 ≤ 4`), and multi-offset coset-BZ over
the same κ = 88 windows (light strata). Every candidate is re-verified
(solution check, slice identity, non-boundary) before it could count.

Results: 2,203 / 2,371 / 2,371 classes **ALL PASS** in ~7 s each. Validation:
f2a6:y 113/113 PASS in 1.7 s — **reproducing the Lean theorem's instance**,
with the lane split (1+7+36+69) matching the census histogram exactly — plus a
soundness control that hand-builds a genuine dangerous logical (weight 20,
overflow 7) and confirms the checker FINDS it at an inflated target.

Sector-completeness: every nonzero cover class's shadow is either a stabilizer
(`b = 0` by LogicalFloor; light `b` by census completeness; heavy `b` by the
slice identity) or a base logical whose class lies in `im p₁ = im δ₂`, hence in
a Δ-class whose full coset the safe floor certified — so `wt(v) ≥ wt(p(v)) ≥ 20`.

### 3.5 A8 — analytic transfer, one rung open

The `b = 0`, hexagon and D-pair `(M)` lower bounds transfer from gross verbatim
(they depend on `d(base) = 6` plus support geometry `\|hexagon\| = 6`,
`\|D-pair union\| = 11`, both identical). Achieved minima are SAT: 12 / 16 / 20
/ ≥24, so the sector minimum is **exactly 12, achieved uniquely at `b = 0`**.
The one gap: the new 21-element weight-8 orbit is **global** (minimum
decomposition 36 hexagons, hexagon-union 110 of 168 cells), so gross's *local*
coset rung — which needs `U ≤ 9` — is hopeless. It is far off the critical path
(SAT gives `m ≥ 8`), but a presentation-free `m(weight-8) ≥ 2` remains open.
Reproduction trap recorded there: the constrained cover-SAT needs the BB
inversion duality `Φ` applied to `b` first, or it is structurally UNSAT for
every `b ≠ 0`.

### 3.6 Bravyi C — the tower slice calculus (no (R), no dangerous census at all)

The `2d − 2 = 22` census at `n = 180` was the wall (~3.3·10¹³ coset-BZ nodes,
40–60 h; the stabilizer census stalled at 151 of ~20k classes). A32 replaced it
by descending a second rung:

- **Lemma 1** (twist-generic): fold is a chain map, `τ` transports stabilizers,
  `p∘τ = 0`, `τ∘p = 1+σ`, slice identity, and `v` is a cycle ⟺ shadow `b` is a
  cycle and the carry system `E v₀ = R b` holds.
- **Theorem 3** (two-level slice): `\|v\| = \|b\| + 2m₁ = \|β\| + 2(m₁+m₂)`,
  with the **overflow square** `m₁ + m₂ = m_x + m_y′` — total overflow is
  path-independent.
- **Theorem 4** (GB-sector trisection): every `d = 24` violation candidate falls
  into sector A (`[β] ∈ W ∖ 0`, `\|β\| ≥ 14`), B (`β = 0`, forced
  `(\|γ\|, m₁) ∈ {(8,≤3),(10,≤1)}`) or C (`β` a GB stabilizer).
- **Theorem 5**: the flat-22 residue's `β = 0` branch is **parity-dead** —
  `\|γ\| = 11` is odd. Zero compute.

Then: censuses at `n = 90` (seconds), bounded-overflow lift fibers per censused
shadow, top rungs per lift. **49,855 rungs / 274k fibers, ALL PASS, ~9.4 min**
— ~250× under the enumeration baseline, and ≥ 93 % of fibers turn out
carry-infeasible, which is the structural reason the tower wins. Bonuses:
`d(BY) = 12` re-derived solver-free in 0.04 s, `d(GB) = 8` in < 1 s, and the BY
stabilizer census re-derived as the fiber union — **one machine serves the
distance floor, the census, and the sector floors at every rung**.

Falsified in-session and worth not re-proposing: the calculus is *not* a
compression statement (banked census compressed only 2.2×); the win is that
census *generation* moves to `n = 90`.

**No SAT is load-bearing here.** The `(M)@24` floors at bands ≤ 20 — the one
leg that reads banked A19 data — are consumed as **8,461 solver-free rung
passes** (`data/a19_scope/scope_results.json`, `non_pass: []`, 735 s, lanes
w6…w22), not as the 8,310 SAT UNSATs they replaced; and `d(BY) = 12`, called
out in A32 §5 as "the last SAT-tier input on the `d(C) = 24` critical path",
fell to 45 fibers in 0.04 s. Full leg-by-leg trust table: §1 Tier F.

**Erratum (found by the A33 port's audit, patched at source).** The
`a32_gb_census.census()` helper never emitted the empty-window coset-base
element. Effect on A32 §3: logicals ≤ 10 — 15 weight-10 vectors missed
(1,623 → 1,638), **all 15 inside stored orbits**, so every per-orbit consumer
(sector B, the dby floor, `d(GB) = 8`) is unaffected; W-cosets ≤ 22 —
2 weight-22 elements missed, both orbits enumerable via translates, both cap-0
fibers empty. **The `d = 24` assembly is unaffected**; raw vector counts carry
`+15 / +2` corrections and the censuses were re-run post-fix.

### 3.7 IBM class Y — census-by-four + deterministic rungs

The A33 dangerous half is the A30 architecture at Y4, but with the census
established **four independent ways**, which is the strongest completeness
evidence in the program for any single object:

1. banked SAT census (historical), 2. the V7 analytic engine (historical),
3. **direct 3-offset coset-BZ** — 118,932 vectors ≤ 18 collapsing to exactly
1,655 G-canonical classes, obtained *free* as the stabilizer offset of the same
kernel pass that did H5, and 4. **the tower fiber union** — Y2 stabilizer census
≤ 18 (165,517 vectors) → 4,605 orbit reps → bounded-overflow fibers → 8,302
lifts + the τ₁-family, class-key union `==` the banked 1,655 exactly (3.2 s).

Floors: **1,655/1,655 deterministic rung PASSes at target 20 in 1.9 s**, in
agreement with all 1,655 banked SAT UNSATs — a solver-tier → certificate-tier
upgrade *and* a 9× speedup. Engine soundness is anchored by a planted control:
a hand-built genuine weight-26 dangerous logical (over a weight-6 shadow,
overflow 10) is **FOUND** by the BZ lane at `M = 11`, with minimum overflow
exactly 10.

Anti-pattern recorded here: the Prop-10 weight-8 stabilizer gap (which held at
gross, Y4, BY, GB) **fails at Y2** — 36 weight-8 stabilizers exist. Do not cite
that gap as universal.

### 3.8 The comparison

| instance | dangerous method | what made it possible |
|---|---|---|
| pair72 | finite dispatch, single-shape rung | base picked for zero seam-hostile classes |
| gross | classification ≤ 11 + coset rungs | flat F₄ dictionary; surveyable by hand |
| A8 | rung transfer from gross + SAT minima | identical support geometry; new class non-binding |
| f2a6 | 113-class classification + S/P/W rungs | σ-correlation `B = xy⁶σ(A)`; window rungs at t=1,2,3 |
| Z₁₅×Z₆ ×3 | affine carry system + 2 complete lanes | `E` class-independent; MITM + coset-BZ |
| IBM Y8 | census-by-four + 1,655 deterministic rungs | a second deck **and** (R): fiber union re-derives the census |
| Bravyi C | trisection on the *second* shadow | a second deck existed; parity killed a whole branch |

---

## 4. The pattern: which half is hard, and why it moves

- **At `d ≤ 6` the dangerous half is nearly free** (pair72's 24 classes, gross's
  hand table) **and the safe half is the discriminator.** A11 measured this
  directly: certificates (R0-sq…R2, linchpin, k, tight witness) are satisfied
  uniformly and are *never* the discriminator — **every non-doubling cell breaks
  in the safe sector**, and a base-side coset-min probe (no cover SAT) separates
  the doubling presentations from the non-doubling ones exactly
  (`{≥12: 63}` vs `{6:12, 8:45, ≥12:6}`).
- **At `d ≥ 8` the dangerous half becomes the harder one**, because the light-
  boundary set stops being support-bounded (f2a6's near-kernel stratum) and the
  census depth `2d − 2` grows faster than any hand table. This is why Tier D has
  eighteen safe floors and one dangerous floor.
- **The safe half is presentation-sensitive; the dangerous half is not.**
  Literal-lift safe floors depend on the chosen presentation (A11); descent-space
  and code-level verdicts do not (A10 L1).
- **Everything even-weight breaks.** The parity lemma is what makes `floor − 2`
  decisive, makes the deficit wall land at exactly `2d − 2`, halves stratum
  tables, and kills entire branches for free (A32 Theorem 5). Its scope is
  provably necessary — §5.

### 4.1 What (R) actually buys — measured on two towers side by side

A33 §6 ran the A32 calculus on a structurally opposite tower and measured the
difference. This is the cleanest statement in the program of what condition (R)
is worth *operationally*, as opposed to what it means:

| quantity | A33 tower (R) holds, both rungs | A32 tower (R) fails, every deck |
|---|---|---|
| `rank τ*` | `k/2` (drop) | `= rank p*`, no drop |
| exactness | both ways (`im τ* = ker p*`, `ker τ* = im p*`) | fails |
| `σ*` on H1 | `= id` | `≠ id` |
| `im p*` | the `k/2`-dim SEAM, with the `ker ∂₂`/seamC dictionary | no such structure |
| descent trisection | **collapses to one branch** — `SEAM ∩ ker p₁* = 0` and `p₁*` injective on SEAM pin one target coset per seam class | all three sectors A/B/C live |
| safe-sector rungs | flat-top dominated (94.6 % at `M = 1`) | ran to `M = 5` |

So (R) buys **exactness, `σ* = id`, the seam dictionary, and a one-branch
descent** — two of three sectors die to a single rank computation. What it does
*not* buy: the sector space is `2^k` on **both** towers, so A32's "2⁸" was
k-genericity, not an (R) effect.

Three further portability facts, both towers agreeing:

- **Twisted lifts are invisible to the calculus.** All four rungs across the two
  towers are twisted; Lemma 1, Theorem 3, parity and the rung engine held
  verbatim, with the transports asserted rather than assumed.
- **Iterated Z₂ suffices.** The composite Z₄ deck appears only in a consistency
  assert; nothing consumed it. Likewise the overflow square degenerates on a
  same-axis tower (one descent order) and nothing needed it.
- **One machine, three jobs, every rung.** The middle code's distance floor and
  its stabilizer census both re-derive from the same fibers at a smaller budget
  — confirmed now on two structurally different towers (`d(BY) = 12` / `d₁ = 10`;
  the BY census / H1).

---

## 5. Negative results that pin the boundary

| result | what it kills |
|---|---|
| **A11 Entry 3** — `(A,B) = (x²y+x³+x³y³+x³y⁴, x⁴(1+y+y²+y⁴))` on Z₅×Z₅, `[[50,2,5]] → [[100,2,8]]` | C-safe (tight witness ∧ safe minima ≥ 2d) holds two-sidedly at *enumeration* grade, `(R)` holds, and `d(cover) = 8 < 10`. The failing weight-8 logical pushes to a **weight-6 stabilizer** — the first observed dangerous-sector bind. Both safe minima are 11, **odd** — impossible under parity. So "C-safe ⟹ doubling" is false weight-agnostically, and any (M)-robustness statement **must consume parity** |
| **A10** — 13 small-frame bases (Z₃×Z₃/Z₃×Z₄ `d=4`, Z₃×Z₅ `d=6`) | their entire 256-cover descent space fails; 3,328 witness rows, numpy-re-verified. Code-level hard negatives — no equivalent presentation lifts to a doubling literal cover. The universal-doubling claim is FALSE |
| **A14 §13** — all five same-axis rung-2 re-doubles | every rung-2 safe floor refuted **by the free tier alone**: the raw seam minimum equals `d(rung-1)` exactly. Mechanism: re-doubling x leaves the y-loop as the bottleneck, `min(2L, L) = L`. Towers are flat. Distance > 12 must come from a larger-d base or a mixed lift, not iteration |
| **A14 §§14–16** — bb_90 / bb_108 / bb_288 | all four axes fail; the bb_108 2,808-cell orbit sweep finds **no** SF-passing presentation (ceiling 18 < 20); Z₁₈×Z₃ re-decomposition equally dead ⟹ the obstruction is code-level. bb_288 SF-refuted on both axes with certificates |
| **A17 §6 / A17-P3** — the deficit wall | across 2,497 refuted cells the reached-weight deficit-2 bucket dominates (28 %) — and `2d − 2` is the *parity-maximal* failing value, not a mystery. **Measurement discipline: recorded reject weights are first-found witness weights = upper bounds on `d_safe`, never minima.** The SS15/SS16 "orbit ceiling 18/32–34" readings were retracted on exact ladders |
| **A20 rung 1** — Y2 `[[72,8,6]]` →y→ Y4 `[[144,8,10]]` | the first **measured instance attaining** the wall: `10 = 2·6 − 2` exactly, on a rung where (R) *holds*. One tower carries both extremes — a maximal-deficit rung and a perfect doubling rung — stacked. Whatever sets the deficit acts **per rung, not per code** |
| **A33 §9** — wall refinement | on Y8 the safe minima sit at `2d₁ − 2 = 18` on the live seam orbit while the other orbit is **empty**. So the wall phenomenon is **per-orbit**, not per-code either. Feed to the A17-P3 taxonomy |
| **A29 §5.4** — by90 negative control | the engine must not certify a known-false floor 16, and does not |

---

## 6. What is owed

1. **Tier D's dangerous half.** Seventeen SF-certified cells with no `(M)` run.
   The A30 rung engine is built, validated, and ~7 s per cell at `n = 180`; the
   `[[126,8,8]]` cells are `n = 126`. This is the cheapest large win on the board.
2. **PR #61 merge.** The three f2a6 floor theorems are sorry-free and audited but
   live on an open PR; `main` still carries the conditional `cover300` statement.
   The teaching doc's instance table should be re-read as branch-state until then.
3. **Lean packaging of a Tier-C cell** — the natural first unconditional `d = 20`
   BB theorem. Route known (A15 pattern: `{logicalFloor, lightClassification,
   seamCosetFloor}` instance layers + re-instantiating
   `dangerousFloorNZ_of_lightClassification` at `(n,d) = (180,10)`). The
   dangerous-sector Lean theorem exists today **only** for the `[[300,8,16]]`
   instance — its architecture is portable, its artifact is not.
4. **Proof-grade upgrade for Tier D.** drat-trim / cake_lpr on the one banked
   6.85 GB DRAT proof; or better, re-decide those cells under coset-BZ (2.6),
   which produces an auditable counting certificate instead of solver trust.
   The Z₂₁×Z₃ cross-check already ran in **0.6 s**.
5. **A29's ε-recursion chapter** — the named residue. A δ-config's exact cost is
   relaxed + a weighted coset-min in the ε-quotient pair map; at `δ = 0` this is
   exactly `q ×` (the safe-floor problem of the ε-quotient code). Charging
   partial ε-costs inside the DFS is the designed route to the frames the engine
   currently refuses.
6. ~~**A32 port to `[[288,8,20]]`**~~ — **DONE, A33 (2026-08-11).** Every A32 §7
   portability claim held; see §4.1. The updated next-instance ranking is:
   (i) the A30 `[[360,4,20]]` **rung-2 re-doubles** — three-level, all-(R),
   same-axis towers, exactly A33's shape, with A14 §13's "freeze at the toric
   bottleneck" verdicts now re-examinable *per element* by §2.8's machinery
   (that is the sharpest open question the new method creates: the freeze was
   established against the *naive* coset floor, which is precisely the
   obligation A33 showed can be false without killing the doubling);
   (ii) `[[756,16,≤34]]` y-rung (n = 378 base — a tower census is the only
   plausible route); (iii) gross over bb72 over (3,6), retrospectively, as a
   teaching-doc chapter.
7. **A8's weight-8 rung** — presentation-free `m ≥ 2`, off the critical path.
8. ~~**Reconcile the "first certified `d = 20`" claim**~~ — **DONE in this
   pass**; see Tier C, A33 §8, and the A33 note's priority paragraph.
9. **Force-add a per-class digest of the A19 `(M)@24` census.** Tier F's trust
   tier is fine (§1 Tier F), but its floors leg reads `data/a19/`, which
   `.gitignore` excludes; only the 8,461-rung summary is committed. A digest
   keyed by class would close the reproduction path without shipping the full
   census. Same audit applies to `data/a20/` for Tier C-ii and `data/a24/`.

---

## Appendix: reproduction map

| claim | check |
|---|---|
| Tier A/B/C pair list, polynomials | `data/a17/docket_decision.jsonl`, `data/a30/decide_*.json`, `data/certify_runs/*/verdict.json` |
| Tier C end-to-end, cover as sole input | `bb-lab certify` / `bb_lab.doubling_certify`; `tests/test_doubling_certify.py` (5/5) |
| Tier D safe floors | `data/a17/docket_decision.jsonl` (28 rows incl. the UNKNOWN→CERTIFIED retries) |
| Tier D engine second opinions | `data/a29/docket_engine_final.jsonl`, `data/a29/extra_targets.jsonl` |
| Tier E hits + ladders | `data/a9/t2_presentation_hits.json`, `data/a9/t2_cover_ladders.json` |
| the 152 T1 pairs + the 111/41 split | `A9_lean_target_screen.md` table (col `safe`), `data/a9/t1_hunt.jsonl` |
| A30 safe floors + rung pass | `uv run python scripts/a30_coset_bz.py validate` (~4 min), `... decide --threads 10` (~15 min); `data/a30/rungs_*.json` |
| A32 tower | `scripts/a32_*.py`, `data/a32/`; §8 verification map (counts per the A33 §8 erratum) |
| Tier C-ii: H5 both ways | `scripts/a33_h5_close.py` (66 s) and `a33_h5_descent.py` (28 s); `data/a33/{h5_direct,h5_descent}.json`, `seam_census.jsonl`, `seam_rungs.jsonl` (1,680 rows each) |
| Tier C-ii: H1/H2 + solver-free upgrade | `a33_validate_banked.py`, `a33_solver_free.py`; `data/a33/{h2_rungs.jsonl, banked_validation.json, solver_free.json}` |
| Tier C-ii: rank lattice, (R) contrast | `a33_tower_cells.py` Part 2; `data/a33/tower_cells.json` |
| Lean capstones | QECLean `main`: `Gross/Distance.lean`, `Z3Z6/Distance.lean`, `Z5Z15F2A6/Distance.lean` |
