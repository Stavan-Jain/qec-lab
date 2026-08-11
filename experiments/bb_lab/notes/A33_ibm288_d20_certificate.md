# A33 — The tower slice calculus on IBM class Y: H5 closed, d([[288,8,20]]) = 20 at certificate tier

**Session 2026-08-11** (worktree `worktree-bravyi360-fibering-fit`). Mission:
port the A32 tower slice calculus to the A20 line — IBM's class-Y
`[[288,8,20]]` (Z₁₈×Z₈, A = 1+xy⁴+x¹⁴y, B = 1+xy²+x²y⁷, d = 20 MILP-exact
per arXiv:2606.02418 Table II) — close the one remaining hypothesis (H5,
the lift-aware seam floor) of its d = 20 theorem, upgrade the whole theorem
to certificate tier, and produce the portability scorecard A32 §7 promised
("next instance #1"). All three landed:

    d([[288,8,20]]) = 20 — certificate tier, end-to-end, no SAT on the
    critical path (X side; Z side by the BB transpose duality, spot-
    asserted at Y8).  Critical-path compute: ~105 s deterministic.

This is the calculus's first test on an **all-(R), same-axis** tower (A32
proved d = 24 on a tower with (R) failing everywhere and mixed axes; here
both rungs are y-decks with Bezout witnesses). §6 documents exactly which
parts of the calculus care. Falsify-first throughout: every banked A20
artifact was re-verified under this port's conventions before use (§7), and
the port's census audit found and patched a latent enumeration edge in the
A32 tooling itself (§8). NOT Lean-checked (follow-on 1, §9).

Scripts `a33_*.py`; data `data/a33/`. Read-only inputs: MAIN checkout
`data/a20/`; worktree `data/a20_ro/`, `data/a32/`.

## §0 Starting state (A20 §7 recap) and the open item

Tower (A20 §1, kernel-verified): Y2 (18,2) `[[72,8,6]]` →y→ Y4 (18,4)
`[[144,8,10]]` →y→ Y8 (18,8) `[[288,8,20]]`; both upper rungs twisted
lifts; (R) holds on both y-rungs (Bezout witnesses banked); k = 8
everywhere; rung 1 attains the deficit wall (10 = 2·6−2), rung 2 doubles
perfectly (20 = 2·10). A20's d = 20 ledger: H1 (dangerous census, 1,655
classes) DONE at V7 grade; H2 (per-class floors at 20) DONE solver-grade;
H3 (LogicalFloor 10) DONE solver-grade; H4 (DeckTrivialOnH1) DONE (Bezout);
H6 (weight-20 witness) DONE. **H5 open**: the naive SeamCosetFloor 20 is
FALSE (weight-18 element in class 0x1's coset), and the per-element
replacement — census both seam-orbit cosets to 19 and certify every
element's lift fiber — had banked partials only.

**Banked-state correction** (the session brief's numbers were stale): the
SAT census banked **278** stab-orbit elements for 0x1 (all lift-UNSAT,
≈ 12.7 h across two runs, census not exhausted) and **0** for 0x3 (its
probe: 3.8 h BUDGET, inconclusive). The dying session's last artifact — a
V7-style seam census (`v7_seam_*.jsonl`, not in the A20 note) — had
already found 280 classes for 0x1 and measured 0x3 EMPTY at ≤ 19; both
verdicts are independently confirmed by BZ here (§2), retroactively
validating that engine run.

## §1 The calculus port (`a33_tower_cells.py`, 0.2 s; `a33_rung_cell.py`)

Frames in the lab convention (H_X = [M_A|M_B] rows; X-logicals =
ker H_Z ∖ rowspace H_X); folds p₂: C₁(Y8) → C₁(Y4) (y mod 4), p₁: C₁(Y4)
→ C₁(Y2) (y mod 2); sheet embeddings, τ = e₀+e₁ per rung. A32's Lemma 1
(deck transport, twist-generic), Lemma 2 (parity; trinomials), Theorem 3
(two-level slice |v| = |β| + 2(m₁+m₂)) all asserted numerically: chain
maps, stabilizer transport, twist-invariance (im S ⊆ ker E), sections,
200 random two-level slice/carry checks, converse lifts. **Same-axis
degeneracy**: the overflow square is vacuous — there is only one descent
order (Y8 → Y4 → Y2); asserted instead: the composite fold equals the
unique Z₄-deck fold. The composite Z₄ structure is never otherwise used —
**iterated Z₂ suffices for everything** (mission question answered).

**The H1 rank lattice under (R)** (all measured, the port's structural
payoff):

- rank p\* = rank τ\* = **4 = k/2 at both rungs**; exactness both ways
  (im τ\* = ker p\*, ker τ\* = im p\*); σ\* = id on H1 at both rungs
  (= H4, re-measured directly; the banked Bezout witnesses 1+y⁴ ∈ (A₈,B₈)
  and 1+y² ∈ (A₄,B₄) also re-verified as convolution identities).
- SEAM := im p₂\* ⊂ H1(Y4): dim 4, 15 nonzero classes in **2 G-orbits
  (12 + 3), class-stabilizer orders 6 / 24** — A20 §5 reproduced. The
  seamC dictionary: the 15 classes of ι(seamC ζ) (A20's chain-frame
  ker-∂₂ parametrization, mapped by the antipode ι) equal SEAM ∖ 0
  exactly; chain-mask 0x1 ↦ class 0xe1 (12-orbit), 0x3 ↦ 0x73 (3-orbit).
- Rung-1 lattice: SEAM1 := im p₁\* also dim 4, orbits 12+3;
  **W₂ := p₁\*(SEAM) has full dim 4** (p₁\* injective on SEAM) and
  **SEAM ∩ ker p₁\* = 0** — one rank computation that collapses the
  descent trisection (§3).
- Contrast datum (A32's non-(R) top rung, measured side-by-side):
  rank τ_y\* = 6 (= rank p_y\*, no drop), ker τ_y\* ≠ im p_y\*, σ_y\* ≠ id.
  Exactness-at-the-cover and deck-triviality are exactly the (R)-bound
  parts of the lattice.

**Rung engine** (generic same-axis port of the a30/A32 architecture):
sector space ker E / im S ≅ H1(base) = 2⁸ at both rungs (κ = 68 top / 32
bottom); dangerous rung with full sector dispatch + restricted MITM lanes
to size 6 + coset-BZ lane; **seam rung** = feasibility-only (every cover
cycle over a nonzero-class element is automatically a nontrivial logical
by stabilizer transport — A32's "pure feasibility" observation held
verbatim; the non-stab property of every enumerated chain is asserted as
a convention tripwire). Engine validations: linear-reduction sector scan
== direct reduction 256/256; planted control — v = τ₂(u₁₀) + H_X(Y8)row,
a genuine weight-26 nontrivial dangerous logical over the weight-6 shadow
with overflow 10 — **FOUND by the BZ lane at M = 11** (min overflow
exactly 10, every violation ≥ 20; 501 s, one-time; the production path
never uses this lane).

## §2 Validation against banked ground truth + H5 closed DIRECT

**Part 3 validations** (`a33_validate_banked.py`, 2.1 s): all 1,655
census rows load as lab Y4-stabilizers (hard assert, full file);
v7_complete_classes == SAT census as G-canonical class sets; rung-1
decomposition of the census (transport + slice asserts all pass; 602
distinct nonzero Y2-shadow orbits, 2.7× compression — the A32 "census
compression is modest" finding repeats; 2 β=0 records = the τ₁-family);
**H2 re-derived: 1,655/1,655 deterministic rung PASSes at target 20 in
1.9 s** (lanes: hexagon M=7 0.3 s, w10 M=5 0.03 s, the rest sub-ms;
banked SAT agreement 1,655/1,655); H6 witness re-verified (tower profile:
b = 0, τ-diagonal over a weight-10 Y4 logical whose own rung-1 shadow is
flat |β| = 10, m₂ = 0); y144/y72 ladder witnesses re-verified.

**H5 direct** (`a33_h5_close.py`, 66 s): one 3-offset coset-BZ pass at
n = 144 (κ = 68, W = 18 by parity, r-pair (9,8), 6.62e10 nodes, exact
node-count asserts, 61 s):

| offset | result |
|---|---|
| stabilizers (t₀ = 0) | 118,932 vectors ≤ 18, weights {6: 72, 10: 432, 12: 2,268, 14: 3,888, 16: 26,928, 18: 85,344} → G-canonical classes = **1,655 == banked EXACTLY** (H1 re-derived DIRECT; the vector count equals A20 §6's V7 accounting 118,932) |
| 12-orbit seam coset (class 0x6) | **1,680 elements** {14: 6, 16: 84, 18: 1,590}, every one verified a class-0x6 logical; 280 G-canonical classes |
| 3-orbit seam coset (class 0x73) | **EMPTY** — the naive SeamCosetFloor 20 is TRUE on this orbit (H5 vacuous there; explains the banked 0x3 SAT BUDGET) |

Banked cross-checks: the 278 SAT-census elements (chain frame, verified
in the chain-frame coset, then ι-mapped) are all present — the dying
census was **2 classes short of 280**; the v7_seam 280 classes equal the
BZ census exactly; the SAT@18 seam_floor witness is present.

**The seam rungs: 1,680/1,680 PASS in 1.8 s** (lanes: 1,590 flat-top
M=1, 84 at M=2, 6 at M=3; zero violations). Banked agreement: all 6
class-stab translates of each of the 278 banked lift-UNSAT elements PASS
(1,668 elements). Covariance spot-check: translated elements give
identical verdicts. All 15 classes follow by G-transport (the fold
G(Y8) ↠ G(Y4) is onto, so every class translates to an orbit rep).

Tightness probes (find-side control + slack measurement): over a probed
w18 element the lightest cover logical weighs 22 (overflow 2, first hit
at M = 3); over a probed w14 element also 22 (overflow 4). Both ≥ 20;
the seam sector does not attain d — the weight-20 witness lives in the
τ-diagonal sector.

**H5 is closed.** With H1–H4, H6 already done, the A20 §7 assembly is
complete: for any nontrivial X-logical v of Y8 with |v| ≤ 18 (even by
parity), b = p₂(v) is (i) 0 → v = τ₂(u), [u] ≠ 0, |v| = 2|u| ≥ 20 by H3;
(ii) a nonzero stabilizer → excluded by H1+H2 (|b| ≤ 18 censused classes,
rungs at M = (20−|b|)/2; |b| ≥ 20 direct by slice); (iii) of nonzero
class → [b] ∈ SEAM ∖ 0 (= p₂\*[v]), excluded by the H5 census + rungs
(12-orbit) or by coset emptiness (3-orbit). Hence d ≥ 20; = 20 by H6. ∎

## §3 H5 closed by DESCENT (`a33_h5_descent.py`, 28 s)

The A32 Part-5 pattern one rung down — and the (R)-collapse makes it
one-branched: SEAM ∩ ker p₁\* = 0 kills the β = 0 (τ-diagonal) and
β-stabilizer branches outright, and p₁\*-injectivity on SEAM pins the
shadows of class-0x6 elements to the **single Y2-class 0x36** (0x73 →
0xf3). This is the sharpened analogue of A24-style reachability pruning:
one rank computation, measured in §1, answers "which Y2-classes can seam
shadows occupy" — exactly one per seam class, and never the trivial one.

- Y2 coset censuses (n = 72, κ = 32, W = 18, r-pair (9,8)): class 0x36 —
  105,328 elements {10: 18, 12: 162, 14: 1,470, 16: 11,910, 18: 91,768};
  class 0xf3 — 81,108 elements. (A live instance of the empty-window
  census edge, §8: the 0x36 window-1 coset base has weight 18 and is a
  census member — handled explicitly.)
- Lift fibers over every censused shadow (a32 MITM enumerators reused,
  deep lane gate-validated == size-4 lane): 12-orbit — 105,328 fibers,
  **90.6% carry-infeasible**, 38,714 lifts, exactly **1,680 in SEAM ==
  the direct census EXACTLY** (sheet-flip expansion via the deck
  translate (0,2); p₁\*-injectivity asserted per lift: any SEAM lift is
  class 0x6). 3-orbit — 81,108 fibers, 8,460 lifts, **0 in SEAM**: the
  direct emptiness re-derived from n = 72 data.
- Identical element set ⟹ identical rung verdicts (§2's 1,680/1,680);
  50-element re-rung sample 50/50 PASS.

## §4 The full solver-free upgrade (`a33_solver_free.py`, 6.6 s)

Every remaining SAT-tier input on the critical path re-derived:

- **d₀ = 6 direct BZ** (n = 72, all 255 classes to weight 6, 6.0e3
  nodes): 84 weight-6 logicals in 24 classes, none lighter (parity kills
  7); banked y72 witness present. 0 attaining classes in SEAM1.
- **H3 = LogicalFloor 10 direct BZ** (n = 144, all 255 classes to weight
  10, 1.2e7 nodes, 0.6 s): no element below 10; 72 weight-10 logicals in
  12 classes; banked u10 and the H6 witness's sheet both present. **0
  attaining classes in SEAM** — "light logicals concentrate outside
  im p\*" repeats at both rungs of this tower (5th and 6th instances of
  the A24/A32 pattern).
- **d₁ = 10 by the tower** (the a32_dby pattern at budget 8): the
  seam1-branch is *census-empty* (the 15 seam1-cosets have no element
  ≤ 8 at all — 0 fibers to run); the β = 0 branch dies by parity + d₀;
  the stabilizer branch closes inside the fiber-union sweep below (min
  non-stab lift weight = 12 ≥ 10). Witness = (d₁ ≤ 10) from the direct
  weight-10 census.
- **H1 census-completeness by the tower** (the census bonus): Y2
  stabilizer census ≤ 18 (165,517 vectors, {6: 36, 8: 36, 10: 252,
  12: 1,368, 14: 5,436, 16: 26,721, 18: 131,668}; µ(Y2) = 6 — note the
  weight-8 gap does NOT hold at Y2, an anti-instance of the Prop-10
  pattern that held at Y4/gross/BY/GB) → 4,605 orbit reps → bounded-
  overflow fibers (caps (18−|β|)/2) → 8,302 lifts (3,328 stabilizer,
  4,974 logical all ≥ 12) + the τ₁-family (τ₁ of the 72 Y2-stabs ≤ 8;
  the seam1 γ-family is empty) → class-key union = **the banked 1,655
  EXACTLY** (3.2 s). H1 now has four independent derivations: SAT
  census, V7 analytic, direct BZ (§2), tower fiber union.

**Hypothesis ledger after this session** (all deterministic-certificate
tier, no SAT):

| # | statement | certificate |
|---|---|---|
| H1 | dangerous census exhaustive | direct BZ (118,932 vectors → 1,655 classes) + tower fiber union + [SAT, V7 as historical cross-checks] |
| H2 | per-class dangerous floors at 20 | 1,655 deterministic rungs, restricted lanes |
| H3 | LogicalFloor 10 at Y4 | direct BZ 255-class census + tower re-derivation |
| H4 | DeckTrivialOnH1 | Bezout convolution identities re-verified + σ\* = id measured on H1 |
| H5 | lift-aware seam floor at 20 | 12-orbit: BZ census (1,680) + 1,680 rungs, direct AND by descent; 3-orbit: coset empty |
| H6 | weight-20 witness | re-verified nontrivial X-logical, tower profile τ-diagonal |

Claim tier: identical to A32/A30 — BZ censuses with exact node-count
invariants, MITM fibers complete by the exact-subset-sum argument with
every solution re-verified, rungs with end-to-end re-verification, plus
external agreement with every banked SAT verdict (1,655 H2 + 278 H5 + 3
ladders). Z side by transpose duality (antipode + block swap,
spot-asserted at Y8). NOT kernel-checked.

## §5 Cost model: before → after (measured)

| obligation | before (A20 architecture) | after (this session) | factor |
|---|---|---|---|
| H5 element census + lift floors | SAT: ≈ 12.7 h for 278/280 classes of ONE orbit (incomplete) + 3.8 h inconclusive 0x3 probe ⇒ ≈ 16.5 h partial | 66 s direct (both orbits, complete, incl. the H1 census) + 28 s descent cross-derivation | ≈ **600×**, and partial → complete |
| H2 dangerous floors | 1,655 SAT UNSATs, ≈ 17 s (solver tier) | 1,655 deterministic rungs, 1.9 s | tier upgrade + 9× |
| H1 census | SAT enumeration (overnight) / V7 analytic 1,105 s | direct BZ 61 s (shared pass) + fiber union 3.2 s | ≈ 18× vs V7 |
| H3 | SAT UNSAT@9 ladders, 180 s | direct BZ 0.6 s + tower | 300× |
| d₀ | SAT ladder 0.5 s | BZ ≈ 0.1 s | — |
| whole d = 20 critical path | mixed SAT + analytic, H5 open | **≈ 105 s deterministic, end-to-end** | — |

One-time engine validation (planted control through the BZ lane): 501 s.

## §6 Portability scorecard (the mission deliverable)

**Same-axis deltas.**
1. The overflow square (A32 Thm 3's path-independence) **degenerates**:
   one descent order exists. Nothing else changes — no calculus step
   consumed the square.
2. Iterated-Z₂ vs composite-Z₄: **iterated Z₂ sufficed for everything**;
   the composite fold appears only in a consistency assert. The
   trisection, fibers, and rungs never needed the Z₄ deck.

**(R) vs non-(R) — what actually changes** (all measured, §1):
3. Sector space = 2^k at every rung on both towers (A32's "2⁸" is
   k-genericity, **not** an (R) effect — both towers have k = 8 bases).
4. (R) shows up in the H1 lattice: rank τ\* drops to k/2 with exactness
   both ways and σ\* = id; im p\* is the k/2-dim SEAM with the
   ker-∂₂/seamC dictionary. On A32's non-(R) top rung all three fail
   (rank τ_y\* = 6, no exactness, σ\* ≠ id).
5. Consequence for the descent: the trisection **collapses to one
   branch** (SEAM ∩ ker p₁\* = 0; p₁\*|SEAM injective ⟹ one target
   coset per seam class). A32's three sectors A/B/C are the non-(R)
   general position; under (R) two of three die by one rank computation.
6. Consequence for the safe sector: the per-element obligations are
   flat-top-dominated (94.6% of the 1,680 rungs at M = 1) because seam
   minima sit at 2d₁−2; A32's safe strata ran to M = 5.

**A32 §7 "fully portable" claims, verdict.**
7. Lemma 1 + Theorem 3: **held verbatim** (both rungs twisted; the
   calculus never noticed — transports asserted).
8. Trisection shape: held; under (R) it *simplifies* (above).
9. Lemma 2 (parity): held (trinomials; every distance even).
10. The iterated-rung engine + both bonus patterns: **held verbatim** —
    the middle code's distance floor re-derived by the same fibers at a
    smaller budget (d₁ = 10, §4), and the middle code's stabilizer
    census re-derived as the fiber union (H1, §4). "One machine serves
    the distance floor, the census, and the sector floors at every
    rung" — confirmed on a second, structurally different tower.
11. Touch-up found: the `a32_gb_census.census()` helper's empty-window
    edge (§8) — my census wrappers handle coset-base elements
    explicitly; **the edge was live in this port** (a Y2 window base at
    weight 18 was a real census member).
12. New species this port adds: the feasibility-only seam rung
    (per-element, no sector dispatch) and the 3-offset combined BZ pass
    (stab census + seam cosets in one kernel run).

**Updated next-instance ranking**: (1) the A30 `[[360,4,20]]` doubled
codes' rung-2 re-doubles — three-level towers, all-(R), same-axis: this
session's exact shape, with the A14 §13 "freeze at the toric bottleneck"
now re-examinable per-element; (2) `[[756,16,≤34]]` y-rung (n = 378
base; tower census the only plausible route); (3) gross over bb72 over
(3,6) retrospective (teaching-doc chapter); (4) revisit A24's odd-|G|
instances where no second deck exists (the calculus cannot apply —
Φ-transfer/λ-certificates remain the only lane there).

## §7 Verification map

| claim | check |
|---|---|
| Lemma 1 both rungs (twisted) | `a33_tower_cells.py` Part 0 + TS.Deck constructor asserts (chain map, stab transport, twist-invariance im S ⊆ ker E, sections) |
| parity exhaustive | Part 1 (kernel bases even, Y2/Y4/Y8) |
| two-level slice/carry | Part 1: 200 random cycles + composite-shadow assert + converse lifts |
| rank lattice / exactness / σ\* = id / SEAM orbits / dictionary / Bezout | Part 2 (all hard asserts; A20 §5 orbit structure reproduced; A32-tower contrast measured side-by-side) |
| banked census loads + rung-1 decomposition | `a33_validate_banked.py` (1,655 hard asserts; v7 == SAT class sets) |
| H2 floors | 1,655 rungs PASS + 1,655/1,655 banked-SAT agreement |
| rung engine soundness | linear sector trick 256/256; planted weight-26 dangerous logical FOUND at its exact overflow (BZ lane); every violation-candidate re-verified in-line (slice identity + non-triviality asserts) |
| H5 census completeness | 3-offset BZ, exact node counts (6.62e10), r-pair (9,8) complete to 18, empty-window bases handled; 118,932-vector H1 side re-derives banked 1,655 as an in-pass control |
| H5 census vs banked | 278/278 SAT elements present (chain-frame coset membership asserted before ι); v7_seam 280 == exact; SAT@18 witness present |
| H5 rungs | 1,680/1,680 PASS; per-chain non-stab asserts; covariance spot-check; G-transport via fold surjectivity |
| H5 descent | element-set equality with direct (1,680 == 1,680; 0 == 0); deep-enumerator gate (== size-4 lane); 50-element re-rung |
| d₀, H3 | direct BZ censuses (255 classes each); banked witnesses + H6-sheet membership asserts |
| d₁ tower re-derivation | seam1-census empty + fiber-sweep min non-stab lift 12 + parity/d₀ kill of β = 0 |
| H1 fiber union | class-key equality with banked 1,655 (orbit-rep fibers, fold-surjectivity completeness) |
| H6 | witness re-verified + tower coordinates |
| a32 census-edge audit | §8: per-species empty-window audit + orbit-membership checks + cap-0 fiber patches |

## §8 Falsified-claims ledger + errata (session-internal)

- **Mission-brief stale numbers**: "banked 260 (0x1) + 270 (0x3), both
  censused" — FALSE; actual banked state was 278 (0x1, incomplete, all
  UNSAT) + 0 (0x3, probe inconclusive), plus the unrecorded v7_seam
  engine run. No banked *verdict* was wrong.
- **"Element-census completeness was never settled" — now settled three
  ways** for 0x1 (BZ direct, descent, v7_seam agreement) and the 0x3
  coset is EMPTY ≤ 18 — the naive floor A20 §5 declared FALSE-as-stated
  is actually TRUE on the 3-orbit; it fails only on the 12-orbit.
- **A32 tooling erratum (found by this port's audit, patched)**: the
  `census()` helper never emits the empty-window coset-base element.
  Audit of every A32 GB census species: logical ≤ 10 — 15 weight-10
  vectors missed, but **all 15 lie in stored orbits** ⟹ every per-orbit
  consumer (sector B, dby floor, d(GB) = 8) unaffected; vector-count
  erratum only (1,623 → 1,638). W-coset ≤ 16 — clean. W-coset ≤ 22 —
  2 weight-22 elements missed, both orbits enumerable via translates
  (orbit tables unaffected), and both cap-0 fibers are EMPTY ⟹ **A32's
  d = 24 assembly is unaffected** (patched directly regardless). Raw
  vector counts in A32 §3's table carry +15/+2 corrections.
- "The seam sector might attain d = 20" — measured NO on probed
  elements: lightest logicals over w18/w14 elements weigh 22; the
  witness sector is τ-diagonal. (Not a general theorem; two probes.)
- The Prop-10 weight-8 stabilizer gap does **not** extend to Y2
  (w8 stabilizers exist, 36 of them) — the gap held at Y4/gross/BY/GB
  only; do not cite it as universal.
- (Inherited, respected): no spectral kill engines, no witness-weights-
  as-floors (probe weights reported as *found logicals*, not minima —
  except where the census IS complete), no upper-as-lower.

## §9 Residue and next steps

1. **Lean packaging** (nothing here is kernel-checked). All species have
   shipped analogues (A28 §6.3 BZ counting certificates, KernelCert
   pivot certificates, A15/A30 rung dispatches); the new species — the
   3-offset BZ pass and the feasibility-only seam rung — are strictly
   simpler than their dangerous-sector counterparts. The 1,680-element
   seam layer is small enough to data-carry.
2. **A20 note bookkeeping**: §7 ledger updated in place (H5 → DONE,
   pointer here); the "two verdicts pending" draft marker resolved by
   the dated addendum.
3. **Deficit-wall postscript**: the safe-sector minima sit at
   2d₁ − 2 = 18 on the live orbit (the wall value, again) while the
   other orbit is empty — the wall phenomenon is per-orbit, not
   per-code. Feed to the A17-P3 taxonomy.
4. **Next instance** (§6 ranking): the A30 `[[360,4,20]]` rung-2
   re-doubles — the freeze-at-2d̄−2 verdicts (A14 §13) deserve a
   per-element re-examination with this session's machinery.
5. **Witness polish**: bank one explicit weight-20 logical + its
   tower profile as the canonical d ≤ 20 certificate (currently: the
   re-verified `y8_weight20_witness.npy`).
