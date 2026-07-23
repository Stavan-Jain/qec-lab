# A19 — Bravyi [[360,12,≤24]] as a deck-nontrivial tower: instance study and path to certification

**Session 2026-07-20.** Origin: methodology review of IBM's LLM-guided BB-code
discovery paper (arXiv:2606.02418) surfaced Bravyi et al.'s [[360,12,≤24]]
(arXiv:2308.07915, Table 3: `(ℓ,m) = (30,6)`, `A = x⁹+y+y²`,
`B = y³+x²⁵+x²⁶`) as the natural gross-playbook target — highest published
FOM, distance unresolved (IBM's deep MILP: incumbent only; "≤" in Bravyi's
table = upper bound known only). This note records the instance study, the
corrected structural picture, and the staged plan toward certified distance.

Scripts: `scripts/a19_*.py` (run from `experiments/bb_lab` via `uv run`).
Data: `data/a19/` (ISD per-class minima npz, SAT round logs).

## 1. The tower and its exact base distances

```
GB (15,3)  --x-->  BY (30,3)          A = x⁹+y+y² throughout
   |                  |               B descends by exponent reduction:
   y                  y               BX: y³+x¹⁰+x¹¹   BY: 1+x²⁵+x²⁶
   v                  v               GB: 1+x¹⁰+x¹¹
BX (15,6)  --x-->  C  (30,6) = the Bravyi code
```

| code | frame | (n,k) | d | provenance |
|---|---|---|---|---|
| GB | Z₁₅×Z₃ | [[90,8]] | **8 exact** | coset-SAT ladder, 35 orbit reps, 187 s (`a15_coset_distance`); witness support `[12,14,27,29,69,71,84,86]` |
| BY | Z₃₀×Z₃ | [[180,8]] | **12 exact** | floor@11 all-UNSAT (35 reps, 3.1 h, `data/a19/by_floor11.log`) + weight-12 witness (`by_floor13.log`, class 0x10) |
| BX | Z₁₅×Z₆ | [[180,8]] | **12 exact** (2026-07-22) | floor@11 all-UNSAT: 14/35 banked (`bx_floor11_partial.log`) + 21 resumed (`a19_bx_resume.py`, 4.1 h, `bx_floor11_resume.log`) + weight-12 ISD witness (`a19_lifts.py`) |
| C | Z₃₀×Z₆ | [[360,12]] | **≥ 12 CERTIFIED (M12, 2026-07-22)**; ≤ 24 constructive | §3, §4, §8 |

- GB is **not** Bravyi's bb_90 ([[90,8,10]] shares the frame and A; different
  B, different d). BX ≁ BY under the full Aut+monomial orbit
  (`a19_bx_by_equiv.py`: distinct canonical pairs).
- Both lower rungs are k-preserving (8→8) with **deficit 4**:
  d = 12 = 2·8 − 4 on *both* axes — two new natural deficit cells for the A17
  docket, again strictly below the wall value 2d−2, and the first
  same-base/different-axis pair where both axes are deficit-equal.
- All three bases pass `a15_class_certify` in full (D1, D2, (iii), FRM, ANN,
  FLR; one-sided engine floors 20/20/10) ⇒ certified d ≥ 6 each, and by the
  A16 transfer clause **the Bravyi cover inherits certified d ≥ 6** — its
  first certified lower bound by any method.

## 2. Deck-nontriviality (the headline structural verdict)

k jumps 8 → 12 into the cover. Checked exhaustively (`a19_orbit_descents.py`,
`a19_diagonal_descent.py`): all **three** decks (x¹⁵, y³, and the diagonal
x¹⁵y³ via `φ(u,v) = (u+15v mod 30, v mod 3)`), each across the full monomial
orbit (units a ∈ (Z/30)*, b ∈ (Z/6)*) — **no k-preserving descent exists**.
By A12 ((R) ⟺ k-preservation ⟺ Bezout) and A10-L1 (descent verdicts are
code-level), the Bravyi flagship is **genuinely deck-nontrivial on every
deck**: the first natural flagship instance outside the doubling template's
condition 2, after 157/157 historically checked covers satisfied it.

Quantitative structure (`a19_cover_census.py`, `a19_seam_sector.py`, all
kernel-level F₂ linear algebra):

- `dim (1+σ)H₁ = 4 = k̃ − k` for **all three decks** — the A12 open
  quantitative conjecture verified on a natural instance, three ways.
- Pushforward class maps to BX and BY each have rank 6 (kernel 6); the
  doubly-new sector `ker p_x* ∩ ker p_y*` has dim exactly 4, and it **equals
  im(1+σ)₊ for every deck** — the A13 D-module summand, deck-independently.
  Sector census over 155 translation-orbit reps: 138 old / 6 new-x / 6 new-y /
  5 new-xy.

## 3. The weight table is transfer-lift arithmetic (d ≤ 24 constructive)

Deep ISD over ker H_Z (54k iterations, all 4095 classes, vectors saved in
`data/a19/isd_class_minima.npz`):

| sector | min wt seen | witness class | mechanism (verified) |
|---|---|---|---|
| old | 24 | 0x6e2 | section-type: projects to nontrivial base logicals wt 24 (BX) / 22 (BY) |
| new-x | 24 | 0x41d | **τ_x(BX weight-12 logical)** — verified explicitly (0x570) |
| new-y | 24 | 0x002 | **τ_y(BY weight-12 logical)** — verified explicitly (0x380) |
| new-xy | 32 | 0x80c | **τ_x∘τ_y(GB weight-8 logical)** — verified explicitly (0xc04) |

(`a19_lifts.py`: each τ-lift checked to be a cycle, nontrivial, of weight
exactly 2·(base weight), landing in the predicted sector.) The transfer lift
τ = full-preimage projects to **zero** (sheets cancel mod 2) — τ-lifts live in
the *new* sector of their deck, `im τ_* ⊆ ker p_*`. So the cover's light
operators are the tower's arithmetic: 24 = 2·12 (both decks), 32 = 4·8.
**d(C) ≤ 24 is now constructive** (explicit verified operator), independent of
solver incumbents. A budgeted beat-24 SAT gate (w ≤ 22 per-coset over the ten
lightest orbit reps, `a19_beat24.py`) died with its host session having
returned zero results — retired for good per the epistemics below; 54k ISD
iterations never beat 24. **Epistemics of the gate:** n = 360 SAT is witness-side only — if the
true coset minimum is 24, a w ≤ 22 query is secretly UNSAT and times out
uninformatively (per-coset UNSAT is priced out beyond ~w = 13 at n = 180,
~w = 14 at n = 300; hopeless at n = 360). Hits are decisive, timeouts are NOT
evidence. The *real* beat-24 instrument is at base scale and coincides with
the certification program: a dangerous-sector operator below 24 forces a
base stabilizer with |b| + 2m(b) < 24 (slice identity, (R)-free) — so the
(M)@24 census is simultaneously the lower-bound machine and the localized
sub-24 hunt; old-sector candidates reduce to merges of SAT-enumerable
weight-12–16 base logicals (n = 180, minutes). Blind n = 360 SAT is retired
from the program after this gate's budget expires.

## 4. Gap analysis for certified lower bounds (the corrected picture)

Audit of `docs/gross-distance-proof.md`: **(R) is consumed only by the
safe-sector confinement (Theorem D)**; the dangerous-sector machinery
(slice identity + light-stabilizer inequality (M), Theorem C / Prop 10) never
invokes it. Consequences for this (R)-failing cover:

- **Portable now:** for any cover logical v with `p_y*[v] = 0` (the new-y and
  new-xy sectors), the slice identity over BY applies verbatim:
  `|v| ≥ |b| + 2m(b)`, `b = p_y(v) ∈ Stab_Z(BY) ∪ {0}`, plus the zero-rung
  for b = 0. Symmetrically over BX.
- **Milestone M12 (near-term):** (M)@12 for BY — every nonzero
  `b ∈ Stab_Z(BY)` with `|b| ≤ 11` has `|b| + 2m(b) ≥ 12` — is *gross Prop-10
  depth* (classification to weight 11: hexagons |A|+|B| = 6, D-pairs, …), at
  n = 180. Combined with the projection floor min(d(BX), d(BY)) = 12 on the
  other three sectors, it yields **certified d(C) ≥ 12** — double the A16
  floor, and stronger than any published bound for this code.
- **Full-24 program:** (i) (M)@24 for BY/BX (classification to weight 23 —
  substantially deeper than gross's 11; SAT-assisted census first — now via
  the shipped A17 (M)-kernel pipeline, see §6);
  (ii) zero-rung at 24 = 2·d(base) (τ-tightness: the dangerous floor is
  *attained* by the §3 lifts, exactly the gross pattern); (iii) the
  **doubly-old (safe) sector ≥ 24 without (R)** — *the route now exists*
  (update 2026-07-22): the seamC↔δ₂ transport closed on 2026-07-20
  (`BBBocksteinTransport.lean`, axiom-clean —
  `bocksteinVanishes_of_elementForm`; `BocksteinElementForm` discharged
  **unconditionally for every ZMod doubling deck, twists included**, which
  covers all three decks of C: the diagonal x¹⁵y³ is the twisted x-deck
  (15,3)). The rank scaffolding is now theorems — `E = k̃ − k`,
  `dim ker ε₊ = k`, `ker τ₊ ≤ range p₊` — so §2's numerics (E = 4, ×3 decks)
  are instances of proven statements. The remaining *structural* piece is the
  iso `H₁ ≅ D^{k̃−k} ⊕ F₂^{2k−k̃}` (f.g. `F₂[ε]/(ε²)`-module classification —
  bounded algebra; all rank inputs proven). Sequencing: A20's
  `SeamCosetFloor 20` is the transport machinery's **first live trial** at
  half scale (odd part Z₉ ⇒ F₆₄ components there, vs F₁₆ here); the cheap
  A19 pilot — doubly-old classes have both projections nontrivial (≥ 12
  each), test whether a joint-support argument clears 16–18 — proceeds in
  parallel, then A20's port pattern comes up-tower.

Record context: the published solver-exact BB record is Bravyi's
[[288,12,18]] (MIP-exact per their Table 3 note; independently E-status in
IBM's Table V). Full certification at [[360,12,24]] would exceed it by 25% in
n and 33% in d; the kernel-verified record (gross, n = 144) by 2.5× in n.

## 5. Falsified-claims ledger (session-internal; doctrine: record them)

- "The 360 code is a doubling tower with d = 4·d(GB)" — dead: (R) fails at
  the top rung, all decks, all presentations.
- "d(BY) = 14" (first ladder witness) — dead: true d(BY) = 12; the witness
  was non-minimal. Same error nearly repeated at BX (MC said 16; ISD found 12
  in 60 s). **Protocol lesson: exhaust cheap witness hunts (ISD with
  early-exit) before buying UNSAT rounds; MC-phase-A of `a15_coset_distance`
  is far weaker than 60 s of ISD.**
- "The cover minimum is an old-sector diagonal lift at 2·d(BY)" — mis-framed:
  τ-lifts project to zero and live in the *new* sectors (§3); the old-sector
  24s are section-type. (Conflated section lift s(u) — weight-preserving,
  p∘s = id, not a cycle in general — with transfer lift τ(u) — weight-2|u|,
  p∘τ = 0, always a cycle.)
- "[[288,12,18]] sits at the wall (its (6,12) base has d = 10)" — dead
  (2026-07-20 addendum): floor@9 + floor@11 all-UNSAT over all 155 orbit
  reps certify d(base) ≥ 12 (§7 verdict).
- "IBM's largest exact is n = 288 at d = 12" — wrong: that is their largest
  *own-discovered* exact; the published record incl. Bravyi's codes is
  [[288,12,18]].

## 6. State and next steps — RE-RANKED 2026-07-22 (post-transport)

Program shifts since the opening session: (a) the A13 seamC↔δ₂ transport
closed — the §4(iii) safe-sector route is scaffolded by theorems, not a
conjecture; (b) the A17 (M)-kernel census pipeline **shipped**
(pivot-certificate sweeps: 2^L enumeration → Gaussian-elimination
certificates, 53 min → 3.9 s) and A20 already reuses it (1,655 classes
floor-certified in ~17 s, 0 SAT hits); (c) A20 stands as the deliberate
half-scale pathfinder — 4/5 rung-2 obligations discharged, `SeamCosetFloor
20` the sole opening and the transport's first live consumer; (d) the
beat-24 gate died with its session, zero hits — retired (§3). Re-ranked
order:

1. **Close BX** (cheapest, hours): targeted resume of the 21 outstanding
   floor@11 orbit reps (14/35 banked all-UNSAT) → both bases exact
   [[180,8,12]]; projection floor 12 on every sector except new-xy. Log the
   two deficit-4 cells in the A17 docket.
2. **M12 → certified d(C) ≥ 12 — DONE (2026-07-22, same session as this
   re-rank; §8).** The pipeline retarget took one census run (227 s) plus
   seven floor queries (0.4 s). Note: the route runs entirely through the
   y-deck — BX (step 1) is redundancy, not a dependency.
3. **(M)@24 census** at BY and BX (depth 23): same pipeline, deeper bands
   — with the **A22 cross-link (2026-07-22)**: A22's (ε,δ) CRT fibering
   re-derived a 113-class enumeration in ~2 min that SAT took 9.6 h, and
   BY carries the same ℤ₅ fiber (z = x⁶, `F₂[z]/(z⁵−1) ≅ F₂ × GF(16)`,
   identical per-site weight table) — the census redesign A20's porting
   lesson called for. Port caveat: A22's base ℤ₁₅ is odd/semisimple; BY's
   base ℤ₆×ℤ₃ carries the ℤ₂ radical, so the ε-side needs the (1+s)-layer
   engine (A5 E2 pattern), not A22's invertible-substitution shortcut.
   Port opened same day (`a19_fibering.py`, §9). A20's band data
   calibrates the SAT fallback (orbit-blocking clause load is its known
   bottleneck).
4. **Doubly-old sector ≥ 24** (moonshot, §4(iii)): sequence *behind* A20's
   SeamCosetFloor-20 trial (F₆₄ there, F₁₆ here); run the joint-support
   pilot (≥ 16–18) in parallel; port the pattern up-tower. The
   module-classification iso, when it lands, supplies the D⁴ ⊕ F₂⁴
   statement for the Lean capstone. **A23 cross-links (2026-07-22)**:
   (a) its instance's hard block and BY's both live over **F₁₆** (ℤ₅
   fiber there, e = 15 here) — the idempotent/character-strata machinery
   is being built in our field; (b) its reduction pattern — a fully
   quantified SeamCosetFloor collapsed to one self-contained convolution
   inequality `|A⋆f + e₀| + |B⋆f| ≥ 16` anchored on the kernel idempotent
   — is the target shape for our sector floors; (c) **negative knowledge,
   binding**: disjoint/fractional packing and region-parity certificate
   families are REFUTED there with a sharp LP-gap (≈13.8 < 15) — do not
   spend budget on packing certificates for A19's floors.
5. **Lean staging** (parallel, engineering): BaseFloors certificates for
   GB/BX/BY (generator-ready; frames floor-bearing); named-hypothesis
   packaging of d = 8/12/12 (Z5Z15F2A6 model); **new post-transport item:**
   instantiate the transport theorems on C — the unconditional ZMod-family
   element form covers all three decks (diagonal = twisted x-deck), making
   E = 4 and the sector accounting Lean theorems for the flagship. (The
   former item "verify the dangerous rungs are (R)-free" is done — §7.)

Open question worth naming: both rungs of this tower are deficit-4, yet the
cover recovers 24 = 2·12 via τ-lifts of the *deficit* value while gaining
k̃ − k = 4 logical qubits that sit heavier (32 = 4·8). Whether
"deficit rungs + deck-nontrivial top" is *why* this code's FOM leads the
catalog — i.e., a design principle (hunt covers that break (R) over certified
bases) — is testable against the IBM catalog's other high-FOM codes.

## 7. Same-session addendum: deck survey of the full Bravyi table

`a19_deck_survey.py` (quotients built by coset convolution — the φ-map trick
does not generalize to even m/2): literal-descent k on every order-2 deck of
every published Bravyi code, orbit-swept where a jump appeared.

| code | decks | verdict |
|---|---|---|
| [[72,12,6]] (6,6) | x, y, xy | **JUMP 8→12 on all three — code-level** (full orbit swept: no k-preserving descent exists) |
| [[90,8,10]] (15,3) | — | \|G\| odd: not a Z₂-cover at all |
| [[108,8,10]] (9,6) | y | R-holds (8→8) |
| [[144,12,12]] (12,6) | x / y / xy | x: R-holds (bb72, the known descent); y and xy: JUMP 8→12 (literal) |
| [[288,12,18]] (12,12) | x, y, xy | **R-holds on all three** (12→12) |
| [[360,12,≤24]] (30,6) | x, y, xy | JUMP 8→12 all three, code-level (§2) |
| [[756,16,≤34]] (21,18) | y | R-holds (16→16; base = [[378,16]]) |

Consequences:

- **The strong design-principle hypothesis is falsified**: [[288,12,18]]
  (FOM 13.5, #2 in the table) is deck-trivial on every deck. But the refined
  pattern holds and is sharper: **every k = 12 flagship acquires its k at
  exactly one deck-nontrivial rung and propagates it through (R)-rungs.**
  Gross and 288 both root at bb72 — the (6,6) quarter of [[288,12,18]] is
  literally bb72 (y⁷ ≡ y mod 6) — and *bb72's own k = 12 is deck-born* from
  k = 8 quarter-codes at (3,6). The 360 code's k = 12 is born at its top rung.
- **Wall-attainment prediction (new A17-docket candidate).** [[288,12,18]]
  is an all-(R) iterated bb72-cover; PAR forces its bases' distances even;
  18 = 2·10 − 2. If its (6,12) base `A = x³+y²+y⁷, B = y³+x+x²` (n = 144,
  k = 12) has d = 10, then [[288,12,18]] sits *exactly at the deficit wall
  2d − 2* — which no measured instance has attained (A17 §corrected picture).
  Ladder is cheap at n = 144; queued behind the in-flight SAT rounds.
  **Verdict (2026-07-20, same day): REFUTED.** The ladder ran — floor@9 and
  floor@11 over all 155 orbit reps, both all-UNSAT (191 s / 21.7 min,
  `data/a19/e288_base_floor{9,11}.log`): **d(base) ≥ 12 certified**, so the
  wall value `2d − 2` is out of reach for this cell. **Closed same day:
  floor@13 returned nine weight-12 witnesses (4.0 h) — d(base) = 12
  EXACT, the cell sits at `18 = 2d − 6`, deficit-6 certified 2-for-2**
  (with bb108-y; A17 §8b). The base is [[144,12,12]] — gross parameters,
  the perfect y-double of its bb72 half — so [[288,12,18]] is bb72
  doubled twice: perfect rung, then deficit-6 rung (mirror of A20's
  wall-then-perfect tower). A catalog-wide deck-birth survey
  (`scripts/a19_deck_survey_catalog.py`, IBM supplemental Tables I–III,
  218 k-validated rows) simultaneously confirmed the refined design
  principle for every k = 12 class (literal) and surfaced two follow-up
  targets: **[[288,12,16]]** (IBM class a, A = 1+y+x, B = y³+x⁵+x¹⁰ at
  (12,12)) — a second deck-nontrivial-top flagship, d = 16 MILP-exact —
  and **[[288,8,20]]** (IBM class Y, (18,8), d = 20 MILP-exact), born at
  (18,2) and tower-structured with bases at n = 144/72: the cheapest
  known d ≥ 20 certification target.
- **Lean landing zone verified**: in `BBDoubling.lean`, `dangerous_zero_rung`
  (l. 441) and the shape rungs consume `StrongBaseFloor`/`LogicalFloor` only;
  `DeckTrivialOnH1` first appears in the assembly theorem (l. 1280). The
  dangerous-side machinery is already (R)-free as formalized — the M12/M24
  program needs no refactor to land. m(b): `gross-distance-proof.md` l. 1059
  (off-support minimum), hexagon rung Lemma 11 l. 1477 — the census port
  starts there.

## 8. M12 CERTIFIED: d([[360,12,≤24]]) ≥ 12 (2026-07-22)

The §6.2 milestone, executed via the A17/A20 (M)-kernel pipeline retargeted
at the y-deck (`a19_m_census.py`, `a19_m_floors.py`; data in `data/a19/`):

1. **Census** (227 s, XOR-native + orbit blocking, bands 2–10, exhaustive):
   BY has exactly **7** translation-orbit classes of nonzero X-stabilizers
   with |b| ≤ 11 — one hexagon (w = 6) and six at w = 10, **no weight-8
   stabilizers** (the gross Prop-10 gap, third instance: gross, A20's Y4,
   now BY). µ(BY) = 6.
2. **Fiber-pinned floors** (7 queries, 0.4 s total, all UNSAT): every
   nontrivial cover X-logical v with p_y(v) = b has
   |v| = |b| + 2·(off-fiber occupancy) ≥ 12 — m_req = 3 over the hexagon,
   1 over each w = 10 class.
3. **b = 0 stratum, analytic**: p_y(v) = 0 chain-level ⟺ v = τ(u) with u a
   BY cycle and |v| = 2|u|; |v| ≤ 10 would need a nonzero BY cycle of
   weight ≤ 5 — dead by the A16 class certificate (µ ≥ 6). No SAT needed.
4. **Assembly**: p_y*[v] ≠ 0 ⇒ |v| ≥ |p_y(v)| ≥ d(BY) = 12 (exact, §1);
   p_y*[v] = 0 ⇒ (2) or (3). **Certified d(C) ≥ 12**, entirely through the
   y-deck — solver-grade (CryptoMiniSat UNSATs) + the A16 analytic
   certificate + the d(BY) ladder. With §3's τ-lift: **12 ≤ d(C) ≤ 24, both
   ends certified/verified** — the strongest bounds this code has had.

Two forward-looking bonuses: (a) the b = 0 stratum discharges at the FULL
24 target too — τ(u) nontrivial forces [u] ≠ 0, hence |u| ≥ 12, |v| ≥ 24 —
so after the deep census (§6.3), only the doubly-old sector separates the
program from full 24; (b) the census weight-gap (no w = 8) is now a
three-instance pattern worth a lemma hunt (why do weight-3 mirrored pairs
never produce weight-8 stabilizers?).

## 9. (M)@24 census redesign: the A22 fibering port, session 1 (2026-07-22)

`a19_fibering.py`. Setup: Z₃₀ = Z₅⟨z = x⁶⟩ × Z₆⟨v = x²⁵⟩ (CRT), base
H = Z₆×Z₃ (18 sites), fiber Z₅; F₂[z]/(z⁵−1) ≅ F₂ × GF(16), A22's exact
per-site weight table reused verbatim. Fibered polynomials:
A = z⁴v³ + y + y², B = 1 + v + zv². Results:

1. **Weight calculus verified**: 2,000/2,000 random stabilizers satisfy
   |·| = Σ_sites W(ε,δ) on both blocks.
2. **Operator components — inverted landscape vs A22** (their ε-side was
   invertible, δ-side carried the kernel): here A_ε, B_ε have rank 12/18
   (kernel dim 6 each — the Z₂ radical, as §6.3's caveat predicted), while
   **A_δ and B_δ are BOTH invertible over GF(16)[H]**. Consequence: u_δ is
   a free GF(16)[H] coordinate and v_δ = Φu_δ for the single fixed
   transfer operator Φ := B̃Ã⁻¹ — the entire δ-side of the census reduces
   to Φ's sparse-to-sparse behavior (the A22 m-site-theorem shape), and
   the radical lives only in the small ε-side (18-dim F₂, kernels dim 6:
   enumerable outright).
3. **Census-through-the-lens calibration** (the 7 M12 classes + band-12/14
   partials): light stabilizers are monomial-heavy — per-block site counts
   ≈ block weight (hexagon 3+3; w10s split {4,6}; w14 strata 5–7 δ-sites
   per block) — i.e. mostly cost-1 μ₅-monomial fibers, the rigidity a
   classification tree wants.

SAT-census status at retirement: bands ≤ 12 complete (7 + 42 classes),
band 14 partial (89+), JSONL resumable as fallback
(`data/a19/m12_census_classes.jsonl`). Next session: the Φ-sparsity
classification tree (δ-active-site bounds from W ≥ 1 per site, ε-kernel
enumeration for the radical layer) → analytic census to weight 23.

### §9 session 2 (2026-07-22, same day): structural inputs + census redesign

- **Φ-sparsity spectrum** (`a19_phi_eps.py` (A)): |supp φ| = 15/18; min
  joint sparsity |supp U| + |supp ΦU| = 16 / 11 / **6** at |supp U| = 1/2/3
  — the 3↔3 minimum is the hexagon's δ-shadow (matches §8's census-lens
  data). Pure-δ stabilizers (a ≤ 3) have weight ≥ 2·6 = 12; the a ∈ {4,5,6}
  completion is queued for the tree session.
- **ε-layer census, exhaustive** (`a19_phi_eps.py` (B), 2¹⁸ enumeration):
  pair map (A_ε f, B_ε f) has kernel dim 4 — the ε-pair code is [36,14]
  with min nonzero weight 6; full weight histogram banked through w = 22
  (288 at w6, 720 at w8, 1728 at w10, …). The ε-backbone of the
  classification tree is now data, not conjecture.
- **Census redesign — lex-leader symmetry breaking** (`a19_m_census2.py`):
  replaces per-class orbit-blocking (90 clauses/class, the A20 bottleneck)
  with 89 static lex-leader chains (every model is its orbit's lex-min rep;
  one blocking clause per class). Built-in completeness harness: bands
  2–12 must reproduce the a19_m_census ground truth or the run aborts.
  **Validated: {6:1, 10:6} in 15 s (was 227 s) and band 12's 42/42 in 35 s
  (was 1,351 s) — 39×.** Bands 14–22 launched
  (`data/a19/m24_census_classes.jsonl`).

The analytic classification tree (the full A22-style proof-grade census)
remains the moonshot follow-on; its two structural inputs (Φ-spectrum,
ε-layer) are banked above, and the redesigned SAT census now provides the
ground truth it must reproduce at every depth.

**Census progress (2026-07-22, redesigned engine):** bands 2–20 ALL
COMPLETE — per-band classes {6: 1, 10: 6, 12: 42, 14: 54, 16: 487,
18: 1,649, 20: 6,071} (band 20 alone: 2.9 h; the old engine's projection
for the same depth was multi-day). Band 22 in flight (raised cap; ~20k
classes projected from the ~3.5×/band growth). **(M)@24 floors launched**
over the 8,310 banked classes (`a19_m_floors.py --target 24`, m_req
strata 2…9; the hexagon at m_req = 9 / bound 16 is the expected hard
tail). Also closed same day: **d(BX) = 12 exact** (§1) — both n = 180
bases now solver-certified at 12, making the §4 projection floor
symmetric across decks.
