# A36 — constructed doubling pairs: engineering q = k·d²/n > 13.5 covers by design

**Question.** A35 (sibling session, same day) searched the *stored* corpus
presentations of every exact-d base with `q_base > 6.75` and got a wall of
certified negatives: 42 certify runs, 0 hits, with a sharp anatomy — 9 runs
froze exactly at `d_base` and 11 stalled at `d_base + 2/+4`; **zero** even
reached the 2d − 2 deficit wall. This line attacks the same FOM target
(`q = k·d²/n > 13.5` via a certified doubling rung `[[n,k,d]] → [[2n,k,2d]]`,
which doubles q) from the *construction* side: doubling is
presentation-sensitive (A11 Entry 1 — every stored-form negative on the
Z₆×Z₆ engine frame flipped to d = 12 in an anchorable presentation), so
instead of testing the presentations the corpus happens to store, engineer
the presentation cell so the doubling template holds by design.

**Budget rule (hard, user-set, same as A35).** ≤ 5 min compute per code
(one certify() target), one fallback retry ≤ 15 min when 5 is infeasible,
nothing beyond 15. Cheap screens (rank checks, seam arithmetic, base-side
coset probes ≤ ~10 s) are screens, not certification, and are exempt.
Refusals are ledger verdicts (`REFUSED-BUDGET`), never silent.

**CPU contention note.** The A35 sibling ran certification batteries on
this same 10-core machine throughout this session; certify() calls here use
`threads = 4` (vs the front-end default 8) and wall times should be read
with that handicap in mind.

**Status: session 1, 2026-08-11.**

---

## 0. The design theory (what "by design" means)

A certified literal-lift doubling needs four template conditions
(extensibility doc §3, A11 §1); on odd-weight pairs the discriminating one
is condition 3's **safe floor** (A11 S2: `safe_floor_ok` sufficient-shaped
at 0/465 on the corpus, and the A35 refutation anatomy is 20/20 safe-floor
kills). The safe floor is a *base-side* quantity: for every kernel class
ζ ∈ ker ∂₂ ∖ 0, the seam-carry coset `t(ζ) + rowspace(S)` must have min
weight ≥ 2d (A14 Prop A14.1, A29 seam offsets, A30 coset-BZ).

The construction levers, in decreasing directness:

1. **Seam re-routing (Mode 1 — presentation orbits).** The seam offset
   `t(ζ) = carry(Ã·ζ̃, B̃·ζ̃)` depends on which monomial products cross the
   doubled-axis seam — exactly what doubled-axis translations of A and B
   change. The seam-relevant equivalence set per (code, axis) is
   **doubled-axis translations of A and B independently × diagonal unit
   automorphisms (x,y) → (x^u, y^v) × A↔B swap** (undoubled-axis
   translations are exact cover symmetries and cannot move seams —
   A14 §15). d(base), k, and ker ∂₂ are orbit-invariant; the seams and (R)
   are not. A11 Entry 1 (hit2/3/5 flips) and A14 §15 (bb_108, no rescue)
   are the two known outcomes of sweeping this orbit: the sweep decides
   which kind of code we have.
2. **(R) by construction.** k-preservation ⟺ (R) ⟺ `1+δ ∈ (A,B)` in the
   cover ring (A12). Checked per cell (lazy, on S0 survivors — an S0
   reject is a genuine coset element regardless of (R), A14 §16).
3. **Fresh enumeration (Mode 2, fallback).** If a code's whole orbit has a
   safe-floor ceiling < 2d (bb_108-style code-level obstruction), enumerate
   *new* (A,B) pairs on the same frame: k-gate by rank, screen seams first,
   establish d(base) exactly only for survivors. D1∧D2 difference-set
   structure (A5 Entry 15's discovery-sieve surprise) as an enumeration
   prior where compatible.
4. **Descent twists (Mode 3, second shot).** ~2^(wA+wB) sheet-twisted
   covers per extension class (A10); rescued hit2/hit5 historically.
   Solver-tier certificates only — noted as a tier difference.

Screens are tiered exactly like A14 §15, kills are sound at every tier
(each tier exhibits a genuine coset element below the floor):

- **T0 (free):** S0 raw seam minimum over all 2^κ − 1 kernel classes
  (κ = k/2); reject iff min |t| < 2d. Pure numpy seam arithmetic.
- **T0.5 (lazy):** k-gate (cover rank) on S0 survivors only.
- **T1 (cheap):** exact per-class coset minima by SAT ladder
  (CaDiCaL, cap 2d − 1, parity-stepped) — the A11 s3 probe at the
  right cap. Reject iff any class min < 2d.
- **T2 (the verdict oracle):** `bb_lab.doubling_certify.certify()` on the
  cover spec, budget-capped — certificate-tier CERTIFIED /
  DOUBLING-REFUTED / etc. Claim tier: **certified computational data, NOT
  kernel-checked Lean**; SAT/DRAT artifacts (Modes 3) are solver-tier.

## 1. Targets (inside the measured budget envelope)

A35 measured the generic front-end envelope on this machine:
(d_base ≤ 12 AND n ≲ 126 — but n = 126 at W = 22 already BUDGET-FINAL)
∪ (d_base ≤ 10 AND n ≲ 168). k = 16 refused (2^k cosets); k = 14 heavy.
Construction docket, ranked by q_cover × iteration cost (certify wall
times from the A35 ledger for the same points):

| # | base point | group | q_cover target | certify-refute wall | orbit cells/axis |
|---|---|---|---|---|---|
| T1 | [[126,12,10]] | Z₇×Z₉ | **[[252,12,20]] 19.05** | ~40 s | 3,528 |
| T2 | [[120,8,12]] | Z₆×Z₁₀ | **[[240,8,24]] 19.20** | ~410 s | 576/1600 |
| T2b | [[120,8,12]] | Z₅×Z₁₂ | **[[240,8,24]] 19.20** | ~420 s (y only) | 400/2304 |
| T3 | [[90,8,10]] | Z₁₅×Z₃ | [[180,8,20]] 17.78 | ~4 s | 7,200/144 |
| T4 | [[98,6,12]] | Z₇×Z₇ | [[196,6,24]] 17.63 | ~180 s | 3,528 |
| T5 | [[112,6,12]] | Z₇×Z₈ | [[224,6,24]] 15.43 | ~280 s | mid |
| T6 | [[84,6,10]] | Z₆×Z₇ | [[168,6,20]] 14.29 | ~4 s | small |

Excluded as already closed or out of envelope: bb_108 [[108,8,10]] Z₉×Z₆
(A14 §15/§16: presentation orbit + Z₁₈×Z₃ re-decomposition exhausted, no
SF-passing cell — code-level obstruction); [[144,8,12]], [[126,6,12]],
[[150,8,12]], [[162,8,12]], [[168,8/12/14,·]] d = 12 points (census
W = 22 at n ≥ 144, or BUDGET-FINAL measured by A35); [[180,16,10]]
(k = 16 refused); [[126,12,12–14]] Lane-D dreams (n = 126 at W = 22 is
past the envelope — theory route only). bb_90 [[90,8,10]] Z₁₅×Z₃ stored
forms are A14 §14-refuted but its presentation orbit was never swept —
T3 does it.

## 2. Controls (falsify-first)

Before any screen verdict is trusted:

- **Known-positive control:** the A30-certified [[180,4,10]] Z₁₅×Z₆ code
  (`A = 1+y+x`, `B = y⁴+x+x¹¹y²`, x-axis) must pass T0/T1 at floor 20.
- **Known-negative control:** bb_108 stored-y (exact d_safe = 14 per
  A17 §8) must be killed by T1 (and ideally T0).
- **A35 consistency:** the stored cells refuted by the sibling at my
  target points must be killed by T0/T1.

## 3. Results ledger

Every certify()/T1.5 verdict lands here. Tier: certify() and T1.5 rows
are **certified computational data** (complete-enumeration BZ with node
count invariants) — NOT kernel-checked Lean; T0/T1 kills are sound coset
elements (screen tier); T1 passes are enrichment only.

### 3.1 certify() calls

| tag | cover | budget | wall | verdict | detail |
|---|---|---|---|---|---|
| T1p0x_c2_s052 | Z₁₄×Z₉ `y²+x⁶+x⁶y⁴` / `y⁶+x⁵y⁴+x⁶` | 300 | 95.3 s | **DOUBLING-REFUTED** | safe-class coset of weight **14** < 20 (T1 UNDET false-pass, exposed) |

### 3.2 Screen-level closures (sound kills, no certify needed)

| point | axis | orbit size | outcome |
|---|---|---|---|
| T3 [[90,8,10]] Z₁₅×Z₃ (bb_90) | y | 72–144 cells/pres ×3 pres | **entire orbit raw-seam-frozen: S0 = 10 = d_base on every cell** (S0 kill = genuine weight-10 coset element; matches A14 §14 stored-y) |
| T3 [[90,8,10]] Z₁₅×Z₃ (bb_90) | x | 1,800 cells/pres ×3 pres (identical S0 histograms — one code) | 12 k-gate fails/pres; **all 24 top-stratum T1 probes freeze-kill at 10 in ~0.02 s each**; finalists 0 (strata below S0 = 30 unprobed — see §4.4) |

**bb_90 verdict: the v1 presentation orbit (doubled-axis translations ×
diagonal units × swap) is dead on both axes** — the orbit-sweep question
A14 §14 left open for bb_90 now has the bb_108 answer. The freeze is at
d_base exactly, not the 2d−2 wall.

## 4. Session log

### 4.1 Screen validation (controls, all green)

`scripts/a36_orbit_screen.py controls`:

| control | expectation | outcome |
|---|---|---|
| A30-certified [[180,4,10]] Z₁₅×Z₆:x (known positive) | T0/T1 must not kill | S0 = 48 ≥ 20; T1 pass (1 orbit class, UNDET at budget — n = 360 UNSAT is certify()'s job) |
| bb_108 stored-y (known negative, exact d_safe = 14 per A17 §8) | T1 kill + exact pin | S0 = 18 (< 20, T0 already kills); T1 kill, first-found 18; pin: SAT@14 / UNDET@12 ✓ |
| [[126,12,10]] Z₇×Z₉ stored p0:x (A35 DOUBLING-REFUTED) | screen must kill | S0 = 28 (survives T0!); T1 kill at first-found 16–18 |
| seam vectors | must equal `bb_lab.fibering.seam_offsets` exactly | asserted per (code, axis), green |

Two design lessons the controls caught, both now in the tool:

1. **UNSAT is the expensive direction.** A naive exact T1 (prove every
   class ≥ 2d by SAT) hangs on the known-POSITIVE control — the A30
   code's floors took 378–499 s via the compiled BZ kernel, and SAT is
   worse. T1 is therefore a conflict-budgeted KILL screen (kills are
   certificates; passes are enrichment), and certify() is the decider.
2. **Single-orbit kernels need diversification.** Z₇×Z₉ k = 12 codes
   have all 63 nonzero kernel classes in ONE G-orbit; one UNDET on the
   single rep masked a weight-16 kill that the all-classes probe found.
   T1 now probes several *translated members* per class (isomorphic
   instances — one SAT kills, one UNSAT decides the class, extras are
   solver diversification) plus a freeze-probe at cap d first (the A35
   anatomy says most cells die AT d_base).

### 4.2 The swap-twin lemma (found by the first sweep's own data)

`BB(B, A)` is `BB(A, B)` with the two qubit blocks swapped — the SAME
code, and the literal lift commutes with the block swap, so **the A↔B
swap orbit direction produces the identical cover code**. In the first
T1 sweep several swap-twin pairs got (pass, kill@18) — the pass side is
an UNDET false-pass. Consequences, both implemented in the runner:
kills cross-apply to the swap twin, and certify targets are deduped by
swap class. (The swap stays in the sweep — probing both orientations is
free diversification that caught exactly these false-passes.)

### 4.3 T1 point sweeps ([[126,12,10]] Z₇×Z₉, floor 20)

Stored-presentation baseline: S0 = 28, T1-killed (§4.1 control 3), A35
certify DOUBLING-REFUTED. The x-axis orbit of pres 0 (3,528 unique
cells): S0 histogram {10: 720, 24: 360, 26: 288, 28: 624, 30: 432,
32: 216, 34: 336, 36: 96, 38: 144, 42: 264, **52: 48**} — the orbit has
a 48-cell top stratum at S0 = 52, nearly twice the stored cell's 28,
with x-supports concentrated at the top of the half-window (wrap-heavy
"anchored" shape). T1 on the top 24: mix of SF-PASS and kill@18 cells;
swap-consistent double-passes go to certify().

*(results ledger §3 accumulates the certify verdicts)*

### 4.4 T1 point closed by exact T1.5 gates; T2/T2b/T4 screen-dead

The full T1.5 (exact BZ safe-floor) gate over every swap-consistent
T1-finalist of the [[126,12,10]] point, both axes: **27 exact
refutations, 0 certifications** — x: {14: 8, 10: 2},
y: {16: 10, 14: 4, 12: 3}. Combined with the sweeps' sound kills, all
three docket presentations are dead on both axes at the probed strata,
with exact safe minima ∈ {10, 12, 14, 16} — never even reaching the
2d − 2 = 18 wall value. T2 (Z₆×Z₁₀) and T2b (Z₅×Z₁₂) x/y sweeps at
floor 24: 0 finalists (kills 20–22 = at/below the wall); T4 (Z₇×Z₇)
both axes: 0 finalists (kills 18–22).

### 4.5 The pivot: light-class censuses (safe floors are CODE data)

The uniform refutation weights forced the right abstraction into view:
a seam-coset minimum is the min weight of a *homology class* of
H₁(base) — a code invariant — and a presentation cell passes the safe
floor iff its transfer image im Δ (the span of its κ = k/2 basis-seam
classes) avoids every light class (class-min ≤ 2d − 2; all class minima
are even, A17 parity). So ONE complete BZ census over all 2^k − 1
logical cosets per code (labels = the symplectic pairing
φ(v) = (v·zrep_j)ⱼ, presentation-canonical) turns the whole orbit sweep
into exact linear algebra (~ms/cell, no SAT, no UNDET):
`scripts/a36_light_census.py` — with the presentation transport
(translations/units/swap as coordinate isomorphisms) validated
empirically against per-cell ground-truth censuses before any scan
(`verify` mode fixes the translation sign by hit-multiset equality).

**Epistemic caveat (recorded for all wall statistics here and in
A35):** a safe-floor refutation refutes the CERTIFICATE ROUTE, not
necessarily the cover's true distance — the A9 rows-112–152 /
A11 41-row overlap-rescued doubles have SF-false yet d(cover) = 2d,
invisible to the whole A30 certificate machinery. This hunt targets
*certifiable* doubles; "dead" always means certificate-dead.

### 4.6 Supply widening + code identities

- S0 histograms identify code classes cheaply: the T1 docket's pres 1
  and pres 2 share identical histograms on both axes (one code), pres 0
  differs — the point had TWO distinct codes, both now T1.5-dead.
- Supply mined for the census line (`data/a36/triage_input.jsonl`,
  501 rows): merit CSVs + store rows across
  Z₆×Z₇ (252) / Z₇×Z₈ (94) / Z₉×Z₆ (72 — bb_108's parameter point, but
  possibly other CODES than bb_108) / Z₇×Z₉ (35) / Z₆×Z₁₀ (18) /
  Z₅×Z₁₂ (14) / Z₇×Z₇ (7) / Z₁₅×Z₃ (store rows).
- Sibling A35 Lane-D observation (their lane, logged here for docket
  awareness): [[126,12,≥14]] floors landing at q_floor = 18.67 on the
  w = 5 T4.2 family — |A| = |B| = 5 odd, but d ≥ 14 exceeds the
  front-end's DBASE_CAP, so those are not doubling bases for this line.

### 4.6b Orbit dimensions not covered by v1 (queued extensions)

The v1 equivalence set (A14 §15) uses only DIAGONAL unit automorphisms.
Mixing automorphisms exist wherever the group has repeated prime parts
and are absent from every sweep in program history: for Z₇×Z₇,
Aut = GL₂(F₇) has order 2,016 vs the 36 diagonal units — a 56× larger
orbit. The census line prices a mixing-extended scan at ~ms/cell (same
one census per code), so it is cheap to add IF a diagonal scan shows
near-misses. Z₉×Z₆ / Z₆×Z₁₀ / Z₁₅×Z₃ have smaller mixing families
(shared Z₃/Z₂ parts). Descent twists (Mode 3) remain the second
uncovered dimension — the twist changes the extension class and hence
which subspaces of H₁ are reachable; the light census directly
identifies which codes have abstractly-avoidable light sets (the
necessary condition for ANY cover of ANY class to pass).

### 4.7 Ops lessons (recorded for reuse)

- `run_window` .mat files are tag-named: concurrent censuses MUST use
  per-instance workdirs (fixed; the original shared dir was a
  collision).
- Wall-killing a python runner does NOT kill its `cosetbz` C
  subprocess — a killed control run left a 4-thread zombie grinding a
  ~10¹² node window (load 48 on 10 cores) until manually killed. Kill
  the process GROUP or `pkill -f cosetbz` after any wall kill.
- `tail`-piped background jobs buffer: write logs to files and poll
  the files, not the pipe.
