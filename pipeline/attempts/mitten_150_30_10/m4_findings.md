# M4 findings — the floor architecture (design + de-risk session)

**Date:** 2026-08-11 · **Tools:** `experiments/bb_lab/scripts/a32_m4_{scoping,collect,probe_gen}.py`
· **Data:** `experiments/bb_lab/data-a32/m4_triples.json` · **Lean rehearsal:**
scratch probe, QECLean toolchain, one `native_decide`.

## G1 — The triple-equation reduction (the architecture)

Writing a `ker H_X` vector in blocks `v = (v₀, v₁, v₂, v₃, v₄)` (grid
`m = 2α+β`, shared `m = 4`), the two check blocks are, in `F₂[G]`
(left-mults for A, right for B; `G = C₅×S₃`):

```
(E0)  a₀·v₀ + a₁·v₂ = v₄·b₀          (E1)  a₀·v₁ + a₁·v₃ = v₄·b₁
```

— coupled ONLY through `v₄`.  Mirror for `ker H_Z` (antipode forms,
right-mults `b₀~ , b₁~` on the pair, left `a_α~` on the shared block).
Hence the floor reduces to:

1. **classify light triples**: all `(u, w, t)` with the instance
   equation and `|u|+|w|+|t| ≤ 9` — four instances: X0, X1 (b₀/b₁ with
   (a₀,a₁)) and Z0, Z1 (a₀~/a₁~ with (b₀~,b₁~));
2. **join** the two instances of a side along the shared `t` with total
   weight `(pair₀) + (pair₁) + |t| ≤ 9`.

Determinacy (`rank L(a₁) = R(b₀) = R(b₁) = 30`, antipode mirrors; only
`L(a₀)/L(a₀~)` are singular, rank 28): any two of `(u,w,t)` determine
the third — uniquely except through `a₀` (a 4-element `Ann(a₀)` coset).
The sweeps therefore enumerate the two lightest components and DERIVE
the third; nothing exponential in the heavy component is ever touched.

In-triple ε-parity (augmentation hom, all base sets odd):
`|u|+|w| ≡ |t| (mod 2)` — triple weight is always even; the split table
is all `(p,q,r)` with `p+q+r ∈ {0,2,4,6,8}` (95 splits), and the
would-be-worst `(3,3,3)` split does not exist.

`Ann(a₀)`, `Ann(a₀~)`: all three nonzero elements have weight **20** —
pure-annihilator configurations are dead far above the cap.

## G2 — Falsify-first catch: the families are NOT the whole story

Expected classification "= the 60 row fragments per instance" is
**false**: the probe simulation fired at split (1,2,3) (15 solutions =
3 C₅-orbits off-family).  The complete exhaustive censuses (vectorized
split sweeps; `t` derived, never enumerated; parity kills `p+q = 9`):

| instance | light triples (≤ 9) | of which row fragments |
|---|---:|---:|
| X0 | 1,996 (incl. 0) | 60 |
| X1 | 3,531 | 60 |
| Z0 | 2,326 | 60 |
| Z1 | 3,251 | 60 |

The stray triples are harmless BY THE JOIN: **joining the two lists of
each side along `t` (weights re-summed, ≤ 9) yields EXACTLY the 60
generator rows of that side — machine-verified for both sides.**  This
is the entire M4 floor statement, checked end-to-end offline
(`a32_m4_collect.py`, 67 s; hard asserts).  Consequence for Lean: the
classification lists are emitted data (~2-3.5k packed 90-bit Nats per
instance), and the join is one more finite `native_decide` per side
(~12M pair checks).

Independent cross-check: the blocking-clause SAT enumeration of the
same four censuses (`a32_m4_scoping.py`) was still running at session
close (thousands of solver rounds); when it lands its lists must equal
`m4_triples.json` — the falsify-first two-tool agreement gate (G5).

## G3 — Lean budget rehearsal: 7.3 s per instance, one `native_decide`

Scratch probe (import-free, QECLean toolchain v4.30): the FULL X0 sweep
— all 95 splits, 8,569,930 swept classes, three derivation modes,
classification checked against the full 1,995-triple list — as ONE
`native_decide`:

| driver version | wall |
|---|---:|
| naive (per-pair `xorFold`, fueled popcount) | 69 s |
| fold-hoisted (`foldedPairs` per element) + SWAR popcount | **7.3 s** |

The 9.4× lesson: hoist the table-fold to the per-ELEMENT lists
(`List (Nat × Nat)` of (mask, folded)); the quadratic body must be xor +
SWAR-popcount + compare only.  `masksOfWt` recursion + 2k-entry packed
list literals need `set_option maxRecDepth 65536` (list literals are
nested `cons`).

Projection: 4 instances ≈ 25-30 s + 2 join natives (~2-4 s) + data
elaboration + chain plumbing → **M4 ≈ 40-70 s against the 180 s
allocation**; total certification trending ≈ 80-110 s of the 300 s cap.
C₅ transport is NOT needed (skipping its equivariance machinery
entirely); revisit only if the real files overshoot.

## G4 — Lean formalization plan (next sessions)

Files under `QEC/Stabilizer/Codes/Mitten/M150/` (leaf-file-per-side so
lake parallelism makes the two sides' sweeps concurrent):

1. `FloorCore.lean` — generic mask layer: `masksOfWt` + completeness
   (`m < 2^30 ∧ popCnt m = k → m ∈ masksOfWt 30 k`), `xorFold`-testBit
   linearity, SWAR-popcount = card-filter bridge, `foldedPairs`, the
   three-mode sweep driver + ONE generic soundness lemma per mode
   (consumes the driver `Bool` via `List.all` completeness; no per-split
   dispatch — the split list is data).
2. `FloorData.lean` (GENERATED; emitter mode `m4` in
   `m150_gen_lean_data.py`) — per instance: derivation-map columns
   (`A₁⁻¹A₀`, `A₁⁻¹R_b`, `R_b⁻¹A₀`, `R_b⁻¹A₁`, pseudo-inverse `P`,
   left-null filters, `Ann` lists — mirrors for Z), split/mode lists,
   packed classification lists, join row tables.  All numpy-validated
   before emission per GENERATORS.md.
3. `FloorZSide.lean` / `FloorXSide.lean` — chain plumbing per side:
   block split of `chainWeight`, `dualBfn`/`lpBoundary1Fn` → triple-form
   bridges (basis expansion + pointwise `native_decide`, the M2 §2
   pattern), chain↔mask bridges (`maskOf` testBit/weight lemmas), the
   per-instance classification theorems, the join, and
   `floorZ : ∀ c, dualBoundary c = 0 → chainWeight c ≤ 9 → c ∈ dualBoundaries`
   (mirror `floorX` on cycles/boundaries).
4. `Distance.lean` (M5) — `chainWeight_lower_bound_transfers` at K = 10
   + the M3 witness → `HasCodeDistance m150StabilizerCode 10` →
   `StabilizerCodeWithDistance 150 30 10`.

Soundness-lemma inventory (all patterned on existing repo proofs):
mode-ut/uw need table-defines-map pointwise natives + `A₁⁻¹A₁ = id`-type
inverse certs; mode-wt needs the gross-§2-style corrected decoder for
the singular `a₀` (dropSet = the 2 free columns) — direct port of
`decoder_identity_X` + `face_kernel_trivial`-with-drops.

## G5 — SAT cross-check: CONFIRMED (landed 2026-08-11, after G6)

The blocking-clause SAT enumerations exhausted (final UNSAT) on all four
systems: X0/X1/Z0/Z1 = 1995/3530/2325/3250 nonzero solutions in
1204/564/1578/614 s.  Counts equal the collector's (minus the zero
triple), `n_missing = 0` (all 60 family members found per instance),
and the per-split histograms are EQUAL dict-for-dict between
`m4_scoping.json` and `m4_triples.json` on every instance.  **Two
independent tools, one answer** — the classification lists the Lean
sweeps certify against are cross-confirmed.  (Original gate text kept
below for the record.)

### (original gate note, superseded)

`a32_m4_scoping.py` (blocking-clause SAT over the same four systems)
was launched this session and left running (>60 min CPU at write time;
thousands of blocking rounds).  **Gate re-scoping, with reasoning:**
the original note said to hold `FloorData.lean` for this comparison.
That was over-cautious: the Lean sweeps are *self-certifying* against
the classification lists — a missing entry makes the corresponding
`checkAll` fire off-list and `native_decide` fails the build (this is
exactly what the probe demonstrated when run against the fam-only
list), and a bogus extra entry either fails the join native or is
sound-harmless (the lists are consumed as upper bounds; the join
re-verifies every pair it admits).  So the *in-build* obligations
subsume the census-correctness question, and `FloorData.lean` was
emitted and verified green without waiting.  The SAT comparison remains
worth recording here when the run lands (expected: per-instance
solution sets equal `m4_triples.json`, i.e. `n_extras` = stray counts,
`n_missing = 0`); record it, but nothing downstream is blocked on it.

## G6 — M4 sweep layer GREEN in QECLean (2026-08-11)

`FloorCore.lean` (generic mask/sweep layer, fully symbolic: fueled/table
popcount + agreement, `xorFold` linearity + basis + `linear_ext`,
`masksOfWt` completeness, the 4-mode driver — mode 3 = coset-on-`t`,
added when Z-side t-heavy splits turned out to need it — with one
soundness lemma per mode, packed-triple extraction, `checkJoin`),
`FloorData.lean` (generated, mode m4), `FloorSweep{X,Z}.lean` (six
`native_decide` obligations).  Measured: FloorCore 2.4 s, FloorData
14 s, SweepX 15 s, SweepZ 40 s — the sweeps and joins of BOTH floors
now verify inside the build at 71.4 s total.  The Z side pays a
~30-step identity-`tP` fold per class on its unique-`u` mode-2 splits
(worth optimizing only if M5 squeezes the budget).
