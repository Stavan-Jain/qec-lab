# Research status and assessment

Audit date: 2026-10-07. This is a current reading guide, not a replacement
for the dated research logs. The strongest coherent result is a method for
proving and certifying BB distances through smaller quotient codes. Its
most useful next step is a reusable descent engine with independently
checkable certificates.

The audit read the canonical handoffs, later corrections, theorem notes,
benchmark records, and relevant formal source in the sibling QECLean
checkouts. It inspected previous compilation and axiom-audit records; it
**did not rebuild QECLean**. The fractal seed checker was rerun through
`s=4`, reproducing classical distances `2, 8, 26, 80`. Larger enumerations
and reported timings below are historical results, not fresh measurements.

## How to read the evidence

- **Formal theorem:** a concrete Lean endpoint and its recorded validation.
  State the branch and trust assumptions when they differ; an old progress
  note is not a fresh build result.
- **Analytic theorem:** an argument for all parameters under explicit
  hypotheses. This audit does not turn a paper proof into a Lean theorem.
- **Computational certificate:** complete finite enumeration or checked
  proof artifacts plus a mathematical reduction. Python/C execution remains
  outside the Lean kernel unless a separate formal endpoint is cited.
- **Evidence or conjecture:** finite tests or plausible extrapolations.
  A logical witness alone establishes an upper bound on distance.

## 1. Even-cover distance theory and the gross code

**Result.** The program proves the gross code's exact distance
`[[144,12,12]]` and abstracts its argument into a sufficient distance-doubling
theorem. The cover must satisfy explicit homological, weight-floor, and
attaining-witness conditions. The structural equivalence

`deck action trivial on H1 ⇔ k(cover)=k(base) ⇔ 1+deck ∈ (A,B)`

turns a homological hypothesis into an algebraic certificate.

**Evidence.** Gross has an unconditional formal distance endpoint in
QECLean `Gross/Distance.lean`, in addition to the analytic argument in
[A4, Theorem D](../experiments/bb_lab/notes/A4_writeup.md#14-theorem-d-assembly-and-status).
The equivalence and its proof are in
[A12 §3](../experiments/bb_lab/notes/A12_deck_homotopy_R.md)
(lines 99–141). The [extensibility guide](gross-distance-extensibility.md)
explains the reusable layer and instance obligations.

**Why researchers should care.** Cover constructions can cheaply generate
many larger codes, but their distance is the expensive missing parameter.
This theorem identifies conditions that actually force doubling and makes
the obstructions inspectable on a smaller code. It addresses the even-cover
case left conjectural in the distance bounds of
[Symons–Rajput–Browne, v3](https://arxiv.org/abs/2511.13560).
The value is the reusable sufficient theorem and its verified instance;
gross's numerical distance and its cover interpretation were already known.

**Limits.** This is not a universal doubling law. The
[descent-cover screen](../experiments/bb_lab/notes/A10_descent_twist_screen.md)
contains bases whose entire enumerated descent space fails to double.
Presentation choices and twists matter. The
[small-cycle theorem](../experiments/bb_lab/notes/A16_class_theorem_writeup.md)
supplies a class-wide `d≥6` result, but explicitly consumes a one-sided
annihilator floor and excludes frames divisible by four; its hypotheses
must not be shortened to “all weight-six BB codes.”

## 2. Descent certificates that reach difficult finite codes

**Result and evidence.** The tower slice calculus decomposes a low-weight
cover cycle into a smaller shadow plus bounded overlap. Complete censuses
and lift checks certify, without SAT on these critical paths:

| Code | Exact distance | Recorded computation | Source |
|---|---:|---:|---|
| Bravyi `[[360,12,d]]` | 24 | about 9.4 min | [A32](../experiments/bb_lab/notes/A32_tower_slice_calculus.md), lines 1–22 |
| IBM class Y `[[288,8,d]]` | 20 | about 105 s | [A33](../experiments/bb_lab/notes/A33_ibm288_d20_certificate.md), lines 1–34 |
| Tour-de-gross `[[432,12,d]]` | 24 | about 20 min | [A40 §3.2](https://github.com/Stavan-Jain/qec-lab/blob/50fef445515748f2f291474a5a1c838f8070fef0/experiments/bb_lab/notes/A40_tour_de_gross_column.md), lines 285–305 |

These are computational certificates, not claims that those three full
results have been checked by Lean. A32 closes the `≤24` entry discussed in
the [original BB paper](https://arxiv.org/abs/2308.07915); A33 independently
reproduces a previously solver-exact value. A40 checks one further member
of the original IBM family, not its all-size distance conjecture.

**Why researchers should care.** A search that produces only a low-weight
logical cannot rule out a still smaller one. These methods provide that
missing lower-bound evidence and can turn a candidate into a code with
known parameters. The method has also been tested beyond its motivating
instances: [A39](../experiments/bb_lab/notes/A39_descent_theory_validation.md)
records 97 new exact distances and 17 additional floors.

The [frozen-cohort scorecard](../experiments/bb_lab/data/descent_theory_test/PHASE4_SCORECARD.md)
(lines 211–234) reports 128/154 rows resolved by descent versus 47/154 by a
CryptoMiniSat ladder, with no SAT-only resolution. **Resolved means a floor
or an exact distance** for the descent arm; the 47 SAT outcomes were exact.
Seventeen independent ground-truth solves agreed, with three timeouts.

**Limits.** This is a restricted, structured cohort and a specific SAT
baseline, not a universal speed comparison. The timing columns use descent
CPU and uninterrupted SAT wall time; SAT received a 60-second minimum
budget and did not include the tandem MaxSAT stage. The preregistered
“GREEN reaches its predicted floor” criterion failed: 54/138, not 90%.
Seventy-four cases had smaller exact distances; ten hit a memory limit.
Those failures are useful design constraints, not results to suppress.

The current bundle's provenance check matches 15/16 original manifest
entries; the remaining source hash is the documented amendment. The
phase-two manifest matches 16/18 entries: `watch_sweep.sh` and
`analysis_final.json` are absent. Main result tables, the scorecard, and
compute drivers match their hashes. Restore or document these ancillary
omissions before claiming a fully reproduced release bundle.

## 3. A formal weight-six family with unbounded actual k d²/n

**Result.** For `q=2^s`, `N=q²−1`, the separable BB checks
`A=1+x+x^q`, `B=1+y+y^q` give

`[[2(4^s−1)², 2·4^s, d_s]]`, with `3^s ≤ 2s d_s`.

Here `d_s` is an attained classical minimum, proved equal to the quantum
distance. A selected subsequence has strictly increasing, unbounded actual
`k d_s²/n`, with weight-six stabilizer generators.

**Evidence.** The formal development is in QECLean branch
`codex/bb-growing-fom`, inspected at commit `1e7f444`. The
[theorem note §7](https://github.com/Stavan-Jain/qec-lab/blob/12156e6f463d5d226a373c2c8ac005202f7b2f59/experiments/bb_lab/notes/BB_FRACTAL_UNBOUNDED_FOM_THEOREM.md#7-exact-finite-values-and-the-formal-endpoint)
links the precise declarations; the
[validation record](https://github.com/Stavan-Jain/qec-lab/blob/12156e6f463d5d226a373c2c8ac005202f7b2f59/pipeline/attempts/bb_growing_fom/validation.md)
reports compilation and only `propext`, `Classical.choice`, and `Quot.sound`
in the axiom audit. This is branch work, not a claim about its merge status.

**Why researchers should care.** It provides an explicit sparse family
whose claimed asymptotic figure of merit follows from a physical-weight
lower bound, rather than extrapolated finite distances. The useful proof
step handles cyclic boundary cancellation by averaging open fractal
histories over translations; that is a concrete technique to reuse.

**Limits.** The seeds and product construction connect to
[Yoshida's fractal codes](https://arxiv.org/abs/1111.3275) and
[hypergraph products](https://arxiv.org/abs/1202.0928); no novelty over those
constructions is established here. `d_s=3^s−1` remains conjectural.
The first four exact checks reach `[[130050,512,80]]` at `s=4`, so this is
not evidence of a small hardware implementation. The offsets grow with
`s`; fixed check weight does not mean fixed-range two-dimensional locality.

## 4. The non-abelian mitten code reaches a formal distance endpoint

**Result.** The `C5×S3` lifted-product code has a packaged
`StabilizerCodeWithDistance 150 30 10`, with both CSS lower bounds and an
attaining logical. Its proof is not inherited from a BB transpose symmetry.

**Evidence.** [Attempt state](../pipeline/attempts/mitten_150_30_10/state.yaml),
[build ledger](../pipeline/attempts/mitten_150_30_10/build_budget.md), and
[independent solver verification](../experiments/bb_lab/notes/mitten150_tandem_verification.md)
record the formal and computational routes. The inspected
`claude/a32-lifted-product` branch still uses `native_decide`; this is not
a current QECLean `main` theorem at its three-axiom standard. The completion
record reported 67.5 seconds of wall time for the 14-module rebuild; that
is a historical measurement, not a fresh build or a timing guarantee.

**Why researchers should care.** This checks whether the verification
architecture survives noncommuting group algebra and asymmetric CSS
sectors. It does so on a non-abelian, rate-one-fifth code, supplying a useful
test case for other lifted-product constructions. It verifies a known
distance; it is not a new code or an improved distance claim.

## 5. An all-size correction and sector classification for original IBM

**Result.** For the original Laurent pair on periods `ell,m` divisible by
six, the analytic module calculation gives

`k = 12 + 14·[127 divides both ell and m]`.

Thus the `b=1` tour-de-gross column has `k=12`, while the square `b=0`
column has `k=26` when `127|r`. The same argument identifies every compact
cylinder-image subspace and the index-two pullback kernels.

**Evidence and significance.** The self-contained
[parity/Tor proof](https://github.com/Stavan-Jain/qec-lab/blob/50fef445515748f2f291474a5a1c838f8070fef0/experiments/bb_lab/notes/A40_PARITY_TOR_PROOF.md)
(lines 21–60 and 143–151) is an analytic theorem. It exposes a parameter
exception that small-size tables would miss and tells a distance proof
exactly which logical sectors it must cover: normally 189 nonzero
cylinder-liftable classes and 3,906 mixed classes. This makes it useful to
researchers studying the original family even before its metric is solved.

**Limits.** The homological coordinates are not Hamming isometries. This
does not prove an IBM distance floor. The
[current program](https://github.com/Stavan-Jain/qec-lab/blob/50fef445515748f2f291474a5a1c838f8070fef0/experiments/bb_lab/notes/IBM_EXACT_DISTANCE_PROGRAM.md)
and [corrected state](https://github.com/Stavan-Jain/qec-lab/blob/50fef445515748f2f291474a5a1c838f8070fef0/experiments/bb_lab/notes/TOUR_DE_GROSS_STATE.md)
record the failed crossing lemma and other invalid metric shortcuts.

## 6. Exact modified families with physical contraction proofs

**Result.** Changing one original IBM monomial gives the analytic family
`[[98r²,12,4r]]` when `3∤r`, and `[[98r²,20,10r/3]]` when `3|r`, with
weight-six checks. A related primitive-trinomial theorem produces
`[[2(2^t−1)²r²,4t,2^(t−1)r]]`, conditional on a specified primitive
trinomial of degree `t`.

**Evidence.** [The one-monomial theorem](https://github.com/Stavan-Jain/qec-lab/blob/50fef445515748f2f291474a5a1c838f8070fef0/experiments/bb_lab/notes/IBM_ONE_MONOMIAL_K12_FAMILY.md)
(lines 8–48) and [two-contraction theorem](https://github.com/Stavan-Jain/qec-lab/blob/50fef445515748f2f291474a5a1c838f8070fef0/experiments/bb_lab/notes/IBM_TWO_CONTRACTIONS_EXACT_THEOREM.md)
(lines 7–57) give analytic proofs, not Lean certification. A
[weight-eight modification](https://github.com/Stavan-Jain/qec-lab/blob/50fef445515748f2f291474a5a1c838f8070fef0/experiments/bb_lab/notes/IBM_MODIFIED_EXACT_DISTANCE_THEOREM.md)
also has exact `[[72r²,12,4r]]` parameters.

**Why researchers should care.** The argument separates two tasks that are
often conflated: algebraic normal forms detect every logical class, while
physical contraction maps control its weight. Together with an empty-seam
lift and an attaining stripe, this gives exact distance for all repetition
parameters, including repeated-root periods. These are useful controlled
examples for developing methods, rather than claims of superior parameters.
They change the original checks; none solves the unchanged IBM conjecture
or establishes novelty or inequivalence to known product constructions.

## Recommended next research direction

Prioritize **a reusable descent certificate format and small independent
checker**, supported by the existing tower engine. It has the strongest
combination of demonstrated reach, measurable benefit, and bounded next
steps. A new code search should consume this machinery, not replace its
validation with more upper-bound tables.

1. Specify the complete low-weight census, lift-fiber coverage, quotient
   transport, and nontrivial attaining witness as separate certificate
   obligations. Make each failure return UNKNOWN or a verified witness.
2. Reproduce one mixed-axis and one same-axis flagship certificate with a
   checker independent of the discovery code; then formalize a reusable
   checker theorem in QECLean. This is proposed work, not an existing result.
3. Repair the documented empty-preimage edge and price candidate memory as
   well as enumeration nodes. Recheck current source before assuming a
   historical “owed” task remains undone.
4. Freeze a fresh cohort and compare against both CMS and tandem MaxSAT
   under matched accounting. Report floors, exact distances, timeouts,
   memory stops, and previously known anchors separately.

The fractal exact-distance formula is a worthwhile secondary theory
problem. The original IBM all-size distance law is a higher-risk program:
several tempting contraction, mirror, and crossing shortcuts have already
been disproved.

The next soundness theorem can be stated precisely: for a free two-sheet
cover, a complete census of relevant base shadows plus complete carry-fiber
checks that classify every lift of weight below `D` as a stabilizer implies
`d(cover) >= D`. The identity `|v| = |shadow| + 2*overflow` bounds the fibers
that must be checked. Start with a self-contained A33 replay, then A32;
require the verifier to reject omitted or duplicated census records. Exact
node counts and hashes are useful checks on a producer, but do not alone
prove to an independent checker that every combination was visited.

### New investigation in this audit

The investigation found and fixed two fail-closed edges in
[`cosetbz.py`](../experiments/bb_lab/src/bb_lab/cosetbz.py): zero-radius
searches incorrectly visited singleton subsets, and rank-zero codes
incorrectly rejected empty information sets. These caused assertions,
not false certificates. The new regression suite passed **3 tests**;
independent checks matched exhaustive results in **100 small coset cases**
and **1,120 bounded carry-fiber comparisons**. These checks do not constitute
a fresh rerun of the large flagship certificates.

The research subagent also replayed the two documented counterexamples:
the [A43 rectangle](https://github.com/Stavan-Jain/qec-lab/blob/50fef445515748f2f291474a5a1c838f8070fef0/experiments/bb_lab/notes/A43_S2_codex_audit.md) has
weight 48 on `(54,6)`, below the proposed bound 54; the
[A42 S19 quotient](https://github.com/Stavan-Jain/qec-lab/blob/50fef445515748f2f291474a5a1c838f8070fef0/experiments/bb_lab/notes/A42_S19_transverse_quotient.md)
has weight 6 below the proposed floor 12. These disprove the indicated
proof shortcuts, not the original IBM distance conjecture.

## Mathematical scope checklist

- [x] Gross exact distance and a sufficient doubling framework.
- [x] Specific tower certificates and the recorded A39 cohort outcomes.
- [x] Fractal uniform lower bound, attained distance, and formal unbounded
  subsequence; finite exact checks through `s=5` are separate evidence.
- [ ] Fractal `d_s=3^s−1` for all `s`. The formal growth subsequence is
  `s=16·4^t`; the sharper analytic `s=16·2^t` also uses an unformalized
  upper witness. Do not conflate them.
- [ ] Original IBM `d=6(2r+b−1)` for every member. Its exceptional `k`
  values must be included in any eventual family statement.
- [ ] A universal doubling or same-axis iteration law. The repo contains
  counterexamples; a narrower theorem needs explicit hypotheses.
- [ ] Lean verification of every Python/C tower certificate, a universal
  runtime guarantee, or a decoding/circuit-performance theorem.
- [ ] Global novelty or “first ever” claims. The source comparisons here
  establish relevant context, not an exhaustive priority determination.

## Repository and reproduction

Keep theorem sources, finite witnesses, compact certificate summaries,
protocols, manifests, generators, and pinned third-party provenance. Large
regenerable outputs need explicit regeneration or archival instructions.
See [repository maintenance](repository-maintenance.md) for this audit's
cleanup record and the commands for reproducing local state.
