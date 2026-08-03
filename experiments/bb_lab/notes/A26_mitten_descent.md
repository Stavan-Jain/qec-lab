# A26 — Mitten-code descent: cover-detect-and-descend beyond bivariate

**Date:** 2026-08-03 · **Branch:** `claude/priceless-merkle-4060c1` ·
**Status:** experiment complete (ISD-estimate tier); exact-d and Lean tiers = follow-ups.
**Paper:** arXiv:2607.28795 "High-rate qLDPC processors" (mitten codes, Caltech/Oratomic).
**Artifacts:** `scripts/a26_mitten_descent.py`, `scripts/a26_export_groups.g`,
`instances/mitten_groups/*.txt`. Shipped-matrix cross-checks need a checkout of
`github.com/a7b/yarn` (`--yarn PATH`) and `pynauty` (not a bb_lab dependency; install ad hoc).

## §0 TL;DR

First port of the lab's descent move to non-abelian LP codes. All eight published
mitten codes rebuilt from Table XIII + GAP `SmallGroup` element order and validated
(5/8 byte-identical to the shipped artifacts, [[150]]/[[200]] equivalent up to
relabeling — and **[[300,60,14]]'s Table XIII row is an erratum**, see §4). Deck
census: central involutions exist for exactly [[200,40,12]], [[300,60,14]],
[[540,108,18]]; no other translation decks anywhere. All three descend cleanly
(k halves; entries stay weight-3; pushforward intertwining verified exactly):

| cover (published) | base = cover/⟨ι⟩ | base group | d ratio |
|---|---|---|---|
| [[200,40,12]] | [[100,20,**≤8**]] | C₂×D₁₀ | 12/8 (2·d_base−4) |
| ↳ [[100,20,8]] | [[50,10,**≤4**]] (entries degenerate to wt 3,1,1,3) | D₁₀ | 8/4 (**2·d_base**) |
| [[300,60,14]] (shipped) | [[150,30,**≤8**]] — **NOT the published [[150,30,10]]** (nauty-inequivalent) | C₅×S₃ | 14/8 (**2·d_base−2**) |
| [[540,108,18]] | [[270,54,**≤10**]] | (C₉⋊C₃)⋊C₂ | 18/10 (**2·d_base−2**) |

Headline findings: (1) **the tower hypothesis is false** — the published family is
not a descent chain; [[300]]'s central quotient is a strictly weaker (estimated d = 8) *sibling* of
[[150,30,10]] on the same group C₅×S₃. (2) Mitten 2-covers are **k-doubling**
(k(cover)=2k(base) in every case) — the complementary regime to the BB gross tower
(k-preserving, distance-doubling); accordingly descent *loses* distance here.
(3) **Deficit-wall echo**: two of four descents sit exactly at d = 2·d_base−2 — the same
wall value A17 found across unrelated bivariate codes — now appearing in a
non-abelian, k-doubling family. (4) Real-code minimum-weight logicals are all
**dangerous-sector** (pushforward = a genuine base logical); only the erratum
variant's spurious wt-6 logicals are deck-odd.

All ≤w distance entries are ISD upper bounds believed tight: thousands of independent hits at the quoted
weight and nothing lighter in 20k trials/side — *estimates, not floors* (lab rule).
Exact certification belongs to the Tandem lane (the [[150,30,10]] pattern of
qec-lab PR #10 applies verbatim; note the `-cost-step` parity premise fails at odd
check weight 9, use naive OPTIMUM + CMS UNSAT ladder).

## §1 Setup and provenance

- Groups: `a26_export_groups.g` exports multiplication tables in GAP's
  `Elements(G)` order for the eight Table VII/XIII `SmallGroup` IDs, plus
  `DirectProduct(C10,S3)` order (needed for the shipped [[300]], §4).
  GAP 4.15.1 via conda-forge; tables vendored under `instances/mitten_groups/` so the
  experiment reruns without GAP.
- Construction: Eq. (J1) conventions (L(g)b(h)=b(gh), R(g)b(h)=b(hg⁻¹),
  a* = sum of inverses). Verification against `yarn/processor_codes`:
  [[500]], [[540]], [[630]], [[780]], [[975]] **byte-identical**; [[150]], [[200]]
  Tanner-isomorphic (shipped files use unnormalized set representatives, extracted
  in §4); [[300]] see §4. CSS condition, k = |G|, full ranks: all eight ✓.
- Deck census (`translation_deck_census`): left decks need t·aᵢ·t⁻¹ = aᵢ, right
  decks s⁻¹·bᵢ·s = bᵢ, central involutions satisfy both unconditionally. Result:
  the ONLY translation decks on any published mitten code are the central
  involutions of [[200]] (ι=#3 ∈ C₄), [[300]] (ι = the C₂ ⊂ C₁₀), [[540]]
  (ι=#3 = y⁶, the group's unique involution). [[150]]/[[500]]/[[630]]/[[780]]
  have involutions but none survive the conjugation conditions; [[975]] has none
  (odd order). Centers: C₅, C₄, C₁₀, C₅, C₂, C₃, C₁₃, C₅ respectively — the three
  even centers are exactly the descendable codes, as predicted from Table VII.

## §2 Descent mechanics (all verified exactly, not numerically)

For each central ι: the qubit/check permutation L(ι) (same on all 5 data and 4
check blocks) satisfies `H[perm rows, perm cols] == H` exactly and is free; the
quotient code is the mitten code on G/⟨ι⟩ with pushed supports (F₂ cancellation
tracked); and the pushforward intertwining `H_base · p = p_check · H_cover` holds
exactly for both H_X and H_Z. k halves in every descent (60→30, 40→20, 20→10,
108→54) — **mitten towers are k-doubling** (deck acts nontrivially on H₁; the
A12 (R)-condition FAILS by construction since k=|G| is forced by shape). Contrast:
gross tower is k-preserving with d = 2·d_base; here d < 2·d_base in 3 of 4 descents.

Support collisions (two set elements in one ⟨ι⟩-coset → F₂ cancellation) never
occur at the first level; the [[100]]→[[50]] level collapses a₁ and b₀ to the
identity alone (entry weights 3,1,1,3), so the bottom [[50,10,≤4]] leaves strict
mitten form but stays LP with the canonical structure intact (a₁ = e is trivially
full-rank).

## §3 The tower hypothesis: REFUTED

(C₁₀×S₃)/⟨ι⟩ ≅ C₅×S₃ made it natural to ask whether [[300,60,14]] descends onto
the published [[150,30,10]]. It does not: the quotient is a [[150,30]] mitten code
on C₅×S₃ with sets a0=(4,9,19) a1=(15,21,28) b0=(3,10,22) b1=(0,8,15) (in the
quotient labeling), estimated d = 8 (11.9k/12.2k hits per side at 8, nothing lighter in
20k trials/side), and it is **not Tanner-isomorphic** to [[150,30,10]] (nauty
certificates differ; the published code's d=10 is now independently certified by
the PR #10 Tandem run). So the published mitten instances are genuinely separate
discoveries, not one chain — and this weaker sibling shows the C₅×S₃ design space
contains strictly worse codes that descent reaches but their pipeline (rightly)
skipped. Same story one level down: [[270,54,≤10]] is dominated by [[150,30,10]]
(same estimated distance, 1.8× the qubits).

## §4 Erratum: Table XIII's [[300,60,14]] row

Three independent facts show the paper's Table XIII row for [[300,60,14]]
(a0={38,51,54}, a1={0,6,45}, b0={25,33,48}, b1={0,16,58} over SmallGroup(60,11))
is not the code of Tables I/VI/VIII:

1. its distance is **6**, not 14 (wt-6 logicals on both sides, found instantly;
   they are deck-odd/fiber-antisymmetric, so they even vanish under pushforward);
2. R(b₁) has rank 58/60 — violating Definition 4's full-rank requirement
   (square-invertibility survives only via the b₀ column);
3. its canonical-basis weights are {16,18}, not the published 22/22.

The real [[300,60,14]] is the shipped artifact: a mitten code on the same group
in `DirectProduct(C10,S3)` element order with sets a0=(4,9,49), a1=(15,21,58),
b0=(3,10,52), b1=(15,30,38) — matches the published d = 14 (ISD minima at 14) and canonical weights 22/22, and is
Tanner-inequivalent to the Table XIII build. (The two codes' quotients also differ:
ISD minima 3 vs 8.) **Action: report upstream to the yarn/paper authors.** Until then, any
lab work on "[[300,60,14]]" must use the shipped sets above, not Table XIII.

## §5 Sector structure and the deficit-wall echo

Pushforward classification of minimum-weight cover logicals (one witness/side):

| cover | lane | wt | pushes to |
|---|---|---|---|
| [[200,40,12]] | Z / X | 12 / 12 | base logical wt 12 / base logical wt **8 = est. d(base)** |
| [[100,20,8]] | Z / X | 8 / 8 | base logical wt **4 = est. d(base)** / base logical wt 8 |
| [[300,60,14]] (shipped) | Z / X | 14 / 14 | base logical wt 10 / base logical wt 12 |
| [[540,108,18]] | Z / X | 18 / 18 | base logical wt 16 / base logical wt 18 |

Every real-code minimum sits in the **dangerous sector** (image is a genuine base
logical), which only enforces d ≥ d_base — the binding constraint is therefore the
deck-odd (safe) sector, exactly the two-tier structure of A14/A17. And the
observed values 14 = 2·8−2 and 18 = 2·10−2 put two of four descents exactly ON
the deficit wall (A17: 2·d_base−2 as the recurring orbit-ceiling value across unrelated
bivariate codes). One descent ([[100]]→[[50]]) is exact doubling, one ([[200]]→
[[100]]) is 2·d_base−4. n=2 is anecdote, not theorem — but the wall value appearing in
a non-abelian, k-doubling family (both regime knobs flipped vs. the bivariate
towers where it was found) upgrades the deficit wall from "bivariate curiosity"
toward "generic 2-cover phenomenon", and any wall theory should now be tested
against mitten towers too.

## §6 Follow-ups (ranked)

1. **Exact d for the four descent codes** via the PR #10 Tandem pattern
   ([[100,20,8]] and [[150,30,8]]-sibling are cheap; [[270,54,10]] moderate).
   Turns §0's ≤ column into certified values and makes the wall echo exact.
2. **Safe-sector floors for mitten towers**: port the A14 safe-floor SAT method
   (augmented H_X membership rows) off bivariate — the descent data says that's
   where mitten cover distance lives. Cheap prototype: [[200]]/ι with base d=8.
3. **Erratum report** to the authors (Table XIII [[300]] row + the labeling-
   convention footnote that processor_codes for [[150]]/[[200]]/[[300]] are not in
   `SmallGroup` Elements order).
4. **A12/A13 non-abelian generalization**: mitten towers give the first natural
   k-doubling test family; the Bockstein/deck-tower layer currently assumes
   abelian G — the ε-map and transfer arguments should be re-derived for central
   C₂ ≤ Z(G) with G non-abelian (BBEpsFreeGroupAlgebra already covers "k[G] free
   over k[⟨σ⟩]" for central σ, abelian-G hypothesis likely removable).
5. **Odd-center fiberings** (not run): [[150]]/[[500]]/[[780]]/[[975]] have odd
   central subgroups (C₅/C₅/C₁₃/C₅), so F₂[G] Maschke-splits along them — the
   A22-style CRT fibering applies with NO cover subtleties. Relevant to any Lean
   floor attempt on [[150,30,10]] (2-part of C₅×S₃ is only C₂).

## §7 Reproduction

```
uv run --project experiments/bb_lab python experiments/bb_lab/scripts/a26_mitten_descent.py \
    --yarn /path/to/yarn --trials 20000
```
Without `--yarn`/pynauty: Table XIII rebuilds, deck census, and all descents run;
shipped-artifact comparisons and Tanner-iso checks are skipped. Regenerate group
tables with `gap -q --nointeract a26_export_groups.g` (GAP ≥ 4.12; conda-forge
`gap-defaults` works). Runtime: ~15 min single-core at default trials.
