# bb_lab notes — approach-number registry

**Rule: approach numbers are claimed HERE, in the same commit that creates the
first `A<N>_*.md` note.** Parallel sessions must not pick a number without
updating this file — two independent "A15" lines were created in July 2026 and
had to be renumbered after the fact (the d≥7 hunt is now A17). Branch names do
not reserve numbers (the branch `claude/a8-literal-lift-criterion` hosts the
A11 line, not A8).

**Next free number: A31.**

| # | Line | Note file(s) | Status |
|---|------|--------------|--------|
| A0 | Analytic baseline scoreboard (Bravyi instances) | `A0_baseline.md` | done |
| A1 | Literature sweep, four lanes + synthesis | `A1_literature_L{1-4}.md`, `A1_synthesis.md` | done |
| A2 | Phase-2 scouting synthesis (track ranking) | `A2_scouting.md` | done |
| A3 | Track 1.1 deep push (Smith h=2 cover transfer) — the goal-1 lab log | `A3_track1p1_log.md` | done (d(gross)=12 analytic) |
| A4 | The write-up of record: d ≥ 6 and d = 12 fully analytic | `A4_writeup.md` | done |
| A5 | Goal 2: analytic bounds for a class of BB codes — running log | `A5_goal2_log.md` | active (Entries 8+ = A15 line) |
| A6 | Lean d=12 finish plan | `A6_lean_d12_finish_plan.md` | done |
| A7 | d=12 finish plan (Tier-3 analytic replacement) | `A7_d12_finish_plan.md` | active (Props 30–31 open) |
| A8 | Doubling extension beyond gross: [[336,12,12]] over Z₃×Z₇ | `A8_doubling_extension_writeup.md` | done |
| A9 | Lean-target screen (stage 1 of the doubling-layer plan) | `A9_lean_target_screen.md` | done |
| A10 | Descent-cover twist screen (doubling-existence question) | `A10_descent_twist_screen.md` | merged 2026-07-18 (annotated) |
| A11 | Literal-lift doubling criterion (C-safe ∧ C-danger) | `A11_*` | merged 2026-07-18 (annotated) |
| A12 | Deck homotopy (R) ⟺ k-preservation ⟺ Bezout | `A12_deck_homotopy_R.md` | done |
| A13 | Deck towers (OQ1) + Bockstein equality (OQ2) + L2 Lean wiring | `A13_deck_tower_plan.md`, `A13_bockstein_equality_plan.md`, `A13_result.md`, `A13_L2_formalization_plan.md` | done (seamC↔δ₂ transport CLOSED 2026-07-20; only the `D`-module structure iso remains) |
| A14 | Safe-floor criterion (OQ4): screens S0–S4, tower bottleneck | `A14_safe_floor_criterion_plan.md` | done |
| A15 | Base-floor generalization: the class small-cycle theorem | `A15_base_floor_class_plan.md` | active (w=5 kills in flight) |
| A16 | Class small-cycle theorem — statement and proof of record | `A16_class_theorem_writeup.md` | done (unconditional) |
| A17 | The d ≥ 7 doubling hunt (ex-"A15"): [[252,8,16]] doubles, [[300,8,16]] two-tier, deficit wall (P3) | `A17_d7plus_doubling_hunt_plan.md`, `A17_deficit_wall.md` | active (near-kernel stratum classified) |
| A18 | Corpus breadth sweep (ex-"a16"): μ_e barrenness, 41 group shapes | `A18_breadth_sweep.md` | done |
| A19 | Bravyi [[360,12,≤24]] tower: deck-nontrivial on all decks, bases [[90,8,8]]/[[180,8,12]]², τ-lift weight table, **M12 CERTIFIED: 12 ≤ d ≤ 24** | `A19_bravyi360_tower.md` | active (M12 done 2026-07-22 §8; next: (M)@24 census + BX close-out; doubly-old behind A20's seam trial) |
| A20 | IBM [[288,8,20]] (class Y) as an all-(R) tower: k born (18,2), Bezout witnesses both y-rungs, d₀ = 6 exact, cheapest d ≥ 20 target | `A20_ibm288_tower.md` | active (d₁ ladder in flight; template obligations next) |
| A21 | Analytic base floor for f2a6f17e: LogicalFloor 8 via the class small-cycle machinery | `A21_analytic_base_floor.md` | active (spawned 2026-07-21) |
| A22 | Analytic light-boundary completeness: LightClassification via σ-correlation/Fourier structure | `A22_analytic_classification.md` | **DONE** (2026-07-22): `lightClassification` sorry-free via CRT fibering (z=y³,w=y⁵), branch `claude/a22-light-classification` |
| A23 | Analytic seam-coset floor: SeamCosetFloor 16 via parity + transport + difference sets | `A23_analytic_seam_floor.md` | **DONE** (2026-07-22): `seamCosetFloor_16` sorry-free via A22-fibering site sweep, branch `claude/a23-seam-transfer` |
| A24 | Bravyi-360 y-deck safe-sector floor at 24 (p_y*≠0 sector; subsumes new-x + doubly-old), falsify-first + A23-fibering port | `A24_y_safe_floor.md` | active (spawned 2026-07-23; note lives on branch `claude/bravyi-y-deck-floor-bbf8cf`, not yet on main) |
| A25 | gross d_circ via SBB Thm-1 extended-code reduction (d_ext^Z ≤ 9 verified; exact-d_ext UNSAT probes) | note on unmerged branch — row added defensively 2026-08-03 to prevent number collision | active |
| A26 | Mitten-code descent (arXiv:2607.28795): decks/quotients/towers on non-abelian LP codes; tower hypothesis refuted, Table XIII [[300]] erratum, 2d̄−2 wall echo ×2 | `A26_mitten_descent.md` | done (ISD tier; exact-d + safe-floor follow-ups queued) |
| A27 | Safe-floor deletion: generality tiering + fibering feasibility probe for the Z₁₅×Z₆ [[180,4,10]] docket UNKNOWNs | `A27_safe_floor_generality.md` | done (synthesis + probe, 2026-08-06; §3.3 P5 sweep-size estimate corrected by A29 §3 erratum) |
| A28 | LSC theory: certified BZ census lane (docket [[180,4,10]] censuses NEW), shift-bound verdict (measured negative w/ certified gap cell), ε-trisection theorem (validated 4,574/4,574) | `A28_light_classification_theory.md` | done (2026-08-06 incl. depth-10 dichotomy curve; §6 gates data-backed for both docket cells) |
| A29 | General fibering engine (`bb_lab.fibering`): portable safe-floor certification; 3 docket UNKNOWNs certified (2 pinned exact), A8 §4.3 open core CLOSED ([[168,12,6]] safe floor 12 in 10 s); ε-recursion chapter = the named residue | `A29_general_fibering.md` | active (session 1, 2026-08-06; engine validated bit-for-bit vs A22/A23; **renumbered from A28 at merge** — claimed A28 in parallel with the LSC line, LSC claim was 31 min earlier) |
| A30 | Coset-BZ (A28 windows × A29 seam offsets): safe floors as coset weights; **all 3 remaining docket UNKNOWNs certified** — both Z₁₅×Z₆ [[180,4,10]] codes double on every axis ([[360,4,20]], first d=10→20, deficit wall cleared) + Z₂₁×Z₃ CMS cross-check in 0.6 s (A29 §5.3 "4 UNKNOWNs" = stale-row miscount); docket CLOSED, ≤ 8.4 min/code | `A30_coset_bz_doubling_certificates.md` | done (2026-08-07; certificate tier; Lean packaging = follow-on 1) |

Other prefixes: `S<NN>` = engineering/packaging sessions (e.g. `S39` = pair72
packaging); `T<N>*` = the May-era Tier-2/Tier-3 classifier and conjecture-mill
clusters (see `HANDOFF.md` for that program's structure); `A_HANDOFF.md` = the
Phase-A program handoff (read first).
