# lp_lab — non-BB lifted-product codes with the BB hardware package

Research lab (uv project, depends on `../bb_lab`) for the LP1 program: weight-6 quantum LDPC
codes built as lifted products `LP(A, b⊗I)` with a 2×2 base matrix `A = [[p,1],[1,q]]` over
F₂[G] that are provably **not** two-block group-algebra / bivariate-bicycle codes but inherit
the thickness-2 Tanner graph, the depth-7 syndrome cycle and (for suitable parameters) a toric
layout of arXiv:2308.07915.

Start with `notes/LP1_results.md` (findings and certificate status), `notes/LP1_theory.md`
(statements and proofs), `notes/LP1_plan.md` (chronological log incl. what did not work), and the
self-contained report page `notes/LP1_report.html`.

Headline codes (code cards with permutation parts, toric coordinates and certificates in
`data/codes/`, Stim circuits alongside):

| code | G | b, p, q | parameters | certified | toric layout |
|---|---|---|---|---|---|
| C224 | Z₂₈×Z₂ | 1+x+x¹²y, 1+x¹⁷, x³y+x⁴ | [[224,12,≤16]] | 14 ≤ d ≤ 16 (SAT) | Z₁₆×Z₂₈ |
| C288 | Z₇₂ | 1+x²+x¹⁰, 1+x⁵, x²¹+x⁵⁸ | [[288,8,≤22]] | d ≥ 14 (SAT), d ≤ 22 (ISD) | Z₃₂×Z₁₈ |

## Layout
* `src/lp_lab/` — `lp.py` (LP construction over abelian group algebras), `fastlp.py` (bit-packed
  builders/k), `design.py` (ideal reduction, kernel codes, maximal-k designs), `distance.py`
  (ISD upper bounds, bb_lab SAT wrapper), `satdist.py` (CryptoMiniSat ladder with translation
  symmetry breaking), `structure.py` (nauty orbits, wheel layers + planarity, toric layouts),
  `circuit.py` (Stim depth-7 memory circuit), `circdist.py` (circuit-level distance: SAT and BP-OSD).
* `scripts/` — `s*_` searches (blind sweeps, design-based), `a*_` analysis/re-screens, `c*_`
  certification (structure, exact d, circuit, BP-OSD d_circ, code cards), `t*_` validations against
  known BB codes, `f1_toric_svg.py` figures, `r1_build_report.py` report page.
* `data/` is gitignored at the repo root; the committed subset is `certs/`, `codes/`, `figs/*.svg`,
  `logs/`, and the design tables for the two headline groups. Sweep dumps are regenerable.

## Running
```bash
cd experiments/lp_lab && uv sync
uv run python scripts/c1_certify.py --orders 28,2 --p 0,34 --q 7,8 --b 0,2,25 --skip-exact --no-circuit
uv run python scripts/c2_satsym.py --orders 28,2 --p 0,34 --q 7,8 --b 0,2,25 --lower 10 --tag C224
uv run python scripts/c4_dcirc_bposd.py --code lp --orders 28,2 --p 0,34 --q 7,8 --b 0,2,25 --cycles 2 --trials 600 --tag C224
```
Indices are group elements in row-major order over the cyclic factors; `--p/--q/--b` give supports.
