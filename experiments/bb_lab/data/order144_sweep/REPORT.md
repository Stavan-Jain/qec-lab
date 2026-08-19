# Order-144 BB breadth sweep — n = 288 (2026-08-17)

Weight-3 x weight-3 BB codes over abelian groups of order 144, sampled enumeration (`bb_samp_` provenance, a18 conventions), canonical Aut x translation x swap dedup, k >= 2, HARD 120 s per-code compute cap. This is a data errand (no approach number); all artifacts live in this directory.

## Headline

- codes processed: **58** (Z12xZ12: 34, Z18xZ8: 6, Z24xZ6: 6, Z36xZ4: 4, Z16xZ9: 4, Z48xZ3: 4)
- exact d: **24** (solver-exact 13, certificate 11); bounded-only: 34
- exact-d histogram: d=6: 2, d=8: 14, d=12: 5, d=16: 3
- bounded-only rows by d_ub: ub=12: 3, ub=18: 3, ub=22: 1, ub=24: 5, ub=26: 2, ub=30: 1, ub=32: 2, ub=34: 3, ub=36: 4, ub=44: 1, ub=48: 5, ub=52: 1, ub=54: 1, ub=62: 1, ub=64: 1
- total per-code wall: 78.5 min (mean 81 s/code)

## Counts by group x outcome

| group | codes | solver-exact | certificate | bounded-only | exact d values |
|---|---|---|---|---|---|
| Z12xZ12 | 34 | 9 | 9 | 16 | [6, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 12, 12, 12, 12, 12, 16] |
| Z18xZ8 | 6 | 1 | 0 | 5 | [8] |
| Z24xZ6 | 6 | 2 | 1 | 3 | [6, 8, 16] |
| Z36xZ4 | 4 | 1 | 0 | 3 | [8] |
| Z16xZ9 | 4 | 0 | 0 | 4 | [] |
| Z48xZ3 | 4 | 0 | 1 | 3 | [16] |

## Notable rows

### High-distance (exact d >= 12 or floor >= 10)

- `bb_samp_Z12xZ12_420d1627` [[288,4]] d = 16 (certificate); A = `1 + y^2 + x`, B = `x*y^5 + x*y^7 + x^3*y^4`; d_ub = 24; floor = None; 21.7 s
- `bb_samp_Z24xZ6_95cbdb6b` [[288,4]] d = 16 (certificate); A = `y + x + x^2`, B = `x^10*y^4 + x^14 + x^15*y^5`; d_ub = 16; floor = None; 11.0 s
- `bb_samp_Z48xZ3_73c33832` [[288,4]] d = 16 (certificate); A = `y + x + x^2`, B = `x^12 + x^26*y^2 + x^34*y`; d_ub = 16; floor = None; 19.9 s
- `bb_samp_Z12xZ12_026cb33a` [[288,8]] d = 12 (certificate); A = `y^3 + x + x^2`, B = `x^2*y^2 + x^4 + x^6*y`; d_ub = 12; floor = None; 15.5 s
- `bb_samp_Z12xZ12_54d10efa` [[288,8]] d = 12 (certificate); A = `1 + y + x`, B = `x^6*y^11 + x^7*y^10 + x^8*y^6`; d_ub = 16; floor = None; 17.5 s
- `bb_samp_Z12xZ12_91d8adba` [[288,8]] d = 12 (certificate); A = `y + y^3 + x^2`, B = `x*y^4 + x*y^5 + x^7`; d_ub = 12; floor = None; 9.3 s
- `bb_samp_Z12xZ12_aa84325d` [[288,4]] d = 12 (certificate); A = `1 + y + x`, B = `y^4 + y^6 + x^3*y^2`; d_ub = 12; floor = None; 23.4 s
- `bb_samp_Z12xZ12_f88ac006` [[288,4]] d = 12 (certificate); A = `1 + y^2 + x`, B = `x*y^2 + x^3*y + x^11*y^6`; d_ub = 12; floor = None; 16.7 s
- `bb_samp_Z48xZ3_0d33e3ce` [[288,4]] UNKNOWN in [10, 32]; A = `y + x + x^2`, B = `x^2*y + x^3*y^2 + x^46`; d_ub = 32; floor = 10; 120.0 s

### Open tail (bounded-only, d_ub >= 14) — the interesting frontier

- `bb_samp_Z24xZ6_d6396252` [[288,4]] UNKNOWN in [8, 64]; A = `1 + y + x^2`, B = `x^3*y^2 + x^18*y^4 + x^20*y`; d_ub = 64; floor = 8; 120.0 s
- `bb_samp_Z18xZ8_f3d6801a` [[288,4]] UNKNOWN in [8, 62]; A = `y^2 + x + x^2`, B = `x^5 + x^6*y^2 + x^7*y^3`; d_ub = 62; floor = 8; 120.0 s
- `bb_samp_Z12xZ12_cfbc9578` [[288,4]] UNKNOWN in [8, 54]; A = `1 + y + x`, B = `x*y^10 + x^6*y^7 + x^6*y^8`; d_ub = 54; floor = 8; 120.1 s
- `bb_samp_Z36xZ4_45c18b57` [[288,4]] UNKNOWN in [9, 52]; A = `1 + x + x^2`, B = `x^5*y + x^7 + x^18`; d_ub = 52; floor = 9; 120.0 s
- `bb_samp_Z18xZ8_a8f8a9cd` [[288,4]] UNKNOWN in [8, 48]; A = `y + x + x^2`, B = `x^13*y^7 + x^14*y^2 + x^15*y^3`; d_ub = 48; floor = 8; 120.0 s
- `bb_samp_Z18xZ8_feee3165` [[288,4]] UNKNOWN in [8, 48]; A = `y^3 + x*y + x^2`, B = `x^7*y^2 + x^11*y + x^15*y`; d_ub = 48; floor = 8; 120.0 s
- `bb_samp_Z16xZ9_2f686b84` [[288,4]] UNKNOWN in [8, 48]; A = `y^8 + x*y + x^3`, B = `y + x^6*y^3 + x^12*y^2`; d_ub = 48; floor = 8; 120.0 s
- `bb_samp_Z16xZ9_8db79545` [[288,4]] UNKNOWN in [9, 48]; A = `y^2 + x*y + x^2`, B = `1 + x^8*y^4 + x^10*y^2`; d_ub = 48; floor = 9; 120.0 s
- `bb_samp_Z48xZ3_d411dd81` [[288,4]] UNKNOWN in [8, 48]; A = `1 + x*y + x^4`, B = `x^18*y + x^30*y^2 + x^44*y`; d_ub = 48; floor = 8; 120.1 s
- `bb_samp_Z12xZ12_918c97a8` [[288,12]] UNKNOWN in [8, 44]; A = `1 + y + x`, B = `x*y^7 + x^4*y^2 + x^5*y^4`; d_ub = 44; floor = 8; 120.0 s
- `bb_samp_Z12xZ12_3b47338f` [[288,4]] UNKNOWN in [9, 36]; A = `1 + y + x`, B = `y^4 + x^9*y^11 + x^10*y^4`; d_ub = 36; floor = 9; 120.0 s
- `bb_samp_Z12xZ12_57bd841d` [[288,4]] UNKNOWN in [8, 36]; A = `1 + y + x`, B = `x^2*y^6 + x^3*y^8 + x^5*y^2`; d_ub = 36; floor = 8; 120.1 s
- `bb_samp_Z12xZ12_9033bcff` [[288,4]] UNKNOWN in [8, 36]; A = `1 + y^2 + x`, B = `x + x*y^5 + x^3*y^2`; d_ub = 36; floor = 8; 120.0 s
- `bb_samp_Z12xZ12_c8fc555a` [[288,4]] UNKNOWN in [8, 36]; A = `1 + y^2 + x`, B = `y^5 + x^3 + x^5*y^8`; d_ub = 36; floor = 8; 120.0 s
- `bb_samp_Z18xZ8_e872aa7c` [[288,8]] UNKNOWN in [7, 34]; A = `y + x^2 + x^4`, B = `x*y^6 + x^11*y + x^15*y^3`; d_ub = 34; floor = 7; 120.0 s
- `bb_samp_Z16xZ9_b2928324` [[288,4]] UNKNOWN in [9, 34]; A = `y + y^2 + x^2`, B = `y^6 + x*y^7 + x^2*y^8`; d_ub = 34; floor = 9; 120.0 s
- `bb_samp_Z48xZ3_d35a61e9` [[288,4]] UNKNOWN in [9, 34]; A = `1 + x + x^2`, B = `x^5*y + x^25*y + x^27`; d_ub = 34; floor = 9; 120.0 s
- `bb_samp_Z36xZ4_241c221e` [[288,8]] UNKNOWN in [8, 32]; A = `y + x^2 + x^4`, B = `x^7*y^3 + x^17*y^2 + x^27`; d_ub = 32; floor = 8; 120.0 s
- `bb_samp_Z48xZ3_0d33e3ce` [[288,4]] UNKNOWN in [10, 32]; A = `y + x + x^2`, B = `x^2*y + x^3*y^2 + x^46`; d_ub = 32; floor = 10; 120.0 s
- `bb_samp_Z24xZ6_35c111f7` [[288,12]] UNKNOWN in [9, 30]; A = `y^3 + x + x^2`, B = `x*y^2 + x^2*y + x^6`; d_ub = 30; floor = 9; 120.0 s
- `bb_samp_Z12xZ12_9d9a89f0` [[288,4]] UNKNOWN in [6, 26]; A = `y + y^3 + x^2`, B = `x^4*y^7 + x^5*y^11 + x^9`; d_ub = 26; floor = 6; 120.0 s
- `bb_samp_Z12xZ12_ad8120f1` [[288,4]] UNKNOWN in [7, 26]; A = `y + y^3 + x^2`, B = `x^5*y^3 + x^6*y + x^10*y^2`; d_ub = 26; floor = 7; 120.0 s
- `bb_samp_Z12xZ12_295acff6` [[288,4]] UNKNOWN in [7, 24]; A = `1 + y + x`, B = `x*y + x^2*y^10 + x^8`; d_ub = 24; floor = 7; 120.0 s
- `bb_samp_Z12xZ12_2f1371cc` [[288,4]] UNKNOWN in [8, 24]; A = `1 + y + x`, B = `y^5 + x^5*y^11 + x^11*y^3`; d_ub = 24; floor = 8; 120.0 s
- `bb_samp_Z12xZ12_a1930091` [[288,4]] UNKNOWN in [8, 24]; A = `1 + y^2 + x`, B = `x^2*y^11 + x^6*y^3 + x^8`; d_ub = 24; floor = 8; 120.0 s
- `bb_samp_Z18xZ8_83a40cd3` [[288,4]] UNKNOWN in [7, 24]; A = `y + x*y + x^5`, B = `x^2*y + x^15*y^3 + x^16*y^5`; d_ub = 24; floor = 7; 120.0 s
- `bb_samp_Z24xZ6_a2de1fcb` [[288,4]] UNKNOWN in [9, 24]; A = `y + x + x^2`, B = `x^3*y^3 + x^4 + x^20*y^3`; d_ub = 24; floor = 9; 120.0 s
- `bb_samp_Z12xZ12_171360b8` [[288,4]] UNKNOWN in [8, 22]; A = `1 + y + x`, B = `y^3 + x^2 + x^6*y^8`; d_ub = 22; floor = 8; 120.0 s
- `bb_samp_Z12xZ12_50a94227` [[288,4]] UNKNOWN in [8, 18]; A = `1 + y^2 + x`, B = `x^3*y^7 + x^6*y^2 + x^6*y^3`; d_ub = 18; floor = 8; 120.0 s
- `bb_samp_Z12xZ12_83410680` [[288,4]] UNKNOWN in [8, 18]; A = `1 + y^4 + x`, B = `x^2*y^4 + x^2*y^5 + x^5*y^3`; d_ub = 18; floor = 8; 120.0 s
- `bb_samp_Z12xZ12_ecf09d8e` [[288,4]] UNKNOWN in [8, 18]; A = `1 + y^4 + x`, B = `x^3*y^7 + x^4*y + x^5*y^10`; d_ub = 18; floor = 8; 120.0 s

### Certified doubles

- `bb_samp_Z12xZ12_026cb33a` [[288,8]] d = 12 (certificate); A = `y^3 + x + x^2`, B = `x^2*y^2 + x^4 + x^6*y`; d_ub = 12; floor = None; 15.5 s
  - base: {'group': [12, 6], 'axis': 1, 'n': 144, 'k': 8, 'kappa': 68}, d_base: 6
- `bb_samp_Z12xZ12_54d10efa` [[288,8]] d = 12 (certificate); A = `1 + y + x`, B = `x^6*y^11 + x^7*y^10 + x^8*y^6`; d_ub = 16; floor = None; 17.5 s
  - base: {'group': [6, 12], 'axis': 0, 'n': 144, 'k': 8, 'kappa': 68}, d_base: 6
- `bb_samp_Z12xZ12_91d8adba` [[288,8]] d = 12 (certificate); A = `y + y^3 + x^2`, B = `x*y^4 + x*y^5 + x^7`; d_ub = 12; floor = None; 9.3 s
  - base: {'group': [12, 6], 'axis': 1, 'n': 144, 'k': 8, 'kappa': 68}, d_base: 6
- `bb_samp_Z12xZ12_420d1627` [[288,4]] d = 16 (certificate); A = `1 + y^2 + x`, B = `x*y^5 + x*y^7 + x^3*y^4`; d_ub = 24; floor = None; 21.7 s
  - base: {'group': [12, 6], 'axis': 1, 'n': 144, 'k': 4, 'kappa': 70}, d_base: 8
- `bb_samp_Z12xZ12_aa84325d` [[288,4]] d = 12 (certificate); A = `1 + y + x`, B = `y^4 + y^6 + x^3*y^2`; d_ub = 12; floor = None; 23.4 s
  - base: {'group': [6, 12], 'axis': 0, 'n': 144, 'k': 4, 'kappa': 70}, d_base: 6
- `bb_samp_Z12xZ12_b824336c` [[288,4]] d = 8 (certificate); A = `1 + y + y^5`, B = `y^11 + x*y + x^6*y`; d_ub = 24; floor = None; 8.7 s
  - base: {'group': [12, 6], 'axis': 1, 'n': 144, 'k': 4, 'kappa': 70}, d_base: 4
- `bb_samp_Z12xZ12_d5c70349` [[288,4]] d = 8 (certificate); A = `1 + y + y^5`, B = `x^2*y^10 + x^3*y + x^8*y^3`; d_ub = 48; floor = None; 6.2 s
  - base: {'group': [12, 6], 'axis': 1, 'n': 144, 'k': 4, 'kappa': 70}, d_base: 4
- `bb_samp_Z12xZ12_f1fb28fb` [[288,4]] d = 8 (certificate); A = `1 + y + y^2`, B = `x^3 + x^4*y + x^6*y`; d_ub = 32; floor = None; 17.6 s
  - base: {'group': [12, 6], 'axis': 1, 'n': 144, 'k': 4, 'kappa': 70}, d_base: 4
- `bb_samp_Z12xZ12_f88ac006` [[288,4]] d = 12 (certificate); A = `1 + y^2 + x`, B = `x*y^2 + x^3*y + x^11*y^6`; d_ub = 12; floor = None; 16.7 s
  - base: {'group': [12, 6], 'axis': 1, 'n': 144, 'k': 4, 'kappa': 70}, d_base: 6
- `bb_samp_Z24xZ6_95cbdb6b` [[288,4]] d = 16 (certificate); A = `y + x + x^2`, B = `x^10*y^4 + x^14 + x^15*y^5`; d_ub = 16; floor = None; 11.0 s
  - base: {'group': [12, 6], 'axis': 0, 'n': 144, 'k': 4, 'kappa': 70}, d_base: 8
- `bb_samp_Z48xZ3_73c33832` [[288,4]] d = 16 (certificate); A = `y + x + x^2`, B = `x^12 + x^26*y^2 + x^34*y`; d_ub = 16; floor = None; 19.9 s
  - base: {'group': [24, 3], 'axis': 0, 'n': 144, 'k': 4, 'kappa': 70}, d_base: 8

### Anomalies

(none)

## Certify-lane outcome mix

- base out of scope (d_ub > 15): 22
- no (R) candidate: 12
- CERTIFIED: 11
- census too large at required W (front-end scope): 7
- DOUBLING-REFUTED: 3
- front-end k cap (k > 14): 2
- error: 1

Notes: the 3 DOUBLING-REFUTED verdicts are certified negatives ('safe-class coset of weight d_base exists', so the safe-floor doubling fails on that axis); two of those rows were independently closed by SAT at exactly d = 8, consistent with the refutation's 'd < 2*d_base on the safe sector'. The 1 error row is an infrastructure timeout of certify's inner `cosetbz` native subprocess (bb_samp_Z36xZ4_23206eb7); the row fell through to the SAT lane normally.

## Base cache (shared quotient ladder, order 72 / 36 / 18)

| quotient group | entries | with d_exact | from corpus | compute wall s |
|---|---|---|---|---|
| Z12xZ3 | 44 | 44 | 44 | 0.0 |
| Z12xZ6 | 34 | 7 | 1 | 74.7 |
| Z18xZ1 | 9 | 9 | 0 | 1.7 |
| Z18xZ2 | 6 | 6 | 0 | 6.0 |
| Z18xZ4 | 6 | 1 | 0 | 46.5 |
| Z24xZ3 | 10 | 1 | 0 | 15.0 |
| Z2xZ9 | 4 | 4 | 0 | 0.7 |
| Z36xZ1 | 4 | 4 | 0 | 2.6 |
| Z36xZ2 | 4 | 3 | 0 | 23.6 |
| Z4xZ9 | 4 | 4 | 4 | 0.0 |
| Z6xZ3 | 38 | 38 | 38 | 0.0 |
| Z8xZ9 | 4 | 1 | 0 | 2.1 |

Cache: 167 distinct quotient codes, 172.9 s total compute (amortized across covers; corpus hits free). Ladder d values are context, NOT floors/bounds for the cover (no bound transport claimed).

## Method / trust-tier legend

- **solver-exact** — CMS SAT ladder in `bb_lab.sat_distance`: witness at d + solver-proved UNSAT at every w < d (solver trust, same tier as corpus CaDiCaL d_exact). Rows closed as 'UNSAT floor + verified L1 witness' are the same two sides established by different engines.
- **certificate** — `bb_lab.doubling_certify` (A30 front-end): d = 2 * d_base by the doubling theorem + counting-invariant enumeration; certificate tier, not kernel-checked.
- **bounded-only** — d in [floor, d_ub]: floor only when the solver genuinely refuted all weights < floor (contiguous UNSAT prefix); d_ub from L1 sampling with a verified logical witness. Witness weights are never reported as floors or exact. UNKNOWN stays UNKNOWN.

## Smoke tests (pre-batch)

- fast canonicalization == `bb_lab.canonical.canonical_bits` on 180 random pairs (Z6xZ6, Z12xZ3, Z6xZ4): 0 mismatches; all 13 corpus Z12xZ12 exact rows round-trip to their instance_ids.
- SAT lane: gross base bb72 [[72,12,6]] -> DISTANCE 6, clean UNSAT ladder, 0.5 s.
- certify lane: gross [[144,12,12]] re-derived CERTIFIED d = 12 (= 2 * 6) in 3.7 s.
- L1 at n = 288: 100k samples ~5 s; finds d_ub = 18 on [[288,12,18]] (its true d).

## Caveats

- Sampled (not exhaustive) enumeration: `bb_samp_` provenance; orbit-size-biased uniform sampling per a18 conventions.
- Per-code 120 s cap means the interesting tail (d_ub >= 12) is mostly bounded-only; floors reflect only solver-completed UNSAT rounds at n = 288 (w <= ~9 typically).
- L1 d_ub at n = 288 with 100k samples is loose for low-d codes (SAT closed those anyway); for the open tail treat d_ub as an upper bound, not an estimate.
- Certify-lane REFUSED/FALLBACK outcomes are scope refusals, not refutations, EXCEPT explicit DOUBLING-REFUTED verdicts (certified negatives for that axis).
- Quotient-ladder d values are informational context; no cover-distance bound is claimed from them.
- Machine was contended by other sessions' solvers (load ~5-6); wall times are upper estimates of uncontended cost.

## Reproduction

```bash
cd experiments/bb_lab
uv run python data/order144_sweep/smoke.py          # lane smokes
uv run python data/order144_sweep/sampler.py        # sampled enumeration -> sweep.duckdb
caffeinate -ims uv run python data/order144_sweep/driver.py   # batch (resumable; skips done rows)
uv run python data/order144_sweep/report.py         # this report
```

Main corpus DB was opened READ-ONLY for dedup + quotient lookups; nothing under the main checkout was written. Skipped as already-exact: 13 corpus Z12xZ12 rows (incl. [[288,12,18]]) + published IBM [[288,8,20]] (Z18xZ8) by canonical id.
