# A38 S1 — corpus envelope census (the empirical burden map)

Generated 2026-08-18 by `scripts/a38_envelope_census.py` (pricing only — no census was run; verdicts are COST verdicts, never distance claims).

Corpus: 58350 rows / 47 group shapes (A18 recorded 41; grown since). Representatives: best exact-d per group + frontier d_ub per group where it out-ranks exact + 5 zoo BB instances + 11 A35-docket calibration rows (all reproduce the banked verdicts exactly).

Verdicts: **AMBER**: 3, **GREEN**: 41, **RED**: 6, **no-deck**: 23

## Headlines

1. **Every exact-d corpus representative with a deck prices GREEN** (27/27): re-certifying the corpus's known distances is entirely within the demonstrated envelope wherever 2 | |G|.
2. **The open frontier is CAP-bound, not node-bound**: of the 8 open d_ub questions with a deck, 8 are AMBER/RED with bottom censuses *within* the 2e11 node envelope but caps 9-17 beyond the demonstrated 8; 0 are node-bound. The corpus frontier wall is W3 (fiber-cap growth — the F2b/F2c target), not W2 (census blowup): [[756]]-style census walls are atypical in-corpus.
3. **W1 (odd |G|) locks out 22 of 57 corpus rows**, including the highest-d_ub open questions (Z9xZ9 ub 40, Z7xZ9 ub 30): F1 demand is real.
4. cap_max < 0 on tiny questions (W below the lightest stabilizer weight) means the dangerous sector is empty at that budget — a pricing artifact, GREEN by construction.

| code | src | group | [[n,k,d]] | v2/axis | odd part | depth | k-chain | W | log10 nodes/level | cap | 2^k(base) | verdict |
|---|---|---|---|---|---|---|---|---|---|---|---|
| Z3xZ3:best-exact | corpus | 3x3 | d exact 4 | [0, 0] | 9 | 0 | — | — | — | — | no-deck |
| Z3xZ4:best-exact | corpus | 3x4 | [[24,4,exact 4]] | [0, 2] | 3 | 2 | [4, 4, 4] | 2 | [1.0, 0.6, 0.0] | -2 | 2^4 | GREEN |
| Z3xZ5:best-exact | corpus | 3x5 | d exact 6 | [0, 0] | 15 | 0 | — | — | — | — | no-deck |
| Z3xZ6:best-exact | corpus | 3x6 | [[36,4,exact 6]] | [0, 1] | 9 | 1 | [4, 4] | 4 | [2.2, 1.5] | -1 | 2^4 | GREEN |
| Z3xZ7:best-exact | corpus | 3x7 | d exact 6 | [0, 0] | 21 | 0 | — | — | — | — | no-deck |
| Z3xZ8:best-exact | corpus | 3x8 | [[48,4,exact 8]] | [0, 3] | 3 | 3 | [4, 4, 4, 4] | 6 | [3.3, 2.4, 1.4, 0.3] | 0 | 2^4 | GREEN |
| Z4xZ6:best-exact | corpus | 4x6 | [[48,8,exact 4]] | [2, 1] | 3 | 3 | [8, 8, 4, 4] | 2 | [1.3, 0.9, 0.6, 0.0] | -2 | 2^4 | GREEN |
| Z3xZ9:best-exact | corpus | 3x9 | d exact 8 | [0, 0] | 27 | 0 | — | — | — | — | no-deck |
| Z4xZ7:best-exact | corpus | 4x7 | [[56,6,exact 8]] | [2, 0] | 7 | 2 | [6, 6, 6] | 6 | [3.5, 2.5, 1.4] | 0 | 2^6 | GREEN |
| Z5xZ6:best-exact | corpus | 5x6 | [[60,4,exact 8]] | [0, 1] | 15 | 1 | [4, 4] | 6 | [3.6, 2.7] | 0 | 2^4 | GREEN |
| Z5xZ7:best-exact | corpus | 5x7 | d exact 8 | [0, 0] | 35 | 0 | — | — | — | — | no-deck |
| Z3xZ12:best-exact | corpus | 3x12 | [[72,8,exact 8]] | [0, 2] | 9 | 2 | [8, 8, 8] | 6 | [3.8, 2.8, 1.6] | 0 | 2^8 | GREEN |
| Z4xZ9:best-exact | corpus | 4x9 | [[72,4,exact 10]] | [2, 0] | 9 | 2 | [4, 4, 4] | 8 | [4.8, 3.5, 2.2] | 1 | 2^4 | GREEN |
| Z6xZ6:best-exact | corpus | 6x6 | [[72,4,exact 8]] | [1, 1] | 9 | 2 | [4, 4, 4] | 6 | [3.9, 2.9, 2.0] | 0 | 2^4 | GREEN |
| Z3xZ13:best-exact | corpus | 3x13 | d exact 10 | [0, 0] | 39 | 0 | — | — | — | — | no-deck |
| Z6xZ7:best-exact | corpus | 6x7 | [[84,6,exact 10]] | [1, 0] | 21 | 1 | [6, 6] | 8 | [5.0, 3.7] | 1 | 2^6 | GREEN |
| Z15xZ3:best-exact | corpus | 15x3 | d exact 10 | [0, 0] | 45 | 0 | — | — | — | — | no-deck |
| Z3xZ15:best-exact | corpus | 3x15 | d exact 8 | [0, 0] | 45 | 0 | — | — | — | — | no-deck |
| Z5xZ9:best-exact | corpus | 5x9 | d exact 12 | [0, 0] | 45 | 0 | — | — | — | — | no-deck |
| Z4xZ12:best-exact | corpus | 4x12 | [[96,4,exact 12]] | [2, 2] | 3 | 4 | [4, 4, 4, 4, 4] | 10 | [6.2, 4.6, 3.0, 1.5, 0.3] | 2 | 2^4 | GREEN |
| Z7xZ7:best-exact | corpus | 7x7 | d exact 12 | [0, 0] | 49 | 0 | — | — | — | — | no-deck |
| Z9xZ6:best-exact | corpus | 9x6 | [[108,4,exact 12]] | [0, 1] | 27 | 1 | [4, 4] | 10 | [6.5, 4.9] | 2 | 2^4 | GREEN |
| Z7xZ8:best-exact | corpus | 7x8 | [[112,6,exact 12]] | [0, 3] | 7 | 3 | [6, 6, 6, 6] | 10 | [6.5, 4.9, 3.2, 1.5] | 2 | 2^6 | GREEN |
| Z5xZ12:best-exact | corpus | 5x12 | [[120,8,exact 12]] | [0, 2] | 15 | 2 | [8, 8, 8] | 10 | [6.7, 5.0, 3.2] | 2 | 2^8 | GREEN |
| Z6xZ10:best-exact | corpus | 6x10 | [[120,8,exact 12]] | [1, 1] | 15 | 2 | [8, 8, 8] | 10 | [6.7, 5.0, 3.2] | 2 | 2^8 | GREEN |
| Z21xZ3:best-exact | corpus | 21x3 | d exact 10 | [0, 0] | 63 | 0 | — | — | — | — | no-deck |
| Z3xZ21:best-exact | corpus | 3x21 | d exact 8 | [0, 0] | 63 | 0 | — | — | — | — | no-deck |
| Z7xZ9:best-exact | corpus | 7x9 | d exact 14 | [0, 0] | 63 | 0 | — | — | — | — | no-deck |
| Z12xZ6:best-exact | corpus | 12x6 | [[144,12,exact 12]] | [2, 1] | 9 | 3 | [12, 12, 8, 8] | 10 | [7.0, 5.3, 3.7, 1.8] | 2 | 2^8 | GREEN |
| Z8xZ9:best-exact | corpus | 8x9 | [[144,8,exact 10]] | [3, 0] | 9 | 3 | [8, 8, 8, 4] | 8 | [6.0, 4.7, 3.3, 2.2] | 1 | 2^4 | GREEN |
| Z15xZ5:best-exact | corpus | 15x5 | d exact 10 | [0, 0] | 75 | 0 | — | — | — | — | no-deck |
| Z5xZ15:best-exact | corpus | 5x15 | d exact 12 | [0, 0] | 75 | 0 | — | — | — | — | no-deck |
| Z7xZ11:best-exact | corpus | 7x11 | d exact 16 | [0, 0] | 77 | 0 | — | — | — | — | no-deck |
| Z6xZ13:best-exact | corpus | 6x13 | [[156,8,exact 10]] | [1, 0] | 39 | 1 | [8, 4] | 8 | [6.1, 4.9] | 1 | 2^4 | GREEN |
| Z9xZ9:best-exact | corpus | 9x9 | d exact 12 | [0, 0] | 81 | 0 | — | — | — | — | no-deck |
| Z6xZ14:best-exact | corpus | 6x14 | [[168,6,exact 16]] | [1, 1] | 21 | 2 | [6, 6, 6] | 14 | [9.6, 7.4, 5.0] | 4 | 2^6 | GREEN |
| Z7xZ12:best-exact | corpus | 7x12 | [[168,6,exact 16]] | [0, 2] | 21 | 2 | [6, 6, 6] | 14 | [9.6, 7.4, 5.0] | 4 | 2^6 | GREEN |
| Z15xZ6:best-exact | corpus | 15x6 | [[180,16,exact 10]] | [0, 1] | 45 | 1 | [16, 8] | 8 | [6.3, 5.1] | 1 | 2^8 | GREEN |
| Z6xZ15:best-exact | corpus | 6x15 | [[180,40,exact 4]] | [1, 0] | 45 | 1 | [40, 40] | 2 | [1.8, 1.4] | -2 | 2^40 | GREEN |
| Z9xZ10:best-exact | corpus | 9x10 | [[180,4,exact 18]] | [0, 1] | 45 | 1 | [4, 4] | 16 | [10.9, 8.3] | 5 | 2^4 | GREEN |
| Z7xZ14:best-exact | corpus | 7x14 | [[196,6,exact 16]] | [0, 1] | 49 | 1 | [6, 6] | 14 | [10.1, 7.9] | 4 | 2^6 | GREEN |
| Z3xZ35:best-exact | corpus | 3x35 | d exact 2 | [0, 0] | 105 | 0 | — | — | — | — | no-deck |
| Z5xZ21:best-exact | corpus | 5x21 | d exact 18 | [0, 0] | 105 | 0 | — | — | — | — | no-deck |
| Z7xZ15:best-exact | corpus | 7x15 | d exact 18 | [0, 0] | 105 | 0 | — | — | — | — | no-deck |
| Z6xZ18:best-exact | corpus | 6x18 | [[216,4,exact 18]] | [1, 1] | 27 | 2 | [4, 4, 4] | 16 | [11.5, 9.0, 6.4] | 5 | 2^4 | GREEN |
| Z3xZ42:best-exact | corpus | 3x42 | [[252,6,exact 18]] | [0, 1] | 63 | 1 | [6, 6] | 16 | [12.1, 9.5] | 5 | 2^6 | GREEN |
| Z12xZ12:best-exact | corpus | 12x12 | [[288,12,exact 18]] | [2, 2] | 9 | 4 | [12, 12, 8, 8, 8] | 16 | [12.5, 9.9, 7.3, 4.4, 1.8] | 5 | 2^8 | GREEN |
| Z9xZ6:frontier-ub | corpus | 9x6 | [[108,4,ub 26]] | [0, 1] | 27 | 1 | [4, 4] | 24 | [11.6, 7.5] | 9 | 2^4 | AMBER |
| Z5xZ12:frontier-ub | corpus | 5x12 | [[120,4,ub 30]] | [0, 2] | 15 | 2 | [4, 4, 4] | 28 | [13.3, 8.4, 4.2] | 11 | 2^4 | AMBER |
| Z6xZ10:frontier-ub | corpus | 6x10 | [[120,4,ub 30]] | [1, 1] | 15 | 2 | [4, 4, 4] | 28 | [13.3, 8.4, 4.2] | 11 | 2^4 | AMBER |
| Z7xZ9:frontier-ub | corpus | 7x9 | d ub 30 | [0, 0] | 63 | 0 | — | — | — | — | no-deck |
| Z12xZ6:frontier-ub | corpus | 12x6 | [[144,4,ub 36]] | [2, 1] | 9 | 3 | [4, 4, 4, 4] | 34 | [16.2, 10.2, 5.1, 2.4] | 14 | 2^4 | RED |
| Z8xZ9:frontier-ub | corpus | 8x9 | [[144,4,ub 36]] | [3, 0] | 9 | 3 | [4, 4, 4, 4] | 34 | [16.2, 10.2, 5.1, 2.4] | 14 | 2^4 | RED |
| Z6xZ13:frontier-ub | corpus | 6x13 | [[156,4,ub 38]] | [1, 0] | 39 | 1 | [4, 4] | 36 | [17.4, 11.1] | 15 | 2^4 | RED |
| Z9xZ9:frontier-ub | corpus | 9x9 | d ub 40 | [0, 0] | 81 | 0 | — | — | — | — | no-deck |
| Z6xZ14:frontier-ub | corpus | 6x14 | [[168,4,ub 36]] | [1, 1] | 21 | 2 | [4, 4, 4] | 34 | [17.4, 11.6, 6.0] | 14 | 2^4 | RED |
| Z7xZ12:frontier-ub | corpus | 7x12 | [[168,4,ub 42]] | [0, 2] | 21 | 2 | [4, 4, 4] | 40 | [19.1, 12.0, 6.0] | 17 | 2^4 | RED |
| zoo:bb72 | zoo | 6x6 | [[72,12,exact 6]] | [1, 1] | 9 | 2 | [12, 8, 8] | 4 | [2.7, 2.1, 1.3] | -1 | 2^8 | GREEN |
| zoo:bb90 | zoo | 15x3 | d exact 10 | [0, 0] | 45 | 0 | — | — | — | — | no-deck |
| zoo:bb108 | zoo | 9x6 | [[108,8,exact 10]] | [0, 1] | 27 | 1 | [8, 8] | 8 | [5.4, 4.1] | 1 | 2^8 | GREEN |
| zoo:gross | zoo | 12x6 | [[144,12,exact 12]] | [2, 1] | 9 | 3 | [12, 12, 8, 8] | 10 | [7.0, 5.3, 3.7, 1.8] | 2 | 2^8 | GREEN |
| zoo:bb288 | zoo | 12x12 | [[288,12,exact 18]] | [2, 2] | 9 | 4 | [12, 12, 8, 8, 8] | 16 | [12.5, 9.9, 7.3, 4.4, 1.8] | 5 | 2^8 | GREEN |
| docket:bravyi360 | a35-docket | 30x6 | [[360,12,exact 24]] | [1, 1] | 45 | 2 | [12, 8, 8] | 22 | [17.0, 13.5, 9.8] | 8 | 2^8 | GREEN |
| docket:ibm288Y | a35-docket | 18x8 | [[288,8,exact 20]] | [1, 3] | 9 | 4 | [8, 8, 8] | 18 | [13.7, 10.8, 7.8] | 6 | 2^8 | GREEN |
| docket:gross_xx | a35-docket | 12x6 | [[144,12,exact 12]] | [2, 1] | 9 | 3 | [12, 12, 8] | 10 | [7.0, 5.3, 3.7] | 2 | 2^8 | GREEN |
| docket:bb288_yxx | a35-docket | 12x12 | [[288,12,exact 18]] | [2, 2] | 9 | 4 | [12, 12, 12, 8] | 16 | [12.5, 9.9, 7.1, 4.4] | 5 | 2^8 | GREEN |
| docket:c37x_360420 | a35-docket | 30x6 | [[360,4,exact 20]] | [1, 1] | 45 | 2 | [4, 4, 4] | 18 | [14.7, 11.9, 9.0] | 6 | 2^4 | GREEN |
| docket:e5e50yy_360420 | a35-docket | 15x12 | [[360,4,exact 20]] | [0, 2] | 45 | 2 | [4, 4, 4] | 18 | [14.7, 11.9, 9.0] | 6 | 2^4 | GREEN |
| docket:c37xx_720 | a35-docket | 60x6 | [[720,4,question 20]] | [2, 1] | 45 | 3 | [4, 4, 4, 4] | 18 | [17.4, 14.7, 11.9, 9.0] | 6 | 2^4 | GREEN |
| docket:a8_336 | a35-docket | 12x14 | [[336,12,exact 12]] | [2, 1] | 21 | 3 | [12, 12, 12] | 10 | [9.0, 7.4, 5.7] | 2 | 2^12 | GREEN |
| docket:bravyi756 | a35-docket | 21x18 | [[756,16,exact 34]] | [0, 1] | 189 | 1 | [16, 16] | 32 | [27.7, 22.6] | 13 | 2^16 | RED |
| docket:cover300 | a35-docket | 5x30 | [[300,8,exact 16]] | [0, 1] | 75 | 1 | [8, 8] | 14 | [11.4, 9.2] | 4 | 2^8 | GREEN |
| docket:pair72 | a35-docket | 6x6 | [[72,4,exact 8]] | [1, 1] | 9 | 2 | [4, 4, 4] | 6 | [3.9, 2.9, 2.0] | 0 | 2^4 | GREEN |

Family-level zoo entries with no fixed instance (not priced): qcga (BB family), bb5 (weight-5 family), bicycle, generalized_bicycle.

G5 (the tau-branch ceiling d <= 2 d(mid) per rung) is unknown for corpus rows — per-level distances are not in the corpus; it is part of any eventual closure, not of this pricing.
