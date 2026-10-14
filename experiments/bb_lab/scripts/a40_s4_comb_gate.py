#!/usr/bin/env python3
"""A40 S4 gate 0 — the comb refutation + row-recurrence conventions.

Claim under test (falsify-first, against the S3 note §8.2 formulation):
the "defect graph min-mean-cycle" formulation is VACUOUS as stated,
because TRIVIAL (stabilizer) y-spanning cycles exist at mean cost far
below the target rate 2/row.  The witness family: "combs" = sums of
single Z-checks tiled with y-period 6 — weight m on an m-row torus
(rate 1), all cyclic y-gaps <= 1.  If the combs check out, §8.2 needs
the S4 reformulation (slab-amortized costs + light core + local-
reduction pruning); the compression lemma was not the missing piece.

Gates:
  C1: combs at (2,0), (2,1), (3,0) members — stabilizer, y-spanning,
      rate 1; period-8 comb at the non-member (12,16) frame — rate 3/4.
  C2: the row recurrence E_j (conv: circulant M[g,h] = P(g-h)):
        E_j : (1+x^-1) v1[j] + x v1[j-3] + v2[j] + v2[j+1]
              + x^-3 v2[j-1] = 0
      identity on random kernel cycles, and the v2-DETERMINISM
      reconstruction (block-2 rows forced by history: the system is a
      convolutional code in y with block-1 rows free).
  C3: the slab-telescope identity sum_j W_j = 4|v| (exact, cyclic).
  C4: small-tooth boundary sweep — tile d(t·delta) at y-period g for
      tooth polynomials t with support in a 2x2 window, g in {4,..,9}:
      the empirical trivial-rate landscape the pruning catalog must
      kill (min rate, y-spanning only).
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np

LAB = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(LAB / "src"))

from bb_lab.tower import TowerCode, validate_banked  # noqa: E402

DATA = LAB / "data" / "a40"
A_L = [(0, 0), (0, 1), (3, -1)]
B_L = [(0, 0), (1, 0), (-1, -3)]


def red(supp, lm):
    return frozenset((e[0] % lm[0], e[1] % lm[1]) for e in supp)


def code_at(lm, name=None):
    return TowerCode(name or f"tdg{lm}", lm, red(A_L, lm), red(B_L, lm))


def rows_of(code, v):
    """Per-row weights and occupied-row set of a qubit vector."""
    l, m = code.G.orders
    w = np.zeros(m, dtype=np.int64)
    for i in np.nonzero(v)[0]:
        g = code.G.from_index(int(i) % code.ng)
        w[g[1]] += 1
    return w


def max_cyclic_gap(wrow):
    m = len(wrow)
    occ = [j for j in range(m) if wrow[j]]
    if not occ:
        return m
    if len(occ) == m:
        return 0
    gaps = [(occ[(i + 1) % len(occ)] - occ[i]) % m - 1
            for i in range(len(occ))]
    return max(gaps)


def stab_from_checks(code, checks):
    """Sum of H_X rows at the given (x,y) check positions."""
    v = np.zeros(code.n, dtype=np.uint8)
    for g in checks:
        v ^= code.HX[code.G.index((g[0] % code.G.orders[0],
                                   g[1] % code.G.orders[1]))]
    return v


def row_polys(code, v):
    """v as (v1_rows, v2_rows): lists of length-m arrays of length-l
    F2 row contents."""
    l, m = code.G.orders
    out = [np.zeros((m, l), dtype=np.uint8) for _ in range(2)]
    for i in np.nonzero(v)[0]:
        blk, gi = divmod(int(i), code.ng)
        gx, gy = code.G.from_index(gi)
        out[blk][gy, gx] = 1
    return out


def xs(row, s, l):
    """Multiply a length-l F2 row by x^s (cyclic shift by +s)."""
    return np.roll(row, s % l)


def main():
    t0 = time.time()
    validate_banked(LAB / "data")
    print("validate_banked: PASS")
    out = {}

    # ---- C1: the combs -------------------------------------------------
    comb_rows = []
    for lm, period in [((12, 12), 6), ((18, 12), 6), ((18, 18), 6),
                       ((12, 16), 8)]:
        c = code_at(lm)
        l, m = lm
        assert m % period == 0
        v = stab_from_checks(c, [(0, period * j) for j in range(m // period)])
        w = int(v.sum())
        wrow = rows_of(c, v)
        gap = max_cyclic_gap(wrow)
        assert c.is_cycle(v) and c.is_stab(v), (lm, "comb not a stabilizer")
        assert gap <= 3, (lm, "comb not y-spanning", gap)
        # slab weights (4-row windows) and the telescope
        W = np.array([wrow[(j - 3) % m] + wrow[(j - 2) % m]
                      + wrow[(j - 1) % m] + wrow[j] for j in range(m)])
        assert W.sum() == 4 * w
        rate = w / m
        comb_rows.append(dict(lm=lm, period=period, weight=w, rate=rate,
                              max_y_gap=int(gap), max_slab=int(W.max())))
        print(f"C1 comb {lm} period {period}: weight {w} = rate "
              f"{rate:g}/row (target rate 2), max y-gap {gap}, "
              f"max slab {W.max()} -> TRIVIAL y-spanning cheap cycle")
    assert any(r["rate"] <= 1.0 for r in comb_rows)
    assert any(r["rate"] < 1.0 for r in comb_rows)   # (12,16) 3/4
    out["C1_combs"] = comb_rows

    # ---- C2: recurrence + determinism ---------------------------------
    c = code_at((18, 12))
    l, m = c.G.orders
    rng = np.random.default_rng(40)
    n_checked = 0
    for _ in range(20):
        v = c.random_cycle(rng)
        v1, v2 = row_polys(c, v)
        # E_j: (1+x^-1) v1[j] + x v1[j-3] + v2[j] + v2[j+1] + x^-3 v2[j-1]
        for j in range(m):
            e = ((v1[j] ^ xs(v1[j], -1, l)) ^ xs(v1[(j - 3) % m], 1, l)
                 ^ v2[j] ^ v2[(j + 1) % m] ^ xs(v2[(j - 1) % m], -3, l))
            assert not e.any(), f"E_{j} fails"
        # determinism: reconstruct v2[2..] from v1[.] + v2[0], v2[1]
        w2 = np.zeros_like(v2)
        w2[0], w2[1] = v2[0], v2[1]
        for j in range(1, m - 1):
            w2[j + 1] = (v1[j] ^ xs(v1[j], -1, l)
                         ^ xs(v1[(j - 3) % m], 1, l)
                         ^ w2[j] ^ xs(w2[(j - 1) % m], -3, l))
        assert (w2 == v2).all(), "v2-determinism reconstruction fails"
        n_checked += 1
    print(f"C2 recurrence + v2-determinism: PASS on {n_checked} random "
          f"cycles at (18,12)")
    out["C2_recurrence"] = dict(frame=[18, 12], n_random_cycles=n_checked)

    # ---- C3 is asserted inside C1 (telescope) -------------------------
    out["C3_telescope"] = "asserted in C1 (sum_j W_j == 4|v|)"

    # ---- C4: small-tooth tiled-boundary sweep -------------------------
    # tooth s = t(x,y)*delta, supp(t) subset {0,1}x{0,1}; tile at period g.
    best = []
    for lm in [(12, 12), (12, 24), (18, 12)]:
        c = code_at(lm)
        l, m = lm
        rows_best = None
        for tmask in range(1, 16):
            t = [(dx, dy) for k, (dx, dy) in enumerate(
                [(0, 0), (1, 0), (0, 1), (1, 1)]) if tmask >> k & 1]
            for g in range(4, 10):
                if m % g:
                    continue
                checks = [(dx, dy + g * j) for j in range(m // g)
                          for (dx, dy) in t]
                v = stab_from_checks(c, checks)
                if not v.any():
                    continue
                wrow = rows_of(c, v)
                if max_cyclic_gap(wrow) > 3:
                    continue
                rate = v.sum() / m
                if rows_best is None or rate < rows_best["rate"]:
                    rows_best = dict(lm=lm, tooth=t, period=g,
                                     weight=int(v.sum()),
                                     rate=float(rate))
        best.append(rows_best)
        print(f"C4 {lm}: cheapest tiled trivial y-spanning boundary: "
              f"rate {rows_best['rate']:g}/row  (tooth {rows_best['tooth']}, "
              f"period {rows_best['period']}, weight {rows_best['weight']})")
    out["C4_tooth_sweep"] = best

    out["wall_s"] = round(time.time() - t0, 1)
    DATA.mkdir(parents=True, exist_ok=True)
    (DATA / "s4_comb_gate.json").write_text(json.dumps(out, indent=1))
    print(f"\nVERDICT: trivial y-spanning cycles at rate <= "
          f"{min(r['rate'] for r in comb_rows):g}/row exist (stabilizer "
          f"combs) -> the bare min-mean-cycle formulation of §8.2 is "
          f"VACUOUS; the S4 reformulation is required.")
    print(f"wrote {DATA/'s4_comb_gate.json'} ({out['wall_s']} s)")


if __name__ == "__main__":
    main()
