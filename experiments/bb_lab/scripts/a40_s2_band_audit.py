#!/usr/bin/env python3
"""A40 S2 / P6 — per-cell cost accounting on certified logicals
(CORRECTED: the ckpt reps are (18,6)-frame vectors; the first cut of
this script mis-read them in the (18,12) frame — caught and fixed;
also note (2,1) itself has NO nontrivial logicals <= 22 at all by the
session-1 sweep, so the audit population lives at the tower's lower
levels).

Audit: the 5,727 + 12 certified nontrivial (18,6)-logicals <= 18 on the
x-cell grid over bb72 (3 cells, d(bb72) = 6), plus the a36 two-gross
witness on its 2x2 grid (archaeology already recorded it: 4 cells,
weights (4,4,8,2), weight 18 < 6*4 — the per-cell lemma's refutation).

Questions: (a) min weight by crossed-cell count t; (b) does ANY
certified logical violate the weaker "total weight >= 6 t"; (c) can a
cell carry < 6 (per-cell-each form)?
"""
from __future__ import annotations

import json
import sys
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


def main():
    validate_banked(LAB / "data")
    print("validate_banked: PASS")
    code = TowerCode("b2", (18, 6), red(A_L, (18, 6)),
                     red(B_L, (18, 6)))
    ng = code.ng

    def xcells(v):
        cells: dict[int, int] = {}
        for i in np.nonzero(v)[0]:
            gx = code.G.from_index(int(i) % ng)[0]
            cells[gx // 6] = cells.get(gx // 6, 0) + 1
        return cells

    stats: dict[tuple, int] = {}
    minw_by_t: dict[int, int] = {}
    viol_6t = []
    min_cell_load = 99
    n = 0
    for line in (DATA / "tdg432" / "ckpt_W22_ntrv1.jsonl").open():
        r = json.loads(line)
        v = np.zeros(code.n, dtype=np.uint8)
        v[r["support"]] = 1
        assert code.is_cycle(v) and not code.is_stab(v)
        cp = xcells(v)
        t = len(cp)
        w = r["w"]
        key = (t, tuple(sorted(cp.values())))
        stats[key] = stats.get(key, 0) + 1
        minw_by_t[t] = min(minw_by_t.get(t, 999), w)
        min_cell_load = min(min_cell_load, min(cp.values()))
        if w < 6 * t:
            viol_6t.append((key, w))
        n += 1
    print(f"(18,6) nontrivial logicals <= 18 audited: {n} "
          f"(x-cell grid over bb72, d(bb72) = 6)")
    print(f"  min weight by crossed-cell count t: "
          f"{dict(sorted(minw_by_t.items()))}")
    print(f"  'total >= 6t' violations: {len(viol_6t)} {viol_6t[:6]}")
    print(f"  min single-cell load seen: {min_cell_load} "
          f"(per-cell-each >= 6 {'FAILS' if min_cell_load < 6 else 'holds'})")
    top = sorted(stats.items(), key=lambda kv: -kv[1])[:12]
    for k, c in top:
        print(f"    {k}: {c}")
    json.dump({"minw_by_t": {str(k): v for k, v in minw_by_t.items()},
               "viol_6t": len(viol_6t),
               "min_cell_load": min_cell_load,
               "hist": {str(k): v for k, v in stats.items()}},
              (DATA / "s2_band_audit.json").open("w"), indent=1)
    print("wrote data/a40/s2_band_audit.json")


if __name__ == "__main__":
    main()
