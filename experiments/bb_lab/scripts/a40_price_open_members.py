#!/usr/bin/env python3
"""A40 P2 — price the open tour-de-gross members through the A35/A38
screen (`bb_lab.tower.screen_tower`).  COST verdicts only, never distance
claims.

Members and their maximal iterated-Z2 towers (P1 inventory):

  (2,1) [[432,12]]  d_conj 24: (18,12) -x-> (9,12) -y-> (9,6) -y-> (9,3)
                    and the (18,6)-mid ordering; bottom n = 54.
  (3,0) [[648,12]]  d_conj 30: (18,18) -x-> (9,18) -y-> (9,9); bottom 162.
  (3,1) [[864,12]]  d_conj 36: (24,18) x,x,x,y -> (3,9); bottom 54.
  (4,0) [[1152,12]] d_conj 42: (24,24) full depth 6 -> (3,3); bottom 18.
  (4,1) [[1440,12]] d_conj 48: (30,24) y,y,y,x -> (15,3); bottom 90.

W per member: the decisive floor question for the conjectured d is
W = d_conj - 2 (all cycles are even: |A|,|B| odd => parity scope).
Smaller W rows price the partial questions (d >= W+2).

S2 refinements carried: the cap gate is n-blind (we also report
C(n_bottom - mu, ceil(cap/2)) per fiber as the honest MITM cost) and the
kernel-shift lane converts deep caps on LIGHT shadows into level-(r-1)
census windows — so cap-RED rows with light-shadow structure may still be
executable; the census/window population is the storage axis (A39).
"""
from __future__ import annotations

import json
import math
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from bb_lab.tower import (  # noqa: E402
    screen_tower, validate_banked, census_nodes,
)

DATA = ROOT / "data" / "a40"
DATA.mkdir(parents=True, exist_ok=True)


def laur(lm):
    A = frozenset(((e[0] % lm[0]), (e[1] % lm[1]))
                  for e in [(0, 0), (0, 1), (3, -1)])
    B = frozenset(((e[0] % lm[0]), (e[1] % lm[1]))
                  for e in [(0, 0), (1, 0), (-1, -3)])
    return A, B


SPECS = []
for name, lm, folds, d_conj, W_list, fib in [
    ("tdg432_xyy", (18, 12), [(0, 9), (1, 6), (1, 3)], 24,
     [10, 14, 16, 18, 22], True),
    ("tdg432_yxy", (18, 12), [(1, 6), (0, 9), (1, 3)], 24,
     [16, 22], True),
    ("tdg648_xy", (18, 18), [(0, 9), (1, 9)], 30, [16, 22, 28], False),
    ("tdg864_xxxy", (24, 18), [(0, 12), (0, 6), (0, 3), (1, 9)], 36,
     [16, 22, 28, 34], False),
    ("tdg1152_full", (24, 24), [(0, 12), (1, 12), (0, 6), (1, 6),
                                (0, 3), (1, 3)], 42,
     [22, 28, 34, 40], False),
    ("tdg1440_yyyx", (30, 24), [(1, 12), (1, 6), (1, 3), (0, 15)], 48,
     [22, 34, 46], False),
]:
    A, B = laur(lm)
    SPECS.append(dict(name=name, top=(lm, A, B), folds=folds,
                      d_top=None, W_eff=d_conj - 2, W_list=W_list,
                      fibers=fib, tag="A40",
                      notes=f"tour-de-gross member, d_conj={d_conj}"))


def main():
    t0 = time.time()
    print("== gate: validate_banked ==")
    validate_banked(ROOT / "data")
    print("   PASS")
    rng = np.random.default_rng(20260825)
    results = []
    for spec in SPECS:
        res = screen_tower(spec, rng=rng, log=lambda s: print(s))
        # honest MITM per-fiber cost at the top rung (S2: cap gate is
        # n-blind): C(n_mid - mu, ceil(cap/2)) for the W_list rows
        n_mid = res["levels"][1]["n"]
        mu = res["costs"][0]["mu"] if res["costs"] else 6
        for c in res["costs"]:
            half = (c["cap_max"] + 1) // 2
            c["mitm_half_subsets_top_rung"] = float(
                math.comb(max(n_mid - c["mu"], 0), max(half, 0)))
        results.append(res)
        print()
    out = {"when": time.strftime("%Y-%m-%d %H:%M:%S"),
           "screens": results,
           "wall_s": round(time.time() - t0, 1)}
    (DATA / "pricing.json").write_text(json.dumps(out, indent=1))
    print(f"wrote {DATA/'pricing.json'} ({out['wall_s']} s)")

    # compact verdict table
    print("\n== pricing table (COST verdicts, never distance claims) ==")
    for res in results:
        ks = "/".join(str(lv["k"]) for lv in res["levels"])
        ns = "/".join(str(lv["n"]) for lv in res["levels"])
        print(f"{res['name']:>14}  n {ns}  k {ks}")
        for c in res["costs"]:
            print(f"    W={c['W']:>2} (d>= {c['W']+2}): nodes/level "
                  f"{c['log10_nodes_per_level']} cap {c['cap_max']} "
                  f"-> {c['verdict']}  "
                  f"[MITM half-subsets top rung "
                  f"{c['mitm_half_subsets_top_rung']:.1e}]")


if __name__ == "__main__":
    main()
