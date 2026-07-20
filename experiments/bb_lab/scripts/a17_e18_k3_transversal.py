#!/usr/bin/env python3
"""A17 E18 — K3: the dA-K₄ transversal census on (2,4) candidates.

Lemma K (all K₄ ⊆ ±B + g) is REFUTED (K1: exotic K₄s abound), so the
(2,4) kill moves to the σ-form level. A match σ_L = σ_R forces (doc
§8): σ splits as (A \ {a₂}) ⊔ ((A \ {a₁}) + δ) — TWO complementary
grid transversals, each a K₄ of Cay(G, dA) — and N_A(σ) ≥ 30 of the
48 cross-translate ordered pairs.

Censused, over every member and every a = 6 K₄ class of its B with
the viable multiplicity profile (n₂ = 6, |σ| = 8):
  N_A(σ) histogram (demand: ≥ 30), and
  #grids admitting ANY transversal that is a dA-K₄ (demand: a
  complementary PAIR of them).
Expect: max N_A well under 30 and zero transversal grids.

Usage: uv run python scripts/a17_e18_k3_transversal.py \
    --members data/a17/members_7x9.jsonl,data/a17/members_6x9_6x10.jsonl \
    --out data/a17/e18_k3_transversal.json
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from collections import Counter
from itertools import combinations, product
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from a17_e18_k4_census import diffs, k4_classes
from bb_lab.group import AbelianGroup
from bb_lab.poly import Poly


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--members", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    t0 = time.time()

    na_hist: Counter = Counter()
    out = {"members": 0, "k4_classes": 0, "viable_grids": 0,
           "transversal_grids": 0, "transversal_examples": []}
    for path in args.members.split(","):
        with open(path) as f:
            rows = [json.loads(l) for l in f if '"A"' in l]
        for r in rows:
            G = AbelianGroup(tuple(r["frame"]))
            A = frozenset(Poly.from_string(r["A"], G).support)
            B = frozenset(Poly.from_string(r["B"], G).support)
            dAs = diffs(A, G)
            out["members"] += 1
            for T in k4_classes(B, G):
                out["k4_classes"] += 1
                mult = Counter(G.add(b, t) for t in T for b in B)
                sigma = [y for y, c in mult.items() if c & 1]
                if len(sigma) != 8:
                    continue
                out["viable_grids"] += 1
                na = sum(1 for y, yp in product(sigma, sigma)
                         if y != yp and G.sub(y, yp) in dAs)
                na_hist[na] += 1
                # the 4x2 grid: two sigma-cells per translate
                grid = []
                for t in T:
                    grid.append([y for y in sigma
                                 if G.sub(y, t) in B])
                if any(len(gcol) != 2 for gcol in grid):
                    continue  # profile violation (none expected)
                for choice in product(*grid):
                    if all(G.sub(v, u) in dAs
                           for u, v in combinations(choice, 2)):
                        out["transversal_grids"] += 1
                        if len(out["transversal_examples"]) < 5:
                            out["transversal_examples"].append({
                                "frame": r["frame"], "A": r["A"],
                                "B": r["B"],
                                "T": [list(t) for t in T],
                                "cells": [list(c) for c in choice]})
                        break
    out["na_hist"] = dict(sorted(na_hist.items()))
    out["na_max"] = max(na_hist) if na_hist else None
    out["secs"] = round(time.time() - t0, 1)
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    Path(args.out).write_text(json.dumps(out, indent=1))
    print(json.dumps(out))
    ok = (out["transversal_grids"] == 0
          and (out["na_max"] is None or out["na_max"] < 30))
    print(f"\nK3: {'CLEAN — cap holds' if ok else 'RESIDUALS'}"
          f" (na_max={out['na_max']}, transversals="
          f"{out['transversal_grids']}/{out['viable_grids']})")


if __name__ == "__main__":
    main()
