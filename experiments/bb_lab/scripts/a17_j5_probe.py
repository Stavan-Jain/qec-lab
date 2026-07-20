#!/usr/bin/env python3
"""A17 E17b — (J5) census probe: can dA = 2·dB happen at all under
(iii) on floor-bearing frames?

The (1,5) Frobenius class (T = B + c) matches only if dA = 2·dB (the
w = 5 Theorem-J analog). This probe tests the UNIVERSAL shape-level
claim over whole pools: for every mono-y Sidon 5-set B and every
mono-x Sidon 5-set A on the frame (no pairing, no D2 — D1 ∧ (iii)
alone, matching A16 Theorem J's hypotheses): dA ≠ 2·dB.

Also recorded: |2·dB| collapse statistics (on frames with a Z₂ part,
doubling can collide on dB, auto-killing the match against |dA| = 20)
and the near-miss distribution max |dA ∩ 2·dB| (the separating
invariant for the hand proof).

Usage (from experiments/bb_lab):
    uv run python scripts/a17_j5_probe.py --frames 7x9,6x9,6x10 \
        --out data/a17/e17_j5_probe.json
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from a15_t42_w5_sweep import enumerate_pool, translation_classes
from bb_lab.group import AbelianGroup


def diffs(S, G) -> frozenset:
    return frozenset(G.sub(a, b) for a in S for b in S if a != b)


def probe_frame(ell: int, m: int) -> dict:
    t0 = time.time()
    G = AbelianGroup((ell, m))
    elems = list(G)
    idx = {g: i for i, g in enumerate(elems)}
    poolA = translation_classes(G, enumerate_pool(G, 5, 0, idx))  # mono-x
    poolB = translation_classes(G, enumerate_pool(G, 5, 1, idx))  # mono-y

    dA_index: dict[frozenset, int] = {}
    for supp, _ in poolA:
        dA_index.setdefault(diffs(supp, G), 0)
        dA_index[diffs(supp, G)] += 1

    dbl_sizes: Counter = Counter()
    hits = []
    best_overlap = -1
    best_overlap_count = 0
    for supp, _ in poolB:
        dB = diffs(supp, G)
        two_dB = frozenset(G.add(d, d) for d in dB)
        dbl_sizes[len(two_dB)] += 1
        if len(two_dB) == 20 and two_dB in dA_index:
            hits.append({"B": [list(b) for b in supp],
                         "n_A_shapes": dA_index[two_dB]})
        ov = max((len(two_dB & dA) for dA in dA_index), default=0)
        if ov > best_overlap:
            best_overlap, best_overlap_count = ov, 1
        elif ov == best_overlap:
            best_overlap_count += 1
    return {
        "frame": f"Z{ell}xZ{m}",
        "poolA_classes": len(poolA), "poolB_classes": len(poolB),
        "distinct_dA": len(dA_index),
        "two_dB_size_hist": {str(k): v for k, v in sorted(dbl_sizes.items())},
        "J5_HITS": hits,                      # expect []
        "max_overlap_dA_2dB": best_overlap,   # separating-invariant data
        "n_B_at_max_overlap": best_overlap_count,
        "secs": round(time.time() - t0, 1),
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--frames", type=str, default="7x9,6x9,6x10")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    results = []
    for tok in args.frames.split(","):
        ell, m = (int(x) for x in tok.lower().split("x"))
        r = probe_frame(ell, m)
        results.append(r)
        print(json.dumps(r), flush=True)
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    Path(args.out).write_text(json.dumps(results, indent=1))
    n_hits = sum(len(r["J5_HITS"]) for r in results)
    print(f"\nJ5 universal-claim violations: {n_hits} (expect 0)")


if __name__ == "__main__":
    main()
