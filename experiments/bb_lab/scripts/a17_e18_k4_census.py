#!/usr/bin/env python3
"""A17 E18 — (2,4) falsify-first: K₄ census in Cay(G, dB) + full
weight-4 end-to-end match sweep on members.

Theorem E18.2 (A17 doc §8) reduces any (2,4) match to: u_R a K₄ in
Cay(G, dB) with all six pairwise collision cells distinct and off σ,
|σ| = 8, δ_L ∈ dA, and a pinned 4×2 incidence grid. The endgames kill
u_R ⊆ B + g (fifth-translate cap) and u_R ⊆ −B + g (an m = 4 cell).
What remains is **Lemma K**: every 4-set with all six differences in
dB lies in a translate of B or −B.

Censused here:
  K1: exhaustive K₄ enumeration per B-class over full mono-y pools,
      classified sub-B / sub-(−B) / both / EXOTIC (falsify Lemma K);
      plus the c* = |B² ∩ (B + b*)| histogram (c* = 1 is the size-
      viable sub-B stratum — killed by the cap — and c* ≥ 2 is
      size-dead).
  K2: every weight-4 u_R on every member, end to end: σ_R =
      odd-support of B·u_R matched (up to translation) against every
      σ_L = A Δ (A + δ). Expect 0 — this verifies Theorem E18.2 +
      endgames + Lemma K jointly at the conclusion level.

Usage: uv run python scripts/a17_e18_k4_census.py \
    --frames 7x9,6x9,6x10 \
    --members data/a17/members_7x9.jsonl,data/a17/members_6x9_6x10.jsonl \
    --out data/a17/e18_k4_census.json
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from collections import Counter
from itertools import combinations
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from a15_t42_w5_sweep import enumerate_pool, translation_classes
from bb_lab.group import AbelianGroup
from bb_lab.poly import Poly


def diffs(S, G) -> frozenset:
    return frozenset(G.sub(a, b) for a in S for b in S if a != b)


def canon_translate(S: frozenset, G) -> tuple:
    """Lexicographic-min translate normal form of a finite subset."""
    best = None
    for s0 in S:
        form = tuple(sorted(G.sub(x, s0) for x in S))
        if best is None or form < best:
            best = form
    return best


def k4_classes(B: frozenset, G) -> dict:
    """All K₄s of Cay(G, dB) up to translation, classified."""
    dB = sorted(diffs(B, G))
    dBs = frozenset(dB)
    Bs = frozenset(B)
    zero = tuple(0 for _ in range(G.rank))
    seen: dict = {}
    for i, d1 in enumerate(dB):
        n1 = [d for d in dB[i + 1:] if G.sub(d, d1) in dBs]
        for a_i in range(len(n1)):
            for b_i in range(a_i + 1, len(n1)):
                d2, d3 = n1[a_i], n1[b_i]
                if G.sub(d3, d2) not in dBs:
                    continue
                T = frozenset((zero, d1, d2, d3))
                c = canon_translate(T, G)
                if c in seen:
                    continue
                sub_b = any(all(G.add(t, b) in Bs for t in T) for b in B)
                sub_nb = any(all(G.sub(b, t) in Bs for t in T) for b in B)
                # six pairwise collision cells (unique by Sidon)
                cells = set()
                for t, tp in combinations(sorted(T), 2):
                    delta = G.sub(tp, t)
                    for b in B:
                        if G.sub(b, delta) in Bs:
                            cells.add(G.add(b, t))
                            break
                seen[c] = (sub_b, sub_nb, len(cells) == 6)
    return seen


def census_frame(ell: int, m: int) -> dict:
    t0 = time.time()
    G = AbelianGroup((ell, m))
    elems = list(G)
    idx = {g: i for i, g in enumerate(elems)}
    poolB = translation_classes(G, enumerate_pool(G, 5, 1, idx))
    counts = Counter()
    k4_per_class = Counter()
    cstar_hist = Counter()
    exotics = []
    for supp, _ in poolB:
        B = frozenset(supp)
        Bs = B
        cls = k4_classes(B, G)
        k4_per_class[len(cls)] += 1
        for c, (sb, snb, dist6) in cls.items():
            key = ("both" if sb and snb else "sub_B" if sb
                   else "sub_negB" if snb else "EXOTIC")
            counts[key] += 1
            if key == "EXOTIC" and len(exotics) < 5:
                exotics.append({"B": [list(b) for b in sorted(B)],
                                "T": [list(t) for t in c],
                                "cells_distinct": dist6})
        for bst in B:
            cstar = sum(
                1 for b in B
                if G.sub(G.add(b, b), bst) in Bs)
            cstar_hist[cstar] += 1
    return {
        "frame": f"Z{ell}xZ{m}", "poolB": len(poolB),
        "k4_per_class_hist": dict(sorted(k4_per_class.items())),
        "class_counts": dict(counts),
        "cstar_hist": dict(sorted(cstar_hist.items())),
        "exotic_examples": exotics,
        "secs": round(time.time() - t0, 1),
    }


def member_sweep(paths: list[str]) -> dict:
    """K2: every weight-4 u_R against every weight-2 u_L, per member."""
    t0 = time.time()
    out = {"members": 0, "triples": 0, "ann_hits": 0,
           "near_size8": 0, "near_size10": 0, "k4_near": 0,
           "matches": 0, "match_examples": []}
    for path in paths:
        with open(path) as f:
            rows = [json.loads(l) for l in f if '"A"' in l]
        for r in rows:
            G = AbelianGroup(tuple(r["frame"]))
            A = frozenset(Poly.from_string(r["A"], G).support)
            B = sorted(Poly.from_string(r["B"], G).support)
            dBs = diffs(frozenset(B), G)
            zero = tuple(0 for _ in range(G.rank))
            nonzero = [g for g in G if g != zero]
            formsL = set()
            for d in nonzero:
                sig = A.symmetric_difference(
                    frozenset(G.add(a, d) for a in A))
                if sig:
                    formsL.add(canon_translate(sig, G))
            for tri in combinations(nonzero, 3):
                out["triples"] += 1
                T = (zero,) + tri
                mult = Counter(G.add(b, t) for t in T for b in B)
                sigma = frozenset(
                    y for y, c in mult.items() if c & 1)
                if not sigma:
                    out["ann_hits"] += 1
                    continue
                if len(sigma) not in (8, 10):
                    continue
                out["near_size8" if len(sigma) == 8
                    else "near_size10"] += 1
                aT = sum(1 for u, v in combinations(T, 2)
                         if G.sub(v, u) in dBs)
                if aT == 6:
                    out["k4_near"] += 1
                if canon_translate(sigma, G) in formsL:
                    out["matches"] += 1
                    if len(out["match_examples"]) < 5:
                        out["match_examples"].append({
                            "frame": r["frame"], "A": r["A"],
                            "B": r["B"], "u_R": [list(t) for t in T],
                            "a": aT, "size": len(sigma)})
            out["members"] += 1
    out["secs"] = round(time.time() - t0, 1)
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--frames", default="7x9,6x9,6x10")
    ap.add_argument("--members", default=None)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    results = []
    for tok in args.frames.split(","):
        ell, m = (int(x) for x in tok.lower().split("x"))
        r = census_frame(ell, m)
        results.append(r)
        print(json.dumps(r), flush=True)
    k2 = None
    if args.members:
        k2 = member_sweep(args.members.split(","))
        print(json.dumps({"K2": k2}), flush=True)
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    Path(args.out).write_text(json.dumps(
        {"frames": results, "K2": k2}, indent=1))
    n_exotic = sum(r["class_counts"].get("EXOTIC", 0) for r in results)
    n_match = (k2 or {}).get("matches", 0)
    print(f"\nK1: {'EXOTIC-FREE' if n_exotic == 0 else f'{n_exotic} EXOTICS'}"
          f" | K2: {'CLEAN' if n_match == 0 else f'{n_match} MATCHES'}")


if __name__ == "__main__":
    main()
