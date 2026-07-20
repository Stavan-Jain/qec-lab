#!/usr/bin/env python3
"""A17 E18 — even-split census: (2,2) atoms + the translate-class
coincidence, over full pools (D1 ∧ (iii), no pairing).

The w = 5 (2,2) analysis (A17 doc §7):
  atom identity   mult_{d(σ_L)}(δ_L) = 5 + [2δ_L ∈ dA] (+5 if 2δ_L=0),
  off-branch cap  mult_{d(σ_R)}(δ_L) ≤ 2·[δ_L∈dB] + [δ_L−δ_R∈dB] +
                  [δ_L+δ_R∈dB] ≤ 4 < 5 when δ_L ≠ ±δ_R,
  ⟹ any (2,2) match has δ_L = ±δ_R (branch 2b is VACUOUS at w = 5),
  ⟹ dB = dA + δ = dA − δ (branch 2a, counting over dA).

Censused here:
  E1: the atom identity, verified over sampled (A, δ) pairs;
  E2: the translate-class coincidence — is dB EVER a translate of dA
      across full mono-x / mono-y pools? (the 2a precondition; expect
      0 — this is the (2,2) analog of the J5 universal probe);
  E3: direct (2,2) match census on members (σ_L ~ σ_R + t): expect 0
      (members have d ≥ 10 > 4).

Usage: uv run python scripts/a17_e18_even_census.py \
    --frames 7x9,6x9,6x10 --members data/a17/members_7x9.jsonl \
    --out data/a17/e18_census.json
"""

from __future__ import annotations

import argparse
import json
import random
import sys
import time
from collections import Counter
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


def sym_diff_translate(S: frozenset, delta, G) -> frozenset:
    return S.symmetric_difference(frozenset(G.add(x, delta) for x in S))


def diff_multiset(S, G) -> Counter:
    c: Counter = Counter()
    for a in S:
        for b in S:
            if a != b:
                c[G.sub(a, b)] += 1
    return c


def census_frame(ell: int, m: int, rng: random.Random,
                 atom_samples: int = 400) -> dict:
    t0 = time.time()
    G = AbelianGroup((ell, m))
    elems = list(G)
    idx = {g: i for i, g in enumerate(elems)}
    zero = tuple(0 for _ in range(G.rank))
    nonzero = [g for g in elems if g != zero]
    poolA = translation_classes(G, enumerate_pool(G, 5, 0, idx))
    poolB = translation_classes(G, enumerate_pool(G, 5, 1, idx))

    # E1: atom identity on random (A, δ)
    e1_bad = 0
    for _ in range(atom_samples):
        supp, _ = rng.choice(poolA)
        dA = diffs(supp, G)
        delta = rng.choice(nonzero)
        if delta in dA:
            continue  # size-8 stratum has its own formula; sample on
        sigma = sym_diff_translate(frozenset(supp), delta, G)
        dm = diff_multiset(sigma, G)
        two = G.add(delta, delta)
        expect = 5 + (1 if two in dA else 0) + (5 if two == zero else 0)
        if dm.get(delta, 0) != expect:
            e1_bad += 1

    # E2: translate-class coincidence dB ~ dA + δ
    classA = {canon_translate(diffs(supp, G), G) for supp, _ in poolA}
    e2_hits = []
    for supp, _ in poolB:
        if canon_translate(diffs(supp, G), G) in classA:
            e2_hits.append([list(b) for b in supp])

    return {
        "frame": f"Z{ell}xZ{m}",
        "poolA": len(poolA), "poolB": len(poolB),
        "E1_atom_violations": e1_bad, "E1_samples": atom_samples,
        "E2_translate_class_hits": e2_hits,
        "secs": round(time.time() - t0, 1),
    }


def match_census_members(path: str) -> dict:
    """E3: direct (2,2) match search per member: σ_L ~ σ_R + t."""
    out = {"members": 0, "matches": 0}
    with open(path) as f:
        rows = [json.loads(l) for l in f if '"A"' in l]
    for r in rows:
        G = AbelianGroup(tuple(r["frame"]))
        A = frozenset(Poly.from_string(r["A"], G).support)
        B = frozenset(Poly.from_string(r["B"], G).support)
        elems = [g for g in G]
        zero = tuple(0 for _ in range(G.rank))
        nonzero = [g for g in elems if g != zero]
        formsR = {}
        for dR in nonzero:
            formsR.setdefault(
                canon_translate(sym_diff_translate(B, dR, G), G), dR)
        for dL in nonzero:
            f_ = canon_translate(sym_diff_translate(A, dL, G), G)
            if f_ in formsR:
                out["matches"] += 1
        out["members"] += 1
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--frames", default="7x9,6x9,6x10")
    ap.add_argument("--members", default=None)
    ap.add_argument("--out", required=True)
    ap.add_argument("--seed", type=int, default=20260719)
    args = ap.parse_args()
    rng = random.Random(args.seed)

    results = []
    for tok in args.frames.split(","):
        ell, m = (int(x) for x in tok.lower().split("x"))
        r = census_frame(ell, m, rng)
        results.append(r)
        print(json.dumps(r), flush=True)
    e3 = None
    if args.members:
        e3 = match_census_members(args.members)
        print(json.dumps({"E3": e3}), flush=True)
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    Path(args.out).write_text(json.dumps(
        {"frames": results, "E3": e3}, indent=1))
    bad = sum(r["E1_atom_violations"] for r in results) \
        + sum(len(r["E2_translate_class_hits"]) for r in results) \
        + ((e3 or {}).get("matches", 0))
    print(f"\nE18 census: {'CLEAN' if bad == 0 else f'{bad} ANOMALIES'}")


if __name__ == "__main__":
    main()
