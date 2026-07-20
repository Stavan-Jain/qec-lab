#!/usr/bin/env python3
"""A17 E18c — (2,6) end-to-end census: dense 6-sets vs σ_L forms.

The (2,6) reduction (doc §8 forward note + profile counting): C2 caps
Σ_j |(B+r_j) ∩ σ| ≤ 12, which with (†) forces a(u_R) ≥ 10 — u_R is a
6-vertex near-clique of Cay(G, dB) (≥ 10 of 15 differences in dB),
n₅ = n₆-profiles pinned per size — and no full ±B-translate fits
inside u_R survivably (fifth-translate / squares argument).

Enumeration soundness: a ≥ 10 ⟹ Σ complement-degrees = 2(15−a) ≤ 10
⟹ some vertex has complement-degree ≤ 1, i.e. is dB-adjacent to ≥ 4
of the other 5. Translate it to 0: T = {0} ∪ S₄ ∪ {v}, S₄ ⊆ dB
(4 of 0's neighbors), v ∈ G arbitrary. All viable T's arise this
way; duplicates removed by translation-canon.

Censused per member: every candidate's (a, |σ|) profile and the full
match test canon(σ_R) ∈ {canon(A Δ (A+δ))}. Expect 0 matches.

Usage: uv run python scripts/a17_e18c_26_census.py \
    --members data/a17/members_7x9.jsonl,data/a17/members_6x9_6x10.jsonl \
    --out data/a17/e18c_26_census.json
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

from a17_e18_k4_census import canon_translate, diffs
from bb_lab.group import AbelianGroup
from bb_lab.poly import Poly


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--members", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    t0 = time.time()

    out = {"members": 0, "candidates": 0, "matches": 0,
           "match_examples": []}
    prof_hist: Counter = Counter()
    for path in args.members.split(","):
        with open(path) as f:
            rows = [json.loads(l) for l in f if '"A"' in l]
        for r in rows:
            G = AbelianGroup(tuple(r["frame"]))
            A = frozenset(Poly.from_string(r["A"], G).support)
            B = sorted(Poly.from_string(r["B"], G).support)
            Bset = frozenset(B)
            dBs = diffs(Bset, G)
            dB = sorted(dBs)
            zero = tuple(0 for _ in range(G.rank))
            elems = list(G)
            formsL = set()
            for d in elems:
                if d == zero:
                    continue
                sig = A.symmetric_difference(
                    frozenset(G.add(a, d) for a in A))
                if sig:
                    formsL.add(canon_translate(sig, G))
            seen: set = set()
            for s4 in combinations(dB, 4):
                p4 = sum(1 for u, v in combinations(s4, 2)
                         if G.sub(v, u) in dBs)
                need = 6 - p4
                if need > 5:
                    continue
                base = (zero,) + s4
                bset = frozenset(base)
                for v6 in elems:
                    if v6 in bset:
                        continue
                    deg = sum(1 for u in base if G.sub(v6, u) in dBs)
                    if deg < need:
                        continue
                    T = bset | {v6}
                    c = canon_translate(T, G)
                    if c in seen:
                        continue
                    seen.add(c)
                    out["candidates"] += 1
                    aT = 4 + p4 + deg
                    mult = Counter(G.add(b, t) for t in T for b in B)
                    sigma = frozenset(
                        y for y, cnt in mult.items() if cnt & 1)
                    prof_hist[(aT, len(sigma))] += 1
                    if len(sigma) not in (8, 10):
                        continue
                    if canon_translate(sigma, G) in formsL:
                        out["matches"] += 1
                        if len(out["match_examples"]) < 5:
                            out["match_examples"].append({
                                "frame": r["frame"], "A": r["A"],
                                "B": r["B"],
                                "u_R": [list(t) for t in sorted(T)],
                                "a": aT, "size": len(sigma)})
            out["members"] += 1
    out["profile_hist"] = {f"a={a},sz={s}": c for (a, s), c
                           in sorted(prof_hist.items())}
    out["secs"] = round(time.time() - t0, 1)
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    Path(args.out).write_text(json.dumps(out, indent=1))
    print(json.dumps(out))
    print(f"\nE18c (2,6): {'CLEAN' if out['matches'] == 0 else 'MATCHES!'}"
          f" ({out['candidates']} candidates, {out['members']} members)")


if __name__ == "__main__":
    main()
