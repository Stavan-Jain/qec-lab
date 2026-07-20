#!/usr/bin/env python3
"""A17 E18d — (4,4) end-to-end census: every weight-4 u_L form vs
every weight-4 u_R form, per member, up to translation.

Opening counts (doc §10): with |u_L| = |u_R| = 4, cells of σ carry
odd A-multiplicity (1 or 3) and odd B-multiplicity; two cells of one
B-translate sharing an A-translate die by D2, so the A-translate
families over (B+r_j) ∩ σ are pairwise disjoint ⟹
Σ_j |(B+r_j) ∩ σ| ≤ 16 ⟹ a_B(u_R) ≥ 2, and symmetrically
a_A(u_L) ≥ 2 — much weaker forcing than (2,4)'s equality chain: the
(4,4) analytic layer is genuinely new. Falsify-first: this census
checks the conclusion (no matches) exhaustively on members.

Usage: uv run python scripts/a17_e18d_44_census.py \
    --members data/a17/members_7x9.jsonl,data/a17/members_6x9_6x10.jsonl \
    --out data/a17/e18d_44_census.json
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

from a17_e18_k4_census import canon_translate
from bb_lab.group import AbelianGroup
from bb_lab.poly import Poly


def forms4(S: list, G, collect) -> dict:
    """canon(supp(S·u)) for every weight-4 u ∋ 0, keyed by canon."""
    zero = tuple(0 for _ in range(G.rank))
    nonzero = [g for g in G if g != zero]
    out: dict = {}
    n_ann = 0
    for tri in combinations(nonzero, 3):
        u = (zero,) + tri
        mult = Counter(G.add(s, t) for t in u for s in S)
        sigma = frozenset(y for y, c in mult.items() if c & 1)
        if not sigma:
            n_ann += 1
            continue
        if collect:
            out.setdefault(canon_translate(sigma, G), u)
    return {"forms": out, "ann": n_ann}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--members", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    t0 = time.time()

    out = {"members": 0, "pairs_tested": 0, "ann_hits": 0,
           "matches": 0, "match_examples": []}
    for path in args.members.split(","):
        with open(path) as f:
            rows = [json.loads(l) for l in f if '"A"' in l]
        for r in rows:
            G = AbelianGroup(tuple(r["frame"]))
            A = sorted(Poly.from_string(r["A"], G).support)
            B = sorted(Poly.from_string(r["B"], G).support)
            fl = forms4(A, G, collect=True)
            formsL = fl["forms"]
            out["ann_hits"] += fl["ann"]
            zero = tuple(0 for _ in range(G.rank))
            nonzero = [g for g in G if g != zero]
            for tri in combinations(nonzero, 3):
                u = (zero,) + tri
                mult = Counter(G.add(b, t) for t in u for b in B)
                sigma = frozenset(y for y, c in mult.items() if c & 1)
                out["pairs_tested"] += 1
                if not sigma:
                    out["ann_hits"] += 1
                    continue
                hit = formsL.get(canon_translate(sigma, G))
                if hit is not None:
                    out["matches"] += 1
                    if len(out["match_examples"]) < 5:
                        out["match_examples"].append({
                            "frame": r["frame"], "A": r["A"],
                            "B": r["B"],
                            "u_L": [list(t) for t in hit],
                            "u_R": [list(t) for t in u],
                            "size": len(sigma)})
            out["members"] += 1
    out["secs"] = round(time.time() - t0, 1)
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    Path(args.out).write_text(json.dumps(out, indent=1))
    print(json.dumps(out))
    print(f"\nE18d (4,4): {'CLEAN' if out['matches'] == 0 else 'MATCHES!'}"
          f" ({out['pairs_tested']} u_R forms x full u_L form-set, "
          f"{out['members']} members, ann={out['ann_hits']})")


if __name__ == "__main__":
    main()
