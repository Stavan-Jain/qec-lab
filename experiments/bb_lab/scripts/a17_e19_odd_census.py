#!/usr/bin/env python3
"""A17 E19 — odd–odd splits: (3,3) / (3,5) / (1,7) census battery.

Derived reductions (doc §11):
  (3,3): the cap (A-families of one B-translate's σ-cells pairwise
    disjoint ⟹ ≤ 3 cells) + (†) force a_B = 3 = a_A — BOTH u's are
    triangles (Cay(G,dA) / Cay(G,dB)) — ñ₃ = 0, |σ| = 9, and σ is a
    bijective 3×3 grid y_ij ∈ (A+g_i) ∩ (B+r_j). D2 makes
    (a,b) ↦ b−a injective (|B−A| = 25 exactly), so grid witnesses
    are FORCED: the surviving demand is a rank-1 shift grid
    {r_j − g_i} of 9 distinct elements inside S = B−A with
    row/column-injective witnesses.
  (3,5): same cap ⟹ a_B(u_R) ≥ 5 (dense 5-set), |σ| ∈ {9,11,13,15}.
  (1,7): σ = A+g and every B-translate meets σ in ≤ 1 cell ⟹
    ñ₃ ≤ 1, ñ₅ = ñ₇ = 0, a_B(u_R) ∈ {15,17,19,21} — a ≥ 15-of-21
    super-dense 7-set (complete enumeration via the anchor argument:
    Σ complement-degrees ≤ 12 ⟹ some vertex adjacent to ≥ 5 of 6).

Censused per member:
  O1 (3,3) end-to-end: all 1,891 weight-3 forms vs forms. Expect 0.
  O2 (3,5) end-to-end: weight-3 forms vs all weight-5 u_R with
     a_B ≥ 5 (valid prune: |σ| ≥ 25 − 2a_B). Expect 0.
  O3 (1,7): dense-7 enumeration, match σ_R = A + g. Expect 0.
  O4 (3,3) grid demand: dA-triangles × dB-triangles with all 9
     shifts in S + forced-witness injectivity. Expect 0 full passes.

Usage: uv run python scripts/a17_e19_odd_census.py \
    --members data/a17/members_7x9.jsonl,data/a17/members_6x9_6x10.jsonl \
    --out data/a17/e19_odd_census.json
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


def forms3(S, G, nonzero, zero):
    """canon(supp(S*u)) for every weight-3 u containing 0."""
    out = {}
    n_ann = 0
    sizes = Counter()
    for pair in combinations(nonzero, 2):
        u = (zero,) + pair
        mult = Counter(G.add(s, t) for t in u for s in S)
        sigma = frozenset(y for y, c in mult.items() if c & 1)
        if not sigma:
            n_ann += 1
            continue
        sizes[len(sigma)] += 1
        out.setdefault(canon_translate(sigma, G), u)
    return out, n_ann, sizes


def triangles(dX: frozenset, G) -> list:
    """All triangles of Cay(G, dX) up to translation, anchored at 0."""
    dl = sorted(dX)
    seen = set()
    out = []
    zero = tuple(0 for _ in range(G.rank))
    for i, d1 in enumerate(dl):
        for d2 in dl[i + 1:]:
            if G.sub(d2, d1) not in dX:
                continue
            T = frozenset((zero, d1, d2))
            c = canon_translate(T, G)
            if c not in seen:
                seen.add(c)
                out.append(c)
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--members", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    t0 = time.time()

    out = {"members": 0,
           "O1": {"pairs": 0, "matches": 0, "ann_L": 0, "ann_R": 0},
           "O2": {"uR_tested": 0, "uR_dense": 0, "matches": 0},
           "O3": {"candidates": 0, "matches": 0},
           "O4": {"tri_A": 0, "tri_B": 0, "grid_pairs": 0,
                  "shift_pass": 0, "witness_pass": 0, "examples": []},
           "S_size_ok": 0}
    for path in args.members.split(","):
        with open(path) as f:
            rows = [json.loads(l) for l in f if '"A"' in l]
        for r in rows:
            G = AbelianGroup(tuple(r["frame"]))
            A = sorted(Poly.from_string(r["A"], G).support)
            B = sorted(Poly.from_string(r["B"], G).support)
            Aset, Bset = frozenset(A), frozenset(B)
            dAs, dBs = diffs(Aset, G), diffs(Bset, G)
            zero = tuple(0 for _ in range(G.rank))
            nonzero = [g for g in G if g != zero]
            canonA = canon_translate(Aset, G)

            fL, annL, _ = forms3(A, G, nonzero, zero)
            out["O1"]["ann_L"] += annL

            # O1: (3,3)
            for pair in combinations(nonzero, 2):
                u = (zero,) + pair
                mult = Counter(G.add(b, t) for t in u for b in B)
                sigma = frozenset(y for y, c in mult.items() if c & 1)
                out["O1"]["pairs"] += 1
                if not sigma:
                    out["O1"]["ann_R"] += 1
                    continue
                if canon_translate(sigma, G) in fL:
                    out["O1"]["matches"] += 1

            # O2: (3,5) — u_R weight 5, prune a_B >= 5
            for quad in combinations(nonzero, 4):
                u = (zero,) + quad
                out["O2"]["uR_tested"] += 1
                aB = sum(1 for x, y in combinations(u, 2)
                         if G.sub(y, x) in dBs)
                if aB < 5:
                    continue
                out["O2"]["uR_dense"] += 1
                mult = Counter(G.add(b, t) for t in u for b in B)
                sigma = frozenset(y for y, c in mult.items() if c & 1)
                if len(sigma) not in (9, 11, 13, 15):
                    continue
                if canon_translate(sigma, G) in fL:
                    out["O2"]["matches"] += 1

            # O3: (1,7) — super-dense 7-sets, a >= 15
            seen7 = set()
            dBl = sorted(dBs)
            for s5 in combinations(dBl, 5):
                p5 = sum(1 for x, y in combinations(s5, 2)
                         if G.sub(y, x) in dBs)
                if p5 + 5 + 6 < 15:
                    continue
                base = (zero,) + s5
                bset = frozenset(base)
                for v7 in nonzero:
                    if v7 in bset:
                        continue
                    deg = sum(1 for x in base if G.sub(v7, x) in dBs)
                    if p5 + 5 + deg < 15:
                        continue
                    T = bset | {v7}
                    c = canon_translate(T, G)
                    if c in seen7:
                        continue
                    seen7.add(c)
                    out["O3"]["candidates"] += 1
                    mult = Counter(G.add(b, t) for t in T for b in B)
                    sigma = frozenset(
                        y for y, cnt in mult.items() if cnt & 1)
                    if len(sigma) == 5 and \
                            canon_translate(sigma, G) == canonA:
                        out["O3"]["matches"] += 1

            # O4: (3,3) grid demand on triangles
            S_map = {}
            for a in A:
                for b in B:
                    S_map[G.sub(b, a)] = (a, b)
            if len(S_map) == 25:
                out["S_size_ok"] += 1
            triA = triangles(dAs, G)
            triB = triangles(dBs, G)
            out["O4"]["tri_A"] += len(triA)
            out["O4"]["tri_B"] += len(triB)
            for Tg in triA:
                for Tr in triB:
                    out["O4"]["grid_pairs"] += 1
                    wit = {}
                    ok = True
                    for i, g in enumerate(Tg):
                        for j, rr in enumerate(Tr):
                            s = G.sub(rr, g)
                            w = S_map.get(s)
                            if w is None:
                                ok = False
                                break
                            wit[(i, j)] = w
                        if not ok:
                            break
                    if not ok:
                        continue
                    out["O4"]["shift_pass"] += 1
                    rows_ok = all(
                        len({wit[(i, j)][0] for j in range(3)}) == 3
                        and len({wit[(i, j)][1] for j in range(3)}) == 3
                        for i in range(3))
                    cols_ok = all(
                        len({wit[(i, j)][0] for i in range(3)}) == 3
                        and len({wit[(i, j)][1] for i in range(3)}) == 3
                        for j in range(3))
                    if rows_ok and cols_ok:
                        out["O4"]["witness_pass"] += 1
                        if len(out["O4"]["examples"]) < 5:
                            out["O4"]["examples"].append({
                                "frame": r["frame"], "A": r["A"],
                                "B": r["B"],
                                "Tg": [list(t) for t in Tg],
                                "Tr": [list(t) for t in Tr]})
            out["members"] += 1
    out["secs"] = round(time.time() - t0, 1)
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    Path(args.out).write_text(json.dumps(out, indent=1))
    print(json.dumps(out))
    bad = (out["O1"]["matches"] + out["O2"]["matches"]
           + out["O3"]["matches"] + out["O4"]["witness_pass"])
    print(f"\nE19: {'CLEAN' if bad == 0 else f'{bad} RESIDUALS'} "
          f"(S-bijection {out['S_size_ok']}/{out['members']})")


if __name__ == "__main__":
    main()
