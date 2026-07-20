#!/usr/bin/env python3
"""A17 E17b — verification for Theorem J5-odd and Prop J5-even
(confirmation of hand proofs, never an ingredient).

W1: over full pools on given frames: every mono-x Sidon A has
    Z_y(dA) ∈ {0, 2, 6} (ordered y-free differences) and every
    mono-y Sidon B has Z_y(dB) ∈ {4, 8, 12, 20} — the two value sets
    of the J5-odd disjointness argument.
W2 (odd frames): |2·dB ∩ {y = 0}| = Z_y(dB) (doubling is an
    automorphism preserving y-free) — the theorem's transport step.
W3 (even-ℓ frames, m odd): no B-(4,1) has a fully τ-closed 4-row
    difference set (Prop J5-even(b): full closure is D1-dead).

Usage: uv run python scripts/a17_j5_shape_checks.py \
    --frames 7x9,6x9,6x10 --out data/a17/e17_j5_shapes.json
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


def zcomp(S, G, comp: int) -> int:
    """Ordered differences whose `comp` component is zero (comp=1:
    y-free, from same-y pairs; comp=0: x-free, from same-x pairs)."""
    c = 0
    for a in S:
        for b in S:
            if a != b and (G.sub(a, b))[comp] == 0:
                c += 1
    return c


def check_frame(ell: int, m: int) -> dict:
    t0 = time.time()
    G = AbelianGroup((ell, m))
    elems = list(G)
    idx = {g: i for i, g in enumerate(elems)}
    poolA = translation_classes(G, enumerate_pool(G, 5, 0, idx))
    poolB = translation_classes(G, enumerate_pool(G, 5, 1, idx))
    odd_frame = (ell % 2 == 1) and (m % 2 == 1)

    zyA: Counter = Counter()   # y-free diffs of mono-x A: expect {0,2,6}
    zxA: Counter = Counter()   # x-free diffs of A (mirror check)
    for supp, _ in poolA:
        zyA[zcomp(supp, G, 1)] += 1
        zxA[zcomp(supp, G, 0)] += 1
    zyB: Counter = Counter()   # y-free diffs of mono-y B: expect {4,8,12,20}
    zxB: Counter = Counter()   # x-free diffs of B (mirror): expect {0,2,6}
    w2_bad = 0
    for supp, _ in poolB:
        z = zcomp(supp, G, 1)
        zyB[z] += 1
        zxB[zcomp(supp, G, 0)] += 1
        if odd_frame:
            dB = {G.sub(a, b) for a in supp for b in supp if a != b}
            two_y0 = {G.add(d, d) for d in dB if d[1] == 0}
            if len(two_y0) != z:
                w2_bad += 1

    w1_ok = (set(zyA) <= {0, 2, 6} and set(zyB) <= {4, 8, 12, 20}
             and set(zxB) <= {0, 2, 6})

    # W3: B-(4,1) full τ-closure absence on even-ℓ frames
    w3 = None
    if ell % 2 == 0:
        tau = (ell // 2, 0)
        closed = 0
        n41 = 0
        for supp, _ in poolB:
            rows: Counter = Counter(p[1] for p in supp)
            if sorted(rows.values()) != [1, 4]:
                continue
            n41 += 1
            y4 = next(y for y, c in rows.items() if c == 4)
            xs = [p[0] for p in supp if p[1] == y4]
            D = {(a - b) % ell for a in xs for b in xs if a != b}
            if all(((u + tau[0]) % ell) in D for u in D):
                closed += 1
        w3 = {"n_B41": n41, "fully_tau_closed": closed}

    return {
        "frame": f"Z{ell}xZ{m}", "odd_frame": odd_frame,
        "zyA_hist": {str(k): v for k, v in sorted(zyA.items())},
        "zyB_hist": {str(k): v for k, v in sorted(zyB.items())},
        "zxA_hist": {str(k): v for k, v in sorted(zxA.items())},
        "zxB_hist": {str(k): v for k, v in sorted(zxB.items())},
        "W1_ok": w1_ok, "W2_bad": w2_bad if odd_frame else None,
        "W3": w3, "secs": round(time.time() - t0, 1),
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--frames", default="7x9,6x9,6x10")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    results = []
    for tok in args.frames.split(","):
        ell, m = (int(x) for x in tok.lower().split("x"))
        r = check_frame(ell, m)
        results.append(r)
        print(json.dumps(r), flush=True)
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    Path(args.out).write_text(json.dumps(results, indent=1))
    ok = all(r["W1_ok"] and not r["W2_bad"]
             and (r["W3"] is None or r["W3"]["fully_tau_closed"] == 0)
             for r in results)
    print(f"\nALL CHECKS {'PASS' if ok else 'FAIL'}")


if __name__ == "__main__":
    main()
