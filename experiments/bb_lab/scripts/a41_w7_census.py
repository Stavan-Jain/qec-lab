#!/usr/bin/env python3
"""A41 G1b — the w = 7 structural census (first sweep of the class at
the next odd weight).

Same gate stack as `a15_t42_w5_sweep.py` (translation-normalized
mono-axis Sidon pools -> per-side zd filter -> D2 pairing -> k > 0 ->
canonical dedupe), with one engineering replacement: the pool
enumerator.  At w = 7 the sweep's `itertools.combinations` pool is
C(|G|-1, 6) ~ 4e8 per axis per frame — infeasible; here Sidon supports
are enumerated by DFS with incremental difference-set pruning (a
repeated difference kills the whole subtree), and the parity-mono
(iii) filter runs at the leaves.

Scope (recorded, not hidden): frames with BOTH axes >= 5, class frame
rule (4 does not divide either order), |G| >= 85 (the D1&D2 counting
bound 2*2*C(7,2) = 84 <= |G|-1).  Axes of order <= 4 are excluded by
cost, not by a proven kill — the w = 5 Lemma E small-axis argument is
weight-specific and has NOT been re-derived at w = 7.

Structural only (no SAT).  Output rows match the w = 5 census schema
(d: null), so the MC falsifier and the G3 screen consume them as-is.

Usage (from experiments/bb_lab):
    uv run python scripts/a41_w7_census.py --selftest
    uv run python scripts/a41_w7_census.py --frames 5x17,9x10
Output: data/a41/w7_census.jsonl + per-frame summary lines.
"""

from __future__ import annotations

import argparse
import itertools
import json
import sys
import time
from collections import Counter
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from a15_t42_w5_sweep import (canon_key, diff_index_table, mono_axes,
                              diff_multiset_free)
from bb_lab.checks import bb_check_matrices
from bb_lab.codeparams import code_params
from bb_lab.diffset_predicates import frobenius_square, is_translate
from bb_lab.group import AbelianGroup
from bb_lab.linalg import rank_f2
from bb_lab.poly import Poly

W = 7


def rank_f2_packed(rows_bits: np.ndarray, nbits: int) -> int:
    """GF(2) rank of a bit-packed matrix (rows of uint64 words) by
    vectorized elimination.  ~100x the uint8 rref at these sizes."""
    M = rows_bits.copy()
    nrows = M.shape[0]
    rank = 0
    for col in range(nbits):
        wrd, bit = col >> 6, np.uint64(1) << np.uint64(col & 63)
        colvals = (M[rank:, wrd] & bit) != 0
        hits = np.flatnonzero(colvals)
        if hits.size == 0:
            continue
        piv = rank + hits[0]
        if piv != rank:
            M[[rank, piv]] = M[[piv, rank]]
        rest = rank + 1 + np.flatnonzero((M[rank + 1:, wrd] & bit) != 0)
        if rest.size:
            M[rest] ^= M[rank]
        rank += 1
        if rank == nrows:
            break
    return rank


def pack_rows(mat: np.ndarray) -> np.ndarray:
    """uint8 (r, c) 0/1 matrix -> (r, ceil(c/64)) uint64 bit rows."""
    r, c = mat.shape
    words = (c + 63) // 64
    padded = np.zeros((r, words * 64), dtype=np.uint8)
    padded[:, :c] = mat
    return np.packbits(padded, axis=1, bitorder="little").reshape(
        r, words * 8).view(np.uint64).reshape(r, words)


def sidon_pools_both(G, idx: dict):
    """0-normalized Sidon W-supports by one DFS, split at the leaves
    into the mono-x and mono-y pools (each support is parity-mono in
    at most one axis, so the split is disjoint)."""
    elems = list(G)
    zero = (0,) * G.rank
    n = len(elems)
    pools = {0: [], 1: []}
    chosen = [zero]
    used: set = set()

    def leaf_check():
        supp = tuple(chosen)
        ax = mono_axes(supp, G.rank)
        if len(ax) != 1:
            return
        mask = 0
        for d in used:
            mask |= 1 << idx[d]
        pools[ax[0]].append((supp, mask))

    def dfs(start: int):
        if len(chosen) == W:
            leaf_check()
            return
        # prune: not enough elements left
        need = W - len(chosen)
        for i in range(start, n - need + 1):
            h = elems[i]
            new = []
            ok = True
            for s in chosen:
                d1, d2 = G.sub(h, s), G.sub(s, h)
                if d1 == d2:  # order-2 difference: ordered multiset
                    ok = False  # holds it TWICE -> never Sidon
                    break
                if d1 in used or d2 in used:
                    ok = False
                    break
                new.append(d1)
                new.append(d2)
            if ok and len(new) != len(set(new)):
                ok = False  # internal collision within the new batch
            if not ok:
                continue
            chosen.append(h)
            for d in new:
                used.add(d)
            dfs(i + 1)
            chosen.pop()
            for d in new:
                used.discard(d)

    dfs(1)
    return pools[0], pools[1]


def sweep_frame_w7(ell: int, m: int, out, log) -> dict:
    t0 = time.time()
    G = AbelianGroup((ell, m))
    n = ell * m
    elems = list(G)
    idx = {g: i for i, g in enumerate(elems)}

    from a15_t42_w5_sweep import translation_classes
    poolA_raw, poolB_raw = sidon_pools_both(G, idx)
    poolA = translation_classes(G, poolA_raw)
    poolB = translation_classes(G, poolB_raw)
    log(f"  pools (one DFS): A {len(poolA_raw)}->{len(poolA)} classes, "
        f"B {len(poolB_raw)}->{len(poolB)} ({time.time() - t0:.0f}s)")
    del poolA_raw, poolB_raw

    have_both = bool(poolA) and bool(poolB)
    dt = diff_index_table(G, elems, idx) if have_both else None

    def zd_filter(pool):
        keep = []
        vec = np.zeros(n, dtype=np.uint8)
        for supp, mask in pool:
            vec[:] = 0
            vec[[idx[s] for s in supp]] = 1
            circ = vec[dt]
            packed = pack_rows(circ)
            if rank_f2_packed(packed, n) < n:
                keep.append((supp, mask, packed))
        return keep

    # B first: the mono-y pool is the small side on tall frames, and an
    # empty zdB kills the frame before the (huge) A-side filter runs.
    zdB = zd_filter(poolB) if have_both else []
    zdA = zd_filter(poolA) if have_both and zdB else []
    log(f"  zd: A {len(zdA)}, B {len(zdB)} ({time.time() - t0:.0f}s)")

    n_pairs = n_kpos = 0
    reps: dict[tuple, tuple] = {}
    for sa, ma, ca in zdA:
        for sb, mb, cb in zdB:
            if ma & mb:
                continue
            n_pairs += 1
            # packed hstack = [A | zero-pad | B | zero-pad]: the pad
            # columns are zero, so the GF(2) rank is that of [A | B];
            # scan every word-bit.
            stacked = np.hstack([ca, cb])
            k = 2 * (n - rank_f2_packed(stacked, stacked.shape[1] * 64))
            if k == 0:
                continue
            n_kpos += 1
            key = canon_key(G, sa, sb)
            if key not in reps:
                reps[key] = (sa, sb, k)

    rows = 0
    for sa, sb, k in reps.values():
        A = Poly.from_support(sa, G)
        B = Poly.from_support(sb, G)
        k_rec = code_params(bb_check_matrices(A, B)).k
        A2, B2 = frobenius_square(A), frobenius_square(B)
        row = {
            "frame": [ell, m], "n": 2 * n, "k": int(k_rec),
            "A": A.canonical_string(), "B": B.canonical_string(),
            "flags": {
                "frob2": is_translate(A, B2) or is_translate(B, A2),
                "frob4": is_translate(A, frobenius_square(B2))
                         or is_translate(B, frobenius_square(A2)),
                "k_fast_vs_record": int(k) == int(k_rec),
            },
            "d": None, "w": W,
        }
        out.write(json.dumps(row) + "\n")
        out.flush()
        rows += 1

    summary = {"frame": f"Z{ell}xZ{m}", "G": n, "w": W,
               "poolA": len(poolA), "poolB": len(poolB),
               "zdA": len(zdA), "zdB": len(zdB),
               "d2_pairs": n_pairs, "kpos_pairs": n_kpos,
               "members": len(reps), "rows": rows,
               "secs": round(time.time() - t0, 1)}
    log(json.dumps(summary))
    return summary


def selftest() -> None:
    """DFS pools must equal the combinations pools at w = 5 on a small
    frame (both axes), and the packed rank must agree with rank_f2."""
    import a15_t42_w5_sweep as w5
    G = AbelianGroup((6, 9))
    elems = list(G)
    idx = {g: i for i, g in enumerate(elems)}
    global W
    W_saved = W
    W = 5
    try:
        pa, pb = sidon_pools_both(G, idx)
        dfsA = {frozenset(s) for s, _ in pa}
        dfsB = {frozenset(s) for s, _ in pb}
        refA = {frozenset(s) for s, _ in w5.enumerate_pool(G, 5, 0, idx)}
        refB = {frozenset(s) for s, _ in w5.enumerate_pool(G, 5, 1, idx)}
        assert dfsA == refA, (len(dfsA), len(refA))
        assert dfsB == refB, (len(dfsB), len(refB))
    finally:
        W = W_saved
    # packed rank vs reference on 300 random 0/1 matrices (incl. the
    # hstack-with-pad shape)
    rng = np.random.default_rng(20260825)
    for _ in range(300):
        r, c = rng.integers(3, 60), rng.integers(3, 90)
        M = rng.integers(0, 2, size=(r, c)).astype(np.uint8)
        assert rank_f2_packed(pack_rows(M), int(c)) == rank_f2(M)
        two = np.hstack([pack_rows(M), pack_rows(M[:, ::-1])])
        assert rank_f2_packed(two, two.shape[1] * 64) == rank_f2(
            np.hstack([M, M[:, ::-1]]))
    print(f"SELFTEST PASS (w=5 DFS pools == combinations on Z6xZ9: "
          f"A {len(dfsA)}, B {len(dfsB)} classes*5; packed rank == "
          f"rank_f2 on 300 random + padded-hstack cases)")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--frames", type=str,
                    default="5x17,5x18,6x15,9x10")
    ap.add_argument("--out", type=str, default="data/a41/w7_census.jsonl")
    ap.add_argument("--selftest", action="store_true")
    args = ap.parse_args()
    if args.selftest:
        selftest()
        return
    frames = []
    for tok in args.frames.split(","):
        a, b = tok.lower().split("x")
        frames.append((int(a), int(b)))
    outp = Path(args.out)
    outp.parent.mkdir(parents=True, exist_ok=True)
    log = lambda s: print(s, flush=True)
    log(f"w=7 frames: {[f'{a}x{b}' for a, b in frames]}")
    summaries = []
    with outp.open("a") as out:
        for ell, m in frames:
            log(f"== Z{ell}xZ{m}")
            summaries.append(sweep_frame_w7(ell, m, out, log))
    total = sum(s["members"] for s in summaries)
    log(f"TOTAL w=7 members: {total} across {len(frames)} frames")


if __name__ == "__main__":
    main()
