"""A23 phase 2: disjoint dual-detector packing for the seam-coset floor.

Goal: for the seam class c = [seamC zeta] (one G-orbit, rep = mask-1 kernel
element), find N >= 15 dual cycles z_1..z_N with
  (i)   D2^T z_i = 0            (z_i kills every boundary in the pairing),
  (ii)  <seamC zeta, z_i> = 1   (odd pairing with the fixed coset offset),
  (iii) supp z_i pairwise disjoint.
Then every w = seamC zeta + d2 f has |w| >= N; parity (all coset weights
even) lifts N = 15 to the tight floor 16.

Search strategy:
  pool  = translates of harvested light dual cycles (reflect-swap of the
          weight-8 ustar + annealed extra shapes), filtered by (ii);
  pack  = greedy + randomized restarts + 2-opt local search on the
          disjointness conflict graph.

The order-5 stabilizer S = <(1,3)> of the kernel rep acts freely: any
S-clean shape contributes 5 mutually disjoint odd-pairing translates.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from bb_lab.linalg import nullspace_f2, rank_f2, rref_f2  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
from a23_seam_calibration import (  # noqa: E402
    A_SUPP,
    B_SUPP,
    LX,
    LY,
    conv_matrix,
    flat,
    seam_maps,
    translate1,
)


def build() -> dict:
    MA = conv_matrix(A_SUPP)
    MB = conv_matrix(B_SUPP)
    D2 = np.vstack([MA, MB])
    D1 = np.hstack([MB, MA])
    K = nullspace_f2(D2)
    # canonical rep: mask 1 = first basis row
    zeta = K[0].copy()
    _, c1 = seam_maps(zeta)
    return {"D2": D2, "D1": D1, "K": K, "zeta": zeta, "c1": c1}


def harvest_light_dual_cycles(
    ctx: dict, rng: np.random.Generator, n_rounds: int = 4000, wmax: int = 10
) -> list[np.ndarray]:
    """Anneal within dual cycle space for light codewords."""
    D2 = ctx["D2"]
    Zdual = nullspace_f2(D2.T)  # rows: 79 basis vectors
    nb = Zdual.shape[0]
    pool: dict[bytes, np.ndarray] = {}

    # seed: reflect-swap of ustar translates handled separately; here anneal
    for _ in range(n_rounds):
        w = np.zeros(150, dtype=np.uint8)
        # random sparse combo
        for i in rng.choice(nb, size=rng.integers(1, 4), replace=False):
            w ^= Zdual[i]
        # greedy descent over basis rows
        improved = True
        while improved:
            improved = False
            order = rng.permutation(nb)
            cw = int(w.sum())
            for i in order:
                cand = w ^ Zdual[i]
                if int(cand.sum()) < cw:
                    w = cand
                    cw = int(cand.sum())
                    improved = True
        if 0 < int(w.sum()) <= wmax:
            pool[w.tobytes()] = w.copy()
    return list(pool.values())


def reflect75(v75: np.ndarray) -> np.ndarray:
    arr = v75.reshape(LX, LY)
    out = np.zeros_like(arr)
    for a in range(LX):
        for b in range(LY):
            out[(-a) % LX, (-b) % LY] = arr[a, b]
    return out.reshape(-1)


def reflect_swap(w: np.ndarray) -> np.ndarray:
    return np.concatenate([reflect75(w[75:]), reflect75(w[:75])])


def max_disjoint(cands: list[np.ndarray], rng: np.random.Generator,
                 n_restarts: int = 4000, target: int = 15) -> list[int]:
    """Randomized greedy max independent set in the overlap-conflict graph."""
    n = len(cands)
    sup = [frozenset(np.flatnonzero(c)) for c in cands]
    # adjacency: overlap
    adj = [set() for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            if sup[i] & sup[j]:
                adj[i].add(j)
                adj[j].add(i)
    best: list[int] = []
    wts = np.array([len(s) for s in sup])
    for it in range(n_restarts):
        if it % 3 == 0:
            order = np.argsort(wts + rng.random(n))  # light-first
        else:
            order = rng.permutation(n)
        chosen: list[int] = []
        blocked: set[int] = set()
        for i in order:
            if i in blocked:
                continue
            chosen.append(i)
            blocked.add(i)
            blocked |= adj[i]
        if len(chosen) > len(best):
            best = chosen
            if len(best) >= target:
                break
    return best


def main() -> None:
    rng = np.random.default_rng(2323)
    ctx = build()
    D2, D1, c1 = ctx["D2"], ctx["D1"], ctx["c1"]

    import json

    data_p = Path(__file__).resolve().parents[1] / "data/a17/f2a6_z5z30_lean_data.json"
    data = json.loads(data_p.read_text())
    w = np.zeros(150, dtype=np.uint8)
    for (a, b) in data["ustar_left"]:
        w[flat(a, b)] = 1
    for (a, b) in data["ustar_right"]:
        w[75 + flat(a, b)] = 1
    zstar = reflect_swap(w)

    # shape pool: annealed light dual cycles + zstar
    print("harvesting light dual cycles (annealing)...")
    shapes = harvest_light_dual_cycles(ctx, rng, n_rounds=3000, wmax=10)
    shapes.append(zstar)
    wt_hist: dict[int, int] = {}
    for s in shapes:
        wt_hist[int(s.sum())] = wt_hist.get(int(s.sum()), 0) + 1
    print(f"  harvested {len(shapes)} distinct shapes, weights {dict(sorted(wt_hist.items()))}")

    # candidates: all 75 translates of each shape, filtered by odd pairing
    cands: dict[bytes, np.ndarray] = {}
    for s in shapes:
        for ga in range(LX):
            for gb in range(LY):
                t = translate1(s, (ga, gb))
                if int((t & c1).sum()) % 2 == 1:
                    cands[t.tobytes()] = t
    cand_list = list(cands.values())
    print(f"  candidate detectors (odd pairing with c1): {len(cand_list)}")

    # sanity: all candidates are dual cycles
    for t in cand_list[:20]:
        assert not ((D2.T @ t) % 2).any()

    best = max_disjoint(cand_list, rng, n_restarts=6000, target=16)
    print(f"packing: best disjoint family size = {len(best)}")
    total = sum(int(cand_list[i].sum()) for i in best)
    print(f"  total support = {total}/150")
    if len(best) >= 15:
        print("  SUCCESS: >= 15 disjoint odd-pairing detectors exist")
        out = {
            "class_rep_mask": 1,
            "detectors": [sorted(map(int, np.flatnonzero(cand_list[i])))
                          for i in best],
        }
        outp = Path(__file__).resolve().parents[1] / "data/a23/packing_rep1.json"
        outp.parent.mkdir(exist_ok=True)
        outp.write_text(json.dumps(out, indent=1))
        print(f"  wrote {outp}")
    else:
        print("  below target; need richer pool or different mechanism")


if __name__ == "__main__":
    main()
