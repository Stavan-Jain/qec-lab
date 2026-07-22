"""A23 phase 5: EXACT packing analysis over the true weight-8 detector universe.

Census facts (a23_w8_cycle_census.py): the weight-8 cycles of the base are
exactly 75 = one G-orbit of logicals; 40 pair odd with c1; all 40 are
S-clean, forming 8 S-orbits of size 5.

Here:
  1. S-orbit triple search: 3 pairwise-disjoint S-orbits => 15 detectors.
  2. Exact max independent set (branch & bound, bitsets) over the 40.
  3. If < 15: mix in annealed weight-10/12 shapes and re-run exact/greedy.
  4. Fractional (PB) retry on the enriched pool if still short.
"""

from __future__ import annotations

import json
import sys
import time
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
    seam_maps,
    translate1,
)
from a23_targeted_packing import S_GENS, anneal_class  # noqa: E402


def to_mask(sup: list[int]) -> int:
    m = 0
    for j in sup:
        m |= 1 << j
    return m


def exact_mis(masks: list[int], ub_target: int | None = None,
              time_cap: float = 120.0) -> list[int]:
    """Exact max independent set on the overlap-conflict graph (bitset B&B)."""
    n = len(masks)
    adj = [0] * n
    for i in range(n):
        for j in range(i + 1, n):
            if masks[i] & masks[j]:
                adj[i] |= 1 << j
                adj[j] |= 1 << i
    best: list[int] = []
    t0 = time.time()
    aborted = [False]

    def bb(cand: int, cur: list[int]) -> None:
        nonlocal best
        if aborted[0]:
            return
        if time.time() - t0 > time_cap:
            aborted[0] = True
            return
        if not cand:
            if len(cur) > len(best):
                best = cur[:]
            return
        # bound: |cur| + popcount(cand)
        if len(cur) + bin(cand).count("1") <= len(best):
            return
        if ub_target and len(best) >= ub_target:
            return
        # pick vertex with max degree within cand
        c = cand
        v = (c & -c).bit_length() - 1
        # branch: v in / v out
        bb(cand & ~(1 << v) & ~adj[v], cur + [v])
        bb(cand & ~(1 << v), cur)

    bb((1 << n) - 1, [])
    if aborted[0]:
        print(f"    (B&B time-capped at {time_cap}s; result is a lower bound)")
    return best


def main() -> None:
    rng = np.random.default_rng(99)
    MA = conv_matrix(A_SUPP)
    MB = conv_matrix(B_SUPP)
    D2 = np.vstack([MA, MB])
    D1 = np.hstack([MB, MA])
    K = nullspace_f2(D2)
    zeta = K[0]
    _, c1 = seam_maps(zeta)

    census = json.loads(
        (Path(__file__).resolve().parents[1] / "data/a23/w8_census.json").read_text()
    )
    dets8 = [np.zeros(150, dtype=np.uint8) for _ in census["detector_supports_odd_c1"]]
    for v, sup in zip(dets8, census["detector_supports_odd_c1"]):
        v[sup] = 1
    print(f"w8 odd-pairing detectors: {len(dets8)}")

    # --- 1. S-orbit structure ---------------------------------------------
    orb_of: dict[bytes, int] = {}
    orbit_reps: list[np.ndarray] = []
    orbit_members: list[list[int]] = []
    for i, v in enumerate(dets8):
        b = v.tobytes()
        if b in orb_of:
            orbit_members[orb_of[b]].append(i)
            continue
        oid = len(orbit_reps)
        for s in [(0, 0)] + S_GENS:
            orb_of[translate1(v, s).tobytes()] = oid
        orbit_reps.append(v)
        orbit_members.append([i])
    print(f"S-orbits among the 40: {len(orbit_reps)}")
    sats = []
    for rep in orbit_reps:
        sat = rep.copy()
        for s in S_GENS:
            sat |= translate1(rep, s)
        sats.append(sat)
        assert int(sat.sum()) == 40, "S-orbit not free/disjoint?"
    n_orb = len(orbit_reps)
    disj = np.zeros((n_orb, n_orb), dtype=bool)
    for i in range(n_orb):
        for j in range(i + 1, n_orb):
            disj[i, j] = disj[j, i] = not (sats[i] & sats[j]).any()
    print("S-orbit disjointness matrix (upper):")
    for i in range(n_orb):
        print("   ", "".join("." if i == j else ("1" if disj[i, j] else "0")
                             for j in range(n_orb)))
    triples = [
        (i, j, k)
        for i in range(n_orb)
        for j in range(i + 1, n_orb)
        for k in range(j + 1, n_orb)
        if disj[i, j] and disj[i, k] and disj[j, k]
    ]
    print(f"disjoint S-orbit triples: {len(triples)}")
    if triples:
        i, j, k = triples[0]
        print(f"  FOUND: orbits {i},{j},{k} -> 15 disjoint w8 detectors!")

    # --- 2. exact MIS over the 40 -----------------------------------------
    masks8 = [to_mask(sorted(map(int, np.flatnonzero(v)))) for v in dets8]
    best8 = exact_mis(masks8, time_cap=60)
    print(f"exact max disjoint family within w8 universe: {len(best8)}")

    if len(best8) >= 15 or triples:
        chosen = ([m for t in triples[:1] for o in t for m in orbit_members[o]]
                  if triples else best8)
        dets = [dets8[i] for i in chosen] if not triples else \
               [dets8[m] for t in triples[:1] for o in t for m in orbit_members[o]]
        _save(dets, c1, D2)
        return

    # --- 3. enrich with w10/w12 anneal ------------------------------------
    print("enriching with annealed w10/w12 shapes...")
    Zdual = nullspace_f2(D2.T)
    Rd, pivd = rref_f2(D1)
    dual_bnd = Rd[: len(pivd)]
    cur = dual_bnd.copy()
    dual_log_basis = []
    for row in Zdual:
        if rank_f2(np.vstack([cur, row])) > rank_f2(cur):
            cur = np.vstack([cur, row])
            dual_log_basis.append(row)
    dual_log_basis = np.array(dual_log_basis, dtype=np.uint8)
    prof = np.array([int((c1 & z).sum() % 2) for z in dual_log_basis])
    odd_classes = [
        eps for eps in range(1, 256)
        if (sum(prof[j] for j in range(8) if (eps >> j) & 1) % 2) == 1
    ]
    gen_rows = D1.astype(np.uint8)
    shapes: dict[bytes, np.ndarray] = {}
    for eps in odd_classes:
        rep = np.zeros(150, dtype=np.uint8)
        for j in range(8):
            if (eps >> j) & 1:
                rep ^= dual_log_basis[j]
        shapes.update(anneal_class(rep, gen_rows, rng, steps=6000, wmax=12))
    # translates, odd pairing
    pool: dict[bytes, np.ndarray] = {v.tobytes(): v for v in dets8}
    for s in shapes.values():
        for ga in range(LX):
            for gb in range(LY):
                t = translate1(s, (ga, gb))
                if int((t & c1).sum()) % 2 == 1:
                    pool[t.tobytes()] = t
    pool_list = list(pool.values())
    wt_hist: dict[int, int] = {}
    for v in pool_list:
        wt_hist[int(v.sum())] = wt_hist.get(int(v.sum()), 0) + 1
    print(f"enriched pool: {len(pool_list)}, weights {dict(sorted(wt_hist.items()))}")
    masks = [to_mask(sorted(map(int, np.flatnonzero(v)))) for v in pool_list]
    best = exact_mis(masks, ub_target=15, time_cap=240)
    print(f"max disjoint family over enriched pool: {len(best)}")
    if len(best) >= 15:
        _save([pool_list[i] for i in best], c1, D2)


def _save(dets: list[np.ndarray], c1: np.ndarray, D2: np.ndarray) -> None:
    cover = np.zeros(150, dtype=int)
    for i, z in enumerate(dets):
        assert not ((D2.T @ z) % 2).any()
        assert int((z & c1).sum()) % 2 == 1
        cover += z
    assert cover.max() <= 1, "not disjoint!"
    print(f"VERIFIED {len(dets)} disjoint odd-pairing detectors "
          f"(support {int(cover.sum())}/150)")
    outp = Path(__file__).resolve().parents[1] / "data/a23/disjoint15_rep1.json"
    outp.write_text(json.dumps({
        "n": len(dets),
        "detectors": [sorted(map(int, np.flatnonzero(z))) for z in dets],
    }, indent=1))
    print(f"wrote {outp}")


if __name__ == "__main__":
    main()
