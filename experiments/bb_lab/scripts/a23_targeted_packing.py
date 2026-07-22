"""A23 phase 2b: targeted per-dual-class harvest + S-structured packing.

Upgrades over a23_detector_packing.py:
  * harvest = per-DUAL-CLASS annealing (fix the logical part, walk the
    dual-boundary coset with weight-6 generator moves), targeting the 128
    dual classes pairing odd with c1;
  * S-structure: the stabilizer S = <(1,3)> of the kernel rep acts freely
    and fixes the seam class, so an S-clean shape (support differences
    avoid S\\0) contributes 5 mutually disjoint odd-pairing translates;
    three shapes with disjoint S-saturations give 15 detectors at once;
  * fallback: randomized packing over the enriched translate pool.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from bb_lab.linalg import nullspace_f2, rref_f2  # noqa: E402

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
from a23_detector_packing import max_disjoint, reflect_swap  # noqa: E402

S_GENS = [(1, 3), (2, 6), (3, 9), (4, 12)]  # S \ 0


def s_translates(v: np.ndarray) -> list[np.ndarray]:
    return [v] + [translate1(v, s) for s in S_GENS]


def is_s_clean(v: np.ndarray) -> bool:
    sup = np.flatnonzero(v)
    base = v
    for s in S_GENS:
        if (translate1(base, s) & v).any():
            return False
    _ = sup
    return True


def s_saturate(v: np.ndarray) -> np.ndarray:
    out = v.copy()
    for s in S_GENS:
        out |= translate1(v, s)
    return out


def anneal_class(
    rep: np.ndarray,
    gen_rows: np.ndarray,
    rng: np.random.Generator,
    steps: int = 6000,
    wmax: int = 10,
) -> dict[bytes, np.ndarray]:
    """Random walk in rep + <gen_rows>, collecting all weight <= wmax states."""
    found: dict[bytes, np.ndarray] = {}
    w = rep.copy()
    cw = int(w.sum())
    ngen = gen_rows.shape[0]
    best = cw
    for t in range(steps):
        i = rng.integers(ngen)
        cand = w ^ gen_rows[i]
        nw = int(cand.sum())
        # accept downhill/plateau always; uphill with decaying probability
        if nw <= cw or rng.random() < 0.02:
            w, cw = cand, nw
            if cw <= wmax and cw > 0:
                found[w.tobytes()] = w.copy()
            if cw < best:
                best = cw
        # occasional restart to the best-so-far basin
        if t % 1500 == 1499 and found:
            k = rng.integers(len(found))
            w = list(found.values())[k].copy()
            cw = int(w.sum())
    return found


def main() -> None:
    rng = np.random.default_rng(4242)
    MA = conv_matrix(A_SUPP)
    MB = conv_matrix(B_SUPP)
    D2 = np.vstack([MA, MB])
    D1 = np.hstack([MB, MA])
    K = nullspace_f2(D2)
    zeta = K[0]
    _, c1 = seam_maps(zeta)

    # sanity: S fixes zeta
    from a23_seam_calibration import translate2

    for s in S_GENS:
        assert np.array_equal(translate2(zeta, s), zeta)

    Zdual = nullspace_f2(D2.T)
    Rd, pivd = rref_f2(D1)
    dual_bnd = Rd[: len(pivd)]

    # dual logical basis (complement of dual_bnd in Zdual)
    from bb_lab.linalg import rank_f2

    cur = dual_bnd.copy()
    dual_log_basis = []
    for row in Zdual:
        if rank_f2(np.vstack([cur, row])) > rank_f2(cur):
            cur = np.vstack([cur, row])
            dual_log_basis.append(row)
    dual_log_basis = np.array(dual_log_basis, dtype=np.uint8)
    assert dual_log_basis.shape[0] == 8

    # 128 odd-pairing dual classes: eps with sum eps_i <c1,z_i> odd
    prof = np.array([int((c1 & z).sum() % 2) for z in dual_log_basis])
    odd_classes = [
        eps
        for eps in range(1, 256)
        if (sum(prof[j] for j in range(8) if (eps >> j) & 1) % 2) == 1
    ]
    print(f"odd-pairing dual classes: {len(odd_classes)}")

    gen_rows = D1.astype(np.uint8)  # 75 weight-6 dual-boundary generators

    # include the ustar-derived shape's class first (known weight-8)
    data_p = Path(__file__).resolve().parents[1] / "data/a17/f2a6_z5z30_lean_data.json"
    data = json.loads(data_p.read_text())
    from a23_seam_calibration import flat

    wst = np.zeros(150, dtype=np.uint8)
    for (a, b) in data["ustar_left"]:
        wst[flat(a, b)] = 1
    for (a, b) in data["ustar_right"]:
        wst[75 + flat(a, b)] = 1
    zstar = reflect_swap(wst)

    shapes: dict[bytes, np.ndarray] = {}
    if int((zstar & c1).sum()) % 2 == 1:
        shapes[zstar.tobytes()] = zstar

    print("targeted per-class annealing over odd classes...")
    for idx, eps in enumerate(odd_classes):
        rep = np.zeros(150, dtype=np.uint8)
        for j in range(8):
            if (eps >> j) & 1:
                rep ^= dual_log_basis[j]
        found = anneal_class(rep, gen_rows, rng, steps=5000, wmax=10)
        shapes.update(found)
        if idx % 16 == 15:
            print(f"  {idx+1}/{len(odd_classes)} classes, pool {len(shapes)}")
    wt_hist: dict[int, int] = {}
    for s in shapes.values():
        wt_hist[int(s.sum())] = wt_hist.get(int(s.sum()), 0) + 1
    print(f"pool: {len(shapes)} shapes, weights {dict(sorted(wt_hist.items()))}")

    # every harvested shape pairs odd with c1 by construction (class-level);
    # verify on a sample
    for s in list(shapes.values())[:10]:
        assert int((s & c1).sum()) % 2 == 1
        assert not ((D2.T @ s) % 2).any()

    # --- S-structured route ------------------------------------------------
    s_clean = [s for s in shapes.values() if is_s_clean(s)]
    print(f"S-clean shapes: {len(s_clean)}")
    sats = [s_saturate(s) for s in s_clean]
    hit = None
    order = np.argsort([int(s.sum()) for s in s_clean])
    for ii in range(len(s_clean)):
        i = order[ii]
        for jj in range(ii + 1, len(s_clean)):
            j = order[jj]
            if (sats[i] & sats[j]).any():
                continue
            for kk in range(jj + 1, len(s_clean)):
                k = order[kk]
                if (sats[i] & sats[k]).any() or (sats[j] & sats[k]).any():
                    continue
                hit = (i, j, k)
                break
            if hit:
                break
        if hit:
            break
    if hit:
        i, j, k = hit
        print(f"S-TRIPLE FOUND: weights "
              f"{[int(s_clean[t].sum()) for t in (i, j, k)]} -> 15 disjoint detectors")
        detectors = []
        for t in (i, j, k):
            detectors += s_translates(s_clean[t])
    else:
        print("no S-triple; falling back to general packing")
        # general pool: all odd-pairing translates
        cands: dict[bytes, np.ndarray] = {}
        for s in shapes.values():
            for ga in range(LX):
                for gb in range(LY):
                    t = translate1(s, (ga, gb))
                    if int((t & c1).sum()) % 2 == 1:
                        cands[t.tobytes()] = t
        cand_list = list(cands.values())
        print(f"translate pool: {len(cand_list)}")
        best = max_disjoint(cand_list, rng, n_restarts=8000, target=15)
        print(f"best disjoint family: {len(best)}")
        detectors = [cand_list[i] for i in best]

    if len(detectors) >= 15:
        # final verification
        ok = True
        for i, z in enumerate(detectors):
            assert not ((D2.T @ z) % 2).any(), "not a dual cycle"
            assert int((z & c1).sum()) % 2 == 1, "even pairing"
            for j in range(i):
                if (z & detectors[j]).any():
                    ok = False
        print(f"verified {len(detectors)} detectors, pairwise disjoint: {ok}, "
              f"total support {int(np.sum([z.sum() for z in detectors]))}/150")
        if ok:
            out = {
                "class_rep": "K[0] (mask 1)",
                "n_detectors": len(detectors),
                "detectors": [sorted(map(int, np.flatnonzero(z)))
                              for z in detectors],
            }
            outp = Path(__file__).resolve().parents[1] / "data/a23/packing_rep1.json"
            outp.parent.mkdir(exist_ok=True)
            outp.write_text(json.dumps(out, indent=1))
            print(f"wrote {outp}")


if __name__ == "__main__":
    main()
