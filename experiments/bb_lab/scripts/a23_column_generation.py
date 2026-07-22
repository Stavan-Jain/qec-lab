"""A23 phase 7: column generation -- is the TRUE packing LP below 14?

Loop:
  1. solve the packing LP on the current pool (HiGHS), get dual cover mu;
  2. hunt (annealing, mu-mass objective) for an odd-pairing dual cycle z
     with mu-mass < 1 (a violated dual constraint = improving column);
  3. if found, add all its odd-pairing translates and re-solve; else stop.

Convergence with no violated column ==> the dual mu is feasible for the
FULL universe (up to annealing power), i.e. the true LP <= sum(mu).

Run:  uv run --with scipy --project experiments/bb_lab python \
          experiments/bb_lab/scripts/a23_column_generation.py
"""

from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np
from scipy.optimize import linprog

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
from a23_lp_estimate import build_or_load_pool  # noqa: E402


def solve_lp(P: np.ndarray) -> tuple[float, np.ndarray]:
    res = linprog(
        c=-np.ones(P.shape[0]),
        A_ub=P.T.astype(float),
        b_ub=np.ones(P.shape[1]),
        bounds=(0, None),
        method="highs",
    )
    assert res.status == 0
    return -res.fun, -res.ineqlin.marginals


def mu_anneal(
    rep: np.ndarray,
    gen_rows: np.ndarray,
    mu: np.ndarray,
    rng: np.random.Generator,
    steps: int = 8000,
) -> tuple[float, np.ndarray]:
    """Anneal within the dual class coset of rep minimizing mu-mass."""
    w = rep.copy()

    def mass(v: np.ndarray) -> float:
        return float(mu[v > 0].sum())

    cw = mass(w)
    best_w, best_m = w.copy(), cw
    ngen = gen_rows.shape[0]
    for t in range(steps):
        i = rng.integers(ngen)
        cand = w ^ gen_rows[i]
        nm = mass(cand)
        if nm <= cw + 1e-12 or rng.random() < 0.03:
            w, cw = cand, nm
            if cw < best_m:
                best_w, best_m = w.copy(), cw
    return best_m, best_w


def main() -> None:
    rng = np.random.default_rng(31415)
    P, c1 = build_or_load_pool()
    pool: dict[bytes, np.ndarray] = {bytes(r.tobytes()): r for r in P}

    MA = conv_matrix(A_SUPP)
    MB = conv_matrix(B_SUPP)
    D2 = np.vstack([MA, MB])
    D1 = np.hstack([MB, MA])
    gen_rows = D1.astype(np.uint8)

    Zdual = nullspace_f2(D2.T)
    Rd, pivd = rref_f2(D1)
    dual_bnd = Rd[: len(pivd)]
    cur = dual_bnd.copy()
    dlb = []
    for row in Zdual:
        if rank_f2(np.vstack([cur, row])) > rank_f2(cur):
            cur = np.vstack([cur, row])
            dlb.append(row)
    dlb = np.array(dlb, dtype=np.uint8)
    prof = np.array([int((c1 & z).sum() % 2) for z in dlb])
    odd_classes = [
        eps for eps in range(1, 256)
        if (sum(prof[j] for j in range(8) if (eps >> j) & 1) % 2) == 1
    ]
    reps = []
    for eps in odd_classes:
        rep = np.zeros(150, dtype=np.uint8)
        for j in range(8):
            if (eps >> j) & 1:
                rep ^= dlb[j]
        reps.append(rep)

    for it in range(12):
        Pm = np.array(list(pool.values()), dtype=np.uint8)
        lp, mu = solve_lp(Pm)
        print(f"[iter {it}] pool {Pm.shape[0]}, LP = {lp:.4f}, "
              f"dual sum = {mu.sum():.4f}")
        # column hunt: best mu-mass over all odd classes
        t0 = time.time()
        best = (1e9, None)
        for rep in reps:
            m, w = mu_anneal(rep, gen_rows, mu, rng,
                             steps=4000 if it else 8000)
            if m < best[0]:
                best = (m, w)
        print(f"  best column mass found: {best[0]:.4f} "
              f"(weight {int(best[1].sum())}) [{time.time()-t0:.0f}s]")
        if best[0] >= 1.0 - 1e-6:
            print("  NO violated column: dual cover is (empirically) feasible "
                  f"for the full universe ==> true LP <= {mu.sum():.4f}")
            break
        # add all odd-pairing translates of the found column
        added = 0
        for ga in range(LX):
            for gb in range(LY):
                t = translate1(best[1], (ga, gb))
                if int((t & c1).sum()) % 2 == 1 and t.tobytes() not in pool:
                    pool[t.tobytes()] = t
                    added += 1
        print(f"  added {added} translates")
        if added == 0:
            print("  column already present (annealing stuck); stopping")
            break

    Pm = np.array(list(pool.values()), dtype=np.uint8)
    lp, mu = solve_lp(Pm)
    print(f"FINAL: pool {Pm.shape[0]}, LP = {lp:.4f}")
    np.savez_compressed(
        Path(__file__).resolve().parents[1] / "data/a23/detector_pool_cg.npz",
        P=Pm)


if __name__ == "__main__":
    main()
