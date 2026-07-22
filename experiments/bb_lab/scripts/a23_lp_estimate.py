"""A23 phase 6: estimate the fractional packing LP value over the detector pool.

LP:  max sum lambda_z  s.t.  per-qubit coverage <= 1, lambda >= 0
Dual: min sum mu_j     s.t.  sum_{j in z} mu_j >= 1 for all z, mu >= 0.

MWU on the dual: maintain mu, repeatedly find the most-violated candidate
(min mu-mass), boost its qubits; the running normalized average gives a
feasible-in-the-limit dual estimate.  We report both a primal lower bound
(best integral/greedy fractional found) and the dual upper bound
(scaled feasible mu), bracketing the LP.

Pool = complete w8 odd-pairing universe (40) + annealed w10/w12 translates.
Saved to data/a23/detector_pool.npz for reuse.
"""

from __future__ import annotations

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
from a23_targeted_packing import anneal_class  # noqa: E402


def build_or_load_pool() -> tuple[np.ndarray, np.ndarray]:
    """Returns (P, c1): P is (n_cand, 150) uint8."""
    import json

    poolp = Path(__file__).resolve().parents[1] / "data/a23/detector_pool.npz"
    MA = conv_matrix(A_SUPP)
    MB = conv_matrix(B_SUPP)
    D2 = np.vstack([MA, MB])
    D1 = np.hstack([MB, MA])
    K = nullspace_f2(D2)
    _, c1 = seam_maps(K[0])
    if poolp.exists():
        P = np.load(poolp)["P"]
        print(f"loaded pool {P.shape}")
        return P, c1

    rng = np.random.default_rng(1234)
    census = json.loads(
        (Path(__file__).resolve().parents[1] / "data/a23/w8_census.json").read_text()
    )
    pool: dict[bytes, np.ndarray] = {}
    for sup in census["detector_supports_odd_c1"]:
        v = np.zeros(150, dtype=np.uint8)
        v[sup] = 1
        pool[v.tobytes()] = v

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
    gen_rows = D1.astype(np.uint8)
    shapes: dict[bytes, np.ndarray] = {}
    t0 = time.time()
    for eps in odd_classes:
        rep = np.zeros(150, dtype=np.uint8)
        for j in range(8):
            if (eps >> j) & 1:
                rep ^= dlb[j]
        shapes.update(anneal_class(rep, gen_rows, rng, steps=7000, wmax=12))
    print(f"annealed {len(shapes)} shapes [{time.time()-t0:.0f}s]")
    for s in shapes.values():
        for ga in range(LX):
            for gb in range(LY):
                t = translate1(s, (ga, gb))
                if int((t & c1).sum()) % 2 == 1:
                    pool[t.tobytes()] = t
    P = np.array(list(pool.values()), dtype=np.uint8)
    poolp.parent.mkdir(exist_ok=True)
    np.savez_compressed(poolp, P=P)
    print(f"pool {P.shape} saved")
    return P, c1


def mwu_dual(P: np.ndarray, iters: int = 20000) -> float:
    """Upper bound on LP via a feasible dual point from MWU averaging."""
    n, m = P.shape
    Pf = P.astype(np.float64)
    mu = np.ones(m) / 8.0
    avg = np.zeros(m)
    eta = 0.05
    for t in range(iters):
        loads = Pf @ mu  # mass per candidate
        i = int(np.argmin(loads))
        if loads[i] >= 1.0 - 1e-9:
            # mu/max(...) is feasible now
            return float(mu.sum() / loads[i]) if loads[i] < 1 else float(mu.sum())
        # boost qubits of the violated candidate
        mu[P[i] > 0] *= (1 + eta)
        mu /= (1 + eta * 8 / m)  # gentle renormalization to control growth
        avg += mu
    avg /= iters
    scale = (Pf @ avg).min()
    return float(avg.sum() / scale)


def greedy_fractional(P: np.ndarray, rounds: int = 400,
                      rng: np.random.Generator | None = None) -> float:
    """Primal lower bound: randomized greedy fractional packing
    (uniform on maximal disjoint families, averaged)."""
    rng = rng or np.random.default_rng(0)
    n, m = P.shape
    lam = np.zeros(n)
    cover = np.zeros(m)
    Pf = P.astype(np.float64)
    # simple FW-ish: repeatedly add eps of the candidate with most headroom
    eps = 0.02
    for _ in range(20000):
        head = 1.0 - (cover[None, :] * (Pf > 0)).max(axis=1) if False else None
        maxcov = np.max(np.where(Pf > 0, cover[None, :], 0.0), axis=1)
        room = 1.0 - maxcov
        i = int(np.argmax(room))
        if room[i] <= 1e-9:
            break
        step = min(eps, room[i])
        lam[i] += step
        cover += step * Pf[i]
    return float(lam.sum())


def main() -> None:
    P, c1 = build_or_load_pool()
    wt = P.sum(axis=1)
    hist: dict[int, int] = {}
    for w in wt:
        hist[int(w)] = hist.get(int(w), 0) + 1
    print(f"pool weights: {dict(sorted(hist.items()))}")

    t0 = time.time()
    lb = greedy_fractional(P)
    print(f"primal (greedy fractional) lower bound: {lb:.2f}  [{time.time()-t0:.0f}s]")
    t0 = time.time()
    ub = mwu_dual(P, iters=30000)
    print(f"dual (feasible mu) upper bound:        {ub:.2f}  [{time.time()-t0:.0f}s]")
    print(f"LP value in [{lb:.2f}, {ub:.2f}]; need > 14 for the "
          f"fractional certificate route")


if __name__ == "__main__":
    main()
