"""A23 phase 8: region-confined parity structure (beating the packing LP).

Partition the 150 qubits into regions; for each region R the constraints on
w /\\ R are the parities <w, z> = <c1, z> for every dual cycle z with
supp z inside R.  The local floor = min-size subset of R satisfying all
confined parities; disjoint regions sum.  Parity information (odd
intersections, even-pairing cycles) is exactly what the packing LP ignores.

Probes:
  1. the disjoint S-orbit saturation pairs (40+40 points) + the 70-point rest;
  2. per region: dimension of the confined dual-cycle space, pairing values,
     and the local floor (annealed upper bound + exact-ish lower probes).
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
    conv_matrix,
    seam_maps,
    translate1,
)
from a23_targeted_packing import S_GENS  # noqa: E402


def confined_space(D2T: np.ndarray, region: np.ndarray) -> np.ndarray:
    """Basis (rows) of dual cycles supported inside `region` (bool mask)."""
    from bb_lab.linalg import nullspace_f2

    cols = np.flatnonzero(region)
    sub = D2T[:, cols]  # (75, |R|)
    ns = nullspace_f2(sub)  # rows in R-coordinates
    out = np.zeros((ns.shape[0], 150), dtype=np.uint8)
    out[:, cols] = ns
    return out


def local_floor_anneal(
    confined: np.ndarray,
    c1: np.ndarray,
    region: np.ndarray,
    rng: np.random.Generator,
    tries: int = 40,
    steps: int = 4000,
) -> int:
    """Upper bound on the local floor: min |T| over T inside region with
    <T, z> = <c1, z> for all confined z.  Solve the affine system and
    minimize weight by annealing over the solution space."""
    cols = np.flatnonzero(region)
    nR = cols.size
    if confined.shape[0] == 0:
        return 0
    Asys = confined[:, cols]  # (m, nR)
    b = np.array([int((c1 & z).sum() % 2) for z in confined], dtype=np.uint8)
    # solve Asys T = b
    aug = np.hstack([Asys, b[:, None]])
    R, piv = rref_f2(aug)
    if any(p == nR for p in piv):
        return -1  # inconsistent: impossible (c1 restricted not compatible?)
    # particular solution
    T0 = np.zeros(nR, dtype=np.uint8)
    for r, p in enumerate(piv):
        T0[p] = R[r, nR]
    ker = nullspace_f2(Asys)
    if ker.shape[0] == 0:
        return int(T0.sum())
    best = int(T0.sum())
    for _ in range(tries):
        w = T0.copy()
        if rng.random() < 0.8:
            for i in rng.choice(ker.shape[0], size=rng.integers(0, min(4, ker.shape[0]) + 1), replace=False):
                w ^= ker[i]
        cw = int(w.sum())
        for t in range(steps):
            i = rng.integers(ker.shape[0])
            cand = w ^ ker[i]
            nw = int(cand.sum())
            if nw <= cw or rng.random() < 0.02:
                w, cw = cand, nw
                best = min(best, cw)
    return best


def main() -> None:
    rng = np.random.default_rng(555)
    MA = conv_matrix(A_SUPP)
    MB = conv_matrix(B_SUPP)
    D2 = np.vstack([MA, MB])
    D2T = (D2.T % 2).astype(np.uint8)
    K = nullspace_f2(D2)
    zeta = K[0]
    _, c1 = seam_maps(zeta)

    import json

    census = json.loads(
        (Path(__file__).resolve().parents[1] / "data/a23/w8_census.json").read_text()
    )
    dets8 = []
    for sup in census["detector_supports_odd_c1"]:
        v = np.zeros(150, dtype=np.uint8)
        v[sup] = 1
        dets8.append(v)

    # rebuild the 8 S-orbits and their saturations
    orb_of: dict[bytes, int] = {}
    orbit_reps: list[np.ndarray] = []
    for v in dets8:
        if v.tobytes() in orb_of:
            continue
        oid = len(orbit_reps)
        for s in [(0, 0)] + S_GENS:
            orb_of[translate1(v, s).tobytes()] = oid
        orbit_reps.append(v)
    sats = []
    for rep in orbit_reps:
        sat = rep.copy()
        for s in S_GENS:
            sat |= translate1(rep, s)
        sats.append(sat)
    n_orb = len(orbit_reps)
    print(f"S-orbits: {n_orb}")
    disj_pairs = [
        (i, j)
        for i in range(n_orb)
        for j in range(i + 1, n_orb)
        if not (sats[i] & sats[j]).any()
    ]
    print(f"disjoint saturation pairs: {disj_pairs}")

    for (i, j) in disj_pairs:
        Ri = sats[i].astype(bool)
        Rj = sats[j].astype(bool)
        Rrest = ~(Ri | Rj)
        print(f"--- partition via orbits ({i},{j}): |R|={Ri.sum()},{Rj.sum()},"
              f"{Rrest.sum()} ---")
        floors = []
        for name, R in (("sat_i", Ri), ("sat_j", Rj), ("rest", Rrest)):
            C = confined_space(D2T, R)
            odd = [z for z in C if int((z & c1).sum()) % 2 == 1]
            fl = local_floor_anneal(C, c1, R, rng)
            floors.append(fl)
            print(f"  {name}: confined dim {C.shape[0]}, "
                  f"odd-pairing among basis {len(odd)}, local floor <= {fl}")
        print(f"  SUM of local floors (upper bounds) = {sum(floors)} "
              f"(need >= 15 for the certificate; >= 16 without parity)")


if __name__ == "__main__":
    main()
