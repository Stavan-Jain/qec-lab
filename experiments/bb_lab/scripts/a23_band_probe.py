"""A23 phase 8b: y-band region partitions + seam support geometry.

Bands are x-complete strips in y (both blocks).  B's y-degrees {6,10,12}
mean confined dual cycles need tall bands; probe several partitions and
report confined dims + local floors.  Also print the per-block seam
supports (the affine targets for the difference-set route).
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from bb_lab.linalg import nullspace_f2  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
from a23_seam_calibration import (  # noqa: E402
    A_SUPP,
    B_SUPP,
    LX,
    LY,
    conv_matrix,
    seam_maps,
)
from a23_region_probe import confined_space, local_floor_anneal  # noqa: E402


def band_mask(rows: list[int], blocks=(0, 1)) -> np.ndarray:
    m = np.zeros(150, dtype=bool)
    for blk in blocks:
        for a in range(LX):
            for b in rows:
                m[blk * 75 + a * LY + b] = True
    return m


def main() -> None:
    rng = np.random.default_rng(808)
    MA = conv_matrix(A_SUPP)
    MB = conv_matrix(B_SUPP)
    D2 = np.vstack([MA, MB])
    D2T = (D2.T % 2).astype(np.uint8)
    K = nullspace_f2(D2)

    print("=== seam geometry of the 4 kernel basis elements ===")
    for i in range(K.shape[0]):
        _, sC = seam_maps(K[i])
        supL = [(j // LY, j % LY) for j in np.flatnonzero(sC[:75])]
        supR = [(j // LY, j % LY) for j in np.flatnonzero(sC[75:])]
        rowsL = sorted(set(b for _, b in supL))
        rowsR = sorted(set(b for _, b in supR))
        print(f"  K[{i}]: |sC_L|={len(supL)} rows {rowsL}; "
              f"|sC_R|={len(supR)} rows {rowsR}")

    zeta = K[0]
    _, c1 = seam_maps(zeta)

    partitions = {
        "7+7+1 rows (0-6 | 7-13 | 14)": [list(range(0, 7)), list(range(7, 14)), [14]],
        "8+7 rows (0-7 | 8-14)": [list(range(0, 8)), list(range(8, 15))],
        "5+5+5 rows": [list(range(0, 5)), list(range(5, 10)), list(range(10, 15))],
        "10+5 rows": [list(range(0, 10)), list(range(10, 15))],
        "12+3 rows": [list(range(0, 12)), list(range(12, 15))],
        "13+2 rows": [list(range(0, 13)), list(range(13, 15))],
        "one 13-row band only (0-12)": [list(range(0, 13))],
        "one 12-row band only (0-11)": [list(range(0, 12))],
    }
    for name, rowsets in partitions.items():
        floors = []
        descr = []
        for rows in rowsets:
            R = band_mask(rows)
            C = confined_space(D2T, R)
            fl = local_floor_anneal(C, c1, R, rng, tries=25, steps=3000)
            floors.append(max(fl, 0))
            descr.append(f"rows{rows[0]}-{rows[-1]}: dim {C.shape[0]}, floor<={fl}")
        print(f"[{name}] " + " | ".join(descr) + f"  SUM<= {sum(floors)}")


if __name__ == "__main__":
    main()
