"""A23 phase 10: the final-form inequality and its Pareto frontier.

Chain of reductions (verified in a23_transfer_structure.py):
  seam cosets  =  C_graph \\ boundaries,   C_graph = boundaries + {(b,0)},
  b in ker d2 \\ 0 = the 15 translates of the block idempotent e0.
So

  SeamCosetFloor 16  <=>  forall f, b != 0:  |A*f + b| + |B*f| >= 16
                     <=>  forall f:          |A*f + e0| + |B*f| >= 16
                                            (b-orbit = G-orbit of e0).

Here:
  1. verify the finite dictionary: for each kernel basis element K[i],
     seamC K[i] = (b_i + A*f_i | B*f_i) with explicit b_i, f_i  (the
     Lean-side bridge: 4 decide-checks + linearity);
  2. the tight witness split (a*, b*) with a* + b* = 16;
  3. anneal the Pareto frontier: for each a = |A*f + e0|, min |B*f|
     (upper bounds; the SAT UNSAT@14 says the true frontier obeys
     a + b >= 16).
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
    conv_matrix,
    seam_maps,
)


def main() -> None:
    rng = np.random.default_rng(161616)
    MA = conv_matrix(A_SUPP)
    MB = conv_matrix(B_SUPP)
    D2 = np.vstack([MA, MB])
    K = nullspace_f2(D2)

    # e0: recompute (unique idempotent acting as id on block)
    from itertools import product

    e0 = None
    for coeffs in product((0, 1), repeat=4):
        if not any(coeffs):
            continue
        e = np.zeros(75, dtype=np.uint8)
        for c, row in zip(coeffs, K):
            if c:
                e ^= row
        Me = conv_matrix([(j // 15, j % 15) for j in np.flatnonzero(e)])
        if all(np.array_equal((Me @ K[i]) % 2, K[i]) for i in range(4)):
            e0 = e
            break
    assert e0 is not None
    print(f"|e0| = {int(e0.sum())}")

    # --- 1. dictionary: seamC K[i] = (b_i + A f_i | B f_i) ------------------
    print("dictionary checks:")
    dict_out = []
    for i in range(4):
        _, sC = seam_maps(K[i])
        sL, sR = sC[:75].copy(), sC[75:].copy()
        # solve B f = sR ; then b := sL + A f must be in ker d2
        aug = np.hstack([MB, sR[:, None]])
        R, piv = rref_f2(aug)
        assert 75 not in piv, f"B f = sC_R unsolvable for K[{i}]"
        f = np.zeros(75, dtype=np.uint8)
        for r, p in enumerate(piv):
            f[p] = R[r, 75]
        assert np.array_equal((MB @ f) % 2, sR)
        b = (sL + MA @ f) % 2
        in_ker = not ((D2 @ b) % 2).any()
        nz = b.any()
        print(f"  K[{i}]: |f_i| = {int(f.sum())}, b_i in ker d2: {in_ker}, "
              f"b_i != 0: {nz}")
        dict_out.append({
            "kernel_basis_index": i,
            "f_support": sorted(map(int, np.flatnonzero(f))),
            "b_support": sorted(map(int, np.flatnonzero(b))),
        })

    # --- 2. tight witness: anneal min |A f + e0| + |B f| -------------------
    def total(f: np.ndarray) -> tuple[int, int]:
        return int(((MA @ f + e0) % 2).sum()), int(((MB @ f) % 2).sum())

    best = (10**9, None, None)
    frontier: dict[int, int] = {}
    f = np.zeros(75, dtype=np.uint8)
    cw = sum(total(f))
    for t in range(120000):
        j = rng.integers(75)
        cand = f.copy()
        cand[j] ^= 1
        a, b = total(cand)
        nw = a + b
        if nw <= cw or rng.random() < 0.03:
            f, cw = cand, nw
            if nw < best[0]:
                best = (nw, a, b)
            if a <= 40 and (a not in frontier or b < frontier[a]):
                frontier[a] = b
        if t % 40000 == 39999:
            f = rng.integers(0, 2, 75).astype(np.uint8)
            cw = sum(total(f))
    print(f"min |A f + e0| + |B f| found: {best[0]} at split "
          f"(a,b) = ({best[1]},{best[2]})  [SAT says true min = 16]")
    front = sorted(frontier.items())
    print("Pareto frontier samples (a = |Af+e0|, min-found b = |Bf|):")
    print("   ", [(a, b) for a, b in front if a + b <= 22][:30])

    outp = Path(__file__).resolve().parents[1] / "data/a23/final_form.json"
    outp.write_text(json.dumps({
        "e0_support": sorted(map(int, np.flatnonzero(e0))),
        "dictionary": dict_out,
        "best_found": {"total": best[0], "a": best[1], "b": best[2]},
        "frontier_found": {str(a): b for a, b in front},
    }, indent=1))
    print(f"wrote {outp}")


if __name__ == "__main__":
    main()
