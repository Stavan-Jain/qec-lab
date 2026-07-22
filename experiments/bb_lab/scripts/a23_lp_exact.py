"""A23 phase 6b: exact fractional packing LP via scipy/HiGHS.

Run:  uv run --with scipy --project experiments/bb_lab python \
          experiments/bb_lab/scripts/a23_lp_exact.py

max 1.lambda  s.t.  P^T lambda <= 1  (coverage), lambda >= 0.
Also reports the optimal dual (fractional cover) and the support structure
of both optima.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
from scipy.optimize import linprog

sys.path.insert(0, str(Path(__file__).resolve().parent))
from a23_lp_estimate import build_or_load_pool  # noqa: E402


def main() -> None:
    P, c1 = build_or_load_pool()
    n, m = P.shape
    # linprog minimizes: min -1.lambda s.t. A_ub lambda <= b_ub
    res = linprog(
        c=-np.ones(n),
        A_ub=P.T.astype(float),  # (150, n)
        b_ub=np.ones(m),
        bounds=(0, None),
        method="highs",
    )
    assert res.status == 0, res.message
    lp = -res.fun
    lam = res.x
    print(f"LP (fractional packing) over pool {P.shape}: {lp:.4f}")
    used = np.flatnonzero(lam > 1e-6)
    print(f"support size {used.size}; weights used: "
          f"{sorted(set(int(P[i].sum()) for i in used))}")
    vals = sorted(set(np.round(lam[used], 4)))
    print(f"distinct lambda values: {vals[:12]}")
    # dual: fractional cover
    dual = res.ineqlin.marginals  # <= 0
    mu = -dual
    print(f"dual cover value: {mu.sum():.4f}; mu support {np.count_nonzero(mu > 1e-9)}"
          f"; distinct mu vals {sorted(set(np.round(mu[mu > 1e-9], 4)))[:12]}")


if __name__ == "__main__":
    main()
