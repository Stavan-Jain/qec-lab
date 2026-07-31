"""Probe: per-block parity rows on the coarse strengthened instance.

The K = G endpoint of the quotient tower: summing all H_Z rows gives
c_L·|v_L| + c_R·|v_R| ≡ 0 where c_L/c_R are the (constant) column
sums of the two blocks — for |A|, |B| odd this is par(v_L) = par(v_R),
and together with the coset-parity theorem (all classes even) each
block's weight is even *individually*. That conjunction is implied
linear algebra but not a single row of the instance, so we hand it to
CMS as two extra XOR rows and measure the refutation round.

Premises verified per code; applies to every BB code with odd |A|,|B|
— including odd-order groups (bb_90) where no axis deck exists.

Usage: uv run python scripts/block_parity_probe.py [--reps 3]
"""

from __future__ import annotations

import argparse
import time
from pathlib import Path

import numpy as np
import yaml

LAB_ROOT = Path(__file__).resolve().parent.parent

import pycryptosat  # noqa: E402
from pysat.card import CardEnc, EncType  # noqa: E402
from pysat.formula import IDPool  # noqa: E402

from bb_lab.checks import bb_check_matrices  # noqa: E402
from bb_lab.group import ZmZn  # noqa: E402
from bb_lab.poly import Poly  # noqa: E402
from bb_lab.shard_distance import (  # noqa: E402
    _int_to_bits,
    compute_class_action,
)


def refute_round(
    checks, action, w: int, block_parity: bool
) -> float:
    """One coarse case-L refutation at ≤ w (the round that costs
    everything); returns solver seconds. Asserts UNSAT."""
    H_Z, L_Z = checks.H_Z, action.L_Z
    n = H_Z.shape[1]
    N = n // 2
    k = action.k
    reps = action.orbit_reps
    pool = IDPool()
    qv = [pool.id() for _ in range(n)]
    s = pycryptosat.Solver()
    for row in H_Z:
        idx = np.flatnonzero(row)
        if idx.size:
            s.add_xor_clause([qv[i] for i in idx], False)
    a = []
    for L in L_Z:
        aj = pool.id()
        s.add_xor_clause([qv[i] for i in np.flatnonzero(L)] + [aj], False)
        a.append(aj)
    sel = [pool.id() for _ in reps]
    s.add_clause(sel)
    for sl, rep in zip(sel, reps):
        cbits = _int_to_bits(rep, k)
        for j in range(k):
            s.add_clause([-sl, a[j] if cbits[j] else -a[j]])
        T = action.transversal[rep]
        s.add_clause([-sl] + [qv[p] for p in T])
    if block_parity:
        # premises: constant odd column sums per block; all classes even
        cs = H_Z.sum(axis=0) % 2
        assert cs[:N].all() and cs[N:].all(), "block colsums not all odd"
        assert not any(int(r.sum()) % 2 for r in checks.H_X)
        assert not any(int(v.sum()) % 2 for v in action.V.sum(axis=1) % 2)
        s.add_xor_clause([qv[i] for i in range(N)], False)
        s.add_xor_clause([qv[i] for i in range(N, n)], False)
    card = CardEnc.atmost(
        lits=qv, bound=w, vpool=pool, encoding=EncType.seqcounter
    )
    for cl in card.clauses:
        s.add_clause(cl)
    t0 = time.perf_counter()
    sat, _ = s.solve()
    dt = time.perf_counter() - t0
    assert sat is False, "expected UNSAT at w = d(tightened) - stayed SAT?"
    return dt


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--reps", type=int, default=3)
    ap.add_argument(
        "--codes", nargs="*", default=["bb_90_8_10", "bb_108_8_10"]
    )
    args = ap.parse_args()
    table = {
        r["code_id"]: r
        for r in yaml.safe_load(
            (LAB_ROOT / "instances" / "bravyi_table.yaml").read_text()
        )["instances"]
    }
    for name in args.codes:
        row = table[name]
        G = ZmZn(row["group"]["ell"], row["group"]["m"])
        A = Poly.from_string(row["polynomials"]["A"], G)
        B = Poly.from_string(row["polynomials"]["B"], G)
        checks = bb_check_matrices(A, B)
        action = compute_class_action(checks)
        d = row["parameters"]["d"]
        w = d - 1 - ((d - 1) % 2)  # parity-tightened (all classes even)
        for bp in (False, True):
            ts = [
                refute_round(checks, action, w, bp)
                for _ in range(args.reps)
            ]
            print(
                f"{name:<12} w≤{w} block_parity={bp!s:<5} "
                + " ".join(f"{t:6.2f}s" for t in ts),
                flush=True,
            )


if __name__ == "__main__":
    main()
