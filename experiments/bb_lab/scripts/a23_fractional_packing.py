"""A23 phase 3: FRACTIONAL detector packing via PB-SAT.

Mechanism upgrade: detectors may overlap.  If z_1..z_N are dual cycles
(with multiplicities m_i <= M), each pairing odd with the seam offset c1,
and every qubit lies in at most M detectors counted with multiplicity,
then for every coset element w:

    sum_i m_i  <=  sum_i m_i |w /\ z_i|  =  sum_{j in w} cover(j)  <=  M |w|

so |w| >= (sum m_i)/M.  Target sum m_i >= 14M + 1  ==>  |w| >= 15, and
parity (coset weights all even) lifts to the tight floor 16.

Encoding: M binary copies per candidate; per-qubit atmost-M; global
atleast-(14M+1).  Solved with pysat (short runs only).
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
    flat,
    seam_maps,
    translate1,
)
from a23_detector_packing import reflect_swap  # noqa: E402
from a23_targeted_packing import anneal_class  # noqa: E402


def build_pool(rng: np.random.Generator, wmax: int = 12,
               steps: int = 9000) -> tuple[list[np.ndarray], np.ndarray, dict]:
    MA = conv_matrix(A_SUPP)
    MB = conv_matrix(B_SUPP)
    D2 = np.vstack([MA, MB])
    D1 = np.hstack([MB, MA])
    K = nullspace_f2(D2)
    zeta = K[0]
    _, c1 = seam_maps(zeta)

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

    data_p = Path(__file__).resolve().parents[1] / "data/a17/f2a6_z5z30_lean_data.json"
    data = json.loads(data_p.read_text())
    wst = np.zeros(150, dtype=np.uint8)
    for (a, b) in data["ustar_left"]:
        wst[flat(a, b)] = 1
    for (a, b) in data["ustar_right"]:
        wst[75 + flat(a, b)] = 1
    zstar = reflect_swap(wst)
    if int((zstar & c1).sum()) % 2 == 1:
        shapes[zstar.tobytes()] = zstar

    print(f"annealing {len(odd_classes)} odd classes to wmax={wmax}...")
    t0 = time.time()
    for idx, eps in enumerate(odd_classes):
        rep = np.zeros(150, dtype=np.uint8)
        for j in range(8):
            if (eps >> j) & 1:
                rep ^= dual_log_basis[j]
        shapes.update(anneal_class(rep, gen_rows, rng, steps=steps, wmax=wmax))
    print(f"  {len(shapes)} shapes in {time.time()-t0:.0f}s")

    # all odd-pairing translates
    cands: dict[bytes, np.ndarray] = {}
    for s in shapes.values():
        for ga in range(LX):
            for gb in range(LY):
                t = translate1(s, (ga, gb))
                if int((t & c1).sum()) % 2 == 1:
                    cands[t.tobytes()] = t
    cand_list = list(cands.values())
    wt_hist: dict[int, int] = {}
    for c in cand_list:
        wt_hist[int(c.sum())] = wt_hist.get(int(c.sum()), 0) + 1
    print(f"candidate translates: {len(cand_list)}, weights {dict(sorted(wt_hist.items()))}")
    return cand_list, c1, {"D2": D2, "D1": D1, "K": K}


def try_pb(cand_list: list[np.ndarray], M: int, target: int,
           timeout_s: int = 180) -> list[int] | None:
    """Feasibility: multiplicities m_z <= M, coverage <= M, sum m_z >= target.
    Returns chosen list of candidate indices with repetition, or None."""
    from pysat.card import CardEnc, EncType
    from pysat.formula import IDPool
    from pysat.solvers import Cadical195

    pool = IDPool()
    n = len(cand_list)
    # binary copy variables x[z][c], c < M, with symmetry chain
    xv = [[pool.id(f"x_{z}_{c}") for c in range(M)] for z in range(n)]
    clauses: list[list[int]] = []
    for z in range(n):
        for c in range(1, M):
            clauses.append([-xv[z][c], xv[z][c - 1]])  # x_c -> x_{c-1}

    incid: list[list[int]] = [[] for _ in range(150)]
    for z, v in enumerate(cand_list):
        for j in np.flatnonzero(v):
            incid[j].extend(xv[z])

    cnf_all = list(clauses)
    for j in range(150):
        if len(incid[j]) > M:
            enc = CardEnc.atmost(lits=incid[j], bound=M, vpool=pool,
                                 encoding=EncType.seqcounter)
            cnf_all.extend(enc.clauses)
    all_x = [x for row in xv for x in row]
    enc = CardEnc.atleast(lits=all_x, bound=target, vpool=pool,
                          encoding=EncType.seqcounter)
    cnf_all.extend(enc.clauses)

    print(f"  PB model: {pool.top} vars, {len(cnf_all)} clauses; "
          f"M={M}, target={target}")
    with Cadical195(bootstrap_with=cnf_all) as s:
        t0 = time.time()
        ok = s.solve()
        dt = time.time() - t0
        print(f"  solve: {'SAT' if ok else 'UNSAT'} in {dt:.1f}s")
        if not ok:
            return None
        model = set(l for l in s.get_model() if l > 0)
        chosen: list[int] = []
        for z in range(n):
            m_z = sum(1 for c in range(M) if xv[z][c] in model)
            chosen.extend([z] * m_z)
        return chosen


def main() -> None:
    rng = np.random.default_rng(777)
    cand_list, c1, ctx = build_pool(rng, wmax=12, steps=8000)
    if len(cand_list) > 2500:
        order = np.argsort([int(c.sum()) for c in cand_list])
        cand_list = [cand_list[i] for i in order[:2500]]
        print(f"capped pool to {len(cand_list)} lightest")

    for M in (2, 3, 4):
        target = 14 * M + 1
        print(f"--- fractional packing M={M}, need sum m >= {target} ---")
        chosen = try_pb(cand_list, M, target)
        if chosen is None:
            continue
        # verify
        cover = np.zeros(150, dtype=int)
        total = 0
        for z in chosen:
            cover += cand_list[z]
            total += 1
            assert int((cand_list[z] & c1).sum()) % 2 == 1
        assert cover.max() <= M, f"coverage {cover.max()} > {M}"
        print(f"  VERIFIED: {total} detector-copies (>= {target}), "
              f"max coverage {cover.max()} <= {M}  ==> floor >= "
              f"ceil({total}/{M}) = {-(-total // M)} (+parity -> 16)")
        out = {
            "M": M,
            "total": total,
            "class_rep": "K[0]",
            "detector_indices_with_multiplicity": sorted(chosen),
            "detectors": [sorted(map(int, np.flatnonzero(cand_list[z])))
                          for z in sorted(set(chosen))],
            "multiplicity": {int(z): chosen.count(z) for z in sorted(set(chosen))},
        }
        outp = Path(__file__).resolve().parents[1] / "data/a23/fractional_rep1.json"
        outp.parent.mkdir(exist_ok=True)
        outp.write_text(json.dumps(out, indent=1))
        print(f"  wrote {outp}")
        return
    print("all M in (2,3,4) UNSAT on this pool -- fractional mechanism "
          "insufficient at these weights")


if __name__ == "__main__":
    main()
