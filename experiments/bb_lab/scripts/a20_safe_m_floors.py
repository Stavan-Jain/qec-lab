"""A20: the lift-aware seam floor — restoring >= 20 after the SAT@18 finding.

`SeamCosetFloor 20` is FALSE as stated (class 0x1 contains a weight-18
element). But the safe-sector conclusion only needs: every cover logical v
with p(v) in a seam coset has |v| >= 20. By the slice identity
|v| = |p(v)| + 2*(off-support count), it suffices that every coset element
w with |w| = 20 - 2j (j >= 1) admits no cover cycle over it with fewer than
j off-support doubly-occupied fibers. Since the cosets avoid boundaries
(delta2-injectivity, calibrated), ANY cover cycle over such w is
automatically a nontrivial logical — no signature clause needed.

Per seam class (orbit rep zeta):
  Phase A (coset census): enumerate ALL elements w = seamC zeta + d2 f with
    |w| <= 19 (SAT + blocking; weights are even by parity, so strata are
    18, 16, ...). Stabilizer-orbit blocking (translations fixing the class).
  Phase B (lift floors): for each found w with |w| = 20 - 2j, the
    fiber-pinned query: cover cycle v (ker H_Z(Y8) XORs), p(v) = w
    (144 fiber XORs), off-support cardinality <= 2(j-1).
    UNSAT  ==> every logical over w weighs >= 20.
    SAT    ==> an explicit sub-20 logical of Y8 (would refute d = 20).
  Verdict per class: census exhausted + all lift queries UNSAT
    ==> lift-aware seam floor at 20 CERTIFIED for that class.

Usage (from experiments/bb_lab):
    uv run scripts/a20_safe_m_floors.py [--class-mask 1] [--census-bound 19]
Results: data/a20/safe_m_floors_<mask>.jsonl; log stdout.
"""
import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np

LAB = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(LAB / "src"))

import pycryptosat
from pysat.card import CardEnc, EncType
from pysat.formula import IDPool

from bb_lab.checks import bb_check_matrices
from bb_lab.group import AbelianGroup
from bb_lab.linalg import nullspace_f2
from bb_lab.poly import Poly

OUT = LAB / "data" / "a20"
LX, LY, CLY = 18, 4, 8
N2, N1 = 72, 144

A4S = [(0, 0), (1, 0), (14, 1)]
B4S = [(0, 0), (1, 2), (2, 3)]
A8S = [(0, 0), (1, 4), (14, 1)]
B8S = [(0, 0), (1, 2), (2, 7)]
Y8 = {"frame": (18, 8), "A": "1 + x*y^4 + x^14*y", "B": "1 + x*y^2 + x^2*y^7"}


def conv_matrix(supp):
    M = np.zeros((N2, N2), dtype=np.uint8)
    for ga in range(LX):
        for gb in range(LY):
            g = (ga % LX) * LY + (gb % LY)
            for (sa, sb) in supp:
                M[((ga + sa) % LX) * LY + ((gb + sb) % LY), g] ^= 1
    return M


def conv_cover(supp, f):
    out = np.zeros_like(f)
    for (sa, sb) in supp:
        out ^= np.roll(np.roll(f, sa, axis=0), sb, axis=1)
    return out


def seam_c(zeta):
    lift = np.zeros((LX, CLY), dtype=np.uint8)
    lift[:, :LY] = zeta.reshape(LX, LY)
    cA = conv_cover(A8S, lift)
    cB = conv_cover(B8S, lift)
    return np.concatenate([cA[:, LY:].reshape(-1), cB[:, LY:].reshape(-1)])


def translate1(v, g):
    out = np.empty_like(v)
    for blk in range(2):
        arr = v[blk * N2:(blk + 1) * N2].reshape(LX, LY)
        out[blk * N2:(blk + 1) * N2] = np.roll(
            np.roll(arr, g[0], axis=0), g[1], axis=1).reshape(-1)
    return out


def translate2(z, g):
    return np.roll(np.roll(z.reshape(LX, LY), g[0], axis=0),
                   g[1], axis=1).reshape(-1)


def lift_query(w144, j, HZ8, fiber):
    """Cover cycle over w with off-support count <= j-1? -> (verdict, data)."""
    pool = IDPool()
    v = [pool.id() for _ in range(2 * 144)]
    solver = pycryptosat.Solver()
    for row in HZ8:
        idx = np.flatnonzero(row)
        if idx.size:
            solver.add_xor_clause([v[i] for i in idx], False)
    off_vars = []
    for q4, (q0, q1) in enumerate(fiber):
        solver.add_xor_clause([v[q0], v[q1]], bool(w144[q4]))
        if not w144[q4]:
            off_vars.extend((v[q0], v[q1]))
    bound = 2 * (j - 1)
    if bound == 0:
        for lit in off_vars:
            solver.add_clause([-lit])
    else:
        card = CardEnc.atmost(lits=off_vars, bound=bound, vpool=pool,
                              encoding=EncType.seqcounter)
        for cl in card.clauses:
            solver.add_clause(cl)
    sat, model = solver.solve()
    if not sat:
        return "UNSAT", None
    vv = np.array([1 if model[x] else 0 for x in v], dtype=np.uint8)
    return "SAT", np.flatnonzero(vv).tolist()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--class-mask", type=int, default=1)
    ap.add_argument("--census-bound", type=int, default=19)
    args = ap.parse_args()

    MA, MB = conv_matrix(A4S), conv_matrix(B4S)
    D2 = np.vstack([MA, MB])
    K = nullspace_f2(D2)
    zeta = np.zeros(N2, dtype=np.uint8)
    for i in range(K.shape[0]):
        if (args.class_mask >> i) & 1:
            zeta ^= K[i]
    sc = seam_c(zeta)

    # stabilizer of the class: translations g with g.zeta = zeta
    stab = [(ta, tb) for ta in range(LX) for tb in range(LY)
            if (translate2(zeta, (ta, tb)) == zeta).all()]
    print(f"class {args.class_mask:#x}: |seamC| = {int(sc.sum())}, "
          f"class-stabilizer order {len(stab)}", flush=True)

    # cover context for lift queries
    G8 = AbelianGroup(Y8["frame"])
    c8 = bb_check_matrices(Poly.from_string(Y8["A"], G8),
                           Poly.from_string(Y8["B"], G8))
    HZ8 = c8.H_Z % 2
    N8 = 144
    fiber = []
    for blk in (0, 1):
        for a in range(LX):
            for b in range(LY):
                fiber.append((blk * N8 + a * CLY + b,
                              blk * N8 + a * CLY + b + 4))

    # Phase A + B interleaved: census with blocking; lift-check each element
    res_path = OUT / f"safe_m_floors_{args.class_mask:#x}.jsonl"
    log = res_path.open("a")
    pool = IDPool()
    w = [pool.id() for _ in range(N1)]
    f = [pool.id() for _ in range(N2)]
    solver = pycryptosat.Solver()
    for jrow in range(N1):
        idx = np.flatnonzero(D2[jrow])
        solver.add_xor_clause([w[jrow]] + [f[i] for i in idx], bool(sc[jrow]))
    card = CardEnc.atmost(lits=w, bound=args.census_bound, vpool=pool,
                          encoding=EncType.seqcounter)
    for cl in card.clauses:
        solver.add_clause(cl)

    # resume: reblock previously-found elements (and count them)
    n_elts, n_lift_unsat, sat_hits, t0 = 0, 0, [], time.time()
    if res_path.exists():
        for line in res_path.read_text().splitlines():
            row = json.loads(line)
            wv = np.zeros(N1, dtype=np.uint8)
            wv[row["w_support"]] = 1
            for g in stab:
                tv = translate1(wv, g)
                solver.add_clause(
                    [-w[q] if tv[q] else w[q] for q in range(N1)])
            n_elts += 1
            if row["lift"] == "UNSAT":
                n_lift_unsat += 1
            else:
                sat_hits.append(row)
        print(f"resumed: {n_elts} elements reblocked "
              f"({n_lift_unsat} lift-UNSAT)", flush=True)
    while True:
        sat, model = solver.solve()
        if not sat:
            break
        wv = np.array([1 if model[x] else 0 for x in w], dtype=np.uint8)
        wt = int(wv.sum())
        assert wt % 2 == 0 and wt <= args.census_bound
        j = (20 - wt + 1) // 2
        verdict, witness = lift_query(wv, j, HZ8, fiber)
        n_elts += 1
        row = {"elt": n_elts, "w": wt, "j": j, "lift": verdict,
               "w_support": np.flatnonzero(wv).tolist()}
        if verdict == "SAT":
            row["sub20_logical"] = witness
            sat_hits.append(row)
            print(f"!! elt {n_elts} wt {wt}: LIFT SAT — sub-20 logical?!",
                  flush=True)
        else:
            n_lift_unsat += 1
        log.write(json.dumps(row) + "\n")
        log.flush()
        # block the stabilizer orbit of this element
        for g in stab:
            tv = translate1(wv, g)
            solver.add_clause([-w[q] if tv[q] else w[q] for q in range(N1)])
        if n_elts % 10 == 0:
            print(f"  {n_elts} elements censused "
                  f"({n_lift_unsat} lift-UNSAT, {time.time()-t0:.0f}s)",
                  flush=True)
    dt = time.time() - t0
    print(f"census EXHAUSTED: {n_elts} stabilizer-orbit elements with "
          f"|w| <= {args.census_bound} ({dt:.0f}s)", flush=True)
    if not sat_hits:
        print(f"LIFT-AWARE SEAM FLOOR AT 20 CERTIFIED for class "
              f"{args.class_mask:#x}: every coset element is >= 20 or all "
              f"its cover logicals weigh >= 20.", flush=True)
    log.close()


if __name__ == "__main__":
    main()
