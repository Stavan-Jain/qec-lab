"""A20: SeamCosetFloor 20 — calibration + per-coset probes for the Y-tower.

Port of the A23 seam calibration to Y4 = (18,4) -> Y8 = (18,8), followed by
the floor probes. Convention (QECLean BBChainComplex/BBCover, as in
a23_seam_calibration.py):

  conv P f (h) = sum_x P(x) f(h-x);   d2 f = (A*f | B*f);  d1 c = B*c_L + A*c_R
  sec (a,b) = (a,b) with y in [0,4);  deckS = (0,4)
  seamC zeta (h,j) = (P_j^cover * lift zeta)(h + (0,4))   [sheet-1 slab]

NOTE the non-literal lift: Y8's polynomials are NOT Y4's with the same
exponents (A8 = 1 + x y^4 + x^14 y reduces to A4 = 1 + x + x^14 y mod y^4),
so seamC must convolve with the ACTUAL Y8 supports on (18,8) and split
sheets — the general XDoubleCoverData shape.

Calibration outputs: ker d2(Y4) structure (dim, weights, G-orbits),
|seamC zeta| per class, seamC is a d1-cycle (assert), delta2-injectivity
(seamC zeta not in im d2 for zeta != 0).

Probe per orbit rep: min weight over the coset seamC zeta + im d2(Y4)
via CryptoMiniSat (w = seamC + D2 f as XORs, CardEnc atmost 19).
UNSAT ==> every coset element has weight >= 20: the floor certified for
that class. Budget-limited (confl_limit); budget exhaustion reported as
BUDGET, never as evidence.

Usage (from experiments/bb_lab):
    uv run scripts/a20_seam_floor.py [--probe-bound 19] [--confl-limit 60000000]
        [--jobs 5] [--calibrate-only]
Results: data/a20/seam_floor_results.jsonl (resumable), log to stdout.
"""
import argparse
import json
import multiprocessing as mp
import sys
import time
from pathlib import Path

import numpy as np

LAB = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(LAB / "src"))

import pycryptosat
from pysat.card import CardEnc, EncType
from pysat.formula import IDPool

from bb_lab.linalg import nullspace_f2, rank_f2

OUT = LAB / "data" / "a20"
RESULTS = OUT / "seam_floor_results.jsonl"

LX, LY, CLY = 18, 4, 8
N2, N1 = LX * LY, 2 * LX * LY               # 72, 144

A4 = [(0, 0), (1, 0), (14, 1)]              # 1 + x + x^14 y      (base)
B4 = [(0, 0), (1, 2), (2, 3)]               # 1 + x y^2 + x^2 y^3 (base)
A8 = [(0, 0), (1, 4), (14, 1)]              # 1 + x y^4 + x^14 y  (cover)
B8 = [(0, 0), (1, 2), (2, 7)]               # 1 + x y^2 + x^2 y^7 (cover)


def conv_matrix(supp):
    M = np.zeros((N2, N2), dtype=np.uint8)
    for ga in range(LX):
        for gb in range(LY):
            g = ga * LY + gb
            for (sa, sb) in supp:
                M[((ga + sa) % LX) * LY + ((gb + sb) % LY), g] ^= 1
    return M


def conv_cover(supp, f):
    out = np.zeros_like(f)
    for (sa, sb) in supp:
        out ^= np.roll(np.roll(f, sa, axis=0), sb, axis=1)
    return out


def translate2(v, g):
    return np.roll(np.roll(v.reshape(LX, LY), g[0], axis=0),
                   g[1], axis=1).reshape(-1)


def seam_c(zeta):
    """seamC zeta as a C1(base) vector (length 144), Y8 polynomials."""
    lift = np.zeros((LX, CLY), dtype=np.uint8)
    lift[:, :LY] = zeta.reshape(LX, LY)
    convA = conv_cover(A8, lift)
    convB = conv_cover(B8, lift)
    return np.concatenate([convA[:, LY:].reshape(-1), convB[:, LY:].reshape(-1)])


def probe(task):
    """(mask, seamC_vec, D2, bound, confl) -> result dict."""
    mask, sc, D2, bound, confl = task
    t0 = time.time()
    pool = IDPool()
    w = [pool.id() for _ in range(N1)]
    f = [pool.id() for _ in range(N2)]
    solver = pycryptosat.Solver(confl_limit=confl)
    for j in range(N1):
        idx = np.flatnonzero(D2[j])
        solver.add_xor_clause([w[j]] + [f[i] for i in idx], bool(sc[j]))
    card = CardEnc.atmost(lits=w, bound=bound, vpool=pool,
                          encoding=EncType.seqcounter)
    for cl in card.clauses:
        solver.add_clause(cl)
    sat, model = solver.solve()
    out = {"class_mask": mask, "seam_wt": int(sc.sum()), "bound": bound,
           "secs": round(time.time() - t0, 1)}
    if sat is None:
        out["verdict"] = "BUDGET"
    elif not sat:
        out["verdict"] = "UNSAT"            # floor >= bound+1 certified
    else:
        vv = np.array([1 if model[x] else 0 for x in w], dtype=np.uint8)
        out["verdict"] = "SAT"
        out["witness_weight"] = int(vv.sum())
        out["witness"] = np.flatnonzero(vv).tolist()
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--probe-bound", type=int, default=19)
    ap.add_argument("--confl-limit", type=int, default=60_000_000)
    ap.add_argument("--jobs", type=int, default=5)
    ap.add_argument("--calibrate-only", action="store_true")
    args = ap.parse_args()

    MA, MB = conv_matrix(A4), conv_matrix(B4)
    D2 = np.vstack([MA, MB])                # (144, 72)
    D1 = np.hstack([MB, MA])                # (72, 144)
    assert not ((D1 @ D2) % 2).any(), "d1 . d2 != 0 — convention bug"

    K = nullspace_f2(D2)
    dimk = K.shape[0]
    print(f"[1] dim ker d2(Y4) = {dimk} ({2**dimk - 1} nonzero classes)")

    elts = {}
    for mask in range(1, 2**dimk):
        z = np.zeros(N2, dtype=np.uint8)
        for i in range(dimk):
            if (mask >> i) & 1:
                z ^= K[i]
        elts[mask] = z

    # G-orbit collapse on the kernel classes
    key = {mask: min(tuple(translate2(z, (ta, tb)))
                     for ta in range(LX) for tb in range(LY))
           for mask, z in elts.items()}
    reps, orbit_of_key = [], {}
    for mask in sorted(elts):
        if key[mask] not in orbit_of_key:
            orbit_of_key[key[mask]] = mask
            reps.append(mask)
    print(f"[2] G-orbits of nonzero kernel classes: {len(reps)} reps "
          f"{reps} (orbit sizes "
          f"{[sum(1 for m in elts if key[m] == key[r]) for r in reps]})")

    rank_d2 = rank_f2(D2)
    aug = {}
    for r in reps:
        sc = seam_c(elts[r])
        assert not ((D1 @ sc) % 2).any(), "seamC not a d1-cycle?!"
        inj = rank_f2(np.vstack([D2.T, sc])) > rank_d2
        aug[r] = (sc, inj)
        print(f"[3] class {r:#x}: zeta wt {int(elts[r].sum())}, "
              f"|seamC| = {int(sc.sum())}, "
              f"delta2-injective: {inj}")
    if args.calibrate_only:
        return

    done = set()
    if RESULTS.exists():
        for line in RESULTS.read_text().splitlines():
            done.add(json.loads(line)["class_mask"])
    tasks = [(r, aug[r][0], D2, args.probe_bound, args.confl_limit)
             for r in reps if r not in done]
    print(f"\n[probe] {len(tasks)} coset queries at w <= {args.probe_bound}, "
          f"confl_limit {args.confl_limit:.0e}", flush=True)
    log = RESULTS.open("a")
    ctx = mp.get_context("spawn")
    with ctx.Pool(min(args.jobs, max(len(tasks), 1))) as pool_:
        for out in pool_.imap_unordered(probe, tasks):
            log.write(json.dumps(out) + "\n")
            log.flush()
            print(f"  class {out['class_mask']:#x} |seamC|={out['seam_wt']}: "
                  f"{out['verdict']}"
                  + (f" wt {out['witness_weight']}" if out["verdict"] == "SAT"
                     else "") + f" ({out['secs']}s)", flush=True)
    log.close()
    print("done — a class's floor is certified ONLY if its verdict above "
          "is UNSAT (SAT = light element found; BUDGET = no evidence).",
          flush=True)


if __name__ == "__main__":
    main()
