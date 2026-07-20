"""Cover-class weight census for Bravyi (30,6) [[360,12,<=24]] — MC-first, budgeted SAT.

Phase 1 (fast): ISD-style Monte Carlo over ker H_Z, binning best witness weight
per logical class signature (12 bits) and per deck-sector (old / new-x / new-y
/ new-xy, via pushforward kernels to the two k=8 bases).
Phase 2 (budgeted): per-coset CryptoMiniSat queries trying to beat the MC minima
on selected classes; hard wall-clock deadline, pool terminated at budget.
"""
import sys, time, importlib.util, json
import multiprocessing as mp
from collections import Counter
from pathlib import Path

LAB = Path.home() / "Code/qec-lab/experiments/bb_lab"
sys.path.insert(0, str(LAB / "src"))
sys.path.insert(0, str(LAB / "scripts"))

import numpy as np
from bb_lab.group import ZmZn
from bb_lab.poly import Poly
from bb_lab.checks import bb_check_matrices
from bb_lab.linalg import rank_f2, nullspace_f2, rref_f2
from bb_lab.sat_distance import find_logical_z


def load_acd():
    spec = importlib.util.spec_from_file_location(
        "a15_coset_distance", LAB / "scripts" / "a15_coset_distance.py")
    acd = importlib.util.module_from_spec(spec)
    sys.modules["a15_coset_distance"] = acd
    spec.loader.exec_module(acd)
    return acd


def main():
    acd = load_acd()
    GC = ZmZn(30, 6)
    chC = bb_check_matrices(Poly.from_string("x^9 + y + y^2", GC),
                            Poly.from_string("y^3 + x^25 + x^26", GC))
    GBX = ZmZn(15, 6)
    chBX = bb_check_matrices(Poly.from_string("x^9 + y + y^2", GBX),
                             Poly.from_string("y^3 + x^10 + x^11", GBX))
    GBY = ZmZn(30, 3)
    chBY = bb_check_matrices(Poly.from_string("x^9 + y + y^2", GBY),
                             Poly.from_string("1 + x^25 + x^26", GBY))
    nC = GC.cardinality

    L_ZC = find_logical_z(chC)
    SC = acd.x_class_reps(chC)
    k = SC.shape[0]
    Pinv = acd.inv_f2((L_ZC @ SC.T) % 2)

    def proj_matrix(Gcov, Gbase, red):
        ncov, nbase = Gcov.cardinality, Gbase.cardinality
        P = np.zeros((2 * nbase, 2 * ncov), dtype=np.uint8)
        for i, g in enumerate(Gcov):
            j = Gbase.index(red(g))
            P[j, i] = 1
            P[nbase + j, ncov + i] = 1
        return P

    Px = proj_matrix(GC, GBX, lambda g: (g[0] % 15, g[1]))
    Py = proj_matrix(GC, GBY, lambda g: (g[0], g[1] % 3))
    Nx = (find_logical_z(chBX) @ (Px @ SC.T % 2)) % 2
    Ny = (find_logical_z(chBY) @ (Py @ SC.T % 2)) % 2

    def sector_of_sig(sig):
        c = (Pinv @ np.fromiter(((sig >> i) & 1 for i in range(k)),
                                dtype=np.uint8, count=k)) % 2
        inx = not np.any((Nx @ c) % 2)
        iny = not np.any((Ny @ c) % 2)
        return ("new-xy" if (inx and iny) else "new-x" if inx
                else "new-y" if iny else "old")

    # ---------- phase 1: ISD Monte Carlo census ----------
    K = nullspace_f2(chC.H_Z)
    print(f"ker H_Z dim = {K.shape[0]} over n2 = {K.shape[1]}", flush=True)
    rng = np.random.default_rng(20260720)
    best_by_sig = {}
    t0 = time.time()
    iters = 0
    while time.time() - t0 < 150:
        iters += 1
        perm = rng.permutation(K.shape[1])
        R, _ = rref_f2(K[:, perm].copy())
        R = R[R.any(axis=1)]
        sigs = (L_ZC[:, perm] @ R.T) % 2
        wts = R.sum(axis=1)
        for j in range(R.shape[0]):
            sig = 0
            col = sigs[:, j]
            for i in range(k):
                if col[i]:
                    sig |= 1 << i
            if sig == 0:
                continue
            w = int(wts[j])
            if w < best_by_sig.get(sig, 10**9):
                v = np.zeros(K.shape[1], dtype=np.uint8)
                v[perm] = R[j]
                best_by_sig[sig] = w
    print(f"MC census: {iters} ISD iterations, {len(best_by_sig)} classes seen",
          flush=True)

    sector_best = {}
    for sig, w in best_by_sig.items():
        s = sector_of_sig(sig)
        if w < sector_best.get(s, (10**9, 0))[0]:
            sector_best[s] = (w, sig)
    print("per-sector MC minima:", flush=True)
    for s, (w, sig) in sorted(sector_best.items()):
        print(f"  {s:>7}: wt {w}  (class {sig:#05x})", flush=True)
    gmin = min(w for w, _ in sector_best.values())
    print(f"MC upper bound on d(cover): {gmin}", flush=True)

    # ---------- phase 2: budgeted SAT tightening ----------
    # targets: per-sector argmin classes at (mc_best - 2), plus new-xy orbit reps
    mats, _ = acd.signature_action(GC, chC, L_ZC, SC)
    reps = acd.orbit_reps(mats, list(GC.orders), k)
    newxy_reps = [sig for sig, _ in reps if sector_of_sig(sig) == "new-xy"]
    targets = []
    for s, (w, sig) in sector_best.items():
        targets.append((sig, max(2, w - 2), f"tighten-{s}"))
    for sig in newxy_reps:
        w = best_by_sig.get(sig)
        bound = (w - 2) if w else 20
        targets.append((sig, max(2, bound), "new-xy-rep"))
    seen = set()
    targets = [t for t in targets if not (t[0] in seen or seen.add(t[0]))]

    print(f"SAT phase: {len(targets)} queries, budget 300 s", flush=True)
    ctx = mp.get_context("spawn")
    deadline = time.time() + 300
    results = []
    with ctx.Pool(8) as pool:
        async_res = [(sig, W, tag,
                      pool.apply_async(acd._query_worker, ((chC, L_ZC, sig, W),)))
                     for sig, W, tag in targets]
        for sig, W, tag, ar in async_res:
            left = deadline - time.time()
            if left <= 0:
                results.append((sig, W, tag, "TIMEOUT"))
                continue
            try:
                _, v, secs = ar.get(timeout=left)
                out = f"SAT wt {int(v.sum())}" if v is not None else f"UNSAT@{W}"
                results.append((sig, W, tag, out))
                print(f"  class {sig:#05x} [{tag}] w<={W}: {out} ({secs:.0f}s)",
                      flush=True)
            except mp.TimeoutError:
                results.append((sig, W, tag, "TIMEOUT"))
                print(f"  class {sig:#05x} [{tag}] w<={W}: TIMEOUT", flush=True)
        pool.terminate()

    print(json.dumps({
        "mc_iters": iters,
        "sector_mc_min": {s: w for s, (w, _) in sector_best.items()},
        "d_cover_upper_bound": gmin,
        "sat_phase": [[f"{sig:#05x}", W, tag, out]
                      for sig, W, tag, out in results],
    }), flush=True)


if __name__ == "__main__":
    main()
