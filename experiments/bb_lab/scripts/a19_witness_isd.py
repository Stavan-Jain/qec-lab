"""Phase 0.2: recover a weight-24 old-sector witness for the (30,6) cover and
test the diagonal-lift hypothesis; save per-class ISD minima + vectors for the
beat-24 hunt (step8).

Outputs: witness_tau_report (stdout), isd_class_minima.npz (sig -> weight,vector).
"""
import sys, time, importlib.util, json
from pathlib import Path

LAB = Path.home() / "Code/qec-lab/experiments/bb_lab"
SCRATCH = Path(__file__).resolve().parent.parent / "data" / "a19"
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


def build():
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
    L_ZC = find_logical_z(chC)
    SC = acd.x_class_reps(chC)
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
    return acd, GC, chC, chBX, chBY, L_ZC, SC, Pinv, Px, Py, Nx, Ny


def main():
    acd, GC, chC, chBX, chBY, L_ZC, SC, Pinv, Px, Py, Nx, Ny = build()
    k = SC.shape[0]

    def sig_coords(sig):
        return (Pinv @ np.fromiter(((sig >> i) & 1 for i in range(k)),
                                   dtype=np.uint8, count=k)) % 2

    def sector_of_sig(sig):
        c = sig_coords(sig)
        inx = not np.any((Nx @ c) % 2)
        iny = not np.any((Ny @ c) % 2)
        return ("new-xy" if (inx and iny) else "new-x" if inx
                else "new-y" if iny else "old")

    # ---- ISD with vector saving (300 s) ----
    K = nullspace_f2(chC.H_Z)
    rng = np.random.default_rng(7)
    best = {}          # sig -> (weight, vector)
    t0 = time.time()
    iters = 0
    while time.time() - t0 < 300:
        iters += 1
        perm = rng.permutation(K.shape[1])
        R, _ = rref_f2(K[:, perm].copy())
        R = R[R.any(axis=1)]
        sigs = (L_ZC[:, perm] @ R.T) % 2
        wts = R.sum(axis=1)
        for j in range(R.shape[0]):
            col = sigs[:, j]
            sig = 0
            for i in range(k):
                if col[i]:
                    sig |= 1 << i
            if sig == 0:
                continue
            w = int(wts[j])
            if w < best.get(sig, (10**9, None))[0]:
                v = np.zeros(K.shape[1], dtype=np.uint8)
                v[perm] = R[j]
                best[sig] = (w, v)
    print(f"ISD: {iters} iterations, {len(best)} classes", flush=True)

    np.savez_compressed(
        SCRATCH / "isd_class_minima.npz",
        sigs=np.array(sorted(best.keys()), dtype=np.int64),
        weights=np.array([best[s][0] for s in sorted(best.keys())], dtype=np.int64),
        vectors=np.array([best[s][1] for s in sorted(best.keys())], dtype=np.uint8))

    sector_min = {}
    for sig, (w, v) in best.items():
        s = sector_of_sig(sig)
        if w < sector_min.get(s, (10**9, None, None))[0]:
            sector_min[s] = (w, sig, v)
    for s, (w, sig, _) in sorted(sector_min.items()):
        print(f"  sector {s:>7}: min wt {w} (class {sig:#05x})", flush=True)

    # ---- diagonal-lift test on the global minimum witness ----
    w_glob, sig_glob, v = min(((w, sig, v) for sig, (w, v) in best.items()))
    print(f"\nglobal ISD minimum: wt {w_glob}, class {sig_glob:#05x} "
          f"[{sector_of_sig(sig_glob)}]", flush=True)

    def tau(P, u):
        return (P.T @ u) % 2

    for name, P, chB in [("p_x -> BX", Px, chBX), ("p_y -> BY", Py, chBY)]:
        pv = (P @ v) % 2
        wpv = int(pv.sum())
        in_ker = not np.any((chB.H_X @ pv) % 2)
        RZ, _ = rref_f2(chB.H_Z.copy())
        RZ = RZ[RZ.any(axis=1)]
        nontriv = rank_f2(np.vstack([RZ, pv[None, :]])) > RZ.shape[0] if in_ker else False
        print(f"  {name}: |p(v)| = {wpv}; cycle: {in_ker}; nontrivial logical: "
              f"{nontriv}; |v| == 2|p(v)|: {w_glob == 2 * wpv}", flush=True)
        if in_ker and nontriv and w_glob == 2 * wpv:
            lift = tau(P, pv)
            same_supp = bool(np.array_equal(lift, v))
            # v equals full preimage lift up to a cover stabilizer?
            RZc, _ = rref_f2(chC.H_Z.copy())
            RZc = RZc[RZc.any(axis=1)]
            diff = (lift + v) % 2
            triv_diff = rank_f2(np.vstack([RZc, diff[None, :]])) == RZc.shape[0]
            print(f"    exact preimage lift: {same_supp}; "
                  f"v = lift + stabilizer: {triv_diff}", flush=True)


if __name__ == "__main__":
    main()
