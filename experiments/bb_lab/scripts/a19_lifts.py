"""Corrected lift analysis (X-side conventions).

(a) ISD at BY and BX for weight-12 X-logicals (u in ker H_Z(base), nontrivial
    vs rowspan H_X(base)).
(b) Transfer lifts tau(u) = P.T @ u at the cover: verify cycle, weight 2|u|,
    nontrivial, and sector (predict: killed by that deck's pushforward).
(c) Double lift of a GB weight-8 logical -> predict new-xy at 32.
(d) Old-sector minimum vector (class 0x6e2 from saved npz): projection weights
    and side-correct verdicts.
"""
import sys, time, importlib.util
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


def isd_min(ch, budget_s, rng):
    """(best_weight, vector) over nontrivial X-logicals via ISD on ker H_Z."""
    K = nullspace_f2(ch.H_Z)
    RX, _ = rref_f2(ch.H_X.copy())
    RX = RX[RX.any(axis=1)]
    rX = RX.shape[0]
    best_w, best_v = 10**9, None
    t0 = time.time()
    while time.time() - t0 < budget_s:
        perm = rng.permutation(K.shape[1])
        R, _ = rref_f2(K[:, perm].copy())
        R = R[R.any(axis=1)]
        wts = R.sum(axis=1)
        for j in np.argsort(wts):
            w = int(wts[j])
            if w >= best_w:
                break
            v = np.zeros(K.shape[1], dtype=np.uint8)
            v[perm] = R[j]
            if rank_f2(np.vstack([RX, v[None, :]])) > rX:  # nontrivial
                best_w, best_v = w, v
    return best_w, best_v


def main():
    acd = load_acd()
    rng = np.random.default_rng(11)
    GC = ZmZn(30, 6)
    chC = bb_check_matrices(Poly.from_string("x^9 + y + y^2", GC),
                            Poly.from_string("y^3 + x^25 + x^26", GC))
    GBX = ZmZn(15, 6)
    chBX = bb_check_matrices(Poly.from_string("x^9 + y + y^2", GBX),
                             Poly.from_string("y^3 + x^10 + x^11", GBX))
    GBY = ZmZn(30, 3)
    chBY = bb_check_matrices(Poly.from_string("x^9 + y + y^2", GBY),
                             Poly.from_string("1 + x^25 + x^26", GBY))
    GGB = ZmZn(15, 3)
    chGB = bb_check_matrices(Poly.from_string("x^9 + y + y^2", GGB),
                             Poly.from_string("1 + x^10 + x^11", GGB))

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

    RXc, _ = rref_f2(chC.H_X.copy())
    RXc = RXc[RXc.any(axis=1)]
    rXc = RXc.shape[0]

    def cover_verdict(v, label):
        cyc = not np.any((chC.H_Z @ v) % 2)
        nontriv = cyc and rank_f2(np.vstack([RXc, v[None, :]])) > rXc
        sig = 0
        if nontriv:
            s = (L_ZC @ v) % 2
            for i in range(k):
                if s[i]:
                    sig |= 1 << i
        c = (Pinv @ np.fromiter(((sig >> i) & 1 for i in range(k)),
                                dtype=np.uint8, count=k)) % 2
        inx = not np.any((Nx @ c) % 2)
        iny = not np.any((Ny @ c) % 2)
        sector = ("new-xy" if (inx and iny) else "new-x" if inx
                  else "new-y" if iny else "old")
        print(f"  {label}: wt {int(v.sum())}; cycle {cyc}; nontrivial {nontriv}; "
              f"class {sig:#05x} [{sector}]", flush=True)
        return nontriv

    # (a)+(b): base minima and their lifts
    for name, chB, P in [("BY", chBY, Py), ("BX", chBX, Px)]:
        w, u = isd_min(chB, 60, rng)
        print(f"{name}: ISD min nontrivial X-logical wt = {w}", flush=True)
        lift = (P.T @ u) % 2
        cover_verdict(lift, f"tau_{name}(u) lift")

    # (c): double lift from GB
    wg, ug = isd_min(chGB, 30, rng)
    print(f"GB: ISD min nontrivial X-logical wt = {wg}", flush=True)
    PxG = proj_matrix(GBX, GGB, lambda g: (g[0] % 15, g[1] % 3))  # BX -> GB? no:
    # build GB -> BX lift then BX -> C lift: qubits of BX project to GB by y mod 3
    PyG = proj_matrix(GBX, GGB, lambda g: (g[0], g[1] % 3))
    lift1 = (PyG.T @ ug) % 2          # GB logical lifted into BX chain
    lift2 = (Px.T @ lift1) % 2        # then into the cover
    cover_verdict(lift2, "tau_x(tau_y(u_GB)) double lift")

    # (d): old-sector minimum from saved ISD data
    data = np.load(SCRATCH / "isd_class_minima.npz")
    sigs, wts, vecs = data["sigs"], data["weights"], data["vectors"]
    idx = {int(s): i for i, s in enumerate(sigs)}
    i = idx[0x6e2]
    v = vecs[i]
    print(f"old-sector min class 0x6e2: wt {int(wts[i])}", flush=True)
    for nm, P, chB in [("p_x", Px, chBX), ("p_y", Py, chBY)]:
        pv = (P @ v) % 2
        cyc = not np.any((chB.H_Z @ pv) % 2)
        RB, _ = rref_f2(chB.H_X.copy())
        RB = RB[RB.any(axis=1)]
        nt = cyc and rank_f2(np.vstack([RB, pv[None, :]])) > RB.shape[0]
        print(f"  {nm}(v): wt {int(pv.sum())}; cycle {cyc}; nontrivial {nt}",
              flush=True)


if __name__ == "__main__":
    main()
