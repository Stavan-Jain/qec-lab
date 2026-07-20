"""Is the doubly-new sector of the (30,6) cover exactly im(1+sigma) on H1?

Checks, for each deck sigma in {x^15, y^3, x^15*y^3}: the class-space image
of (1+sigma) (dim 4 = k~-k, verified earlier) vs the doubly-new sector
ker Nx  cap  ker Ny (dim 4). Equality identifies the deck-born qubits with
the A13 D-module summand across all decks at once.
"""
import sys, importlib.util
from pathlib import Path

LAB = Path.home() / "Code/qec-lab/experiments/bb_lab"
sys.path.insert(0, str(LAB / "src"))
sys.path.insert(0, str(LAB / "scripts"))

import numpy as np
from bb_lab.group import ZmZn
from bb_lab.poly import Poly
from bb_lab.checks import bb_check_matrices
from bb_lab.linalg import rank_f2, nullspace_f2
from bb_lab.sat_distance import find_logical_z

spec = importlib.util.spec_from_file_location(
    "a15_coset_distance", LAB / "scripts" / "a15_coset_distance.py")
acd = importlib.util.module_from_spec(spec)
sys.modules["a15_coset_distance"] = acd
spec.loader.exec_module(acd)

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


Nx = (find_logical_z(chBX) @ (proj_matrix(GC, GBX, lambda g: (g[0] % 15, g[1])) @ SC.T % 2)) % 2
Ny = (find_logical_z(chBY) @ (proj_matrix(GC, GBY, lambda g: (g[0], g[1] % 3)) @ SC.T % 2)) % 2

# doubly-new sector in class coordinates: ker(Nx) cap ker(Ny)
DN = nullspace_f2(np.vstack([Nx, Ny]))   # rows: class-coord vectors, dim 4
print(f"doubly-new sector dim = {DN.shape[0]}")

for name, t in [("sigma_x  (x^15)", (15, 0)), ("sigma_y  (y^3)", (0, 3)),
                ("sigma_xy (x^15 y^3)", (15, 3))]:
    perm = acd.translation_perm(GC, t, nC)
    moved = (SC[:, perm] + SC) % 2                    # (1+sigma) on class reps
    # class coordinates of (1+sigma) images: signature -> coords via Pinv
    sigs = (L_ZC @ moved.T) % 2                       # k x k
    coords = (Pinv @ sigs) % 2                        # columns = images in class coords
    img = coords.T                                    # rows
    img = img[img.any(axis=1)]
    dim_img = rank_f2(img)
    # containment tests
    in_DN = rank_f2(np.vstack([DN, img])) == rank_f2(DN)
    contains_DN = rank_f2(np.vstack([img, DN])) == dim_img
    print(f"{name}: dim im(1+sigma)_* = {dim_img};  im subset of doubly-new: {in_DN};"
          f"  equal: {in_DN and dim_img == DN.shape[0]}")
