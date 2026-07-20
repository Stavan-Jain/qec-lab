"""Orbit-sweep the cover's monomial-equivalence presentations for a k-preserving descent.

Cover: (30,6), A=x^9+y+y^2, B=y^3+x^25+x^26, k=12.
For each unit a in (Z/30)*, b in (Z/6)*: transform exponents (x->x^a, y->y^b)
(an automorphism of the cover code), then descend along x (exponents mod 15)
and along y (exponents mod 3). A k=12 descent <=> (R) holds for that presentation
(A12: deck-trivial iff k-preserving).
"""
import sys
from pathlib import Path

LAB = Path.home() / "Code/qec-lab/experiments/bb_lab"
sys.path.insert(0, str(LAB / "src"))

from math import gcd
from bb_lab.group import ZmZn
from bb_lab.poly import Poly
from bb_lab.checks import bb_check_matrices
from bb_lab.codeparams import code_params

ELL, M = 30, 6
A_SUPP = [(9, 0), (0, 1), (0, 2)]
B_SUPP = [(0, 3), (25, 0), (26, 0)]


def descend_k(supp_a, supp_b, ell, m):
    G = ZmZn(ell, m)
    A = Poly.from_support([(e % ell, f % m) for (e, f) in supp_a], G)
    B = Poly.from_support([(e % ell, f % m) for (e, f) in supp_b], G)
    if A.weight() < 3 or B.weight() < 3:
        return None  # support collision under reduction: not a weight-3 descent
    ch = bb_check_matrices(A, B)
    return code_params(ch).k


units30 = [a for a in range(1, 30) if gcd(a, 30) == 1]
units6 = [b for b in range(1, 6) if gcd(b, 6) == 1]

print(f"{'a':>3} {'b':>2} | k(x-descent (15,6)) | k(y-descent (30,3))")
hits = []
for a in units30:
    for b in units6:
        sa = [((e * a) % ELL, (f * b) % M) for (e, f) in A_SUPP]
        sb = [((e * a) % ELL, (f * b) % M) for (e, f) in B_SUPP]
        kx = descend_k(sa, sb, 15, 6)
        ky = descend_k(sa, sb, 30, 3)
        tag = ""
        if kx == 12 or ky == 12:
            tag = "  <-- k-preserving descent!"
            hits.append((a, b, kx, ky))
        print(f"{a:>3} {b:>2} | {str(kx):>19} | {str(ky):>19}{tag}")

print()
print(f"k-preserving descents found: {len(hits)}" + (f" -> {hits}" if hits else ""))
