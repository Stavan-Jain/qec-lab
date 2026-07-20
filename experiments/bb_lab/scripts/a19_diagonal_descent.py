"""Diagonal-deck descent of the (30,6) Bravyi cover: quotient by <x^15 y^3>.

Iso phi: Z30xZ6 / <(15,3)> ~= Z30xZ3 via phi(u,v) = (u+15v mod 30, v mod 3).
Sweep the cover's monomial orbit (units a mod 30, b mod 6) composed with phi.
k(descent)=12 <=> (R) holds on the diagonal deck for that presentation.
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

A_SUPP = [(9, 0), (0, 1), (0, 2)]
B_SUPP = [(0, 3), (25, 0), (26, 0)]


def diag_descend_k(supp_a, supp_b):
    G = ZmZn(30, 3)
    def phi(e, f):
        return ((e + 15 * f) % 30, f % 3)
    A = Poly.from_support([phi(e, f) for (e, f) in supp_a], G)
    B = Poly.from_support([phi(e, f) for (e, f) in supp_b], G)
    if A.weight() < 3 or B.weight() < 3:
        return None, None, None
    ch = bb_check_matrices(A, B)
    return code_params(ch).k, A.canonical_string(), B.canonical_string()


units30 = [a for a in range(1, 30) if gcd(a, 30) == 1]
units6 = [b for b in range(1, 6) if gcd(b, 6) == 1]

hits = []
for a in units30:
    for b in units6:
        sa = [((e * a) % 30, (f * b) % 6) for (e, f) in A_SUPP]
        sb = [((e * a) % 30, (f * b) % 6) for (e, f) in B_SUPP]
        k, As, Bs = diag_descend_k(sa, sb)
        mark = "  <-- k-preserving!" if k == 12 else ""
        print(f"a={a:>2} b={b} | k(diag descent (30,3)) = {k}{mark}")
        if k == 12:
            hits.append((a, b, As, Bs))

print()
if hits:
    print(f"{len(hits)} k-preserving diagonal descents; representatives:")
    for a, b, As, Bs in hits[:4]:
        print(f"  a={a} b={b}:  A = {As}   B = {Bs}")
else:
    print("no k-preserving diagonal descent")
