"""Are BX=(15,6) and BY=(30,3) the same code up to Aut+monomial equivalence?

Both are quotients of the (30,6) cover by the two axis decks; groups are
abstractly isomorphic (Z2 x Z3^2 x Z5). Re-coordinate BX into Z30xZ3 via CRT
iso and compare canonical pairs under the full Aut(G) x translation x swap
orbit (bb_lab.canonical).
"""
import sys
from pathlib import Path

LAB = Path.home() / "Code/qec-lab/experiments/bb_lab"
sys.path.insert(0, str(LAB / "src"))

from bb_lab.group import ZmZn
from bb_lab.canonical import canonical_pair, build_perm_table

G = ZmZn(30, 3)

# BX on Z15xZ6: A = x^9+y+y^2, B = y^3+x^10+x^11
# iso phi: Z15 x Z6 -> Z30 x Z3, phi(u,v) = (s, t): s = CRT(u mod 15, v mod 2), t = v mod 3
def phi(u, v):
    for s in range(30):
        if s % 15 == u % 15 and s % 2 == v % 2:
            return (s, v % 3)
    raise AssertionError

BX_A = {phi(9, 0), phi(0, 1), phi(0, 2)}
BX_B = {phi(0, 3), phi(10, 0), phi(11, 0)}

# BY on Z30xZ3 (native): A = x^9+y+y^2, B = 1+x^25+x^26
BY_A = {(9, 0), (0, 1), (0, 2)}
BY_B = {(0, 0), (25, 0), (26, 0)}

perms = build_perm_table(G)
cx = canonical_pair(BX_A, BX_B, G, perms=perms)
cy = canonical_pair(BY_A, BY_B, G, perms=perms)
print("BX canonical:", sorted(cx.A_support), sorted(cx.B_support))
print("BY canonical:", sorted(cy.A_support), sorted(cy.B_support))
print("EQUIVALENT" if cx.key == cy.key else "NOT equivalent (distinct codes)")
