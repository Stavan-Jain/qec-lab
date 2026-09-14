import numpy as np, time
from lp_lab.lp import Group, lp_square
from lp_lab.structure import tanner_automorphisms, two_block_perms, thickness2_certificate
G = Group((12, 6))
A = {(0, 0): G.parse("x^3 + y + y^2")}; B = G.parse("y^3 + x + x^2")
code = lp_square(G, A, 1, G.conj(B))
P = two_block_perms(G, A, 1, G.conj(B), code)
print("perm parts verified against H_X, H_Z")
cert = thickness2_certificate(P)
print("planar layer choices:", sum(c['planar_A'] and c['planar_B'] for c in cert), "/ 9;", cert[0])
# a 2x2 example (non-BB) to exercise the general path
G = Group((6, 6))
A = {(0,0): G.parse("1 + x*y^2"), (0,1): G.mono(0,0), (1,0): G.mono(0,0), (1,1): G.parse("x^2 + y^5")}
b = G.parse("1 + x + y^3")
code = lp_square(G, A, 2, b); code.check_css()
P = two_block_perms(G, A, 2, b, code)
cert = thickness2_certificate(P)
print("2x2 example: planar layer choices:", sum(c['planar_A'] and c['planar_B'] for c in cert), "/ 9; wheelsA", cert[0]['wheelsA'], "wheelsB", cert[0]['wheelsB'])
aut = tanner_automorphisms(code.H_X, code.H_Z); print("2x2 Aut:", {k: v for k, v in aut.items() if k != 'orbit_labels'})
