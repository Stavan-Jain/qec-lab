from lp_lab.lp import Group, lp_square
from lp_lab.structure import two_block_perms, toric_layouts
# gross: expect toric layout with torus (2m, 2ℓ)-ish = Z_12 x Z_24 (mu=m=6, lam=ell=12)
G = Group((12, 6)); A = {(0, 0): G.parse("x^3 + y + y^2")}; B = G.parse("y^3 + x + x^2")
code = lp_square(G, A, 1, G.conj(B)); P = two_block_perms(G, A, 1, G.conj(B), code)
lay = toric_layouts(P); print("gross toric layouts:", len(lay), lay[:2])
# a 2x2 designed example over Z56: b = 1+x^2+x^6, p = 1 + x, q with pq ≡ 1 mod (b*) -- just test the machinery on any code
G = Group((56,)); one = G.mono(0)
A = {(0,0): G.parse("1 + x^8"), (0,1): one, (1,0): one, (1,1): G.parse("x^3 + x^21")}
b = G.parse("1 + x^2 + x^6")
code = lp_square(G, A, 2, b); P = two_block_perms(G, A, 2, b, code)
lay = toric_layouts(P); print("2x2 Z56 example toric layouts:", len(lay), lay[:3])
