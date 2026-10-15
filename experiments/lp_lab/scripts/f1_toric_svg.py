"""Draw the toric layout of an a=2 code as SVG: grid edges (4 per vertex) plus the two
long-range edges of a few highlighted vertices, in the style of Fig. 1B of arXiv:2308.07915."""
import argparse
from lp_lab.lp import Group, lp_square
from lp_lab.structure import two_block_perms, toric_layouts
from lp_lab.distance import code_k

ap = argparse.ArgumentParser()
ap.add_argument("--orders", required=True); ap.add_argument("--p", required=True); ap.add_argument("--q", required=True); ap.add_argument("--b", required=True)
ap.add_argument("--out", required=True); ap.add_argument("--layout", type=int, default=0); ap.add_argument("--highlight", type=int, default=3)
args = ap.parse_args()
orders = tuple(int(x) for x in args.orders.split(",")); G = Group(orders); one = G.mono(*([0] * len(orders)))
p = frozenset(int(x) for x in args.p.split(",")); q = frozenset(int(x) for x in args.q.split(",")); b = frozenset(int(x) for x in args.b.split(","))
A = {(0, 0): p, (0, 1): one, (1, 0): one, (1, 1): q}
code = lp_square(G, A, 2, b); P = two_block_perms(G, A, 2, b, code); m = P.m
lay = toric_layouts(P)[args.layout]
i, j, g, h = lay["i"], lay["j"], lay["g"], lay["h"]; mu, lam = lay["mu"], lay["lam"]
inv = P.inv
sigma = P.A[j][inv(P.A[i])]; tau = P.B[g][inv(P.B[h])]
pos = {}; cur = 0
for a in range(mu):
    x = cur
    for bb in range(lam):
        pos[x] = (a, bb); x = tau[x]
    cur = sigma[cur]
Ai_inv = inv(P.A[i])
coord = {}   # vertex id -> (col, row); ids: X c, L m+l, R 2m+r, Z 3m+z
for ell, (a, bb) in pos.items():
    coord[m + ell] = (2 * a, 2 * bb)
    c = Ai_inv[ell]; coord[c] = (2 * a + 1, 2 * bb)
    coord[2 * m + P.B[g][c]] = (2 * a + 1, 2 * bb + 1)
    coord[3 * m + P.B[g][ell]] = (2 * a, 2 * bb + 1)
M1, M2 = 2 * mu, 2 * lam
# all Tanner edges
edges = []
for c in range(m):
    for pp in range(3):
        edges.append((c, m + P.A[pp][c])); edges.append((c, 2 * m + P.B[pp][c]))
for z in range(m):
    for pp in range(3):
        edges.append((3 * m + z, 2 * m + inv(P.A[pp])[z])); edges.append((3 * m + z, m + inv(P.B[pp])[z]))
def is_grid(u, v):
    (a1, b1), (a2, b2) = coord[u], coord[v]
    da = (a1 - a2) % M1; db = (b1 - b2) % M2
    return (da in (1, M1 - 1) and db == 0) or (db in (1, M2 - 1) and da == 0)
S = 22; W = M1 * S + 60; H = M2 * S + 80
def xy(v):
    a, bb = coord[v]; return 30 + a * S + S / 2, 40 + bb * S + S / 2
svg = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" viewBox="0 0 {W} {H}" font-family="Helvetica, Arial, sans-serif">',
       f'<rect width="{W}" height="{H}" fill="white"/>',
       f'<text x="30" y="24" font-size="14" fill="#222">Toric layout of the [[{code.n},{code_k(code)},·]] code on Z_{M1} × Z_{M2}: 4 grid couplers per vertex; long-range edges shown for {2*args.highlight} vertices</text>']
# grid edges
for u, v in edges:
    if is_grid(u, v):
        (x1, y1), (x2, y2) = xy(u), xy(v)
        if abs(x1 - x2) > S * 1.5 or abs(y1 - y2) > S * 1.5:  # wrap-around: draw short stubs
            continue
        svg.append(f'<line x1="{x1:.1f}" y1="{y1:.1f}" x2="{x2:.1f}" y2="{y2:.1f}" stroke="#999" stroke-width="1.2"/>')
# long-range edges for a few L qubits and a few X checks
colors = ["#d62728", "#1f77b4", "#2ca02c", "#9467bd", "#ff7f0e", "#8c564b"]
hl = [m + ell for ell in list(pos.keys())[:: max(1, m // args.highlight)][: args.highlight]]
hl += [Ai_inv[ell] for ell in list(pos.keys())[1:: max(1, m // args.highlight)][: args.highlight]]
ci = 0
for v in hl:
    col = colors[ci % len(colors)]; ci += 1
    for u, w in edges:
        if v in (u, w) and not is_grid(u, w):
            (x1, y1), (x2, y2) = xy(u), xy(w)
            mx, my = (x1 + x2) / 2, (y1 + y2) / 2 - 0.25 * abs(x2 - x1) - 8
            svg.append(f'<path d="M {x1:.1f} {y1:.1f} Q {mx:.1f} {my:.1f} {x2:.1f} {y2:.1f}" fill="none" stroke="{col}" stroke-width="1.8" stroke-dasharray="4 3"/>')
# vertices
style = {"X": ("#e8a33d", "rect"), "Z": ("#5b8def", "rect"), "L": ("#4c9a2a", "circle"), "R": ("#c0392b", "circle")}
for v, (a, bb) in coord.items():
    t = "X" if v < m else "L" if v < 2 * m else "R" if v < 3 * m else "Z"
    x, y = xy(v); col, shape = style[t]
    if shape == "rect":
        svg.append(f'<rect x="{x-5:.1f}" y="{y-5:.1f}" width="10" height="10" fill="{col}" stroke="#333" stroke-width="0.6"/>')
    else:
        svg.append(f'<circle cx="{x:.1f}" cy="{y:.1f}" r="5.5" fill="{col}" stroke="#333" stroke-width="0.6"/>')
ly = H - 22
svg.append(f'<text x="30" y="{ly}" font-size="12" fill="#333">'
           f'<tspan fill="#4c9a2a">●</tspan> L data  <tspan fill="#c0392b">●</tspan> R data  '
           f'<tspan fill="#e8a33d">■</tspan> X check  <tspan fill="#5b8def">■</tspan> Z check   (torus {M1}×{M2}; wrap-around grid edges omitted)</text>')
svg.append("</svg>")
open(args.out, "w").write("\n".join(svg)); print("wrote", args.out, "torus", (M1, M2), "vertices", len(coord), "edges", len(edges))
