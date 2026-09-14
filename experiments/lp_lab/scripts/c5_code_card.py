"""Export a self-contained 'code card' JSON for an a=2 code: definition, parameters, permutation
parts for the depth-7 circuit, wheel layers, toric layout coordinates, and available certificates."""
import argparse, json, glob
from pathlib import Path
import numpy as np
from lp_lab.lp import Group, lp_square
from lp_lab.distance import code_k
from lp_lab.structure import two_block_perms, thickness2_certificate, toric_layouts, tanner_automorphisms

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--orders", required=True); ap.add_argument("--p", required=True); ap.add_argument("--q", required=True); ap.add_argument("--b", required=True)
    ap.add_argument("--name", required=True)
    args = ap.parse_args()
    orders = tuple(int(x) for x in args.orders.split(",")); G = Group(orders); one = G.mono(*([0] * len(orders)))
    p = frozenset(int(x) for x in args.p.split(",")); q = frozenset(int(x) for x in args.q.split(",")); b = frozenset(int(x) for x in args.b.split(","))
    A = {(0, 0): p, (0, 1): one, (1, 0): one, (1, 1): q}
    code = lp_square(G, A, 2, b); code.check_css(); n = code.n; k = code_k(code)
    P = two_block_perms(G, A, 2, b, code)
    aut = tanner_automorphisms(code.H_X, code.H_Z); aut.pop("orbit_labels")
    cert = thickness2_certificate(P); good = [c for c in cert if c["planar_A"] and c["planar_B"]]
    lay = toric_layouts(P)
    card = {
        "name": args.name,
        "family": "monomial lifted product LP(A, b⊗I), a = 2 (experiments/lp_lab, LP1)",
        "group": {"orders": orders, "variables": "xyz"[:len(orders)], "element_index": "row-major over the cyclic factors"},
        "A": {"(0,0)": G.poly_str(p), "(0,1)": "1", "(1,0)": "1", "(1,1)": G.poly_str(q)},
        "b": G.poly_str(b),
        "check_matrices": "H_X = [A | b* I_2], H_Z = [b I_2 | A†] with M(f)[g,h] = [g−h ∈ supp f]; qubit order L0,L1,R0,R1 blocks of |G|",
        "n": n, "k": k, "weights": code.weights(),
        "aut": aut, "not_2BGA_certificate": aut["orbits_X"] >= 2,
        "thickness2": {"planar_layer_choices_of_9": len(good), "layer_example": good[0] if good else None,
                       "layers": "GA = {A2, A3, B3} edges, GB = {A1, B1, B2} (indices per layer_example)"},
        "toric_layouts": lay,
        "circuit": {"schedule": "arXiv:2308.07915 Table 5 with the permutation parts below",
                    "A_parts (X-check c -> L qubit)": [P.A[i].tolist() for i in range(3)],
                    "B_parts (X-check c -> R qubit)": [P.B[i].tolist() for i in range(3)]},
        "certificates": {},
    }
    for f in glob.glob(f"data/certs/*{args.name}*.json"):
        try:
            card["certificates"][Path(f).name] = json.loads(Path(f).read_text())
        except Exception as e:
            card["certificates"][Path(f).name] = f"unreadable: {e}"
    Path("data/codes").mkdir(exist_ok=True)
    out = Path(f"data/codes/{args.name}.json"); out.write_text(json.dumps(card, indent=1, default=str))
    print("wrote", out, f"[[{n},{k},?]]", "toric layouts:", len(lay), "planar choices:", len(good), "X-orbits:", aut["orbits_X"])

if __name__ == "__main__":
    main()
