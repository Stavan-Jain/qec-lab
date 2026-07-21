"""A20 instance study: IBM's [[288,8,20]] (class Y, arXiv:2606.02418) as a tower.

The code: (ell,m) = (18,8), A = 1+xy^4+x^14y, B = 1+xy^2+x^2y^7 (weight-3,
weight-6 stabilizers), k = 8, d = 20 MILP-exact per IBM's Table II — found by
the A19 catalog deck-birth survey to be born at (18,2) and to propagate
through two y-doublings: the cheapest known d >= 20 certification target.

This script computes, kernel-level (F2 linear algebra only):
  * k at every rung of the y-tower (18,8) -> (18,4) -> (18,2), plus every
    order-2 deck quotient at each rung (x / y / diagonal, coset convolution);
  * the (R)/Bezout verdict for each y-rung: is 1 + y^{m/2} in the ideal
    (A,B) of F2[Z_ell x Z_m]?  If yes, saves an explicit Bezout witness
    (P, Q) with P*A + Q*B = 1 + y^{m/2} to data/a20/ — the input for the
    Lean `deckTrivial_of_bezout` certificate;
  * the reduced polynomials at each rung (for ladder JSON rows).

Usage (from experiments/bb_lab):
    uv run scripts/a20_ibm288_tower.py

Follow-ups (A20 note SS3): exact base distances via
    uv run scripts/a15_coset_distance.py '<row json printed below>' ...
"""
import json
from pathlib import Path

import numpy as np

ELL, M = 18, 8
A_SUPP = {(0, 0), (1, 4), (14, 1)}
B_SUPP = {(0, 0), (1, 2), (2, 7)}
OUT = Path(__file__).resolve().parent.parent / "data" / "a20"


def rank_f2(Mx):
    Mx = Mx.copy().astype(np.uint8)
    r = 0
    for c in range(Mx.shape[1]):
        piv = np.nonzero(Mx[r:, c])[0]
        if piv.size == 0:
            continue
        p = r + piv[0]
        Mx[[r, p]] = Mx[[p, r]]
        rows = np.nonzero(Mx[:, c])[0]
        Mx[rows[rows != r]] ^= Mx[r]
        r += 1
        if r == Mx.shape[0]:
            break
    return r


def solve_f2(Mx, b):
    """One solution of Mx z = b over F2, or None."""
    n_unk = Mx.shape[1]
    Aug = np.concatenate([Mx.astype(np.uint8), b.reshape(-1, 1)], axis=1)
    r, pivots = 0, []
    for c in range(n_unk):
        piv = np.nonzero(Aug[r:, c])[0]
        if piv.size == 0:
            continue
        p = r + piv[0]
        Aug[[r, p]] = Aug[[p, r]]
        rows = np.nonzero(Aug[:, c])[0]
        Aug[rows[rows != r]] ^= Aug[r]
        pivots.append(c)
        r += 1
        if r == Aug.shape[0]:
            break
    if any(Aug[i, n_unk] and not Aug[i, :n_unk].any() for i in range(Aug.shape[0])):
        return None
    z = np.zeros(n_unk, dtype=np.uint8)
    for i, c in enumerate(pivots):
        z[c] = Aug[i, n_unk]
    return z


def conv_matrix(ell, m, supp):
    elems = [(a, b) for a in range(ell) for b in range(m)]
    idx = {e: i for i, e in enumerate(elems)}
    Mx = np.zeros((len(elems),) * 2, dtype=np.uint8)
    for j, (gx, gy) in enumerate(elems):
        for (ax, ay) in supp:
            Mx[idx[((gx + ax) % ell, (gy + ay) % m)], j] ^= 1
    return Mx


def code_k(ell, m, sA, sB):
    MA, MB = conv_matrix(ell, m, sA), conv_matrix(ell, m, sB)
    HX = np.concatenate([MA, MB], axis=1)
    HZ = np.concatenate([MB.T, MA.T], axis=1)
    return 2 * ell * m - rank_f2(HX) - rank_f2(HZ)


def quotient_k(ell, m, sA, sB, t):
    tx, ty = t
    seen, reps, ridx = set(), [], {}
    for a in range(ell):
        for b in range(m):
            c = frozenset({(a, b), ((a + tx) % ell, (b + ty) % m)})
            if c not in seen:
                seen.add(c)
                reps.append(min(c))
    for i, rr in enumerate(reps):
        ridx[rr] = i
        ridx[((rr[0] + tx) % ell, (rr[1] + ty) % m)] = i

    def qconv(supp):
        Mx = np.zeros((len(reps),) * 2, dtype=np.uint8)
        for j, (gx, gy) in enumerate(reps):
            for (ax, ay) in supp:
                Mx[ridx[((gx + ax) % ell, (gy + ay) % m)], j] ^= 1
        return Mx

    MA, MB = qconv(sA), qconv(sB)
    HX = np.concatenate([MA, MB], axis=1)
    HZ = np.concatenate([MB.T, MA.T], axis=1)
    return 2 * len(reps) - rank_f2(HX) - rank_f2(HZ)


def half_y(supp, half):
    out = {}
    for (ax, ay) in supp:
        e = (ax, ay % half)
        out[e] = out.get(e, 0) ^ 1
    return {e for e, c in out.items() if c}


def poly_str(supp):
    """bb_lab Poly.from_string format: explicit '*' products, e.g. 'x^2*y^3'."""
    def term(e):
        ex, ey = e
        fs = []
        if ex:
            fs.append("x" if ex == 1 else f"x^{ex}")
        if ey:
            fs.append("y" if ey == 1 else f"y^{ey}")
        return "*".join(fs) if fs else "1"
    return " + ".join(term(e) for e in sorted(supp))


def bezout_deck_y(ell, m, sA, sB):
    """Witness (P,Q) with P*A + Q*B = 1 + y^{m/2} in F2[Z_ell x Z_m], or None."""
    MA, MB = conv_matrix(ell, m, sA), conv_matrix(ell, m, sB)
    elems = [(a, b) for a in range(ell) for b in range(m)]
    idx = {e: i for i, e in enumerate(elems)}
    b = np.zeros(len(elems), dtype=np.uint8)
    b[idx[(0, 0)]] ^= 1
    b[idx[(0, m // 2)]] ^= 1
    z = solve_f2(np.concatenate([MA, MB], axis=1), b)
    if z is None:
        return None
    P, Q = z[:len(elems)], z[len(elems):]
    assert ((MA @ P + MB @ Q) % 2 == b).all()
    return {"P": [elems[i] for i in np.nonzero(P)[0]],
            "Q": [elems[i] for i in np.nonzero(Q)[0]]}


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    report = {"code": "IBM class Y [[288,8,20]] (arXiv:2606.02418 Table II)",
              "rungs": []}
    ell, m, sA, sB = ELL, M, A_SUPP, B_SUPP
    while True:
        k = code_k(ell, m, sA, sB)
        decks = {}
        if ell % 2 == 0:
            decks["x"] = quotient_k(ell, m, sA, sB, (ell // 2, 0))
        if m % 2 == 0:
            decks["y"] = quotient_k(ell, m, sA, sB, (0, m // 2))
        if ell % 2 == 0 and m % 2 == 0:
            decks["xy"] = quotient_k(ell, m, sA, sB, (ell // 2, m // 2))
        rung = {"frame": [ell, m], "n": 2 * ell * m, "k": int(k),
                "deck_k": {kk: int(v) for kk, v in decks.items()},
                "A": poly_str(sA), "B": poly_str(sB),
                "ladder_row": json.dumps({"frame": [ell, m],
                                          "A": poly_str(sA),
                                          "B": poly_str(sB)})}
        if m % 2 == 0:
            wit = bezout_deck_y(ell, m, sA, sB)
            rung["R_y_deck"] = ("HOLDS (Bezout witness saved)"
                               if wit is not None else "FAILS (1+y^{m/2} not in (A,B))")
            if wit is not None:
                (OUT / f"bezout_y_{ell}x{m}.json").write_text(json.dumps(wit))
        report["rungs"].append(rung)
        print(f"({ell},{m})  n={2*ell*m}  k={k}  decks={decks}  "
              f"R(y)={rung.get('R_y_deck', '-')}")
        print(f"   A = {rung['A']}\n   B = {rung['B']}")
        if m % 2 == 0 and decks.get("y") == k:
            m //= 2
            sA, sB = half_y(sA, m), half_y(sB, m)
        else:
            break
    (OUT / "tower_report.json").write_text(json.dumps(report, indent=1))
    print(f"\nwrote {OUT / 'tower_report.json'}")
    print("\nladder rows (feed to a15_coset_distance.py):")
    for rung in report["rungs"]:
        print(f"  n={rung['n']}: {rung['ladder_row']}")


if __name__ == "__main__":
    main()
