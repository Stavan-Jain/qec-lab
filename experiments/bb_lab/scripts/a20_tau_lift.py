"""A20: tau-lift the Y4 weight-10 witness to a weight-20 X-logical of Y8.

Reads the exact-ladder witness for Y4 = [[144,8,10]] (data/a20/y144_ladder.log,
final JSON line), verifies it is a weight-10 X-logical of Y4 under bb_lab's own
conventions, transfer-lifts it (both sheets of the y-fiber) to Y8 = (18,8), and
verifies the lift is a weight-20 nontrivial X-logical: **d(Y8) <= 20
constructive**, independent of IBM's MILP incumbent — and the tau-tightness
witness for the rung-2 doubling template (the dangerous floor is attained).

Usage (from experiments/bb_lab):  uv run scripts/a20_tau_lift.py
"""
import json
import sys
from pathlib import Path

import numpy as np

LAB = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(LAB / "src"))

from bb_lab.checks import bb_check_matrices
from bb_lab.group import AbelianGroup
from bb_lab.linalg import rank_f2
from bb_lab.poly import Poly

OUT = LAB / "data" / "a20"

Y4 = {"frame": (18, 4), "A": "1 + x + x^14*y", "B": "1 + x*y^2 + x^2*y^3"}
Y8 = {"frame": (18, 8), "A": "1 + x*y^4 + x^14*y", "B": "1 + x*y^2 + x^2*y^7"}


def checks_of(spec):
    G = AbelianGroup(spec["frame"])
    A = Poly.from_string(spec["A"], G)
    B = Poly.from_string(spec["B"], G)
    return G, bb_check_matrices(A, B)


def is_x_logical(checks, v):
    """X-operator v: commutes with Z-checks (ker H_Z) and is not an
    X-stabilizer (rowspace H_X). The ladder's witnesses are X-side;
    d_X = d_Z by inversion duality."""
    if (checks.H_Z @ v % 2).any():
        return False, "not in ker H_Z"
    hx = checks.H_X % 2
    if rank_f2(np.vstack([hx, v])) == rank_f2(hx):
        return False, "in rowspace(H_X) (stabilizer)"
    return True, "nontrivial X-logical"


def main():
    witness = None
    for line in (OUT / "y144_ladder.log").read_text().splitlines():
        try:
            row = json.loads(line)
            witness = row["witness"]
        except (json.JSONDecodeError, KeyError):
            continue
    if witness is None:
        sys.exit("no witness JSON in data/a20/y144_ladder.log")

    G4, c4 = checks_of(Y4)
    els4 = list(G4)
    n4 = c4.H_X.shape[1]
    v4 = np.zeros(n4, dtype=np.uint8)
    v4[witness] = 1
    okay, why = is_x_logical(c4, v4)
    print(f"Y4 witness weight {int(v4.sum())}: {why}")
    if not okay:
        sys.exit("Y4 witness failed verification — indexing convention?")

    # tau-lift: qubit (block, (a,b)) of Y4 -> both fiber qubits (block, (a,b)),
    # (block, (a,b+4)) of Y8. bb_lab qubit order: [left block | right block],
    # group elements enumerated by AbelianGroup iteration order.
    G8, c8 = checks_of(Y8)
    els8 = list(G8)
    e8 = {g: i for i, g in enumerate(els8)}
    half4 = len(els4)
    half8 = len(els8)
    v8 = np.zeros(c8.H_X.shape[1], dtype=np.uint8)
    for q in np.nonzero(v4)[0]:
        blk, gi = divmod(int(q), half4)
        a, b = els4[gi]
        for bb in (b, b + 4):
            v8[blk * half8 + e8[(a, bb % 8)]] ^= 1
    okay, why = is_x_logical(c8, v8)
    print(f"tau-lift weight {int(v8.sum())}: {why}")
    if okay and int(v8.sum()) == 20:
        print("VERIFIED: d(Y8) <= 20 constructive; tau-tightness witness "
              "for the rung-2 template (dangerous floor attained at 2*d1).")
        np.save(OUT / "y8_weight20_witness.npy", v8)
        print(f"saved {OUT / 'y8_weight20_witness.npy'}")
    else:
        sys.exit("lift failed")


if __name__ == "__main__":
    main()
