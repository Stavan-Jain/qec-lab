"""A28 diagnostics — (a) does the syzygy joint game beat block-additivity?
(b) certify the (0,1,5) technique-gap cell rigorously (2-set BZ on its ideal).

Run: uv run --project experiments/bb_lab python experiments/bb_lab/scripts/a28_diag_joint.py
"""

import itertools
import math
import subprocess
import sys
import time
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from a28_lsc_lib import REGISTRY, rref
from a28_spectral import Spectral
from a28_shift_engine import ShiftGame, ideal_basis
from a28_joint_game import JointGame
from a28_bz_census import C_SRC, BUILD, write_mat


def mask_of(idxs):
    m = 0
    for i in idxs:
        m |= 1 << i
    return m


def main():
    inst = REGISTRY["f2a6"]
    sp = Spectral(inst.G)
    game = ShiftGame(sp)
    joint = JointGame(sp, sp.fourier(inst.A), sp.fourier(inst.B))
    N = inst.G.N
    ZA = mask_of(sp.zero_set(inst.A))
    ZB = mask_of(sp.zero_set(inst.B))
    kernel = ZA & ZB
    orbits = [mask_of(o) for o in sp.galois_orbits()]
    free_orbits = [o for o in orbits if not (o & kernel)]
    all_free = 0
    for o in free_orbits:
        all_free |= o

    # (a) joint vs union on sample cells (P_in = orbit unions, P_out = rest)
    print("cell (P_in orbits | P_out = rest): union vs joint")
    samples = [(), (0,), (5,), (0, 1), (0, 5), (0, 1, 2), (0, 1, 5),
               (0, 1, 2, 3), (2, 7, 11)]
    for sel in samples:
        P_in = 0
        for j in sel:
            P_in |= free_orbits[j]
        P_out = all_free & ~P_in
        t0 = time.time()
        bu, _ = game.best_bound(ZA | P_in, free=P_out, beam=500, max_level=15)
        bv, _ = game.best_bound(ZB | P_in, free=P_out, beam=500, max_level=15)
        bj, _ = joint.best_bound(P_in, P_out, target=15, beam=400)
        print(f"  P_in={sel!s:14} union={bu}+{bv}={bu+bv:2d} joint={bj:2d} "
              f"delta={bj - bu - bv:+d}  [{time.time()-t0:.0f}s]")

    # (b) rigorous emptiness for the (0,1,5) subtree: BZ with two disjoint
    # info sets on the ideal's boundary image (certified: no nonzero boundary
    # of weight <= 14 with fhat vanishing on orbits {0,1,5} + kernel)
    P_in = free_orbits[0] | free_orbits[1] | free_orbits[5]
    kb = ideal_basis(sp, P_in | kernel)
    rows = [inst.pack(inst.boundary(f)) for f in kb]
    ncols = 2 * N
    b1, p1 = rref(rows, ncols)
    kappa = len(b1)
    comp = [c for c in range(ncols) if c not in set(p1)]
    b2, p2 = rref(rows, ncols, col_order=comp)
    print(f"(b) ideal(orbits 0,1,5): dim={len(kb)} rank(boundary)={kappa} "
          f"second-set rank={len(b2)}")
    assert len(b2) == kappa, "need fallback info sets"
    W, r = 14, 7
    exp_nodes = sum(math.comb(kappa, s) for s in range(1, r + 1))
    BUILD.mkdir(parents=True, exist_ok=True)
    binp = BUILD / "bzkern"
    if not binp.exists():
        (BUILD / "bzkern.c").write_text(C_SRC)
        subprocess.run(["cc", "-O3", "-o", str(binp), str(BUILD / "bzkern.c"),
                        "-lpthread"], check=True)
    tot_hits = 0
    for j, basis in enumerate((b1, b2), 1):
        mat = BUILD / f"gapcell_G{j}.mat"
        write_mat(mat, basis)
        res = subprocess.run([str(binp), str(mat), str(r), str(W), "8",
                              str(BUILD / f"gapcell_G{j}")],
                             capture_output=True, text=True, check=True)
        parts = dict(kv.split("=") for kv in res.stdout.strip().split())
        assert int(parts["nodes"]) == exp_nodes
        tot_hits += int(parts["hits"])
        print(f"    set {j}: {res.stdout.strip()}")
    print(f"    CERTIFIED: subtree ideal has {'NO' if tot_hits == 0 else tot_hits} "
          f"boundaries of weight <= {W} (game bound there was 8 — "
          f"technique gap {'CONFIRMED' if tot_hits == 0 else 'REFUTED'})")


if __name__ == "__main__":
    main()
