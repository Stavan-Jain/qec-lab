"""A28 diagnostics — is the shift-game weakness intrinsic or a search artifact?

D1: exhaustive BFS game value (no beam) vs beam value, small oracles.
D2: TRUE cell minima: for sample cells (P_in, P_out), enumerate the ideal
    {f : fhat|_(P_in u kernel) = 0} with the BZ C kernel at W = 14 and check
    which cells actually contain light boundaries.  Gap anatomy:
      true cell min > 14, game bound <= 14  -> technique gap (kill was
                                               possible for SOME method)
      exhaustive value > beam value         -> search gap only.
"""

import subprocess
import sys
import time
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from a28_lsc_lib import DATA, REGISTRY, rref
from a28_spectral import Spectral
from a28_shift_engine import ShiftGame, ideal_basis
from a28_bz_census import C_SRC, BUILD, write_mat


def mask_of(idxs):
    m = 0
    for i in idxs:
        m |= 1 << i
    return m


def exhaustive_value(game, oracle, max_level=8, cap_states=250_000):
    """True BFS with full dedupe by translation-canonical fingerprint."""
    states = {(): 0}
    best = 0
    N = game.N
    for level in range(max_level):
        nxt = {}
        for fp in states:
            # reconstruct one representative mask from the fingerprint is
            # impossible; store masks instead
            pass
        break
    # simpler: BFS storing canonical masks
    seen = {0}
    frontier = [0]
    closer_pool = ~oracle & ((1 << N) - 1)
    for level in range(max_level):
        nfront = set()
        for mask in frontier:
            for mu in range(N):
                tm = game.translate_mask(mask, mu) if mask else 0
                if tm & ~oracle:
                    continue
                pool = closer_pool & ~tm
                while pool:
                    low = pool & -pool
                    pool ^= low
                    nm = tm | low
                    # canonicalize
                    cbest = None
                    for m2 in range(N):
                        t2 = game.translate_mask(nm, m2)
                        if cbest is None or t2 < cbest:
                            cbest = t2
                    if cbest not in nfront:
                        nfront.add(cbest)
                if not mask:
                    break
        if not nfront:
            return best
        best = level + 1
        if len(nfront) > cap_states:
            return best, "state-cap-hit"
        frontier = list(nfront)
    return best


def main():
    inst = REGISTRY["f2a6"]
    sp = Spectral(inst.G)
    game = ShiftGame(sp)
    N = inst.G.N
    ZA = mask_of(sp.zero_set(inst.A))
    ZB = mask_of(sp.zero_set(inst.B))
    kernel = ZA & ZB
    orbits = [mask_of(o) for o in sp.galois_orbits()]
    free_orbits = [o for o in orbits if not (o & kernel)]

    # ---- D1: exhaustive vs beam on Z_A and Z_A + one orbit ----
    for tag, oracle in [("Z_A", ZA)] + [
            (f"Z_A+orb{j}", ZA | free_orbits[j]) for j in (0, 5, 11)]:
        t0 = time.time()
        ex = exhaustive_value(game, oracle, max_level=9)
        bm, _ = game.best_bound(oracle, beam=600)
        print(f"D1 {tag}: exhaustive={ex} beam={bm} [{time.time()-t0:.0f}s]")

    # ---- D2: true cell minima for sample 2-orbit and 3-orbit P_in cells ----
    src = BUILD / "bzkern.c"
    binp = BUILD / "bzkern"
    BUILD.mkdir(parents=True, exist_ok=True)
    if not binp.exists():
        src.write_text(C_SRC)
        subprocess.run(["cc", "-O3", "-o", str(binp), str(src), "-lpthread"],
                       check=True)

    import itertools
    samples = list(itertools.combinations(range(len(free_orbits)), 2))[:6] + \
        list(itertools.combinations(range(len(free_orbits)), 3))[:4]
    for orbsel in samples:
        P_in = 0
        for j in orbsel:
            P_in |= free_orbits[j]
        # ideal mod K: fhat = 0 on P_in and on the kernel
        kb = ideal_basis(sp, P_in | kernel)
        # boundary images as packed rows
        rows = [inst.pack(inst.boundary(f)) for f in kb]
        basis, piv = rref(rows, 2 * N)
        kappa = len(basis)
        mat = BUILD / "diag.mat"
        write_mat(mat, basis)
        # enumerate: any codeword of the SUBcode with weight <= 14; single
        # info set: coverage r = ... we need min over subcode: escalate r
        # until floor r+1 > 14 (r = 14) — too big; instead r=4 probe:
        # presence of light words is what we need, plus exactness via the
        # 2-disjoint-set trick is overkill here; use r s.t. C(kappa,r) ~ 1e8.
        r = 4
        res = subprocess.run([str(binp), str(mat), str(r), "14", "8",
                              str(BUILD / "diag")],
                             capture_output=True, text=True, check=True)
        hits = 0
        light = Counter()
        for t in range(8):
            for line in open(str(BUILD / f"diag_t{t:02d}.hits")):
                w = bin(int(line, 16)).count("1")
                light[w] += 1
                hits += 1
        bu, _ = game.best_bound(ZA | P_in, beam=400)
        bv, _ = game.best_bound(ZB | P_in, beam=400)
        print(f"D2 P_in=orbits{orbsel} dim_ideal(modK)={len(kb)} "
              f"game(union)={bu}+{bv}={bu+bv} "
              f"light@r{r}: {dict(sorted(light.items())) if hits else 'NONE'}")


if __name__ == "__main__":
    main()
