"""A28 M-series — f2a6 spectral measurements (falsify-first).

M0: spectral layer self-checks (Fourier multiplicativity, Galois closure,
    |Z_A cap Z_B| = k/2).
M1: empirical patterns of the 113 census classes: exact zero sets of both
    blocks, the visible-Z(f) pattern, orbit structure.
M2: game soundness on real words: bound(Z(u)) <= |u| for every class block.
M3: kill power: game values for Z_A u (small orbit unions) — how fast do
    bounds grow with the hypothesized pattern?
M4: cell-exact minima vs game bounds on the empirical patterns.

Run: uv run --project experiments/bb_lab python experiments/bb_lab/scripts/a28_m_f2a6.py
"""

import sys
import time
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from a28_lsc_lib import REGISTRY, load_f2a6_census
from a28_spectral import Spectral, gf_mul, solve_boundary
from a28_shift_engine import ShiftGame, ideal_basis, min_weight_bz


def mask_of(idxs):
    m = 0
    for i in idxs:
        m |= 1 << i
    return m


def orbits_in(mask, orbit_list):
    return [o for o in orbit_list if mask_of(o) & mask == mask_of(o)]


def main():
    inst = REGISTRY["f2a6"]
    G = inst.G
    sp = Spectral(G)
    game = ShiftGame(sp)
    N = G.N

    # ---- M0 ----
    import random
    rng = random.Random(1)
    for _ in range(20):
        f = rng.getrandbits(N) & ((1 << N) - 1)
        g = rng.getrandbits(N) & ((1 << N) - 1)
        fg = G.mul(f, g)
        Ff, Fg, Ffg = sp.fourier(f), sp.fourier(g), sp.fourier(fg)
        assert all(Ffg[i] == gf_mul(Ff[i], Fg[i]) for i in range(N)), "Fourier not multiplicative"
    ZA = sp.zero_set(inst.A)
    ZB = sp.zero_set(inst.B)
    orbit_list = sp.galois_orbits()
    orb_of = {}
    for o in orbit_list:
        for ci in o:
            orb_of[ci] = o
    for Z in (ZA, ZB):
        assert all(orb_of[ci] <= Z for ci in Z), "zero set not Galois-closed"
    assert len(ZA & ZB) == inst.expect_k // 2, (len(ZA & ZB), inst.expect_k)
    print(f"M0: Fourier ✓  |Z_A|={len(ZA)} |Z_B|={len(ZB)} "
          f"|Z_A∩Z_B|={len(ZA & ZB)}=k/2 ✓  orbits={len(orbit_list)} "
          f"sizes={sorted(Counter(len(o) for o in orbit_list).items())}")

    ZA_m, ZB_m = mask_of(ZA), mask_of(ZB)
    shared = ~(ZA_m | ZB_m) & ((1 << N) - 1)

    # ---- M1 + M2 ----
    census = load_f2a6_census()
    pats = Counter()
    pat_examples = {}
    t0 = time.time()
    viol = 0
    rows_out = []
    for row in census:
        u = G.from_support([(x, y) for blk, x, y in row["b_support"] if blk == 0])
        v = G.from_support([(x, y) for blk, x, y in row["b_support"] if blk == 1])
        Zu, Zv = sp.zero_set(u), sp.zero_set(v)
        Zu_m, Zv_m = mask_of(Zu), mask_of(Zv)
        # consistency: on the shared region, Z(u) and Z(v) agree (both = Z_f)
        assert (Zu_m & shared) == (Zv_m & shared), "shared-region mismatch"
        # visible Z(f): shared part + (Z_v inside Z_A\Z_B) + (Z_u inside Z_B\Z_A)
        vis = (Zu_m & shared) | (Zv_m & ZA_m & ~ZB_m) | (Zu_m & ZB_m & ~ZA_m)
        # Galois-closed?
        vi = {i for i in range(N) if (vis >> i) & 1}
        assert all(orb_of[ci] <= vi for ci in vi), "visible pattern not Galois-closed"
        n_orb = len(orbits_in(vis, orbit_list))
        # M2 game soundness: oracle = exact zero set, bound <= |u|
        bu, _ = game.best_bound(Zu_m, beam=250, max_level=row["u_weight"] + 2)
        bv, _ = game.best_bound(Zv_m, beam=250, max_level=row["v_weight"] + 2)
        if bu > row["u_weight"] or bv > row["v_weight"]:
            viol += 1
            print("  SOUNDNESS VIOLATION", row, bu, bv)
        pats[(bin(vis).count('1'), n_orb)] += 1
        pat_examples.setdefault(vis, []).append(row["b_weight"])
        rows_out.append((row["b_weight"], bin(vis).count("1"), n_orb,
                         row["u_weight"], bu, row["v_weight"], bv))
    assert viol == 0
    print(f"M1: 113 classes, visible-pattern (|Z_f_vis|, #orbits) histogram: "
          f"{dict(sorted(pats.items()))}")
    print(f"    distinct visible patterns: {len(pat_examples)}")
    print(f"M2: game soundness on all 226 block words ✓ (no bound exceeded truth) "
          f"[{time.time()-t0:.0f}s]")
    tight_u = sum(1 for w, _, _, uw, bu, vw, bv in rows_out if bu == uw)
    tight_v = sum(1 for w, _, _, uw, bu, vw, bv in rows_out if bv == vw)
    gaps = [(uw - bu) + (vw - bv) for w, _, _, uw, bu, vw, bv in rows_out]
    print(f"    tight blocks: u {tight_u}/113, v {tight_v}/113; "
          f"pair-gap histogram {sorted(Counter(gaps).items())}")

    # ---- M3: kill power on small hypothesized patterns ----
    print("M3: game values for O = Z_A (+ orbit unions):")
    base_u, _ = game.best_bound(ZA_m, beam=600)
    base_v, _ = game.best_bound(ZB_m, beam=600)
    print(f"    l(Z_A)={base_u}  l(Z_B)={base_v}  pair floor={base_u + base_v} "
          f"(true d(census)=6)")
    free_orbits = [o for o in orbit_list
                   if not (mask_of(o) & (ZA_m | ZB_m))]
    print(f"    free orbits (outside Z_A∪Z_B): "
          f"{sorted(Counter(len(o) for o in free_orbits).items())}")
    t0 = time.time()
    grow = []
    for o in free_orbits:
        om = mask_of(o)
        bu, _ = game.best_bound(ZA_m | om, beam=400)
        bv, _ = game.best_bound(ZB_m | om, beam=400)
        grow.append((len(o), bu, bv, bu + bv))
    print(f"    single-orbit additions: (|orbit|, l_u, l_v, sum) -> "
          f"{sorted(Counter(grow).items())} [{time.time()-t0:.0f}s]")

    # ---- M4: cell-exact minima vs bounds on empirical patterns ----
    print("M4: per empirical pattern: cell-exact min |b| vs game pair bound")
    done = 0
    for vis, weights in sorted(pat_examples.items(),
                               key=lambda kv: bin(kv[0]).count("1"))[:8]:
        Ou = ZA_m | vis
        Ov = ZB_m | vis
        bu, _ = game.best_bound(Ou, beam=600)
        bv, _ = game.best_bound(Ov, beam=600)
        # cell-exact: enumerate the ideal {f : fhat|_vis = 0}, min |del f| over
        # members whose visible pattern is exactly vis
        kb = ideal_basis(sp, vis)
        dim = len(kb)
        best = None
        if dim <= 22:
            import itertools
            for bits in range(1, 1 << dim):
                f = 0
                bb = bits
                i = 0
                while bb:
                    if bb & 1:
                        f ^= kb[i]
                    bb >>= 1
                    i += 1
                uu, vv = inst.boundary(f)
                if uu == 0 and vv == 0:
                    continue
                w = G.weight(uu) + G.weight(vv)
                if best is None or w < best:
                    best = w
        print(f"    |vis|={bin(vis).count('1'):2d} dim(ideal)={dim:2d} "
              f"census weights {sorted(set(weights))} bound {bu}+{bv}={bu+bv} "
              f"cell-min {best}")
        done += 1
        if done >= 8:
            break


if __name__ == "__main__":
    main()
