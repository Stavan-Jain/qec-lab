"""A36 Parts 0-2: the tower slice calculus ported to [[288,12,18]] (bb288).

Tower (A35 screen row `bb288_yxx`, all folds literal descent):

  B36 (3,6) [[36,8,4]] <--x3-- B72 (6,6) [[72,12,6]] <--x6-- GR (12,6)
  [[144,12,12]] (THE GROSS CODE) <--y6-- G8 (12,12) [[288,12,18]]

  G8:  A = x^3 + y^2 + y^7,  B = y^3 + x + x^2   (d = 18 SAT-established,
       arXiv:2308.07915 Table 3; the published solver-exact record)
  GR:  A = x^3 + y + y^2,    B = y^3 + x + x^2   (y^7 -> y mod 6; d = 12
       kernel-checked in QECLean — consumed at Lean tier by the assembly)
  B72: same polys read mod (6,6);  B36: A = 1 + y + y^2 (x^3 -> 1 mod 3)

Rung facts (A35, re-asserted here): top y-rung TWISTED (y^7 carries
sigma = y^6), (R) holds (k 12 -> 12); mid x-rung untwisted (R); bottom
x-rung twisted, (R) FAILS (k 12 -> 8).  Pair (top, mid) is regime R3:
SEAM = im p_y* equals ker p_x* exactly (dim 6), so W = p_x*(SEAM) = 0.

Parts (falsify-first, hard asserts):
  0  frames + all three decks (Lemma 1 constructor asserts); the
     overflow square via the OTHER descent order G8 -> (6,12) -> (6,6)
     (mixed-axis: the square is live on this tower, unlike A33)
  1  parity (Lemma 2, exhaustive via kernel bases, all 4 levels);
     two-level slice/carry on random cycles; converse lifts
  2  H1 rank lattice per rung (rank p* = rank tau* = 6/6/6, exactness
     both ways on the two (R) rungs, cover-side only on the bottom
     rung); SEAM = im p_y*: dim 6, 63 classes, translation-orbit split,
     reps, class-stab orders; the R3 pair identity (SEAM = ker p_x*);
     lab X<->Z transpose-duality spot-check at G8.

Also the shared MODULE for a36_direct_close / a36_descent
(build_tower(), seam_info(), rep helpers re-exported).

Output: data/a36/tower_cells.json
"""

# Provenance: copied verbatim 2026-08-18 (A38 S1) from the unmerged
# branch claude/tower-slice-calculus-generalize-410ed1 (the A35/A36
# session). That branch stays the source of truth until it merges;
# library-grade ports live in bb_lab.tower, not in edits here.


from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np

LAB = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(LAB / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from a30_rung_pass import rref_ints, v2i  # noqa: E402
import a32_tower_slice as TS  # noqa: E402
from a33_tower_cells import h1_map, iota_perm, rep_for, translation_mat  # noqa: E402

DATA = LAB / "data" / "a36"
DATA.mkdir(parents=True, exist_ok=True)

RNG = np.random.default_rng(20260811)

SPECS = {
    "G8": ((12, 12), "x^3 + y^2 + y^7", "y^3 + x + x^2"),
    "GR": ((12, 6), "x^3 + y + y^2", "y^3 + x + x^2"),
    "B72": ((6, 6), "x^3 + y + y^2", "y^3 + x + x^2"),
    "B36": ((3, 6), "1 + y + y^2", "y^3 + x + x^2"),
}


def build_tower():
    """The four codes + three decks (constructor asserts = Lemma 1)."""
    codes = {n: TS.BBCode(n, *s) for n, s in SPECS.items()}
    G8, GR, B72, B36 = (codes[k] for k in ("G8", "GR", "B72", "B36"))
    deck_y = TS.Deck(G8, GR, lambda e: (e[0], e[1] % 6),
                     lambda e, s: (e[0], e[1] + 6 * s))
    deck_x = TS.Deck(GR, B72, lambda e: (e[0] % 6, e[1]),
                     lambda e, s: (e[0] + 6 * s, e[1]))
    deck_b = TS.Deck(B72, B36, lambda e: (e[0] % 3, e[1]),
                     lambda e, s: (e[0] + 3 * s, e[1]))
    return G8, GR, B72, B36, deck_y, deck_x, deck_b


def seam_info(GR: TS.BBCode, My: np.ndarray):
    """SEAM = im p_y* on H1(GR): basis, 63 classes, orbit split, reps,
    class-stab orders (generic version of a33's seam_data — no
    instance-specific orbit-size assert; the caller asserts)."""
    Sb = TS._colspace(My)
    Sbb, Sbp = rref_ints(list(Sb))
    pts = TS._span_points(Sbb, Sbp) - {0}
    mats = TS._translation_action(GR)
    orbs = sorted((sorted(o) for o in TS._orbits(pts, mats)),
                  key=lambda o: (len(o), o[0]))
    reps = [min(o) for o in orbs]
    stab_orders = []
    for rep in reps:
        cnt = sum(1 for M in mats if TS._apply(M, rep) == rep)
        stab_orders.append(cnt)
    return {"basis": Sbb, "piv": Sbp, "classes": sorted(pts),
            "orbits": orbs, "reps": reps, "stab_orders": stab_orders}


def main():
    t0 = time.monotonic()
    out: dict = {}

    # ------------------------------------------------------------- Part 0
    G8, GR, B72, B36, deck_y, deck_x, deck_b = build_tower()
    assert (G8.k, GR.k, B72.k, B36.k) == (12, 12, 12, 8), \
        (G8.k, GR.k, B72.k, B36.k)
    print(f"[{time.monotonic()-t0:5.1f}s] Part 0: codes built, "
          f"k = 12/12/12/8; Lemma 1 constructor asserts PASS on all "
          f"three rungs")
    # the OTHER descent order (x first): G8 -> (6,12) -> (6,6) — the
    # overflow square is live (mixed-axis tower)
    GX = TS.BBCode("GX", (6, 12), "x^3 + y^2 + y^7", "y^3 + x + x^2")
    deck_xc = TS.Deck(G8, GX, lambda e: (e[0] % 6, e[1]),
                      lambda e, s: (e[0] + 6 * s, e[1]))
    deck_yb = TS.Deck(GX, B72, lambda e: (e[0], e[1] % 6),
                      lambda e, s: (e[0], e[1] + 6 * s))
    lhs = (deck_x.P @ deck_y.P) % 2
    rhs = (deck_yb.P @ deck_xc.P) % 2
    assert (lhs == rhs).all(), "folds do not commute"
    assert GX.k == 12, GX.k
    print(f"[{time.monotonic()-t0:5.1f}s]   other route G8->(6,12)->(6,6)"
          f" built (k(GX) = 12); folds commute")
    out["part0"] = {"k": [12, 12, 12, 8], "k_GX": 12,
                    "folds_commute": True}

    # ------------------------------------------------------------- Part 1
    for code in (G8, GR, B72, B36, GX):
        odd = [int(kv.sum()) % 2 for kv in code.kerHZ]
        assert not any(odd), f"odd kernel basis vector in {code.name}"
    print(f"[{time.monotonic()-t0:5.1f}s] Part 1: parity exhaustive "
          f"(kernel bases even, all 5 frames) => every cycle even")

    n_sq = 0
    for _ in range(200):
        v = G8.random_cycle()
        b, m1, _ = deck_y.slice_data(v)
        assert GR.is_cycle(b)
        beta, m2, _ = deck_x.slice_data(b)
        assert B72.is_cycle(beta)
        assert int(v.sum()) == int(beta.sum()) + 2 * (m1 + m2)
        bx, mx, _ = deck_xc.slice_data(v)
        beta2, my2, _ = deck_yb.slice_data(bx)
        assert (beta2 == beta).all(), "composite shadows differ"
        assert m1 + m2 == mx + my2, "overflow square fails"
        n_sq += 1
    # converse lifts: random GROSS cycles are almost never liftable here
    # (liftable classes = im p_y*, density 2^-6 — the A35 delta law), so
    # draw shadows OF cover cycles instead: those must always reconstruct.
    n_conv = n_unlift = 0
    for _ in range(60):
        b = (deck_y.P @ G8.random_cycle()) % 2
        rhs_v = (deck_y.RHS @ b) % 2
        aug = np.concatenate([deck_y.E, rhs_v[:, None]], axis=1)
        R, piv = TS._rref_np(aug)
        assert not any(p == deck_y.E.shape[1] for p in piv), \
            "shadow of a cover cycle has inconsistent carry system?!"
        v0 = np.zeros(deck_y.E.shape[1], dtype=np.uint8)
        for i, p in enumerate(piv):
            v0[p] = R[i, -1]
        v = deck_y.lift(v0, b)
        assert G8.is_cycle(v), "reconstructed lift not a cycle"
        n_conv += 1
    for _ in range(60):   # and the delta-obstruction side: random cycles
        b = GR.random_cycle()
        rhs_v = (deck_y.RHS @ b) % 2
        aug = np.concatenate([deck_y.E, rhs_v[:, None]], axis=1)
        _, piv = TS._rref_np(aug)
        if any(p == deck_y.E.shape[1] for p in piv):
            n_unlift += 1
    print(f"[{time.monotonic()-t0:5.1f}s]   {n_sq} two-level slice/carry/"
          f"square checks; converse lifts {n_conv}/60 on true shadows; "
          f"{n_unlift}/60 random gross cycles obstructed (delta law: "
          f"expect ~ 60*(1 - 2^-6) = 59)")
    out["part1"] = {"square_checks": n_sq, "converse_lifts": n_conv,
                    "random_cycles_obstructed": n_unlift}

    # ------------------------------------------------------------- Part 2
    My = h1_map(deck_y)      # p_y*: H1(G8) -> H1(GR), 12x12
    Mx = h1_map(deck_x)      # p_x*: H1(GR) -> H1(B72), 12x12
    Mb = h1_map(deck_b)      # p_b*: H1(B72) -> H1(B36), 8x12
    Ty = h1_map(deck_y, tau=True)
    Tx = h1_map(deck_x, tau=True)
    Tb = h1_map(deck_b, tau=True)
    ranks = {"py": TS.gf2_rank([v2i(c) for c in My.T]),
             "px": TS.gf2_rank([v2i(c) for c in Mx.T]),
             "pb": TS.gf2_rank([v2i(c) for c in Mb.T]),
             "ty": TS.gf2_rank([v2i(c) for c in Ty.T]),
             "tx": TS.gf2_rank([v2i(c) for c in Tx.T]),
             "tb": TS.gf2_rank([v2i(c) for c in Tb.T])}
    assert ranks == {"py": 6, "px": 6, "pb": 6,
                     "ty": 6, "tx": 6, "tb": 6}, ranks
    # exactness: cover side on all three rungs; base side on (R) rungs only
    for name, Mp, Mt, base_exact in [("y", My, Ty, True),
                                     ("x", Mx, Tx, True),
                                     ("b", Mb, Tb, False)]:
        assert TS._span_eq(TS._colspace(Mt), TS._kernel_ints(Mp)), \
            f"im tau* != ker p* on rung {name}"
        got = TS._span_eq(TS._kernel_ints(Mt), TS._colspace(Mp))
        assert got == base_exact, f"base exactness on rung {name}: {got}"
    # deck-triviality: sigma* = id on the two (R) rungs, != id on bottom
    assert (translation_mat(G8, (0, 6))
            == np.eye(12, dtype=np.uint8)).all(), "sigma_y* != id"
    assert (translation_mat(GR, (6, 0))
            == np.eye(12, dtype=np.uint8)).all(), "sigma_x* != id"
    assert not (translation_mat(B72, (3, 0))
                == np.eye(12, dtype=np.uint8)).all(), "sigma_b* = id?!"
    print(f"[{time.monotonic()-t0:5.1f}s] Part 2: rank p* = rank tau* = 6 "
          f"all rungs; exactness cover-side 3/3, base-side on (R) rungs "
          f"y,x only; sigma* = id on y,x; != id on bottom  [A35 row "
          f"REPRODUCED]")
    out["part2_ranks"] = ranks

    # SEAM = im p_y*: 63 classes, orbit split under the 72 translations
    sd = seam_info(GR, My)
    assert len(sd["basis"]) == 6 and len(sd["classes"]) == 63
    osz = [len(o) for o in sd["orbits"]]
    print(f"[{time.monotonic()-t0:5.1f}s]   SEAM = im p_y*: dim 6, 63 "
          f"classes in {len(osz)} translation orbits, sizes {osz}, "
          f"class-stab orders {sd['stab_orders']}")
    out["seam"] = {"dim": 6, "classes": 63, "orbit_sizes": osz,
                   "reps": [f"{r:#x}" for r in sd["reps"]],
                   "stab_orders": sd["stab_orders"]}

    # the R3 pair identity: SEAM == ker p_x* exactly (=> W = p_x*(SEAM) = 0)
    Kx = TS._kernel_ints(Mx)
    assert TS._span_eq(sd["basis"], Kx), "SEAM != ker p_x* (R3 fails?)"
    Wimg = sorted({TS._apply(Mx, s) for s in
                   TS._span_points(sd["basis"], sd["piv"])} - {0})
    assert not Wimg, f"W nonzero: {len(Wimg)} classes"
    print(f"[{time.monotonic()-t0:5.1f}s]   R3 pair identity: SEAM = "
          f"ker p_x* (dim 6); W = p_x*(SEAM) = 0 — sector A EMPTY one "
          f"rung down  [A35 pair(0,1) REPRODUCED]")
    out["pair_R3"] = {"seam_eq_ker_px": True, "dim_W": 0}

    # X<->Z transpose duality spot-check at G8 (antipode + block swap)
    iota8 = iota_perm(G8)
    swap8 = np.concatenate([np.arange(G8.ng, 2 * G8.ng),
                            np.arange(0, G8.ng)])
    dual = lambda v: v[iota8][swap8]  # noqa: E731
    for kv in G8.kerHZ[:20]:
        assert not ((G8.HX @ dual(kv)) % 2).any(), "duality: cycle fails"
    for row in G8.HX[::12]:
        assert TS.in_span(v2i(dual(row)),
                          *rref_ints([v2i(r) for r in G8.HZ])), \
            "duality: stab row fails"
    print(f"[{time.monotonic()-t0:5.1f}s]   X<->Z duality spot-check OK "
          f"=> Z-side floors follow from X-side")
    out["duality"] = "spot-checked at G8"

    out["wall_s"] = round(time.monotonic() - t0, 1)
    (DATA / "tower_cells.json").write_text(json.dumps(out, indent=1))
    print(f"total {out['wall_s']}s -> {DATA / 'tower_cells.json'}")


if __name__ == "__main__":
    main()
