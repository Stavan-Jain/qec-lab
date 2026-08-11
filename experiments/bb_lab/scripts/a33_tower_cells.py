"""A33 Parts 0-2: the A32 tower slice calculus ported to IBM class Y.

Tower (A20 SS1, all-(R), SAME-AXIS — both rungs are y-decks):

  Y2 (18,2) [[72,8,6]] --y--> Y4 (18,4) [[144,8,10]] --y--> Y8 (18,8) [[288,8,20]]
  Y2: A = 1+x+x^14 y,      B = 1+x+x^2 y
  Y4: A = 1+x+x^14 y,      B = 1+x y^2+x^2 y^3
  Y8: A = 1+x y^4+x^14 y,  B = 1+x y^2+x^2 y^7

Both upper rungs are TWISTED lifts (Y8's x y^4 and x^2 y^7 terms carry the
deck sigma = y^4 of Z8 -> Z4; Y4's x y^2 / x^2 y^3 carry y^2 of Z4 -> Z2).
Lemma 1 of A32 is twist-generic; every transport is asserted numerically
here exactly as A32 did (TS.Deck's constructor asserts).

Parts (falsify-first, hard asserts):
  0  frames + decks (chain maps, stab transport, twist-invariance,
     sections); composite-fold consistency (SAME-AXIS: the overflow
     square of A32 Thm 3 DEGENERATES — there is only one descent order,
     Y8 -> Y4 -> Y2; asserted as composite-shadow consistency instead)
  1  parity (Lemma 2, exhaustive via kernel bases), two-level
     slice/carry identities on random cycles, converse lifts
  2  the H1 rank lattice under (R): p*/tau* ranks (4 = k/2 each rung),
     exactness (im tau* = ker p*, ker tau* = im p*), deck-triviality
     (sigma* = id, = H4 re-measured), SEAM := im p2* (15 classes, 2
     G-orbits 12+3, class-stab orders 6/24 — A20 SS5 reproduced), the
     rung-1 lattice (SEAM1, K1, W2 := p1*(SEAM) — the descent
     reachability datum), the banked Bezout witnesses re-verified as
     convolution identities, the seamC dictionary (A20's ker-d2
     parametrization == im p2* classes, via the antipode iota), and the
     lab X<->Z transpose-duality check.

Also doubles as the shared MODULE for the other a33 scripts
(build_tower(), h1_map(), seam_data(), iota_perm(), ...).

Output: data/a33/tower_cells.json
Read-only inputs: MAIN checkout data/a20 (never written).
"""

from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np

LAB = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(LAB / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from bb_lab.linalg import nullspace_f2  # noqa: E402

from a30_rung_pass import i2v, reduce_int, rref_ints, v2i  # noqa: E402
import a32_tower_slice as TS  # noqa: E402

MAIN = Path("/Users/stavanjain/Code/qec-lab/experiments/bb_lab/data")
DATA = LAB / "data" / "a33"
DATA.mkdir(parents=True, exist_ok=True)

RNG = np.random.default_rng(20260811)

SPECS = {
    "Y2": ((18, 2), "1 + x + x^14*y", "1 + x + x^2*y"),
    "Y4": ((18, 4), "1 + x + x^14*y", "1 + x*y^2 + x^2*y^3"),
    "Y8": ((18, 8), "1 + x*y^4 + x^14*y", "1 + x*y^2 + x^2*y^7"),
}

# chain-frame (a20_seam_floor / QECLean BBChainComplex) polynomial supports
LX, LY, CLY = 18, 4, 8
A4S = [(0, 0), (1, 0), (14, 1)]
B4S = [(0, 0), (1, 2), (2, 3)]
A8S = [(0, 0), (1, 4), (14, 1)]
B8S = [(0, 0), (1, 2), (2, 7)]


def build_tower():
    """Build the three codes + both decks (constructor asserts = Lemma 1)."""
    codes = {name: TS.BBCode(name, *spec) for name, spec in SPECS.items()}
    Y2, Y4, Y8 = codes["Y2"], codes["Y4"], codes["Y8"]
    deck_top = TS.Deck(Y8, Y4, lambda e: (e[0], e[1] % 4),
                       lambda e, s: (e[0], e[1] + 4 * s))
    deck_bot = TS.Deck(Y4, Y2, lambda e: (e[0], e[1] % 2),
                       lambda e, s: (e[0], e[1] + 2 * s))
    return Y2, Y4, Y8, deck_top, deck_bot


def h1_map(deck: TS.Deck, tau: bool = False) -> np.ndarray:
    """Matrix of p_* (or tau_*): H1(cover) -> H1(base) (or reverse)."""
    src = deck.cover if not tau else deck.base
    dst = deck.base if not tau else deck.cover
    S = np.array([src.sig(r) for r in src.xreps], dtype=np.uint8)
    op = (deck.P if not tau else deck.TAU)
    D = np.array([dst.sig((op @ r) % 2) for r in src.xreps], dtype=np.uint8)
    return (D.T @ TS._gf2_inv(S.T)) % 2


def translation_mat(code: TS.BBCode, t) -> np.ndarray:
    """Sig-space matrix of the translation by t acting on H1."""
    S = np.array([code.sig(r) for r in code.xreps], dtype=np.uint8)
    Sinv = TS._gf2_inv(S.T)
    perm = TS._perm_for(code, t)
    D = np.array([code.sig(r[perm]) for r in code.xreps], dtype=np.uint8)
    return (D.T @ Sinv) % 2


def rep_for(code: TS.BBCode, sig_int: int) -> np.ndarray:
    """A cycle with the prescribed H1 signature."""
    S = np.array([code.sig(r) for r in code.xreps], dtype=np.uint8)
    SinvT = TS._gf2_inv(S.T)
    tvec = i2v(sig_int, code.k)
    coeff = (SinvT @ tvec) % 2
    v = np.zeros(code.n, dtype=np.uint8)
    for i in range(code.k):
        if coeff[i]:
            v ^= code.xreps[i]
    assert v2i(code.sig(v)) == sig_int
    return v


def iota_perm(code: TS.BBCode) -> np.ndarray:
    """Antipode permutation on C1 (per block, g -> -g): the chain-frame
    <-> lab-frame conjugation."""
    ng = code.ng
    perm = np.zeros(2 * ng, dtype=np.int64)
    for i, e in enumerate(code.G):
        j = code.G.index(code.G.neg(e))
        perm[i] = j
        perm[ng + i] = ng + j
    return perm


# ------------------------- chain-frame seamC (a20_seam_floor, verbatim port)
def conv_matrix(supp):
    M = np.zeros((LX * LY, LX * LY), dtype=np.uint8)
    for ga in range(LX):
        for gb in range(LY):
            g = ga * LY + gb
            for (sa, sb) in supp:
                M[((ga + sa) % LX) * LY + ((gb + sb) % LY), g] ^= 1
    return M


def conv_cover(supp, f):
    out = np.zeros_like(f)
    for (sa, sb) in supp:
        out ^= np.roll(np.roll(f, sa, axis=0), sb, axis=1)
    return out


def seam_c(zeta):
    """seamC zeta as a C1(Y4) chain-frame vector (length 144)."""
    lift = np.zeros((LX, CLY), dtype=np.uint8)
    lift[:, :LY] = zeta.reshape(LX, LY)
    cA = conv_cover(A8S, lift)
    cB = conv_cover(B8S, lift)
    return np.concatenate([cA[:, LY:].reshape(-1), cB[:, LY:].reshape(-1)])


def chain_kernel_classes():
    """The 15 nonzero ker d2(Y4) elements in the chain frame, indexed by
    the a20 basis masks (same nullspace_f2 basis order as the banked run)."""
    MA, MB = conv_matrix(A4S), conv_matrix(B4S)
    D2 = np.vstack([MA, MB])
    K = nullspace_f2(D2)
    assert K.shape[0] == 4, f"dim ker d2 = {K.shape[0]} != 4"
    elts = {}
    for mask in range(1, 16):
        z = np.zeros(LX * LY, dtype=np.uint8)
        for i in range(4):
            if (mask >> i) & 1:
                z ^= K[i]
        elts[mask] = z
    return D2, elts


def seam_data(Y4: TS.BBCode, M2: np.ndarray):
    """SEAM = im p2* on H1(Y4): basis, the 15 classes, orbit split, reps.

    Returns dict with basis ints, class list, orbits (as sorted lists),
    orbit rep class ints (min per orbit), class-stab orders.
    """
    Sb = TS._colspace(M2)
    Sbb, Sbp = rref_ints(list(Sb))
    pts = TS._span_points(Sbb, Sbp) - {0}
    assert len(Sbb) == 4 and len(pts) == 15
    mats = TS._translation_action(Y4)
    orbs = TS._orbits(pts, mats)
    orbs = sorted((sorted(o) for o in orbs), key=len, reverse=True)
    assert sorted(len(o) for o in orbs) == [3, 12], \
        f"seam orbit sizes {[len(o) for o in orbs]} != [12, 3]"
    reps = [min(o) for o in orbs]
    stab_orders = []
    for rep in reps:
        cnt = 0
        for t, M in zip(Y4.G, mats):
            if TS._apply(M, rep) == rep:
                cnt += 1
        stab_orders.append(cnt)
    return {"basis": Sbb, "piv": Sbp, "classes": sorted(pts),
            "orbits": orbs, "reps": reps, "stab_orders": stab_orders,
            "mats": mats}


def main():
    t0 = time.monotonic()
    out: dict = {}

    # ------------------------------------------------------------- Part 0
    Y2, Y4, Y8, deck_top, deck_bot = build_tower()
    assert (Y2.k, Y4.k, Y8.k) == (8, 8, 8), (Y2.k, Y4.k, Y8.k)
    print(f"[{time.monotonic()-t0:5.1f}s] Part 0: frames built, k = 8/8/8; "
          f"Deck constructor asserts (chain map, stab transport, "
          f"twist-invariance, sections) PASS on both rungs")
    # composite fold Y8 -> Y2 (y mod 2 from Z8) factors through Y4: the
    # SAME-AXIS degeneracy — there is no second descent order, so A32's
    # overflow square is vacuous here; assert the composite instead.
    P_comp = (deck_bot.P @ deck_top.P) % 2
    Pc_direct = TS.fold_matrix(Y8.G, Y2.G, lambda e: (e[0], e[1] % 2))
    assert (P_comp == Pc_direct).all(), "composite fold mismatch"
    print(f"[{time.monotonic()-t0:5.1f}s]   composite fold = the unique "
          f"Z4-deck fold Y8->Y2 (overflow square DEGENERATE: one path)")
    out["part0"] = {"k": [8, 8, 8], "composite_fold": "consistent",
                    "overflow_square": "degenerate (same-axis, one path)"}

    # ------------------------------------------------------------- Part 1
    for code in (Y2, Y4, Y8):
        odd = [int(kv.sum()) % 2 for kv in code.kerHZ]
        assert not any(odd), f"odd kernel basis vector in {code.name}"
    print(f"[{time.monotonic()-t0:5.1f}s] Part 1: parity lemma exhaustive "
          f"(kernel bases even, Y2/Y4/Y8) => every cycle weight is even")

    n_checks = 0
    for _ in range(200):
        v = Y8.random_cycle()
        b, m1, v0 = deck_top.slice_data(v)
        assert Y4.is_cycle(b)
        assert (((deck_top.E @ v0) + (deck_top.RHS @ b)) % 2 == 0).all()
        beta, m2, b0 = deck_bot.slice_data(b)
        assert Y2.is_cycle(beta)
        assert (((deck_bot.E @ b0) + (deck_bot.RHS @ beta)) % 2 == 0).all()
        assert int(v.sum()) == int(beta.sum()) + 2 * (m1 + m2)
        assert (beta == (P_comp @ v) % 2).all()
        n_checks += 1
    n_conv = 0
    for _ in range(60):
        b = Y4.random_cycle()
        rhs = (deck_top.RHS @ b) % 2
        aug = np.concatenate([deck_top.E, rhs[:, None]], axis=1)
        R, piv = TS._rref_np(aug)
        if any(p == deck_top.E.shape[1] for p in piv):
            continue
        v0 = np.zeros(deck_top.E.shape[1], dtype=np.uint8)
        for i, p in enumerate(piv):
            v0[p] = R[i, -1]
        v = deck_top.lift(v0, b)
        assert Y8.is_cycle(v), "reconstructed lift not a cycle"
        n_conv += 1
    print(f"[{time.monotonic()-t0:5.1f}s]   {n_checks} two-level "
          f"slice/carry checks, {n_conv} converse lifts OK")
    out["part1"] = {"random_cycles": n_checks, "converse_lifts": n_conv}

    # ------------------------------------------------------------- Part 2
    M2 = h1_map(deck_top)             # p2*: H1(Y8) -> H1(Y4)
    M1 = h1_map(deck_bot)             # p1*: H1(Y4) -> H1(Y2)
    T2 = h1_map(deck_top, tau=True)   # tau2*: H1(Y4) -> H1(Y8)
    T1 = h1_map(deck_bot, tau=True)   # tau1*: H1(Y2) -> H1(Y4)
    ranks = {n: TS.gf2_rank([v2i(c) for c in M.T])
             for n, M in [("p2", M2), ("p1", M1), ("tau2", T2),
                          ("tau1", T1)]}
    assert ranks == {"p2": 4, "p1": 4, "tau2": 4, "tau1": 4}, ranks
    # exactness at both rungs
    K2 = TS._kernel_ints(M2)          # ker p2* in H1(Y8)
    K1 = TS._kernel_ints(M1)          # ker p1* in H1(Y4)
    imT2 = TS._colspace(T2)           # im tau2* in H1(Y8)
    imT1 = TS._colspace(T1)           # im tau1* in H1(Y4)
    imM2 = TS._colspace(M2)           # im p2* in H1(Y4)  (= SEAM)
    imM1 = TS._colspace(M1)           # im p1* in H1(Y2)  (= SEAM1)
    kerT2 = TS._kernel_ints(T2)       # ker tau2* in H1(Y4)
    kerT1 = TS._kernel_ints(T1)       # ker tau1* in H1(Y2)
    assert TS._span_eq(imT2, K2), "im tau2* != ker p2*"
    assert TS._span_eq(imT1, K1), "im tau1* != ker p1*"
    assert TS._span_eq(kerT2, imM2), "ker tau2* != im p2*"
    assert TS._span_eq(kerT1, imM1), "ker tau1* != im p1*"
    # deck-triviality (H4 re-measured directly): sigma* = id on H1(cover)
    sig2 = translation_mat(Y8, (0, 4))
    sig1 = translation_mat(Y4, (0, 2))
    assert (sig2 == np.eye(8, dtype=np.uint8)).all(), "sigma2* != id"
    assert (sig1 == np.eye(8, dtype=np.uint8)).all(), "sigma1* != id"
    print(f"[{time.monotonic()-t0:5.1f}s] Part 2: rank p* = rank tau* = 4 "
          f"= k/2 at BOTH rungs; exactness im tau* = ker p*, "
          f"ker tau* = im p*; sigma* = id (DeckTrivialOnH1, H4) both rungs")
    out["part2_ranks"] = ranks

    # SEAM = im p2*: 15 classes, orbits 12+3, class-stab orders
    sd = seam_data(Y4, M2)
    print(f"[{time.monotonic()-t0:5.1f}s]   SEAM = im p2*: dim 4, 15 "
          f"classes, orbits {[len(o) for o in sd['orbits']]}, class-stab "
          f"orders {sd['stab_orders']}  [A20 SS5 12+3 REPRODUCED]")
    out["seam"] = {"classes": sd["classes"],
                   "orbit_sizes": [len(o) for o in sd["orbits"]],
                   "reps": sd["reps"], "stab_orders": sd["stab_orders"]}

    # rung-1 lattice: SEAM1 (im p1*), its orbit split, and the descent
    # reachability datum W2 = p1*(SEAM)
    S1b, S1p = rref_ints(list(imM1))
    pts1 = TS._span_points(S1b, S1p) - {0}
    assert len(S1b) == 4 and len(pts1) == 15
    mats2 = TS._translation_action(Y2)
    orbs1 = sorted((sorted(o) for o in TS._orbits(pts1, mats2)),
                   key=len, reverse=True)
    W2 = sorted({TS._apply(M1, s)
                 for s in TS._span_points(*rref_ints(list(imM2)))} - {0})
    W2b, W2p = rref_ints(list(W2))
    # per-orbit-rep pushforward class (which Y2 class the shadows occupy)
    rep_push = {f"{rep:#x}": f"{TS._apply(M1, rep):#x}" for rep in sd["reps"]}
    # SEAM cap K1 (seam classes with [beta] = 0 shadows possible)
    Sbb, Sbp = sd["basis"], sd["piv"]
    K1b, K1p = rref_ints(list(K1))
    seam_cap_k1 = [s for s in sd["classes"] if TS.in_span(s, K1b, K1p)]
    print(f"[{time.monotonic()-t0:5.1f}s]   rung-1 lattice: SEAM1 = im p1* "
          f"dim 4, orbits {[len(o) for o in orbs1]}; "
          f"W2 = p1*(SEAM): dim {len(W2b)} ({len(W2)} classes); "
          f"rep pushforwards {rep_push}; |SEAM cap ker p1*| = "
          f"{len(seam_cap_k1)} classes")
    out["rung1_lattice"] = {
        "seam1_orbit_sizes": [len(o) for o in orbs1],
        "W2_dim": len(W2b), "W2_classes": [f"{w:#x}" for w in W2],
        "rep_pushforwards": rep_push,
        "seam_cap_K1": [f"{s:#x}" for s in seam_cap_k1],
        "K1_dim": len(K1),
    }

    # Bezout witnesses (banked, re-verified as convolution identities)
    from bb_lab.checks import circulant  # noqa: E402
    for cid, code, mhalf in [("bezout_y_18x8.json", Y8, 4),
                             ("bezout_y_18x4.json", Y4, 2)]:
        bez = json.loads((MAIN / "a20" / cid).read_text())
        vec = np.zeros(code.ng, dtype=np.uint8)
        for a, b in bez["P"]:
            vec[code.G.index((a, b))] ^= 1
        acc = (circulant(code.A).astype(np.uint8) @ vec) % 2
        vecq = np.zeros(code.ng, dtype=np.uint8)
        for a, b in bez.get("Q", []):
            vecq[code.G.index((a, b))] ^= 1
        acc = (acc + circulant(code.B).astype(np.uint8) @ vecq) % 2
        target = np.zeros(code.ng, dtype=np.uint8)
        target[code.G.index((0, 0))] = 1
        target[code.G.index((0, mhalf))] = 1
        assert (acc == target).all(), f"Bezout witness {cid} FAILS"
    print(f"[{time.monotonic()-t0:5.1f}s]   Bezout witnesses verified: "
          f"1+y^4 in (A8,B8), 1+y^2 in (A4,B4)  [(R) both rungs, "
          f"kernel-level]")
    out["bezout"] = "both witnesses re-verified as convolution identities"

    # seamC dictionary: A20's ker-d2 parametrization == SEAM via iota
    D2c, elts = chain_kernel_classes()
    iota4 = iota_perm(Y4)
    dict_map = {}
    sig_set = set()
    for mask, zeta in sorted(elts.items()):
        sc = seam_c(zeta)
        w_lab = sc[iota4]
        assert Y4.is_cycle(w_lab), f"iota(seamC {mask:#x}) not a cycle"
        assert not Y4.is_stab(w_lab), \
            f"iota(seamC {mask:#x}) is a stabilizer (delta2-inj fails?)"
        s = v2i(Y4.sig(w_lab))
        assert TS.in_span(s, Sbb, Sbp) and s != 0
        dict_map[mask] = s
        sig_set.add(s)
    assert sig_set == set(sd["classes"]), "seamC classes != SEAM classes"
    orb_of = {}
    for oi, o in enumerate(sd["orbits"]):
        for c in o:
            orb_of[c] = oi
    print(f"[{time.monotonic()-t0:5.1f}s]   seamC dictionary: the 15 "
          f"iota(seamC zeta) classes == SEAM \\ 0 EXACTLY; chain-mask "
          f"0x1 -> class {dict_map[1]:#x} (orbit "
          f"{['12','3'][orb_of[dict_map[1]]] if len(sd['orbits'][0])==12 else '?'}-sized"
          f" idx {orb_of[dict_map[1]]}), 0x3 -> class {dict_map[3]:#x} "
          f"(orbit idx {orb_of[dict_map[3]]})")
    out["seamC_dictionary"] = {
        "mask_to_class": {f"{m:#x}": f"{s:#x}" for m, s in dict_map.items()},
        "mask1_orbit_size": len(sd["orbits"][orb_of[dict_map[1]]]),
        "mask3_orbit_size": len(sd["orbits"][orb_of[dict_map[3]]]),
    }
    assert len(sd["orbits"][orb_of[dict_map[1]]]) == 12
    assert len(sd["orbits"][orb_of[dict_map[3]]]) == 3

    # X<->Z transpose duality spot-check at Y8 (antipode + block swap)
    iota8 = iota_perm(Y8)
    swap8 = np.concatenate([np.arange(Y8.ng, 2 * Y8.ng),
                            np.arange(0, Y8.ng)])
    dual = lambda v: v[iota8][swap8]  # noqa: E731
    for kv in Y8.kerHZ[:20]:
        assert not ((Y8.HX @ dual(kv)) % 2).any(), "duality: cycle fails"
    for row in Y8.HX[::12]:
        d = dual(row)
        assert TS.in_span(v2i(d), *rref_ints([v2i(r) for r in Y8.HZ])), \
            "duality: stab row fails"
    print(f"[{time.monotonic()-t0:5.1f}s]   X<->Z duality spot-check "
          f"(antipode + block swap) OK => Z-side floors follow from X-side")
    out["duality"] = "spot-checked at Y8"

    # comparison datum: the A32 (non-(R)) tower's top-rung tau_y* rank
    C = TS.BBCode("C", (30, 6), "x^9 + y + y^2", "y^3 + x^25 + x^26")
    BY = TS.BBCode("BY", (30, 3), "x^9 + y + y^2", "1 + x^25 + x^26")
    deck_y32 = TS.Deck(C, BY, lambda e: (e[0], e[1] % 3),
                       lambda e, s: (e[0], e[1] + 3 * s))
    Ty32 = h1_map(deck_y32, tau=True)   # H1(BY) -> H1(C), 12x8
    My32 = h1_map(deck_y32)             # H1(C) -> H1(BY), 8x12
    r_ty = TS.gf2_rank([v2i(c) for c in Ty32.T])
    imM32 = TS._colspace(My32)
    kerT32 = TS._kernel_ints(Ty32)
    ker_eq_im = TS._span_eq(kerT32, imM32)
    sig_y32 = translation_mat(C, (0, 3))
    deck_triv_32 = (sig_y32 == np.eye(12, dtype=np.uint8)).all()
    print(f"[{time.monotonic()-t0:5.1f}s]   A32-tower comparison (top "
          f"y-rung, (R) fails): rank tau_y* = {r_ty}, ker tau_y* = im "
          f"p_y*: {ker_eq_im}, sigma_y* = id: {deck_triv_32}")
    out["a32_comparison"] = {"rank_tau_top": int(r_ty),
                             "rank_p_top": 6,
                             "ker_tau_eq_im_p": bool(ker_eq_im),
                             "deck_trivial_top": bool(deck_triv_32)}

    out["wall_s"] = round(time.monotonic() - t0, 1)
    (DATA / "tower_cells.json").write_text(json.dumps(out, indent=1))
    print(f"total {out['wall_s']}s -> {DATA / 'tower_cells.json'}")


if __name__ == "__main__":
    main()
