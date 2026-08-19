"""A38 S1: the R4-regime methods closure — d(gross) = 12 re-derived
through the gross_xx tower's lower pair, completing A35's pair-regime
coverage 4/4 (A32 = R1, A33 = R2, A36 = R3, THIS = R4).

Tower (A35 screen row `gross_xx`): GR (12,6) [[144,12,12]] --x6-->
B72 (6,6) [[72,12,6]] --x3--> B36 (3,6) [[36,8,4]].  The pair lattice is
R4 = partial overlap (banked A35: dim S = dim K = 6, dim S^K = 2,
dim W = 4, K !<= S, preimage != S): no automatic pruning — the FULL
trisection is live, including sector A, which R3 (A36) killed outright.

Values here are NOT new: d(gross) = 12 is kernel-checked in QECLean
(the stronger tier), d(B72) = 6 is census-complete (A36 §4), d(B36) = 4
is full-enumeration exact (A35 [1b]).  The deliverable is the REGIME
EXECUTION: what the calculus's obligations look like when the pair
lattice offers no collapse, measured end-to-end, plus the R4-specific
descent lane (below).  Claim tier of everything computed here:
deterministic certificate (exact node counts, complete-by-construction
enumerations, in-line re-verification) — consistent with, and weaker
than, the Lean theorem it re-derives.

Parts:
  0  frames + decks through bb_lab.tower (Lemma-1 constructor asserts);
     the R4 pair lattice hard-asserted against the banked A35 numbers;
     SEAM = im p_x* (63 classes) split into B72-translation orbits;
     the CLASS-DETERMINED trisection of the seam species: [b] in S^K\\0
     (3 classes) => shadow beta is stab-or-0 (sectors B/C); [b] in
     S\\(S^K) (60 classes) => [beta] in W\\0 (sector A live).
  1  d(B72) = 6 census-complete re-derivation (the b = 0 branch input):
     full cycle census <= 8 over the dim-42 kernel; banked A36 numbers
     hard-asserted (1,110 cycles {6: 120, 8: 990}, 84 weight-6
     nontrivial logicals).
  2  the direct closure at n = 72: one multi-offset coset-BZ pass
     (kappa = 30, W = 10, r-pair (5,4), exact node asserts) delivering
     the stabilizer census <= 10 AND the seam orbit-rep coset censuses
     <= 10; dangerous rungs (V1 sector-scan validation, 2^12 dispatch)
     + seam rungs at M = (12 - |b|)/2; covariance spot-checks;
     G-transport via fold surjectivity.
  3  the witness (d <= 12): tau_x(gamma) for a weight-6 B72 logical
     gamma with tau_x*[gamma] != 0, verified end-to-end (weight 12,
     nontrivial gross logical, shadow 0 — the diagonal sector, as A33's
     H6).
  4  the R4 DESCENT LANE (the regime deliverable): both n = 72 census
     species re-derived from n = 36 data via the full trisection —
     sector C fibers over B36 stab orbit reps, sector B tau_b-family
     over the light cycle census, and the R4-specific SECTOR A fibers
     over W-coset element orbit reps, with the preimage LEAKAGE
     measured (lifts over W-shadows whose class falls outside SEAM —
     the quantitative form of "no automatic pruning": in R1/A32 the
     preimage identity made reachability decidable one rung below; here
     it provably cannot).  Base-exactness is NOT assumed on the bottom
     rung ((R) fails there): the tau_b family is classified by ker
     tau_b* directly, and the ker tau_b* vs im p_b* gap is measured.
     Gate: G-canonical key-set EQUALITY with the direct censuses.

Output: data/a38/r4_close.json (+ census jsonl files)
Run:    cd experiments/bb_lab && uv run python scripts/a38_r4_close.py
"""

from __future__ import annotations

import json
import math
import sys
import time
from pathlib import Path

import numpy as np

LAB = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(LAB / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from bb_lab import tower as tw  # noqa: E402
from bb_lab.cosetbz import (  # noqa: E402
    build_kernel, coset_base, disjoint_info_sets, rref, run_window, unpack3,
)

from a33_rung_cell import YRungCell  # noqa: E402

DATA = LAB / "data" / "a38"
DATA.mkdir(parents=True, exist_ok=True)

TARGET = 12
W = TARGET - 2          # parity: every cycle weight is even

T0 = time.monotonic()


def log(msg: str) -> None:
    print(f"[{time.monotonic()-T0:6.1f}s] {msg}", flush=True)


def main() -> None:
    out: dict = {}

    # ------------------------------------------------------------- Part 0
    GR = tw.TowerCode("GR", (12, 6), "x^3 + y + y^2", "y^3 + x + x^2")
    B72 = tw.TowerCode("B72", (6, 6), "x^3 + y + y^2", "y^3 + x + x^2")
    B36 = tw.TowerCode("B36", (3, 6), "1 + y + y^2", "y^3 + x + x^2")
    deck_x = tw.AxisDeck(GR, B72, 0)     # untwisted, (R)
    deck_b = tw.AxisDeck(B72, B36, 0)    # twisted, non-(R)
    assert (GR.k, B72.k, B36.k) == (12, 12, 8)
    assert not deck_x.twisted() and deck_b.twisted()
    for code in (GR, B72, B36):          # parity exhaustive
        assert not any(int(kv.sum()) % 2 for kv in code.kerHZ)
    log("Part 0: codes built k = 12/12/8; Lemma-1 asserts PASS both "
        "rungs; parity exhaustive (W = 10)")

    Mx = tw.h1_map(deck_x)               # p_x*: H1(GR) -> H1(B72)
    Tx = tw.h1_map(deck_x, tau=True)
    Mb = tw.h1_map(deck_b)               # p_b*: H1(B72) -> H1(B36)
    Tb = tw.h1_map(deck_b, tau=True)
    S = tw.colspace(Mx)                  # SEAM = im p_x*  (dim 6)
    K = tw.kernel_ints(Mb)               # ker p_b*        (dim 6)
    dS, dK = tw.gf2_rank(S), tw.gf2_rank(K)
    dSK = dS + dK - tw.gf2_rank(list(S) + list(K))
    Sb, Sp = tw.rref_ints(list(S))
    Wimg = [tw.apply_map(Mb, s) for s in tw.span_points(Sb)]
    dW = tw.gf2_rank(Wimg)
    K_in_S = all(tw.in_span(x, Sb, Sp) for x in K)
    # the banked A35 R4 lattice, hard-asserted
    assert (dS, dK, dSK, dW, K_in_S) == (6, 6, 2, 4, False), \
        (dS, dK, dSK, dW, K_in_S)
    # (R) rung lattice on the x-rung; non-(R) on the b-rung
    assert tw.span_eq(tw.colspace(Tx), tw.kernel_ints(Mx))
    assert tw.span_eq(tw.kernel_ints(Tx), S)          # base-exact (R)
    assert not tw.span_eq(tw.kernel_ints(Tb), tw.colspace(Mb)), \
        "bottom rung base-exactness should FAIL (non-(R))"
    log(f"  R4 pair lattice REPRODUCED: dim S = {dS}, K = {dK}, "
        f"S^K = {dSK}, W = {dW}, K<=S: {K_in_S} — full trisection live")

    # SEAM classes and their B72-translation orbit split
    seam_classes = tw.span_points(Sb) - {0}
    assert len(seam_classes) == 63
    mats72 = tw.translation_action(B72)
    seam_orbs = sorted((sorted(o) for o in tw.orbits(seam_classes, mats72)),
                       key=lambda o: (len(o), o[0]))
    seam_reps = [min(o) for o in seam_orbs]
    # class-determined trisection: W-image zero <=> class in S^K
    SKb, SKp = tw.rref_ints([x for x in Wimg if False] or
                            [x for x in K if tw.in_span(x, Sb, Sp)])
    n_bc = sum(1 for c in seam_classes if tw.apply_map(Mb, c) == 0)
    assert n_bc == 3, f"S^K nonzero classes {n_bc} != 3 (dim 2)"
    orb_secA = [i for i, r in enumerate(seam_reps)
                if tw.apply_map(Mb, r) != 0]
    log(f"  SEAM: 63 classes in {len(seam_orbs)} B72-orbits, sizes "
        f"{[len(o) for o in seam_orbs]}; class-determined trisection: "
        f"3 classes (S^K) -> sectors B/C, 60 -> sector A")
    out["part0"] = {"k": [12, 12, 8],
                    "pair": {"dim_S": dS, "dim_K": dK, "dim_SK": dSK,
                             "dim_W": dW, "K_in_S": K_in_S},
                    "seam_orbit_sizes": [len(o) for o in seam_orbs],
                    "seam_classes_SK": n_bc, "seam_classes_A": 60}

    # ------------------------------------------------------------- Part 1
    # d(B72) = 6 census-complete (single-window BZ over the cycle space)
    binp = build_kernel()
    Kk = np.array(B72.kerHZ, dtype=np.uint8)
    assert Kk.shape[0] == 42
    Rn, piv = rref(Kk)
    assert len(piv) == 42
    res = run_window(binp, "a38_b72cyc", Rn,
                     [np.zeros(B72.n, dtype=np.uint8)], 8, 8,
                     time.monotonic() + 600)
    assert res["nodes"] == sum(math.comb(42, s) for s in range(1, 9))
    gammas = []
    for _, hx in res.pop("hit_rows"):
        g = unpack3(hx, B72.n)
        if g.any():
            gammas.append(g)
    gws = np.array([int(g.sum()) for g in gammas])
    ghist = {int(w): int((gws == w).sum()) for w in sorted(set(gws))}
    ntl = [(int(w), g) for w, g in zip(gws, gammas) if not B72.is_stab(g)]
    d72 = min(w for w, _ in ntl)
    n_d72 = sum(1 for w, _ in ntl if w == d72)
    # banked A36 §4 numbers, hard-asserted
    assert len(gammas) == 1110 and ghist == {6: 120, 8: 990}, ghist
    assert d72 == 6 and n_d72 == 84, (d72, n_d72)
    log(f"Part 1: bb72 cycle census <= 8 == banked (1,110 cycles "
        f"{ghist}); d(B72) = 6 CENSUS-COMPLETE ({n_d72} weight-6 "
        f"logicals) — the b = 0 branch input, solver-free")
    out["part1"] = {"cycles": len(gammas), "weight_hist": ghist,
                    "d_B72": d72, "n_min": n_d72,
                    "nodes": res["nodes"]}

    # ------------------------------------------------------------- Part 2
    # direct closure at n = 72: one multi-offset BZ pass
    I1, G1, I2, G2, kappa = disjoint_info_sets(B72.HX)
    assert kappa == 30
    rep_vecs = [tw.rep_for(B72, c) for c in seam_reps]
    offsets = [("S", np.zeros(B72.n, dtype=np.uint8))] + \
        [(f"R{i}", rv) for i, rv in enumerate(rep_vecs)]
    r1, r2 = 5, 4
    hits: dict[str, set[int]] = {lab: set() for lab, _ in offsets}
    nodes_total = 0
    for wi, (window, Gs, r) in enumerate([(I1, G1, r1), (I2, G2, r2)]):
        bases = []
        for lab, tv in offsets:
            cb = coset_base(Gs, window, tv)
            wcb = int(cb.sum())
            if 0 < wcb <= W:          # the empty-window coset-base edge
                hits[lab].add(tw.v2i(cb))
            bases.append(cb)
        res = run_window(binp, f"a38_r4_w{wi}", Gs, bases, r, W,
                         time.monotonic() + 600)
        nodes_total += res["nodes"]
        for j, hx in res.pop("hit_rows"):
            v = unpack3(hx, B72.n)
            if v.any():
                hits[offsets[j][0]].add(tw.v2i(v))
    exp = sum(math.comb(30, s) for s in range(1, r1 + 1)) + \
        sum(math.comb(30, s) for s in range(1, r2 + 1))
    assert nodes_total == exp
    log(f"Part 2: 6-window BZ pass done, nodes {nodes_total:.3e} "
        f"(exact assert OK), {len(offsets)} offsets")
    out["bz_pass"] = {"kappa": kappa, "W": W, "r_pair": [r1, r2],
                      "nodes": nodes_total, "offsets": len(offsets)}

    # stab census <= 10
    perms72 = tw.translation_perms(B72)
    stab_vecs = np.array([tw.i2v(h, B72.n) for h in sorted(hits["S"])],
                         dtype=np.uint8)
    ws = stab_vecs.sum(axis=1)
    assert (ws % 2 == 0).all() and (ws <= W).all()
    for v in stab_vecs[:: max(1, len(stab_vecs) // 40)]:
        assert B72.is_stab(v)
    whist = {int(w): int((ws == w).sum()) for w in sorted(set(ws))}
    # cross-assert vs the banked A36 <=16 census truncated to <=10
    assert whist == {6: 36, 10: 216}, whist
    keys = tw.batch_keys(stab_vecs, perms72)
    orb_rep: dict[bytes, int] = {}
    for i, kk in enumerate(keys):
        orb_rep.setdefault(bytes(kk), i)
    orb_whist: dict[int, int] = {}
    for i in orb_rep.values():
        orb_whist[int(ws[i])] = orb_whist.get(int(ws[i]), 0) + 1
    log(f"  stab census <= 10: {len(stab_vecs)} vectors {whist} == "
        f"banked truncation; {len(orb_rep)} orbits {orb_whist}")
    out["stab_census"] = {"vectors": int(len(stab_vecs)),
                          "weight_hist": whist, "orbits": len(orb_rep),
                          "orbit_weight_hist": orb_whist}

    # seam coset censuses <= 10
    seam_rows = []
    seam_whist: dict[int, int] = {}
    for oi, (lab, _) in enumerate(offsets[1:]):
        els = sorted(hits[lab])
        for h in els:
            v = tw.i2v(h, B72.n)
            wv = int(v.sum())
            assert B72.is_cycle(v) and not B72.is_stab(v)
            assert tw.v2i(B72.sig(v)) == seam_reps[oi], "wrong class"
            assert wv >= d72, \
                f"seam element weight {wv} < d(B72) = 6 ?!"
            seam_rows.append((oi, v, wv))
            seam_whist[wv] = seam_whist.get(wv, 0) + 1
    log(f"  seam censuses <= 10: {len(seam_rows)} elements "
        f"{dict(sorted(seam_whist.items()))} across "
        f"{len(seam_reps)} orbit-rep cosets  [minima >= 6 = d(B72)]")
    out["seam_census"] = {"elements": len(seam_rows),
                          "weight_hist": dict(sorted(seam_whist.items())),
                          "per_orbit": {f"R{i}": len(hits[f"R{i}"])
                                        for i in range(len(seam_reps))}}
    with (DATA / "r4_stab_census_orbits.jsonl").open("w") as f:
        for key, i in orb_rep.items():
            f.write(json.dumps({
                "w": int(ws[i]),
                "b_support": sorted(int(j) for j in
                                    np.nonzero(stab_vecs[i])[0])}) + "\n")
    with (DATA / "r4_seam_census.jsonl").open("w") as f:
        for oi, v, wv in seam_rows:
            f.write(json.dumps({
                "orbit": oi, "class": f"{seam_reps[oi]:#x}", "w": wv,
                "w_support": sorted(int(j) for j in
                                    np.nonzero(v)[0])}) + "\n")

    # rung engine (a33 YRungCell, generic over the deck)
    cell = YRungCell("r4", B72, GR, deck_x)
    assert len(cell.sector_basis) == 12
    i6 = next(i for i in orb_rep.values() if ws[i] == 6)
    r6 = cell.rung(stab_vecs[i6], (TARGET - 6) // 2,
                   time.monotonic() + 1200, validate_sectors=True)
    assert r6["verdict"] == "PASS", r6
    log(f"  V1: 4096/4096 sector linear trick OK + weight-6 rung M=3 "
        f"PASS (lane {r6['lane']})")

    verd: dict[str, int] = {}
    lanes: dict[str, int] = {}
    tR = time.monotonic()
    for key, i in orb_rep.items():
        b = stab_vecs[i]
        M = (TARGET - int(ws[i])) // 2
        r = (r6 if i == i6 else cell.rung(b, M, time.monotonic() + 600))
        verd[r["verdict"]] = verd.get(r["verdict"], 0) + 1
        lanes[r["lane"]] = lanes.get(r["lane"], 0) + 1
        assert r["verdict"] == "PASS", r
    log(f"  dangerous rungs: {len(orb_rep)}/{len(orb_rep)} PASS "
        f"({time.monotonic()-tR:.1f}s, lanes {lanes})")
    out["dangerous_rungs"] = {"rungs": len(orb_rep), "verdicts": verd,
                              "lanes": lanes}

    verd2: dict[str, int] = {}
    lanes2: dict[str, int] = {}
    tS = time.monotonic()
    for oi, v, wv in seam_rows:
        M = (TARGET - wv) // 2
        r = cell.seam_rung(v, M)
        verd2[r["verdict"]] = verd2.get(r["verdict"], 0) + 1
        lanes2[r["lane"]] = lanes2.get(r["lane"], 0) + 1
        assert r["verdict"] == "PASS", r
    log(f"  seam rungs: {len(seam_rows)}/{len(seam_rows)} PASS "
        f"({time.monotonic()-tS:.1f}s, lanes {lanes2})")
    out["seam_rungs"] = {"rungs": len(seam_rows), "verdicts": verd2,
                         "lanes": lanes2}

    # covariance spot-checks + G-transport
    g = (2, 1)
    perm_g = tw.perm_for(B72, g)
    for i in list(orb_rep.values())[:3]:
        rt = cell.rung(stab_vecs[i][perm_g],
                       (TARGET - int(ws[i])) // 2, time.monotonic() + 600)
        assert rt["verdict"] == "PASS"
    for oi, v, wv in seam_rows[:3]:
        rt = cell.seam_rung(v[perm_g], (TARGET - wv) // 2)
        assert rt["verdict"] == "PASS"
    # fold surjectivity on groups: G(GR) ->> G(B72)
    img = {deck_x.fold(e) for e in GR.G}
    assert len(img) == B72.G.cardinality
    log("  covariance 3+3 OK; fold G(GR) ->> G(B72) onto — verdicts "
        "transport to every translate")

    # ------------------------------------------------------------- Part 3
    # witness: tau_x(gamma), gamma a weight-6 B72 logical with
    # tau_x*[gamma] != 0  (tau-lift nontrivial <=> class outside ker
    # tau_x* = im p_x* — base-exactness HOLDS on the (R) x-rung)
    n_in = n_out = 0
    wit = None
    imb, imp = tw.rref_ints(list(tw.colspace(Mx)))
    for wv, g in ntl:
        if wv != 6:
            continue
        inside = tw.in_span(tw.v2i(B72.sig(g)), imb, imp)
        n_in += int(inside)
        n_out += int(not inside)
        if not inside and wit is None:
            wit = g
    assert wit is not None
    v_wit = (deck_x.TAU @ wit) % 2
    assert int(v_wit.sum()) == TARGET
    assert GR.is_cycle(v_wit) and not GR.is_stab(v_wit)
    bchk, m1, _ = deck_x.slice_data(v_wit)
    assert not bchk.any() and m1 == 6      # diagonal sector, |v| = 2|u|
    log(f"Part 3: WITNESS tau_x(gamma) verified end-to-end: weight 12 "
        f"nontrivial gross logical, shadow 0 (diagonal sector); the 84 "
        f"weight-6 logicals split {n_in} inside / {n_out} outside "
        f"im p_x*")
    out["witness"] = {"weight": 12, "sector": "diagonal",
                      "gamma_w": 6, "m1": int(m1),
                      "w6_classes_in_imp": n_in,
                      "w6_classes_out_imp": n_out,
                      "v_support": sorted(int(j) for j in
                                          np.nonzero(v_wit)[0])}

    # ------------------------------------------------------------- Part 4
    # the R4 descent lane: both n = 72 species from n = 36 data
    log("Part 4: R4 descent lane (full trisection, sector A live)")
    perms36 = tw.translation_perms(B36)

    def keyset(vlist):
        if not vlist:
            return set()
        return {bytes(kk) for kk in
                tw.batch_keys(np.array(vlist, dtype=np.uint8), perms72)}

    k_stab_direct = keyset(list(stab_vecs))
    k_seam_direct = keyset([v for _, v, _ in seam_rows])

    # --- B36 ground sets by full enumeration (stabs 2^14, cycles 2^22)
    stab36_ints = list(tw.span_points(B36.rsHX_b))
    assert len(stab36_ints) == 1 << 14
    stab36 = [tw.i2v(x, B36.n) for x in stab36_ints
              if 0 < bin(x).count("1") <= W]
    sw36 = {}
    for v in stab36:
        sw36[int(v.sum())] = sw36.get(int(v.sum()), 0) + 1
    assert sw36 == {6: 18, 8: 45, 10: 108}, sw36   # A35 demo truncation
    s36keys = tw.batch_keys(np.array(stab36, dtype=np.uint8), perms36)
    s36_reps: dict[bytes, int] = {}
    for i, kk in enumerate(s36keys):
        s36_reps.setdefault(bytes(kk), i)
    # W-space membership test on H1(B36)
    Wb2, Wp2 = tw.rref_ints([x for x in Wimg if x])
    Kk36 = np.array(B36.kerHZ, dtype=np.uint8)
    dimk36 = Kk36.shape[0]
    assert dimk36 == 22
    wcos36 = []          # nontrivial logicals, class in W \ 0, w <= 10
    light36 = []         # nonzero cycles of weight <= 4 (tau_b family)
    for lo in range(0, 1 << dimk36, 1 << 18):
        idx = np.arange(lo, lo + (1 << 18), dtype=np.int64)
        bits = ((idx[:, None] >> np.arange(dimk36)) & 1).astype(np.uint8)
        V = (bits @ Kk36) % 2
        wts = V.sum(axis=1)
        sigs = (V @ B36.zreps.T) % 2
        svals = sigs @ (1 << np.arange(B36.k, dtype=np.int64))
        sel_w = (wts <= W) & (svals != 0)
        for v, s in zip(V[sel_w], svals[sel_w].tolist()):
            if tw.in_span(int(s), Wb2, Wp2):
                wcos36.append(v.copy())
        sel_l = (wts > 0) & (wts <= 4)
        for v in V[sel_l]:
            light36.append(v.copy())
    ww36 = {}
    for v in wcos36:
        ww36[int(v.sum())] = ww36.get(int(v.sum()), 0) + 1
    lw36 = {}
    for v in light36:
        lw36[int(v.sum())] = lw36.get(int(v.sum()), 0) + 1
    w36keys = tw.batch_keys(np.array(wcos36, dtype=np.uint8), perms36)
    w36_reps: dict[bytes, int] = {}
    for i, kk in enumerate(w36keys):
        w36_reps.setdefault(bytes(kk), i)
    log(f"  B36 ground sets: stabs <= 10 {sw36} ({len(s36_reps)} "
        f"orbits); W-coset elements <= 10 {dict(sorted(ww36.items()))} "
        f"({len(w36_reps)} orbits); light cycles <= 4 "
        f"{dict(sorted(lw36.items()))}")
    out["b36_ground"] = {"stab_whist": sw36, "stab_orbits": len(s36_reps),
                         "wcoset_whist": dict(sorted(ww36.items())),
                         "wcoset_orbits": len(w36_reps),
                         "light_cycles": dict(sorted(lw36.items()))}

    # --- sector B: tau_b family over the light cycles (beta = 0).
    # (R) FAILS on the bottom rung, so base-exactness may not classify
    # tau-lifts; classify by ker tau_b* directly and MEASURE the gap.
    kerTb_b, kerTb_p = tw.rref_ints(tw.kernel_ints(Tb))
    imPb_b, imPb_p = tw.rref_ints(list(tw.colspace(Mb)))
    tau_stab, tau_seam, tau_other = [], [], []
    n_exact_agree = n_exact_disagree = 0
    for g in light36:
        b = (deck_b.TAU @ g) % 2
        assert int(b.sum()) == 2 * int(g.sum())
        if int(b.sum()) > W:
            continue
        sg = tw.v2i(B36.sig(g))
        pred_ker = tw.in_span(sg, kerTb_b, kerTb_p)   # tau-lift trivial?
        pred_imp = tw.in_span(sg, imPb_b, imPb_p)
        is_st = B72.is_stab(b)
        assert is_st == pred_ker, "tau* kernel misclassified a lift"
        n_exact_agree += int(pred_ker == pred_imp)
        n_exact_disagree += int(pred_ker != pred_imp)
        if is_st:
            tau_stab.append(b)
        elif tw.in_span(tw.v2i(B72.sig(b)), Sb, Sp):
            tau_seam.append(b)
        else:
            tau_other.append(b)
    log(f"  sector B (tau_b family): {len(tau_stab)} stab / "
        f"{len(tau_seam)} seam / {len(tau_other)} other-class lifts; "
        f"ker tau_b* vs im p_b* class agreement {n_exact_agree}/"
        f"{n_exact_agree + n_exact_disagree} (base-exactness FAILS on "
        f"this rung — the disagreements are the non-(R) signature)")
    out["sectorB"] = {"tau_stab": len(tau_stab),
                      "tau_seam": len(tau_seam),
                      "tau_other": len(tau_other),
                      "ker_vs_imp_agree": n_exact_agree,
                      "ker_vs_imp_disagree": n_exact_disagree}

    # --- sectors C and A: bounded-overflow fibers over B36 orbit reps
    def run_fibers(rep_dict, vecs, sector_name):
        fib_stab, fib_seam, fib_other = [], [], []
        n_fib = n_lift = n_empty = 0
        for kk, i in rep_dict.items():
            beta = vecs[i]
            wb = int(beta.sum())
            cap = (W - wb) // 2
            if cap < 0:
                continue
            lifts = tw.enumerate_lifts(deck_b, beta, cap)
            n_fib += 1
            if not lifts:
                n_empty += 1
            for v0_int, m2 in lifts.items():
                b = deck_b.lift(tw.i2v(v0_int, B36.n), beta)
                assert int(b.sum()) == wb + 2 * m2 <= W
                n_lift += 1
                if B72.is_stab(b):
                    fib_stab.append(b)
                elif tw.in_span(tw.v2i(B72.sig(b)), Sb, Sp):
                    fib_seam.append(b)
                else:
                    fib_other.append(b)
        log(f"  sector {sector_name}: {n_fib} fibers, {n_lift} lifts "
            f"({len(fib_stab)} stab / {len(fib_seam)} seam / "
            f"{len(fib_other)} other-class = leakage), {n_empty} empty "
            f"({(n_empty/n_fib if n_fib else 0):.0%})")
        return fib_stab, fib_seam, fib_other, n_fib, n_lift, n_empty

    c_stab, c_seam, c_other, cf, cl, ce = run_fibers(
        s36_reps, stab36, "C (stab shadows)")
    a_stab, a_seam, a_other, af, al, ae = run_fibers(
        w36_reps, wcos36, "A (W-coset shadows)")
    assert not a_stab, "a W-shadow lift cannot be a B72 stabilizer"
    out["sectorC"] = {"fibers": cf, "lifts": cl, "empty": ce,
                      "stab": len(c_stab), "seam": len(c_seam),
                      "other": len(c_other)}
    out["sectorA"] = {"fibers": af, "lifts": al, "empty": ae,
                      "seam": len(a_seam), "other_leakage": len(a_other)}

    # --- key-set equality with the direct censuses
    k_stab_desc = keyset(tau_stab) | keyset(c_stab)
    k_seam_desc = keyset(tau_seam) | keyset(c_seam) | keyset(a_seam)
    assert k_stab_desc == k_stab_direct, \
        (len(k_stab_desc), len(k_stab_direct))
    assert k_seam_desc == k_seam_direct, \
        (len(k_seam_desc), len(k_seam_direct))
    log(f"  KEY-SET EQUALITY: stab {len(k_stab_desc)} == direct; seam "
        f"{len(k_seam_desc)} == direct  [both n = 72 species re-derived "
        f"from n = 36 data through the FULL R4 trisection]")
    out["equality"] = {"stab_keys": len(k_stab_desc),
                       "seam_keys": len(k_seam_desc),
                       "stab_equal": True, "seam_equal": True}

    # ------------------------------------------------------ the assembly
    log(f"""
ASSEMBLY — d([[144,12,12]]) = 12 re-derived at certificate tier:
  |v| <= 10 (parity), b = p_x(v) for nontrivial X-logical v of GR:
  1. b = 0:      v = tau_x(u), [u] != 0, |v| = 2|u| >= 2 d(B72) = 12
                 [d(B72) = 6 census-complete, Part 1 — boundary case of
                 the G5 ceiling: 12 = 2*6 exactly]
  2. b stab:     census <= 10 complete ({out['stab_census']['vectors']}
                 vectors, {out['stab_census']['orbits']} orbits) +
                 {out['dangerous_rungs']['rungs']} rungs ALL PASS +
                 G-transport
  3. [b] != 0:   [b] in SEAM (63 classes, {len(seam_orbs)} orbits);
                 coset censuses <= 10 + {out['seam_rungs']['rungs']}
                 seam rungs ALL PASS + G-transport
  => d >= 12;  the tau-diagonal witness gives d <= 12.  d = 12.
  Z side by BB transpose duality.  CONSISTENT with the kernel-checked
  QECLean theorem (the stronger tier); R4 regime coverage COMPLETE.
""")
    out["wall_s"] = round(time.monotonic() - T0, 1)
    (DATA / "r4_close.json").write_text(json.dumps(out, indent=1))
    log(f"total {out['wall_s']}s -> {DATA / 'r4_close.json'}")


if __name__ == "__main__":
    main()
