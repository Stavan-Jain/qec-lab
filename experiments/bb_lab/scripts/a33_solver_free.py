"""A33 Part 5: the full solver-free upgrade — every remaining SAT-tier
input on the d(Y8) = 20 critical path re-derived at certificate tier.

  (c) d0 = 6 by direct BZ at n = 72: all 255 Y2 logical-class cosets
      censused to weight 6 (complete; parity kills odd weights) —
      floor AND witnesses in one pass; banked y72 witness in the set.
  (d) H3 = LogicalFloor 10 at Y4 by direct BZ at n = 144: all 255
      classes to weight 10 — no element below 10 (floor), the weight-10
      element census (witnesses; banked u10 + the H6 witness's sheet
      both present); seam classes do NOT attain (minima 14 / >= 20).
  (b) d1 = 10 re-derived by the TOWER (the a32_dby pattern one rung
      down, budget 8): a weight <= 8 Y4-logical u would have shadow
      beta = p1(u) with
        [beta] in SEAM1 \\ 0  -> beta in the 15-coset census <= 8,
                                 lift fiber cap (8-|beta|)/2: any lift
                                 is a counterexample (fibers EMPTY);
        beta = 0             -> u = tau1(gamma), |gamma| <= 4 < 6 = d0:
                                 dead by (c);
        beta in Stab \\ 0    -> non-stab lifts of the stab-orbit fiber
                                 sweep below; all weigh >= 10 (asserted).
  (a) H1 census-completeness re-derived by the tower (fiber union):
      every Y4-stabilizer b, |b| <= 18, has beta = p1(b) in Stab(Y2);
      beta != 0: b is a bounded-overflow lift of a censused Y2-stab
      orbit rep (caps (18-|beta|)/2, complete at class level since the
      fold G(Y4) ->> G(Y2) is onto); beta = 0: b = tau1(gamma) with
      [gamma] in ker tau1* = SEAM1 u {0}, |gamma| <= 8 (parity), i.e.
      tau1(Stab(Y2) <= 8) u tau1(seam1-census <= 8).  The class-key
      union must equal the banked 1,655 EXACTLY.

Output: data/a33/solver_free.json
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

from a30_coset_bz import (  # noqa: E402
    build_kernel, coset_base, disjoint_info_sets, run_window, unpack3,
)
from a30_rung_pass import i2v, v2i  # noqa: E402
import a32_tower_slice as TS  # noqa: E402
from a32_subclosures import enumerate_lifts  # noqa: E402
from a32_deep_fibers import enumerate_lifts_deep  # noqa: E402
from a32_sectorAC_full import batch_keys  # noqa: E402
from a33_tower_cells import build_tower, h1_map, rep_for, seam_data  # noqa: E402

MAIN = Path("/Users/stavanjain/Code/qec-lab/experiments/bb_lab/data")
DATA = LAB / "data" / "a33"


def census_multi(binp, tag, HX, offsets, W, rpair, n, deadline):
    """Multi-offset coset census (both windows, deduped per offset,
    empty-window bases handled).  offsets: list[(label, t0_vec)]."""
    I1, G1, I2, G2, kappa = disjoint_info_sets(HX)
    hits: dict[str, set[int]] = {lab: set() for lab, _ in offsets}
    nodes = 0
    empties = {}
    for wi, (window, Gs, r) in enumerate(
            [(I1, G1, rpair[0]), (I2, G2, rpair[1])]):
        bases = []
        for lab, tv in offsets:
            cb = coset_base(Gs, window, tv)
            wcb = int(cb.sum())
            empties.setdefault(lab, []).append(wcb)
            if 0 < wcb <= W:
                hits[lab].add(v2i(cb))
            bases.append(cb)
        res = run_window(binp, f"{tag}_w{wi}", Gs, bases, r, W, deadline)
        nodes += res["nodes"]
        for j, hx in res.pop("hit_rows"):
            v = unpack3(hx, n)
            if v.any():
                hits[offsets[j][0]].add(v2i(v))
    return hits, kappa, nodes, empties


def main():
    t0 = time.monotonic()
    out: dict = {}
    Y2, Y4, Y8, deck_top, deck_bot = build_tower()
    M2 = h1_map(deck_top)
    M1 = h1_map(deck_bot)
    sd = seam_data(Y4, M2)
    Sbb, Sbp = sd["basis"], sd["piv"]
    from a30_rung_pass import rref_ints  # noqa: E402
    S1b, S1p = rref_ints(list(TS._colspace(M1)))   # SEAM1 basis
    seam1_classes = sorted(TS._span_points(S1b, S1p) - {0})
    assert len(seam1_classes) == 15
    perms_y2 = TS._translation_perms(Y2)
    perms_y4 = TS._translation_perms(Y4)
    binp = build_kernel()
    deadline = time.monotonic() + 3600

    # ---------------------------------------------- (c) d0 = 6 direct BZ
    offs = [(f"L{s}", rep_for(Y2, s)) for s in range(1, 256)]
    hits, kappa2, nodes, _ = census_multi(
        binp, "a33_d0", Y2.HX, offs, 6, (3, 2), Y2.n, deadline)
    assert kappa2 == 32
    all6 = set()
    for lab, hs in hits.items():
        for h in hs:
            w = bin(h).count("1")
            assert w == 6, f"logical of weight {w} < 6 found?!"
            all6.add(h)
    wit6 = None
    for line in (MAIN / "a20" / "y72_ladder.log").read_text().splitlines():
        try:
            wit6 = json.loads(line)["witness"]
        except (json.JSONDecodeError, KeyError):
            continue
    u6 = np.zeros(Y2.n, dtype=np.uint8)
    u6[wit6] = 1
    assert v2i(u6) in all6, "banked y72 witness not in the census"
    cls6 = {lab for lab, hs in hits.items() if hs}
    seam1_attain = sum(1 for lab in cls6
                       if int(lab[1:]) in set(seam1_classes))
    print(f"[{time.monotonic()-t0:5.1f}s] (c) d0 = 6 DIRECT: 255-class "
          f"census to 6 — {len(all6)} weight-6 logicals in {len(cls6)} "
          f"classes (none lighter; parity kills 7); banked witness "
          f"present; seam1 classes attaining: {seam1_attain}  "
          f"[{nodes:.2e} nodes]")
    out["d0"] = {"weight6_vectors": len(all6), "classes_attaining":
                 len(cls6), "seam1_classes_attaining": seam1_attain,
                 "nodes": nodes}

    # ---------------------------------------------- (d) H3 direct BZ
    offs4 = [(f"L{s}", rep_for(Y4, s)) for s in range(1, 256)]
    hits4, kappa4, nodes4, _ = census_multi(
        binp, "a33_h3", Y4.HX, offs4, 10, (5, 4), Y4.n, deadline)
    assert kappa4 == 68
    all10 = set()
    for lab, hs in hits4.items():
        for h in hs:
            w = bin(h).count("1")
            assert w == 10, f"Y4 logical of weight {w} < 10?!"
            all10.add(h)
    wit10 = None
    for line in (MAIN / "a20" / "y144_ladder.log").read_text().splitlines():
        try:
            wit10 = json.loads(line)["witness"]
        except (json.JSONDecodeError, KeyError):
            continue
    u10 = np.zeros(Y4.n, dtype=np.uint8)
    u10[wit10] = 1
    assert v2i(u10) in all10, "banked u10 not in the census"
    w20 = np.load(MAIN / "a20" / "y8_weight20_witness.npy").astype(np.uint8)
    sheet_u = (deck_top.EMB[0].T @ w20) % 2
    assert v2i(sheet_u) in all10, "H6 witness sheet not in the census"
    cls10 = {int(lab[1:]) for lab, hs in hits4.items() if hs}
    seam_attain = sum(1 for s in cls10
                      if TS.in_span(s, Sbb, Sbp))
    assert seam_attain == 0, "a seam class attains d1 = 10?!"
    print(f"[{time.monotonic()-t0:5.1f}s] (d) H3 = LogicalFloor 10 "
          f"DIRECT: 255-class census to 10 — no element below 10; "
          f"{len(all10)} weight-10 logicals in {len(cls10)} classes "
          f"(banked u10 + H6-sheet present; 0 in SEAM: light logicals "
          f"concentrate outside im p2*)  [{nodes4:.2e} nodes]")
    out["h3"] = {"weight10_vectors": len(all10),
                 "classes_attaining": len(cls10),
                 "seam_classes_attaining": 0, "nodes": nodes4}

    # ------------------------------- Y2 censuses feeding (a) and (b)
    # stabilizers <= 18 (base 0)
    hits_s, _, nodes_s, _ = census_multi(
        binp, "a33_y2stab", Y2.HX,
        [("S", np.zeros(Y2.n, dtype=np.uint8))], 18, (9, 8), Y2.n,
        deadline)
    y2stabs = sorted(hits_s["S"])
    whs: dict[int, int] = {}
    for h in y2stabs:
        whs[bin(h).count("1")] = whs.get(bin(h).count("1"), 0) + 1
    mu2 = min(whs)
    print(f"[{time.monotonic()-t0:5.1f}s] Y2 stab census <= 18: "
          f"{len(y2stabs)} vectors {dict(sorted(whs.items()))} "
          f"(mu(Y2) = {mu2})  [{nodes_s:.2e} nodes]")
    out["y2_stabs"] = {"vectors": len(y2stabs),
                       "whist": {str(k): v for k, v in sorted(whs.items())},
                       "mu": mu2}

    # seam1-class cosets <= 8 (the gamma census for (a)+(b))
    offs1 = [(f"Z{s}", rep_for(Y2, s)) for s in seam1_classes]
    hits1, _, nodes1, _ = census_multi(
        binp, "a33_seam1", Y2.HX, offs1, 8, (4, 3), Y2.n, deadline)
    gammas = []
    for lab, hs in hits1.items():
        for h in hs:
            gammas.append(h)
    whg: dict[int, int] = {}
    for h in gammas:
        whg[bin(h).count("1")] = whg.get(bin(h).count("1"), 0) + 1
    print(f"[{time.monotonic()-t0:5.1f}s] seam1-coset census <= 8: "
          f"{len(gammas)} logicals {dict(sorted(whg.items()))} across "
          f"the 15 classes  [{nodes1:.2e} nodes]")
    out["seam1_census8"] = {"vectors": len(gammas),
                            "whist": {str(k): v
                                      for k, v in sorted(whg.items())}}

    # ---------------------------------------------- (b) d1 seam1 fibers
    cex = 0
    for h in gammas:
        g = i2v(h, Y2.n)
        cap = (8 - int(g.sum())) // 2
        lifts = enumerate_lifts(deck_bot, g, cap=cap)
        if lifts:
            cex += len(lifts)
            print(f"  !! weight <= 8 Y4-logical over seam1 shadow "
                  f"(|gamma| = {int(g.sum())}) — d1 < 10 ?!")
    assert cex == 0
    print(f"[{time.monotonic()-t0:5.1f}s] (b) d1 seam1 branch: all "
          f"{len(gammas)} lift fibers at cap (8-|g|)/2 EMPTY "
          f"[carry-infeasible] — no sub-10 logical with nonzero shadow "
          f"class")

    # ---------------------------------------------- (a) the fiber union
    # orbit compression: fibers over Y2-stab orbit reps (class-level
    # completeness: the fold G(Y4) ->> G(Y2) is onto)
    svecs = np.array([i2v(h, Y2.n) for h in y2stabs], dtype=np.uint8)
    keys = batch_keys(svecs, perms_y2)
    _, idx = np.unique(keys, axis=0, return_index=True)
    reps = svecs[idx]
    print(f"[{time.monotonic()-t0:5.1f}s] (a) fiber union: "
          f"{len(y2stabs)} Y2-stabs -> {len(reps)} orbit reps; "
          f"running caps (18-|beta|)/2 ...")
    union_keys: set[bytes] = set()
    n_lifts = 0
    n_stab_lifts = 0
    n_logical_lifts = 0
    min_logical_w = 99
    tF = time.monotonic()
    for ri, beta in enumerate(reps):
        wb = int(beta.sum())
        cap = (18 - wb) // 2
        lifts = (enumerate_lifts(deck_bot, beta, cap=cap) if cap <= 4
                 else enumerate_lifts_deep(deck_bot, beta, cap=cap))
        h = v2i(beta)
        for v0c, m2 in lifts.items():
            for v0i in (v0c, v0c ^ h):
                b = deck_bot.lift(i2v(v0i, Y2.n), beta)
                n_lifts += 1
                if Y4.is_stab(b):
                    n_stab_lifts += 1
                    union_keys.add(bytes(batch_keys(
                        b[None, :], perms_y4)[0]))
                else:
                    n_logical_lifts += 1
                    min_logical_w = min(min_logical_w, int(b.sum()))
        if ri % 500 == 0 and ri:
            print(f"    ... {ri}/{len(reps)} reps "
                  f"({time.monotonic()-t0:.0f}s)")
    assert min_logical_w >= 10, \
        f"non-stab lift of weight {min_logical_w} < 10 — d1 < 10?!"
    # beta = 0 family: tau1 of stabs <= 8 and of the seam1 census <= 8
    n_tau = 0
    for h in [x for x in y2stabs if bin(x).count("1") <= 8] + gammas:
        g = i2v(h, Y2.n)
        b = (deck_bot.TAU @ g) % 2
        assert Y4.is_stab(b), "tau1 of ker-tau1* cycle not a stabilizer"
        if 0 < int(b.sum()) <= 18:
            union_keys.add(bytes(batch_keys(b[None, :], perms_y4)[0]))
            n_tau += 1
    dtF = time.monotonic() - tF
    # compare with the banked class set
    banked = []
    for line in (MAIN / "a20" / "m_census_classes.jsonl").open():
        r = json.loads(line)
        if "b_support" in r:
            banked.append(r)
    bvecs = np.zeros((len(banked), Y4.n), dtype=np.uint8)
    for i, e in enumerate(banked):
        bvecs[i, e["b_support"]] = 1
    k_banked = {bytes(k) for k in batch_keys(bvecs, perms_y4)}
    assert union_keys == k_banked, (
        f"fiber union {len(union_keys)} classes != banked "
        f"{len(k_banked)}; diff {len(union_keys ^ k_banked)}")
    print(f"[{time.monotonic()-t0:5.1f}s] (a) fiber union == banked "
          f"1,655 classes EXACTLY ({dtF:.0f}s): {n_lifts} lifts "
          f"({n_stab_lifts} stab, {n_logical_lifts} logical, min "
          f"logical weight {min_logical_w} >= 10 => (b) stab branch "
          f"closed), tau1-family {n_tau}")
    out["h1_union"] = {
        "y2_orbit_reps": int(len(reps)), "lifts": n_lifts,
        "stab_lifts": n_stab_lifts, "logical_lifts": n_logical_lifts,
        "min_logical_lift_w": min_logical_w, "tau1_family": n_tau,
        "classes": len(union_keys), "equals_banked": True,
        "wall_s": round(dtF, 1),
    }
    print(f"[{time.monotonic()-t0:5.1f}s] (b) d1 = 10 TOWER-DERIVED: "
          f"seam1 fibers empty + stab-branch lifts >= 10 + beta=0 dead "
          f"by d0; witness = (d)'s weight-10 census")
    out["d1_tower"] = {"seam1_fibers_empty": True,
                       "stab_branch_min": min_logical_w,
                       "beta0": "dead by d0 = 6 (parity + tau)"}

    out["wall_s"] = round(time.monotonic() - t0, 1)
    (DATA / "solver_free.json").write_text(json.dumps(out, indent=1))
    print(f"total {out['wall_s']}s -> {DATA / 'solver_free.json'}")
    print("\nASSEMBLY (see A33 note): H1 (4 derivations) + H2 (1,655 "
          "deterministic rungs) + H3 (direct BZ + tower) + H4 (Bezout + "
          "sigma* = id) + H5 (direct BZ census + 1,680 rungs + descent) "
          "+ H6 (verified witness) => d(Y8) = 20, certificate tier, "
          "no SAT on the critical path.")


if __name__ == "__main__":
    main()
