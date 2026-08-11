"""A33 Part 4: H5 closed DIRECT — the theorem-completing step.

One 3-offset coset-BZ pass at n = 144 (kappa = 68, W = 19 -> 18 by
parity, r-pair (9, 8), complete to r1+r2+1 = 18):

  offset S   t0 = 0                : the Y4 stabilizer census <= 18
             -> H1 DIRECTLY re-derived at n = 144 (fourth independent
                derivation: SAT, V7-analytic, [tower fibers in Part 5],
                and now direct BZ) — G-canonical class set must equal
                the banked 1,655 EXACTLY
  offset R12 t0 = rep(class 0x6)   : the 12-orbit seam coset census
             -> every element w, |w| <= 18; expect 1,680 elements
                (280 class-stab-orbit classes x 6), weights 14/16/18
  offset R3  t0 = rep(class 0x73)  : the 3-orbit seam coset census
             -> expected EMPTY (v7_seam engine measured 0; the naive
                SeamCosetFloor 20 is TRUE on this orbit)

Then the lift-aware floor: for EVERY censused seam element w, the seam
rung at M = (20 - |w|)/2 — PASS = every cover cycle over w (all of which
are automatically nontrivial logicals, stabilizer transport) has
|v| = |w| + 2*overflow >= 20.  All 15 classes follow by G-transport
(covariance spot-checked).

Banked cross-checks (falsify-first):
  - the 278 safe_m_floors_0x1 elements (chain frame -> iota) must all
    appear in the census and their banked lift-UNSAT == my rung PASS
  - the 280 v7_seam_0x1 classes (lab frame) == my census's G-canonical
    key set exactly
  - the seam_floor SAT@18 witness appears in the census
  - the empty-window coset-base weights are checked (the latent
    coset-base edge of the census() helper — elements with empty window
    restriction are handled explicitly here)

Tightness probe: escalate M on one w18 and one w14 element until the
rung FINDS a logical (every find must weigh >= 20 — else d < 20).

Output: data/a33/h5_direct.json + seam_census.jsonl + seam_rungs.jsonl
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

from a30_coset_bz import (  # noqa: E402
    build_kernel, coset_base, disjoint_info_sets, run_window, unpack3,
)
from a30_rung_pass import i2v, rref_ints, v2i  # noqa: E402
import a32_tower_slice as TS  # noqa: E402
from a32_sectorAC_full import batch_keys  # noqa: E402
from a33_tower_cells import (  # noqa: E402
    build_tower, chain_kernel_classes, h1_map, iota_perm, rep_for, seam_c,
    seam_data,
)
from a33_rung_cell import YRungCell  # noqa: E402

MAIN = Path("/Users/stavanjain/Code/qec-lab/experiments/bb_lab/data")
DATA = LAB / "data" / "a33"


def main():
    t0 = time.monotonic()
    out: dict = {}
    Y2, Y4, Y8, deck_top, deck_bot = build_tower()
    M2 = h1_map(deck_top)
    sd = seam_data(Y4, M2)
    rep12_cls, rep3_cls = sd["reps"]
    assert len(sd["orbits"][0]) == 12 and len(sd["orbits"][1]) == 3
    rep12 = rep_for(Y4, rep12_cls)
    rep3 = rep_for(Y4, rep3_cls)
    perms_y4 = TS._translation_perms(Y4)
    print(f"[{time.monotonic()-t0:5.1f}s] seam reps: 12-orbit class "
          f"{rep12_cls:#x}, 3-orbit class {rep3_cls:#x}")

    # ------------------------------------------------- the 3-offset BZ pass
    binp = build_kernel()
    I1, G1, I2, G2, kappa = disjoint_info_sets(Y4.HX)
    assert kappa == 68 and not (set(I1) & set(I2))
    offsets = [("S", np.zeros(Y4.n, dtype=np.uint8)),
               ("R12", rep12), ("R3", rep3)]
    W, r1, r2 = 18, 9, 8
    hits: dict[str, set[int]] = {lab: set() for lab, _ in offsets}
    nodes_total = 0
    # empty-window elements (the coset-base edge): handled explicitly
    empty_info = []
    for wi, (window, Gs, r) in enumerate([(I1, G1, r1), (I2, G2, r2)]):
        bases = []
        for lab, tv in offsets:
            cb = coset_base(Gs, window, tv)
            wcb = int(cb.sum())
            empty_info.append({"window": wi, "offset": lab, "w": wcb})
            if 0 < wcb <= W:
                hits[lab].add(v2i(cb))
            bases.append(cb)
        res = run_window(binp, f"a33_h5_w{wi}", Gs, bases, r, W,
                         time.monotonic() + 3600)
        nodes_total += res["nodes"]
        for j, hx in res.pop("hit_rows"):
            v = unpack3(hx, Y4.n)
            if v.any():
                hits[offsets[j][0]].add(v2i(v))
        print(f"[{time.monotonic()-t0:5.1f}s]   window {wi}: r={r}, "
              f"nodes {res['nodes']:.3e} (exact-count assert OK), "
              f"wall {res['wall_s']}s")
    exp_nodes = sum(math.comb(68, s) for s in range(1, r1 + 1)) + \
        sum(math.comb(68, s) for s in range(1, r2 + 1))
    assert nodes_total == exp_nodes
    out["bz_pass"] = {"kappa": kappa, "W": W, "r_pair": [r1, r2],
                      "nodes": nodes_total,
                      "empty_window_bases": empty_info}
    print(f"[{time.monotonic()-t0:5.1f}s] BZ pass complete: "
          f"{nodes_total:.3e} nodes; empty-window base weights "
          f"{[e['w'] for e in empty_info]} (all > {W} except handled)")

    # ------------------------------------------------- offset S: H1 direct
    stab_vecs = np.array([i2v(h, Y4.n) for h in sorted(hits["S"])],
                         dtype=np.uint8)
    ws = stab_vecs.sum(axis=1)
    assert (ws % 2 == 0).all() and (ws <= 18).all()
    for v in stab_vecs[:: max(1, len(stab_vecs) // 40)]:
        assert Y4.is_stab(v)
    k_direct = {bytes(k) for k in batch_keys(stab_vecs, perms_y4)}
    census = []
    for line in (MAIN / "a20" / "m_census_classes.jsonl").open():
        r = json.loads(line)
        if "b_support" in r:
            census.append(r)
    cvecs = np.zeros((len(census), Y4.n), dtype=np.uint8)
    for i, e in enumerate(census):
        cvecs[i, e["b_support"]] = 1
    k_banked = {bytes(k) for k in batch_keys(cvecs, perms_y4)}
    assert k_direct == k_banked and len(k_direct) == 1655, \
        f"H1 direct census: {len(k_direct)} classes vs banked 1655"
    whist = {}
    for w in ws:
        whist[int(w)] = whist.get(int(w), 0) + 1
    print(f"[{time.monotonic()-t0:5.1f}s] offset S: {len(stab_vecs)} "
          f"stabilizer vectors <= 18, weight hist {whist}; G-canonical "
          f"classes = 1,655 == banked EXACTLY  [H1 re-derived DIRECT]")
    out["h1_direct"] = {"vectors": int(len(stab_vecs)),
                        "weight_hist": whist, "classes": 1655,
                        "equals_banked": True}

    # ------------------------------------------------- offset R3: EMPTY
    assert not hits["R3"], \
        f"3-orbit coset has {len(hits['R3'])} elements <= 18 (expected 0)"
    print(f"[{time.monotonic()-t0:5.1f}s] offset R3: coset EMPTY at "
          f"<= 18 => naive SeamCosetFloor 20 HOLDS on the 3-orbit "
          f"(H5 vacuous there; matches v7_seam engine + explains the "
          f"banked 0x3 SAT BUDGET)")
    out["r3"] = {"elements": 0, "verdict": "coset minimum >= 20 outright"}

    # ------------------------------------------------- offset R12: census
    seam_els = sorted(hits["R12"])
    svecs = np.array([i2v(h, Y4.n) for h in seam_els], dtype=np.uint8)
    sws = svecs.sum(axis=1)
    swh = {}
    for w in sws:
        swh[int(w)] = swh.get(int(w), 0) + 1
    assert len(seam_els) == 1680, f"{len(seam_els)} elements != 1680"
    assert swh == {14: 6, 16: 84, 18: 1590}, swh
    Sbb, Sbp = sd["basis"], sd["piv"]
    for v in svecs:
        assert Y4.is_cycle(v) and not Y4.is_stab(v)
        assert v2i(Y4.sig(v)) == rep12_cls, "census element in wrong class"
    keys12 = batch_keys(svecs, perms_y4)
    kset12 = {bytes(k) for k in keys12}
    assert len(kset12) == 280, f"{len(kset12)} G-canonical classes != 280"
    print(f"[{time.monotonic()-t0:5.1f}s] offset R12: 1,680 elements "
          f"(weights {swh}), all verified class-{rep12_cls:#x} logicals; "
          f"280 G-canonical classes  [v7_seam counts REPRODUCED]")
    out["r12_census"] = {"elements": 1680, "weight_hist": swh,
                         "g_classes": 280}
    with (DATA / "seam_census.jsonl").open("w") as f:
        for v, w in zip(svecs, sws):
            f.write(json.dumps({
                "w": int(w),
                "w_support": sorted(int(j) for j in np.nonzero(v)[0]),
            }) + "\n")

    # banked cross-checks
    iota4 = iota_perm(Y4)
    D2c, elts = chain_kernel_classes()
    sc1 = seam_c(elts[1])
    # (i) safe_m_floors_0x1: chain frame; assert coset membership there,
    # then iota-map and require census membership
    D2cols = [v2i(D2c[:, j]) for j in range(D2c.shape[1])]
    D2b, D2p = rref_ints(D2cols)
    from a30_rung_pass import reduce_int  # noqa: E402
    banked278 = []
    for line in (MAIN / "a20" / "safe_m_floors_0x1.jsonl").open():
        r = json.loads(line)
        wv = np.zeros(Y4.n, dtype=np.uint8)
        wv[r["w_support"]] = 1
        assert reduce_int(v2i((wv + sc1) % 2), D2b, D2p) == 0, \
            "banked element not in the chain-frame 0x1 coset"
        banked278.append((wv, r["lift"]))
    assert len(banked278) == 278
    bvecs = np.array([wv[iota4] for wv, _ in banked278], dtype=np.uint8)
    for v in bvecs:
        assert Y4.is_cycle(v) and not Y4.is_stab(v)
        assert TS.in_span(v2i(Y4.sig(v)), Sbb, Sbp)
    kb = {bytes(k) for k in batch_keys(bvecs, perms_y4)}
    assert len(kb) == 278 and kb <= kset12, \
        "banked 278 not a subset of the census classes"
    assert all(lift == "UNSAT" for _, lift in banked278)
    missing = len(kset12 - kb)
    # (ii) v7_seam_0x1 (lab frame): exact class-set equality
    v7rows = [json.loads(x) for x in
              (MAIN / "a20" / "v7_seam_0x1.jsonl").open()]
    v7v = np.zeros((len(v7rows), Y4.n), dtype=np.uint8)
    for i, r in enumerate(v7rows):
        v7v[i, r["w_support"]] = 1
    for v in v7v:
        assert Y4.is_cycle(v) and not Y4.is_stab(v)
        assert TS.in_span(v2i(Y4.sig(v)), Sbb, Sbp)
    kv7 = {bytes(k) for k in batch_keys(v7v, perms_y4)}
    assert kv7 == kset12, "v7_seam classes != BZ census classes"
    # (iii) the seam_floor SAT@18 witness
    sf = [json.loads(x) for x in
          (MAIN / "a20" / "seam_floor_results.jsonl").open()]
    wit18 = np.zeros(Y4.n, dtype=np.uint8)
    wit18[[r for r in sf if r["class_mask"] == 1][0]["witness"]] = 1
    kwit = bytes(batch_keys(wit18[None, iota4], perms_y4)[0])
    assert kwit in kset12
    print(f"[{time.monotonic()-t0:5.1f}s] banked cross-checks: 278/278 "
          f"safe_m_floors elements in census (dying census was "
          f"{missing} classes short of 280); v7_seam 280 classes == "
          f"census EXACTLY; SAT@18 witness present")
    out["banked_crosschecks"] = {
        "safe_m_floors_278": "all present, all lift-UNSAT",
        "census_classes_missing_from_banked": missing,
        "v7_seam_280": "exact class-set equality",
        "sat18_witness": "present",
    }

    # ------------------------------------------------- the seam rungs
    cell_top = YRungCell("top", Y4, Y8, deck_top)
    verd = {}
    lane_hist = {}
    rows = []
    tR = time.monotonic()
    banked_keys = kb
    agree = 0
    for v, w in zip(svecs, sws):
        M = (21 - int(w)) // 2
        r = cell_top.seam_rung(v, M)
        verd[r["verdict"]] = verd.get(r["verdict"], 0) + 1
        lane_hist[r["lane"]] = lane_hist.get(r["lane"], 0) + 1
        if r["verdict"] != "PASS":
            print("  !! SEAM VIOLATION:", json.dumps(r)[:300])
        key = bytes(batch_keys(v[None, :], perms_y4)[0])
        if key in banked_keys and r["verdict"] == "PASS":
            agree += 1
        rows.append({"w": int(w), "M": M, "verdict": r["verdict"],
                     "lane": r["lane"]})
    dt = time.monotonic() - tR
    assert verd == {"PASS": 1680}, verd
    print(f"[{time.monotonic()-t0:5.1f}s] seam rungs: 1,680/1,680 PASS "
          f"({dt:.1f}s, lanes {lane_hist}); banked-agreement elements "
          f"{agree} (all 6 translates of each banked class PASS)")
    out["seam_rungs"] = {"rungs": 1680, "verdicts": verd,
                         "lanes": lane_hist, "wall_s": round(dt, 1),
                         "banked_agreeing_elements": agree}
    with (DATA / "seam_rungs.jsonl").open("w") as f:
        for r in rows:
            f.write(json.dumps(r) + "\n")

    # covariance spot-check: verdicts transport along G
    g = (3, 1)
    perm_g = TS._perm_for(Y4, g)
    for v, w in list(zip(svecs, sws))[:3]:
        vt = v[perm_g]
        M = (21 - int(w)) // 2
        rt = cell_top.seam_rung(vt, M)
        assert rt["verdict"] == "PASS"
    print(f"[{time.monotonic()-t0:5.1f}s] covariance spot-check: "
          f"translated elements give identical verdicts (3/3)")

    # ------------------------------------------------- tightness probes
    probes = []
    for wtarget in (18, 14):
        idx = int(np.nonzero(sws == wtarget)[0][0])
        v = svecs[idx]
        Mstart = (21 - wtarget) // 2
        hit = None
        for M in range(Mstart + 1, 8):
            r = cell_top.seam_rung(v, M)
            if r["verdict"] == "VIOLATION":
                ovs = [x["overflow"] for x in r["violations"]]
                wts = [x["weight"] for x in r["violations"]]
                assert min(wts) >= 20, "sub-20 logical in probe?!"
                hit = {"w": wtarget, "first_hit_M": M,
                       "min_overflow": min(ovs),
                       "lightest_logical_over_w": min(wts)}
                break
        if hit is None:
            hit = {"w": wtarget, "first_hit_M": None,
                   "note": f"no logical over w at overflow <= 6"}
        probes.append(hit)
        print(f"[{time.monotonic()-t0:5.1f}s] tightness probe w{wtarget}: "
              f"{hit}")
    out["tightness_probes"] = probes

    out["wall_s"] = round(time.monotonic() - t0, 1)
    (DATA / "h5_direct.json").write_text(json.dumps(out, indent=1))
    print(f"total {out['wall_s']}s -> {DATA / 'h5_direct.json'}")
    print("H5 CLOSED (direct): 12-orbit = 1,680/1,680 element rungs PASS; "
          "3-orbit = coset empty (floor >= 20 outright); all 15 classes "
          "by G-transport.")


if __name__ == "__main__":
    main()
