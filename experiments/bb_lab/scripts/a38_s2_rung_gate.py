"""A38 S2 gate: the promoted rung engine + deep fiber lane vs EVERY banked
rung battery (charter §6.0 falsify-first — nothing new runs until this is
green).

What is promoted (bb_lab.tower): `enumerate_lifts_deep` (ordered-split
MITM, cap <= 8 — the a32_deep_fibers lane) and `RungCell` (the a30/a33
YRungCell architecture: dangerous rung with full sector dispatch + seam
feasibility rung; restricted lanes extended 7-8 via the ordered-split;
BZ lane via bb_lab.cosetbz, lazy).

Batteries (every banked number is a hard assert):

[A] a32 deep fiber lane
    A1  caps<=4 equivalence: enumerate_lifts_deep == enumerate_lifts on
        all 397 banked sector-C orbit fibers (bit-level, incl. m2 values)
    A2  the 28 banked deep fibers (|beta| in {6,10,12}, caps 8/6/5):
        per-fiber lift weight-hist == banked deep_fibers.jsonl EXACT
        (2,371 lifts total), then the banked trichotomy dispatch:
        stab lifts <= 20 in the banked A19 m24 census (by G-canonical
        key), light logicals == the banked A24 band-16 reachable set,
        flat-22 + logical rungs re-run through the LIBRARY RungCell on
        the BY->C top deck — verdict counts == banked (2,030 rungs PASS)

[B] a33 ibm288Y batteries (tower via library TowerCode/AxisDeck)
    B1  1,655 dangerous rungs at M = (21-w)//2 (banked a20 census, M
        asserted == banked m_req): per-row verdict+lane == banked
        h2_rungs.jsonl EXACT
    B2  1,680 seam rungs at M = (21-w)//2 over the banked seam census:
        per-row verdict+lane == banked seam_rungs.jsonl EXACT
    B3  hexagon rung at M = 7 (validate_sectors=True: V1 sector-linear
        trick 256/256): PASS, lane restricted<=6, 240 nontrivial sectors
        == banked rung_validation.json
    B4  the banked planted control (v = tau(u10) + HXc[0], |v| = 26,
        overflow 10 over the w6 stab shadow), TWICE:
        - M = 9: the NEW restricted<=8 ordered-split lane must PASS
          (banked BZ found min overflow exactly 10 => overflow <= 8
          empty — the deep join validated against the banked boundary)
        - M = 11: the library BZ lane must re-find it — VIOLATION with
          min overflow == 10, n_viol == 540, every weight >= 20
          (== banked rung_validation.json planted block)

[C] a36 bb288 direct-close reproduced END-TO-END through the library
    (TowerCode + AxisDeck + cosetbz + RungCell + batch_keys): 6-offset
    BZ pass 7,502,279,774 nodes exact; stab census 33,588 vectors /
    469 orbits, hists exact, orbit key set == banked file; seam census
    395 elements {12:3, 14:18, 16:374}, per-orbit {6,18,70,199,102};
    V1 4,096/4,096; 469/469 dangerous + 395/395 seam PASS with lane
    hists == banked direct_close_banked.json; witness re-verified
    end-to-end from banked support; covariance 3+3.

Output: data/a38/s2_rung_gate.json
"""

from __future__ import annotations

import gzip
import json
import math
import sys
import time
from pathlib import Path

import numpy as np

LAB = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(LAB / "src"))

from bb_lab import cosetbz  # noqa: E402
from bb_lab.tower import (  # noqa: E402
    AxisDeck, RungCell, TowerCode, batch_keys, enumerate_lifts,
    enumerate_lifts_deep, h1_map, i2v, orbits, perm_for, rep_for,
    rref_ints, span_points, translation_action, translation_perms, v2i,
)

MAIN = Path("/Users/stavanjain/Code/qec-lab/experiments/bb_lab/data")
DATA = LAB / "data"
OUT = DATA / "a38"


def battery_a32(out: dict, t0: float) -> None:
    GB = TowerCode("GB", (15, 3), "x^9 + y + y^2", "1 + x^10 + x^11")
    BY = TowerCode("BY", (30, 3), "x^9 + y + y^2", "1 + x^25 + x^26")
    C = TowerCode("C", (30, 6), "x^9 + y + y^2", "y^3 + x^25 + x^26")
    deck_x = AxisDeck(BY, GB, 0)
    deck_y = AxisDeck(C, BY, 1)
    perms_by = translation_perms(BY)

    stab_recs = [json.loads(x)
                 for x in (DATA / "a32" / "gb_census_stab.jsonl").open()]

    # ---- A1: caps<=4 equivalence on all 397 banked sector-C fibers
    reps: dict[int, np.ndarray] = {}
    for r in stab_recs:
        if r["w"] in (14, 16) and r["canon"] not in reps:
            v = np.zeros(GB.n, dtype=np.uint8)
            v[r["support"]] = 1
            reps[r["canon"]] = v
    assert len(reps) == 397
    n_eq = 0
    for beta in reps.values():
        cap = (22 - int(beta.sum())) // 2
        assert cap <= 4
        l1 = enumerate_lifts(deck_x, beta, cap)
        l2 = enumerate_lifts_deep(deck_x, beta, cap)
        assert l1 == l2, "deep != shallow at cap <= 4"
        n_eq += 1
    print(f"[{time.monotonic()-t0:6.1f}s] A1: deep == shallow enumerator "
          f"on {n_eq}/397 banked fibers (bit-level, caps 4/3)")
    out["A1"] = {"fibers": n_eq, "equal": True}

    # ---- A2: the 28 deep fibers + banked trichotomy dispatch
    banked_rows = [json.loads(x)
                   for x in (DATA / "a32" / "deep_fibers.jsonl").open()]
    assert len(banked_rows) == 28
    banked_by_w: dict[int, set] = {}
    bvs, bws = [], []
    for line in (MAIN / "a19" / "m24_census_classes.jsonl").open():
        r = json.loads(line)
        if "b_support" in r:
            b = np.zeros(BY.n, dtype=np.uint8)
            b[r["b_support"]] = 1
            bvs.append(b)
            bws.append(r["w"])
    BK = batch_keys(np.array(bvs, dtype=np.uint8), perms_by)
    for i, w in enumerate(bws):
        banked_by_w.setdefault(w, set()).add(bytes(BK[i]))
    reach16 = set()
    for x in (MAIN / "a24" /
              "cell_census_reach_band16_mchecks.jsonl").open():
        e = json.loads(x)
        b = np.zeros(BY.n, dtype=np.uint8)
        b[e["b_support"]] = 1
        reach16.add(bytes(batch_keys(b[None, :], perms_by)[0]))

    deep_reps: dict[int, dict] = {}
    for r in stab_recs:
        if r["w"] in (6, 10, 12):
            deep_reps.setdefault(r["canon"], r)
    assert len(deep_reps) == 28
    cell = RungCell("byc", BY, C, deck_y)
    assert len(cell.sector_basis) == 8
    stats = {"fibers": 0, "lifts": 0, "rungs": 0, "violations": 0,
             "flat22_rungs": 0, "banked_checked": 0,
             "light_logical_hits": 0, "verdicts": {}, "lift_hist": {}}
    rows_iter = iter(sorted(banked_rows, key=lambda r: (r["wbeta"],)))
    for rep_i, r in enumerate(sorted(deep_reps.values(),
                                     key=lambda r: (r["w"], r["canon"]))):
        beta = np.zeros(GB.n, dtype=np.uint8)
        beta[r["support"]] = 1
        wbeta = r["w"]
        cap = (22 - wbeta) // 2
        lifts = enumerate_lifts_deep(deck_x, beta, cap=cap)
        wbh: dict[str, int] = {}
        vh: dict[str, int] = {}
        for v0c, m2 in sorted(lifts.items()):
            b = deck_x.lift(i2v(v0c, GB.n), beta)
            wb = int(b.sum())
            assert wb == wbeta + 2 * m2 and BY.is_cycle(b)
            stats["lifts"] += 1
            k = f"{wb}"
            wbh[k] = wbh.get(k, 0) + 1
            stats["lift_hist"][k] = stats["lift_hist"].get(k, 0) + 1
            if BY.is_stab(b):
                if wb <= 20:
                    key = bytes(batch_keys(b[None, :], perms_by)[0])
                    assert key in banked_by_w.get(wb, set()), \
                        f"stab lift w{wb} missing from banked census!"
                    stats["banked_checked"] += 1
                    continue
                res = cell.rung(b, 1, time.monotonic() + 600)
                assert res["verdict"] != "ABORT", res
                stats["flat22_rungs"] += 1
                stats["rungs"] += 1
                vh[res["verdict"]] = vh.get(res["verdict"], 0) + 1
                stats["verdicts"][res["verdict"]] = \
                    stats["verdicts"].get(res["verdict"], 0) + 1
                if res["verdict"] == "VIOLATION":
                    stats["violations"] += 1
                continue
            if wb <= 16:
                stats["light_logical_hits"] += 1
                assert wb == 16, f"NEW reachable logical at w{wb}!"
                key = bytes(batch_keys(b[None, :], perms_by)[0])
                assert key in reach16, "NEW band-16 reachable class!"
                continue
            M = (24 - wb + 1) // 2
            res = cell.rung(b, M, time.monotonic() + 600)
            assert res["verdict"] != "ABORT", res
            stats["rungs"] += 1
            vh[res["verdict"]] = vh.get(res["verdict"], 0) + 1
            stats["verdicts"][res["verdict"]] = \
                stats["verdicts"].get(res["verdict"], 0) + 1
            if res["verdict"] == "VIOLATION":
                stats["violations"] += 1
        bk = next(rows_iter)
        assert bk["wbeta"] == wbeta and bk["cap"] == cap
        assert bk["lift_whist"] == wbh, \
            f"deep fiber {rep_i}: lift hist {wbh} != banked " \
            f"{bk['lift_whist']}"
        assert bk["verdicts"] == vh, \
            f"deep fiber {rep_i}: verdicts {vh} != banked {bk['verdicts']}"
        stats["fibers"] += 1
    assert stats["lifts"] == 2371, stats["lifts"]
    assert stats["violations"] == 0
    assert stats["rungs"] == 2030, stats["rungs"]
    print(f"[{time.monotonic()-t0:6.1f}s] A2: 28/28 deep fibers == banked "
          f"(2,371 lifts, per-fiber hists exact); trichotomy dispatch "
          f"through library RungCell: {stats['rungs']} rungs "
          f"{stats['verdicts']}, {stats['banked_checked']} stab lifts "
          f"banked-checked, {stats['light_logical_hits']} band-16 hits "
          f"reachable-known  == banked EXACT")
    out["A2"] = stats


def battery_a33(out: dict, t0: float) -> None:
    SPECS = {
        "Y2": ((18, 2), "1 + x + x^14*y", "1 + x + x^2*y"),
        "Y4": ((18, 4), "1 + x + x^14*y", "1 + x*y^2 + x^2*y^3"),
        "Y8": ((18, 8), "1 + x*y^4 + x^14*y", "1 + x*y^2 + x^2*y^7"),
    }
    codes = {n: TowerCode(n, *s) for n, s in SPECS.items()}
    Y4, Y8 = codes["Y4"], codes["Y8"]
    deck_top = AxisDeck(Y8, Y4, 1)
    cell = RungCell("top", Y4, Y8, deck_top)
    assert len(cell.sector_basis) == 8

    # ---- B1: 1,655 dangerous rungs vs banked h2_rungs.jsonl
    census = []
    for line in (MAIN / "a20" / "m_census_classes.jsonl").open():
        r = json.loads(line)
        if "b_support" in r:
            census.append(r)
    assert len(census) == 1655
    floors = {}
    for line in (MAIN / "a20" / "m_floors_results.jsonl").open():
        r = json.loads(line)
        floors[tuple(r["b_support"])] = r
    banked_h2 = [json.loads(x)
                 for x in (DATA / "a33" / "h2_rungs.jsonl").open()]
    assert len(banked_h2) == 1655
    verd: dict[str, int] = {}
    for i, e in enumerate(census):
        b = np.zeros(Y4.n, dtype=np.uint8)
        b[e["b_support"]] = 1
        M = (21 - e["w"]) // 2
        banked = floors[tuple(e["b_support"])]
        assert banked["m_req"] == M and banked["verdict"] == "UNSAT"
        r = cell.rung(b, M, time.monotonic() + 900)
        bk = banked_h2[i]
        assert bk["w"] == e["w"] and bk["M"] == M
        assert r["verdict"] == bk["verdict"] and r["lane"] == bk["lane"], \
            f"h2 row {i}: {r['verdict']}/{r['lane']} != banked " \
            f"{bk['verdict']}/{bk['lane']}"
        verd[r["verdict"]] = verd.get(r["verdict"], 0) + 1
    assert verd == {"PASS": 1655}, verd
    print(f"[{time.monotonic()-t0:6.1f}s] B1: 1,655/1,655 dangerous rungs "
          f"verdict+lane == banked h2_rungs.jsonl EXACT (M == banked "
          f"m_req 1,655/1,655)")
    out["B1"] = {"rungs": 1655, "exact_match": True}

    # ---- B2: 1,680 seam rungs vs banked seam_rungs.jsonl
    seam_rows = [json.loads(x)
                 for x in (DATA / "a33" / "seam_census.jsonl").open()]
    banked_sr = [json.loads(x)
                 for x in (DATA / "a33" / "seam_rungs.jsonl").open()]
    assert len(seam_rows) == 1680 and len(banked_sr) == 1680
    verd2: dict[str, int] = {}
    for i, e in enumerate(seam_rows):
        v = np.zeros(Y4.n, dtype=np.uint8)
        v[e["w_support"]] = 1
        M = (21 - e["w"]) // 2
        r = cell.seam_rung(v, M)
        bk = banked_sr[i]
        assert bk["w"] == e["w"] and bk["M"] == M
        assert r["verdict"] == bk["verdict"] and r["lane"] == bk["lane"], \
            f"seam row {i}: {r['verdict']}/{r['lane']} != banked"
        verd2[r["verdict"]] = verd2.get(r["verdict"], 0) + 1
    assert verd2 == {"PASS": 1680}, verd2
    print(f"[{time.monotonic()-t0:6.1f}s] B2: 1,680/1,680 seam rungs "
          f"verdict+lane == banked seam_rungs.jsonl EXACT")
    out["B2"] = {"rungs": 1680, "exact_match": True}

    # ---- B3: hexagon rung, V1 sector trick 256/256
    hexb = next(e for e in census if e["w"] == 6)
    bhex = np.zeros(Y4.n, dtype=np.uint8)
    bhex[hexb["b_support"]] = 1
    r = cell.rung(bhex, 7, time.monotonic() + 600, validate_sectors=True)
    assert r["verdict"] == "PASS" and r["lane"] == "restricted<=6"
    assert r["sectors_nontrivial"] == 240
    print(f"[{time.monotonic()-t0:6.1f}s] B3: V1 sector-linear 256/256 + "
          f"hexagon M=7 PASS lane restricted<=6, 240 nontrivial sectors "
          f"== banked")
    out["B3"] = r

    # ---- B4: the banked planted control, deep lane + BZ lane
    wit = None
    for line in (MAIN / "a20" / "y144_ladder.log").read_text().splitlines():
        try:
            wit = json.loads(line)["witness"]
        except (json.JSONDecodeError, KeyError):
            continue
    u10 = np.zeros(Y4.n, dtype=np.uint8)
    u10[wit] = 1
    assert Y4.is_cycle(u10) and not Y4.is_stab(u10) and u10.sum() == 10
    v = ((deck_top.TAU @ u10) + Y8.HX[0]) % 2
    assert Y8.is_cycle(v) and not Y8.is_stab(v)
    b6 = (deck_top.P @ v) % 2
    assert Y4.is_stab(b6) and int(b6.sum()) == 6
    _, ov, _ = deck_top.slice_data(v)
    wv = int(v.sum())
    assert wv == 6 + 2 * ov and (wv, ov) == (26, 10), (wv, ov)
    # (i) M = 9 -> restricted<=8 (the NEW ordered-split sizes).  Session
    # finding (mismatch resolved): the banked rung_validation.json
    # "found_min_overflow: 10" was computed over the TRUNCATED
    # violations[:5] list (arbitrary append order), not over all 540 —
    # the true fiber has overflow-8 solutions (weight-22 logicals,
    # consistent with d = 20).  Verified three ways below: the deep lane
    # finds them, a one-offset BZ walk over their own sector coset finds
    # them, and the full M = 11 BZ lane finds them; the banked walk's
    # completeness itself stands (n_viol 540 reproduced EXACTLY), and no
    # production verdict anywhere consumed a truncated-min datum.
    tD = time.monotonic()
    r9 = cell.rung(b6, 9, time.monotonic() + 3600, full_viols=True)
    dt9 = time.monotonic() - tD
    assert r9["verdict"] == "VIOLATION" \
        and r9["lane"] == "restricted<=8", r9
    assert r9["min_overflow"] == 8 and r9["ov_hist"] == {"8": 8}, r9
    assert all(x["weight"] >= 20 for x in r9["violations"])
    deep_v0s = {x["v0_hex"] for x in r9["violations"]}
    print(f"[{time.monotonic()-t0:6.1f}s] B4i: planted shadow at M=9 "
          f"(restricted<=8, NEW sizes 7-8): 8 distinct overflow-8 "
          f"solutions, all weight 22 >= 20 ({dt9:.1f}s)")
    # (ii) M = 11 -> BZ lane: must re-find the banked violation count
    # AND the deep lane's ov<=8 set exactly.
    tB = time.monotonic()
    r11 = cell.rung(b6, 11, time.monotonic() + 3600, full_viols=True)
    dtB = time.monotonic() - tB
    assert r11["verdict"] == "VIOLATION" and r11["lane"] == "bz", \
        {k: r11[k] for k in ("verdict", "lane")}
    assert r11["n_viol"] == 540, r11["n_viol"]
    assert r11["min_overflow"] == 8, r11["min_overflow"]
    assert all(x["weight"] >= 20 for x in r11["violations"])
    bz_v0s = {x["v0_hex"] for x in r11["violations"]
              if x["overflow"] <= 8}
    assert bz_v0s == deep_v0s, \
        f"lane disagreement: deep {len(deep_v0s)} vs bz {len(bz_v0s)}"
    print(f"[{time.monotonic()-t0:6.1f}s] B4ii: planted at M=11 (BZ "
          f"lane): n_viol 540 == banked EXACT; min overflow 8 (the "
          f"banked '10' was the truncated-[:5] reporting artifact); "
          f"ov<=8 distinct-v0 set == deep lane EXACTLY ({dtB:.1f}s)")
    out["B4"] = {"planted": {"pv_weight": wv, "overflow": ov},
                 "M9_deep": {"verdict": r9["verdict"], "lane": r9["lane"],
                             "ov_hist": r9["ov_hist"],
                             "secs": round(dt9, 1)},
                 "M11_bz": {"verdict": r11["verdict"],
                            "n_viol": r11["n_viol"],
                            "min_overflow": r11["min_overflow"],
                            "ov_hist": r11["ov_hist"],
                            "banked_min_overflow_10":
                                "truncation artifact (violations[:5])",
                            "secs": round(dtB, 1)},
                 "lane_agreement_ov_le_8": sorted(deep_v0s)}


def battery_a36(out: dict, t0: float) -> None:
    banked = json.loads(
        (DATA / "a36" / "direct_close_banked.json").read_text())
    SPECS = {
        "G8": ((12, 12), "x^3 + y^2 + y^7", "y^3 + x + x^2"),
        "GR": ((12, 6), "x^3 + y + y^2", "y^3 + x + x^2"),
    }
    G8 = TowerCode("G8", *SPECS["G8"])
    GR = TowerCode("GR", *SPECS["GR"])
    deck_y = AxisDeck(G8, GR, 1)
    My = h1_map(deck_y)
    from bb_lab.tower import colspace
    Sb = colspace(My)
    Sbb, Sbp = rref_ints(list(Sb))
    pts = span_points(Sbb) - {0}
    mats = translation_action(GR)
    orbs = sorted((sorted(o) for o in orbits(pts, mats)),
                  key=lambda o: (len(o), o[0]))
    reps_cls = [min(o) for o in orbs]
    assert [len(o) for o in orbs] == [3, 3, 9, 12, 36]
    rep_vecs = [rep_for(GR, c) for c in reps_cls]
    perms = translation_perms(GR)

    # ---- 6-offset BZ pass through bb_lab.cosetbz
    W, r1, r2 = 16, 8, 7
    binp = cosetbz.build_kernel()
    I1, G1, I2, G2, kappa = cosetbz.disjoint_info_sets(GR.HX)
    assert kappa == 66 and not (set(I1) & set(I2))
    offsets = [("S", np.zeros(GR.n, dtype=np.uint8))] + \
        [(f"R{i}", rv) for i, rv in enumerate(rep_vecs)]
    hits: dict[str, set[int]] = {lab: set() for lab, _ in offsets}
    nodes_total = 0
    for wi, (window, Gs, r) in enumerate([(I1, G1, r1), (I2, G2, r2)]):
        bases = []
        for lab, tv in offsets:
            cb = cosetbz.coset_base(Gs, window, tv)
            wcb = int(cb.sum())
            if 0 < wcb <= W:
                hits[lab].add(v2i(cb))
            bases.append(cb)
        res = cosetbz.run_window(binp, f"a38gate_w{wi}", Gs, bases, r, W,
                                 time.monotonic() + 3600)
        nodes_total += res["nodes"]
        for j, hx in res.pop("hit_rows"):
            v = cosetbz.unpack3(hx, GR.n)
            if v.any():
                hits[offsets[j][0]].add(v2i(v))
    assert nodes_total == banked["bz_pass"]["nodes"] == 7502279774
    print(f"[{time.monotonic()-t0:6.1f}s] C: BZ pass {nodes_total:,} nodes "
          f"== banked EXACT")

    # ---- stab census
    stab_vecs = np.array([i2v(h, GR.n) for h in sorted(hits["S"])],
                         dtype=np.uint8)
    ws = stab_vecs.sum(axis=1)
    whist = {str(int(w)): int((ws == w).sum()) for w in sorted(set(ws))}
    assert len(stab_vecs) == 33588
    assert whist == banked["stab_census"]["weight_hist"]
    for v in stab_vecs[:: max(1, len(stab_vecs) // 60)]:
        assert GR.is_stab(v)
    keys = batch_keys(stab_vecs, perms)
    orb_rep: dict[bytes, int] = {}
    for i, k in enumerate(keys):
        orb_rep.setdefault(bytes(k), i)
    assert len(orb_rep) == 469
    orb_whist: dict[str, int] = {}
    for i in orb_rep.values():
        orb_whist[str(int(ws[i]))] = orb_whist.get(str(int(ws[i])), 0) + 1
    assert orb_whist == banked["stab_census"]["orbit_weight_hist"]
    # orbit key set == banked orbit file
    bvecs = []
    for line in (DATA / "a36" / "stab_census_orbits_banked.jsonl").open():
        r = json.loads(line)
        b = np.zeros(GR.n, dtype=np.uint8)
        b[r["b_support"]] = 1
        bvecs.append(b)
    assert len(bvecs) == 469
    bkeys = {bytes(k) for k in
             batch_keys(np.array(bvecs, dtype=np.uint8), perms)}
    assert bkeys == set(orb_rep.keys()), "orbit key sets differ"
    print(f"[{time.monotonic()-t0:6.1f}s] C: stab census 33,588 vectors / "
          f"469 orbits, hists + orbit key set == banked EXACT")

    # ---- seam census
    seam_rows = []
    seam_whist: dict[str, int] = {}
    per_orbit = {}
    for oi, (lab, _) in enumerate(offsets[1:]):
        els = sorted(hits[lab])
        per_orbit[lab] = len(els)
        for h in els:
            v = i2v(h, GR.n)
            wv = int(v.sum())
            assert GR.is_cycle(v) and not GR.is_stab(v)
            assert v2i(GR.sig(v)) == reps_cls[oi]
            assert wv >= 12
            seam_rows.append((oi, v, wv))
            seam_whist[str(wv)] = seam_whist.get(str(wv), 0) + 1
    assert len(seam_rows) == 395
    assert seam_whist == banked["seam_census"]["weight_hist"]
    assert per_orbit == banked["seam_census"]["per_orbit"]
    print(f"[{time.monotonic()-t0:6.1f}s] C: seam census 395 elements "
          f"{seam_whist} per-orbit {per_orbit} == banked EXACT "
          f"(all >= 12 = d(gross))")

    # ---- rungs
    cell = RungCell("g8", GR, G8, deck_y)
    assert len(cell.sector_basis) == 12
    i6 = next(i for i in orb_rep.values() if ws[i] == 6)
    r6 = cell.rung(stab_vecs[i6], 6, time.monotonic() + 1200,
                   validate_sectors=True)
    assert r6["verdict"] == "PASS"
    print(f"[{time.monotonic()-t0:6.1f}s] C: V1 4,096/4,096 sector trick + "
          f"w6 rung PASS")
    verd: dict[str, int] = {}
    lanes: dict[str, int] = {}
    for key, i in orb_rep.items():
        b = stab_vecs[i]
        M = (18 - int(ws[i])) // 2
        r = (r6 if i == i6 else cell.rung(b, M, time.monotonic() + 3600))
        verd[r["verdict"]] = verd.get(r["verdict"], 0) + 1
        lanes[r["lane"]] = lanes.get(r["lane"], 0) + 1
    assert verd == {"PASS": 469}
    assert lanes == banked["dangerous_rungs"]["lanes"], lanes
    print(f"[{time.monotonic()-t0:6.1f}s] C: dangerous 469/469 PASS, "
          f"lanes {lanes} == banked EXACT")
    verd2: dict[str, int] = {}
    lanes2: dict[str, int] = {}
    for oi, v, wv in seam_rows:
        M = (18 - wv) // 2
        r = cell.seam_rung(v, M)
        verd2[r["verdict"]] = verd2.get(r["verdict"], 0) + 1
        lanes2[r["lane"]] = lanes2.get(r["lane"], 0) + 1
    assert verd2 == {"PASS": 395}
    assert lanes2 == banked["seam_rungs"]["lanes"], lanes2
    print(f"[{time.monotonic()-t0:6.1f}s] C: seam 395/395 PASS, lanes "
          f"{lanes2} == banked EXACT")

    # ---- covariance spot-checks
    g = (5, 1)
    perm_g = perm_for(GR, g)
    for i in list(orb_rep.values())[:3]:
        bt = stab_vecs[i][perm_g]
        rt = cell.rung(bt, (18 - int(ws[i])) // 2, time.monotonic() + 600)
        assert rt["verdict"] == "PASS"
    for oi, v, wv in seam_rows[:3]:
        rt = cell.seam_rung(v[perm_g], (18 - wv) // 2)
        assert rt["verdict"] == "PASS"

    # ---- witness end-to-end
    deck_x = AxisDeck(GR, TowerCode("B72", (6, 6),
                                    "x^3 + y + y^2", "y^3 + x + x^2"), 0)
    wit_row = json.loads(
        (DATA / "a36" / "w18_witness_banked.json").read_text())
    vwit = np.zeros(G8.n, dtype=np.uint8)
    vwit[wit_row["v_support"]] = 1
    assert G8.is_cycle(vwit) and not G8.is_stab(vwit)
    assert int(vwit.sum()) == 18
    bvec, m1, _ = deck_y.slice_data(vwit)
    assert int(bvec.sum()) == wit_row["shadow_w"] and m1 == wit_row["m1"]
    assert GR.is_cycle(bvec) and not GR.is_stab(bvec)
    from bb_lab.tower import in_span
    assert in_span(v2i(GR.sig(bvec)), Sbb, Sbp)
    beta, m2, _ = deck_x.slice_data(bvec)
    assert int(beta.sum()) == wit_row["beta_w"] and m2 == wit_row["m2"]
    print(f"[{time.monotonic()-t0:6.1f}s] C: witness re-verified "
          f"end-to-end (w18 over the w12 seam element, m1 = 3, "
          f"tau-diagonal profile) == banked; covariance 3+3 OK")
    out["C"] = {"bz_nodes": nodes_total, "stab": 33588, "orbits": 469,
                "seam": 395, "dangerous_pass": 469, "seam_pass": 395,
                "lanes_dangerous": lanes, "lanes_seam": lanes2,
                "witness": "re-verified", "exact_match": True}


def main() -> None:
    t0 = time.monotonic()
    out: dict = {}
    battery_a32(out, t0)
    battery_a33(out, t0)
    battery_a36(out, t0)
    out["wall_s"] = round(time.monotonic() - t0, 1)
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "s2_rung_gate.json").write_text(json.dumps(out, indent=1))
    print(f"GATE GREEN in {out['wall_s']}s -> "
          f"{OUT / 's2_rung_gate.json'}")


if __name__ == "__main__":
    main()
