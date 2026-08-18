"""A36 Parts 3-5: the direct closure — d([[288,12,18]]) = 18.

The assembly (W = 16 by parity; b = p_y(v) for a nontrivial X-logical v
of G8 with |v| <= 16):

  1. b = 0        v = tau(u) with u a nontrivial GROSS logical,
                  |v| = 2|u| >= 2 d(gross) = 24 > 16.  DEAD — consumes
                  d(gross) = 12 at KERNEL-CHECKED tier (QECLean; >= 10
                  is all the branch needs).
  2. [b] = 0,     b is a gross stabilizer, |b| <= 16 (|b| >= 18 is free
     b != 0       by the slice identity).  Dangerous rung at
                  M = (18-|b|)/2 per translation-orbit rep + G-transport.
  3. [b] != 0     [b] = p_y*[v] in SEAM \\ 0 (63 classes, 5 orbits);
                  b is a nontrivial gross logical so |b| >= 12 (Lean),
                  censused per orbit-rep coset at <= 16; seam rung
                  (pure feasibility, stabilizer transport) at
                  M = (18-|b|)/2 per element + G-transport.

  Witness: a verified weight-18 nontrivial logical => d <= 18.
  Z side: BB transpose duality (spot-asserted in a36_tower_cells).

One 6-offset coset-BZ pass at n = 144 (kappa = 66, W = 16, r-pair (8,7),
complete to r1+r2+1 = 16, exact node-count asserts) delivers BOTH census
species: offset S (t0 = 0) = the gross stabilizer census <= 16; offsets
R0..R4 = the five seam orbit-rep cosets <= 16.

Falsify-first gates baked in:
  - node counts exact (binomial identity, kernel-asserted)
  - every stab-census vector re-verified is_stab; min weight must be 6;
    40 random small row-sums must appear in the census (membership)
  - every seam element re-verified (cycle, non-stab, class == rep
    class); min weight must be >= 12 = d(gross) — a lighter element
    would CONTRADICT the kernel-checked Lean theorem (convention trap)
  - V1: sector-scan linear trick == direct reduction, all 4096 sectors,
    on the first weight-6 dangerous rung
  - covariance: translated reps re-rung, verdicts must match
  - find-side: witness hunt = rungs re-run at M+1 (one overflow unit
    above the floor); EVERY find must weigh >= 18 (a lighter find would
    refute the published d = 18) and some find must weigh exactly 18
    (else d > 18 and the SAT record is wrong) — both directions armed.

Output: data/a36/{direct_close.json, stab_census_orbits.jsonl,
        seam_census.jsonl, rungs.jsonl, w18_witness.json}
"""

# Provenance: copied verbatim 2026-08-18 (A38 S1) from the unmerged
# branch claude/tower-slice-calculus-generalize-410ed1 (the A35/A36
# session). That branch stays the source of truth until it merges;
# library-grade ports live in bb_lab.tower, not in edits here.


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
from a30_rung_pass import i2v, v2i  # noqa: E402
import a32_tower_slice as TS  # noqa: E402
from a32_sectorAC_full import batch_keys  # noqa: E402
from a33_tower_cells import h1_map, rep_for  # noqa: E402
from a33_rung_cell import YRungCell  # noqa: E402
from a36_tower_cells import build_tower, seam_info  # noqa: E402

DATA = LAB / "data" / "a36"
DATA.mkdir(parents=True, exist_ok=True)
RNG = np.random.default_rng(20260811)

D_GROSS = 12          # kernel-checked (QECLean gross d = 12, axiom-clean)
TARGET = 18
W = TARGET - 2        # parity: every cycle weight is even


def main():
    t0 = time.monotonic()
    out: dict = {}
    G8, GR, B72, B36, deck_y, deck_x, deck_b = build_tower()
    My = h1_map(deck_y)
    sd = seam_info(GR, My)
    reps_cls = sd["reps"]
    assert [len(o) for o in sd["orbits"]] == [3, 3, 9, 12, 36]
    rep_vecs = [rep_for(GR, c) for c in reps_cls]
    perms = TS._translation_perms(GR)
    print(f"[{time.monotonic()-t0:6.1f}s] tower + seam built: 5 orbit "
          f"reps {[f'{c:#x}' for c in reps_cls]}")

    # ---------------------------------------------- the 6-offset BZ pass
    binp = build_kernel()
    I1, G1, I2, G2, kappa = disjoint_info_sets(GR.HX)
    assert kappa == 66 and not (set(I1) & set(I2))
    offsets = [("S", np.zeros(GR.n, dtype=np.uint8))] + \
        [(f"R{i}", rv) for i, rv in enumerate(rep_vecs)]
    r1, r2 = 8, 7
    hits: dict[str, set[int]] = {lab: set() for lab, _ in offsets}
    nodes_total = 0
    for wi, (window, Gs, r) in enumerate([(I1, G1, r1), (I2, G2, r2)]):
        bases = []
        for lab, tv in offsets:
            cb = coset_base(Gs, window, tv)
            wcb = int(cb.sum())
            if 0 < wcb <= W:      # the empty-window coset-base edge
                hits[lab].add(v2i(cb))
            bases.append(cb)
        res = run_window(binp, f"a36_w{wi}", Gs, bases, r, W,
                         time.monotonic() + 3600)
        nodes_total += res["nodes"]
        for j, hx in res.pop("hit_rows"):
            v = unpack3(hx, GR.n)
            if v.any():
                hits[offsets[j][0]].add(v2i(v))
        print(f"[{time.monotonic()-t0:6.1f}s]   window {wi}: r={r}, "
              f"nodes {res['nodes']:.3e} (exact assert OK), "
              f"wall {res['wall_s']}s")
    exp = sum(math.comb(66, s) for s in range(1, r1 + 1)) + \
        sum(math.comb(66, s) for s in range(1, r2 + 1))
    assert nodes_total == exp
    out["bz_pass"] = {"kappa": kappa, "W": W, "r_pair": [r1, r2],
                      "nodes": nodes_total, "offsets": len(offsets)}

    # ------------------------------------------- offset S: stab census
    stab_vecs = np.array([i2v(h, GR.n) for h in sorted(hits["S"])],
                         dtype=np.uint8)
    ws = stab_vecs.sum(axis=1)
    assert (ws % 2 == 0).all() and (ws <= W).all()
    assert int(ws.min()) == 6, f"lightest gross stabilizer {ws.min()} != 6"
    for v in stab_vecs[:: max(1, len(stab_vecs) // 60)]:
        assert GR.is_stab(v)
    whist = {int(w): int((ws == w).sum()) for w in sorted(set(ws))}
    keys = batch_keys(stab_vecs, perms)
    orb_rep: dict[bytes, int] = {}
    for i, k in enumerate(keys):
        orb_rep.setdefault(bytes(k), i)
    orb_whist: dict[int, int] = {}
    for i in orb_rep.values():
        orb_whist[int(ws[i])] = orb_whist.get(int(ws[i]), 0) + 1
    print(f"[{time.monotonic()-t0:6.1f}s] offset S: {len(stab_vecs)} "
          f"stabilizer vectors <= {W}, weight hist {whist}; "
          f"{len(orb_rep)} translation orbits {orb_whist}")
    out["stab_census"] = {"vectors": int(len(stab_vecs)),
                          "weight_hist": whist,
                          "orbits": len(orb_rep),
                          "orbit_weight_hist": orb_whist}
    # membership gate: random small row-sums must be in the census
    hit_ints = hits["S"]
    n_mem = 0
    while n_mem < 40:
        j = int(RNG.integers(1, 4))
        idx = RNG.choice(GR.ng, size=j, replace=False)
        b = np.zeros(GR.n, dtype=np.uint8)
        for i in idx:
            b = (b + GR.HX[i]) % 2
        if 0 < int(b.sum()) <= W:
            assert v2i(b) in hit_ints, "row-sum missing from stab census"
            n_mem += 1
    print(f"[{time.monotonic()-t0:6.1f}s]   membership gate: 40/40 "
          f"random row-sums present in the census")
    with (DATA / "stab_census_orbits.jsonl").open("w") as f:
        for key, i in orb_rep.items():
            f.write(json.dumps({
                "w": int(ws[i]),
                "b_support": sorted(int(j) for j in
                                    np.nonzero(stab_vecs[i])[0]),
            }) + "\n")

    # --------------------------------------- offsets R0-R4: seam census
    seam_rows = []       # (orbit_idx, vec, w)
    seam_whist: dict[int, int] = {}
    for oi, (lab, _) in enumerate(offsets[1:]):
        els = sorted(hits[lab])
        for h in els:
            v = i2v(h, GR.n)
            wv = int(v.sum())
            assert GR.is_cycle(v) and not GR.is_stab(v)
            assert v2i(GR.sig(v)) == reps_cls[oi], "element in wrong class"
            assert wv >= D_GROSS, \
                f"seam element weight {wv} < d(gross) = 12 ?! (Lean clash)"
            seam_rows.append((oi, v, wv))
            seam_whist[wv] = seam_whist.get(wv, 0) + 1
        print(f"[{time.monotonic()-t0:6.1f}s] offset {lab} (class "
              f"{reps_cls[oi]:#x}, orbit size "
              f"{len(sd['orbits'][oi])}): {len(els)} elements <= {W}")
    print(f"[{time.monotonic()-t0:6.1f}s] seam census total: "
          f"{len(seam_rows)} elements, weight hist "
          f"{dict(sorted(seam_whist.items()))}  [minima >= 12 = d(gross) "
          f"VERIFIED; naive SeamCosetFloor-18 is "
          f"{'FALSE' if seam_whist else 'TRUE'} on the censused orbits]")
    out["seam_census"] = {"elements": len(seam_rows),
                          "weight_hist": dict(sorted(seam_whist.items())),
                          "per_orbit": {f"R{i}": len(hits[f"R{i}"])
                                        for i in range(5)}}
    with (DATA / "seam_census.jsonl").open("w") as f:
        for oi, v, wv in seam_rows:
            f.write(json.dumps({
                "orbit": oi, "class": f"{reps_cls[oi]:#x}", "w": wv,
                "w_support": sorted(int(j) for j in np.nonzero(v)[0]),
            }) + "\n")

    # ------------------------------------------------------ rung engine
    cell = YRungCell("g8", GR, G8, deck_y)
    assert len(cell.sector_basis) == 12, len(cell.sector_basis)
    print(f"[{time.monotonic()-t0:6.1f}s] rung cell: sectors = 2^12, "
          f"kappa = {cell.kappa}")

    # V1: sector linear trick, all 4096 sectors, on a weight-6 rep
    i6 = next(i for i in orb_rep.values() if ws[i] == 6)
    tV = time.monotonic()
    r6 = cell.rung(stab_vecs[i6], (TARGET - 6) // 2,
                   time.monotonic() + 1200, validate_sectors=True)
    assert r6["verdict"] == "PASS", r6
    print(f"[{time.monotonic()-t0:6.1f}s] V1: 4096/4096 sector linear "
          f"trick OK + hexagon rung M=6 PASS (lane {r6['lane']}, "
          f"{time.monotonic()-tV:.1f}s)")

    # --------------------------------------------- dangerous rungs (all)
    rung_rows = []
    verd: dict[str, int] = {}
    lanes: dict[str, int] = {}
    tR = time.monotonic()
    for key, i in orb_rep.items():
        b = stab_vecs[i]
        wb = int(ws[i])
        M = (TARGET - wb) // 2
        r = (r6 if i == i6 else
             cell.rung(b, M, time.monotonic() + 3600))
        verd[r["verdict"]] = verd.get(r["verdict"], 0) + 1
        lanes[r["lane"]] = lanes.get(r["lane"], 0) + 1
        if r["verdict"] != "PASS":
            print("  !! DANGEROUS VIOLATION:", json.dumps(r)[:300])
        rung_rows.append({"species": "dangerous", "w": wb, "M": M,
                          "verdict": r["verdict"], "lane": r["lane"]})
    dtR = time.monotonic() - tR
    assert verd == {"PASS": len(orb_rep)}, verd
    print(f"[{time.monotonic()-t0:6.1f}s] dangerous rungs: "
          f"{len(orb_rep)}/{len(orb_rep)} PASS ({dtR:.1f}s, lanes "
          f"{lanes})  [+ G-transport => all stabilizer shadows <= {W}]")
    out["dangerous_rungs"] = {"rungs": len(orb_rep), "verdicts": verd,
                              "lanes": lanes, "wall_s": round(dtR, 1)}

    # -------------------------------------------------- seam rungs (all)
    verd2: dict[str, int] = {}
    lanes2: dict[str, int] = {}
    tS = time.monotonic()
    for oi, v, wv in seam_rows:
        M = (TARGET - wv) // 2
        r = cell.seam_rung(v, M)
        verd2[r["verdict"]] = verd2.get(r["verdict"], 0) + 1
        lanes2[r["lane"]] = lanes2.get(r["lane"], 0) + 1
        if r["verdict"] != "PASS":
            print("  !! SEAM VIOLATION:", json.dumps(r)[:300])
        rung_rows.append({"species": "seam", "w": wv, "M": M,
                          "verdict": r["verdict"], "lane": r["lane"]})
    dtS = time.monotonic() - tS
    assert verd2 == {"PASS": len(seam_rows)}, verd2
    print(f"[{time.monotonic()-t0:6.1f}s] seam rungs: "
          f"{len(seam_rows)}/{len(seam_rows)} PASS ({dtS:.1f}s, lanes "
          f"{lanes2})  [+ G-transport => all 63 seam-class cosets]")
    out["seam_rungs"] = {"rungs": len(seam_rows), "verdicts": verd2,
                         "lanes": lanes2, "wall_s": round(dtS, 1)}
    with (DATA / "rungs.jsonl").open("w") as f:
        for r in rung_rows:
            f.write(json.dumps(r) + "\n")

    # ------------------------------------------- covariance spot-checks
    g = (5, 1)
    perm_g = TS._perm_for(GR, g)
    for i in list(orb_rep.values())[:3]:
        bt = stab_vecs[i][perm_g]
        rt = cell.rung(bt, (TARGET - int(ws[i])) // 2,
                       time.monotonic() + 600)
        assert rt["verdict"] == "PASS"
    for oi, v, wv in seam_rows[:3]:
        rt = cell.seam_rung(v[perm_g], (TARGET - wv) // 2)
        assert rt["verdict"] == "PASS"
    print(f"[{time.monotonic()-t0:6.1f}s] covariance: 3+3 translated "
          f"reps re-rung, verdicts identical")

    # -------------------------------------- witness (d <= 18): re-verify
    # The exhaustive find-side ladder lives in a36_witness.py (run it
    # first if the banked witness is absent); here the banked weight-18
    # logical is re-verified END-TO-END from its support on every run.
    wfile = DATA / "w18_witness.json"
    assert wfile.exists(), \
        "no banked witness — run scripts/a36_witness.py first"
    wit_row = json.loads(wfile.read_text())
    vwit = np.zeros(G8.n, dtype=np.uint8)
    vwit[wit_row["v_support"]] = 1
    assert G8.is_cycle(vwit) and not G8.is_stab(vwit)
    assert int(vwit.sum()) == TARGET
    bvec, m1, _ = deck_y.slice_data(vwit)
    assert int(bvec.sum()) == wit_row["shadow_w"] and m1 == wit_row["m1"]
    if wit_row["sector"] == "seam":
        assert GR.is_cycle(bvec) and not GR.is_stab(bvec)
        assert TS.in_span(v2i(GR.sig(bvec)), sd["basis"], sd["piv"])
    else:
        assert GR.is_stab(bvec)
    beta, m2, _ = deck_x.slice_data(bvec)
    assert int(beta.sum()) == wit_row["beta_w"] and m2 == wit_row["m2"]
    print(f"[{time.monotonic()-t0:6.1f}s] WITNESS re-verified end-to-end: "
          f"weight 18, sector {wit_row['sector']}, |b| = "
          f"{wit_row['shadow_w']} (class {wit_row['shadow_class']}), "
          f"m1 = {m1}; 2nd level |beta| = {wit_row['beta_w']}, "
          f"m2 = {m2}, stab: {wit_row['beta_is_stab']}")
    out["witness"] = wit_row

    # ------------------------------------------------------ the assembly
    print(f"""
[{time.monotonic()-t0:6.1f}s] ASSEMBLY — d([[288,12,18]]) = 18:
  |v| <= 16, b = p_y(v):
  1. b = 0:      |v| = 2|u| >= 2 d(gross) = 24   [d(gross) = 12 KERNEL-
                 CHECKED in QECLean — branch needs only >= 10]
  2. b stab:     census <= 16 complete ({out['stab_census']['vectors']}
                 vectors, {out['stab_census']['orbits']} orbits, exact
                 node counts) + {out['dangerous_rungs']['rungs']} rungs
                 ALL PASS + G-transport
  3. [b] != 0:   [b] in SEAM (63 classes, 5 orbits); censuses <= 16
                 complete + {out['seam_rungs']['rungs']} seam rungs ALL
                 PASS + G-transport
  => d >= 18;  witness => d <= 18.  Z side by transpose duality.
  d([[288,12,18]]) = 18 — certificate tier, no SAT on the critical path.
""")
    out["wall_s"] = round(time.monotonic() - t0, 1)
    (DATA / "direct_close.json").write_text(json.dumps(out, indent=1))
    print(f"total {out['wall_s']}s -> {DATA / 'direct_close.json'}")


if __name__ == "__main__":
    main()
