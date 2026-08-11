"""A33 Part 3: falsify-first validation against every banked A20 artifact.

  (a) the 1,655-class dangerous census (H1): loads as lab Y4-stabilizers
      (hard assert, all rows), decomposes through the rung-1 fold
      (transport + slice asserts; shadow-orbit compression measured);
      v7_complete_classes == m_census as G-canonical class sets
  (b) H2 re-derivation: all 1,655 per-class dangerous floors at target
      20 as DETERMINISTIC rung certificates (restricted lanes), replacing
      the banked SAT UNSATs; M values asserted == banked m_req, verdicts
      asserted == banked UNSAT 1:1
  (c) witnesses: the weight-20 Y8 witness (H6) re-verified + decomposed
      in tower coordinates; the weight-10 Y4 and weight-6 Y2 ladder
      witnesses re-verified

Output: data/a33/banked_validation.json + h2_rungs.jsonl
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

from a30_rung_pass import v2i  # noqa: E402
import a32_tower_slice as TS  # noqa: E402
from a32_sectorAC_full import batch_keys  # noqa: E402
from a33_tower_cells import build_tower  # noqa: E402
from a33_rung_cell import YRungCell  # noqa: E402

MAIN = Path("/Users/stavanjain/Code/qec-lab/experiments/bb_lab/data")
DATA = LAB / "data" / "a33"


def main():
    t0 = time.monotonic()
    out: dict = {}
    Y2, Y4, Y8, deck_top, deck_bot = build_tower()
    perms_y4 = TS._translation_perms(Y4)
    perms_y2 = TS._translation_perms(Y2)

    # ------------------------------------------------- (a) census loads
    census = []
    for line in (MAIN / "a20" / "m_census_classes.jsonl").open():
        r = json.loads(line)
        if "b_support" in r:
            census.append(r)
    assert len(census) == 1655
    vecs = np.zeros((1655, Y4.n), dtype=np.uint8)
    for i, e in enumerate(census):
        vecs[i, e["b_support"]] = 1
        assert int(vecs[i].sum()) == e["w"]
    # hard assert: every row is a lab Y4-stabilizer
    for i in range(1655):
        assert Y4.is_stab(vecs[i]), f"census row {i} not a Y4 stabilizer"
    wh = {}
    for e in census:
        wh[e["w"]] = wh.get(e["w"], 0) + 1
    assert wh == {6: 1, 10: 6, 12: 33, 14: 54, 16: 375, 18: 1186}, wh
    print(f"[{time.monotonic()-t0:5.1f}s] (a) census: 1,655/1,655 rows are "
          f"lab Y4-stabilizers; bands {wh}")

    # v7 == m_census as G-canonical class sets
    v7 = []
    for line in (MAIN / "a20" / "v7_complete_classes.jsonl").open():
        r = json.loads(line)
        if "b_support" in r:
            v7.append(r)
    assert len(v7) == 1655
    v7v = np.zeros((1655, Y4.n), dtype=np.uint8)
    for i, e in enumerate(v7):
        v7v[i, e["b_support"]] = 1
    k_sat = {bytes(k) for k in batch_keys(vecs, perms_y4)}
    k_v7 = {bytes(k) for k in batch_keys(v7v, perms_y4)}
    assert len(k_sat) == 1655 and k_sat == k_v7, \
        "v7 classes != SAT census classes"
    print(f"[{time.monotonic()-t0:5.1f}s]     v7 analytic census == SAT "
          f"census as G-canonical class sets (1,655 == 1,655)")

    # decomposition through the rung-1 fold
    hist: dict[tuple[int, int], int] = {}
    beta_zero = 0
    betas = np.zeros((1655, Y2.n), dtype=np.uint8)
    for i in range(1655):
        beta, m2, _ = deck_bot.slice_data(vecs[i])
        assert Y2.is_stab(beta), "shadow of a Y4-stab not a Y2-stab"
        betas[i] = beta
        wb = int(beta.sum())
        hist[(int(vecs[i].sum()), wb)] = hist.get(
            (int(vecs[i].sum()), wb), 0) + 1
        if wb == 0:
            beta_zero += 1
    keys_b = batch_keys(betas, perms_y2)
    nz = [i for i in range(1655) if betas[i].any()]
    orb_beta = {bytes(keys_b[i]) for i in nz}
    print(f"[{time.monotonic()-t0:5.1f}s]     rung-1 decomposition: all "
          f"transport+slice asserts pass; beta=0 records {beta_zero}; "
          f"distinct nonzero Y2-shadow orbits {len(orb_beta)} "
          f"(compression {len(nz)/max(len(orb_beta),1):.1f}x)")
    out["census"] = {
        "records": 1655, "bands": wh, "v7_equals_sat": True,
        "beta_zero_records": beta_zero,
        "beta_orbits_nonzero": len(orb_beta),
        "joint_hist": {f"{w},{wb}": c for (w, wb), c in sorted(hist.items())},
    }

    # ------------------------------------------------- (b) H2 rungs
    floors = {}
    for line in (MAIN / "a20" / "m_floors_results.jsonl").open():
        r = json.loads(line)
        floors[tuple(r["b_support"])] = r
    assert len(floors) == 1655
    cell_top = YRungCell("top", Y4, Y8, deck_top)
    lanes: dict[str, list[float]] = {}
    verd: dict[str, int] = {}
    rows = []
    n_match = 0
    for i, e in enumerate(census):
        b = vecs[i]
        M = (21 - e["w"]) // 2
        banked = floors[tuple(e["b_support"])]
        assert banked["m_req"] == M, \
            f"M mismatch: banked {banked['m_req']} vs {M}"
        assert banked["verdict"] == "UNSAT"
        tc = time.monotonic()
        r = cell_top.rung(b, M, time.monotonic() + 900)
        dt = time.monotonic() - tc
        lanes.setdefault(f"w{e['w']}", []).append(dt)
        verd[r["verdict"]] = verd.get(r["verdict"], 0) + 1
        if r["verdict"] == "PASS":
            n_match += 1
        else:
            print("  !! NON-PASS (banked said UNSAT):", json.dumps(r)[:300])
        rows.append({"w": e["w"], "M": M, "verdict": r["verdict"],
                     "lane": r.get("lane"), "secs": round(dt, 4)})
        if i % 400 == 0 and i:
            print(f"    ... {i}/1655 ({time.monotonic()-t0:.0f}s)")
    assert verd == {"PASS": 1655}, verd
    print(f"[{time.monotonic()-t0:5.1f}s] (b) H2 re-derived: 1,655/1,655 "
          f"rungs PASS at target 20 (deterministic restricted lanes) — "
          f"banked SAT agreement {n_match}/1655")
    for k in sorted(lanes, key=lambda s: int(s[1:])):
        v = lanes[k]
        print(f"      {k}: n={len(v)} tot={sum(v):.1f}s max={max(v):.3f}s")
    out["h2"] = {"rungs": 1655, "verdicts": verd,
                 "banked_agreement": n_match,
                 "lanes": {k: {"n": len(v), "tot_s": round(sum(v), 2)}
                           for k, v in lanes.items()}}
    with (DATA / "h2_rungs.jsonl").open("w") as f:
        for r in rows:
            f.write(json.dumps(r) + "\n")

    # ------------------------------------------------- (c) witnesses
    w20 = np.load(MAIN / "a20" / "y8_weight20_witness.npy").astype(np.uint8)
    assert w20.shape == (Y8.n,) and int(w20.sum()) == 20
    assert Y8.is_cycle(w20) and not Y8.is_stab(w20)
    b, m1, v0 = deck_top.slice_data(w20)
    prof = {"shadow_w": int(b.sum()), "m1": m1}
    assert int(b.sum()) == 0 and m1 == 10, prof
    u = (deck_top.EMB[0].T @ w20) % 2  # sheet 0 = the diagonal preimage
    assert Y4.is_cycle(u) and not Y4.is_stab(u) and int(u.sum()) == 10
    beta_u, m2u, _ = deck_bot.slice_data(u)
    prof["sheet_u_w"] = 10
    prof["u_shadow_w"] = int(beta_u.sum())
    prof["u_m2"] = m2u
    print(f"[{time.monotonic()-t0:5.1f}s] (c) weight-20 Y8 witness: "
          f"verified nontrivial X-logical; tower profile b=0 (tau-diagonal)"
          f", sheet u = weight-10 Y4 logical, u-shadow |beta|="
          f"{prof['u_shadow_w']} m2={prof['u_m2']}")
    wit10 = None
    for line in (MAIN / "a20" / "y144_ladder.log").read_text().splitlines():
        try:
            wit10 = json.loads(line)["witness"]
        except (json.JSONDecodeError, KeyError):
            continue
    u10 = np.zeros(Y4.n, dtype=np.uint8)
    u10[wit10] = 1
    assert Y4.is_cycle(u10) and not Y4.is_stab(u10) and u10.sum() == 10
    wit6 = None
    for line in (MAIN / "a20" / "y72_ladder.log").read_text().splitlines():
        try:
            wit6 = json.loads(line)["witness"]
        except (json.JSONDecodeError, KeyError):
            continue
    u6 = np.zeros(Y2.n, dtype=np.uint8)
    u6[wit6] = 1
    assert Y2.is_cycle(u6) and not Y2.is_stab(u6) and u6.sum() == 6
    print(f"[{time.monotonic()-t0:5.1f}s]     ladder witnesses re-verified: "
          f"weight-10 Y4 logical, weight-6 Y2 logical")
    out["witnesses"] = {"y8_w20": "verified, tau-diagonal over u10",
                        "y8_w20_profile": prof,
                        "y4_w10": "verified", "y2_w6": "verified"}

    out["wall_s"] = round(time.monotonic() - t0, 1)
    (DATA / "banked_validation.json").write_text(json.dumps(out, indent=1))
    print(f"total {out['wall_s']}s -> {DATA / 'banked_validation.json'}")


if __name__ == "__main__":
    main()
