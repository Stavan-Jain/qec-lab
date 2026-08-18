"""A36 Part 6: the weight-18 witness, tower-native (d <= 18 half).

The direct-close hunt established: NO weight-18 logical exists at
overflow 1 over any weight-16 shadow (374 seam elements + 375 stab
orbit reps, transport-complete), nor at overflow 2 over the 18 w14 seam
elements.  d = 18 (published SAT-exact) forces a witness at exactly one
of the remaining strata, so the ladder here is EXHAUSTIVE:

  L1  |b| = 14 stab (54 reps)   m1 = 2   rung at M = 3
      |b| = 12 seam (3 els)     m1 = 3   seam rung at M = 4
      |b| = 12 stab (33 reps)   m1 = 3   rung at M = 4
      |b| = 10 stab (6 reps)    m1 = 4   rung at M = 5
      |b| = 6  stab (1 rep)     m1 = 6   rung at M = 7 (cap-6 lane)
  L2  |b| = 18 flat (m1 = 0): W = 18 extension of the 6-offset BZ pass
      (r-pair (9,8), exact node asserts) -> weight-18 seam elements
      (seam rung at M = 1: any flat cycle IS a witness, transport-
      nontrivial) and weight-18 stab orbit reps (rung at M = 1: flat
      nontrivial logicals).

Every find is asserted >= 18 (a lighter one would refute the published
record — and the d >= 18 floor already closed by a36_direct_close);
finding NOTHING anywhere would refute d = 18 itself.  Either way the
session cannot end agnostic.

Output: data/a36/w18_witness.json (+ witness_hunt.json trace)
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
TARGET = 18


def main():
    t0 = time.monotonic()
    out: dict = {"ladder": []}
    G8, GR, B72, B36, deck_y, deck_x, deck_b = build_tower()
    My = h1_map(deck_y)
    sd = seam_info(GR, My)
    reps_cls = sd["reps"]
    rep_vecs = [rep_for(GR, c) for c in reps_cls]
    perms = TS._translation_perms(GR)
    cell = YRungCell("g8w", GR, G8, deck_y)

    stab_reps = []
    for line in (DATA / "stab_census_orbits.jsonl").open():
        r = json.loads(line)
        b = np.zeros(GR.n, dtype=np.uint8)
        b[r["b_support"]] = 1
        stab_reps.append((r["w"], b))
    seam_els = []
    for line in (DATA / "seam_census.jsonl").open():
        r = json.loads(line)
        v = np.zeros(GR.n, dtype=np.uint8)
        v[r["w_support"]] = 1
        seam_els.append((r["w"], v))
    print(f"[{time.monotonic()-t0:6.1f}s] loaded {len(stab_reps)} stab "
          f"orbit reps + {len(seam_els)} seam elements")

    witness = None
    finds: list[int] = []

    def consider(vio, sector, bvec):
        nonlocal witness
        for x in vio:
            assert x["weight"] >= TARGET, \
                f"sub-{TARGET} find {x} — contradicts the closed floor!"
            finds.append(x["weight"])
            if x["weight"] == TARGET and witness is None:
                witness = {"sector": sector,
                           "shadow_w": int(bvec.sum()),
                           "b_int": v2i(bvec),
                           "v0_hex": x["v0_hex"],
                           "overflow": x["overflow"]}

    # ------------------------------------------------ L1: shadows <= 16
    ladder = [("stab", 14, 3), ("seam", 12, 4), ("stab", 12, 4),
              ("stab", 10, 5), ("stab", 6, 7)]
    for kind, wb, M in ladder:
        if witness:
            break
        tL = time.monotonic()
        n_run = 0
        if kind == "stab":
            for w, b in stab_reps:
                if w != wb or witness:
                    continue
                r = cell.rung(b, M, time.monotonic() + 1800)
                n_run += 1
                if r["verdict"] == "VIOLATION":
                    consider(r["violations"], "dangerous", b)
        else:
            for w, v in seam_els:
                if w != wb or witness:
                    continue
                r = cell.seam_rung(v, M)
                n_run += 1
                if r["verdict"] == "VIOLATION":
                    consider(r["violations"], "seam", v)
        row = {"stratum": f"{kind}|b|={wb}@M={M}", "rungs": n_run,
               "finds": len(finds), "witness": witness is not None,
               "wall_s": round(time.monotonic() - tL, 1)}
        out["ladder"].append(row)
        print(f"[{time.monotonic()-t0:6.1f}s] L1 {row['stratum']}: "
              f"{n_run} rungs, finds so far {len(finds)}, "
              f"witness: {witness is not None} ({row['wall_s']}s)")

    # ------------------------------------------------ L2: |b| = 18 flat
    if witness is None:
        print(f"[{time.monotonic()-t0:6.1f}s] L2: extending censuses to "
              f"W = 18 (r-pair (9,8)) for the flat stratum ...")
        binp = build_kernel()
        I1, G1, I2, G2, kappa = disjoint_info_sets(GR.HX)
        assert kappa == 66
        offsets = [("S", np.zeros(GR.n, dtype=np.uint8))] + \
            [(f"R{i}", rv) for i, rv in enumerate(rep_vecs)]
        hits18: dict[str, set[int]] = {lab: set() for lab, _ in offsets}
        nodes = 0
        for wi, (window, Gs, r) in enumerate([(I1, G1, 9), (I2, G2, 8)]):
            bases = []
            for lab, tv in offsets:
                cb = coset_base(Gs, window, tv)
                if int(cb.sum()) == 18:
                    hits18[lab].add(v2i(cb))
                bases.append(cb)
            res = run_window(binp, f"a36_wit_w{wi}", Gs, bases, r, 18,
                             time.monotonic() + 3600)
            nodes += res["nodes"]
            for j, hx in res.pop("hit_rows"):
                v = unpack3(hx, GR.n)
                if int(v.sum()) == 18:
                    hits18[offsets[j][0]].add(v2i(v))
            print(f"[{time.monotonic()-t0:6.1f}s]   window {wi}: nodes "
                  f"{res['nodes']:.3e}, wall {res['wall_s']}s")
        exp = sum(math.comb(66, s) for s in range(1, 10)) + \
            sum(math.comb(66, s) for s in range(1, 9))
        assert nodes == exp
        # seam-18 first: any flat cycle is a witness (transport-nontrivial)
        n18_seam = sum(len(hits18[f"R{i}"]) for i in range(5))
        n_flat = 0
        for i in range(5):
            if witness:
                break
            for h in sorted(hits18[f"R{i}"]):
                v = i2v(h, GR.n)
                assert GR.is_cycle(v) and not GR.is_stab(v)
                r = cell.seam_rung(v, 1)
                n_flat += 1
                if r["verdict"] == "VIOLATION":
                    consider(r["violations"], "seam", v)
                if witness:
                    break
        print(f"[{time.monotonic()-t0:6.1f}s] L2 seam-18: {n18_seam} "
              f"elements, {n_flat} flat rungs probed, witness: "
              f"{witness is not None}")
        out["ladder"].append({"stratum": "seam|b|=18@M=1",
                              "elements": n18_seam, "rungs": n_flat,
                              "witness": witness is not None})
        if witness is None:
            svecs = np.array([i2v(h, GR.n)
                              for h in sorted(hits18["S"])], dtype=np.uint8)
            keys = batch_keys(svecs, perms)
            orb = {}
            for i, k in enumerate(keys):
                orb.setdefault(bytes(k), i)
            print(f"[{time.monotonic()-t0:6.1f}s] L2 stab-18: "
                  f"{len(svecs)} vectors, {len(orb)} orbit reps; "
                  f"flat rungs ...")
            n_run = 0
            for i in orb.values():
                r = cell.rung(svecs[i], 1, time.monotonic() + 1800)
                n_run += 1
                if r["verdict"] == "VIOLATION":
                    consider(r["violations"], "dangerous", svecs[i])
                if witness:
                    break
            out["ladder"].append({"stratum": "stab|b|=18@M=1",
                                  "orbits": len(orb), "rungs": n_run,
                                  "witness": witness is not None})
            print(f"[{time.monotonic()-t0:6.1f}s] L2 stab-18: {n_run} "
                  f"flat rungs, witness: {witness is not None}")

    assert witness is not None, \
        "EXHAUSTIVE ladder found no weight-18 logical — d > 18?!"

    # --------------------------------------- verify + bank the witness
    bvec = i2v(witness["b_int"], GR.n)
    v0 = i2v(int(witness["v0_hex"], 16), GR.n)
    vwit = deck_y.lift(v0, bvec)
    assert G8.is_cycle(vwit) and not G8.is_stab(vwit)
    assert int(vwit.sum()) == TARGET
    bchk, m1, _ = deck_y.slice_data(vwit)
    assert (bchk == bvec).all() and 2 * m1 + int(bvec.sum()) == TARGET
    beta, m2, _ = deck_x.slice_data(bvec)
    wit_row = {
        "weight": TARGET, "sector": witness["sector"],
        "shadow_w": witness["shadow_w"], "m1": int(m1),
        "shadow_class": f"{v2i(GR.sig(bvec)):#x}",
        "beta_w": int(beta.sum()), "m2": int(m2),
        "beta_is_stab": bool(B72.is_stab(beta)),
        "b_support": sorted(int(j) for j in np.nonzero(bvec)[0]),
        "v_support": sorted(int(j) for j in np.nonzero(vwit)[0]),
        "finds_total": len(finds), "finds_min": min(finds),
    }
    (DATA / "w18_witness.json").write_text(json.dumps(wit_row, indent=1))
    print(f"[{time.monotonic()-t0:6.1f}s] WITNESS verified end-to-end: "
          f"weight 18, sector {wit_row['sector']}, |b| = "
          f"{wit_row['shadow_w']} (class {wit_row['shadow_class']}), "
          f"m1 = {wit_row['m1']}; 2nd level |beta| = {wit_row['beta_w']} "
          f"m2 = {wit_row['m2']} stab: {wit_row['beta_is_stab']}; all "
          f"{len(finds)} finds >= 18")
    out["witness"] = wit_row
    out["wall_s"] = round(time.monotonic() - t0, 1)
    (DATA / "witness_hunt.json").write_text(json.dumps(out, indent=1))
    print(f"total {out['wall_s']}s -> {DATA / 'witness_hunt.json'}")


if __name__ == "__main__":
    main()
