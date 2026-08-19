"""A38 S2: the [[720,4]] freeze — phase 2 (rungs + assembly).

Companion to `a38_c37xx_freeze.py --census-only` (see its docstring):
reloads the checkpointed L1 obligations (orbit-rep supports),
RE-VERIFIES every vector from scratch (cycle membership, stabilizer /
seam-class membership, weight window, the A30 d(L1)/d(L2) tripwires),
rebuilds the kernel-shift window from the loaded stabilizer census,
and runs the dangerous + seam rungs at M = (TARGET - |b|)/2 with the
same lane assignment, cross-validations, covariance spot-checks, and
assembly as the one-shot script.  Split exists only because the
harness kills background tasks at ~1 h wall; the certificate content
is identical (the loaded data re-passes the same in-line asserts).

Usage: python a38_c37xx_rungs.py 22

Output: data/a38/c37xx/freeze_W{W}.json + rungs_W{W}.jsonl
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

from bb_lab.group import AbelianGroup  # noqa: E402
from bb_lab.tower import (  # noqa: E402
    AxisDeck, RungCell, TowerCode, _as_support, colspace, fold_support,
    h1_map, i2v, in_span, perm_for, rref_ints, span_points, v2i,
)
from a38_c37xx_freeze import (  # noqa: E402
    A30, D_L1, D_L2, FOLDS, KernelShift, SPEC, row_lift_v0, whist,
)

DATA = LAB / "data" / "a38" / "c37xx"


def main() -> None:
    W = int(sys.argv[1]) if len(sys.argv) > 1 else 22
    assert W in (18, 22)
    TARGET = W + 2
    t0 = time.monotonic()
    out: dict = {"W": W, "target": TARGET, "phase": "rungs"}

    def log(msg: str) -> None:
        print(f"[{time.monotonic()-t0:7.1f}s] {msg}", flush=True)

    # ---------------- tower + structure (same asserts as phase 1)
    orders, As, Bs = SPEC
    G_top = AbelianGroup(orders)
    levels = [(orders, _as_support(As, G_top), _as_support(Bs, G_top))]
    for axis, newmod in FOLDS:
        plm, pA, pB = levels[-1]
        nlm = tuple(newmod if a == axis else plm[a]
                    for a in range(len(plm)))
        levels.append((nlm, fold_support(pA, axis, newmod),
                       fold_support(pB, axis, newmod)))
    L = [TowerCode(f"L{i}", *lv) for i, lv in enumerate(levels)]
    assert [c.k for c in L] == [4, 4, 4, 4]
    deck0 = AxisDeck(L[0], L[1], 0)
    bank = json.loads((A30 / "decide_37a70e02.json").read_text())
    assert bank["d_base"] == D_L2 and bank["floor"] == D_L1
    Mp0 = h1_map(deck0)
    SEAM = colspace(Mp0)
    Sbb, Sbp = rref_ints(list(SEAM))
    seam_set = span_points(Sbb) - {0}
    assert len(seam_set) == 3
    log("tower rebuilt; structure asserts PASS (phase-1 gates recorded "
        "in freeze_W{}_census.json)".format(W))

    # ---------------- reload + RE-VERIFY the checkpointed obligations
    stab1_reps: list[np.ndarray] = []
    for line in (DATA / f"ckpt_W{W}_stab1.jsonl").open():
        r = json.loads(line)
        b = np.zeros(L[1].n, dtype=np.uint8)
        b[r["support"]] = 1
        assert int(b.sum()) == r["w"] and 0 < r["w"] <= W
        assert L[1].is_cycle(b)
        assert in_span(v2i(b), L[1].rsHX_b, L[1].rsHX_p), \
            "checkpoint row is not an L1 stabilizer"
        stab1_reps.append(b)
    seam1_reps: list[np.ndarray] = []
    for line in (DATA / f"ckpt_W{W}_seam1.jsonl").open():
        r = json.loads(line)
        v = np.zeros(L[1].n, dtype=np.uint8)
        v[r["support"]] = 1
        assert int(v.sum()) == r["w"] and 0 < r["w"] <= W
        assert L[1].is_cycle(v)
        assert not in_span(v2i(v), L[1].rsHX_b, L[1].rsHX_p)
        assert v2i(L[1].sig(v)) in seam_set, "not a SEAM-class element"
        assert r["w"] >= D_L1, \
            f"seam element {r['w']} < d(L1): A30 contradiction!"
        seam1_reps.append(v)
    mu1 = min(int(b.sum()) for b in stab1_reps)
    assert mu1 >= 6
    log(f"checkpoint reloaded + re-verified: {len(stab1_reps)} stab "
        f"orbit reps {whist(stab1_reps)} (mu1 = {mu1}); "
        f"{len(seam1_reps)} seam reps {whist(seam1_reps)} "
        f"(all >= {D_L1})")
    out["L1_stab"] = {"orbits": len(stab1_reps),
                      "whist": whist(stab1_reps), "mu1": mu1}
    out["L1_seam"] = {"orbits": len(seam1_reps),
                      "whist": whist(seam1_reps)}

    # ---------------- rungs (identical logic to the one-shot script)
    cell = RungCell("c37xx_top", L[1], L[0], deck0)
    assert len(cell.sector_basis) == 4
    ks_top = KernelShift(deck0, stab1_reps, complete_to=W)

    def kernel_shift_rung(b: np.ndarray, M: int) -> dict:
        wb = int(b.sum())
        v0p, ovp = row_lift_v0(deck0, b)
        B = wb + (M - 1) + ovp
        assert B < D_L1, \
            f"kernel-shift window {B} >= d(L1): |b|={wb} M={M}"
        rhs = (deck0.RHS @ b) % 2
        viols = []
        seen_v0: set[int] = set()
        bmask = v2i(b)
        for v0i in ks_top.candidates(b, v0p, M - 1):
            canon = min(v0i, v0i ^ bmask)
            if canon in seen_v0:
                continue
            seen_v0.add(canon)
            v0 = i2v(v0i, L[1].n)
            assert not (((deck0.E @ v0) + rhs) % 2).any()
            ch = (deck0.EMB[0] @ v0
                  + deck0.EMB[1] @ ((v0 + b) % 2)) % 2
            if in_span(v2i(ch), L[0].rsHX_b, L[0].rsHX_p):
                continue
            ov = bin(v0i & ~bmask).count("1")
            wt = int(ch.sum())
            assert wt == wb + 2 * ov, "slice identity violated"
            viols.append({"overflow": ov, "weight": wt,
                          "v0_hex": f"{v0i:x}"})
        if viols:
            return {"verdict": "VIOLATION", "lane": "kernel-shift",
                    "w_b": wb, "M": M, "violations": viols[:5],
                    "n_viol": len(viols),
                    "min_weight": min(x["weight"] for x in viols)}
        return {"verdict": "PASS", "lane": "kernel-shift", "w_b": wb,
                "M": M}

    rung_rows = []
    verd: dict[str, int] = {}
    lanes: dict[str, int] = {}
    viol_finds: list[dict] = []
    first3: list[tuple[np.ndarray, int, str]] = []
    n_xval = xval_cheap = xval_deep = 0
    n_done = 0
    tR = time.monotonic()
    for b in sorted(stab1_reps, key=lambda b: -int(b.sum())):
        wb = int(b.sum())
        M = (TARGET - wb) // 2
        if M <= 0:
            continue
        if (M - 1) <= 4:
            r = cell.rung(b, M, time.monotonic() + 3600)
            if xval_cheap < 3:
                v0p, ovp = row_lift_v0(deck0, b)
                if wb + (M - 1) + ovp < D_L1:
                    rks = kernel_shift_rung(b, M)
                    assert rks["verdict"] == r["verdict"], \
                        f"LANE MISMATCH |b|={wb} M={M}"
                    if r["verdict"] == "VIOLATION":
                        assert rks["n_viol"] == r["n_viol"]
                    n_xval += 1
                    xval_cheap += 1
        else:
            v0p, ovp = row_lift_v0(deck0, b)
            B = wb + (M - 1) + ovp
            if B < D_L1 and B <= W:
                r = kernel_shift_rung(b, M)
                if xval_deep < 3 and (M - 1) <= 6:
                    rd = cell.rung(b, M, time.monotonic() + 3600)
                    assert rd["verdict"] == r["verdict"], \
                        f"LANE MISMATCH |b|={wb} M={M}"
                    if r["verdict"] == "VIOLATION":
                        assert rd["n_viol"] == r["n_viol"]
                    n_xval += 1
                    xval_deep += 1
            else:
                assert (M - 1) <= 6, \
                    f"cell |b|={wb} M={M}: no sound lane (B={B})"
                r = cell.rung(b, M, time.monotonic() + 3600)
        verd[r["verdict"]] = verd.get(r["verdict"], 0) + 1
        lanes[r["lane"]] = lanes.get(r["lane"], 0) + 1
        rung_rows.append({"species": "dangerous", "w": wb, "M": M,
                          "verdict": r["verdict"], "lane": r["lane"]})
        if r["verdict"] == "VIOLATION":
            viol_finds.append(r)
        if len(first3) < 3:
            first3.append((b, M, r["verdict"]))
        n_done += 1
        if n_done % 20000 == 0:
            log(f"  ... dangerous rungs {n_done}/{len(stab1_reps)} "
                f"({time.monotonic()-tR:.0f}s)")
    dtR = time.monotonic() - tR
    log(f"dangerous rungs: {sum(verd.values())} at target {TARGET} "
        f"({dtR:.1f}s): verdicts {verd}, lanes {lanes}; xval {n_xval}")
    out["dangerous_rungs"] = {"rungs": sum(verd.values()),
                              "verdicts": verd, "lanes": lanes,
                              "xval_cells": n_xval,
                              "wall_s": round(dtR, 1)}

    verd2: dict[str, int] = {}
    seam_finds: list[dict] = []
    for w_el in sorted(seam1_reps, key=lambda v: -int(v.sum())):
        ww = int(w_el.sum())
        M = (TARGET - ww) // 2
        if M <= 0:
            continue
        r = cell.seam_rung(w_el, M)
        verd2[r["verdict"]] = verd2.get(r["verdict"], 0) + 1
        rung_rows.append({"species": "seam", "w": ww, "M": M,
                          "verdict": r["verdict"], "lane": r["lane"]})
        if r["verdict"] == "VIOLATION":
            viol_finds.append(r)
            seam_finds.append(r)
    log(f"seam rungs: {sum(verd2.values())} at target {TARGET}: {verd2}")
    out["seam_rungs"] = {"rungs": sum(verd2.values()),
                         "verdicts": verd2}
    if len(rung_rows) <= 50000:
        with (DATA / f"rungs_W{W}.jsonl").open("w") as f:
            for rr in rung_rows:
                f.write(json.dumps(rr) + "\n")
    else:
        agg: dict[tuple, int] = {}
        for rr in rung_rows:
            key = (rr["species"], rr["w"], rr["M"], rr["verdict"],
                   rr["lane"])
            agg[key] = agg.get(key, 0) + 1
        with (DATA / f"rungs_W{W}.jsonl").open("w") as f:
            for (sp, w_, M_, vd, ln), cnt in sorted(agg.items()):
                f.write(json.dumps({"species": sp, "w": w_, "M": M_,
                                    "verdict": vd, "lane": ln,
                                    "count": cnt}) + "\n")

    # covariance spot-checks
    g = (7, 1)
    perm_g = perm_for(L[1], g)
    for b, M, v_verd in first3:
        bt = b[perm_g]
        v0p, ovp = row_lift_v0(deck0, bt)
        B = int(bt.sum()) + (M - 1) + ovp
        rt = (kernel_shift_rung(bt, M) if B < D_L1 and B <= W
              else cell.rung(bt, M, time.monotonic() + 600))
        assert rt["verdict"] == v_verd, "covariance broken"
    log("covariance: 3 translated reps re-rung, verdicts transport")

    # ---------------- assembly
    n_viol_cells = len(viol_finds)
    if n_viol_cells == 0:
        log(f"ASSEMBLY (W = {W}): branch b = 0 dead (>= {2*D_L1}); "
            f"dangerous closed ({sum(verd.values())} PASS + "
            f"G-transport); seam closed ({sum(verd2.values())} PASS + "
            f"G-transport) => NO nontrivial X-logical of weight "
            f"<= {W}: d([[720,4]]) >= {TARGET} at certificate tier "
            f"(consuming d(L1) = 20, d(L2) = 10 at A30 certificate "
            f"tier). Z side by BB transpose duality.")
        out["verdict"] = {"floor": TARGET, "all_pass": True}
    else:
        wts = []
        for r in viol_finds:
            wts.extend(x["weight"] for x in r["violations"])
            if "min_weight" in r:
                wts.append(r["min_weight"])
        wmin = min(wts)
        # the lightest find, re-verified end-to-end + banked
        best = None
        for r in viol_finds:
            for x in r["violations"]:
                if best is None or x["weight"] < best["weight"]:
                    best = dict(x)
                    best["cell_w"] = r.get("w_b", r.get("w_w"))
                    best["cell_M"] = r["M"]
                    best["lane"] = r["lane"]
        log(f"ASSEMBLY (W = {W}): VIOLATIONS — explicit nontrivial "
            f"L0-logicals found (every candidate re-verified in-line: "
            f"E-system, non-stab, slice identity), min weight {wmin} "
            f"over {n_viol_cells} cells (seam cells: "
            f"{len(seam_finds)}). d([[720,4]]) <= {wmin}; with the "
            f"Q1 floor d >= 20: "
            + (f"d = 20 EXACT — the A14 SS13 freeze is CONFIRMED at "
               f"certificate tier." if wmin == 20 else
               f"d in [20, {wmin}]."))
        out["verdict"] = {"upper": wmin, "all_pass": False,
                          "viol_cells": n_viol_cells,
                          "seam_viol_cells": len(seam_finds),
                          "witness": best}

    out["wall_s"] = round(time.monotonic() - t0, 1)
    (DATA / f"freeze_W{W}.json").write_text(json.dumps(out, indent=1))
    log(f"total {out['wall_s']}s -> {DATA / f'freeze_W{W}.json'}")


if __name__ == "__main__":
    main()
