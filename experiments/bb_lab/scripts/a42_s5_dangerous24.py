#!/usr/bin/env python3
"""A42 S5 Stage 0b — the W_x sector of the (18,12) weight-24 census
(the dangerous rungs at target 26), completing a42_s5_mixed24.py.

A weight-24 L0 logical v with [v] in W_x = ker p0* has a TRIVIAL fold:
either P0 v = 0 (v = tau0(u), u a weight-12 L1 logical — the 12 banked
orbits, all pure-x and x-windowed, s5_mixed24_W24.json 'controls') or
P0 v = b a nonzero L1 STABILIZER of weight wb <= 24 with overflow
(24 - wb)/2.  This script runs every stab1 orbit rep b (the W = 24 stab
census from the seam phase) through the tdg432 rung dispatch with
TARGET = 26 (overflow < M = (26 - wb)/2, i.e. |v| <= 24), keeping EVERY
violation (= a verified weight-24 nontrivial L0 logical over b), and
classifies the objects: class kind (must be pure-x), gap sector, fibre
profile — the profiles of the period-12 cylinder floor cycles that are
NOT y-pullbacks.

Lanes (soundness as in a40_tdg432_close.py): restricted MITM for
M-1 <= 8 (exact off-support subset sum), kernel-shift for small |b|
with the nontrivial-L1 census complete to WNT = 20 (built here: L2
im-p1* cosets <= 20 + stab2 fibres + tau2), cross-validated on a sample.

Run: cd experiments/bb_lab && caffeinate -i uv run python \
     scripts/a42_s5_dangerous24.py
Needs: data/a42/s5_ckpt_W24_{L2stab,L2all16,L1stab}.jsonl
Output: data/a42/s5_dangerous24.json (+ .log)
"""
from __future__ import annotations

import json
import os
import resource
import sys
import time
from multiprocessing import get_context
from pathlib import Path

import numpy as np

LAB = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(LAB / "src"))
sys.path.insert(0, str(LAB / "scripts"))

import importlib.util  # noqa: E402


def _load(name):
    spec = importlib.util.spec_from_file_location(
        name, Path(__file__).parent / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


M24 = _load("a42_s5_mixed24")          # engine + tower (W from argv: 24)
from bb_lab import cosetbz  # noqa: E402
from bb_lab.tower import (  # noqa: E402
    RungCell, colspace, h1_map, i2v, in_span, rep_for, rref_ints,
    span_points, translation_perms, v2i, validate_banked,
)
from a38_c37xx_freeze import KernelShift, row_lift_v0, whist  # noqa: E402
import a40_s11_compare as C  # noqa: E402

DATA = LAB / "data" / "a42"
W = 24
TARGET = 26
WNT = 20
WALL = 16
T0 = time.monotonic()
LOGF = (DATA / "s5_dangerous24.log").open("a")
M24.LOG = LOGF          # the engine's census logs go to this file too


def log(msg):
    line = f"[{time.monotonic() - T0:8.1f}s] {msg}"
    print(line, flush=True)
    LOGF.write(line + "\n")
    LOGF.flush()


_G: dict = {}


def nontrivial_lifts(args):
    """All NONTRIVIAL L1 lifts of one L2 orbit rep within weight Wb
    (tdg432 lanes: deep for cap <= 6, kernel-shift with the all-class
    <= WALL census otherwise)."""
    from bb_lab.tower import enumerate_lifts_deep
    supp, Wb = args
    d1, L2_, L1_ = M24._G["deck1"], M24._G["L2"], M24._G["L1"]
    b2 = np.zeros(L2_.n, dtype=np.uint8)
    b2[supp] = 1
    wb2 = int(b2.sum())
    cap = (Wb - wb2) // 2
    outl: list[int] = []
    if cap < 0:
        return outl
    if cap <= 6:
        lifts = enumerate_lifts_deep(d1, b2, cap=cap)
        cands = [d1.lift(i2v(v0c, L2_.n), b2) for v0c in lifts]
    else:
        v0p, ovp = row_lift_v0(d1, b2)
        B = wb2 + cap + ovp
        assert B <= WALL, f"no lane |b2|={wb2} cap={cap} B={B}"
        rhs = (d1.RHS @ b2) % 2
        seen: set[int] = set()
        bmask = v2i(b2)
        cands = []
        for v0i in M24._G["ks"].candidates(b2, v0p, cap):
            canon = min(v0i, v0i ^ bmask)
            if canon in seen:
                continue
            seen.add(canon)
            v0 = i2v(v0i, L2_.n)
            assert not (((d1.E @ v0) + rhs) % 2).any()
            cands.append(d1.lift(v0, b2))
    for b1 in cands:
        w_ = int(b1.sum())
        if 0 < w_ <= Wb and not in_span(v2i(b1), L1_.rsHX_b, L1_.rsHX_p):
            outl.append(v2i(b1))
    return outl


def kernel_shift_rung(b: np.ndarray, M: int) -> dict:
    deck0, L0, L1, ks = _G["deck0"], _G["L0"], _G["L1"], _G["ks"]
    wb = int(b.sum())
    v0p, ovp = row_lift_v0(deck0, b)
    B = wb + (M - 1) + ovp
    assert B <= WNT, f"top kernel-shift window {B} > WNT"
    rhs = (deck0.RHS @ b) % 2
    viols = []
    seen: set[int] = set()
    bmask = v2i(b)
    for v0i in ks.candidates(b, v0p, M - 1):
        canon = min(v0i, v0i ^ bmask)
        if canon in seen:
            continue
        seen.add(canon)
        v0 = i2v(v0i, L1.n)
        assert not (((deck0.E @ v0) + rhs) % 2).any()
        ch = (deck0.EMB[0] @ v0 + deck0.EMB[1] @ ((v0 + b) % 2)) % 2
        if in_span(v2i(ch), L0.rsHX_b, L0.rsHX_p):
            continue
        ov = bin(v0i & ~bmask).count("1")
        assert int(ch.sum()) == wb + 2 * ov
        viols.append({"overflow": ov, "weight": int(ch.sum()),
                      "v0_hex": f"{v0i:x}"})
    if viols:
        return {"verdict": "VIOLATION", "lane": "kernel-shift", "w_b": wb,
                "M": M, "n_viol": len(viols), "violations": viols}
    return {"verdict": "PASS", "lane": "kernel-shift", "w_b": wb, "M": M}


def rung_job(args):
    supp, xval = args
    deck0, L1, cell = _G["deck0"], _G["L1"], _G["cell"]
    b = np.zeros(L1.n, dtype=np.uint8)
    b[supp] = 1
    wb = int(b.sum())
    M = (TARGET - wb) // 2
    if M <= 0:
        return None
    lane_note = None
    if (M - 1) <= 4:
        r = cell.rung(b, M, time.monotonic() + 3600, full_viols=True)
        if xval:
            v0p, ovp = row_lift_v0(deck0, b)
            if wb + (M - 1) + ovp <= WNT:
                rks = kernel_shift_rung(b, M)
                assert rks["verdict"] == r["verdict"], "LANE MISMATCH"
                if r["verdict"] == "VIOLATION":
                    assert rks["n_viol"] == r["n_viol"], "LANE COUNT"
                lane_note = "xval"
    else:
        v0p, ovp = row_lift_v0(deck0, b)
        B = wb + (M - 1) + ovp
        if B <= WNT:
            r = kernel_shift_rung(b, M)
            if xval and (M - 1) <= 6:
                rd = cell.rung(b, M, time.monotonic() + 3600,
                               full_viols=True)
                assert rd["verdict"] == r["verdict"], "LANE MISMATCH"
                if r["verdict"] == "VIOLATION":
                    assert rd["n_viol"] == r["n_viol"], "LANE COUNT"
                lane_note = "xval"
        else:
            assert (M - 1) <= 8, f"cell |b|={wb} M={M}: no sound lane"
            r = cell.rung(b, M, time.monotonic() + 3600, full_viols=True)
    chains = []
    if r["verdict"] == "VIOLATION":
        for x in r["violations"]:
            v0i = int(x["v0_hex"], 16)
            chains.append(cell.chain_int(v0i, b))
    return {"w_b": wb, "M": M, "verdict": r["verdict"], "lane": r["lane"],
            "n_viol": r.get("n_viol", 0), "chains": chains,
            "xval": lane_note}


def main():
    log(f"=== a42_s5_dangerous24 W={W} TARGET={TARGET} WNT={WNT} "
        f"pid={os.getpid()}")
    validate_banked(LAB / "data")
    L, (deck0, deck1, deck2), seam_set, S1pset, bases = M24.build_tower()
    L0, L1, L2, L3 = L
    perms0 = translation_perms(L0)
    perms1 = translation_perms(L1)
    perms2 = translation_perms(L2)
    binp = cosetbz.build_kernel()
    ck_stab2 = DATA / f"s5_ckpt_W{W}_L2stab.jsonl"
    ck_all2 = DATA / f"s5_ckpt_W{W}_L2all{WALL}.jsonl"
    ck_stab1 = DATA / f"s5_ckpt_W{W}_L1stab.jsonl"
    ck_ntrv20 = DATA / f"s5_ckpt_W{W}_L1ntrv{WNT}.jsonl"
    for p in (ck_stab2, ck_all2, ck_stab1):
        assert p.exists(), f"missing {p} — run a42_s5_mixed24.py 24 first"
    stab2 = M24.load_reps(ck_stab2, L2.n, lambda v: L2.is_stab(v))
    all2 = M24.load_reps(ck_all2, L2.n)
    stab1 = M24.load_reps(ck_stab1, L1.n, lambda v: L1.is_stab(v))
    log(f"reloaded: stab2 {len(stab2)}, all2<=16 {len(all2)}, stab1 "
        f"{len(stab1)} {whist(stab1)}")
    out: dict = {"W": W, "TARGET": TARGET, "WNT": WNT}

    # ---------------- nontrivial L1 census <= WNT (kernel-shift lane)
    if not ck_ntrv20.exists():
        Mp1 = h1_map(deck1)
        imp1_set = span_points(rref_ints(list(colspace(Mp1)))[0]) - {0}
        assert len(imp1_set) == 63
        imp1 = []
        ims = sorted(imp1_set)
        for lo in range(0, len(ims), 51):
            chunk = ims[lo:lo + 51]
            st, _ = M24.census_stream(
                binp, L2, [(f"C{c}", rep_for(L2, c)) for c in chunk], WNT,
                f"s5_L2imp1_{WNT}_{lo}", perms2)
            for c in chunk:
                for v in st[f"C{c}"].vectors():
                    assert L2.is_cycle(v) and not L2.is_stab(v)
                    imp1.append(v)
        log(f"L2 im-p1* cosets <= {WNT}: {len(imp1)} orbit reps "
            f"{whist(imp1)}")
        M24._init_worker(deck1, L1, L2, seam_set, WALL, all2)
        coll = M24.OrbitStore(L1.n, perms1)
        # tau2 family (nontrivial images only) + fibres over stab2 <= WNT
        # and the im-p1* cosets <= WNT: every nontrivial L1 cycle <= WNT
        # has fold in {0} u stab u im-p1*-cosets (the shadow-class law)
        for u in stab2 + all2:
            if 2 * int(u.sum()) <= WNT:
                b1 = (deck1.TAU @ u) % 2
                assert L1.is_cycle(b1)
                if 0 < int(b1.sum()) <= WNT and not in_span(
                        v2i(b1), L1.rsHX_b, L1.rsHX_p):
                    coll.add(b1)
        plan = ([([int(j) for j in np.nonzero(b)[0]], WNT) for b in stab2
                 if int(b.sum()) <= WNT]
                + [([int(j) for j in np.nonzero(b)[0]], WNT) for b in imp1])
        tF = time.monotonic()
        with get_context("fork").Pool(8) as pool:
            for i, res in enumerate(pool.imap_unordered(nontrivial_lifts,
                                                        plan, chunksize=64)):
                if res:
                    coll.add_vecs(np.array([i2v(x, L1.n) for x in res],
                                           dtype=np.uint8))
        ntrv = coll.vectors()
        for v in ntrv:
            assert L1.is_cycle(v) and not L1.is_stab(v)
        wh = whist(ntrv)
        log(f"L1 nontrivial census <= {WNT}: {len(ntrv)} orbit reps {wh} "
            f"({time.monotonic() - tF:.0f}s)")
        assert wh.get("12") == 12 and min(int(v.sum()) for v in ntrv) == 12
        # regression vs the banked W22 run: <= 18 slice = 5727 orbits
        n18 = sum(c for w_, c in wh.items() if int(w_) <= 18)
        assert n18 == 5727, f"nontrivial <= 18 slice {n18} != banked 5727"
        log("GATE: nontrivial <= 18 slice == banked 5727 orbits")
        M24.save_reps(ck_ntrv20, ntrv)
        out["L1_ntrv20"] = {"orbits": len(ntrv), "whist": wh}
    ntrv = M24.load_reps(ck_ntrv20, L1.n)

    # ------------------------------------------------ the rungs
    cell = RungCell("s5_top", L1, L0, deck0)
    ks_top = KernelShift(deck0, stab1 + ntrv, complete_to=WNT)
    _G.update(deck0=deck0, L0=L0, L1=L1, ks=ks_top, cell=cell)
    plan = []
    nx = 0
    for b in sorted(stab1, key=lambda b: -int(b.sum())):
        wb = int(b.sum())
        xv = (nx < 12 and wb <= 18)
        if xv:
            nx += 1
        plan.append(([int(j) for j in np.nonzero(b)[0]], xv))
    found = M24.OrbitStore(L0.n, perms0)
    verd: dict = {}
    lanes: dict = {}
    per_w: dict = {}
    n_chains = 0
    tR = time.monotonic()
    nproc = max(1, min(8, (os.cpu_count() or 2) - 2))
    log(f"dangerous rungs: {len(plan)} stab1 reps at target {TARGET} on "
        f"{nproc} workers")
    with get_context("fork").Pool(nproc) as pool:
        for i, r in enumerate(pool.imap_unordered(rung_job, plan,
                                                  chunksize=32)):
            if r is None:
                continue
            verd[r["verdict"]] = verd.get(r["verdict"], 0) + 1
            lanes[r["lane"]] = lanes.get(r["lane"], 0) + 1
            key = f"w{r['w_b']}"
            per_w.setdefault(key, {"rungs": 0, "viol": 0, "chains": 0})
            per_w[key]["rungs"] += 1
            if r["verdict"] == "VIOLATION":
                per_w[key]["viol"] += 1
                per_w[key]["chains"] += len(r["chains"])
                vs = []
                for ci in r["chains"]:
                    ch = i2v(ci, L0.n)
                    assert L0.is_cycle(ch) and not L0.is_stab(ch)
                    assert int(ch.sum()) == W
                    vs.append(ch)
                n_chains += len(vs)
                found.add_vecs(np.array(vs, dtype=np.uint8))
            if (i + 1) % 20000 == 0:
                log(f"  ... rungs {i + 1}/{len(plan)} ({time.monotonic() - tR:.0f}s) "
                    f"verdicts {verd} objects {len(found.reps)} RSS "
                    f"{resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 2**30:.2f}")
    log(f"dangerous rungs done: {verd}, lanes {lanes}, "
        f"{n_chains} weight-24 chains -> {len(found.reps)} L0 orbits "
        f"({time.monotonic() - tR:.0f}s)")
    objs = found.vectors()
    hist: dict = {}
    hist_el: dict = {}
    prof_hist: dict = {}
    classes: dict = {}
    recs = []
    for v in objs:
        sig = int(v2i(L0.sig(v)))
        kind = M24.class_kind(sig, bases)
        assert kind == "pure-x", f"non-W_x object with stab fold?! {kind}"
        cells = C.cells_of(L0, v)
        xs = [c[1] for c in cells]
        ys = [c[2] for c in cells]
        nx_, gx = C.gap_structure(xs, 18)
        ny_, gy = C.gap_structure(ys, 12)
        gsec = ("both-gap" if gx >= 4 and gy >= 4 else
                "x-gap" if gx >= 4 else "y-gap" if gy >= 4 else "gap-dense")
        prof = M24.fibre_profile(L0, v)
        osz = M24.orbit_size(v, perms0)
        sh, m, _ = deck0.slice_data(v)
        key = f"{kind}|{gsec}"
        hist[key] = hist.get(key, 0) + 1
        hist_el[key] = hist_el.get(key, 0) + osz
        classes.setdefault(kind, set()).add(sig)
        pk = f"{gsec}|eps{prof['eps']}|n3={prof['n3']}"
        prof_hist[pk] = prof_hist.get(pk, 0) + 1
        recs.append({"kind": kind, "class": sig, "gap_sector": gsec,
                     "gap_x": gx, "gap_y": gy, "fold_w": int(sh.sum()),
                     "overflow": m, "orbit": osz, "profile": prof,
                     "support": [int(j) for j in np.nonzero(v)[0]]})
    log(f"W_x sector at weight 24 (non-pullback, fold = stab): orbits "
        f"{hist}, elements {hist_el}, classes "
        f"{ {k: len(s) for k, s in classes.items()} }")
    log(f"  profiles (gap|eps|n3 -> orbits): {dict(sorted(prof_hist.items()))}")
    out["rungs"] = {"verdicts": verd, "lanes": lanes, "per_w": per_w,
                    "n_chains": n_chains, "L0_orbits": len(objs),
                    "wall_s": round(time.monotonic() - tR, 1)}
    out["census"] = {"hist_orbits": hist, "hist_elements": hist_el,
                     "classes": {k: len(s) for k, s in classes.items()},
                     "profiles": dict(sorted(prof_hist.items())),
                     "objects": recs}
    out["wall_s"] = round(time.monotonic() - T0, 1)
    (DATA / "s5_dangerous24.json").write_text(json.dumps(out, indent=1))
    log(f"DONE {out['wall_s']}s -> {DATA / 's5_dangerous24.json'}")


if __name__ == "__main__":
    main()
