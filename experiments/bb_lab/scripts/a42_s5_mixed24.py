#!/usr/bin/env python3
"""A42 S5 Stage 0 — the (18,12) weight-24 census OUTSIDE W_x, class-kind
and fibre-profile tabulated (the calibration ground truth for HM at q=4).

Setting.  L0 = (18,12) = [[432,12,24]] (member (2,1); d = 24 certified
by the tdg432 W=22 descent sweep, which enumerates NOTHING at weight 24).
H1(L0) = W_x (+) W_y with W_x = ker p0* for the y-fold p0 : L0 -> L1 =
(18,6) (A42 §2.13: pushforward is zero at a = 2; checked here as a span
equality).  So

    [v] not in W_x   <=>   P0 v is a NONTRIVIAL logical of L1 in a SEAM
                           class (SEAM = im p0*, 63 classes),

and |P0 v| = |v| - 2 c with c = #doubled fibres.  The banked seam census
(sweep_W22_census.json: seam1 <= 22 has whist {22: 68}) says the seam
classes have minimum 22, hence a weight-24 non-W_x logical has c <= 1:
its fold is a seam element of weight 22 (one doubled fibre) or 24 (a
section).  This script therefore

  (1) rebuilds the L2 = (9,6) censuses at W = 24 (stab, the 3 S1'
      classes; all-class <= 16 for the kernel-shift lane) by coset-BZ,
      STREAMING the hit files into translation-orbit reps (the W = 22 run
      materialized 3.9M S1' elements; at 24 that is ~25M — RSS-unsafe);
  (2) descends them to the L1 SEAM census <= 24 (tau2 family + fibres,
      the tdg432 v2 lanes verbatim; the weight-22 slice must reproduce
      the banked 68 orbits — a regression gate);
  (3) enumerates EVERY weight-24 L0 cycle over every seam rep (sections
      of the weight-24 elements, overflow-1 lifts of the weight-22 ones;
      RungCell's exact restricted lane, all solutions kept), re-verifies
      each (cycle, non-stab, weight, fold), reduces to L0 translation
      orbits, and classifies each object: class kind (pure-y / pure-d /
      mixed via the S11 transfer images W_x, W_y, W_d — pure-x is
      IMPOSSIBLE here and asserted so), gap sector (Lemma K), and the
      Z_3-fibre profile (n1, n2, n3) -> (|S|, |s|, |S cap s|, eps, n3).

Positive controls: the L12-stack witness is pure-x and x-windowed; every
object found is non-W_x; the seam-22 slice reproduces 68 orbits; the
fibre identity w = 2|S| + 3|s| - 4|S cap s| holds per object.

Run:  cd experiments/bb_lab && caffeinate -i uv run python \
      scripts/a42_s5_mixed24.py [--phase census|seam|rungs|all]
Data: data/a42/s5_mixed24*.json, s5_ckpt_*.jsonl (resumable).
"""
from __future__ import annotations

import gc
import json
import math
import os
import resource
import subprocess
import sys
import time
from multiprocessing import get_context
from pathlib import Path

import numpy as np

LAB = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(LAB / "src"))
sys.path.insert(0, str(LAB / "scripts"))

from bb_lab import cosetbz  # noqa: E402
from bb_lab.tower import (  # noqa: E402
    AxisDeck, RungCell, TowerCode, batch_keys, colspace,
    enumerate_lifts_deep, gf2_rank, h1_map, i2v, in_span, kernel_ints,
    reduce_int, rep_for, rref_ints, span_eq, span_points,
    translation_mat, translation_perms, v2i, validate_banked,
)
from a38_c37xx_freeze import KernelShift, row_lift_v0, whist  # noqa: E402
import a40_s11_compare as C  # noqa: E402

DATA = LAB / "data" / "a42"
DATA.mkdir(parents=True, exist_ok=True)
WORK = DATA / "s5_bzwork"
WORK.mkdir(parents=True, exist_ok=True)
W = int(sys.argv[1]) if len(sys.argv) > 1 and sys.argv[1].isdigit() else 24
assert W in (22, 24)
WALL = min(16, W)
RSS_CAP_GB = 3.0
T0 = time.monotonic()
LOG = (DATA / f"s5_mixed24_W{W}.log").open("a")


def log(msg: str) -> None:
    line = f"[{time.monotonic() - T0:8.1f}s] {msg}"
    print(line, flush=True)
    LOG.write(line + "\n")
    LOG.flush()


def rss_gb() -> float:
    r = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return r / (1024 ** 3)   # darwin reports bytes


def rss_guard(tag: str = "") -> None:
    g = rss_gb()
    if g > RSS_CAP_GB:
        log(f"RSS GUARD {g:.2f} GB > {RSS_CAP_GB} at {tag}: aborting")
        sys.exit(3)


# ----------------------------------------------------------- hex <-> vecs
def hex_rows_to_vecs(hexes: list[str], n: int) -> np.ndarray:
    """(N, n) uint8 from the kernel's 48-hex-char rows (bit i = coord i)."""
    raw = bytes.fromhex("".join(hexes))
    arr = np.frombuffer(raw, dtype=np.uint8).reshape(len(hexes), 24)
    bits = np.unpackbits(arr[:, ::-1], axis=1, bitorder="little")
    return np.ascontiguousarray(bits[:, :n])


def vec_to_int_rows(vecs: np.ndarray) -> list[int]:
    return [v2i(v) for v in vecs]


class OrbitStore:
    """Translation-orbit reps accumulated chunk-wise (batch_keys)."""

    def __init__(self, n: int, perms):
        self.n, self.perms = n, perms
        self.reps: dict[bytes, np.ndarray] = {}

    def add_vecs(self, vecs: np.ndarray) -> None:
        if vecs.shape[0] == 0:
            return
        keys = batch_keys(vecs, self.perms)
        for i in range(vecs.shape[0]):
            k = bytes(keys[i])
            if k not in self.reps:
                self.reps[k] = vecs[i].copy()

    def add(self, v: np.ndarray) -> None:
        self.add_vecs(v.reshape(1, -1))

    def vectors(self) -> list[np.ndarray]:
        return list(self.reps.values())


def orbit_size(v: np.ndarray, perms) -> int:
    vi = v2i(v)
    stab = sum(1 for p in perms if v2i(v[p]) == vi)
    assert len(perms) % stab == 0
    return len(perms) // stab


def census_stream(binp, code: TowerCode, offsets, Wc: int, tag: str,
                  perms, deadline_s: float = 7200.0,
                  chunk: int = 200_000) -> tuple[dict, dict]:
    """Complete two-window coset-BZ census of the cosets `offsets`
    (label, vector) at weight <= Wc, streamed into orbit reps per label.
    Returns ({label: OrbitStore}, stats).  Node counts asserted exactly
    (the certificate); hit files deleted as consumed."""
    I1, G1, I2, G2, kappa = cosetbz.disjoint_info_sets(code.HX)
    assert not (set(I1) & set(I2))
    r1, r2 = Wc // 2, max(Wc - Wc // 2 - 1, 0)
    stores = {lab: OrbitStore(code.n, perms) for lab, _ in offsets}
    stats = {"kappa": kappa, "windows": []}
    threads = 8
    for wi, (window, Gs, r) in enumerate([(I1, G1, r1), (I2, G2, r2)]):
        bases = []
        for lab, tv in offsets:
            cb = cosetbz.coset_base(Gs, window, tv)
            if 0 < int(cb.sum()) <= Wc:
                stores[lab].add(cb)
            bases.append(cb)
        assert len(bases) <= cosetbz.NOFF_MAX
        n = Gs.shape[1]
        assert n <= cosetbz.NMAX
        wtag = f"{tag}_w{wi}"
        mat = WORK / f"{wtag}.mat"
        with mat.open("w") as f:
            f.write(f"{kappa} {r} {Wc} {len(bases)}\n")
            for row in Gs:
                w0, w1, w2 = cosetbz.pack3(row)
                f.write(f"{w2:x} {w1:x} {w0:x}\n")
            for b in bases:
                w0, w1, w2 = cosetbz.pack3(b)
                f.write(f"{w2:x} {w1:x} {w0:x}\n")
        t0 = time.monotonic()
        out = subprocess.run([str(binp), str(mat), str(threads),
                              str(WORK / wtag)],
                             capture_output=True, text=True,
                             timeout=deadline_s)
        assert out.returncode == 0, out.stderr
        parts = dict(p.split("=") for p in out.stdout.strip().split())
        nodes, hits = int(parts["nodes"]), int(parts["hits"])
        expect = sum(math.comb(kappa, s) for s in range(1, r + 1))
        assert nodes == expect, f"node count {nodes} != {expect}"
        tk = time.monotonic() - t0
        log(f"  {wtag}: kappa {kappa} r {r} nodes {nodes:,} (exact) hits "
            f"{hits:,} kernel {tk:.0f}s; streaming ...")
        seen_hits = 0
        for t in range(threads):
            p = WORK / f"{wtag}_t{t:02d}.hits"
            if not p.exists():
                continue
            buf_j: list[int] = []
            buf_h: list[str] = []

            def flush():
                if not buf_h:
                    return
                vecs = hex_rows_to_vecs(buf_h, code.n)
                if seen_hits <= chunk:      # decoder self-test, first chunk
                    for ii in range(0, len(buf_h), max(1, len(buf_h) // 50)):
                        ref = cosetbz.unpack3(buf_h[ii], code.n)
                        assert (ref == vecs[ii]).all(), "hex decoder"
                js = np.array(buf_j)
                for li, (lab, _) in enumerate(offsets):
                    sel = vecs[js == li]
                    sel = sel[sel.sum(axis=1) > 0]
                    stores[lab].add_vecs(sel)
                buf_j.clear()
                buf_h.clear()

            with p.open() as fh:
                for line in fh:
                    j, hx = line.split()
                    buf_j.append(int(j))
                    buf_h.append(hx)
                    seen_hits += 1
                    if len(buf_h) >= chunk:
                        flush()
                        rss_guard(wtag)
            flush()
            p.unlink()
        assert seen_hits == hits, (seen_hits, hits)
        mat.unlink()
        stats["windows"].append({"window": wi, "r": r, "nodes": nodes,
                                 "hits": hits, "kernel_s": round(tk, 1),
                                 "total_s": round(time.monotonic() - t0,
                                                  1)})
        log(f"  {wtag}: streamed {hits:,} hits in "
            f"{time.monotonic() - t0 - tk:.0f}s; reps so far "
            f"{sum(len(s.reps) for s in stores.values())}")
    return stores, stats


# --------------------------------------------------------- checkpoints
def save_reps(path: Path, vecs: list[np.ndarray], extra=None) -> None:
    with path.open("w") as f:
        for i, b in enumerate(vecs):
            rec = {"w": int(b.sum()),
                   "support": [int(j) for j in np.nonzero(b)[0]]}
            if extra is not None:
                rec.update(extra[i])
            f.write(json.dumps(rec) + "\n")


def load_reps(path: Path, n: int, checker=None) -> list[np.ndarray]:
    reps = []
    for line in path.open():
        r = json.loads(line)
        v = np.zeros(n, dtype=np.uint8)
        v[r["support"]] = 1
        assert int(v.sum()) == r["w"]
        if checker:
            checker(v)
        reps.append(v)
    return reps


# --------------------------------------------------------------- tower
def build_tower():
    L0 = C.member_code(18, 12)
    LMS = [(18, 6), (9, 6), (9, 3)]
    L = [L0] + [TowerCode(f"L{i + 1}", lm, C.red(C.A_L, lm),
                          C.red(C.B_L, lm)) for i, lm in enumerate(LMS)]
    assert [c.k for c in L] == [12, 12, 8, 8]
    assert [c.n for c in L] == [432, 216, 108, 54]
    deck0 = AxisDeck(L[0], L[1], 1)
    deck1 = AxisDeck(L[1], L[2], 0)
    deck2 = AxisDeck(L[2], L[3], 1)
    for c in L:
        assert not any(int(kv.sum()) % 2 for kv in c.kerHZ), "parity!"
    Mp0, Mp1 = h1_map(deck0), h1_map(deck1)
    Mt0 = h1_map(deck0, tau=True)
    SEAMb, _ = rref_ints(list(colspace(Mp0)))
    comp = (Mp1 @ Mp0) % 2
    S1pb, _ = rref_ints(list(colspace(comp)))
    assert len(SEAMb) == 6 and len(S1pb) == 2
    seam_set = span_points(SEAMb) - {0}
    S1pset = span_points(S1pb) - {0}
    assert len(seam_set) == 63 and len(S1pset) == 3
    assert span_eq(kernel_ints(Mt0), list(colspace(Mp0)))
    assert (translation_mat(L[0], deck0.sigma)
            == np.eye(L[0].k, dtype=np.uint8)).all(), "sigma* != id"
    assert deck0.twisted() and deck1.twisted() and deck2.twisted()
    # ker p0* == W_x (the S11 transfer image in direction (0,1))
    kerp0 = kernel_ints(Mp0)
    Wx, _ = C.image_classes(L0, (0, 1), K=4)
    Wy, _ = C.image_classes(L0, (1, 0), K=4)
    Wd, _ = C.image_classes(L0, (1, 1), K=4)
    assert len(Wx) == len(Wy) == len(Wd) == 6
    assert span_eq(kerp0, list(Wx)), "ker p0* != W_x"
    log("tower: k 12/12/8/8, SEAM 63 classes, S1' 3 classes, top rung "
        "(R)+exact+sigma*=id; ker p0* == W_x (dim 6) VERIFIED")
    return L, (deck0, deck1, deck2), seam_set, S1pset, (Wx, Wy, Wd)


# ---------------------------------------------------------- fibre worker
_G: dict = {}


def _init_worker(deck1, L1, L2, seam_set, WALL_, ks_reps):
    _G["deck1"] = deck1
    _G["L1"], _G["L2"] = L1, L2
    _G["seam_set"] = seam_set
    _G["ks"] = KernelShift(deck1, ks_reps, complete_to=WALL_)


def _classify_l1(b1: np.ndarray, out: dict) -> None:
    L1 = _G["L1"]
    w_ = int(b1.sum())
    if w_ == 0 or w_ > W:
        return
    if in_span(v2i(b1), L1.rsHX_b, L1.rsHX_p):
        if w_ <= W:
            out["stab"].append(v2i(b1))
        return
    out["dmin"] = min(out["dmin"], w_)
    if v2i(L1.sig(b1)) in _G["seam_set"]:
        out["seam"].append(v2i(b1))
    elif w_ <= 12:
        out["ntrv12"].append(v2i(b1))


def fibre_job(args):
    """Lifts of one L2 orbit rep b2 to L1 within weight W; returns the
    classified lift ints."""
    supp, Wb = args
    deck1, L2 = _G["deck1"], _G["L2"]
    b2 = np.zeros(L2.n, dtype=np.uint8)
    b2[supp] = 1
    wb2 = int(b2.sum())
    out = {"stab": [], "seam": [], "ntrv12": [], "dmin": 10 ** 9,
           "lane": None}
    cap = (Wb - wb2) // 2
    if cap < 0:
        return out
    if cap <= 6:
        out["lane"] = f"deep{cap}"
        lifts = enumerate_lifts_deep(deck1, b2, cap=cap)
        for v0c, m2 in lifts.items():
            _classify_l1(deck1.lift(i2v(v0c, L2.n), b2), out)
    else:
        v0p, ovp = row_lift_v0(deck1, b2)
        B = wb2 + cap + ovp
        if B <= WALL:
            out["lane"] = f"ks{cap}"
            rhs = (deck1.RHS @ b2) % 2
            seen: set[int] = set()
            bmask = v2i(b2)
            for v0i in _G["ks"].candidates(b2, v0p, cap):
                canon = min(v0i, v0i ^ bmask)
                if canon in seen:
                    continue
                seen.add(canon)
                v0 = i2v(v0i, L2.n)
                assert not (((deck1.E @ v0) + rhs) % 2).any()
                _classify_l1(deck1.lift(v0, b2), out)
        else:
            assert cap <= 8, f"no sound lane: |b2|={wb2} cap={cap} B={B}"
            out["lane"] = f"deep{cap}"
            lifts = enumerate_lifts_deep(deck1, b2, cap=cap)
            for v0c, m2 in lifts.items():
                _classify_l1(deck1.lift(i2v(v0c, L2.n), b2), out)
    return out


# ------------------------------------------------------------ profiles
def fibre_profile(L0: TowerCode, v: np.ndarray):
    """Z_3-fibre profile on Z_12 = Z_3 x Z_4: fibre of (blk, x, y) is
    {y, y+4, y+8}."""
    fib: dict = {}
    for blk, x, y in C.cells_of(L0, v):
        key = (blk, x, y % 4)
        fib[key] = fib.get(key, 0) + 1
    n = [0, 0, 0, 0]
    for k in fib.values():
        n[k] += 1
    w = int(v.sum())
    assert n[1] + 2 * n[2] + 3 * n[3] == w
    S, s, Ss = n[1] + n[2], n[1] + n[3], n[1]
    assert 2 * S + 3 * s - 4 * Ss == w
    eps = S - 12
    D = 4 * Ss - 3 * s
    cols = sorted({(blk, x) for (blk, x, _) in fib})
    return {"n1": n[1], "n2": n[2], "n3": n[3], "S": S, "s": s,
            "Ss": Ss, "eps": eps, "D": D, "ncols": len(cols)}


def class_kind(sig: int, bases) -> str:
    Wx, Wy, Wd = bases
    kind = "mixed"
    for tag, b in (("x", Wx), ("y", Wy), ("d", Wd)):
        bb, pp = rref_ints(list(b))
        if in_span(sig, bb, pp):
            kind = "pure-" + tag
    return kind


# ------------------------------------------------------------------ main
def main() -> None:
    phase = "all"
    if "--phase" in sys.argv:
        phase = sys.argv[sys.argv.index("--phase") + 1]
    log(f"=== a42_s5_mixed24 phase={phase} W={W} pid={os.getpid()} "
        f"cpus={os.cpu_count()}")
    validate_banked(LAB / "data")
    log("validate_banked: PASS")
    L, (deck0, deck1, deck2), seam_set, S1pset, bases = build_tower()
    L0, L1, L2, L3 = L
    perms1 = translation_perms(L1)
    perms2 = translation_perms(L2)
    perms0 = translation_perms(L0)
    binp = cosetbz.build_kernel()
    out: dict = {"W": W, "WALL": WALL}

    ck_stab2 = DATA / f"s5_ckpt_W{W}_L2stab.jsonl"
    ck_s1p2 = DATA / f"s5_ckpt_W{W}_L2s1p.jsonl"
    ck_all2 = DATA / f"s5_ckpt_W{W}_L2all{WALL}.jsonl"
    ck_seam1 = DATA / f"s5_ckpt_W{W}_L1seam.jsonl"
    ck_stab1 = DATA / f"s5_ckpt_W{W}_L1stab.jsonl"
    ck_ntrv1 = DATA / f"s5_ckpt_W{W}_L1ntrv12.jsonl"

    # ------------------------------------------------ 1. L2 censuses
    if phase in ("census", "all") and not (ck_stab2.exists()
                                           and ck_s1p2.exists()
                                           and ck_all2.exists()):
        # 1a. stab <= 24
        if not ck_stab2.exists():
            st, stats = census_stream(
                binp, L2, [("S", np.zeros(L2.n, np.uint8))], W,
                "s5_L2stab24", perms2)
            reps = st["S"].vectors()
            for v in reps[:: max(1, len(reps) // 100)]:
                assert L2.is_stab(v)
            n_el = sum(orbit_size(v, perms2) for v in reps)
            log(f"L2 stab <= {W}: {len(reps)} orbit reps {whist(reps)} = "
                f"{n_el:,} elements; mu2 = "
                f"{min(int(v.sum()) for v in reps)}")
            out["L2_stab"] = {"orbits": len(reps), "whist": whist(reps),
                              "elements": n_el, "bz": stats}
            save_reps(ck_stab2, reps)
            del st, reps
            gc.collect()
            rss_guard("after stab2")
        # 1b. S1' cosets <= 24 (3 classes)
        if not ck_s1p2.exists():
            offs = [(f"C{c}", rep_for(L2, c)) for c in sorted(S1pset)]
            st, stats = census_stream(binp, L2, offs, W, "s5_L2s1p24",
                                      perms2)
            reps = []
            for c in sorted(S1pset):
                vs = st[f"C{c}"].vectors()
                for v in vs:
                    assert L2.is_cycle(v) and not L2.is_stab(v)
                reps.extend(vs)
            store = OrbitStore(L2.n, perms2)
            for i in range(0, len(reps), 50_000):
                store.add_vecs(np.array(reps[i:i + 50_000], dtype=np.uint8))
            reps = store.vectors()
            for v in reps:
                assert v2i(L2.sig(v)) in S1pset
            n_el = sum(orbit_size(v, perms2) for v in reps)
            log(f"L2 S1'-cosets <= {W} (3 classes): {len(reps)} orbit reps "
                f"{whist(reps)} = {n_el:,} elements")
            out["L2_s1p"] = {"orbits": len(reps), "whist": whist(reps),
                             "elements": n_el, "bz": stats}
            save_reps(ck_s1p2, reps)
            del st, reps, store
            gc.collect()
            rss_guard("after s1p")
        # 1c. all classes <= WALL (kernel-shift lane + tau2 sources)
        if not ck_all2.exists():
            store = OrbitStore(L2.n, perms2)
            stabreps = load_reps(ck_stab2, L2.n)
            for v in stabreps:
                if int(v.sum()) <= WALL:
                    store.add(v)
            alll2 = sorted(range(1, 1 << L2.k))
            d_l2 = None
            CH = 51
            for lo in range(0, len(alll2), CH):
                chunk = alll2[lo:lo + CH]
                st, _ = census_stream(
                    binp, L2, [(f"C{c}", rep_for(L2, c)) for c in chunk],
                    WALL, f"s5_L2all16_{lo}", perms2)
                for c in chunk:
                    for v in st[f"C{c}"].vectors():
                        assert L2.is_cycle(v) and not L2.is_stab(v)
                        w_ = int(v.sum())
                        d_l2 = w_ if d_l2 is None else min(d_l2, w_)
                        store.add(v)
            assert d_l2 == 10, f"d(L2) = {d_l2} != banked 10"
            reps = store.vectors()
            log(f"L2 all-class <= {WALL}: {len(reps)} orbit reps; "
                f"d(L2) = {d_l2} EXACT (banked: 7780 reps)")
            out["L2_all16"] = {"orbits": len(reps), "d_L2": d_l2}
            save_reps(ck_all2, reps)
        (DATA / f"s5_mixed24_W{W}_census.json").write_text(
            json.dumps(out, indent=1))
        if phase == "census":
            log("census phase complete")
            return

    # ------------------------------------------- 2. L1 seam census <= 24
    if phase in ("seam", "all") and not ck_seam1.exists():
        stab2_reps = load_reps(ck_stab2, L2.n, lambda v: L2.is_stab(v))
        s1p_reps = load_reps(ck_s1p2, L2.n)
        all2_reps = load_reps(ck_all2, L2.n)
        log(f"reloaded L2 reps: stab {len(stab2_reps)}, S1' "
            f"{len(s1p_reps)}, all<=16 {len(all2_reps)}")
        coll = {"stab": OrbitStore(L1.n, perms1),
                "seam": OrbitStore(L1.n, perms1),
                "ntrv12": OrbitStore(L1.n, perms1)}
        dmin = [10 ** 9]

        def absorb(res: dict) -> None:
            for key in ("stab", "seam", "ntrv12"):
                if res[key]:
                    vecs = np.array([i2v(x, L1.n) for x in res[key]],
                                    dtype=np.uint8)
                    coll[key].add_vecs(vecs)
            dmin[0] = min(dmin[0], res["dmin"])

        # tau2 family: L2 elements with 2w <= W (reps suffice: tau is
        # translation-covariant and every collected set is invariant)
        _init_worker(deck1, L1, L2, seam_set, WALL, all2_reps)
        n_tau = 0
        for u in stab2_reps + all2_reps:
            if 2 * int(u.sum()) <= W:
                b1 = (deck1.TAU @ u) % 2
                assert L1.is_cycle(b1)
                res = {"stab": [], "seam": [], "ntrv12": [], "dmin": 10 ** 9}
                _classify_l1(b1, res)
                absorb(res)
                n_tau += 1
        log(f"L1 tau2 family: {n_tau} rep sources (2w <= {W})")
        plan = ([([int(j) for j in np.nonzero(b)[0]], W) for b in stab2_reps]
                + [([int(j) for j in np.nonzero(b)[0]], W) for b in s1p_reps])
        plan.sort(key=lambda t: -len(t[0]))
        nproc = max(1, min(8, (os.cpu_count() or 2) - 2))
        log(f"L2 -> L1 fibres: {len(plan)} (stab2 + S1') on {nproc} "
            f"workers")
        lanes: dict[str, int] = {}
        tF = time.monotonic()
        # fork: children inherit _G (AxisDeck holds closures — unpicklable)
        with get_context("fork").Pool(nproc) as pool:
            for i, res in enumerate(pool.imap_unordered(fibre_job, plan,
                                                        chunksize=64)):
                absorb(res)
                if res["lane"]:
                    lanes[res["lane"]] = lanes.get(res["lane"], 0) + 1
                if (i + 1) % 20000 == 0:
                    log(f"  ... fibres {i + 1}/{len(plan)} "
                        f"({time.monotonic() - tF:.0f}s) seam reps "
                        f"{len(coll['seam'].reps)} stab reps "
                        f"{len(coll['stab'].reps)} RSS {rss_gb():.2f}")
                    rss_guard("fibres")
        # CONTROL: the nontrivial L1 slice <= 12 is complete only with the
        # im-p1* fibres (fold class in im p1*, 63 classes) — a small
        # census at W = 12 (cap <= 1 since d(L2) = 10) closes it and
        # re-derives d(L1) = 12.
        Mp1 = h1_map(deck1)
        imp1_set = span_points(rref_ints(list(colspace(Mp1)))[0]) - {0}
        assert len(imp1_set) == 63
        imp1_reps = []
        ims = sorted(imp1_set)
        for lo in range(0, len(ims), 51):
            chunk = ims[lo:lo + 51]
            st, _ = census_stream(
                binp, L2, [(f"C{c}", rep_for(L2, c)) for c in chunk], 12,
                f"s5_L2imp1_12_{lo}", perms2)
            for c in chunk:
                for v in st[f"C{c}"].vectors():
                    assert L2.is_cycle(v) and not L2.is_stab(v)
                    imp1_reps.append(v)
        log(f"L2 im-p1* cosets <= 12 (control census): {len(imp1_reps)} "
            f"orbit reps {whist(imp1_reps)}")
        for b2 in imp1_reps:
            res = fibre_job(([int(j) for j in np.nonzero(b2)[0]], 12))
            absorb(res)
        seam1 = coll["seam"].vectors()
        stab1 = coll["stab"].vectors()
        ntrv1 = coll["ntrv12"].vectors()
        for v in seam1:
            assert L1.is_cycle(v) and not L1.is_stab(v)
            assert v2i(L1.sig(v)) in seam_set
        assert dmin[0] == 12, f"d(L1) = {dmin[0]} != banked 12"
        assert whist(ntrv1).get("12") == 12, \
            f"weight-12 L1 orbit count {whist(ntrv1)} != banked 12"
        wh = whist(seam1)
        log(f"L1 SEAM census <= {W}: {len(seam1)} orbit reps {wh}; stab "
            f"reps {len(stab1)} {whist(stab1)}; nontrivial <= 12: "
            f"{len(ntrv1)} {whist(ntrv1)}; d(L1) = {dmin[0]}; lanes "
            f"{lanes}; {time.monotonic() - tF:.0f}s")
        assert wh.get("22", 0) == 68, \
            f"REGRESSION: seam-22 slice {wh.get('22')} != banked 68"
        assert min(int(v.sum()) for v in seam1) == 22, "seam min != 22"
        log("GATE: seam weight-22 slice == banked 68 orbits; seam min 22")
        save_reps(ck_seam1, seam1,
                  [{"class": int(v2i(L1.sig(v)))} for v in seam1])
        save_reps(ck_stab1, stab1)
        save_reps(ck_ntrv1, ntrv1)
        out["L1"] = {"seam_orbits": len(seam1), "seam_whist": wh,
                     "stab_orbits": len(stab1), "stab_whist": whist(stab1),
                     "ntrv12_orbits": len(ntrv1), "d_L1": dmin[0],
                     "lanes": lanes,
                     "fibre_s": round(time.monotonic() - tF, 1)}
        (DATA / f"s5_mixed24_W{W}_seam.json").write_text(
            json.dumps(out, indent=1))
        if phase == "seam":
            log("seam phase complete")
            return

    # ------------------------------------------------ 3. rungs at L0
    def chk_seam(v):
        assert L1.is_cycle(v) and not L1.is_stab(v)
        assert v2i(L1.sig(v)) in seam_set
    seam1 = load_reps(ck_seam1, L1.n, chk_seam)
    log(f"seam reps reloaded + re-verified: {len(seam1)} {whist(seam1)}")
    cell = RungCell("s5_top", L1, L0, deck0)
    found = OrbitStore(L0.n, perms0)
    n_raw = 0
    per_seam = []
    tR = time.monotonic()
    for si, w_el in enumerate(sorted(seam1, key=lambda v: -int(v.sum()))):
        ww = int(w_el.sum())
        M = (W + 2 - ww) // 2          # overflow < M  <=>  |v| <= W
        assert M >= 1
        rhs = (cell.RHS_OP @ w_el) % 2
        v0p = cell.solve_E(rhs)
        cnt = 0
        if v0p is not None:
            wmask = v2i(w_el)
            rhs_i = v2i(rhs)
            vs = []
            for X in sorted(cell._hits_X(w_el, rhs, M - 1)):
                for v0_int in cell._expand_X(w_el, X, rhs_i):
                    ov = bin(v0_int & ~wmask).count("1")
                    if ov > M - 1:
                        continue
                    ch = i2v(cell.chain_int(v0_int, w_el), L0.n)
                    assert L0.is_cycle(ch) and not L0.is_stab(ch)
                    assert int(ch.sum()) == ww + 2 * ov == W
                    sh, m, _ = deck0.slice_data(ch)
                    assert (sh == w_el).all() and m == ov
                    vs.append(ch)
                    cnt += 1
            if vs:
                found.add_vecs(np.array(vs, dtype=np.uint8))
        n_raw += cnt
        per_seam.append({"w": ww, "M": M, "n_lifts": cnt})
        if (si + 1) % 200 == 0:
            log(f"  ... seam rungs {si + 1}/{len(seam1)}: {n_raw} weight-24 "
                f"lifts so far, {len(found.reps)} L0 orbits "
                f"({time.monotonic() - tR:.0f}s)")
    objs = found.vectors()
    log(f"L0 weight-{W} cycles over seam reps: {n_raw} lifts (of seam "
        f"orbit reps) -> {len(objs)} L0 translation-orbit reps "
        f"({time.monotonic() - tR:.0f}s)")

    # ------------------------------------------------ 4. classify
    Wx, Wy, Wd = bases
    hist: dict[str, int] = {}
    hist_el: dict[str, int] = {}
    classes: dict[str, set] = {}
    prof_hist: dict[str, int] = {}
    records = []
    for v in objs:
        sig = int(v2i(L0.sig(v)))
        kind = class_kind(sig, bases)
        assert kind != "pure-x", "a W_x class with nonzero fold?!"
        cells = C.cells_of(L0, v)
        xs = [c[1] for c in cells]
        ys = [c[2] for c in cells]
        nx, gx = C.gap_structure(xs, 18)
        ny, gy = C.gap_structure(ys, 12)
        gsec = ("both-gap" if gx >= 4 and gy >= 4 else
                "x-gap" if gx >= 4 else "y-gap" if gy >= 4 else "gap-dense")
        assert gx < 4, "x-windowed object outside W_x?!"
        prof = fibre_profile(L0, v)
        osz = orbit_size(v, perms0)
        ncomp = len(C.components(L0, cells))
        key = f"{kind}|{gsec}"
        hist[key] = hist.get(key, 0) + 1
        hist_el[key] = hist_el.get(key, 0) + osz
        classes.setdefault(kind, set()).add(sig)
        pk = f"{kind}|eps{prof['eps']}|n3={prof['n3']}"
        prof_hist[pk] = prof_hist.get(pk, 0) + 1
        records.append({"kind": kind, "class": sig, "gap_sector": gsec,
                        "gap_x": gx, "gap_y": gy, "nx": nx, "ny": ny,
                        "ncomp": ncomp, "orbit": osz, "profile": prof,
                        "support": [int(j) for j in np.nonzero(v)[0]]})
    log(f"class-kind x gap-sector (orbit reps): {hist}")
    log(f"  elements: {hist_el}")
    log(f"  distinct classes: { {k: len(s) for k, s in classes.items()} }")
    log(f"  fibre profiles (kind|eps|n3 -> orbits): "
        f"{dict(sorted(prof_hist.items()))}")

    # ------------------------------------------------ 5. controls
    wit = C.l12_stack(L0)
    wsig = int(v2i(L0.sig(wit)))
    wkind = class_kind(wsig, bases)
    sh, m, _ = deck0.slice_data(wit)
    assert wkind == "pure-x" and L1.is_stab(sh), "L12 stack control"
    xs = [c[1] for c in C.cells_of(L0, wit)]
    _, gxw = C.gap_structure(xs, 18)
    assert gxw >= 4
    wprof = fibre_profile(L0, wit)
    log(f"CONTROL: L12 stack is pure-x, fold is an L1 stab (weight "
        f"{int(sh.sum())}, overflow {m}), x-gap {gxw}; profile {wprof}")
    # tau0 family: pullbacks of the weight-12 L1 logicals (non-seam)
    tau_recs = []
    if ck_ntrv1.exists():
        for u in load_reps(ck_ntrv1, L1.n):
            if int(u.sum()) != 12 or v2i(L1.sig(u)) in seam_set:
                continue
            v = (deck0.TAU @ u) % 2
            assert L0.is_cycle(v) and not L0.is_stab(v) and v.sum() == 24
            k = class_kind(int(v2i(L0.sig(v))), bases)
            xs = [c[1] for c in C.cells_of(L0, v)]
            _, gx = C.gap_structure(xs, 18)
            tau_recs.append({"kind": k, "gap_x": gx,
                             "profile": fibre_profile(L0, v)})
        kinds = {}
        for r in tau_recs:
            kk = f"{r['kind']}|gx{r['gap_x']}|eps{r['profile']['eps']}"
            kinds[kk] = kinds.get(kk, 0) + 1
        log(f"CONTROL: tau0 pullbacks of the {len(tau_recs)} weight-12 L1 "
            f"orbit reps: {kinds}")
    out["rungs"] = {"n_lifts_of_reps": n_raw, "L0_orbits": len(objs),
                    "per_seam": per_seam,
                    "wall_s": round(time.monotonic() - tR, 1)}
    out["census"] = {"hist_orbits": hist, "hist_elements": hist_el,
                     "classes": {k: len(s) for k, s in classes.items()},
                     "profiles": dict(sorted(prof_hist.items())),
                     "objects": records}
    out["controls"] = {"l12_stack": {"kind": wkind, "gap_x": gxw,
                                     "profile": wprof},
                       "tau0_family": tau_recs}
    out["wall_s"] = round(time.monotonic() - T0, 1)
    out["max_rss_gb"] = round(rss_gb(), 2)
    (DATA / f"s5_mixed24_W{W}.json").write_text(json.dumps(out, indent=1))
    log(f"DONE total {out['wall_s']}s, max RSS {out['max_rss_gb']} GB -> "
        f"{DATA / f's5_mixed24_W{W}.json'}")


if __name__ == "__main__":
    main()
