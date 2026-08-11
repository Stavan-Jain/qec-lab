"""A33 Part 4b: H5 closed by DESCENT — the A32 Part-5 pattern one rung down.

The trisection of the seam-element census through rung 1 COLLAPSES on this
(R)-tower (measured in Part 2, `a33_tower_cells`):

  - SEAM cap ker p1* = 0  =>  no seam element has a homologically trivial
    Y2-shadow: the beta = 0 (tau-diagonal) branch and the beta-stabilizer
    branch are BOTH DEAD analytically (the A24-reachability analogue,
    sharper: one rank computation kills two of three branches);
  - p1* is INJECTIVE on SEAM (W2 = p1*(SEAM) has full dim 4), so the
    shadows of class-0x6 elements occupy EXACTLY the single Y2-class
    0x36 — the census moves to ONE coset at n = 72.

Descent enumeration (complete):  every seam element w (class 0x6,
|w| <= 18) has beta = p1(w) in the 0x36-coset census at <= 18, and w is a
lift of beta with overflow m2 = (|w| - |beta|)/2 <= (18 - |beta|)/2.  So:

  1. BZ census of the Y2 class-0x36 coset at <= 18 (n = 72, kappa = 32);
  2. per censused beta, the bounded-overflow lift fiber (a32 MITM
     enumerators, deep lane for caps 5-6, gate-validated);
  3. filter lifts to class 0x6 (consistency: any lift with class in SEAM
     must be exactly 0x6 — p1*-injectivity, asserted);
  4. the resulting element set must EQUAL the direct (Part 4a) census
     EXACTLY (sheet-flip expansion via the deck translate (0,2)).

Also: 3-orbit descent (class 0x73 -> Y2-class 0xf3 coset) must produce
zero elements <= 18 — the direct emptiness re-derived from n = 72 data.

Appendix: the empty-window coset-base weights of the A32 GB censuses
(the latent census() edge found while porting — verified harmless there).

Output: data/a33/h5_descent.json
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
from a33_tower_cells import (  # noqa: E402
    build_tower, h1_map, rep_for, seam_data,
)
from a33_rung_cell import YRungCell  # noqa: E402

MAIN = Path("/Users/stavanjain/Code/qec-lab/experiments/bb_lab/data")
DATA = LAB / "data" / "a33"


def census_coset(binp, tag, HX, t0v, W, rpair, n, deadline):
    """All elements of t0 + rowspace(HX) with 0 < |c| <= W (both windows,
    deduped; empty-window bases handled explicitly)."""
    I1, G1, I2, G2, kappa = disjoint_info_sets(HX)
    hits: set[int] = set()
    empties = []
    nodes = 0
    for wi, (window, Gs, r) in enumerate(
            [(I1, G1, rpair[0]), (I2, G2, rpair[1])]):
        cb = coset_base(Gs, window, t0v)
        wcb = int(cb.sum())
        empties.append(wcb)
        if 0 < wcb <= W:
            hits.add(v2i(cb))
        res = run_window(binp, f"{tag}_w{wi}", Gs, [cb], r, W, deadline)
        nodes += res["nodes"]
        for j, hx in res.pop("hit_rows"):
            v = unpack3(hx, n)
            if v.any():
                hits.add(v2i(v))
    return hits, kappa, nodes, empties


def main():
    t0 = time.monotonic()
    out: dict = {}
    Y2, Y4, Y8, deck_top, deck_bot = build_tower()
    M2 = h1_map(deck_top)
    M1 = h1_map(deck_bot)
    sd = seam_data(Y4, M2)
    rep12_cls, rep3_cls = sd["reps"]
    push12 = TS._apply(M1, rep12_cls)
    push3 = TS._apply(M1, rep3_cls)
    print(f"[{time.monotonic()-t0:5.1f}s] descent targets: class "
          f"{rep12_cls:#x} -> Y2 class {push12:#x}; {rep3_cls:#x} -> "
          f"{push3:#x}; dead branches: SEAM cap K1 = 0 (measured)")

    # ------------------------------------------- validation gate (deep MITM)
    # deep enumerator == shallow lane at cap <= 4, on two rung-1 fibers
    stab_rows = [Y2.HX[0], (Y2.HX[3] + Y2.HX[17]) % 2]
    for sv in stab_rows:
        l1 = enumerate_lifts(deck_bot, sv, cap=3)
        l2 = enumerate_lifts_deep(deck_bot, sv, cap=3)
        assert l1 == l2, "deep enumerator mismatch on rung-1 fiber"
    print(f"[{time.monotonic()-t0:5.1f}s] validation gate: deep MITM == "
          f"size-4 lane on 2 rung-1 fibers  [EXACT]")

    # ------------------------------------------- Y2 coset censuses (n = 72)
    binp = build_kernel()
    b12 = rep_for(Y2, push12)
    b3 = rep_for(Y2, push3)
    deadline = time.monotonic() + 1200
    els12, kappa2, nodes12, emp12 = census_coset(
        binp, "a33_desc12", Y2.HX, b12, 18, (9, 8), Y2.n, deadline)
    els3, _, nodes3, emp3 = census_coset(
        binp, "a33_desc3", Y2.HX, b3, 18, (9, 8), Y2.n, deadline)
    assert kappa2 == 32
    wh12: dict[int, int] = {}
    for h in els12:
        w = bin(h).count("1")
        wh12[w] = wh12.get(w, 0) + 1
    wh3: dict[int, int] = {}
    for h in els3:
        w = bin(h).count("1")
        wh3[w] = wh3.get(w, 0) + 1
    print(f"[{time.monotonic()-t0:5.1f}s] Y2 censuses: class {push12:#x} "
          f"coset {len(els12)} elements <= 18 {dict(sorted(wh12.items()))}; "
          f"class {push3:#x} coset {len(els3)} elements "
          f"{dict(sorted(wh3.items()))} [empty-window weights "
          f"{emp12}/{emp3}]")
    out["y2_censuses"] = {
        "class12": f"{push12:#x}", "n12": len(els12),
        "whist12": {str(k): v for k, v in sorted(wh12.items())},
        "class3": f"{push3:#x}", "n3": len(els3),
        "whist3": {str(k): v for k, v in sorted(wh3.items())},
        "nodes": nodes12 + nodes3,
    }

    # ------------------------------------------- fibers over each shadow
    Sbb, Sbp = sd["basis"], sd["piv"]
    sigma1 = TS._perm_for(Y4, (0, 2))  # rung-1 deck translate on C1(Y4)

    def run_descent(els: set[int], want_cls: int):
        found: set[int] = set()
        n_lifts = 0
        n_seam = 0
        fiber_stats = {"fibers": 0, "vacuous": 0}
        for h in sorted(els):
            beta = i2v(h, Y2.n)
            wb = int(beta.sum())
            cap = (18 - wb) // 2
            lifts = (enumerate_lifts(deck_bot, beta, cap=cap) if cap <= 4
                     else enumerate_lifts_deep(deck_bot, beta, cap=cap))
            fiber_stats["fibers"] += 1
            if not lifts:
                fiber_stats["vacuous"] += 1
            for v0c, m2 in lifts.items():
                for v0i in (v0c, v0c ^ h):  # both sheet choices
                    w_el = deck_bot.lift(i2v(v0i, Y2.n), beta)
                    ww = int(w_el.sum())
                    assert ww == wb + 2 * m2 <= 18
                    assert Y4.is_cycle(w_el)
                    n_lifts += 1
                    s = v2i(Y4.sig(w_el))
                    if TS.in_span(s, Sbb, Sbp) and s != 0:
                        assert s == want_cls, \
                            "SEAM lift in wrong class (p1* injectivity?!)"
                        n_seam += 1
                        found.add(v2i(w_el))
        return found, n_lifts, n_seam, fiber_stats

    found12, nl12, ns12, fs12 = run_descent(els12, rep12_cls)
    print(f"[{time.monotonic()-t0:5.1f}s] 12-orbit descent: "
          f"{fs12['fibers']} fibers ({fs12['vacuous']} carry-infeasible), "
          f"{nl12} lifts, {ns12} in SEAM -> {len(found12)} distinct "
          f"class-{rep12_cls:#x} elements")
    found3, nl3, ns3, fs3 = run_descent(els3, rep3_cls)
    print(f"[{time.monotonic()-t0:5.1f}s] 3-orbit descent: "
          f"{fs3['fibers']} fibers ({fs3['vacuous']} carry-infeasible), "
          f"{nl3} lifts, {ns3} in SEAM -> {len(found3)} elements "
          f"(expect 0)")
    assert len(found3) == 0, "3-orbit descent produced elements?!"

    # ------------------------------------------- equality with the direct
    direct = set()
    for line in (DATA / "seam_census.jsonl").open():
        r = json.loads(line)
        v = np.zeros(Y4.n, dtype=np.uint8)
        v[r["w_support"]] = 1
        direct.add(v2i(v))
    assert len(direct) == 1680
    assert found12 == direct, (
        f"descent set != direct set: {len(found12)} vs {len(direct)}, "
        f"diff {len(found12 ^ direct)}")
    print(f"[{time.monotonic()-t0:5.1f}s] DESCENT == DIRECT: the 1,680 "
          f"elements re-derived exactly from n = 72 censuses "
          f"(=> identical rung verdicts; 4a's 1,680/1,680 PASS applies)")
    out["descent"] = {
        "n12_fibers": fs12["fibers"], "n12_vacuous": fs12["vacuous"],
        "n12_lifts": nl12, "n12_seam": ns12,
        "equals_direct": True, "elements": 1680,
        "n3_fibers": fs3["fibers"], "n3_vacuous": fs3["vacuous"],
        "n3_lifts": nl3, "n3_elements": 0,
    }

    # belt-and-suspenders: re-rung a sample of the descent list
    cell_top = YRungCell("top", Y4, Y8, deck_top)
    rng = np.random.default_rng(7)
    sample = rng.choice(sorted(found12), size=50, replace=False)
    for h in sample:
        v = i2v(int(h), Y4.n)
        M = (21 - int(v.sum())) // 2
        r = cell_top.seam_rung(v, M)
        assert r["verdict"] == "PASS"
    print(f"[{time.monotonic()-t0:5.1f}s] 50-element re-rung sample: "
          f"50/50 PASS")

    # ------------------------------------------- appendix: a32 census edge
    # The a32_gb_census.census() helper never visits the empty-window
    # coset element c_empty_j (the BZ kernel XORs nonempty generator
    # combinations only).  Found while porting; my census_coset above
    # handles it (and the edge was LIVE here: the Y2 class-0x36 window-1
    # base has weight 18 <= W).  Audit every a32 GB census species:
    # an element is truly missed only if BOTH windows miss it
    # (empty in one window AND restriction > r in the other).
    GB = TS.BBCode("GB", (15, 3), "x^9 + y + y^2", "1 + x^10 + x^11")
    I1g, G1g, I2g, G2g, _ = disjoint_info_sets(GB.HX)
    S = np.array([GB.sig(r) for r in GB.xreps], dtype=np.uint8)
    SinvT = TS._gf2_inv(S.T)
    C = TS.BBCode("C", (30, 6), "x^9 + y + y^2", "y^3 + x^25 + x^26")
    BY = TS.BBCode("BY", (30, 3), "x^9 + y + y^2", "1 + x^25 + x^26")
    deck_y32 = TS.Deck(C, BY, lambda e: (e[0], e[1] % 3),
                       lambda e, s: (e[0], e[1] + 3 * s))
    deck_x32 = TS.Deck(BY, GB, lambda e: (e[0] % 15, e[1]),
                       lambda e, s: (e[0] + 15 * s, e[1]))
    My32 = h1_map(deck_y32)
    Mx32 = h1_map(deck_x32)
    Ry = TS._colspace(My32)
    from a30_rung_pass import rref_ints  # noqa: E402
    Ryb, Ryp = rref_ints(list(Ry))
    W_ints = sorted({TS._apply(Mx32, s)
                     for s in TS._span_points(Ryb, Ryp)} - {0})

    def gb_rep(si: int) -> np.ndarray:
        tvec = i2v(si, 8)
        coeff = (SinvT @ tvec) % 2
        v = np.zeros(GB.n, dtype=np.uint8)
        for i in range(8):
            if coeff[i]:
                v ^= GB.xreps[i]
        return v

    windows_g = [(I1g, G1g), (I2g, G2g)]
    perms_gb = TS._translation_perms(GB)
    audit = {}
    missed: list[tuple[str, np.ndarray]] = []
    # (species, W, r-pair, class list)
    species = [("logical<=10", 10, (5, 4), list(range(1, 256))),
               ("wcoset<=16", 16, (8, 7), W_ints),
               ("wcoset<=22", 22, (11, 10), W_ints)]
    for name, Wg, rp, classes in species:
        risky = 0
        truly = 0
        for si in classes:
            v = gb_rep(si)
            for wi, (window, Gs) in enumerate(windows_g):
                ce = coset_base(Gs, window, v)
                if not 0 < int(ce.sum()) <= Wg:
                    continue
                risky += 1
                other = windows_g[1 - wi][0]
                r_other = rp[1 - wi]
                if int(ce[other].sum()) > r_other:
                    truly += 1
                    missed.append((name, ce))
        audit[name] = {"empty_in_range": risky, "truly_missed": truly}
    print(f"[{time.monotonic()-t0:5.1f}s] appendix: a32 census() "
          f"empty-window audit: {audit}")
    # logical <= 10: the 15 missed vectors are all weight-10; their
    # translation ORBITS must all be in the stored 38-orbit census
    # (then every per-orbit consumer — sector B, dby floor, d(GB) = 8 —
    # is unaffected; erratum is at raw-vector-count level only)
    stored_canons = {json.loads(x)["canon"] for x in
                     (DATA.parent / "a32" / "gb_census_logical.jsonl"
                      ).open()}
    n_log_missed = 0
    for name, ce in missed:
        if name != "logical<=10":
            continue
        assert int(ce.sum()) == 10
        assert TS._canon(ce, perms_gb) in stored_canons, \
            "missed logical in a NEW orbit — a32 sector B has a hole!"
        n_log_missed += 1
    print(f"    logical<=10: {n_log_missed} missed vectors, ALL in "
          f"stored orbits => per-orbit consumers (sector B, dby, "
          f"d(GB)=8) UNAFFECTED; erratum = vector counts only "
          f"(1,623 -> 1,638)")
    # wcoset <= 22: for each missed element, (a) was its orbit
    # enumerable via some translate?  (b) patch regardless: the only
    # band-22 consumer is the sector-A cap-0 fiber (flat lifts, M = 1)
    patch = []
    for name, ce in missed:
        if name != "wcoset<=22":
            continue
        assert int(ce.sum()) == 22
        orbit_found = False
        for t in GB.G:
            tv = ce[TS._perm_for(GB, t)]
            if int(tv[I1g].sum()) <= 11 or int(tv[I2g].sum()) <= 10:
                orbit_found = True
                break
        lifts = enumerate_lifts(deck_x32, ce, cap=0)
        entry = {"species": name, "w": 22, "orbit_enumerable": orbit_found,
                 "flat_lifts": 0, "rungs": []}
        if lifts:
            from scope_bravyi_rung import BravyiRungCell  # noqa: E402
            cell32 = BravyiRungCell()
            for v0c, m2 in lifts.items():
                b = deck_x32.lift(i2v(v0c, GB.n), ce)
                assert BY.is_cycle(b) and int(b.sum()) == 22
                entry["flat_lifts"] += 1
                res = cell32.rung(b, 1, time.monotonic() + 600)
                assert res["verdict"] == "PASS", res
                entry["rungs"].append(res["verdict"])
        patch.append(entry)
        print(f"    missed wcoset<=22 element: orbit enumerable via "
              f"translate = {orbit_found}; {entry['flat_lifts']} flat "
              f"lifts, rungs {entry['rungs'] or 'none needed'}")
    out["a32_census_edge"] = {
        "audit": audit,
        "logical_missed_all_in_stored_orbits": n_log_missed,
        "wcoset22_patch": patch,
        "verdict": "A32 d=24 assembly UNAFFECTED (per-orbit consumers "
                   "covered; missed band-22 elements patched directly)",
    }

    out["wall_s"] = round(time.monotonic() - t0, 1)
    (DATA / "h5_descent.json").write_text(json.dumps(out, indent=1))
    print(f"total {out['wall_s']}s -> {DATA / 'h5_descent.json'}")


if __name__ == "__main__":
    main()
