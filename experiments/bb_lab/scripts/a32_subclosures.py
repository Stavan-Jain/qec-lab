"""A32 Part 6: live sub-closures of the SF24-y strata via the tower.

6a  SECTOR B (diagonal, beta = 0): the only live contribution is
    |gamma| = 10 (parity kills odd |gamma|; |gamma| = 8 gives band 16,
    already closed by A24 but re-run here at M = 4 for slack + engine
    consistency).  For each of the 36 weight-10 GB-logical orbit reps
    with class outside im p_x* (= all of them, Part 5a):
    b = tau'(gamma) (weight 20), top rung at M = 2.
    ALL PASS  =>  the diagonal sub-stratum of (20,1) is CLOSED.

6b  SECTOR A pilot (|beta| = 14): for each of the 6 weight-14 W-coset
    orbit reps: enumerate ALL x-deck lifts b with m2 <= 4 (restricted
    MITM lane, complete), then per lift run the top rung at
    M = 5 - m2.  ALL PASS  =>  sector A at |beta| = 14 contributes
    nothing below 24 (strata (16,3), (18,2), (20,1), (22,0) closed on
    this fiber family).
    Falsify-first cross-checks baked in:
      - m2 = 0 lifts must NOT exist (A24 band-14: no reachable-class
        weight-14 BY logical);
      - m2 = 1 lifts must reproduce EXACTLY the 3 banked band-16
        A-type classes (A24 cell_census_reach_band16).

6c  SECTOR A |beta| = 16 (the 68 W-coset orbit reps at 16, cap m2 <= 3,
    M = 4 - m2): same engine, closing sector A at |beta| = 16.

Transport: rung verdicts are constant on GB-translation orbits of beta
(translations lift through both decks; spot-verified below), and on the
sheet flip v0 -> v0 + beta (the deck translate sigma'), which halves the
lift lists.

Output: data/a32/subclosure_{B,A14,A16}.jsonl + summary.
"""

from __future__ import annotations

import itertools
import json
import sys
import time
from pathlib import Path

import numpy as np

LAB = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(LAB / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from a30_rung_pass import i2v, reduce_int, rref_ints, v2i  # noqa: E402
import a32_tower_slice as TS  # noqa: E402
from scope_bravyi_rung import BravyiRungCell  # noqa: E402

MAIN = Path("/Users/stavanjain/Code/qec-lab/experiments/bb_lab/data")
DATA = LAB / "data" / "a32"
DATA.mkdir(parents=True, exist_ok=True)


def enumerate_lifts(deck: TS.Deck, beta: np.ndarray, cap: int,
                    kernel_cap: int = 16) -> dict[int, int]:
    """All v0 with E v0 = RHS beta and |v0 off supp(beta)| <= cap.

    Returns {canonical v0 int -> m2}, canonical under v0 -> v0 + beta
    (the deck translate).  Complete by the exact-off-support subset-sum
    argument: a solution with off-support X0 forces
    sum_{j in X0} red[j] = rhs_res, |X0| <= cap; every such X is tried
    and the full kernel over supp(beta) u X enumerated.
    """
    assert cap <= 4, "MITM lane implemented to size 4"
    n = deck.base.n
    E_cols = [v2i(deck.E[:, j]) for j in range(n)]
    rhs = (deck.RHS @ beta) % 2
    rhs_i = v2i(rhs)
    bsupp = [int(j) for j in np.nonzero(beta)[0]]
    bmask = v2i(beta)
    # reduce off-support columns modulo the on-support span
    bcols = [E_cols[j] for j in bsupp]
    bb, bp = rref_ints(bcols)
    rhs_res = reduce_int(rhs_i, bb, bp)
    offb = [j for j in range(n) if not (bmask >> j) & 1]
    red = {j: reduce_int(E_cols[j], bb, bp) for j in offb}
    by_val: dict[int, list[int]] = {}
    for j in offb:
        by_val.setdefault(red[j], []).append(j)
    hits_X: set[tuple[int, ...]] = set()
    if rhs_res == 0:
        hits_X.add(())
    if cap >= 1:
        for j in by_val.get(rhs_res, []):
            hits_X.add((j,))
    if cap >= 2:
        for j1 in offb:
            for j2 in by_val.get(rhs_res ^ red[j1], []):
                if j2 > j1:
                    hits_X.add((j1, j2))
    if cap >= 3:
        for j1, j2 in itertools.combinations(offb, 2):
            for j3 in by_val.get(rhs_res ^ red[j1] ^ red[j2], []):
                if j3 > j2:
                    hits_X.add((j1, j2, j3))
    if cap >= 4:
        pair_sum: dict[int, list[tuple[int, int]]] = {}
        for j1, j2 in itertools.combinations(offb, 2):
            pair_sum.setdefault(red[j1] ^ red[j2], []).append((j1, j2))
        for val, prs in pair_sum.items():
            for j3, j4 in pair_sum.get(rhs_res ^ val, []):
                for j1, j2 in prs:
                    if j2 < j3:
                        hits_X.add((j1, j2, j3, j4))
    out: dict[int, int] = {}
    for X in sorted(hits_X):
        cols = bsupp + list(X)
        # rref with combination tracking
        b3: list[int] = []
        p3: list[int] = []
        h3: list[int] = []
        deps: list[int] = []
        for ci, j in enumerate(cols):
            cur, h = E_cols[j], 1 << ci
            for bb3, pp3, hh in zip(b3, p3, h3):
                if (cur >> pp3) & 1:
                    cur ^= bb3
                    h ^= hh
            if cur:
                b3.append(cur)
                p3.append((cur & -cur).bit_length() - 1)
                h3.append(h)
            else:
                deps.append(h)
        cur, hsel = rhs_i, 0
        for bb3, pp3, hh in zip(b3, p3, h3):
            if (cur >> pp3) & 1:
                cur ^= bb3
                hsel ^= hh
        if cur:
            continue
        assert len(deps) <= kernel_cap, f"kernel 2^{len(deps)} at X={X}"
        for kt in range(1 << len(deps)):
            sel = hsel
            for jj in range(len(deps)):
                if (kt >> jj) & 1:
                    sel ^= deps[jj]
            v0_int = 0
            for ci, j in enumerate(cols):
                if (sel >> ci) & 1:
                    v0_int |= 1 << j
            m2 = bin(v0_int & ~bmask).count("1")
            if m2 > cap:
                continue
            # verify solution
            v0 = i2v(v0_int, n)
            assert not (((deck.E @ v0) + rhs) % 2).any(), "not a solution"
            canon = min(v0_int, v0_int ^ bmask)
            prev = out.get(canon)
            if prev is None or m2 < prev:
                out[canon] = m2
    return out


def main():
    t0 = time.monotonic()
    out: dict = {}
    GB = TS.BBCode("GB", (15, 3), "x^9 + y + y^2", "1 + x^10 + x^11")
    BY = TS.BBCode("BY", (30, 3), "x^9 + y + y^2", "1 + x^25 + x^26")
    deck_x = TS.Deck(BY, GB, lambda e: (e[0] % 15, e[1]),
                     lambda e, s: (e[0] + 15 * s, e[1]))
    perms_by = TS._translation_perms(BY)
    cell = BravyiRungCell()
    print(f"[{time.monotonic()-t0:5.1f}s] frames + top-rung cell built")

    # --------------------------------------------------------- 6a sector B
    rows_b = []
    log_recs = [json.loads(x) for x in
                (DATA / "gb_census_logical.jsonl").open()]
    reps: dict[int, dict] = {}
    for r in log_recs:
        if not r["in_im_px"]:
            reps.setdefault(r["canon"], r)
    assert len(reps) == 38  # 2 at w8 + 36 at w10
    viol = 0
    verd_hist: dict[str, int] = {}
    diag16_canons = set()
    for r in sorted(reps.values(), key=lambda r: (r["w"], r["canon"])):
        gam = np.zeros(GB.n, dtype=np.uint8)
        gam[r["support"]] = 1
        b = (deck_x.TAU @ gam) % 2
        assert BY.is_cycle(b) and not BY.is_stab(b), \
            "tau-lift not a BY logical (class filter failed?)"
        wb = int(b.sum())
        assert wb == 2 * r["w"]
        M = (24 - wb + 1) // 2
        res = cell.rung(b, M, time.monotonic() + 600)
        assert res["verdict"] != "ABORT", res
        assert res["verdict"] != "PASS" or res.get("lane") != "all-trivial"
        verd_hist[res["verdict"]] = verd_hist.get(res["verdict"], 0) + 1
        if res["verdict"] == "VIOLATION":
            viol += 1
            print("  !! VIOLATION (sector B):", json.dumps(res)[:300])
        if r["w"] == 8:
            diag16_canons.add(TS._canon(b, perms_by))
        rows_b.append({"w_gamma": r["w"], "gamma_canon": r["canon"],
                       "wb": wb, "M": M, "verdict": res["verdict"],
                       "lane": res.get("lane")})
    # cross-check: the two w8 diagonals are exactly A24's band-16 B-classes
    b16 = [json.loads(x) for x in
           (MAIN / "a24" / "cell_census_reach_band16_mchecks.jsonl").open()]
    banked_diag = set()
    for e in b16:
        b = np.zeros(BY.n, dtype=np.uint8)
        b[e["b_support"]] = 1
        beta = (deck_x.P @ b) % 2
        if not beta.any():
            banked_diag.add(TS._canon(b, perms_by))
    assert banked_diag == diag16_canons, \
        "w8 diagonal orbits != A24 band-16 B-classes"
    print(f"[{time.monotonic()-t0:5.1f}s] 6a sector B: 38 rungs "
          f"{verd_hist}, violations={viol}; w8-diagonals == A24 band-16 "
          f"B-classes [cross-check PASS]")
    out["sectorB"] = {"rungs": len(rows_b), "verdicts": verd_hist,
                      "violations": viol}
    with (DATA / "subclosure_B.jsonl").open("w") as f:
        for r in rows_b:
            f.write(json.dumps(r) + "\n")

    # ------------------------------------------------- 6b sector A, w14
    w_recs = [json.loads(x) for x in
              (DATA / "gb_census_wcoset.jsonl").open()]
    reps14: dict[int, dict] = {}
    for r in w_recs:
        if r["w"] == 14:
            reps14.setdefault(r["canon"], r)
    assert len(reps14) == 6
    banked_A16 = set()
    for e in b16:
        b = np.zeros(BY.n, dtype=np.uint8)
        b[e["b_support"]] = 1
        beta = (deck_x.P @ b) % 2
        if beta.any() and not GB.is_stab(beta):
            banked_A16.add(TS._canon(b, perms_by))
    assert len(banked_A16) == 3
    rows_a = []
    found_m1_lifts = set()
    viol14 = 0
    total_rungs = 0
    for rep_i, r in enumerate(sorted(reps14.values(),
                                     key=lambda r: r["canon"])):
        beta = np.zeros(GB.n, dtype=np.uint8)
        beta[r["support"]] = 1
        tL = time.monotonic()
        lifts = enumerate_lifts(deck_x, beta, cap=4)
        m2h: dict[int, int] = {}
        for m2 in lifts.values():
            m2h[m2] = m2h.get(m2, 0) + 1
        assert m2h.get(0, 0) == 0, \
            "FLAT lift of a W-14 member exists — contradicts A24 band-14!"
        # rungs per lift
        vh: dict[str, int] = {}
        for v0c, m2 in sorted(lifts.items()):
            b = deck_x.lift(i2v(v0c, GB.n), beta)
            assert BY.is_cycle(b) and not BY.is_stab(b)
            wb = int(b.sum())
            assert wb == 14 + 2 * m2
            if m2 == 1:
                found_m1_lifts.add(TS._canon(b, perms_by))
            M = (24 - wb + 1) // 2
            res = cell.rung(b, M, time.monotonic() + 600)
            assert res["verdict"] != "ABORT", res
            vh[res["verdict"]] = vh.get(res["verdict"], 0) + 1
            total_rungs += 1
            if res["verdict"] == "VIOLATION":
                viol14 += 1
                print("  !! VIOLATION (A14):", json.dumps(res)[:300])
        rows_a.append({"beta_canon": r["canon"], "m2_hist": m2h,
                       "lifts": len(lifts), "verdicts": vh,
                       "wall_s": round(time.monotonic() - tL, 2)})
        print(f"    A14 fiber {rep_i}: lifts {m2h} verdicts {vh} "
              f"({rows_a[-1]['wall_s']}s)")
    assert found_m1_lifts == banked_A16, \
        f"m2=1 lifts {len(found_m1_lifts)} orbits != banked 3 A-classes"
    print(f"[{time.monotonic()-t0:5.1f}s] 6b sector A |beta|=14: 6 fibers, "
          f"{total_rungs} rungs, violations={viol14}; m2=0 EMPTY "
          f"[A24 band-14 consistency] and m2=1 == banked band-16 A-classes "
          f"[cross-check PASS]")
    out["sectorA14"] = {"fibers": rows_a, "violations": viol14,
                        "total_rungs": total_rungs}
    with (DATA / "subclosure_A14.jsonl").open("w") as f:
        for r in rows_a:
            f.write(json.dumps(r) + "\n")

    # covariance spot-check: rung verdict transports along GB translation
    r0 = sorted(reps14.values(), key=lambda r: r["canon"])[0]
    beta = np.zeros(GB.n, dtype=np.uint8)
    beta[r0["support"]] = 1
    t = (3, 1)
    perm = TS._perm_for(GB, t)
    beta_t = beta[perm]
    lifts0 = enumerate_lifts(deck_x, beta, cap=2)
    lifts_t = enumerate_lifts(deck_x, beta_t, cap=2)
    h0 = sorted(lifts0.values())
    ht = sorted(lifts_t.values())
    assert h0 == ht, f"lift m2-profiles differ under translation {h0} {ht}"
    print(f"[{time.monotonic()-t0:5.1f}s] covariance spot-check: translated "
          f"fiber has identical m2-profile {h0}")
    out["covariance"] = {"translation": t, "m2_profile": h0}

    # ------------------------------------------------- 6c sector A, w16
    reps16: dict[int, dict] = {}
    for r in w_recs:
        if r["w"] == 16:
            reps16.setdefault(r["canon"], r)
    assert len(reps16) == 68
    rows16 = []
    viol16 = 0
    rungs16 = 0
    lifts16_h: dict[int, int] = {}
    for rep_i, r in enumerate(sorted(reps16.values(),
                                     key=lambda r: r["canon"])):
        beta = np.zeros(GB.n, dtype=np.uint8)
        beta[r["support"]] = 1
        lifts = enumerate_lifts(deck_x, beta, cap=3)
        m2h: dict[int, int] = {}
        for m2 in lifts.values():
            m2h[m2] = m2h.get(m2, 0) + 1
        for m2, c in m2h.items():
            lifts16_h[m2] = lifts16_h.get(m2, 0) + c
        vh: dict[str, int] = {}
        for v0c, m2 in sorted(lifts.items()):
            b = deck_x.lift(i2v(v0c, GB.n), beta)
            wb = int(b.sum())
            assert wb == 16 + 2 * m2
            M = (24 - wb + 1) // 2
            res = cell.rung(b, M, time.monotonic() + 600)
            assert res["verdict"] != "ABORT", res
            vh[res["verdict"]] = vh.get(res["verdict"], 0) + 1
            rungs16 += 1
            if res["verdict"] == "VIOLATION":
                viol16 += 1
                print("  !! VIOLATION (A16):", json.dumps(res)[:300])
        rows16.append({"beta_canon": r["canon"], "m2_hist": m2h,
                       "verdicts": vh})
        if rep_i % 10 == 0:
            print(f"    A16 fiber {rep_i}/68 done "
                  f"({time.monotonic()-t0:.0f}s)")
    print(f"[{time.monotonic()-t0:5.1f}s] 6c sector A |beta|=16: 68 fibers, "
          f"{rungs16} rungs, lift m2-histogram {lifts16_h}, "
          f"violations={viol16}")
    out["sectorA16"] = {"fibers": len(rows16), "total_rungs": rungs16,
                        "lift_m2_hist": lifts16_h, "violations": viol16}
    with (DATA / "subclosure_A16.jsonl").open("w") as f:
        for r in rows16:
            f.write(json.dumps(r) + "\n")

    out["wall_s"] = round(time.monotonic() - t0, 1)
    (DATA / "subclosure_summary.json").write_text(json.dumps(out, indent=1))
    print(f"total {out['wall_s']}s -> {DATA / 'subclosure_summary.json'}")


if __name__ == "__main__":
    main()
