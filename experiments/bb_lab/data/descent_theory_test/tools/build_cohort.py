"""Phase 0 — stratified cohort construction for the descent-theory test.

Strata (see PROTOCOL.md for the frozen definitions):
  frontier      no d_exact, d_ub >= 12, n in [96,168], even |G|
  anchor        corpus-exact d in {6,8,10,12} across 2-adic depths 1..4
  scope-control odd |G| (prediction: machinery must REFUSE)
  parity-strip  freshly sampled codes with an even-weight polynomial

Selection is deterministic (fixed seeds / lexicographic tie-breaks).
Every screened candidate is logged to prescreen_pool.jsonl; the selected
cohort goes to cohort.jsonl; per-row structural screens are cached in
screens_cache.jsonl for Phase 1 to consume unchanged.

NO distance work happens here beyond (i) reading corpus values and
(ii) L1 d_ub sampling for the fresh parity-strip rows (witnesses only).
"""

from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import duckdb
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import dtt_lib as D                                    # noqa: E402
from dtt_lib import OUT                                # noqa: E402

T0 = time.monotonic()
BUDGET_S = 90 * 60
MAX_SCREENED = 1000

FRONTIER_QUOTA = {          # group -> (quota, v2 depth)
    "Z9xZ6": (20, 1), "Z6xZ13": (12, 1),
    "Z6xZ10": (8, 2), "Z5xZ12": (8, 2), "Z6xZ14": (10, 2),
    "Z7xZ12": (10, 2),
    "Z12xZ6": (14, 3), "Z8xZ9": (14, 3),
}
OVERDRAW = 1.6              # screen ~1.6x quota, select the quota

# (group, d_exact, note) anchor cells; one row per cell, min(instance_id),
# except Z12xZ6/d=12 which prefers the gross code itself when present.
ANCHOR_CELLS = [
    ("Z6xZ7", 10, "depth1"), ("Z9xZ6", 12, "depth1"),
    ("Z6xZ13", 8, "depth1"), ("Z15xZ6", 10, "depth1 (mitten host shape)"),
    ("Z6xZ14", 12, "depth2"), ("Z6xZ10", 10, "depth2"),
    ("Z4xZ9", 8, "depth2"), ("Z6xZ6", 6, "depth2 (gross-base shape)"),
    ("Z12xZ6", 12, "depth3 (prefer gross)"), ("Z8xZ9", 10, "depth3"),
    ("Z7xZ8", 12, "depth3"), ("Z3xZ8", 6, "depth3"),
    ("Z4xZ12", 12, "depth4"), ("Z12xZ12", 12, "depth4 (n=288)"),
]

# (group, want_exact_d | None, count) for odd-|G| scope controls
SCOPE_CELLS = [
    ("Z7xZ7", 12, 1), ("Z7xZ7", 10, 1),
    ("Z7xZ11", 16, 1), ("Z7xZ11", 12, 1),
    ("Z9xZ9", 12, 1), ("Z9xZ9", None, 1),      # None = unfilled row
    ("Z7xZ9", 14, 1), ("Z7xZ9", None, 1),
    ("Z15xZ3", 10, 1), ("Z5xZ9", 12, 1),
]

PARITY_GROUPS = [(6, 10), (12, 6)]     # frontier shapes, 5 rows each
PARITY_PER_GROUP = 5
PARITY_TRY_CAP = 30_000

GROSS = ((12, 6), "x^3 + y + y^2", "y^3 + x + x^2")


def log(msg: str) -> None:
    print(f"[{time.monotonic()-T0:7.1f}s] {msg}", flush=True)


def jline(fh, obj) -> None:
    fh.write(json.dumps(obj) + "\n")
    fh.flush()


def screen_summary(s: dict) -> dict:
    return {
        "depth_used": s["depth_used"],
        "depth_available": s["depth_available"],
        "k_chain": [lv["k"] for lv in s["levels"]],
        "regimes": s["regimes"],
        "rank_law_all": bool(all(r.get("rank_law_holds")
                                 for r in s["rungs"])),
        "truncated": bool(s["truncated"]),
    }


def main() -> None:
    con = duckdb.connect(str(D.MAIN_DB), read_only=True)
    pool_f = (OUT / "prescreen_pool.jsonl").open("w")
    cache_f = (OUT / "screens_cache.jsonl").open("w")
    cohort_f = (OUT / "cohort.jsonl").open("w")
    n_screened = 0
    selection_log: list[dict] = []

    def budget_ok() -> bool:
        return (time.monotonic() - T0) < BUDGET_S \
            and n_screened < MAX_SCREENED

    def screen_and_log(row: dict, stratum: str) -> dict | None:
        """Structural screen + pool/cache logging. Returns screen dict."""
        nonlocal n_screened
        n_screened += 1
        try:
            s = D.screen_structure(row["ell"], row["m"], row["A_poly"],
                                   row["B_poly"], row["instance_id"][:8])
        except Exception as e:
            jline(pool_f, {**row, "stratum": stratum,
                           "screen_error": f"{type(e).__name__}: {e}"})
            return None
        jline(pool_f, {**row, "stratum": stratum,
                       "screen": screen_summary(s)})
        jline(cache_f, {"instance_id": row["instance_id"], "screen": s})
        return s

    # ---------------------------------------------------------- frontier
    log("[frontier] drawing candidates (deterministic d_ub round-robin)")
    frontier_selected: list[dict] = []
    frontier_screens: dict[str, dict] = {}
    for grp, (quota, depth) in FRONTIER_QUOTA.items():
        rows = con.execute(
            """select instance_id, code_id, group_struct, ell, m, n, k,
                      A_poly, B_poly, d_lb, d_ub, d_exact
               from bb_instances
               where group_struct = ? and d_exact is null and d_ub >= 12
                     and n between 96 and 168
               order by d_ub, instance_id""", [grp]).fetchall()
        cols = ["instance_id", "code_id", "group_struct", "ell", "m", "n",
                "k", "A_poly", "B_poly", "d_lb", "d_ub", "d_exact"]
        rows = [dict(zip(cols, r)) for r in rows]
        # round-robin across d_ub bands (spread), inst_id order within band
        bands: dict[int, list[dict]] = {}
        for r in rows:
            bands.setdefault(r["d_ub"], []).append(r)
        order: list[dict] = []
        idx = 0
        while len(order) < len(rows):
            for du in sorted(bands):
                if idx < len(bands[du]):
                    order.append(bands[du][idx])
            idx += 1
        n_cand = min(len(order), int(quota * OVERDRAW) + 1)
        cands = order[:n_cand]
        screened = []
        for r in cands:
            if not budget_ok():
                log("  BUDGET HIT during frontier prescreen")
                break
            s = screen_and_log(r, "frontier")
            if s is not None:
                screened.append((r, s))
        # selection: prefer regime diversity, k=6 rows, d_ub spread
        picked: list[tuple[dict, dict]] = []
        seen_reg: set[str] = set()
        seen_dub: set[int] = set()

        def keyfun(t):
            r, s = t
            new_reg = len(set(s["regimes"]) - seen_reg)
            return (-new_reg, -r["k"],
                    0 if r["d_ub"] not in seen_dub else 1,
                    r["instance_id"])
        rem = list(screened)
        while rem and len(picked) < quota:
            rem.sort(key=keyfun)
            r, s = rem.pop(0)
            picked.append((r, s))
            seen_reg |= set(s["regimes"])
            seen_dub.add(r["d_ub"])
        for r, s in picked:
            r2 = {**r, "stratum": "frontier", "v2_depth": depth,
                  "d_ub_provenance": "corpus L1-sampling (a18_fill_ubs)"}
            frontier_selected.append(r2)
            frontier_screens[r["instance_id"]] = s
        log(f"  {grp}: pool {len(rows)}, screened {len(screened)}, "
            f"selected {len(picked)}/{quota}")
        selection_log.append({"stratum": "frontier", "group": grp,
                              "pool": len(rows),
                              "screened": len(screened),
                              "selected": len(picked)})

    # ----------------------------------------------------------- anchors
    log("[anchors] fixed (group, d) cells")
    anchor_selected: list[dict] = []
    gross_iid = D.canonical_id(GROSS[0], D.poly_terms(GROSS[1], GROSS[0]),
                               D.poly_terms(GROSS[2], GROSS[0]))[0]
    for grp, d, note in ANCHOR_CELLS:
        q = """select instance_id, code_id, group_struct, ell, m, n, k,
                      A_poly, B_poly, d_lb, d_ub, d_exact, d_method
               from bb_instances where group_struct = ? and d_exact = ?
               order by (instance_id != ?) , instance_id limit 1"""
        prefer = gross_iid if "gross" in note and grp == "Z12xZ6" else ""
        row = con.execute(q, [grp, d, prefer]).fetchone()
        if row is None:
            log(f"  cell ({grp}, d={d}): EMPTY in corpus -- skipped")
            selection_log.append({"stratum": "anchor", "group": grp,
                                  "d": d, "selected": 0,
                                  "note": "cell empty in corpus"})
            continue
        cols = ["instance_id", "code_id", "group_struct", "ell", "m", "n",
                "k", "A_poly", "B_poly", "d_lb", "d_ub", "d_exact",
                "d_method"]
        r = dict(zip(cols, row))
        if not budget_ok():
            break
        s = screen_and_log(r, "anchor")
        if s is None:
            continue
        r2 = {**r, "stratum": "anchor", "v2_depth": D.v2(r["ell"] * r["m"]),
              "anchor_note": note,
              "d_provenance": f"corpus-exact ({r['d_method']})"}
        anchor_selected.append(r2)
        frontier_screens[r["instance_id"]] = s
        log(f"  ({grp}, d={d}): {r['instance_id']} [[{r['n']},{r['k']}]] "
            f"{'<-- GROSS' if r['instance_id'] == gross_iid else ''}")
        selection_log.append({"stratum": "anchor", "group": grp, "d": d,
                              "selected": 1, "iid": r["instance_id"]})

    # ----------------------------------------------------- scope controls
    log("[scope-controls] odd-|G| rows (prediction: REFUSE)")
    scope_selected: list[dict] = []
    for grp, d, cnt in SCOPE_CELLS:
        if d is None:
            q = ("select instance_id, code_id, group_struct, ell, m, n, k,"
                 " A_poly, B_poly, d_lb, d_ub, d_exact, d_method"
                 " from bb_instances where group_struct = ? and d_exact is"
                 " null order by instance_id limit ?")
            rows = con.execute(q, [grp, cnt]).fetchall()
        else:
            q = ("select instance_id, code_id, group_struct, ell, m, n, k,"
                 " A_poly, B_poly, d_lb, d_ub, d_exact, d_method"
                 " from bb_instances where group_struct = ? and d_exact = ?"
                 " order by instance_id limit ?")
            rows = con.execute(q, [grp, d, cnt]).fetchall()
        cols = ["instance_id", "code_id", "group_struct", "ell", "m", "n",
                "k", "A_poly", "B_poly", "d_lb", "d_ub", "d_exact",
                "d_method"]
        for row in rows:
            r = dict(zip(cols, row))
            assert (r["ell"] * r["m"]) % 2 == 1
            r2 = {**r, "stratum": "scope-control",
                  "v2_depth": 0,
                  "prediction_basis": "odd |G|: no free Z2 deck exists "
                                      "(A35 C1); machinery must refuse"}
            scope_selected.append(r2)
            jline(pool_f, {**r, "stratum": "scope-control",
                           "screen": {"refused": "odd |G|"}})
        selection_log.append({"stratum": "scope-control", "group": grp,
                              "d": d, "selected": len(rows)})
    log(f"  {len(scope_selected)} scope-control rows")

    # ------------------------------------------------------- parity strip
    log("[parity-strip] sampling wt-3 x wt-4 codes (k >= 2)")
    parity_selected: list[dict] = []
    rng = np.random.default_rng(D.L1_SEED)
    for (pl, pm) in PARITY_GROUPS:
        G = D.AbelianGroup((pl, pm))
        N = G.cardinality
        label = G.label()
        got, tried, seen = 0, 0, set()
        while got < PARITY_PER_GROUP and tried < PARITY_TRY_CAP:
            tried += 1
            # A weight 3 (odd), B weight 4 (even) -- "at least one even"
            elemsA = rng.choice(N, size=3, replace=False)
            elemsB = rng.choice(N, size=4, replace=False)
            tA = frozenset(G.from_index(int(i)) for i in elemsA)
            tB = frozenset(G.from_index(int(i)) for i in elemsB)
            k = D.bbcode_k((pl, pm), tA, tB)
            if k < 2:
                continue
            iid, lab, cA, cB = D.canonical_id((pl, pm), tA, tB)
            if iid in seen:
                continue
            seen.add(iid)
            # fresh L1 d_ub (witness only; allowed rule c)
            A = D.Poly(support=tA, group=G)
            B = D.Poly(support=tB, group=G)
            res = D.l1_distance_ub(D.bb_check_matrices(A, B),
                                   n_samples=D.L1_SAMPLES_TARGET,
                                   seed=D.L1_SEED)
            r = {"instance_id": iid,
                 "code_id": f"dtt_parity_{label}_{iid[:8]}",
                 "group_struct": label, "ell": pl, "m": pm,
                 "n": 2 * N, "k": int(k), "A_poly": cA, "B_poly": cB,
                 "d_lb": None, "d_ub": int(res.distance_ub),
                 "d_exact": None}
            s = screen_and_log(r, "parity-strip")
            if s is None:
                continue
            r2 = {**r, "stratum": "parity-strip",
                  "v2_depth": D.v2(pl * pm),
                  "d_ub_provenance": f"fresh L1 sampling "
                                     f"({D.L1_SAMPLES_TARGET} samples, "
                                     f"seed {D.L1_SEED}; witness only)",
                  "weights": [3, 4]}
            parity_selected.append(r2)
            frontier_screens[iid] = s
            got += 1
        log(f"  {label}: {got}/{PARITY_PER_GROUP} rows in {tried} tries")
        selection_log.append({"stratum": "parity-strip", "group": label,
                              "selected": got, "tried": tried})

    # ------------------------------------------- regime coverage + top-up
    log("[regimes] occupancy check over in-scope cohort rows")
    def occupancy() -> dict[str, int]:
        occ: dict[str, int] = {"R1": 0, "R2": 0, "R3": 0, "R4": 0}
        for iid, s in frontier_screens.items():
            for reg in s["regimes"]:
                occ[reg] += 1
        return occ
    occ = occupancy()
    log(f"  initial occupancy (rows containing each regime): {occ}")
    missing = [reg for reg, c in occ.items() if c == 0]
    topup_log = []
    if missing:
        log(f"  missing regimes {missing}: top-up sweep over deeper pool")
        # sweep additional multi-rung candidates (depth >= 2 groups),
        # deterministic order, until filled or pool/budget exhausted
        rows = con.execute(
            """select instance_id, code_id, group_struct, ell, m, n, k,
                      A_poly, B_poly, d_lb, d_ub, d_exact
               from bb_instances
               where d_exact is null and d_ub >= 12
                     and n between 96 and 168
                     and group_struct in
                     ('Z5xZ12','Z6xZ10','Z6xZ14','Z7xZ12','Z12xZ6','Z8xZ9')
               order by instance_id""").fetchall()
        cols = ["instance_id", "code_id", "group_struct", "ell", "m", "n",
                "k", "A_poly", "B_poly", "d_lb", "d_ub", "d_exact"]
        pool2 = [dict(zip(cols, r)) for r in rows
                 if r[0] not in frontier_screens]
        for r in pool2:
            if not missing or not budget_ok():
                break
            s = screen_and_log(r, "frontier-topup")
            if s is None:
                continue
            hit = set(s["regimes"]) & set(missing)
            if hit:
                r2 = {**r, "stratum": "frontier",
                      "v2_depth": D.v2(r["ell"] * r["m"]),
                      "d_ub_provenance":
                          "corpus L1-sampling (a18_fill_ubs)",
                      "topup_for_regime": sorted(hit)}
                frontier_selected.append(r2)
                frontier_screens[r["instance_id"]] = s
                missing = [x for x in missing if x not in hit]
                topup_log.append({"iid": r["instance_id"],
                                  "filled": sorted(hit)})
                log(f"  top-up {r['instance_id']} fills {sorted(hit)}")
        if missing:
            log(f"  UNFILLABLE from pool within budget: {missing}")
    selection_log.append({"stratum": "regime-topup", "topups": topup_log,
                          "unfillable": missing,
                          "final_occupancy": occupancy()})

    # ------------------------------------------------------------- freeze
    cohort = (frontier_selected + anchor_selected + scope_selected
              + parity_selected)
    for r in cohort:
        jline(cohort_f, r)
    for fh in (pool_f, cache_f, cohort_f):
        fh.close()
    (OUT / "selection_log.json").write_text(
        json.dumps({"selection": selection_log,
                    "n_screened": n_screened,
                    "wall_s": round(time.monotonic() - T0, 1)}, indent=1))
    by_stratum: dict[str, int] = {}
    for r in cohort:
        by_stratum[r["stratum"]] = by_stratum.get(r["stratum"], 0) + 1
    log(f"cohort frozen: {len(cohort)} rows {by_stratum}; "
        f"screened {n_screened}; occupancy {occupancy()}")


if __name__ == "__main__":
    main()
