"""Phase 1 — frozen, pre-registered predictions for every cohort row.

Consumes cohort.jsonl + screens_cache.jsonl (Phase 0, unchanged), adds
quotient-level distance info (corpus read-only / n<=40 full enumeration /
L1 d_ub witnesses -- the ONLY allowed distance inputs), prices the cost
gates, and emits one prediction record per row to predictions.jsonl.

The prediction records are the frozen artifact Phase 2+ is scored
against; nothing in here runs a closure, an exact solve, or a floor
certification on any Level-0 target code.
"""

from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import duckdb

sys.path.insert(0, str(Path(__file__).resolve().parent))
import dtt_lib as D                                    # noqa: E402
from dtt_lib import OUT                                # noqa: E402

T0 = time.monotonic()


def log(msg: str) -> None:
    print(f"[{time.monotonic()-T0:7.1f}s] {msg}", flush=True)


def route_str(folds: list[list[int]]) -> str:
    return ", ".join(f"{'xy'[ax]}->{nm}" for ax, nm in folds) or "(none)"


def headline(row: dict, pred: dict) -> str:
    v = pred["verdict"]
    if v == "REFUSED":
        return ("refused: no deck (|G| = "
                f"{row['ell'] * row['m']} is odd; A35 condition C1)")
    if v == "NO-ROUTE":
        return ("no usable route: every available Z2 fold kills all "
                "logical content (k=0 quotient); tower offers no census "
                "win -- predict Phase 2 falls back to direct methods")
    w = pred["g5_window"]
    lo, hi = w["window_lo"], w["window_hi"]
    est = " (window top is bounded-only)" if w["ceiling_is_upper_estimate"] \
        else ""
    if w["ceiling"] is not None and row.get("d_ub") is not None \
            and w["ceiling"] < row["d_ub"] and pred["stratum"] != "anchor":
        est += (f" [G5: if true d > {w['ceiling']}, route certifies the "
                f"FLOOR {w['ceiling']} only]")
    cost = pred["cost_operative"]
    rt = route_str(pred["screen"]["folds"])
    stall = ""
    if w["chain_stalls"]:
        s0 = w["chain_stalls"][0]
        stall = (f"; NOTE route stalls at rung {s0['rung']} "
                 f"(floor caps at {2 * s0['d_base']})")
    parity_note = ""
    if not pred["parity_ok"]:
        parity_note = ("; parity layer lost (even-weight polynomial): "
                       "odd cycle weights possible, W_eff = d-1, no "
                       "beta=0 kills, odd d values allowed (A35 L3)")
    if v == "GREEN":
        t_est = cost["wall_estimate_s"]
        return (f"closes at certificate tier in ~{t_est:.0f} s "
                f"(x3 uncertainty) via route [{rt}] with "
                f"d in [{lo}, {hi}]{est}{stall}{parity_note}")
    if v == "AMBER":
        g = ("G2 cap " + str(cost["cap_max"])
             if cost["cap_max"] > D.GREEN_CAP else
             f"G1 nodes 1e{cost['log10_nodes_per_level'][-1]}")
        return (f"AMBER at {g}: extended budget likely closes it via "
                f"[{rt}], d in [{lo}, {hi}]{est}{stall}{parity_note}")
    binding = []
    if cost["bottom_nodes"] > D.AMBER_NODES:
        binding.append(f"G1 (bottom census 1e"
                       f"{cost['log10_nodes_per_level'][-1]})")
    if cost["cap_max"] > D.AMBER_CAP:
        binding.append(f"G2 (cap {cost['cap_max']})")
    return (f"RED at {' + '.join(binding) or 'gates'}: predict NON-closure "
            f"within certificate-tier budgets via [{rt}]"
            f"{est}{stall}{parity_note}")


def main() -> None:
    con = duckdb.connect(str(D.MAIN_DB), read_only=True)
    screens = {}
    with (OUT / "screens_cache.jsonl").open() as fh:
        for line in fh:
            o = json.loads(line)
            screens[o["instance_id"]] = o["screen"]
    rows = [json.loads(x) for x in (OUT / "cohort.jsonl").open()]
    log(f"{len(rows)} cohort rows, {len(screens)} cached screens")

    pred_f = (OUT / "predictions.jsonl").open("w")
    stats = {"verdicts": {}, "rank_law_violations": 0,
             "lemma1_failures": 0, "parity_consistency_flags": 0,
             "g5_containment": {"pass": 0, "fail": 0, "indeterminate": 0}}
    occupancy: dict[tuple[int, str], int] = {}

    for row in rows:
        iid = row["instance_id"]
        pred: dict = {
            "instance_id": iid, "code_id": row["code_id"],
            "stratum": row["stratum"],
            "group": row["group_struct"], "ell": row["ell"], "m": row["m"],
            "n": row["n"], "k": row["k"],
            "A_poly": row["A_poly"], "B_poly": row["B_poly"],
            "v2_depth": row.get("v2_depth"),
            "d_lb": row.get("d_lb"), "d_ub": row.get("d_ub"),
            "d_exact_corpus": row.get("d_exact"),
            "d_ub_provenance": row.get("d_ub_provenance")
                or row.get("d_provenance"),
            "predicted_at": None,        # stamped at write time below
        }

        # ---------------------------------------------- refusal stratum
        if row["stratum"] == "scope-control":
            pred["verdict"] = "REFUSED"
            pred["refusal_reason"] = ("odd |G|: no free Z2 deck exists "
                                      "(A35 conditions map, item 1)")
            pred["falsifier"] = (
                "any Phase-2 construction of a free Z2 BB deck for this "
                "code, or any tower-calculus closure of it, falsifies "
                "condition C1 as stated")
            pred["headline"] = headline(row, pred)
            stats["verdicts"]["REFUSED"] = \
                stats["verdicts"].get("REFUSED", 0) + 1
            jline(pred_f, pred)
            continue

        s = screens[iid]
        pred["screen"] = s               # full structural screen, frozen

        # ------------------------------------------- structural law audit
        rl_viol = [r["rung"] for r in s["rungs"]
                   if not r.get("rank_law_holds", True)
                   or not r.get("codim_lift_law_holds", True)
                   or r.get("exact_cover") is False]
        if rl_viol:
            stats["rank_law_violations"] += len(rl_viol)
        pred["rank_law_violation_rungs"] = rl_viol
        if s["rungs_lemma1_failures"]:
            stats["lemma1_failures"] += len(s["rungs_lemma1_failures"])
        # (R)-trio prediction per rung (A12): exact_base & sigma_id <=> R
        pred["r_trio_mismatch_rungs"] = [
            r["rung"] for r in s["rungs"] if "LEMMA1_FAILED" not in r
            and r.get("exact_base") is not None
            and ((r["exact_base"] and r["sigma_id"]) != r["R_holds"])]

        # parity consistency (Lemma 2): odd polys => all cycles even
        parity_flags = [lv["level"] for lv in s["levels"]
                        if lv["parity_scope"] and not lv["cycles_all_even"]]
        if parity_flags:
            stats["parity_consistency_flags"] += len(parity_flags)
        pred["lemma2_violation_levels"] = parity_flags
        parity_ok = all(lv["parity_scope"] and lv["cycles_all_even"]
                        for lv in s["levels"])
        pred["parity_ok"] = parity_ok

        for p in s["pairs"]:
            if "regime" in p:
                key = (row.get("v2_depth"), p["regime"])
                occupancy[key] = occupancy.get(key, 0) + 1

        # ------------------------------------------- quotient distances
        dinfo: dict[int, dict] = {}
        for lv in s["levels"][1:]:
            lm = tuple(lv["lm"])
            tA = D.poly_terms(lv["A"], lm)
            tB = D.poly_terms(lv["B"], lm)
            try:
                dinfo[lv["level"]] = D.quotient_distance(
                    con, lm, tA, tB, lv["level"])
            except Exception as e:
                dinfo[lv["level"]] = {
                    "d": None, "d_is_exact": False,
                    "provenance": f"unavailable ({type(e).__name__}: {e})"}
        dinfo = {j: v for j, v in dinfo.items() if v.get("d") is not None}
        pred["quotient_distances"] = {str(j): v for j, v in dinfo.items()}

        # ------------------------------------------------- G5 window
        d_ub_top = row.get("d_ub")
        w = D.g5_window(s["levels"], dinfo, d_ub_top, parity_ok)
        pred["g5_window"] = w

        # ---------------------------------------------------- cost gates
        slack = 2 if parity_ok else 1
        costs = []
        W_op = None
        if row["stratum"] == "anchor" and row.get("d_exact") is not None:
            W_op = row["d_exact"] - slack
            cb = D.cost_block(s["levels"], W_op, s["mu_min"])
            cb["question"] = (f"calibration: certify d = {row['d_exact']} "
                              "(corpus-exact)")
            costs.append(cb)
            if w["window_hi"] is not None and w["window_hi"] != \
                    row["d_exact"]:
                W2 = w["window_hi"] - slack
                if W2 >= 1:
                    cb2 = D.cost_block(s["levels"], W2, s["mu_min"])
                    cb2["question"] = (f"window top d = {w['window_hi']} "
                                       "(from quotient chain, no d_exact "
                                       "consumed)")
                    costs.append(cb2)
            # G5 containment self-test (anchors only; uses corpus values)
            if w["ceiling"] is not None \
                    and not w["ceiling_is_upper_estimate"]:
                ok = row["d_exact"] <= w["ceiling"]
                pred["g5_containment"] = {
                    "d_exact": row["d_exact"], "ceiling": w["ceiling"],
                    "pass": bool(ok),
                    "interpretation": (
                        "d within route ceiling" if ok else
                        "route ceiling below true d: predict FLOOR-ONLY "
                        f"closure at {w['ceiling']} via this route "
                        "(exceeding it would falsify G5)")}
                stats["g5_containment"]["pass" if ok else "fail"] += 1
            else:
                pred["g5_containment"] = {"pass": None,
                                          "note": "chain has bounded-only "
                                                  "terms"}
                stats["g5_containment"]["indeterminate"] += 1
        elif w["window_hi"] is not None:
            W_op = w["window_hi"] - slack
            cb = D.cost_block(s["levels"], W_op, s["mu_min"])
            cb["question"] = (f"certify d = {w['window_hi']} "
                              "(predicted window top)")
            costs.append(cb)
            if d_ub_top is not None and d_ub_top != w["window_hi"]:
                W2 = d_ub_top - slack
                cb2 = D.cost_block(s["levels"], W2, s["mu_min"])
                cb2["question"] = f"certify d = {d_ub_top} (= d_ub)"
                costs.append(cb2)
        pred["costs"] = costs
        pred["cost_operative"] = costs[0] if costs else None
        pred["W_operative"] = W_op
        pred["parity_slack_used"] = slack
        # G3 envelope flag (A35: 2^12 dispatch demonstrated; larger is
        # outside the measured envelope, recorded but not verdict-changing)
        if costs and costs[0]["sector_dispatch_max"] is not None \
                and costs[0]["sector_dispatch_max"] > 4096:
            pred["g3_beyond_envelope"] = costs[0]["sector_dispatch_max"]

        # ------------------------------------------------------ verdict
        if s["depth_used"] == 0:
            pred["verdict"] = "NO-ROUTE"
        elif costs:
            pred["verdict"] = costs[0]["verdict"]
        else:
            pred["verdict"] = "NO-WINDOW"
        stats["verdicts"][pred["verdict"]] = \
            stats["verdicts"].get(pred["verdict"], 0) + 1
        pred["headline"] = headline(row, pred)
        jline(pred_f, pred)

    pred_f.close()
    occ_table = {f"depth{d_}|{reg}": c
                 for (d_, reg), c in sorted(occupancy.items())}
    (OUT / "aggregates.json").write_text(json.dumps(
        {"stats": stats, "depth_regime_occupancy": occ_table,
         "wall_s": round(time.monotonic() - T0, 1)}, indent=1))
    log(f"stats: {json.dumps(stats)}")
    log(f"depth x regime occupancy (pair count): {occ_table}")


def jline(fh, obj) -> None:
    import datetime as _dt
    if isinstance(obj, dict) and obj.get("predicted_at", "x") is None:
        obj["predicted_at"] = _dt.datetime.now(_dt.UTC).isoformat()
    fh.write(json.dumps(obj) + "\n")
    fh.flush()


if __name__ == "__main__":
    main()
