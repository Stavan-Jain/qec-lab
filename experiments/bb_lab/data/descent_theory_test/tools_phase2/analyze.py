"""Phase-4 scorecard builder: computes the five pre-registered criteria
from phase2_results.jsonl / phase3_groundtruth.jsonl / phase3_armB.jsonl
and prints a JSON summary (consumed while drafting PHASE4_SCORECARD.md).

Retrodictions are excluded from (iii)-(iv) by construction (they were
never run in Stage A).  Scope-control rows enter only the refusal count.
"""

from __future__ import annotations

import json
import statistics
import sys
from collections import Counter, defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
DTT = HERE.parent


def load(p: str) -> list[dict]:
    f = DTT / p
    return [json.loads(line) for line in f.open()] if f.exists() else []


def main() -> None:
    preds = {r["instance_id"]: r for r in load("predictions.jsonl")}
    preds2 = {r["instance_id"]: r for r in load("predictions_batch2.jsonl")}
    allpred = {**preds, **preds2}
    res = load("phase2_results.jsonl")
    # last record wins per (row, question) — reruns supersede
    by_key: dict[tuple, dict] = {}
    for r in res:
        if r.get("kind") == "closure":
            by_key[(r["instance_id"], r.get("q"))] = r
    cl = list(by_key.values())
    refusals = [r for r in res if r.get("kind") == "refusal-check"]
    controls = [r for r in res if r.get("kind") == "control"]
    out: dict = {}

    # ---------------- criterion (i): structural laws
    n_rungs = sum(len(r.get("rungs_measured") or []) for r in cl)
    viols = []
    for r in cl:
        for v in r["eval"].get("structural_violations", []):
            viols.append({"instance_id": r["instance_id"], **v})
    # dedupe rungs measured multiple times (per question) — count unique
    uniq_rungs = set()
    rank_fail = exact_fail = 0
    for r in cl:
        for g in r.get("rungs_measured") or []:
            key = (r["instance_id"], g["rung"])
            if key in uniq_rungs:
                continue
            uniq_rungs.add(key)
            if not g.get("rank_law_holds", True):
                rank_fail += 1
            if g.get("exactness_at_cover_holds") is False:
                exact_fail += 1
    out["criterion_i"] = {
        "rung_measurements_total": n_rungs,
        "unique_rungs": len(uniq_rungs),
        "rank_law_failures": rank_fail,
        "exactness_failures": exact_fail,
        "violation_records": viols,
        "verdict": "PASS" if not viols else "FAIL",
    }

    # ---------------- criterion (ii): node formula within x1.1
    ratios = []
    per_q = []
    fallbacks = []
    for r in cl:
        e = r["eval"]
        nr = e.get("node_ratio_vs_prediction")
        if nr is None:
            continue
        fb = any(v.get("single_window_fallback")
                 for k2, v in (r.get("stages") or {}).items()
                 if k2.endswith("_census"))
        if fb:
            # D2 fallback census: different r-pair semantics, excluded
            # from the formula sample by design (see scorecard D2)
            fallbacks.append({"iid": r["instance_id"], "q": r["q"],
                              "vs_pred": nr})
            continue
        ratios.append(nr)
        per_q.append({"iid": r["instance_id"], "q": r["q"], "vs_pred": nr,
                      "vs_formula": e.get("node_ratio_vs_formula")})
    bad = [p for p in per_q if not (1 / 1.1 <= p["vs_pred"] <= 1.1)]
    out["criterion_ii"] = {
        "executed_censuses": len(ratios),
        "ratio_min": min(ratios) if ratios else None,
        "ratio_max": max(ratios) if ratios else None,
        "outside_1p1": bad,
        "fallback_censuses_excluded_D2": fallbacks,
        "verdict": "PASS" if ratios and not bad else
                   ("NO-DATA" if not ratios else "FAIL"),
    }

    # ---------------- criterion (iii): GREEN closure / RED non-closure
    # row-level on OPERATIVE questions, excluding retrodictions (never run)
    def op_record(iid):
        cand = [r for r in cl if r["instance_id"] == iid
                and r.get("operative")]
        return cand[0] if cand else None

    strata = defaultdict(list)
    green_rows = []
    red_rows = []
    amber_rows = []
    for iid, p in allpred.items():
        if p.get("row_class") == "retrodiction":
            continue
        if p.get("stratum") == "scope-control":
            continue
        rec = op_record(iid)
        if rec is None:
            continue
        e = rec["eval"]
        verdict = e["predicted_verdict"]
        row = {"iid": iid, "stratum": p["stratum"], "outcome": e["outcome"],
               "closed": bool(e.get("closed")),
               "closed_any_time": bool(e.get("closed_any_time")),
               "cex_w": e.get("counterexample_weight"),
               "wall_ratio": e.get("wall_ratio"),
               "stop": e.get("stop_reason")}
        strata[p["stratum"]].append(row)
        (green_rows if verdict == "GREEN" else
         red_rows if verdict == "RED" else amber_rows).append(row)
    g_closed = sum(1 for r in green_rows if r["closed"])
    g_cex = sum(1 for r in green_rows if r["outcome"] == "COUNTEREXAMPLE")
    g_stop = sum(1 for r in green_rows
                 if r["outcome"] in ("BLOWUP_STOP", "ENVELOPE_STOP"))
    g_kill = sum(1 for r in green_rows
                 if r["outcome"] in ("BUDGET_KILL", "HARD_KILL"))
    g_other = len(green_rows) - g_closed - g_cex - g_stop - g_kill
    resolved = g_closed + g_cex           # certificate-tier resolution
    out["criterion_iii"] = {
        "green_rows": len(green_rows),
        "green_closed_as_registered": g_closed,
        "green_rate_as_registered": round(g_closed / len(green_rows), 3)
        if green_rows else None,
        "green_counterexample_exact": g_cex,
        "green_engine_stops": g_stop,
        "green_budget_kills": g_kill,
        "green_other": g_other,
        "green_resolved_certificate_tier": resolved,
        "green_resolved_rate": round(resolved / len(green_rows), 3)
        if green_rows else None,
        "red_rows": len(red_rows),
        "red_closures (falsify gate if >0)": sum(
            1 for r in red_rows if r["closed"]),
        "amber_rows": len(amber_rows),
        "amber_closed": sum(1 for r in amber_rows if r["closed"]),
        "refusals_confirmed": sum(bool(r.get("refusal_confirmed"))
                                  for r in refusals),
        "refusals_total": len(refusals),
        "verdict_as_registered": (
            "PASS" if green_rows
            and g_closed / len(green_rows) >= 0.9
            and not any(r["closed"] for r in red_rows) else "FAIL"),
    }
    out["criterion_iii"]["by_stratum"] = {
        s: {"n": len(v),
            "closed": sum(1 for r in v if r["closed"]),
            "cex": sum(1 for r in v if r["outcome"] == "COUNTEREXAMPLE")}
        for s, v in sorted(strata.items())}

    # secondary questions summary
    sec = [r for r in cl if not r.get("operative")]
    sec_sum = Counter((r["eval"]["predicted_verdict"],
                       r["eval"]["outcome"]) for r in sec)
    out["secondary_questions"] = {f"{k[0]}->{k[1]}": v
                                  for k, v in sorted(sec_sum.items())}
    red_sec_closed = [r["instance_id"] for r in sec
                      if r["eval"]["predicted_verdict"] == "RED"
                      and r["eval"].get("closed")]
    out["red_secondary_closures_conservatism"] = red_sec_closed

    # ---------------- criterion (iv): ground truth agreement
    gt = load("phase3_groundtruth.jsonl")
    agree = [r for r in gt if r.get("agreement") is True]
    disagree = [r for r in gt if r.get("agreement") is False]
    out["criterion_iv"] = {
        "rows_run": len(gt),
        "exact_solves": sum(1 for r in gt if r.get("outcome") == "exact"),
        "timeouts": sum(1 for r in gt if r.get("outcome") == "timeout"),
        "agree": len(agree),
        "disagree": len(disagree),
        "disagree_rows": disagree,
        "verdict": ("PASS" if not disagree and agree else
                    "FAIL" if disagree else "NO-DATA"),
    }

    # ---------------- criterion (v): equal-compute head-to-head
    armB_raw = [r for r in load("phase3_armB.jsonl")
                if not r.get("superseded")]
    armB_by = {}
    for r in armB_raw:
        armB_by[r["instance_id"]] = r   # last record wins (re-runs)
    armB = list(armB_by.values())
    if armB:
        # arm A per-row resolution (certificate tier)
        armA = {}
        for r in cl:
            iid = r["instance_id"]
            e = r["eval"]
            a = armA.setdefault(iid, {"resolved": False, "floor": 0,
                                      "exact": None, "cpu": 0.0})
            a["cpu"] += float(e.get("cpu_s") or 0)
            if e["outcome"] == "CERTIFIED_FLOOR":
                a["floor"] = max(a["floor"], e.get("floor") or 0)
                p = allpred.get(iid, {})
                tgt = (p.get("cost_operative") or {}).get("W", 0) + \
                    (2 if p.get("parity_ok") else 1)
                if e.get("floor", 0) >= tgt:
                    a["resolved"] = True
            if e["outcome"] == "COUNTEREXAMPLE":
                a["exact"] = e["counterexample_weight"]
                a["resolved"] = True
        both = [(iid, armA[iid], b) for b in armB
                for iid in [b["instance_id"]] if iid in armA]
        aw = sum(1 for _, a, b in both if a["resolved"])
        bw = sum(1 for _, a, b in both if b.get("closed_v"))
        only_a = [iid for iid, a, b in both
                  if a["resolved"] and not b.get("closed_v")]
        only_b = [iid for iid, a, b in both
                  if not a["resolved"] and b.get("closed_v")]
        out["criterion_v"] = {
            "paired_rows": len(both),
            "armA_resolved": aw,
            "armB_closed": bw,
            "only_armA": only_a,
            "only_armB": only_b,
            "armA_total_cpu_s": round(sum(a["cpu"] for _, a, _ in both), 1),
            "armB_total_wall_s": round(sum(b["wall_s"]
                                           for _, _, b in both), 1),
            "armB_exact": sum(1 for _, _, b in both
                              if b.get("outcome") == "exact"),
            "armB_timeout": sum(1 for _, _, b in both
                                if b.get("outcome") == "timeout"),
        }

    # ---------------- controls, censoring, calibration
    out["controls"] = [
        {k: c.get(k) for k in ("control", "target", "passed", "outcome",
                               "instance_id")}
        for c in controls]
    censored = [{"iid": r["instance_id"], "q": r["q"],
                 "outcome": r["eval"]["outcome"],
                 "budget_s": r["engine"].get("budget_s"),
                 "reason": r["eval"].get("stop_reason")}
                for r in cl
                if r["eval"]["outcome"] in ("BUDGET_KILL", "HARD_KILL")]
    out["censored_rows"] = censored
    stops = Counter(r["eval"]["outcome"] for r in cl)
    out["outcome_histogram_all_questions"] = dict(stops)

    # wall ratios (calibration of the x3 wall model)
    wr = [r["eval"]["wall_ratio"] for r in cl
          if r["eval"].get("wall_ratio") is not None]
    if wr:
        out["wall_ratio"] = {
            "n": len(wr), "median": round(statistics.median(wr), 3),
            "p90": round(sorted(wr)[int(0.9 * len(wr))], 3),
            "max": round(max(wr), 3)}

    # anchor calibration
    anch = [r for r in cl if r.get("stratum") == "anchor"]
    out["anchor_calibration"] = {
        "calibration_floors_ok": sum(
            1 for r in anch if r.get("operative")
            and r["eval"].get("closed_any_time")),
        "window_cex_agree_corpus": sum(
            1 for r in anch if not r.get("operative")
            and r["eval"].get("counterexample_weight")
            == (allpred.get(r["instance_id"], {}).get("d_exact_corpus"))),
        "window_questions": sum(1 for r in anch if not r.get("operative")),
    }

    json.dump(out, sys.stdout, indent=1, default=str)
    print()


if __name__ == "__main__":
    main()
