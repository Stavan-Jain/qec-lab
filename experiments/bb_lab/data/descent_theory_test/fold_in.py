"""Fold-in: ingest the order-144 sweep as cohort BATCH 2 (depth-4
stratum) of the descent-theory test.

DO NOT run this while the sweep is still writing
(`experiments/bb_lab/data/order144_sweep/`). It was written blind to the
sweep's exact schema (the running directory was not read); the FIELD_MAP
below tries the plausible key spellings and FAILS LOUDLY listing the
keys it saw if none match — fix the map, do not guess silently.

Honesty rule (FOLD_IN.md; encoded in `row_class`):
  * sweep rows that already carry an exact/certificate outcome are
    RETRODICTION rows — the outcome predates the prediction, so they
    test the screen's CALIBRATION, not its foresight;
  * only the sweep's UNKNOWN / bounded-only tail yields genuine
    pre-registered predictions for Phase 2.
Both classes get the same blind screen: the row's own d_exact is NEVER
fed to the window/cost machinery; for retrodiction rows it is attached
separately as `calibration_truth` plus an immediate `retro_check`.

Outputs (all in this directory):
  cohort_batch2.jsonl, predictions_batch2.jsonl, aggregates_batch2.json,
  MANIFEST_batch2.sha256

Usage:
  cd experiments/bb_lab && uv run python data/descent_theory_test/fold_in.py
      [--results PATH] [--limit N]
"""

from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import json
import sys
import time
from pathlib import Path

import duckdb

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE / "tools"))
import dtt_lib as D                                    # noqa: E402

T0 = time.monotonic()
DEFAULT_RESULTS = D.LAB / "data" / "order144_sweep" / "results.jsonl"

# candidate key spellings, first hit wins
FIELD_MAP = {
    "instance_id": ["instance_id", "iid", "id", "canonical_id"],
    "ell": ["ell", "l"],
    "m": ["m"],
    "group_struct": ["group_struct", "group", "group_label"],
    "A_poly": ["A_poly", "A", "a_poly", "polyA", "A_str"],
    "B_poly": ["B_poly", "B", "b_poly", "polyB", "B_str"],
    "n": ["n"],
    "k": ["k"],
    "d_exact": ["d_exact", "d"],
    "d_ub": ["d_ub", "dub", "d_upper", "ub"],
    "d_lb": ["d_lb", "lb", "d_lower", "floor"],
    "d_method": ["d_method", "method", "provenance", "outcome"],
}
CERT_METHOD_TOKENS = ("cert", "sat", "exact", "milp", "maxsat", "bz",
                      "census", "kernel")


def log(msg: str) -> None:
    print(f"[{time.monotonic()-T0:7.1f}s] {msg}", flush=True)


def map_row(raw: dict) -> dict:
    out = {}
    for field, keys in FIELD_MAP.items():
        for key in keys:
            if key in raw and raw[key] is not None:
                out[field] = raw[key]
                break
    missing = [f for f in ("A_poly", "B_poly") if f not in out]
    if "ell" not in out or "m" not in out:
        g = out.get("group_struct")
        if isinstance(g, str) and g.startswith("Z") and "xZ" in g:
            try:
                a, b = g[1:].split("xZ")
                out["ell"], out["m"] = int(a), int(b)
            except ValueError:
                pass
    if "ell" not in out or "m" not in out:
        missing.append("ell/m")
    if missing:
        raise KeyError(
            f"cannot map sweep row: missing {missing}; row keys were "
            f"{sorted(raw.keys())} -- edit FIELD_MAP in fold_in.py")
    return out


def classify(row: dict) -> str:
    """RETRODICTION if the sweep already fixed the outcome, else
    PREDICTION (the genuine pre-registered tail)."""
    if row.get("d_exact") is not None:
        return "retrodiction"
    meth = str(row.get("d_method") or "").lower()
    if any(t in meth for t in CERT_METHOD_TOKENS) and row.get("d_lb"):
        return "retrodiction"
    return "prediction"


def predict_one(con, row: dict) -> dict:
    """Mirror of make_predictions.py record assembly (PROTOCOL.md §4),
    blind to the row's own d_exact."""
    ell, m = int(row["ell"]), int(row["m"])
    pred: dict = {
        "instance_id": str(row.get("instance_id") or
                           D.canonical_id((ell, m),
                                          D.poly_terms(row["A_poly"],
                                                       (ell, m)),
                                          D.poly_terms(row["B_poly"],
                                                       (ell, m)))[0]),
        "batch": 2, "stratum": "order144-depth4",
        "group": row.get("group_struct") or f"Z{ell}xZ{m}",
        "ell": ell, "m": m,
        "n": row.get("n") or 2 * ell * m, "k": row.get("k"),
        "A_poly": row["A_poly"], "B_poly": row["B_poly"],
        "v2_depth": D.v2(ell * m),
        "d_ub": row.get("d_ub"),
        "row_class": classify(row),
        "predicted_at": dt.datetime.now(dt.UTC).isoformat(),
    }
    if (ell * m) % 2 == 1:
        pred["verdict"] = "REFUSED"
        pred["refusal_reason"] = "odd |G|: no free Z2 deck (A35 C1)"
        return pred
    s = D.screen_structure(ell, m, row["A_poly"], row["B_poly"],
                           pred["instance_id"][:8])
    pred["screen"] = s
    if not s["levels"]:                  # target itself k=0 / degenerate
        pred["verdict"] = "NO-ROUTE"
        pred["refusal_reason"] = (
            "target level unusable: "
            + (s["truncated"] or {}).get("reason", "no k>0 level"))
        return pred
    if pred.get("k") is None:
        pred["k"] = s["levels"][0]["k"]
    pred["rank_law_violation_rungs"] = [
        r["rung"] for r in s["rungs"]
        if not r.get("rank_law_holds", True)
        or not r.get("codim_lift_law_holds", True)
        or r.get("exact_cover") is False]
    parity_ok = all(lv["parity_scope"] and lv["cycles_all_even"]
                    for lv in s["levels"])
    pred["parity_ok"] = parity_ok
    dinfo = {}
    for lv in s["levels"][1:]:
        lm2 = tuple(lv["lm"])
        tA = D.poly_terms(lv["A"], lm2)
        tB = D.poly_terms(lv["B"], lm2)
        try:
            di = D.quotient_distance(con, lm2, tA, tB, lv["level"])
            if di.get("d") is not None:
                dinfo[lv["level"]] = di
        except Exception as e:
            pred.setdefault("quotient_errors", []).append(
                f"L{lv['level']}: {type(e).__name__}: {e}")
    pred["quotient_distances"] = {str(j): v for j, v in dinfo.items()}
    w = D.g5_window(s["levels"], dinfo, pred["d_ub"], parity_ok)
    pred["g5_window"] = w
    slack = 2 if parity_ok else 1
    costs = []
    if w["window_hi"] is not None and w["window_hi"] - slack >= 1:
        cb = D.cost_block(s["levels"], w["window_hi"] - slack, s["mu_min"])
        cb["question"] = f"certify d = {w['window_hi']} (window top)"
        costs.append(cb)
        if pred["d_ub"] is not None and pred["d_ub"] != w["window_hi"]:
            cb2 = D.cost_block(s["levels"], pred["d_ub"] - slack,
                               s["mu_min"])
            cb2["question"] = f"certify d = {pred['d_ub']} (= d_ub)"
            costs.append(cb2)
    pred["costs"] = costs
    pred["cost_operative"] = costs[0] if costs else None
    if s["depth_used"] == 0:
        pred["verdict"] = "NO-ROUTE"
    elif costs:
        pred["verdict"] = costs[0]["verdict"]
    else:
        pred["verdict"] = "NO-WINDOW"
    # retrodiction rows: attach the known truth SEPARATELY + retro check
    if pred["row_class"] == "retrodiction":
        truth = {"d_exact": row.get("d_exact"),
                 "d_lb": row.get("d_lb"),
                 "d_method": row.get("d_method")}
        pred["calibration_truth"] = truth
        if row.get("d_exact") is not None and w["ceiling"] is not None:
            pred["retro_check"] = {
                "d_in_window": bool(w["window_lo"] <= row["d_exact"]
                                    <= (w["window_hi"] or 10**9)),
                "ceiling_contains_d": bool(row["d_exact"] <= w["ceiling"]),
                "ceiling_is_upper_estimate":
                    w["ceiling_is_upper_estimate"]}
    return pred


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--results", default=str(DEFAULT_RESULTS))
    ap.add_argument("--limit", type=int, default=None)
    args = ap.parse_args()
    res = Path(args.results)
    if not res.exists():
        sys.exit(f"results file not found: {res} -- is the sweep done?")

    known_batch1 = set()
    c1 = HERE / "cohort.jsonl"
    if c1.exists():
        known_batch1 = {json.loads(x)["instance_id"] for x in c1.open()}

    raws, bad = [], 0
    with res.open() as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            try:
                raws.append(json.loads(line))
            except json.JSONDecodeError:
                bad += 1
    if args.limit:
        raws = raws[:args.limit]
    log(f"{len(raws)} sweep rows ({bad} malformed lines skipped)")

    con = duckdb.connect(str(D.MAIN_DB), read_only=True)
    coh_f = (HERE / "cohort_batch2.jsonl").open("w")
    pred_f = (HERE / "predictions_batch2.jsonl").open("w")
    stats = {"verdicts": {}, "row_class": {}, "rank_law_violations": 0,
             "retro_ceiling_fail": 0, "skipped_batch1_dupes": 0}
    occ: dict[str, int] = {}
    for raw in raws:
        try:
            row = map_row(raw)
        except KeyError as e:
            sys.exit(str(e))
        pred = predict_one(con, row)
        if pred["instance_id"] in known_batch1:
            stats["skipped_batch1_dupes"] += 1
            continue
        coh_f.write(json.dumps({**row,
                                "instance_id": pred["instance_id"],
                                "stratum": pred["stratum"],
                                "row_class": pred["row_class"]}) + "\n")
        pred_f.write(json.dumps(pred) + "\n")
        for fh in (coh_f, pred_f):
            fh.flush()
        stats["verdicts"][pred["verdict"]] = \
            stats["verdicts"].get(pred["verdict"], 0) + 1
        stats["row_class"][pred["row_class"]] = \
            stats["row_class"].get(pred["row_class"], 0) + 1
        stats["rank_law_violations"] += \
            len(pred.get("rank_law_violation_rungs", []))
        rc = pred.get("retro_check")
        if rc and not rc["ceiling_contains_d"] \
                and not rc["ceiling_is_upper_estimate"]:
            stats["retro_ceiling_fail"] += 1
        for p in pred.get("screen", {}).get("pairs", []):
            if "regime" in p:
                key = f"depth{pred['v2_depth']}|{p['regime']}"
                occ[key] = occ.get(key, 0) + 1
    coh_f.close()
    pred_f.close()
    (HERE / "aggregates_batch2.json").write_text(json.dumps(
        {"stats": stats, "depth_regime_occupancy": occ,
         "results_file": str(res),
         "wall_s": round(time.monotonic() - T0, 1)}, indent=1))
    # batch-2 manifest
    lines = [f"# MANIFEST_batch2 -- {dt.datetime.now(dt.UTC).isoformat()}"]
    for f in ["cohort_batch2.jsonl", "predictions_batch2.jsonl",
              "aggregates_batch2.json"]:
        h = hashlib.sha256((HERE / f).read_bytes()).hexdigest()
        lines.append(f"{h}  {f}")
    (HERE / "MANIFEST_batch2.sha256").write_text("\n".join(lines) + "\n")
    log(f"batch 2 done: {stats}")
    log(f"occupancy: {occ}")


if __name__ == "__main__":
    main()
