"""Phase-2 orchestrator: budget-scaled closures against the frozen
predictions (PROTOCOL.md criteria (i)-(iii) data collection).

Per in-scope row, every pre-registered cost question is executed via
dtt_close_worker.py in its own nice'd subprocess under a HARD kill at
budget + grace.  Budget per question: min(10 x wall_estimate, 1800 s).
[Deviation note: the tasking said 3x; PROTOCOL criterion (iii) scores
closure "within 10x the recorded wall estimate", so a 3x kill would
censor the pre-registered scoring window — PROTOCOL wins.  Recorded here
and in the scorecard.]

Scope-control rows get the mechanical refusal check (odd |G| -> no
index-2 axis subgroup -> no free Z2 deck -> the frozen REFUSED verdict).

Results append to ../phase2_results.jsonl (crash-safe resume by
(instance_id, q) key).  Run:

  cd experiments/bb_lab && caffeinate -ims nice -n 5 \
      uv run python data/descent_theory_test/tools_phase2/run_phase2.py
"""

from __future__ import annotations

import argparse
import json
import os
import shutil
import signal
import subprocess
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
DTT = HERE.parent
LAB = DTT.parent.parent
RESULTS = DTT / "phase2_results.jsonl"
LOGS = HERE / "logs"
WORK = HERE / "work"
PY = sys.executable

BUDGET_CAP_S = 1800.0
BUDGET_MULT = 10.0          # PROTOCOL (iii) closure window (see docstring)
GRACE_S = 45.0


def v2(x: int) -> int:
    c = 0
    while x % 2 == 0 and x > 0:
        x //= 2
        c += 1
    return c


def free_mb() -> float:
    return shutil.disk_usage(str(LAB)).free / 1e6


def load_jsonl(p: Path) -> list[dict]:
    return [json.loads(line) for line in p.open()] if p.exists() else []


def emit(rec: dict) -> None:
    with RESULTS.open("a") as f:
        f.write(json.dumps(rec) + "\n")


def run_question(row: dict, qidx: int, cost: dict, extra: dict | None = None
                 ) -> dict:
    iid = row["instance_id"]
    tag = f"{iid}:q{qidx}"
    wd = WORK / f"{iid[:12]}_q{qidx}"
    wd.mkdir(parents=True, exist_ok=True)
    budget = min(BUDGET_MULT * float(cost["wall_estimate_s"]), BUDGET_CAP_S)
    spec = {
        "ell": row["ell"], "m": row["m"],
        "A": row["A_poly"], "B": row["B_poly"],
        "folds": row["screen"]["folds"],
        "W": cost["W"],
        "threads": 4,
        "workdir": str(wd),
        "budget_s": budget,
        "tag": tag,
        "expected_levels": [{"k": lv["k"], "kappa": lv["kappa"]}
                            for lv in row["screen"]["levels"]],
    }
    if extra:
        spec.update(extra)
    spec_p = wd / "spec.json"
    spec_p.write_text(json.dumps(spec))
    log_p = LOGS / f"{iid[:12]}_q{qidx}.log"
    t0 = time.time()
    with log_p.open("w") as lf:
        proc = subprocess.Popen(
            ["nice", "-n", "10", PY, str(HERE / "dtt_close_worker.py"),
             str(spec_p)],
            stdout=lf, stderr=subprocess.STDOUT,
            start_new_session=True, cwd=str(LAB))
        killed = False
        try:
            proc.wait(timeout=budget + GRACE_S)
        except subprocess.TimeoutExpired:
            killed = True
            try:
                os.killpg(proc.pid, signal.SIGKILL)
            except ProcessLookupError:
                pass
            proc.wait()
    wall = time.time() - t0
    res: dict = {"outcome": "HARD_KILL"} if killed else {}
    for line in log_p.read_text().splitlines():
        if line.startswith("CLOSE_JSON "):
            res = json.loads(line[len("CLOSE_JSON "):])
            break
    if killed:
        res["outcome"] = "HARD_KILL"
        res.setdefault("reason", f"killed at budget {budget:.0f}s + grace")
        res["cpu_s"] = res.get("cpu_s", budget)  # conservative charge
    res["parent_wall_s"] = round(wall, 2)
    res["budget_s"] = budget
    return res


def eval_question(row: dict, cost: dict, res: dict) -> dict:
    """Prediction-vs-outcome deltas for one executed question."""
    d: dict = {"predicted_verdict": cost["verdict"],
               "predicted_wall_s": cost["wall_estimate_s"],
               "predicted_bottom_nodes": cost["bottom_nodes"],
               "question": cost["question"], "W": cost["W"]}
    out = res.get("outcome")
    d["outcome"] = out
    cen = None
    for k2, v in (res.get("stages") or {}).items():
        if k2.endswith("_census"):
            cen = v
    if cen:
        d["realized_census_nodes"] = cen["nodes_total"]
        d["census_kappa"] = cen["kappa"]
        d["node_ratio_vs_formula"] = cen.get("node_ratio_vs_formula")
        d["node_ratio_vs_prediction"] = (
            round(cen["nodes_total"] / cost["bottom_nodes"], 4)
            if cost.get("bottom_nodes") else None)
    if res.get("wall_s") is not None:
        d["realized_wall_s"] = res["wall_s"]
        d["wall_ratio"] = round(res["wall_s"]
                                / max(cost["wall_estimate_s"], 1e-9), 3)
    d["cpu_s"] = res.get("cpu_s")
    # closure per PROTOCOL (iii): certificate floor at priced W within
    # 10x wall estimate
    q_target = cost["W"] + (2 if row.get("parity_ok") else 1)
    if out == "CERTIFIED_FLOOR":
        d["floor"] = res.get("floor")
        d["closed"] = bool(res.get("floor", 0) >= q_target
                           and res.get("wall_s", 1e18)
                           <= BUDGET_MULT * cost["wall_estimate_s"])
        d["closed_any_time"] = bool(res.get("floor", 0) >= q_target)
    elif out == "COUNTEREXAMPLE":
        d["closed"] = False
        d["counterexample_weight"] = res["counterexample"]["weight"]
        d["counterexample_support"] = res["counterexample"]["support"]
    else:
        d["closed"] = False
        d["stop_reason"] = res.get("reason")
    # criterion (i): structural laws re-measured by the engine
    viol = []
    for r in res.get("rungs", []):
        if not r.get("rank_law_holds", True):
            viol.append({"rung": r["rung"], "law": "rank",
                         "rank_p": r["rank_p"], "rank_tau": r["rank_tau"],
                         "k_cover": r["k_cover"]})
        if r.get("exactness_at_cover_holds") is False:
            viol.append({"rung": r["rung"], "law": "exact_cover"})
    d["structural_violations"] = viol
    return d


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--only", default=None,
                    help="comma list of strata or instance_ids")
    ap.add_argument("--limit", type=int, default=None)
    ap.add_argument("--redo", default=None,
                    help="comma list of instance_ids to re-run even if "
                         "recorded (append-only; the scorecard uses the "
                         "last record per (row, question))")
    args = ap.parse_args()

    LOGS.mkdir(exist_ok=True)
    WORK.mkdir(exist_ok=True)
    preds = load_jsonl(DTT / "predictions.jsonl")
    preds2 = load_jsonl(DTT / "predictions_batch2.jsonl")
    for r in preds:
        r["batch"] = 1
    rows = preds + [r for r in preds2 if r.get("row_class") == "prediction"]

    done = set()
    for r in load_jsonl(RESULTS):
        done.add((r["instance_id"], r.get("q")))
    if args.redo:
        redo = set(args.redo.split(","))
        done = {(i, q) for (i, q) in done if i not in redo}

    # ordering: anchors (burn-in) -> scope-control -> frontier ->
    # parity-strip -> batch-2 predictions
    stratum_order = {"anchor": 0, "scope-control": 1, "frontier": 2,
                     "parity-strip": 3, "order144-depth4": 4}
    rows.sort(key=lambda r: (stratum_order.get(r["stratum"], 9),
                             r["instance_id"]))
    if args.only:
        keys = set(args.only.split(","))
        rows = [r for r in rows
                if r["stratum"] in keys or r["instance_id"] in keys]
    if args.limit:
        rows = rows[: args.limit]

    n_done = 0
    t_batch = time.time()
    for ri, row in enumerate(rows):
        if free_mb() < 250:
            print("DISK GUARD < 250 MB — stopping", flush=True)
            break
        iid = row["instance_id"]
        stratum = row["stratum"]
        if stratum == "scope-control":
            if (iid, "refusal") in done:
                continue
            ok_odd = row["ell"] % 2 == 1 and row["m"] % 2 == 1
            rec = {
                "instance_id": iid, "stratum": stratum, "q": "refusal",
                "kind": "refusal-check",
                "group": row.get("group"),
                "both_axes_odd": ok_odd,
                "free_z2_deck_possible": not ok_odd,
                "refusal_confirmed": ok_odd,
                "note": ("|G| odd: no index-2 subgroup along either axis; "
                         "no free Z2 deck exists — machinery REFUSES "
                         "(frozen prediction: REFUSED)" if ok_odd else
                         "UNEXPECTED: an even axis exists"),
                "ts": time.strftime("%Y-%m-%dT%H:%M:%S"),
            }
            emit(rec)
            n_done += 1
            print(f"[{ri+1}/{len(rows)}] {iid[:12]} scope-control "
                  f"refusal_confirmed={ok_odd}", flush=True)
            continue

        costs = row.get("costs") or []
        # dedupe identical (W) questions
        seen_w = {}
        for qi, cost in enumerate(costs):
            if (iid, f"q{qi}") in done:
                continue
            if cost["W"] in seen_w:
                rec = {"instance_id": iid, "stratum": stratum,
                       "q": f"q{qi}", "kind": "dup-of",
                       "dup_of": seen_w[cost["W"]],
                       "question": cost["question"], "W": cost["W"]}
                emit(rec)
                continue
            seen_w[cost["W"]] = f"q{qi}"
            res = run_question(row, qi, cost)
            ev = eval_question(row, cost, res)
            rec = {
                "instance_id": iid, "stratum": stratum,
                "batch": row.get("batch", 2), "q": f"q{qi}",
                "kind": "closure",
                "group": row.get("group"), "n": row.get("n"),
                "k": row.get("k"), "v2_depth": row.get("v2_depth"),
                "operative": bool(cost is (row.get("cost_operative") or {})
                                  or cost["W"] == (row.get("cost_operative")
                                                   or {}).get("W")),
                "eval": ev,
                "engine": {k2: res.get(k2) for k2 in
                           ("outcome", "floor", "reason", "wall_s",
                            "cpu_s", "budget_s", "maxrss_mb",
                            "top_cycles_all_even", "planted_found")},
                "rungs_measured": res.get("rungs"),
                "stages": res.get("stages"),
                "counterexample": res.get("counterexample"),
                "ts": time.strftime("%Y-%m-%dT%H:%M:%S"),
            }
            emit(rec)
            n_done += 1
            flag = ""
            if ev.get("structural_violations"):
                flag = "  !! STRUCTURAL VIOLATION — CRITERION (i) EVENT"
            if (ev["outcome"] == "COUNTEREXAMPLE"
                    and stratum == "anchor"
                    and row.get("d_exact_corpus") is not None
                    and ev["counterexample_weight"]
                    < row["d_exact_corpus"]):
                flag += ("  !! CEX BELOW CORPUS d_exact — "
                         "STOP-AND-INVESTIGATE")
            print(f"[{ri+1}/{len(rows)}] {iid[:12]} {stratum} q{qi} "
                  f"W={cost['W']} {cost['verdict']} -> {ev['outcome']}"
                  f"{' floor=' + str(ev.get('floor')) if ev.get('floor') else ''}"
                  f" ({res.get('wall_s', '?')}s / est "
                  f"{cost['wall_estimate_s']}s){flag}",
                  flush=True)
        if (ri + 1) % 10 == 0:
            print(f"--- progress {ri+1}/{len(rows)} rows, {n_done} new "
                  f"records, {time.time()-t_batch:.0f}s elapsed",
                  flush=True)
    print(f"PHASE2 BATCH DONE: {n_done} new records in "
          f"{time.time()-t_batch:.0f}s", flush=True)


if __name__ == "__main__":
    main()
