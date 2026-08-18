"""Phase-3 Arm B: the equal-compute SAT-only arm (PROTOCOL criterion (v)).

Paired design on the same cohort rows Arm A (the descent/certificate
pipeline) ran in Stage A.  Per row, Arm B gets the SAME CPU budget Arm A
actually spent on that row's closure questions, floored at 60 s
(PROTOCOL (v) verbatim), on the corpus production SAT lane
(bb_lab.sat_distance CMS ladder = `sat-cms-ladder (pycryptosat)`, the
stack that closed the corpus rows this cohort was drawn from; the
tandem-maxsat second stage is not wired here — recorded as an Arm-B
limitation in the scorecard).  d_ub (corpus/L1, available to both arms)
caps the ladder.  Scheduling: cheapest-first (ascending n, then d_ub).

Score per row: exact d within budget, or honest UNSAT-prefix floor;
closure for (v) = exact d OR floor >= the row's operative window top.

Appends to ../phase3_armB.jsonl (resume-safe).  Two workers, nice'd.
"""

from __future__ import annotations

import argparse
import json
import os
import signal
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

HERE = Path(__file__).resolve().parent
DTT = HERE.parent
LAB = DTT.parent.parent
OUT = DTT / "phase3_armB.jsonl"
PY = sys.executable

MIN_BUDGET_S = 60.0


def load_jsonl(p: Path) -> list[dict]:
    return [json.loads(line) for line in p.open()] if p.exists() else []


def collect_rows() -> list[dict]:
    preds = load_jsonl(DTT / "predictions.jsonl")
    preds2 = [r for r in load_jsonl(DTT / "predictions_batch2.jsonl")
              if r.get("row_class") == "prediction"]
    rows = {r["instance_id"]: r for r in preds + preds2
            if r.get("stratum") != "scope-control"}
    # Arm A realized CPU per row (closure records only)
    cpu: dict[str, float] = {}
    target: dict[str, int] = {}
    for rec in load_jsonl(DTT / "phase2_results.jsonl"):
        if rec.get("kind") != "closure":
            continue
        iid = rec["instance_id"]
        cpu[iid] = cpu.get(iid, 0.0) + float(rec["eval"].get("cpu_s") or 0)
    out = []
    for iid, r in rows.items():
        if iid not in cpu:
            continue                      # arm A never ran it -> unpaired
        co = r.get("cost_operative") or {}
        w_top = (co.get("W") or 0) + (2 if r.get("parity_ok") else 1)
        out.append({
            "instance_id": iid, "stratum": r["stratum"],
            "group": r.get("group"), "n": r.get("n"), "k": r.get("k"),
            "ell": r["ell"], "m": r["m"],
            "A": r["A_poly"], "B": r["B_poly"],
            "d_ub": r.get("d_ub"),
            "armA_cpu_s": round(cpu[iid], 2),
            "budget_s": round(max(cpu[iid], MIN_BUDGET_S), 2),
            "window_top": w_top,
        })
    out.sort(key=lambda x: (x["n"], x.get("d_ub") or 99, x["instance_id"]))
    return out


def run_row(row: dict) -> dict:
    tag = row["instance_id"][:12]
    prog = HERE / "work" / f"armB_{tag}.prog"
    prog.unlink(missing_ok=True)
    cmd = [PY, str(HERE / "sat_worker.py"),  # parent is nice 5; worker renices itself
           "--ell", str(row["ell"]), "--m", str(row["m"]),
           "--A", row["A"], "--B", row["B"], "--progress", str(prog)]
    if row.get("d_ub"):
        cmd += ["--wmax", str(row["d_ub"])]
    budget = row["budget_s"]
    t0 = time.time()
    res: dict = {}
    try:
        p = subprocess.Popen(cmd, stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE, text=True,
                             start_new_session=True, cwd=str(LAB))
        try:
            so, se = p.communicate(timeout=budget * 1.15 + 5)
            if p.returncode == 0 and "DISTANCE" in so:
                res["outcome"] = "exact"
                res["d"] = int(so.split("DISTANCE")[1].split()[0])
            else:
                res["outcome"] = "error"
                res["error"] = (se.strip().splitlines() or ["?"])[-1][:200]
        except subprocess.TimeoutExpired:
            try:
                os.killpg(p.pid, signal.SIGKILL)
            except ProcessLookupError:
                pass
            p.wait()
            res["outcome"] = "timeout"
    finally:
        res["wall_s"] = round(time.time() - t0, 1)
    rungs = []
    if prog.exists():
        for line in prog.read_text().splitlines():
            w, s, dt = line.split(",")
            rungs.append([int(w), s, float(dt)])
        prog.unlink()
    res["rungs"] = rungs
    floor_w = 0
    for w, s, _ in rungs:
        if s == "UNSAT" and w == floor_w + 1:
            floor_w = w
        else:
            break
    if floor_w:
        res["floor"] = floor_w + 1        # solver-proved d >= floor_w + 1
    # (v) closure: exact d, or floor >= operative window top
    res["closed_v"] = bool(
        res.get("outcome") == "exact"
        or (res.get("floor") or 0) >= row["window_top"])
    res.update({k: row[k] for k in
                ("instance_id", "stratum", "group", "n", "k", "d_ub",
                 "armA_cpu_s", "budget_s", "window_top")})
    res["method"] = "sat-cms-ladder (pycryptosat via bb_lab.sat_distance)"
    res["ts"] = time.strftime("%Y-%m-%dT%H:%M:%S")
    return res


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=2)
    ap.add_argument("--limit", type=int, default=None)
    args = ap.parse_args()
    (HERE / "work").mkdir(exist_ok=True)
    done = {r["instance_id"] for r in load_jsonl(OUT)
            if not r.get("superseded")}
    rows = [r for r in collect_rows() if r["instance_id"] not in done]
    if args.limit:
        rows = rows[: args.limit]
    print(f"armB: {len(rows)} rows to run ({len(done)} done); "
          f"total budget {sum(r['budget_s'] for r in rows)/3600:.2f} h",
          flush=True)
    t0 = time.time()
    n = 0
    with ThreadPoolExecutor(max_workers=args.workers) as ex:
        futs = {ex.submit(run_row, r): r for r in rows}
        for fut in as_completed(futs):
            rec = fut.result()
            with OUT.open("a") as f:
                f.write(json.dumps(rec) + "\n")
            n += 1
            if n % 10 == 0 or rec.get("outcome") != "exact":
                print(f"[{n}/{len(rows)}] {rec['instance_id'][:12]} "
                      f"n={rec['n']} {rec['outcome']} "
                      f"d={rec.get('d')} floor={rec.get('floor')} "
                      f"closed_v={rec['closed_v']} "
                      f"({rec['wall_s']}s/{rec['budget_s']}s) "
                      f"[batch {time.time()-t0:.0f}s]", flush=True)
    print(f"ARMB DONE {n} rows in {time.time()-t0:.0f}s", flush=True)


if __name__ == "__main__":
    main()
