"""Phase-3 ground-truth arm (PROTOCOL criterion (iv)).

--select: choose ~N rows stratified over v2-depth x stratum x Stage-A
outcome lane, prioritizing (a) certificate-tier Stage-A results
(CERTIFIED_FLOOR or COUNTEREXAMPLE-exact) with SAT-plausible parameters
(n <= 168, or n = 288 with small W), and (b) a spread of rows Stage A
left bounded-only (ENVELOPE/BLOWUP stops).  Writes phase3_gt_selection.json.

--run: independent SAT to exactness (CMS ladder, no window hints beyond
d_ub) with a generous per-row cap; <= 2 concurrent, nice'd.  Appends to
../phase3_groundtruth.jsonl with the agreement verdict per row:

  certificate floor F  vs solver d:  AGREE iff d >= F  (floor sound)
  counterexample w    vs solver d:  AGREE iff d == w  (census-exactness)

ANY disagreement is a stop-and-investigate event (flagged loudly).
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
OUT = DTT / "phase3_groundtruth.jsonl"
SEL = DTT / "phase3_gt_selection.json"
PY = sys.executable

CAP_S = 3600.0


def load_jsonl(p: Path) -> list[dict]:
    return [json.loads(line) for line in p.open()] if p.exists() else []


def stage_a_rows() -> dict[str, dict]:
    preds = load_jsonl(DTT / "predictions.jsonl")
    preds2 = [r for r in load_jsonl(DTT / "predictions_batch2.jsonl")
              if r.get("row_class") == "prediction"]
    rows = {r["instance_id"]: r for r in preds + preds2
            if r.get("stratum") != "scope-control"}
    per: dict[str, dict] = {}
    for rec in load_jsonl(DTT / "phase2_results.jsonl"):
        if rec.get("kind") != "closure":
            continue
        iid = rec["instance_id"]
        e = per.setdefault(iid, {"floor": None, "cex_w": None,
                                 "outcomes": []})
        ev = rec["eval"]
        e["outcomes"].append(ev["outcome"])
        if ev["outcome"] == "CERTIFIED_FLOOR" and ev.get("floor"):
            e["floor"] = max(e["floor"] or 0, ev["floor"])
        if ev["outcome"] == "COUNTEREXAMPLE":
            w = ev["counterexample_weight"]
            e["cex_w"] = min(e["cex_w"] or 10**9, w)
    for iid, e in per.items():
        r = rows.get(iid)
        if r is None:
            continue
        e.update({"stratum": r["stratum"], "n": r["n"], "k": r["k"],
                  "group": r.get("group"), "v2": r.get("v2_depth"),
                  "ell": r["ell"], "m": r["m"], "A": r["A_poly"],
                  "B": r["B_poly"], "d_ub": r.get("d_ub"),
                  "d_exact_corpus": r.get("d_exact_corpus"),
                  "W_op": (r.get("cost_operative") or {}).get("W")})
    return per


def select(n_target: int = 20) -> list[dict]:
    per = stage_a_rows()
    cands = []
    for iid, e in per.items():
        if "stratum" not in e or e["stratum"] == "anchor":
            continue                      # anchors are corpus-exact already
        lane = ("cex-exact" if e["cex_w"] is not None else
                "cert-floor" if e["floor"] is not None else "bounded-only")
        sat_ok = e["n"] <= 168 or (e["n"] <= 288
                                   and (e["cex_w"] or e["W_op"] or 99) <= 12)
        cands.append({"instance_id": iid, "lane": lane,
                      "sat_plausible": sat_ok, **e})
    # stratify: (v2, stratum, lane) cells, round-robin fill
    cells: dict[tuple, list[dict]] = {}
    for c in cands:
        if not c["sat_plausible"]:
            continue
        cells.setdefault((c["v2"], c["stratum"], c["lane"]), []).append(c)
    for cell in cells.values():
        cell.sort(key=lambda c: (c["n"], c["instance_id"]))
    sel: list[dict] = []
    rounds = 0
    while len(sel) < n_target and rounds < 50:
        for key in sorted(cells):
            if cells[key] and len(sel) < n_target:
                sel.append(cells[key].pop(0))
        rounds += 1
    return sel


def run_row(row: dict) -> dict:
    tag = row["instance_id"][:12]
    prog = HERE / "work" / f"gt_{tag}.prog"
    prog.unlink(missing_ok=True)
    cmd = [PY, str(HERE / "sat_worker.py"),  # parent is nice 5; worker renices itself
           "--ell", str(row["ell"]), "--m", str(row["m"]),
           "--A", row["A"], "--B", row["B"], "--progress", str(prog)]
    if row.get("d_ub"):
        cmd += ["--wmax", str(row["d_ub"])]
    t0 = time.time()
    res: dict = {}
    try:
        p = subprocess.Popen(cmd, stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE, text=True,
                             start_new_session=True, cwd=str(LAB))
        try:
            so, se = p.communicate(timeout=CAP_S)
            if p.returncode == 0 and "DISTANCE" in so:
                res["outcome"] = "exact"
                res["d_sat"] = int(so.split("DISTANCE")[1].split()[0])
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
        res["sat_floor"] = floor_w + 1
    # agreement (criterion iv)
    agree = None
    notes = []
    d = res.get("d_sat")
    if d is not None:
        if row.get("floor") is not None:
            ok = d >= row["floor"]
            agree = ok if agree is None else (agree and ok)
            notes.append(f"cert floor {row['floor']} vs d {d}: "
                         f"{'AGREE' if ok else 'DISAGREE'}")
        if row.get("cex_w") is not None:
            ok = d == row["cex_w"]
            agree = ok if agree is None else (agree and ok)
            notes.append(f"census-exact {row['cex_w']} vs d {d}: "
                         f"{'AGREE' if ok else 'DISAGREE'}")
        if row.get("floor") is None and row.get("cex_w") is None:
            notes.append(f"bounded-only row: solver d = {d} (new info)")
    else:
        if row.get("floor") is not None and res.get("sat_floor"):
            # partial consistency: solver floor and cert floor can't cross
            # a verified witness; only comparable when both bound the
            # same side — record both, no verdict
            notes.append(f"timeout: sat_floor {res.get('sat_floor')} "
                         f"(cert floor {row['floor']}); no exactness")
    res["agreement"] = agree
    res["agreement_notes"] = notes
    res.update({k: row.get(k) for k in
                ("instance_id", "stratum", "group", "n", "k", "lane",
                 "floor", "cex_w", "d_ub", "d_exact_corpus", "W_op")})
    res["method"] = "sat-cms-ladder (pycryptosat via bb_lab.sat_distance)"
    res["ts"] = time.strftime("%Y-%m-%dT%H:%M:%S")
    return res


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--select", action="store_true")
    ap.add_argument("--run", action="store_true")
    ap.add_argument("--n", type=int, default=20)
    ap.add_argument("--workers", type=int, default=2)
    args = ap.parse_args()
    (HERE / "work").mkdir(exist_ok=True)
    if args.select:
        sel = select(args.n)
        SEL.write_text(json.dumps(
            {"selected": sel, "n": len(sel),
             "ts": time.strftime("%Y-%m-%dT%H:%M:%S")}, indent=1))
        for c in sel:
            print(f"  {c['instance_id'][:12]} {c['stratum']:<16} v2="
                  f"{c['v2']} {c['lane']:<12} n={c['n']} floor="
                  f"{c['floor']} cex={c['cex_w']} d_ub={c['d_ub']}")
        print(f"selected {len(sel)} -> {SEL}", flush=True)
    if args.run:
        sel = json.loads(SEL.read_text())["selected"]
        done = {r["instance_id"] for r in load_jsonl(OUT)}
        todo = [r for r in sel if r["instance_id"] not in done]
        print(f"groundtruth: {len(todo)} rows (cap {CAP_S:.0f}s each, "
              f"{args.workers} workers)", flush=True)
        t0 = time.time()
        with ThreadPoolExecutor(max_workers=args.workers) as ex:
            futs = {ex.submit(run_row, r): r for r in todo}
            for fut in as_completed(futs):
                rec = fut.result()
                with OUT.open("a") as f:
                    f.write(json.dumps(rec) + "\n")
                flag = ("  !! DISAGREEMENT — STOP AND INVESTIGATE"
                        if rec.get("agreement") is False else "")
                print(f"  {rec['instance_id'][:12]} {rec['outcome']} "
                      f"d_sat={rec.get('d_sat')} vs floor={rec.get('floor')} "
                      f"cex={rec.get('cex_w')} agree={rec.get('agreement')} "
                      f"({rec['wall_s']}s){flag}", flush=True)
        print(f"GROUNDTRUTH DONE in {time.time()-t0:.0f}s", flush=True)


if __name__ == "__main__":
    main()
