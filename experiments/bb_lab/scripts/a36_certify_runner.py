"""A36 — budget-capped certify() runner for constructed doubling cells.

Same budget discipline as the A35 sibling runner (user-set, hard):
300 s per cover code, one retry at 590 s iff the failure was
budget-shaped (total <= ~900 s = the 15-minute wall); predictably heavy
cells (d_base = 12) take a single 850 s shot.  An outer process kill at
budget + 60 s grace enforces the wall against enumeration-chunk
overshoot.  threads = 4 (NOT the front-end default 8): the A35 sibling
session shares this 10-core machine.

Input cells come from the a36 orbit screen's SF-PASS finalists — each a
JSON object {"orders": [L, M], "A": ..., "B": ..., "tag": ...} giving
the COVER spec.  Verdicts append to data/a36/verdicts.jsonl.

Usage:
    uv run python scripts/a36_certify_runner.py --cell '{"orders": [14, 9],
        "A": "...", "B": "...", "tag": "T1p0x_cell..."}' [--budget 300]
    uv run python scripts/a36_certify_runner.py --from-screen \
        data/a36/T1_x_screen.json [--max-cells 2]
"""

from __future__ import annotations

import argparse
import json
import multiprocessing as mp
import time
from pathlib import Path

LAB = Path(__file__).resolve().parents[1]
DATA = LAB / "data" / "a36"
LEDGER = DATA / "verdicts.jsonl"
RUNDIR = DATA / "runs"

BUDGET_FAST = 300.0
BUDGET_RETRY = 590.0
BUDGET_HEAVY = 850.0
GRACE = 60.0
THREADS = 4


def _child(orders, A, B, budget_s, outfile, logfile):
    import sys
    log = open(logfile, "w")
    sys.stdout = sys.stderr = log

    def prog(stage, **kw):
        print(f"[{time.strftime('%H:%M:%S')}] {stage} {kw}", flush=True)

    from bb_lab.doubling_certify import certify, scrub_json
    v = certify(tuple(orders), A, B, budget_s=budget_s, threads=THREADS,
                workdir=Path(outfile).parent, progress=prog)
    Path(outfile).write_text(json.dumps(_slim(scrub_json(v)), indent=1,
                                        default=str))


def _slim(v: dict) -> dict:
    st = v.get("stages", {})
    keep = {}
    for name in ("detect", "base", "d_base", "census", "safe_floor",
                 "rung_pass", "witness", "abort_reason", "candidate_log"):
        if name in st:
            x = st[name]
            if isinstance(x, dict):
                x = {k: w for k, w in x.items()
                     if k not in ("engine", "rows", "witness", "classes")}
            keep[name] = x
    return {**{k: w for k, w in v.items() if k != "stages"}, "stages": keep}


def classify(v: dict | None, killed: bool) -> tuple[str, str]:
    if v is None:
        return ("BUDGET", "outer wall kill" if killed else "child died")
    s, reason = v.get("status", "?"), str(v.get("reason", ""))
    st = v.get("stages", {})
    ab = str(st.get("abort_reason", ""))
    if s in ("CERTIFIED", "FLOOR-ONLY", "DOUBLING-REFUTED"):
        return (s, ab or reason or str(
            v.get("distance", {}).get("statement", "")))
    if "no literal-lift axis doubling with (R)" in reason:
        return ("NO-R", reason)
    if s == "REFUSED":
        return ("SCOPE/ANOMALY", reason)
    joined = " | ".join(x for x in [ab, reason] if x)
    if "nodes > cap" in joined or "beyond the front-end" in joined:
        return ("SCALE-REFUSED", joined)
    return ("BUDGET", joined or "no candidate completed in budget")


def run_cell(cell: dict, budget_s: float) -> dict:
    tag = cell.get("tag") or f"Z{cell['orders'][0]}x{cell['orders'][1]}"
    wd = RUNDIR / f"{tag}_{int(budget_s)}"
    wd.mkdir(parents=True, exist_ok=True)
    outfile, logfile = wd / "verdict_slim.json", wd / "run.log"
    if outfile.exists():
        outfile.unlink()
    ctx = mp.get_context("spawn")
    p = ctx.Process(target=_child, args=(
        tuple(cell["orders"]), cell["A"], cell["B"], budget_s,
        str(outfile), str(logfile)))
    t0 = time.monotonic()
    p.start()
    p.join(budget_s + GRACE)
    killed = False
    if p.is_alive():
        killed = True
        p.terminate()
        p.join(10)
        if p.is_alive():
            p.kill()
            p.join()
    wall = round(time.monotonic() - t0, 1)
    v = None
    if outfile.exists():
        v = json.loads(outfile.read_text())
    status, note = classify(v, killed)
    rec = {
        "tag": tag, "cell": {k: cell[k] for k in ("orders", "A", "B")},
        "budget_s": budget_s, "wall_s": wall, "status": status,
        "note": note[:400],
        "distance": (v or {}).get("distance"),
        "tier": (v or {}).get("tier"),
        "ts": time.strftime("%Y-%m-%d %H:%M:%S"),
    }
    with LEDGER.open("a") as fh:
        fh.write(json.dumps(rec) + "\n")
    print(f"[{rec['ts']}] {tag} budget={budget_s:.0f}s wall={wall}s "
          f"-> {status}  {note[:150]}", flush=True)
    return rec


def run_with_ladder(cell: dict, heavy: bool,
                    spent_s: float = 0.0) -> dict:
    """The budget rule: 300 s, one 590 s retry iff budget-shaped; heavy
    cells (d_base = 12) take a single 850 s shot.  `spent_s` = budget
    already charged to this cover code (e.g. its T1.5 exact safe-floor
    pre-stage) — subtracted so the per-code wall stays inside 5/15."""
    if heavy:
        return run_cell(cell, max(BUDGET_HEAVY - spent_s, 120.0))
    rec = run_cell(cell, max(BUDGET_FAST - spent_s, 120.0))
    if rec["status"] == "BUDGET":
        rec = run_cell(cell, max(BUDGET_RETRY - spent_s, 120.0))
        if rec["status"] == "BUDGET":
            rec["status"] = "BUDGET-FINAL"
            with LEDGER.open("a") as fh:
                fh.write(json.dumps({**rec, "status": "BUDGET-FINAL"})
                         + "\n")
    return rec


def t15_exact_safe_floor(cover: dict, d_base: int,
                         budget_s: float = 120.0) -> dict:
    """T1.5: the front-end's own exact BZ safe-floor decision, run
    directly at base scale (no d_base recompute, no census).  Only
    exact-certified cells deserve a certify() run.

    `cover` = {"orders", "A", "B"} cover spec; the base is recovered by
    detect() (first (R) candidate)."""
    import time as _t

    from bb_lab.doubling_certify import BaseTools, detect, safe_floor
    from bb_lab.group import AbelianGroup
    from bb_lab.poly import Poly

    G = AbelianGroup(tuple(cover["orders"]))
    A = Poly.from_string(cover["A"], G)
    B = Poly.from_string(cover["B"], G)
    cands = [c for c in detect(G, A, B) if c.R_holds]
    if not cands:
        return {"t15": "NO-R"}
    cand = cands[0]
    Gb = AbelianGroup(cand.base_group)
    bt = BaseTools(Gb, Poly.from_string(cand.base_A, Gb),
                   Poly.from_string(cand.base_B, Gb))
    wd = RUNDIR / "t15"
    wd.mkdir(parents=True, exist_ok=True)
    t0 = _t.monotonic()
    try:
        sf = safe_floor(bt, cand.axis, 2 * d_base,
                        _t.monotonic() + budget_s, THREADS, wd)
    except Exception as exc:  # TimeoutError / TimeoutExpired: no verdict
        return {"t15": "TIMEOUT", "err": str(exc)[:120],
                "wall_s": round(_t.monotonic() - t0, 1)}
    return {"t15": "CERTIFIED" if sf["certified"] else "REFUTED",
            "min_refuted": sf.get("min_refuted"),
            "wall_s": round(_t.monotonic() - t0, 1)}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--cell", type=str, help="JSON cover spec")
    ap.add_argument("--from-screen", type=str,
                    help="screen output JSON; runs its SF-PASS finalists")
    ap.add_argument("--max-cells", type=int, default=3)
    ap.add_argument("--budget", type=float, default=None)
    ap.add_argument("--heavy", action="store_true",
                    help="single 850 s shot (d_base = 12 cells)")
    ap.add_argument("--t15-only", action="store_true",
                    help="run the exact safe-floor pre-stage only")
    ap.add_argument("--d-base", type=int, default=None,
                    help="base distance for --t15-only")
    args = ap.parse_args()
    RUNDIR.mkdir(parents=True, exist_ok=True)
    if args.cell:
        cell = json.loads(args.cell)
        if args.t15_only:
            r = t15_exact_safe_floor(cell, args.d_base)
            print(json.dumps({"cell": cell, **r}), flush=True)
            return
        if args.budget:
            run_cell(cell, args.budget)
        else:
            run_with_ladder(cell, args.heavy)
        return
    if args.d_base is None:
        raise SystemExit("--from-screen requires --d-base")
    src = json.loads(Path(args.from_screen).read_text())
    # swap-twin logic: BB(B, A) is BB(A, B) with the qubit blocks
    # swapped — the SAME cover code.  A kill on either orientation kills
    # both; a certify target is the swap-class, deduped.
    killed_pairs = set()
    for r in src["records"]:
        if r.get("stage") == "t1" and not r.get("pass"):
            killed_pairs.add(frozenset((r["A"], r["B"])))
    finals = []
    seen_class: set = set()
    for r in src["records"]:
        if not (r.get("stage") == "t1" and r.get("pass")):
            continue
        cls = frozenset((r["A"], r["B"]))
        if cls in killed_pairs:
            continue  # swap twin was killed: same cover, cross-applied
        if cls in seen_class:
            continue
        seen_class.add(cls)
        finals.append(r)
    print(f"{len(finals)} swap-consistent SF-PASS cover classes in "
          f"{args.from_screen}; T1.5 gate, then up to {args.max_cells} "
          "certify runs", flush=True)
    d_base = args.d_base
    certified_runs = 0
    t15_ledger = DATA / "t15_results.jsonl"
    for r in finals:
        if certified_runs >= args.max_cells:
            break
        key = "cover_A" if "cover_A" in r else "A"
        keyb = "cover_B" if "cover_B" in r else "B"
        keyo = "cover_orders" if "cover_orders" in r else "orders"
        cell = {"orders": r[keyo], "A": r[key], "B": r[keyb]}
        t15 = t15_exact_safe_floor(cell, d_base, budget_s=90.0)
        with t15_ledger.open("a") as fh:
            fh.write(json.dumps({"cell": cell, "point": src["point"],
                                 "pres": r["pres"], "axis": r["axis"],
                                 **t15}) + "\n")
        print(f"  T1.5 {t15.get('t15')} "
              f"(min_refuted={t15.get('min_refuted')}, "
              f"{t15.get('wall_s')}s) A={cell['A']} B={cell['B']}",
              flush=True)
        if t15.get("t15") != "CERTIFIED":
            continue
        tag = (f"{src['point']}p{r['pres']}{r['axis']}_"
               f"s0{r['s0']}_"
               f"{abs(hash(r[key] + r[keyb])) % 10**8:08d}")
        certified_runs += 1
        run_with_ladder({**cell, "tag": tag}, args.heavy,
                        spent_s=t15.get("wall_s", 0.0))


if __name__ == "__main__":
    main()
