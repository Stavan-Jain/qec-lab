"""A20 V7: fiber-pinned UNSAT floors for band-18 classes the SAT census missed.

Input: data/a20/v7_band18_new.jsonl (from a20_v7_completeness.py) — the
weight-18 boundary classes of Y4 present in the analytic complete list but
absent from the (unterminated) SAT census band.  Runs the identical
per-class query as a20_m_floors.py (imported, read-only): for each class b,
m_req = ceil((20 - 18)/2) = 1, i.e. cover cycles v of Y8 = [[288,8,20]] with
pushforward b and ZERO doubly-occupied off-support fibers; UNSAT ==> every
nontrivial cover logical over b has weight >= 20.

Together with the 1,655 classes already certified in m_floors_results.jsonl,
UNSAT on all rows here completes DangerousFloorNZ 20 over the DEFINITIVE
(V7-complete) class list.

Usage: cd experiments/bb_lab && uv run python scripts/a20_v7_new_floors.py [--jobs 6]
Writes data/a20/v7_floors_results.jsonl (append-only, resumable; kept
separate from m_floors_results.jsonl, which parallel sessions may extend).
"""
import argparse
import json
import multiprocessing as mp
import sys
import time
from pathlib import Path

LAB = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(LAB / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import a20_m_floors as mf                     # noqa: E402  (read-only reuse)

NEW = mf.OUT / "v7_band18_new.jsonl"
RESULTS = mf.OUT / "v7_floors_results.jsonl"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--jobs", type=int, default=6)
    args = ap.parse_args()

    done = set()
    if RESULTS.exists():
        for line in RESULTS.read_text().splitlines():
            done.add(tuple(json.loads(line)["b_support"]))
        print(f"resumed: {len(done)} classes already certified", flush=True)

    tasks = []
    for line in NEW.read_text().splitlines():
        r = json.loads(line)
        if tuple(r["b_support"]) in done:
            continue
        m_req = -(-(20 - r["w"]) // 2)
        tasks.append((r, m_req))
    print(f"{len(tasks)} new-class floor queries queued", flush=True)

    log = RESULTS.open("a")
    t0, ok, sat_hits = time.time(), 0, []
    ctx = mp.get_context("spawn")
    with ctx.Pool(args.jobs, initializer=mf.build_context) as pool:
        for out in pool.imap_unordered(mf.floor_query, tasks):
            log.write(json.dumps(out) + "\n")
            log.flush()
            ok += 1
            if out["verdict"] == "SAT":
                sat_hits.append(out)
                print(f"!! SAT: sub-20 logical weight {out['witness_weight']} "
                      f"— IBM d=20 REFUTED?!", flush=True)
            if ok % 100 == 0:
                print(f"  [{ok}/{len(tasks)}] ({time.time()-t0:.0f}s)",
                      flush=True)
    log.close()
    print(f"\ndone: {ok} queries in {time.time()-t0:.0f}s; "
          f"SAT hits: {len(sat_hits)}", flush=True)
    if not sat_hits:
        print("ALL NEW-CLASS FLOORS CERTIFIED: combined with "
              "m_floors_results.jsonl (1,655 classes), DangerousFloorNZ 20 "
              "holds over the definitive V7-complete class list.", flush=True)


if __name__ == "__main__":
    main()
