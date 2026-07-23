"""A19 M12 Phase 2: per-class dangerous-floor certificates (M)@12 over BY's census.

Adapted from a20_m_floors.py (A20 note SS5) for the Bravyi-360 y-deck:
base BY = (30,3) [[180,8,12]], cover C = (30,6) [[360,12,<=24]].

For each census class b (a light X-stabilizer of BY with |b| <= 11, from
a19_m_census.py), certify by constrained UNSAT that every NONTRIVIAL cover
X-logical v of C with pushforward p_y(v) = b has weight >= 12. Slice
identity: |v| = |b| + 2*(off-support doubly-occupied count), so per class:

    v in ker H_Z(C)                     (cover X-cycle; XOR rows)
    p_y(v) = b                          (180 fiber XORs: v_q0 + v_q1 = b_q)
    v nontrivial                        (OR over signatures <L_i, v>)
    sum of v over OFF-supp(b) fibers <= 2*(m_req - 1)

UNSAT  ==>  m(b) >= m_req = ceil((12 - |b|)/2)  ==>  |v| >= 12 over b.
SAT    ==>  an explicit sub-12 nontrivial logical of C (would beat the ISD
census floor) — recorded, never expected.

The b = 0 stratum needs no query: p_y(v) = 0 chain-level forces v = tau(u)
sigma-symmetric with u a BY cycle, |v| = 2|u|; |v| <= 10 would need a nonzero
BY cycle of weight <= 5 — dead by the A16 class certificate (mu >= 6).
Together with d(BY) = 12 (projection bound on classes with p_y*[v] != 0),
all-UNSAT here completes certified d(C) >= 12.

Usage (from experiments/bb_lab):
    uv run scripts/a19_m_floors.py [--jobs 6] [--max-weight 10]
Resumable: results append to data/a19/m12_floors_results.jsonl.
"""
import argparse
import json
import multiprocessing as mp
import sys
import time
from pathlib import Path

import numpy as np

LAB = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(LAB / "src"))

import pycryptosat
from pysat.card import CardEnc, EncType
from pysat.formula import IDPool

from bb_lab.checks import bb_check_matrices
from bb_lab.group import AbelianGroup
from bb_lab.poly import Poly
from bb_lab.sat_distance import find_logical_z

OUT = LAB / "data" / "a19"
CENSUS = OUT / "m12_census_classes.jsonl"
RESULTS = OUT / "m12_floors_results.jsonl"

BASE = {"frame": (30, 3), "A": "x^9 + y + y^2", "B": "1 + x^25 + x^26"}
COVER = {"frame": (30, 6), "A": "x^9 + y + y^2", "B": "y^3 + x^25 + x^26"}

TARGET = 12

_CTX = {}


def build_context():
    """Shared matrices, built once per worker (spawn-safe via initializer)."""
    Gc = AbelianGroup(COVER["frame"])
    cc = bb_check_matrices(Poly.from_string(COVER["A"], Gc),
                           Poly.from_string(COVER["B"], Gc))
    L_Z = find_logical_z(cc) % 2
    # fiber map: BY qubit (blk, a*3 + b) -> C qubits (blk, a*6 + b), +3
    Nb, Nc = 90, 180
    fiber = []
    for blk in (0, 1):
        for a in range(30):
            for b in range(3):
                fiber.append((blk * Nc + a * 6 + b, blk * Nc + a * 6 + b + 3))
    assert len(fiber) == 2 * Nb
    _CTX.update(HZc=cc.H_Z % 2, L_Z=L_Z, fiber=fiber)


def floor_query(task):
    """(class_row, m_req) -> result dict. UNSAT = floor certified."""
    row, m_req = task
    if not _CTX:
        build_context()
    HZc, L_Z, fiber = _CTX["HZc"], _CTX["L_Z"], _CTX["fiber"]
    nc = HZc.shape[1]
    bb = np.zeros(len(fiber), dtype=np.uint8)
    bb[row["b_support"]] = 1

    t0 = time.time()
    pool = IDPool()
    v = [pool.id() for _ in range(nc)]
    solver = pycryptosat.Solver()
    for hz_row in HZc:                      # cover cycle
        idx = np.flatnonzero(hz_row)
        if idx.size:
            solver.add_xor_clause([v[i] for i in idx], False)
    off_vars = []
    for qb, (q0, q1) in enumerate(fiber):   # pushforward pinned to b
        solver.add_xor_clause([v[q0], v[q1]], bool(bb[qb]))
        if not bb[qb]:
            off_vars.extend((v[q0], v[q1]))
    sig_vars = []
    for L in L_Z:                           # nontriviality
        idx = np.flatnonzero(L)
        s = pool.id()
        solver.add_xor_clause([s] + [v[i] for i in idx], False)
        sig_vars.append(s)
    solver.add_clause(sig_vars)
    bound = 2 * (m_req - 1)
    if bound < len(off_vars):
        card = CardEnc.atmost(lits=off_vars, bound=bound, vpool=pool,
                              encoding=EncType.seqcounter)
        for cl in card.clauses:
            solver.add_clause(cl)
    sat, model = solver.solve()
    out = {"w": row["w"], "b_support": row["b_support"], "m_req": m_req,
           "bound": bound, "secs": round(time.time() - t0, 2)}
    if not sat:
        out["verdict"] = "UNSAT"            # floor certified for this class
    else:
        vv = np.array([1 if model[x] else 0 for x in v], dtype=np.uint8)
        assert not (HZc @ vv % 2).any()
        out["verdict"] = "SAT"              # sub-12 logical?! record it
        out["witness_weight"] = int(vv.sum())
        out["witness"] = np.flatnonzero(vv).tolist()
    return out


def main():
    global TARGET, CENSUS, RESULTS
    ap = argparse.ArgumentParser()
    ap.add_argument("--jobs", type=int, default=6)
    ap.add_argument("--max-weight", type=int, default=10)
    ap.add_argument("--target", type=int, default=12)
    ap.add_argument("--census", default=None)
    ap.add_argument("--results", default=None)
    args = ap.parse_args()
    TARGET = args.target
    if args.census:
        CENSUS = OUT / args.census
    if args.results:
        RESULTS = OUT / args.results

    done = set()
    if RESULTS.exists():
        for line in RESULTS.read_text().splitlines():
            r = json.loads(line)
            done.add(tuple(r["b_support"]))
        print(f"resumed: {len(done)} classes already certified", flush=True)

    tasks = []
    for line in CENSUS.read_text().splitlines():
        r = json.loads(line)
        if "w" not in r or r["w"] > args.max_weight:
            continue
        if tuple(r["b_support"]) in done:
            continue
        m_req = -(-(TARGET - r["w"]) // 2)
        tasks.append((r, m_req))
    tasks.sort(key=lambda t: -t[0]["w"])    # cheapest strata (high w) first
    print(f"{len(tasks)} floor queries queued "
          f"(strata: {sorted({t[0]['w'] for t in tasks}, reverse=True)})",
          flush=True)

    log = RESULTS.open("a")
    t0, ok, sat_hits = time.time(), 0, []
    ctx = mp.get_context("spawn")
    with ctx.Pool(args.jobs, initializer=build_context) as pool_:
        for out in pool_.imap_unordered(floor_query, tasks):
            log.write(json.dumps(out) + "\n")
            log.flush()
            ok += 1
            if out["verdict"] == "SAT":
                sat_hits.append(out)
                print(f"!! SAT at w={out['w']}: sub-{TARGET} logical weight "
                      f"{out['witness_weight']} — ISD floor beaten?!",
                      flush=True)
            print(f"  [{ok}/{len(tasks)}] w={out['w']} bound={out['bound']}"
                  f" {out['verdict']} ({out['secs']}s)", flush=True)
    print(f"\n(M)@{TARGET} floors (w <= {args.max_weight}): {ok} queries in "
          f"{time.time() - t0:.0f}s; SAT hits: {len(sat_hits)}", flush=True)
    if not sat_hits:
        print("ALL FLOORS CERTIFIED: every nontrivial cover logical with "
              f"p_y(v) a census-listed b has weight >= {TARGET}.", flush=True)
    log.close()


if __name__ == "__main__":
    main()
