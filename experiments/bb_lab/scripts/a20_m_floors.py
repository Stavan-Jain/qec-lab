"""A20 Phase 2: per-class dangerous-floor certificates (M)@20 over Y4's census.

For each census class b (a light X-stabilizer of Y4 = [[144,8,10]], from
a20_m_census.py), certify by constrained UNSAT that every NONTRIVIAL cover
X-logical v of Y8 = [[288,8,20]] with pushforward p(v) = b has weight >= 20.
By the slice identity |v| = |b| + 2*(off-support doubly-occupied count), the
query per class is:

    v in ker H_Z(Y8)                    (cover X-cycle; XOR rows)
    p(v) = b                            (144 fiber XORs: v_q0 + v_q1 = b_q)
    v nontrivial                        (OR over signatures <L_i, v> = sig_i,
                                         L_i = Z-logical basis of Y8)
    sum of v over OFF-supp(b) fibers <= 2*(m_req - 1)

UNSAT  ==>  m(b) >= m_req = ceil((20 - |b|)/2)  ==>  |v| >= 20 over b.
SAT    ==>  an explicit sub-20 nontrivial logical of Y8 (would refute IBM's
MILP-exact d = 20) — recorded, never expected.

Strata (completed census bands only, w <= 16), cheapest first:
    w=16: 375 classes, bound 2      w=14: 54, bound 4
    w=12: 33, bound 6               w=10: 6, bound 8
    w=6 (hexagon): 1, bound 12      [w=18 stratum: Phase 2b, uniform m>=1]

Usage (from experiments/bb_lab):
    uv run scripts/a20_m_floors.py [--jobs 6] [--max-weight 16]
Resumable: results append to data/a20/m_floors_results.jsonl.
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

OUT = LAB / "data" / "a20"
CENSUS = OUT / "m_census_classes.jsonl"
RESULTS = OUT / "m_floors_results.jsonl"

Y4 = {"frame": (18, 4), "A": "1 + x + x^14*y", "B": "1 + x*y^2 + x^2*y^3"}
Y8 = {"frame": (18, 8), "A": "1 + x*y^4 + x^14*y", "B": "1 + x*y^2 + x^2*y^7"}

_CTX = {}


def build_context():
    """Shared matrices, built once per worker (spawn-safe via initializer)."""
    G8 = AbelianGroup(Y8["frame"])
    c8 = bb_check_matrices(Poly.from_string(Y8["A"], G8),
                           Poly.from_string(Y8["B"], G8))
    L_Z = find_logical_z(c8) % 2
    # fiber map: Y4 qubit (blk, a*4 + b) -> Y8 qubits (blk, a*8 + b), +4
    N4, N8 = 72, 144
    fiber = []
    for blk in (0, 1):
        for a in range(18):
            for b in range(4):
                fiber.append((blk * N8 + a * 8 + b, blk * N8 + a * 8 + b + 4))
    assert len(fiber) == 2 * N4
    _CTX.update(HZ8=c8.H_Z % 2, HX8=c8.H_X % 2, L_Z=L_Z, fiber=fiber)


def floor_query(task):
    """(class_row, m_req) -> result dict. UNSAT = floor certified."""
    row, m_req = task
    if not _CTX:
        build_context()
    HZ8, L_Z, fiber = _CTX["HZ8"], _CTX["L_Z"], _CTX["fiber"]
    n8 = HZ8.shape[1]
    b4 = np.zeros(len(fiber), dtype=np.uint8)
    b4[row["b_support"]] = 1

    t0 = time.time()
    pool = IDPool()
    v = [pool.id() for _ in range(n8)]
    solver = pycryptosat.Solver()
    for hz_row in HZ8:                      # cover cycle
        idx = np.flatnonzero(hz_row)
        if idx.size:
            solver.add_xor_clause([v[i] for i in idx], False)
    off_vars = []
    for q4, (q0, q1) in enumerate(fiber):   # pushforward pinned to b
        solver.add_xor_clause([v[q0], v[q1]], bool(b4[q4]))
        if not b4[q4]:
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
        assert not (HZ8 @ vv % 2).any()
        out["verdict"] = "SAT"              # sub-20 logical?! record it
        out["witness_weight"] = int(vv.sum())
        out["witness"] = np.flatnonzero(vv).tolist()
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--jobs", type=int, default=6)
    ap.add_argument("--max-weight", type=int, default=16)
    args = ap.parse_args()

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
        m_req = -(-(20 - r["w"]) // 2)
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
                print(f"!! SAT at w={out['w']}: sub-20 logical weight "
                      f"{out['witness_weight']} — IBM d=20 REFUTED?!",
                      flush=True)
            if ok % 25 == 0 or out["w"] <= 12:
                print(f"  [{ok}/{len(tasks)}] w={out['w']} bound={out['bound']}"
                      f" {out['verdict']} ({out['secs']}s)", flush=True)
    print(f"\nPhase 2 (w <= {args.max_weight}): {ok} queries in "
          f"{time.time() - t0:.0f}s; SAT hits: {len(sat_hits)}", flush=True)
    if not sat_hits:
        print("ALL FLOORS CERTIFIED on completed strata: every nontrivial "
              "cover logical over a census-listed b with |b| <= "
              f"{args.max_weight} has weight >= 20.", flush=True)
    log.close()


if __name__ == "__main__":
    main()
