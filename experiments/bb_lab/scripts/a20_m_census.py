"""A20 Phase 1: census of light X-stabilizers of Y4 = [[144,8,10]] to weight 19.

The (M)@20 dangerous-floor obligation for the [[288,8,20]] tower needs
|b| + 2 m(b) >= 20 for every nonzero b in im d2 = rowspace(H_X(Y4)) with
|b| <= 19.  This script enumerates the census table: every G-translation
orbit class of such b, by CryptoMiniSat (XOR-native: b = H_X^T f is a pure
XOR system) with ascending even-weight bands and full-orbit blocking
clauses.  Augmentation parity => stabilizer weights are even, so bands run
W = 2, 4, ..., 18; within band W, all lighter orbits are already blocked
(previous bands terminated UNSAT = complete), so every model has weight
exactly W (asserted).

Output: data/a20/m_census_classes.jsonl (one class per line: weight, qubit
support of b, check support of f, band, secs; band-complete marker lines
{"band_complete": W}).  Resumable: reload + reblock on restart.  A final
"CENSUS COMPLETE" line means every band UNSAT-terminated: the table is
exhaustive up to translation.

Usage (from experiments/bb_lab):
    uv run scripts/a20_m_census.py                # full run (overnight)
    uv run scripts/a20_m_census.py --bands 2,4,6  # smoke test
Phase 2 (per-class m(b) floors) consumes the JSONL; see A20 note SS5.
"""
import argparse
import json
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

OUT = LAB / "data" / "a20"
JSONL = OUT / "m_census_classes.jsonl"

Y4 = {"frame": (18, 4), "A": "1 + x + x^14*y", "B": "1 + x*y^2 + x^2*y^3"}


def build_code():
    G = AbelianGroup(Y4["frame"])
    A = Poly.from_string(Y4["A"], G)
    B = Poly.from_string(Y4["B"], G)
    checks = bb_check_matrices(A, B)
    return G, checks


def orbit_of(b, ell, m):
    """Distinct G-translates of b, as frozensets of flat qubit indices
    (qubit order: block * |G| + (a * m + bb), matching bb_lab / a20_tau_lift)."""
    N = ell * m
    seen, orbit = set(), []
    supp = [(int(q) // N, divmod(int(q) % N, m)) for q in np.flatnonzero(b)]
    for ta in range(ell):
        for tb in range(m):
            key = frozenset(blk * N + ((a + ta) % ell) * m + (bb + tb) % m
                            for blk, (a, bb) in supp)
            if key not in seen:
                seen.add(key)
                orbit.append(key)
    return orbit


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bands", default="2,4,6,8,10,12,14,16,18")
    ap.add_argument("--max-classes", type=int, default=6000)
    args = ap.parse_args()
    bands = [int(w) for w in args.bands.split(",")]

    G, checks = build_code()
    ell, m = Y4["frame"]
    N = ell * m
    n2 = checks.num_qubits
    HX = checks.H_X % 2

    OUT.mkdir(parents=True, exist_ok=True)
    classes, done_bands = [], set()
    if JSONL.exists():
        for line in JSONL.read_text().splitlines():
            row = json.loads(line)
            if "band_complete" in row:
                done_bands.add(row["band_complete"])
            else:
                classes.append(row)
        print(f"resumed: {len(classes)} classes, bands done {sorted(done_bands)}",
              flush=True)

    log = JSONL.open("a")
    total_t0 = time.time()
    for W in bands:
        if W in done_bands:
            continue
        t_band = time.time()
        pool = IDPool()
        b_vars = [pool.id() for _ in range(n2)]
        f_vars = [pool.id() for _ in range(N)]
        solver = pycryptosat.Solver()
        # b_j = XOR of f_i over rows i with HX[i, j] = 1
        for j in range(n2):
            idx = np.flatnonzero(HX[:, j])
            solver.add_xor_clause([b_vars[j]] + [f_vars[i] for i in idx], False)
        solver.add_clause(b_vars)  # b != 0
        card = CardEnc.atmost(lits=b_vars, bound=W, vpool=pool,
                              encoding=EncType.seqcounter)
        for cl in card.clauses:
            solver.add_clause(cl)
        for row in classes:  # reblock everything found so far
            b = np.zeros(n2, dtype=np.uint8)
            b[row["b_support"]] = 1
            for orb in orbit_of(b, ell, m):
                solver.add_clause(
                    [-b_vars[j] if j in orb else b_vars[j] for j in range(n2)])

        found_this_band = 0
        while True:
            if len(classes) >= args.max_classes:
                print(f"CAP HIT at {len(classes)} classes — census INCOMPLETE",
                      flush=True)
                log.close()
                return
            t0 = time.time()
            sat, model = solver.solve()
            dt = time.time() - t0
            if not sat:
                print(f"band W={W}: complete — {found_this_band} classes "
                      f"({time.time() - t_band:.0f}s)", flush=True)
                log.write(json.dumps({"band_complete": W}) + "\n")
                log.flush()
                break
            b = np.array([1 if model[v] else 0 for v in b_vars], dtype=np.uint8)
            f = np.array([1 if model[v] else 0 for v in f_vars], dtype=np.uint8)
            assert (((HX.T @ f) % 2) == b).all(), "model violates b = HX^T f"
            w = int(b.sum())
            assert w == W, f"weight {w} != band {W} (completeness hole?)"
            row = {"w": w, "b_support": np.flatnonzero(b).tolist(),
                   "f_support": np.flatnonzero(f).tolist(),
                   "band": W, "secs": round(dt, 2)}
            classes.append(row)
            log.write(json.dumps(row) + "\n")
            log.flush()
            found_this_band += 1
            for orb in orbit_of(b, ell, m):
                solver.add_clause(
                    [-b_vars[j] if j in orb else b_vars[j] for j in range(n2)])
            if found_this_band % 20 == 0:
                print(f"  W={W}: {found_this_band} classes "
                      f"({len(classes)} total, last solve {dt:.1f}s)",
                      flush=True)

    hist = {}
    for row in classes:
        hist[row["w"]] = hist.get(row["w"], 0) + 1
    print(f"\nCENSUS COMPLETE in {time.time() - total_t0:.0f}s: "
          f"{len(classes)} translation-orbit classes with |b| <= "
          f"{max(bands)}", flush=True)
    print("weight histogram:", json.dumps(hist, sort_keys=True), flush=True)
    log.close()


if __name__ == "__main__":
    main()
