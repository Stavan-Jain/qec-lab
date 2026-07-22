"""A19 (M)@24 census, redesigned: lex-leader symmetry breaking (session 2c).

Same XOR system as a19_m_census.py (b = H_X^T f over BY = (30,3), nonzero,
ascending even bands, atmost-W cardinality), but instead of blocking all 90
translates of every found class (the A20-flagged clause-bloat bottleneck),
the solver carries STATIC lex-leader constraints

    b <=_lex T_g(b)   for every nonidentity translation g (89 chains)

so every model IS its orbit's lex-minimal representative; each found class
is then blocked with ONE clause. Completeness: every orbit contains its
lex-min element, which satisfies all chains; the band-ascending invariant
(models in band W have weight exactly W) is preserved since all lighter
canonical reps are blocked and non-canonical vectors are excluded by the
chains.

Validation harness: bands 2-12 must reproduce the a19_m_census.py ground
truth {6: 1, 10: 6, 12: 42} (0 elsewhere) — the script HARD-FAILS on any
mismatch before entering unexplored bands.

Output: data/a19/m24_census_classes.jsonl (canonical reps).
Usage:  uv run scripts/a19_m_census2.py [--bands 2,...,22] [--max-classes N]
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

OUT = LAB / "data" / "a19"
JSONL = OUT / "m24_census_classes.jsonl"

BASE = {"frame": (30, 3), "A": "x^9 + y + y^2", "B": "1 + x^25 + x^26"}
EXPECTED = {2: 0, 4: 0, 6: 1, 8: 0, 10: 6, 12: 42}   # a19_m_census ground truth


def build_code():
    G = AbelianGroup(BASE["frame"])
    A = Poly.from_string(BASE["A"], G)
    B = Poly.from_string(BASE["B"], G)
    return G, bb_check_matrices(A, B)


def translations(ell, m, n2):
    """Permutation arrays perm[t][i] = index of qubit i under translation t."""
    N = ell * m
    perms = []
    for ta in range(ell):
        for tb in range(m):
            if ta == 0 and tb == 0:
                continue
            p = np.empty(n2, dtype=np.int64)
            for blk in (0, 1):
                for a in range(ell):
                    for b in range(m):
                        src = blk * N + a * m + b
                        dst = blk * N + ((a + ta) % ell) * m + (b + tb) % m
                        p[dst] = src
            perms.append(p)
    return perms


def canonical(b, ell, m):
    """Lex-min translate of b (as uint8 vector), for reblocking stored rows."""
    N = ell * m
    n2 = 2 * N
    best = None
    for ta in range(ell):
        for tb in range(m):
            c = np.empty(n2, dtype=np.uint8)
            for blk in (0, 1):
                blkv = b[blk * N:(blk + 1) * N].reshape(ell, m)
                c[blk * N:(blk + 1) * N] = np.roll(
                    np.roll(blkv, ta, axis=0), tb, axis=1).reshape(-1)
            key = c.tobytes()
            if best is None or key < best[0]:
                best = (key, c)
    return best[1]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bands", default="2,4,6,8,10,12,14,16,18,20,22")
    ap.add_argument("--max-classes", type=int, default=20000)
    args = ap.parse_args()
    bands = [int(w) for w in args.bands.split(",")]

    G, checks = build_code()
    ell, m = BASE["frame"]
    n2 = checks.num_qubits
    HX = checks.H_X % 2
    N = ell * m
    perms = translations(ell, m, n2)

    OUT.mkdir(parents=True, exist_ok=True)
    classes, done_bands = [], set()
    if JSONL.exists():
        for line in JSONL.read_text().splitlines():
            row = json.loads(line)
            if "band_complete" in row:
                done_bands.add(row["band_complete"])
            else:
                classes.append(row)
        print(f"resumed: {len(classes)} classes, bands done "
              f"{sorted(done_bands)}", flush=True)

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
        for j in range(n2):
            idx = np.flatnonzero(HX[:, j])
            solver.add_xor_clause([b_vars[j]] + [f_vars[i] for i in idx],
                                  False)
        solver.add_clause(b_vars)
        card = CardEnc.atmost(lits=b_vars, bound=W, vpool=pool,
                              encoding=EncType.seqcounter)
        for cl in card.clauses:
            solver.add_clause(cl)
        # lex-leader chains: b <=_lex perm(b) for each nonidentity translation
        for p in perms:
            e_prev = None
            for i in range(n2):
                bi, ci = b_vars[i], b_vars[int(p[i])]
                if e_prev is None:
                    solver.add_clause([-bi, ci])          # b_0 <= c_0
                    if i < n2 - 1:
                        e = pool.id()                     # e <-> (b_0 = c_0)
                        solver.add_clause([-e, -bi, ci])
                        solver.add_clause([-e, bi, -ci])
                        solver.add_clause([e, bi, ci])
                        solver.add_clause([e, -bi, -ci])
                        e_prev = e
                else:
                    solver.add_clause([-e_prev, -bi, ci])
                    if i < n2 - 1:
                        e = pool.id()   # e <-> e_prev & (b_i = c_i)
                        solver.add_clause([-e, e_prev])
                        solver.add_clause([-e, -bi, ci])
                        solver.add_clause([-e, bi, -ci])
                        solver.add_clause([e, -e_prev, -bi, -ci])
                        solver.add_clause([e, -e_prev, bi, ci])
                        e_prev = e
        for row in classes:                # block canonical reps, 1 clause each
            b = np.zeros(n2, dtype=np.uint8)
            b[row["b_support"]] = 1
            c = canonical(b, ell, m)
            solver.add_clause([-b_vars[j] if c[j] else b_vars[j]
                               for j in range(n2)])

        found = 0
        while True:
            if len(classes) >= args.max_classes:
                print(f"CAP HIT at {len(classes)} — census INCOMPLETE",
                      flush=True)
                log.close()
                return
            t0 = time.time()
            sat, model = solver.solve()
            dt = time.time() - t0
            if not sat:
                print(f"band W={W}: complete — {found} classes "
                      f"({time.time() - t_band:.0f}s)", flush=True)
                if W in EXPECTED and found != EXPECTED[W]:
                    print(f"VALIDATION FAIL: band {W} found {found}, "
                          f"expected {EXPECTED[W]} — ABORTING", flush=True)
                    log.close()
                    sys.exit(1)
                log.write(json.dumps({"band_complete": W}) + "\n")
                log.flush()
                break
            b = np.array([1 if model[v] else 0 for v in b_vars],
                         dtype=np.uint8)
            f = np.array([1 if model[v] else 0 for v in f_vars],
                         dtype=np.uint8)
            assert (((HX.T @ f) % 2) == b).all()
            w = int(b.sum())
            assert w == W, f"weight {w} != band {W}"
            row = {"w": w, "b_support": np.flatnonzero(b).tolist(),
                   "f_support": np.flatnonzero(f).tolist(), "band": W,
                   "secs": round(dt, 2)}
            classes.append(row)
            log.write(json.dumps(row) + "\n")
            log.flush()
            found += 1
            solver.add_clause([-b_vars[j] if b[j] else b_vars[j]
                               for j in range(n2)])
            if found % 50 == 0:
                print(f"  W={W}: {found} classes ({len(classes)} total, "
                      f"last {dt:.1f}s)", flush=True)

    hist = {}
    for row in classes:
        hist[row["w"]] = hist.get(row["w"], 0) + 1
    print(f"\nCENSUS COMPLETE in {time.time() - total_t0:.0f}s: "
          f"{len(classes)} classes, |b| <= {max(bands)}", flush=True)
    print("weight histogram:", json.dumps(hist, sort_keys=True), flush=True)
    log.close()


if __name__ == "__main__":
    main()
