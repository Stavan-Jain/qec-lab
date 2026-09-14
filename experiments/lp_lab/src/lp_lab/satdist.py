"""Exact CSS distance by a CryptoMiniSat ladder with translation-symmetry breaking.

d_X = min{|v| : H_Z v = 0, v ∉ rowspace(H_X)}.  Encoding per weight w:
  * one Boolean per qubit; each H_Z row is a native XOR clause (= 0);
  * for each Z-logical representative L_j: XOR(v ∩ L_j) = a_j (aux literal), clause ⋁ a_j;
  * Σ v ≤ w via pysat CardEnc (sequential counter);
  * symmetry breaking: if a group G acts on the code by qubit permutations (translations),
    any weight-w logical can be translated so that a chosen representative qubit r of
    some orbit carries a 1; we solve one instance per orbit representative with v_r = 1
    (sound: the union over representatives covers all logicals up to symmetry).
Witnesses are re-verified.  Distances are reported with the UNSAT weights certified.
"""
from __future__ import annotations

import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass

import numpy as np
from pysat.card import CardEnc, EncType
from pysat.formula import IDPool


def _build_and_solve(args):
    HZ_rows, LZ_rows, n, w, fixed_one, seed = args
    import pycryptosat
    s = pycryptosat.Solver(threads=1)
    pool = IDPool(start_from=n + 1)
    # parity rows
    for r in HZ_rows:
        lits = [i + 1 for i in r]
        if lits:
            s.add_xor_clause(lits, False)
    # logical detectors
    a_lits = []
    for L in LZ_rows:
        a = pool.id(("a", len(a_lits)))
        s.add_xor_clause([i + 1 for i in L] + [a], False)   # XOR(v∩L) ⊕ a = 0  ⇒ a = parity
        a_lits.append(a)
    s.add_clause(a_lits)
    # cardinality
    card = CardEnc.atmost(lits=list(range(1, n + 1)), bound=w, vpool=pool, encoding=EncType.seqcounter)
    for cl in card.clauses:
        s.add_clause(cl)
    if fixed_one is not None:
        s.add_clause([fixed_one + 1])
    sat, model = s.solve()
    if not sat:
        return (fixed_one, None)
    v = np.zeros(n, dtype=np.uint8)
    for i in range(n):
        if model[i + 1]:
            v[i] = 1
    return (fixed_one, v)


@dataclass
class SymDistanceResult:
    distance: int
    witness: np.ndarray
    unsat_weights: list
    log: list


def x_distance_sym(H_X: np.ndarray, H_Z: np.ndarray, LZ: np.ndarray, reps: list[int] | None,
                   lower: int, upper: int, workers: int = 4, verbose: bool = True) -> SymDistanceResult:
    """Ladder w = lower..upper.  `reps`: one qubit index per symmetry orbit (None → no
    symmetry breaking, single instance).  Rungs are UNSAT only if all rep-instances are UNSAT."""
    n = H_Z.shape[1]
    HZ_rows = [list(map(int, np.flatnonzero(r))) for r in H_Z]
    LZ_rows = [list(map(int, np.flatnonzero(r))) for r in LZ]
    log = []; unsat = []
    for w in range(lower, upper + 1):
        t0 = time.time()
        tasks = [(HZ_rows, LZ_rows, n, w, r, 0) for r in (reps if reps else [None])]
        witness = None
        with ProcessPoolExecutor(max_workers=min(workers, len(tasks))) as ex:
            futs = [ex.submit(_build_and_solve, t) for t in tasks]
            for f in as_completed(futs):
                r, v = f.result()
                if v is not None:
                    witness = v
                    for g in futs: g.cancel()
                    break
        dt = time.time() - t0
        if witness is not None:
            wt = int(witness.sum())
            assert not ((H_Z.astype(np.int64) @ witness) & 1).any(), "witness not in ker H_Z"
            assert ((LZ.astype(np.int64) @ witness) & 1).any(), "witness is trivial"
            log.append((w, "SAT", dt))
            if verbose: print(f"  [w ≤ {w:3d}] SAT   weight {wt}  ({dt:.1f}s)", flush=True)
            return SymDistanceResult(wt, witness, unsat, log)
        unsat.append(w); log.append((w, "UNSAT", dt))
        if verbose: print(f"  [w ≤ {w:3d}] UNSAT ({dt:.1f}s)", flush=True)
    raise RuntimeError("no logical found up to upper bound")
