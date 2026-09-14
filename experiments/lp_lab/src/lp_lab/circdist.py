"""Circuit-level distance via SAT on the detector error model.

Stim's DEM lists error mechanisms e with detector sets D(e) and observable sets L(e).
An undetectable logical fault set S has ⊕_{e∈S} D(e) = ∅ and ⊕_{e∈S} L(e) ≠ ∅; d_circ is
the minimum |S| (each mechanism = one faulty operation, matching the arXiv:2308.07915 count
since Stim's DEPOLARIZE2 mechanisms are the 15 single-gate Paulis).  We solve the ladder
w = lower.. with CryptoMiniSat (native XOR per detector), seqcounter cardinality.
"""
from __future__ import annotations

import time
from dataclasses import dataclass

import numpy as np
import stim
from pysat.card import CardEnc, EncType
from pysat.formula import IDPool


def dem_incidence(circ: stim.Circuit):
    dem = circ.detector_error_model(flatten_loops=True)
    mechs = []   # (detectors tuple, observables tuple, probability)
    for inst in dem:
        if inst.type == "error":
            dets = []; obs = []
            for t in inst.targets_copy():
                if t.is_relative_detector_id(): dets.append(t.val)
                elif t.is_logical_observable_id(): obs.append(t.val)
            mechs.append((tuple(sorted(dets)), tuple(sorted(obs)), inst.args_copy()[0]))
    return mechs, circ.num_detectors, circ.num_observables


@dataclass
class CircDistResult:
    distance: int
    mechanisms: list
    unsat_weights: list
    log: list


def circuit_distance_sat(circ: stim.Circuit, lower: int = 1, upper: int = 30, verbose: bool = True,
                         fix_first: bool = False) -> CircDistResult:
    import pycryptosat
    mechs, nD, nO = dem_incidence(circ)
    M = len(mechs)
    det_to_m = [[] for _ in range(nD)]; obs_to_m = [[] for _ in range(nO)]
    for i, (dets, obs, _) in enumerate(mechs):
        for d in dets: det_to_m[d].append(i + 1)
        for o in obs: obs_to_m[o].append(i + 1)
    log = []; unsat = []
    for w in range(lower, upper + 1):
        t0 = time.time()
        s = pycryptosat.Solver(threads=4)
        pool = IDPool(start_from=M + 1)
        for d in range(nD):
            if det_to_m[d]:
                s.add_xor_clause(det_to_m[d], False)
        o_lits = []
        for o in range(nO):
            if not obs_to_m[o]: continue
            a = pool.id(("o", o)); s.add_xor_clause(obs_to_m[o] + [a], False); o_lits.append(a)
        s.add_clause(o_lits)
        card = CardEnc.atmost(lits=list(range(1, M + 1)), bound=w, vpool=pool, encoding=EncType.seqcounter)
        for cl in card.clauses: s.add_clause(cl)
        sat, model = s.solve()
        dt = time.time() - t0
        if sat:
            S = [i for i in range(M) if model[i + 1]]
            # verify
            det_par = np.zeros(nD, dtype=np.uint8); obs_par = np.zeros(nO, dtype=np.uint8)
            for i in S:
                for d in mechs[i][0]: det_par[d] ^= 1
                for o in mechs[i][1]: obs_par[o] ^= 1
            assert not det_par.any() and obs_par.any(), "SAT witness failed verification"
            log.append((w, "SAT", dt))
            if verbose: print(f"  [d_circ ≤ {w}] SAT with {len(S)} mechanisms ({dt:.1f}s)", flush=True)
            return CircDistResult(len(S), [mechs[i] for i in S], unsat, log)
        unsat.append(w); log.append((w, "UNSAT", dt))
        if verbose: print(f"  [d_circ ≤ {w}] UNSAT ({dt:.1f}s)", flush=True)
    raise RuntimeError("no undetectable logical fault set found up to the upper bound")


def circuit_distance_bposd(circ: stim.Circuit, trials: int = 100, seed: int = 0, osd_order: int = 10,
                           max_iter: int = 60, verbose: bool = False) -> dict:
    """Upper bound on d_circ by the arXiv:2308.07915 §6 method: for random η = (random
    combination of detector rows) + (random nonzero combination of observable rows), decode the
    syndrome (0,…,0,1) of the extended matrix [D; η] with BP-OSD; every solution e satisfies
    D e = 0, η e = 1, hence is an undetectable logical fault set; return the minimum weight."""
    import numpy as np
    import scipy.sparse as sp
    from ldpc import BpOsdDecoder
    mechs, nD, nO = dem_incidence(circ)
    M = len(mechs)
    rows, cols = [], []
    orow, ocol = [], []
    for i, (dets, obs, _) in enumerate(mechs):
        for d in dets: rows.append(d); cols.append(i)
        for o in obs: orow.append(o); ocol.append(i)
    D = sp.csr_matrix((np.ones(len(rows), dtype=np.uint8), (rows, cols)), shape=(nD, M))
    L = sp.csr_matrix((np.ones(len(orow), dtype=np.uint8), (orow, ocol)), shape=(nO, M))
    rng = np.random.default_rng(seed)
    best = None; best_e = None
    Dd = D.toarray().astype(np.uint8); Ld = L.toarray().astype(np.uint8)
    for t in range(trials):
        cD = rng.integers(0, 2, nD).astype(np.uint8); cL = rng.integers(0, 2, nO).astype(np.uint8)
        if not cL.any(): cL[rng.integers(nO)] = 1
        eta = ((cD @ Dd) + (cL @ Ld)) & 1
        H = sp.vstack([D, sp.csr_matrix(eta.reshape(1, -1))]).tocsr()
        dec = BpOsdDecoder(H, error_rate=0.01, max_iter=max_iter, bp_method="minimum_sum", ms_scaling_factor=0.9,
                           osd_method="osd_cs", osd_order=osd_order)
        s = np.zeros(nD + 1, dtype=np.uint8); s[-1] = 1
        e = dec.decode(s).astype(np.uint8)
        ok = (not ((Dd @ e) & 1).any()) and ((Ld @ e) & 1).any()
        if ok:
            w = int(e.sum())
            if best is None or w < best:
                best, best_e = w, e
                if verbose: print(f"   trial {t}: new best d_circ ≤ {w}", flush=True)
    return {"dcirc_ub": best, "trials": trials, "mechanisms": M, "detectors": nD,
            "witness_mechanisms": [mechs[i] for i in np.flatnonzero(best_e)] if best_e is not None else None}
