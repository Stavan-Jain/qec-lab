"""Depth-7 (8-round) BB-style syndrome-measurement circuit for a two-block
code given by permutation parts (Table 5 of arXiv:2308.07915), built as a
Stim memory experiment, plus circuit-level distance search.

Qubit layout: X-anc c ∈ [0,m), L data [m,2m), R data [2m,3m), Z-anc [3m,4m).
Round schedule (unitary part, Eq. (9) of the paper):
  1: CNOT_{A1}(R→Z)            2: CNOT_{A2}(X→L), CNOT_{A3}(R→Z)
  3: CNOT_{B2}(X→R), CNOT_{B1}(L→Z)    4: CNOT_{B1}(X→R), CNOT_{B2}(L→Z)
  5: CNOT_{B3}(X→R), CNOT_{B3}(L→Z)    6: CNOT_{A1}(X→L), CNOT_{A2}(R→Z)
  7: CNOT_{A3}(X→L), MeasZ(Z)          8: MeasX(X), InitZ(Z) [InitX(X) at round 1]
X→L via A_p: control X(c), target L(A_p[c]).  X→R via B_q: target R(B_q[c]).
R→Z via A_p: control R(A_p^{-1}[z]), target Z(z).  L→Z via B_q: control L(B_q^{-1}[z]).
"""
from __future__ import annotations

import numpy as np
import stim

from .structure import Perms


def _inv(p):
    q = np.empty_like(p); q[p] = np.arange(len(p)); return q


class Schedule:
    """Order of the 7 CNOT layers, as (X-side perm, Z-side perm) per round.
    Default = the paper's Eq. (9)."""

    def __init__(self, x_order=("A2", "B2", "B1", "B3", "A1", "A3"), z_order=("A1", "A3", "B1", "B2", "B3", "A2")):
        # rounds 2..7 for X (6 gates), rounds 1..6 for Z (6 gates)
        self.x_order = x_order; self.z_order = z_order

    def rounds(self):
        # round r (1..7): X gate index r-2 (if 0<=r-2<6), Z gate index r-1 (if r-1<6)
        out = []
        for r in range(1, 8):
            xg = self.x_order[r - 2] if 1 <= r - 1 <= 6 else None
            zg = self.z_order[r - 1] if r - 1 < 6 else None
            out.append((xg, zg))
        return out


def build_memory_circuit(P: Perms, H_X: np.ndarray, H_Z: np.ndarray, logicals: np.ndarray,
                         basis: str, cycles: int, p: float, schedule: Schedule | None = None) -> stim.Circuit:
    """Memory experiment: data prepared in |0⟩^n (basis='Z', detects X errors via Z-checks,
    observables = Z-logicals) or |+⟩^n (basis='X', detects Z errors via X-checks, observables =
    X-logicals). `logicals` rows are length-n binary vectors over data qubits (L then R).
    Noise: standard circuit-level depolarizing model with rate p (paper's model)."""
    sched = schedule or Schedule()
    m = P.m; n = 2 * m
    Xa = list(range(0, m)); L = list(range(m, 2 * m)); R = list(range(2 * m, 3 * m)); Za = list(range(3 * m, 4 * m))
    data = L + R
    perm = {f"A{i+1}": P.A[i] for i in range(3)}; perm.update({f"B{i+1}": P.B[i] for i in range(3)})
    inv = {k: _inv(v) for k, v in perm.items()}

    def cnot_pairs(gate, side):
        if gate is None:
            return []
        pairs = []
        if side == "X":
            if gate[0] == "A":
                for c in range(m): pairs += [Xa[c], L[perm[gate][c]]]
            else:
                for c in range(m): pairs += [Xa[c], R[perm[gate][c]]]
        else:
            if gate[0] == "A":
                for z in range(m): pairs += [R[inv[gate][z]], Za[z]]
            else:
                for z in range(m): pairs += [L[inv[gate][z]], Za[z]]
        return pairs

    def noisy(c: stim.Circuit, targets, gate):
        c.append(gate, targets)

    circ = stim.Circuit()
    # data prep
    if basis == "Z":
        circ.append("R", data)
    else:
        circ.append("RX", data)
    circ.append("R", Za)   # Z ancillas start in |0>
    if p > 0:
        circ.append("X_ERROR" if basis == "Z" else "Z_ERROR", data, p)
        circ.append("X_ERROR", Za, p)
    circ.append("TICK")

    meas_index = 0
    rec = {}  # (cycle, 'Z'|'X', ancilla) -> measurement index
    for t in range(cycles):
        for r, (xg, zg) in enumerate(sched.rounds(), start=1):
            active = set()
            if r == 1:
                circ.append("RX", Xa)
                if p > 0: circ.append("Z_ERROR", Xa, p)
                active |= set(Xa)
            pairs = cnot_pairs(xg, "X") + cnot_pairs(zg, "Z")
            if pairs:
                circ.append("CX", pairs)
                if p > 0: circ.append("DEPOLARIZE2", pairs, p)
                active |= set(pairs)
            if r == 7:
                if p > 0: circ.append("X_ERROR", Za, p)
                circ.append("M", Za)
                for z in range(m):
                    rec[(t, "Z", z)] = meas_index; meas_index += 1
                active |= set(Za)
            idle = [q for q in range(4 * m) if q not in active]
            if p > 0 and idle: circ.append("DEPOLARIZE1", idle, p)
            circ.append("TICK")
        # round 8: MeasX(X), InitZ(Z)
        if p > 0: circ.append("Z_ERROR", Xa, p)
        circ.append("MX", Xa)
        for c in range(m):
            rec[(t, "X", c)] = meas_index; meas_index += 1
        circ.append("R", Za)
        if p > 0:
            circ.append("X_ERROR", Za, p)
            circ.append("DEPOLARIZE1", data, p)
        circ.append("TICK")
        # detectors
        if basis == "Z":
            for z in range(m):
                if t == 0:
                    circ.append("DETECTOR", [stim.target_rec(rec[(0, "Z", z)] - meas_index)])
                else:
                    circ.append("DETECTOR", [stim.target_rec(rec[(t, "Z", z)] - meas_index),
                                             stim.target_rec(rec[(t - 1, "Z", z)] - meas_index)])
            if t >= 1:
                for c in range(m):
                    circ.append("DETECTOR", [stim.target_rec(rec[(t, "X", c)] - meas_index),
                                             stim.target_rec(rec[(t - 1, "X", c)] - meas_index)])
        else:
            for c in range(m):
                if t == 0:
                    circ.append("DETECTOR", [stim.target_rec(rec[(0, "X", c)] - meas_index)])
                else:
                    circ.append("DETECTOR", [stim.target_rec(rec[(t, "X", c)] - meas_index),
                                             stim.target_rec(rec[(t - 1, "X", c)] - meas_index)])
            if t >= 1:
                for z in range(m):
                    circ.append("DETECTOR", [stim.target_rec(rec[(t, "Z", z)] - meas_index),
                                             stim.target_rec(rec[(t - 1, "Z", z)] - meas_index)])
    # final data measurement
    if p > 0:
        circ.append("X_ERROR" if basis == "Z" else "Z_ERROR", data, p)
    circ.append("M" if basis == "Z" else "MX", data)
    data_rec = {q: meas_index + i for i, q in enumerate(data)}
    meas_index += len(data)
    T = cycles - 1
    if basis == "Z":
        for z in range(m):
            targets = [stim.target_rec(rec[(T, "Z", z)] - meas_index)]
            for q in np.flatnonzero(H_Z[z]):
                targets.append(stim.target_rec(data_rec[data[q]] - meas_index))
            circ.append("DETECTOR", targets)
    else:
        for c in range(m):
            targets = [stim.target_rec(rec[(T, "X", c)] - meas_index)]
            for q in np.flatnonzero(H_X[c]):
                targets.append(stim.target_rec(data_rec[data[q]] - meas_index))
            circ.append("DETECTOR", targets)
    for i, lg in enumerate(logicals):
        targets = [stim.target_rec(data_rec[data[q]] - meas_index) for q in np.flatnonzero(lg)]
        circ.append("OBSERVABLE_INCLUDE", targets, i)
    return circ


def circuit_distance_upper_bound(circ: stim.Circuit, max_set_size: int = 5, max_degree: int = 9999):
    """Stim's exhaustive-ish search for a minimum-weight undetectable logical error."""
    errs = circ.search_for_undetectable_logical_errors(
        dont_explore_detection_event_sets_with_size_above=max_set_size,
        dont_explore_edges_with_degree_above=max_degree,
        dont_explore_edges_increasing_symptom_degree=False,
        canonicalize_circuit_errors=True,
    )
    return len(errs), errs
