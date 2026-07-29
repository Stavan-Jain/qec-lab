"""Shard-decomposition distance solver: correctness against the
monolithic SAT oracle and the published Bravyi values, plus the
soundness invariants of the class-action machinery.
"""

from __future__ import annotations

import numpy as np
import pytest

from bb_lab.checks import bb_check_matrices
from bb_lab.group import ZmZn
from bb_lab.poly import Poly
from bb_lab.sat_distance import find_logical_z, x_distance
from bb_lab.shard_distance import (
    _bits_to_int,
    _int_to_bits,
    compute_class_action,
    shard_distance,
)


def _build(ell, m, a_str, b_str):
    G = ZmZn(ell, m)
    A = Poly.from_string(a_str, G)
    B = Poly.from_string(b_str, G)
    return bb_check_matrices(A, B)


def _verify_witness(checks, v: np.ndarray) -> None:
    assert not ((checks.H_Z @ v) % 2).any()
    L_Z = find_logical_z(checks)
    assert ((L_Z @ v) % 2).any()


# --- oracle cross-checks on small codes ------------------------------------


def test_toric_z3z3_matches_monolith():
    """[[18,2,3]] toric code as a BB instance — the ρ-trivial edge case
    (translations act trivially on homology, so every orbit is a
    singleton and anchoring carries the whole reduction), and the
    mixed-parity edge case for the parity-tightened bounds (d odd)."""
    checks = _build(3, 3, "1 + x", "1 + y")
    mono = x_distance(checks)
    res = shard_distance(checks)
    assert res.distance == mono.distance == 3
    _verify_witness(checks, res.witness)
    assert int(res.witness.sum()) == res.distance


@pytest.mark.parametrize("orbits_per_shard", [None, 1, 8, 155])
@pytest.mark.parametrize("encoding", ["pin", "coset"])
def test_bb72_granularity_encoding_matrix(orbits_per_shard, encoding):
    """Every point of the granularity dial and both encodings agree
    with the monolith oracle."""
    checks = _build(6, 6, "x^3 + y + y^2", "y^3 + x + x^2")
    res = shard_distance(
        checks, orbits_per_shard=orbits_per_shard, encoding=encoding
    )
    assert res.distance == 6
    _verify_witness(checks, res.witness)


@pytest.mark.parametrize("ell,m,a,b,d_pub", [
    (15, 3, "x^9 + y + y^2", "1 + x^2 + x^7", 10),   # [[90,8,10]]
    (9, 6, "x^3 + y + y^2", "y^3 + x + x^2", 10),    # [[108,8,10]]
])
def test_published_distances(ell, m, a, b, d_pub):
    checks = _build(ell, m, a, b)
    res = shard_distance(checks)
    assert res.distance == d_pub
    _verify_witness(checks, res.witness)


@pytest.mark.slow
def test_gross_144():
    """[[144,12,12]] — the IBM gross code, coarse default config."""
    checks = _build(12, 6, "x^3 + y + y^2", "y^3 + x + x^2")
    res = shard_distance(checks, jobs=2)
    assert res.distance == 12
    _verify_witness(checks, res.witness)


# --- class-action invariants ------------------------------------------------


def test_class_action_invariants_bb72():
    checks = _build(6, 6, "x^3 + y + y^2", "y^3 + x + x^2")
    action = compute_class_action(checks)
    k = action.k
    assert k == 12
    # orbits partition the nonzero classes
    assert sum(action.orbit_sizes) == 2**k - 1
    # μ-normalization: V rows pair to identity
    assert np.array_equal(
        (action.V @ action.L_Z.T) % 2, np.eye(k, dtype=np.uint8)
    )
    # every V row is a valid nontrivial logical
    for row in action.V:
        _verify_witness(checks, row)
    # transversal sizes match |G| / |Stab|
    N = checks.group.cardinality
    for rep in action.orbit_reps:
        assert len(action.transversal[rep]) == N // action.stab_sizes[rep]


def test_transport_soundness_bb72():
    """End-to-end transport check, independent of ρ's construction:
    translate actual logical vectors and verify (a) they stay
    same-weight codewords, (b) their directly-recomputed class lands in
    the same orbit as the original — i.e. the rep shard really covers
    every translate."""
    from bb_lab.shard_distance import _translation_perm

    checks = _build(6, 6, "x^3 + y + y^2", "y^3 + x + x^2")
    action = compute_class_action(checks)
    G = checks.group
    for v in action.V[:4]:
        c = _bits_to_int((action.L_Z @ v) % 2)
        for t in list(G):
            perm = _translation_perm(G, t)
            w = np.zeros_like(v)
            w[perm] = v
            assert not ((checks.H_Z @ w) % 2).any()
            assert w.sum() == v.sum()
            c_t = _bits_to_int((action.L_Z @ w) % 2)
            assert c_t != 0, "translate of a logical became a stabilizer"
            assert (
                action.orbit_rep_of[c_t] == action.orbit_rep_of[c]
            ), "translated class escaped its orbit"


def test_bits_roundtrip():
    for c in [1, 5, 4095, 2**11]:
        assert _bits_to_int(_int_to_bits(c, 12)) == c
