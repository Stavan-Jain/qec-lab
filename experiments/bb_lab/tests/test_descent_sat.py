"""Tests for the descent-strengthened shard SAT layer.

The Tandem-side orchestration (`maxsat_distance_descent`) needs the
fork binary and is exercised by `scripts/descent_maxsat_ab.py` (its
distances are asserted there); these tests cover the plain-SAT layer,
which is self-contained.
"""

from __future__ import annotations

import numpy as np
import pytest

from bb_lab.checks import bb_check_matrices
from bb_lab.descent_sat import (
    axis_decks,
    base_coset_floors,
    compute_descent,
    descent_shard_distance,
)
from bb_lab.group import ZmZn
from bb_lab.poly import Poly
from bb_lab.sat_distance import x_distance
from bb_lab.shard_distance import compute_class_action


def _bb72():
    G = ZmZn(6, 6)
    A = Poly.from_string("x^3 + y + y^2", G)
    B = Poly.from_string("y^3 + x + x^2", G)
    return A, B, bb_check_matrices(A, B)


def _bb108():
    G = ZmZn(9, 6)
    A = Poly.from_string("x^3 + y + y^2", G)
    B = Poly.from_string("y^3 + x + x^2", G)
    return A, B, bb_check_matrices(A, B)


def test_bb72_both_decks_d6():
    A, B, checks = _bb72()
    decks = axis_decks(checks)
    assert decks == [(3, 0), (0, 3)]
    for sigma in decks:
        r = descent_shard_distance(checks, A, B, sigma=sigma)
        assert r.distance == 6
        # witness is a genuine logical (driver re-verifies; sanity here)
        assert int(r.witness.sum()) == 6


def test_bb108_with_floors_d10():
    A, B, checks = _bb108()
    r = descent_shard_distance(
        checks, A, B, sigma=(0, 3), use_floors=True
    )
    assert r.distance == 10
    assert r.floor_seconds > 0


def test_sector_a_analytic_kills_match_les_prediction():
    """A14 LES: ker p₁ = im τ₁ has dim k/2, so at most 2^(k/2)−1 of the
    2^k−1 classes can be (a)-reachable; with the parity filter the
    reachable orbit-rep count must be well below the total."""
    A, B, checks = _bb72()
    r = descent_shard_distance(checks, A, B, sigma=(0, 3))
    assert r.a_reachable_reps + r.a_killed_reps == r.num_orbits
    assert r.a_killed_reps > r.num_orbits // 2


def test_floors_agree_with_base_distance():
    """min over nonzero-class floors == d_X(base): the floor engine and
    the trusted base SAT must tell the same story."""
    A, B, checks = _bb72()
    action = compute_class_action(checks)
    dd = compute_descent(checks, action, (0, 3), A=A, B=B)
    classes = set(range(1, 2 ** dd.kb))
    floors = base_coset_floors(dd, classes, cap=12)
    d_base = x_distance(dd.base_checks).distance
    assert min(floors.values()) == d_base


def test_non_descending_input_fails_loud():
    """Falsify-first: a code that does NOT descend along σ must raise
    in the verification battery, not return a wrong decomposition.
    Swapping A and B breaks the intertwining (the base is built from
    the swapped images while the cover uses the true polynomials)."""
    A, B, checks = _bb72()
    action = compute_class_action(checks)
    with pytest.raises(AssertionError):
        compute_descent(checks, action, (3, 0), A=B, B=A)


def test_pushforward_witness_expansion():
    """p* of a sector-(a) model must be σ-invariant with doubled
    weight (the dichotomy's arithmetic)."""
    A, B, checks = _bb72()
    action = compute_class_action(checks)
    dd = compute_descent(checks, action, (0, 3), A=A, B=B)
    rng = np.random.default_rng(3)
    vb = rng.integers(0, 2, size=dd.S.shape[0]).astype(np.uint8)
    v = (dd.P_lift @ vb) % 2
    assert int(v.sum()) == 2 * int(vb.sum())
    assert not ((dd.S @ v) % 2).any()  # p₊ p* = 0


def test_fiber_certificate_format(tmp_path):
    """emit_fiber_certificate: parseable, consistent with the solver's
    reader contract, and arithmetically sane (invariant floor even and
    >= 2, moving floors >= 1, pairs partition the qubit ids)."""
    from bb_lab.maxsat_distance import emit_fiber_certificate, write_wcnf

    A, B, checks = _bb72()
    qv, a_lits = write_wcnf(
        checks, tmp_path / "naive.wcnf", mode="naive"
    )
    p = emit_fiber_certificate(
        checks, A, B, (0, 3), qv, a_lits, tmp_path / "c.flb"
    )
    lines = p.read_text().splitlines()
    tag, k, npairs, inv = (
        lines[0].split()[1], int(lines[0].split()[2]),
        int(lines[0].split()[3]), int(lines[0].split()[4]),
    )
    assert (tag, k) == ("fiberlb", 12)
    assert npairs == checks.num_qubits // 2
    assert inv >= 2 and inv % 2 == 0
    avars = [int(x) for x in lines[1].split()[1:]]
    assert avars == [int(x) for x in a_lits]
    tbl = [int(x) for x in lines[2].split()[1:]]
    assert len(tbl) == 1 << k and all(f >= 1 for f in tbl[1:])
    pair_ids = [int(x) for l in lines[3:] for x in l.split()]
    assert sorted(pair_ids) == sorted(int(x) for x in qv)


def test_moving_cost_floor_v2():
    """The v2 cost-floor engine: certified F2 must be >= the v1
    base-coset relaxation everywhere, and on the dangerous fiber
    (lam = 0, nontrivial-pinned) it must reach the code distance on
    bb_72 (d = 6: no lighter nontrivial logical exists, moving or
    not)."""
    from bb_lab.descent_sat import (
        base_coset_floors,
        moving_cost_floor_budgeted,
    )

    A, B, checks = _bb72()
    action = compute_class_action(checks)
    dd = compute_descent(checks, action, (0, 3), A=A, B=B)
    v1 = base_coset_floors(dd, {0, 1}, cap=8)
    f0 = moving_cost_floor_budgeted(
        checks, dd, 0, cap=5, confl_budget=500_000, start=v1[0]
    )
    f1 = moving_cost_floor_budgeted(
        checks, dd, 1, cap=5, confl_budget=500_000, start=v1[1]
    )
    assert f0 >= v1[0] and f1 >= v1[1]
    assert f0 >= 6  # dangerous fiber at the true distance


def test_masked_K_fast_matches_reference():
    """masked_K_value (residual-hash) must agree with the exhaustive
    masked_K_at_least reference on random base cycles (bb_72)."""
    from bb_lab.twist_floors import (
        compute_twist,
        masked_K_at_least,
        masked_K_value,
    )
    from bb_lab.linalg import nullspace_f2

    A, B, checks = _bb72()
    action = compute_class_action(checks)
    dd = compute_descent(checks, action, (0, 3), A=A, B=B)
    td = compute_twist(checks, dd)
    ker = nullspace_f2(td.Hb_Z)
    rng = np.random.default_rng(5)
    for _ in range(15):
        coeff = rng.integers(0, 2, size=ker.shape[0]).astype(np.uint8)
        w = (coeff @ ker) % 2
        if not w.any():
            continue
        fast = masked_K_value(td, w, 3)
        for k in (1, 2, 3):
            assert masked_K_at_least(td, w, k) == (fast >= k), (
                f"disagreement at k={k}: fast={fast}"
            )
