"""Regression tests for the A28 general fibering engine.

Fast ground-truth anchors (seconds):
  * pair72 base [[36,4,4]]: safe floor certify@8 / refute@10 with
    weight-8 witnesses (A14 anchors_screens.jsonl exact minima).
  * f2a6 [[150,8,8]]: seam-floor-16 sweep counts must match A23 §9.2
    bit-for-bit (300 consistent / 6135 inconsistent, kdims {0:120,
    4:180}, 3000 reps); boundary floor 6 certified with 75 weight-6
    violations at 8 (A22's single free orbit).

The 113-class census (65 s) lives in scripts/, not here.
"""

from __future__ import annotations

import numpy as np
import pytest

from bb_lab.fibering import (
    FiberFrame,
    FiberSweep,
    best_frames,
    enumerate_fiber_generators,
    kernel_orbit_reps,
    safe_floor_certify,
    seam_offsets,
)
from bb_lab.group import AbelianGroup
from bb_lab.poly import Poly


@pytest.fixture(scope="module")
def pair72():
    G = AbelianGroup((3, 6))
    A = Poly.from_support([(2, 0), (0, 1), (0, 3)], G)
    B = Poly.from_support([(0, 0), (1, 0), (0, 2)], G)
    return A, B


@pytest.fixture(scope="module")
def f2a6():
    G = AbelianGroup((5, 15))
    A = Poly.from_support([(0, 0), (0, 1), (1, 0)], G)
    B = Poly.from_support([(1, 6), (1, 10), (2, 12)], G)
    return A, B


def test_weight_formula_all_frames(pair72):
    A, B = pair72
    rng = np.random.default_rng(0)
    for z in enumerate_fiber_generators(A.group):
        fr = FiberFrame(A, B, z)
        for _ in range(25):
            v = rng.integers(0, 2, fr.n).astype(np.uint8)
            fr.chain_weight_check(v)


def test_pair72_floor_exactly_8(pair72):
    A, B = pair72
    rep8 = safe_floor_certify(A, B, axis=0, target=8)
    assert rep8.certified and not rep8.refuted
    rep10 = safe_floor_certify(A, B, axis=0, target=10)
    assert rep10.refuted
    wts = {v["weight"] for r in rep10.per_class for v in r.violations}
    assert wts == {8}


def test_pair72_seam_raw_weight(pair72):
    A, B = pair72
    offs = seam_offsets(A, B, axis=0)
    assert len(offs) == 1  # single G-orbit
    _, su, sv = offs[0]
    assert int(su.sum() + sv.sum()) == 12  # A14 s0_raw_weights


def test_f2a6_frame_structure(f2a6):
    A, B = f2a6
    fr = FiberFrame(A, B, (0, 3))
    assert (fr.q, fr.S, fr.r_delta, len(fr.K0)) == (5, 15, 56, 19)
    # the ε-monomial link at the site of x̄
    site_x = int(fr.site_of[list(A.group).index((1, 0))])
    assert fr.link_shift == site_x
    reps = kernel_orbit_reps(A, B)
    assert [int(z.sum()) for z in reps] == [40]  # e₀ orbit


def test_f2a6_seam_floor_16_counts(f2a6):
    A, B = f2a6
    rep = safe_floor_certify(A, B, axis=1, target=16, z=(0, 3))
    assert rep.certified
    (res,) = rep.per_class
    # 300 consistent leaves with kdims {0:120, 4:180} → 3000 reps is the
    # A23 §9.2 ground truth; inconsistent *subtrees* are pruned wholesale
    # by the DFS, so n_incons counts prune events, not the 6135 leaves.
    assert res.n_cons == 300 and res.n_incons > 0
    assert res.kdims == {0: 120, 4: 180}
    assert res.n_reps == 3000


def test_f2a6_boundary_floor(f2a6):
    A, B = f2a6
    sw = FiberSweep(FiberFrame(A, B, (0, 3)))
    z0 = np.zeros(75, dtype=np.uint8)
    r6 = sw.floor_sweep((z0, z0), 6)
    assert r6.certified
    r8 = sw.floor_sweep((z0, z0), 8)
    assert not r8.certified
    assert len(r8.violations) == 75  # the ∂₂(monomial) free orbit
    assert all(v["exact"] == 6 for v in r8.violations)


def test_feasibility_scorer_shapes(f2a6):
    A, B = f2a6
    scores = best_frames(A, B, 16)
    by_z = {s.z: s for s in scores}
    s = by_z[(0, 3)]
    assert s.mode == "paired" and s.feasible and s.n_masks == 6435
