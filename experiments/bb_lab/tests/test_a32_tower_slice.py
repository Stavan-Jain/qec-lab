"""Regression tests for the `a32_tower_slice` linear-algebra helpers.

`_preimage(M, Wb, Wp)` used to raise ValueError when span(Wb) was the
FULL codomain: the annihilator basis is empty, ``np.array([])`` of the
empty kernel list is 1-D, and the ``F @ M`` matmul died.  Mathematically
the preimage of the full space is the whole domain.  Never hit on the
A32/A33/A35 production towers (proper W-subspaces); first hit 2026-08-17
by the descent-theory-test fold-in over the order-144 sweep (towers
where p_bot*(S) spans H1(bottom)) -- see the descent-theory-test
AMENDMENTS.md #1 for the field analysis.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

LAB = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(LAB / "src"))
sys.path.insert(0, str(LAB / "scripts"))

from a32_tower_slice import _preimage, i2v, rref_ints, v2i  # noqa: E402


def _xor_span(basis: list[int]) -> set[int]:
    pts = {0}
    for b in basis:
        pts |= {p ^ b for p in pts}
    return pts


def test_preimage_full_codomain_returns_domain_basis():
    """span(Wb) = full codomain => the whole domain (used to raise)."""
    M = np.array([[1, 0, 1, 1],
                  [0, 1, 1, 0]], dtype=np.uint8)
    Wb, Wp = rref_ints([0b01, 0b11])  # spans all of F_2^2
    pre = _preimage(M, Wb, Wp)
    assert len(pre) == M.shape[1]
    basis, _ = rref_ints(list(pre))
    assert len(basis) == M.shape[1], "not a basis of the full domain"


def test_preimage_matches_enumeration():
    """Agreement with brute force on random small instances, across
    zero, proper and full W-spans (the last is the degenerate case)."""
    rng = np.random.default_rng(20260817)
    saw_full = False
    for _ in range(40):
        n = int(rng.integers(1, 5))        # codomain dim
        dom = int(rng.integers(1, 6))      # domain dim
        M = rng.integers(0, 2, size=(n, dom)).astype(np.uint8)
        nw = int(rng.integers(0, n + 1))
        raw = [v2i(r)
               for r in rng.integers(0, 2, size=(nw, n)).astype(np.uint8)]
        Wb, Wp = rref_ints([x for x in raw if x])
        saw_full = saw_full or len(Wb) == n
        pre = _preimage(M, Wb, Wp)
        span_W = _xor_span(Wb)
        want = {s for s in range(1 << dom)
                if v2i((M @ i2v(s, dom)) % 2) in span_W}
        assert _xor_span(pre) == want, "preimage span != enumeration"
        assert len(rref_ints(list(pre))[0]) == len(pre), "dependent basis"
    assert saw_full, "battery never hit the full-span degenerate case"
