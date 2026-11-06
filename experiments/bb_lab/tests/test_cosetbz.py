"""Check low-radius coset census completeness against exhaustive enumeration."""

from __future__ import annotations

import time

import numpy as np
import pytest

from bb_lab import cosetbz


@pytest.fixture
def census_kernel(tmp_path, monkeypatch):
    monkeypatch.setattr(cosetbz, "BUILD", tmp_path / "kernel")
    return cosetbz.build_kernel()


def _words(generator, offset):
    words = {tuple(int(x) for x in offset)}
    for row in generator:
        words |= {
            tuple(int(x) for x in np.bitwise_xor(word, row))
            for word in list(words)
        }
    return words


def _sweep(kernel, tmp_path, generator, offsets, radius, tag):
    first, g1, second, g2, _ = cosetbz.disjoint_info_sets(generator)
    actual = [set() for _ in offsets]

    def on_hit(index, word):
        actual[index].add(tuple(int(x) for x in word))

    result = cosetbz.coset_sweep(
        kernel, tag, (first, g1, second, g2),
        lambda window, systematic: [
            cosetbz.coset_base(systematic, window, offset)
            for offset in offsets
        ],
        radius, time.monotonic() + 30, threads=2,
        workdir=tmp_path, on_hit=on_hit,
    )
    expected = [
        {word for word in _words(generator, offset) if sum(word) <= radius}
        for offset in offsets
    ]
    assert actual == expected
    return result


def test_low_radius_cosets_match_exhaustive_enumeration(census_kernel, tmp_path):
    rng = np.random.default_rng(640718)
    for dimension in range(1, 5):
        identity = np.eye(dimension, dtype=np.uint8)
        invertible = identity.copy()
        for _ in range(8):
            i, j = rng.integers(0, dimension, size=2)
            if i != j:
                invertible[i] ^= invertible[j]
        extra = rng.integers(0, 2, size=(dimension, 2), dtype=np.uint8)
        generator = np.concatenate([identity, invertible, extra], axis=1)
        n = generator.shape[1]
        offsets = [np.zeros(n, dtype=np.uint8), np.ones(n, dtype=np.uint8)]
        offsets.extend(rng.integers(0, 2, size=(2, n), dtype=np.uint8))
        for radius in range(5):
            result = _sweep(
                census_kernel, tmp_path, generator, offsets, radius,
                f"dimension{dimension}_radius{radius}",
            )
            if radius <= 1:
                assert all(window["nodes"] == 0 for window in result["windows"])


def test_zero_code_uses_empty_information_sets(census_kernel, tmp_path):
    generator = np.zeros((2, 5), dtype=np.uint8)
    first, g1, second, g2, dimension = cosetbz.disjoint_info_sets(generator)
    assert first == second == []
    assert dimension == 0 and g1.shape == g2.shape == (0, 5)
    offsets = [
        np.array([0, 0, 0, 0, 0], dtype=np.uint8),
        np.array([1, 0, 0, 0, 0], dtype=np.uint8),
        np.array([1, 1, 1, 1, 1], dtype=np.uint8),
    ]
    for radius in range(5):
        result = _sweep(
            census_kernel, tmp_path, generator, offsets, radius,
            f"zero_code_radius{radius}",
        )
        assert all(window["nodes"] == 0 for window in result["windows"])


def test_negative_window_radius_is_rejected(tmp_path):
    with pytest.raises(ValueError, match="radius must be nonnegative"):
        cosetbz.run_window(
            tmp_path / "unused", "negative", np.zeros((0, 1), dtype=np.uint8),
            [], -1, 0, time.monotonic() + 30, workdir=tmp_path,
        )
