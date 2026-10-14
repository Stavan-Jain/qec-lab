"""Small convention checks for the exact separable BB comparison family.

The infinite distance proof is in notes/IBM_EXACT_DISTANCE_PROGRAM.md.
This diagnostic checks classical seeds, ranks, and paired witnesses only.
It does not search quantum minimum distances.
"""

from __future__ import annotations

import argparse
import itertools
import json
import resource
import signal
import sys
import time
from pathlib import Path

import numpy as np

LAB = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(LAB / "src"))

from bb_lab.tower import TowerCode, kernel_basis, validate_banked


def guard():
    usage = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    divisor = 1024**3 if sys.platform == "darwin" else 1024**2
    if usage / divisor > 1:
        raise RuntimeError("RSS exceeded 1 GiB")


def timeout_handler(signum, frame):
    raise TimeoutError("Diagnostic exceeded its 60-second budget")


def seed_check(length, offsets, dimension, distance):
    matrix = np.zeros((length, length), dtype=np.uint8)
    for row in range(length):
        for offset in offsets:
            matrix[row, (row - offset) % length] ^= 1
    bases = []
    histograms = []
    for operator in (matrix, matrix.T):
        basis = kernel_basis(operator)
        assert len(basis) == dimension
        histogram = {}
        words = []
        for coefficients in itertools.product((0, 1), repeat=dimension):
            if not any(coefficients):
                continue
            word = np.zeros(length, dtype=np.uint8)
            for coefficient, generator in zip(coefficients, basis):
                if coefficient:
                    word ^= generator
            assert not ((operator @ word) % 2).any()
            weight = int(word.sum())
            assert weight == distance
            histogram[weight] = histogram.get(weight, 0) + 1
            words.append(word)
        bases.append(next(word for word in words if word[0]))
        histograms.append(histogram)
    guard()
    return bases, histograms


def member_check(parameter):
    length_x, length_y = 6 * parameter, 7 * parameter
    distance = 4 * parameter
    words_x, _ = seed_check(length_x, (0, 1, 2), 2, distance)
    words_y, _ = seed_check(length_y, (0, 1, 3), 3, distance)
    code = TowerCode(
        f"separable({parameter})",
        (length_x, length_y),
        frozenset(((0, 0), (1, 0), (2, 0))),
        frozenset(((0, 0), (0, 1), (0, 3))),
    )
    assert code.n == 84 * parameter**2
    assert code.k == 12
    assert np.all(code.HX.sum(axis=1) == 6)
    assert np.all(code.HZ.sum(axis=1) == 6)
    witness_x = np.zeros(code.n, dtype=np.uint8)
    witness_z = np.zeros(code.n, dtype=np.uint8)
    for coordinate in range(length_y):
        witness_x[code.G.index((0, coordinate))] = words_y[1][coordinate]
    for coordinate in range(length_x):
        witness_z[code.G.index((coordinate, 0))] = words_x[0][coordinate]
    assert code.is_cycle(witness_x)
    assert not code.is_stab(witness_x)
    assert not ((code.HX @ witness_z) % 2).any()
    assert int(witness_x @ witness_z) % 2 == 1
    assert int(witness_x.sum()) == int(witness_z.sum()) == distance
    guard()
    return {
        "parameter": parameter,
        "n": code.n,
        "k": code.k,
        "paired_witness_weight": distance,
        "pairing": 1,
        "check_weights": [6, 6],
        "x_support": np.flatnonzero(witness_x).tolist(),
        "z_support": np.flatnonzero(witness_z).tolist(),
        "tier": "rank and witness validation; infinite distance is proved in the note",
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    signal.signal(signal.SIGALRM, timeout_handler)
    signal.alarm(60)
    started = time.monotonic()
    validate_banked(LAB / "data")
    print("validate_banked: PASS", flush=True)
    seeds = []
    for repeats in (1, 2, 3, 4, 8):
        for period, offsets, dimension, ones in (
            (3, (0, 1, 2), 2, 2),
            (7, (0, 1, 3), 3, 4),
        ):
            _, histogram = seed_check(period * repeats, offsets, dimension, ones * repeats)
            seeds.append({"period": period, "repeats": repeats, "spectra": histogram})
    members = []
    for parameter in (1, 2, 3):
        members.append(member_check(parameter))
        print(f"parameter={parameter}: ranks and paired witnesses PASS", flush=True)
    result = {
        "family": "A=1+x+x^2, B=1+y+y^3 on Z_(6t) x Z_(7t)",
        "seeds": seeds,
        "members": members,
        "wall_seconds": round(time.monotonic() - started, 3),
        "distance_claim_source": "analytic proof, not these finite checks",
    }
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))
    signal.alarm(0)


if __name__ == "__main__":
    main()
