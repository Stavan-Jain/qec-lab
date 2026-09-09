"""Tiny exact checks for A43_S2_codex_audit.md; no solvers or searches."""

from __future__ import annotations

import json
import resource
import signal
import time


ORIGINAL_A = ((0, 0), (0, 1), (3, -1))
ORIGINAL_B = ((0, 0), (1, 0), (-1, -3))
SEED = frozenset(
    [(0, 1, 0), (0, 1, 2), (0, 2, 1), (0, 2, 3), (0, 3, 1), (0, 3, 5)]
    + [(1, 3, 1), (1, 3, 5), (1, 4, 0), (1, 4, 1), (1, 4, 2), (1, 4, 5)]
)
DUAL_TILE = frozenset(
    [(0, 0, 1), (0, 0, 2), (0, 0, 3), (0, 1, 2), (0, 1, 3),
     (0, 2, 2), (0, 3, 0), (0, 3, 2), (0, 4, 1),
     (1, 0, 4), (1, 1, 2), (1, 2, 3)]
)


def toggle(support, cell):
    if cell in support:
        support.remove(cell)
    else:
        support.add(cell)


def coordinate(value, period):
    return value if period is None else value % period


def syndrome(cells, periods, first, second, dual=False):
    output = set()
    supports = (first, second) if dual else (second, first)
    direction = 1 if dual else -1
    for block, column, row in cells:
        for shift_x, shift_y in supports[block]:
            check = (
                coordinate(column + direction * shift_x, periods[0]),
                coordinate(row + direction * shift_y, periods[1]),
            )
            toggle(output, check)
    return output


def boundary(column, row, periods):
    return frozenset(
        (block, coordinate(column - shift_x, periods[0]),
         coordinate(row - shift_y, periods[1]))
        for block, support in enumerate((ORIGINAL_A, ORIGINAL_B))
        for shift_x, shift_y in support
    )


def longest_gap(cells, axis, period):
    occupied = sorted({cell[axis + 1] % period for cell in cells})
    return max((occupied[(index + 1) % len(occupied)] - current - 1) % period
               for index, current in enumerate(occupied))


def matrix_rows(periods, first, second, dual=False):
    column_period, row_period = periods
    group_size = column_period * row_period
    rows = [0] * group_size
    for block in (0, 1):
        for column in range(column_period):
            for row in range(row_period):
                qubit = block * group_size + column * row_period + row
                for check_x, check_y in syndrome(
                    {(block, column, row)}, periods, first, second, dual
                ):
                    rows[check_x * row_period + check_y] ^= 1 << qubit
    return rows


def rank(rows):
    pivots = {}
    for row in rows:
        while row:
            pivot = row.bit_length() - 1
            if pivot in pivots:
                row ^= pivots[pivot]
            else:
                pivots[pivot] = row
                break
    return len(pivots)


def vector(cells, periods):
    column_period, row_period = periods
    return sum(1 << (block * column_period * row_period
                     + column * row_period + row)
               for block, column, row in cells)


def field_product(left, right):
    product = 0
    while right:
        if right & 1:
            product ^= left
        right >>= 1
        left <<= 1
        if left & 4:
            left ^= 7
    return product


def homology_moments(cells, periods):
    omega_powers = (1, 2, 3)
    first_moment = [0] * periods[0]
    second_moment = [0] * periods[1]
    for block, column, row in cells:
        first_coefficient = omega_powers[row % 3]
        second_coefficient = omega_powers[column % 3]
        first_shifts = (2,) if block == 0 else (0, 1)
        second_shifts = (1, 2) if block == 0 else (0,)
        if block == 1:
            first_coefficient = field_product(2, first_coefficient)
        if block == 0:
            second_coefficient = field_product(2, second_coefficient)
        for shift in first_shifts:
            first_moment[(column + shift) % periods[0]] ^= first_coefficient
        for shift in second_shifts:
            second_moment[(row + shift) % periods[1]] ^= second_coefficient
    for moment in (first_moment, second_moment):
        assert all(moment[index] ^ moment[index - 1] ^ moment[index - 2] == 0
                   for index in range(len(moment)))
    return first_moment, second_moment


def audit_original():
    assert syndrome(SEED, (None, None), ORIGINAL_A, ORIGINAL_B) == {
        (column, row + shift)
        for column, row in ((0, 0), (1, 0), (3, 0), (4, -1), (4, 2))
        for shift in (0, 6)
    }
    assert syndrome(DUAL_TILE, (None, None), ORIGINAL_A, ORIGINAL_B, dual=True) == {
        (column + shift, row)
        for column, row in ((-1, 1), (0, -1), (0, 1), (1, 0))
        for shift in (0, 6)
    }
    assert not syndrome(SEED, (None, 6), ORIGINAL_A, ORIGINAL_B)
    assert not syndrome(DUAL_TILE, (6, None), ORIGINAL_A, ORIGINAL_B, dual=True)
    assert len(SEED & DUAL_TILE) == 1
    records = []
    for parameter in (2, 3):
        periods = (24 * parameter + 6, 6)
        centers = [11 + 8 * index for index in range(3 * parameter)]
        dressed = set(SEED)
        for center in centers:
            generator = boundary(center, 0, periods)
            assert len(generator) == 6
            assert not (dressed & generator)
            assert not syndrome(generator, (None, 6), ORIGINAL_A, ORIGINAL_B)
            assert not any(map(any, homology_moments(generator, periods)))
            dressed.symmetric_difference_update(generator)
        dual = {
            (block, column + 6 * repetition, row)
            for block, column, row in DUAL_TILE
            for repetition in range(periods[0] // 6)
        }
        assert not syndrome(dressed, periods, ORIGINAL_A, ORIGINAL_B)
        assert not syndrome(dressed, (None, 6), ORIGINAL_A, ORIGINAL_B)
        assert not syndrome(dual, periods, ORIGINAL_A, ORIGINAL_B, dual=True)
        assert len(dressed & dual) % 2 == 1
        assert len(dressed) == 12 + 18 * parameter < sum(periods) - 6
        gaps = [longest_gap(dressed, axis, periods[axis]) for axis in (0, 1)]
        assert gaps == [3, 1]
        rows = matrix_rows(periods, ORIGINAL_A, ORIGINAL_B, dual=True)
        stabilizer_rank = rank(rows)
        assert rank(rows + [vector(dressed, periods)]) == stabilizer_rank + 1
        assert homology_moments(dressed, periods) == homology_moments(SEED, periods)
        records.append({
            "frame": periods, "centers": centers, "weight": len(dressed),
            "claimed_floor": sum(periods) - 6, "gaps": gaps,
            "pairing": 1, "stabilizer_rank": stabilizer_rank,
            "injective_x_cylinder_lift": True,
        })
    rotated = {(block, row, -column) for block, column, row in SEED}
    rotated_a = tuple((shift_y, -shift_x) for shift_x, shift_y in ORIGINAL_A)
    rotated_b = tuple((shift_y, -shift_x) for shift_x, shift_y in ORIGINAL_B)
    assert set(rotated_a) == set(ORIGINAL_B)
    assert set(rotated_b) == {(-shift_x, -shift_y) for shift_x, shift_y in ORIGINAL_A}
    assert not syndrome(rotated, (6, None), rotated_a, rotated_b)
    wrong_syndrome = syndrome(rotated, (6, None), ORIGINAL_A, ORIGINAL_B)
    assert wrong_syndrome
    six_cycle = {(0, 0, 4), (0, 2, 1), (0, 3, 1),
                 (0, 4, 4), (0, 4, 1), (0, 5, 1)}
    assert not syndrome(six_cycle, (6, 6), ORIGINAL_A, ORIGINAL_B)
    moment_weights = [sum(value != 0 for value in moment)
                      for moment in homology_moments(six_cycle, (6, 6))]
    assert moment_weights == [4, 4]
    return {
        "counterexamples": records,
        "rotated_seed_original_pair_syndrome": sorted(wrong_syndrome),
        "rotated_seed_original_pair_syndrome_weight": len(wrong_syndrome),
        "six_cycle_moment_weights": moment_weights,
    }


def audit_separable():
    first = ((0, 0), (1, 0), (2, 0))
    second = ((0, 0), (0, 1), (0, 2))
    records = []
    for column_period, row_period in ((6, 6), (12, 6), (12, 12), (18, 12)):
        periods = column_period, row_period
        check_rows = matrix_rows(periods, first, second)
        stabilizers = matrix_rows(periods, first, second, dual=True)
        qubits = 2 * column_period * row_period
        assert all((check & generator).bit_count() % 2 == 0
                   for check in check_rows for generator in stabilizers)
        logical_dimension = qubits - rank(check_rows) - rank(stabilizers)
        assert logical_dimension == 8
        witness = {(0, 0, row) for row in range(row_period) if row % 3 != 2}
        assert not syndrome(witness, periods, first, second)
        assert len(witness) == 2 * min(periods) // 3
        assert rank(stabilizers + [vector(witness, periods)]) == rank(stabilizers) + 1
        records.append({"frame": periods, "n": qubits, "k": logical_dimension,
                        "verified_witness_weight": len(witness)})
    return records


def main():
    signal.alarm(55)
    resource.setrlimit(resource.RLIMIT_CPU, (50, 50))
    started = time.monotonic()
    output = {"original": audit_original(), "separable": audit_separable()}
    output["wall_seconds"] = round(time.monotonic() - started, 4)
    output["peak_rss_bytes_macos"] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    assert output["peak_rss_bytes_macos"] < 768 * 1024**2
    print(json.dumps(output, indent=2))


if __name__ == "__main__":
    main()
