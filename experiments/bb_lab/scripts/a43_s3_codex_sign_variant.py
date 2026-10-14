"""Bounded exact identities, ranks, and specified witnesses for the sign variant."""

from __future__ import annotations

import json
import resource
import signal
import time


A_STANDARD = frozenset({(3, 0), (0, 1), (0, 2)})
B_STANDARD = frozenset({(0, 3), (1, 0), (2, 0)})
A_LAURENT = frozenset({(0, 0), (0, 1), (3, -1)})
B_PLUS = frozenset({(0, 0), (1, 0), (-1, 3)})
SEED = frozenset({
    (0, 1, 0), (0, 1, 2), (0, 2, 1), (0, 2, 3), (0, 3, 1), (0, 3, 5),
    (1, 3, 1), (1, 3, 5), (1, 4, 0), (1, 4, 1), (1, 4, 2), (1, 4, 5),
})
HELIX_SEED = frozenset({
    (0, 1, 0), (0, 3, 2), (0, 5, 4),
    (1, 0, 1), (1, 2, 3), (1, 4, 5),
})


def polynomial_product(left, right):
    product = set()
    for left_x, left_y in left:
        for right_x, right_y in right:
            cell = (left_x + right_x, left_y + right_y)
            product.symmetric_difference_update({cell})
    return frozenset(product)


def polynomial_power(polynomial, exponent):
    result = frozenset({(0, 0)})
    for _ in range(exponent):
        result = polynomial_product(result, polynomial)
    return result


def polynomial_shift(polynomial, column=0, row=0):
    return frozenset((current_x + column, current_y + row)
                     for current_x, current_y in polynomial)


def coordinate(value, period):
    return value if period is None else value % period


def syndrome(cells, periods, first=A_LAURENT, second=B_PLUS, dual=False):
    result = set()
    supports = (first, second) if dual else (second, first)
    direction = 1 if dual else -1
    for block, column, row in cells:
        for shift_x, shift_y in supports[block]:
            check = (coordinate(column + direction * shift_x, periods[0]),
                     coordinate(row + direction * shift_y, periods[1]))
            result.symmetric_difference_update({check})
    return result


def place(cells, periods, translation=(0, 0)):
    result = set()
    for block, column, row in cells:
        result.symmetric_difference_update({
            (block, coordinate(column + translation[0], periods[0]),
             coordinate(row + translation[1], periods[1]))
        })
    return result


def rows(periods, first=A_LAURENT, second=B_PLUS, dual=False):
    column_period, row_period = periods
    group_size = column_period * row_period
    result = [0] * group_size
    for block in (0, 1):
        for column in range(column_period):
            for row in range(row_period):
                qubit = block * group_size + column * row_period + row
                for check_x, check_y in syndrome(
                    {(block, column, row)}, periods, first, second, dual
                ):
                    result[check_x * row_period + check_y] ^= 1 << qubit
    return result


def rank(vectors):
    pivots = {}
    for vector in vectors:
        while vector:
            pivot = vector.bit_length() - 1
            if pivot in pivots:
                vector ^= pivots[pivot]
            else:
                pivots[pivot] = vector
                break
    return len(pivots)


def packed(cells, periods):
    return sum(1 << (block * periods[0] * periods[1] + column * periods[1] + row)
               for block, column, row in cells)


def verify_ideal():
    cyclotomic = frozenset({(0, 0), (1, 0), (2, 0)})
    linear = frozenset({(1, 1), (0, 1), (1, 0)})
    assert B_STANDARD ^ polynomial_product({(0, 0), (0, 1)}, A_STANDARD) == (
        polynomial_product(cyclotomic, linear)
    )
    coefficient_a = frozenset({(3, 0), (2, 0), (0, 0), (0, 1), (0, 2), (-1, 2)})
    coefficient_b = frozenset({(0, 1), (-1, 1), (-1, 0)})
    assert polynomial_product(coefficient_a, A_STANDARD) ^ polynomial_product(
        coefficient_b, B_STANDARD
    ) == polynomial_power(cyclotomic, 3)
    return {"f_cubed_bezout": {"coefficient_A": sorted(coefficient_a),
                               "coefficient_B": sorted(coefficient_b)}}


def solve_binary(equations):
    pivots = {}
    for coefficients, value in equations:
        while coefficients:
            pivot = coefficients.bit_length() - 1
            if pivot in pivots:
                prior_coefficients, prior_value = pivots[pivot]
                coefficients ^= prior_coefficients
                value ^= prior_value
            else:
                pivots[pivot] = coefficients, value
                break
        if not coefficients and value:
            return None
    solution = 0
    for pivot in sorted(pivots):
        coefficients, value = pivots[pivot]
        if (coefficients & solution).bit_count() % 2 != value:
            solution ^= 1 << pivot
    return solution


def compact_dual(height):
    columns = [(block, column, row) for block in (0, 1)
               for column in range(6) for row in range(height)]
    checks = {}
    pairing = 0
    for index, cell in enumerate(columns):
        for check in syndrome({cell}, (6, None), dual=True):
            checks[check] = checks.get(check, 0) ^ (1 << index)
        if (cell[0], cell[1], cell[2] % 6) in SEED:
            pairing ^= 1 << index
    solution = solve_binary([(coefficients, 0) for coefficients in checks.values()]
                            + [(pairing, 1)])
    if solution is None:
        return None
    support = {cell for index, cell in enumerate(columns) if solution & (1 << index)}
    assert not syndrome(support, (6, None), dual=True)
    return support


def field_product(left, right):
    result = 0
    while right:
        if right & 1:
            result ^= left
        right >>= 1
        left <<= 1
        if left & 4:
            left ^= 7
    return result


def local_product(left, right):
    result = [0, 0, 0]
    for left_index, left_coefficient in enumerate(left):
        for right_index, right_coefficient in enumerate(right):
            if left_index + right_index < 3:
                result[left_index + right_index] ^= field_product(left_coefficient, right_coefficient)
    return tuple(result)


def local_power(value, exponent):
    result = (1, 0, 0)
    for _ in range(exponent):
        result = local_product(result, value)
    return result


def verify_local_algebra():
    local_x = (2, 1, 0)
    local_y = (3, 3, 0)
    for polynomial in (A_STANDARD, B_STANDARD):
        value = [0, 0, 0]
        for column, row in polynomial:
            term = local_product(local_power(local_x, column), local_power(local_y, row))
            value = [current ^ coefficient for current, coefficient in zip(value, term)]
        assert value == [0, 0, 0]
    full_period_residues = []
    for column in range(12):
        for row in range(12):
            diagonal_fixed = (column + row) % 3 == 0
            antidiagonal_fixed = local_product(
                local_power(local_x, column), local_power(local_y, row)
            ) == (1, 0, 0)
            if diagonal_fixed and antidiagonal_fixed:
                full_period_residues.append((column, row))
    assert full_period_residues == [(0, 0)]
    return {"verified_local_length": 3, "full_period_residues_mod_12": full_period_residues}


def main():
    signal.alarm(55)
    resource.setrlimit(resource.RLIMIT_CPU, (50, 50))
    started = time.monotonic()
    output = {"identities": verify_ideal(), "local_algebra": verify_local_algebra(), "frames": []}
    direct_dual = compact_dual(6)
    assert direct_dual is not None
    output["linear_solve_dual"] = None if direct_dual is None else sorted(direct_dual)
    output["seed_plane_syndrome"] = sorted(syndrome(SEED, (None, None)))
    output["dual_plane_syndrome"] = sorted(syndrome(direct_dual, (None, None), dual=True))
    assert syndrome(SEED, (None, None)) == {
        (column, row + offset)
        for column, row in ((0, 0), (1, 0), (2, -3), (2, -1), (3, -2), (4, -2), (4, -1))
        for offset in (0, 6)
    }
    assert syndrome(direct_dual, (None, None), dual=True) == {
        (column + offset, row)
        for column, row in ((-1, 4), (-1, 6), (0, 1), (1, 0), (1, 1), (2, 3), (2, 4))
        for offset in (0, 6)
    }
    assert not syndrome(SEED, (None, 6))
    assert not syndrome(HELIX_SEED, (6, 6))
    dual_tile = {(block, (-row - 1) % 6, -column - 1)
                 for block, column, row in SEED}
    assert not syndrome(dual_tile, (6, None), dual=True)
    dual_tile = place(dual_tile, (6, None), (0, 6))
    output["dual_tile"] = sorted(dual_tile)
    output["dual_seed_pairings"] = []
    for shift in range(6):
        translated_dual = place(dual_tile, (6, 6), (shift, 0))
        output["dual_seed_pairings"].append(len(SEED & translated_dual) % 2)
    for periods in ((3, 3), (6, 6), (12, 6), (12, 12), (24, 12), (36, 18)):
        syndrome_rows = rows(periods)
        boundary_rows = rows(periods, dual=True)
        syndrome_rank = rank(syndrome_rows)
        boundary_rank = rank(boundary_rows)
        assert syndrome_rank == boundary_rank
        expected_k = 8 if min(periods) == 3 else (16 if all(
            period % 12 == 0 for period in periods) else 12)
        dimension = 2 * periods[0] * periods[1] - syndrome_rank - boundary_rank
        assert dimension == expected_k
        record = {"frame": periods, "k": dimension}
        if periods[0] >= 12:
            witness = set()
            for repetition in range(periods[1] // 6):
                witness.symmetric_difference_update(place(SEED, periods, (0, 6 * repetition)))
            assert not syndrome(witness, periods)
            record["stack_weight"] = len(witness)
            record["stack_nontrivial"] = rank(boundary_rows + [packed(witness, periods)]) > boundary_rank
            assert record["stack_nontrivial"]
            direct_tiled = set()
            for repetition in range(periods[0] // 6):
                direct_tiled.symmetric_difference_update(place(direct_dual, periods, (6 * repetition, 0)))
            assert not syndrome(direct_tiled, periods, dual=True)
            assert len(witness & direct_tiled) % 2 == 1
            record["direct_dual_pairing"] = 1
            dual_pairings = []
            for shift in range(6):
                tiled_dual = set()
                for repetition in range(periods[0] // 6):
                    tiled_dual.symmetric_difference_update(place(
                        dual_tile, periods, (shift + 6 * repetition, 0)))
                assert not syndrome(tiled_dual, periods, dual=True)
                dual_pairings.append(len(witness & tiled_dual) % 2)
            record["dual_pairings"] = dual_pairings
        output["frames"].append(record)
    output["wall_seconds"] = round(time.monotonic() - started, 4)
    output["peak_rss_bytes_macos"] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    assert output["peak_rss_bytes_macos"] < 768 * 1024**2
    print(json.dumps(output, indent=2))


if __name__ == "__main__":
    main()
