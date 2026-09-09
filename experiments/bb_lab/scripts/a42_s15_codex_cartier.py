"""Fixed Cartier identities and three-basis actions, not a distance search."""

import resource
import signal
import sys
import time

sys.dont_write_bytecode = True

from a42_s6_codex_identities import add, field_mul, monomial, multiply, pairs, square


CHANNELS = ((0, 0), (1, 0), (0, 1), (1, 1))


def cartier(polynomial, column_parity, row_parity, target_period):
    return add(*(
        monomial((column - column_parity) // 2,
                 (row - row_parity) // 2,
                 target_period, field_mul(coefficient, coefficient))
        for (column, row), coefficient in polynomial.items()
        if column % 2 == column_parity and row % 2 == row_parity
    ))


def cycle_map(vector, channel, target_period):
    source_period = 2 * target_period
    _, second_check, first_check = pairs(source_period)
    return (
        cartier(multiply(second_check, vector[0], source_period),
                *channel, target_period),
        cartier(multiply(first_check, vector[1], source_period),
                *channel, target_period),
    )


def syndrome(vector, period):
    _, second_check, first_check = pairs(period)
    return add(multiply(first_check, vector[0], period),
               multiply(second_check, vector[1], period))


def boundary(generator, period):
    _, second_check, first_check = pairs(period)
    return (multiply(second_check, generator, period),
            multiply(first_check, generator, period))


def vector_weight(vector):
    return len(vector[0]) + len(vector[1])


def label_polynomial(constant, linear):
    return {(column, 0): coefficient for column, coefficient
            in ((0, constant), (1, linear)) if coefficient}


def scaled_sum(polynomials, coefficients, period):
    return add(*(multiply(monomial(0, 0, period, coefficient), polynomial, period)
                 for polynomial, coefficient in zip(polynomials, coefficients)
                 if coefficient))


def height_two_basis():
    repeated = {(0, 0): 1, (0, 1): 1}
    one_plus_column = {(0, 0): 1, (1, 0): 1}
    first_basis = (repeated, multiply(one_plus_column, repeated, 2))
    second_basis = tuple(multiply(monomial(1, 0, 2), block, 2)
                         for block in first_basis)
    third_basis = (
        {(1, 0): 1, (0, 0): 2},
        {(0, 0): 1, (2, 1): 1, (1, 1): 3, (0, 1): 3},
    )
    return first_basis, second_basis, third_basis


def fixed_chain_checks():
    source_period = 4
    target_period = 2
    _, source_second, source_first = pairs(source_period)
    _, target_second, target_first = pairs(target_period)
    generator = {(-3, 3): 2, (2, 0): 3, (0, 1): 1}
    vector = ({(-1, 3): 3, (4, 0): 2}, {(2, 1): 1, (-4, 2): 2})
    source_product = multiply(source_first, source_second, source_period)
    for channel in CHANNELS:
        image = cycle_map(vector, channel, target_period)
        assert syndrome(image, target_period) == cartier(
            multiply(source_product, syndrome(vector, source_period), source_period),
            *channel, target_period)
        assert cycle_map(boundary(generator, source_period), channel, target_period) == (
            multiply(target_second, cartier(generator, *channel, target_period), target_period),
            multiply(target_first, cartier(generator, *channel, target_period), target_period),
        )
        assert cartier(multiply(square(target_first, source_period), generator, source_period),
                       *channel, target_period) == multiply(
            target_first, cartier(generator, *channel, target_period), target_period)
        unwrapped = {(column, row + 8): coefficient
                     for (column, row), coefficient in generator.items()}
        assert cartier(unwrapped, *channel, target_period) == cartier(
            generator, *channel, target_period)
    assert sum(vector_weight(cycle_map(vector, channel, target_period))
               for channel in CHANNELS) == (
        len(multiply(source_second, vector[0], source_period))
        + len(multiply(source_first, vector[1], source_period)))


def fixed_action_checks():
    expected_columns = (
        ((1, 1), (0, 1), (3, 0)),
        ((1, 0), (1, 1), (1, 0)),
        ((1, 1), (0, 1), (0, 2)),
        ((1, 0), (1, 1), (3, 1)),
    )
    kernels = ((2, 2, 1), (1, 0, 1), (0, 3, 1), (3, 1, 1))
    basis = height_two_basis()
    one_plus_column = {(0, 0): 1, (1, 0): 1}
    for vector in basis:
        assert syndrome(vector, 2) == {}
    for channel, columns in zip(CHANNELS, expected_columns):
        for vector, column in zip(basis, columns):
            image = cycle_map(vector, channel, 1)
            expected_first = label_polynomial(*column)
            assert image == (expected_first, multiply(one_plus_column, expected_first, 1))
            transfer = {(0, 0): 1, (0, 2): 1}
            lifted = tuple(multiply(transfer, block, 4) for block in vector)
            assert syndrome(lifted, 4) == {}
            assert cycle_map(lifted, channel, 2) == tuple(
                multiply({(0, 0): 1, (0, 1): 1}, block, 2) for block in image)
        determinant = field_mul(columns[0][0], columns[1][1]) ^ field_mul(
            columns[0][1], columns[1][0])
        assert determinant
    for kernel_index, kernel in enumerate(kernels):
        images = []
        for columns in expected_columns:
            image = scaled_sum([label_polynomial(*column) for column in columns],
                               [field_mul(coefficient, coefficient) for coefficient in kernel], 1)
            images.append(image)
        assert [index for index, image in enumerate(images) if not image] == [kernel_index]
    source_coordinates = ((1, 0, 0), (0, 1, 0), (0, 0, 1))
    for coordinates in source_coordinates:
        coefficient_first, coefficient_second, coefficient_third = coordinates
        for channel_index, columns in enumerate(expected_columns):
            image = scaled_sum([label_polynomial(*column) for column in columns], coordinates, 1)
            at_omega = 0
            at_omega_squared = 0
            for (column, _), coefficient in image.items():
                at_omega ^= field_mul(coefficient, 1 if column == 0 else 2)
                at_omega_squared ^= field_mul(coefficient, 1 if column == 0 else 3)
            first_value = coefficient_first ^ field_mul(3, coefficient_second) ^ coefficient_third
            second_value = coefficient_first ^ field_mul(2, coefficient_second)
            expected_first = field_mul(3, first_value) if channel_index in (0, 2) else first_value
            expected_second = (
                field_mul(2, second_value) ^ field_mul(3, coefficient_third),
                second_value ^ coefficient_third,
                field_mul(2, second_value) ^ coefficient_third,
                second_value,
            )[channel_index]
            assert (at_omega, at_omega_squared) == (expected_first, expected_second)


def fixed_cost_obstructions():
    source_period = 16
    target_period = 8
    _, second_check, first_check = pairs(source_period)
    generator = multiply(first_check, second_check, source_period)
    vector = boundary(generator, source_period)
    assert syndrome(vector, source_period) == {}
    assert tuple(map(len, vector)) == (9, 9)
    assert (len(multiply(second_check, vector[0], source_period)),
            len(multiply(first_check, vector[1], source_period))) == (27, 27)
    assert sum(vector_weight(cycle_map(vector, channel, target_period))
               for channel in CHANNELS) == 3 * vector_weight(vector) == 54
    expected_factors = (
        {(0, 1): 2, (1, 2): 3},
        {(1, 0): 1, (0, 2): 3},
        {(0, 0): 3, (2, 1): 1, (1, 2): 2},
        {(2, 1): 1, (0, 2): 2},
    )
    for channel, expected in zip(CHANNELS, expected_factors):
        assert cartier(generator, *channel, target_period) == expected
    nontrivial_vector = height_two_basis()[2]
    images = [cycle_map(nontrivial_vector, channel, 1) for channel in CHANNELS]
    assert vector_weight(nontrivial_vector) == 6
    assert [vector_weight(image) for image in images] == [3, 3, 3, 5]
    repeated = {(0, row): 1 for row in range(4)}
    canonical = (repeated, multiply({(0, 0): 1, (1, 0): 1}, repeated, 4))
    assert [vector_weight(cycle_map(canonical, channel, 2))
            for channel in CHANNELS] == [8, 6, 8, 6]
    print({"sharp_boundary_cost": "54/18=3", "height_two_nontrivial_output_costs": [3, 3, 3, 5]})


def main():
    started = time.monotonic()
    signal.alarm(55)
    resource.setrlimit(resource.RLIMIT_CPU, (50, 50))
    fixed_chain_checks()
    fixed_action_checks()
    fixed_cost_obstructions()
    peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    peak_bytes = peak if sys.platform == "darwin" else 1024 * peak
    assert peak_bytes <= 768 * 1024 ** 2
    print({"elapsed_seconds": round(time.monotonic() - started, 6), "peak_rss_bytes": peak_bytes})
    print("PASS: fixed identities and basis actions only; no sharp pure floor asserted.")


if __name__ == "__main__":
    main()
