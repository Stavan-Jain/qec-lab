"""Fixed polynomial and explicit-vector certificates, not a distance search."""

import resource
import signal
import sys
import time


def add(*polynomials):
    result = set()
    for polynomial in polynomials:
        result.symmetric_difference_update(polynomial)
    return result


def multiply(left, right):
    result = set()
    for left_column, left_row in left:
        for right_column, right_row in right:
            result.symmetric_difference_update({(left_column + right_column, left_row + right_row)})
    return result


def square(polynomial):
    return {(2 * column, 2 * row) for column, row in polynomial}


def reduce_periods(polynomial, width, height):
    result = set()
    for column, row in polynomial:
        result.symmetric_difference_update({(column % width, row % height)})
    return result


def field_multiply(left, right):
    result = 0
    for shift in range(2):
        if right & (1 << shift):
            result ^= left << shift
    if result & 4:
        result ^= 7
    return result


def dual_multiply(left, right):
    return (field_multiply(left[0], right[0]),
            field_multiply(left[0], right[1]) ^ field_multiply(left[1], right[0]))


def dual_power(value, exponent):
    assert exponent >= 0
    result = (1, 0)
    for _ in range(exponent):
        result = dual_multiply(result, value)
    return result


def evaluate_dual(polynomial, column_value, row_value):
    constant = 0
    nilpotent = 0
    for column, row in polynomial:
        value = dual_multiply(dual_power(column_value, column), dual_power(row_value, row))
        constant ^= value[0]
        nilpotent ^= value[1]
    return constant, nilpotent


def main():
    signal.alarm(55)
    resource.setrlimit(resource.RLIMIT_CPU, (50, 50))
    started = time.monotonic()
    identity = {(0, 0)}
    column = {(1, 0)}
    row = {(0, 1)}
    inverse_column = {(-1, 0)}
    cyclotomic = {(0, 0), (1, 0), (2, 0)}
    first_factor = add(row, column)
    second_factor = add(first_factor, identity)
    first_check = {(3, 0), (0, 1), (0, 2)}
    new_check = {(2, 1), (1, 1), (0, 1), (3, 0), (0, 0)}
    old_check = {(0, 0), (1, 3), (2, 3)}
    neighboring_columns = {(1, 0), (2, 0)}
    assert first_check == add(multiply(first_factor, second_factor), multiply(column, cyclotomic))
    assert new_check == multiply(cyclotomic, second_factor)
    assert add(multiply(cyclotomic, first_check), multiply(first_factor, new_check)) == multiply(column, square(cyclotomic))
    assert add(multiply(second_factor, first_check), multiply(column, new_check)) == multiply(first_factor, square(second_factor))
    original_multiplier = multiply(neighboring_columns, add(row, identity))
    original_linear_quotient = {(1, 0), (3, 0)}
    assert old_check == add(multiply(original_multiplier, first_check),
                            multiply(original_linear_quotient, new_check),
                            multiply(cyclotomic, square(cyclotomic)))
    assert multiply(column, original_linear_quotient) == add(identity, square(cyclotomic))
    assert old_check == add(multiply(original_multiplier, first_check),
                            multiply(square(column), new_check),
                            multiply(square(cyclotomic), add(identity, multiply(column, row))))
    assert add(multiply(add(column, identity), column), {(3, 0)}) == multiply(column, cyclotomic)
    old_constant = multiply(cyclotomic, {(0, 0), (1, 0), (3, 0)})
    old_linear = multiply(cyclotomic, original_linear_quotient)
    exceptional = {(7, 0), (6, 0), (4, 0), (2, 0), (0, 0)}
    old_determinant = add(square(old_constant), multiply(old_constant, old_linear),
                          multiply({(3, 0)}, square(old_linear)))
    assert old_determinant == multiply(multiply(square(cyclotomic), cyclotomic), exceptional)
    assert add(old_constant, old_linear) == cyclotomic
    assert evaluate_dual(exceptional, (2, 0), (1, 0)) == (3, 0)
    horizontal_first = multiply(multiply(inverse_column, square(add(column, identity))), cyclotomic)
    horizontal_second = multiply(multiply(inverse_column, square(add(column, identity))), first_factor)
    horizontal_period = {(6, 0), (0, 0)}
    vertical_period = {(0, 6), (0, 0)}
    assert add(multiply(horizontal_first, first_check), multiply(horizontal_second, new_check)) == horizontal_period
    vertical_first = multiply(square(add(row, identity)), add(first_check, horizontal_first))
    vertical_second = multiply(square(add(row, identity)), horizontal_second)
    assert add(multiply(vertical_first, first_check), multiply(vertical_second, new_check)) == vertical_period
    simple_column, simple_row = (2, 0), (2, 0)
    double_column, double_row = (2, 1), (3, 3)
    for column_value, row_value in ((simple_column, simple_row), (double_column, double_row)):
        assert evaluate_dual(first_check, column_value, row_value) == (0, 0)
        assert evaluate_dual(new_check, column_value, row_value) == (0, 0)
        assert evaluate_dual(old_check, column_value, row_value) == (0, 0)
        assert evaluate_dual(horizontal_period, column_value, row_value) == (0, 0)
        assert evaluate_dual(vertical_period, column_value, row_value) == (0, 0)
    assert evaluate_dual(cyclotomic, double_column, double_row) == (0, 1)
    assert evaluate_dual(second_factor, double_column, double_row) == (0, 2)
    assert evaluate_dual(second_factor, simple_column, simple_row) == (1, 0)
    stripe = {(column_index, 0) for start in (0, 3, 6, 9) for column_index in (start, start + 1)}
    fixed_period = {(12, 0), (0, 0)}
    assert multiply(cyclotomic, stripe) == fixed_period
    assert multiply(new_check, stripe) == multiply(second_factor, fixed_period)
    assert reduce_periods(multiply(new_check, stripe), 12, 6) == set()
    assert len(reduce_periods(stripe, 12, 6)) == 8
    assert len(reduce_periods(first_check, 12, 6)) + len(reduce_periods(new_check, 12, 6)) == 8
    peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    peak_bytes = peak if sys.platform == "darwin" else 1024 * peak
    assert peak_bytes <= 768 * 1024 ** 2
    print({"status": "PASS", "plane_dimension_proved_in_note": 6,
           "fixed_frame": [12, 6], "stripe_weight": 8,
           "stripe_period_coefficient_at_simple_point": 1,
           "wall_seconds": round(time.monotonic() - started, 6), "peak_rss_bytes": peak_bytes,
           "scope": "fixed identities and one explicit vector; no exact full distance asserted"})


if __name__ == "__main__":
    main()
