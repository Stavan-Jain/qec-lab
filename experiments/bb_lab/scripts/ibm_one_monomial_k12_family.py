"""Fixed identities and two prescribed witnesses for the one-monomial family."""

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
    for first_column, first_row in left:
        for second_column, second_row in right:
            result.symmetric_difference_update({(first_column + second_column, first_row + second_row)})
    return result


def periodic(polynomial, circumference, height=None):
    result = set()
    for column, row in polynomial:
        actual_row = row if height is None else row % height
        result.symmetric_difference_update({(column % circumference, actual_row)})
    return result


def substitute_row(polynomial, exponent):
    result = set()
    for column, row in polynomial:
        result.symmetric_difference_update({(column + exponent * row, 0)})
    return result


def field_multiply(left, right):
    result = 0
    while right:
        if right & 1:
            result ^= left
        right >>= 1
        left <<= 1
        if left & 8:
            left ^= 13
    return result


def field_power(value, exponent):
    assert value
    if exponent < 0:
        exponent %= 7
    result = 1
    for _ in range(exponent):
        result = field_multiply(result, value)
    return result


def detector(vector, row_exponent, first_shift, circumference):
    return periodic(add(multiply({(first_shift, 0)}, substitute_row(vector[0], row_exponent)),
                        substitute_row(vector[1], row_exponent)), circumference)


def verify_witness(circumference, stripe, expected_weight, first_check, second_check):
    assert all(0 <= column < circumference and row == 0 for column, row in stripe)
    vector = set(), stripe
    assert len(stripe) == expected_weight
    assert periodic(multiply(second_check, stripe), circumference, circumference) == set()
    assert detector(vector, -1, -2, circumference) == stripe
    assert detector(vector, 4, 3, circumference) == stripe
    assert len(periodic(first_check, circumference, circumference)) == 3
    assert len(periodic(second_check, circumference, circumference)) == 3


def main():
    signal.alarm(55)
    resource.setrlimit(resource.RLIMIT_CPU, (50, 50))
    started = time.monotonic()
    identity = {(0, 0)}
    first_check = {(3, 0), (0, 1), (0, 2)}
    second_check = {(0, 0), (1, 0), (5, 0)}
    cyclotomic = {(0, 0), (1, 0), (2, 0)}
    primitive = {(0, 0), (2, 0), (3, 0)}
    stripe_seven = {(0, 0), (2, 0), (3, 0), (4, 0)}
    assert multiply(cyclotomic, primitive) == second_check
    assert multiply(primitive, stripe_seven) == {(0, 0), (7, 0)}
    assert substitute_row(first_check, -1) == multiply({(-2, 0)}, second_check)
    assert substitute_row(first_check, 4) == multiply({(3, 0)}, second_check)
    assert add({(-1, 0)}, {(4, 0)}) == add(identity, multiply({(-1, 0)}, second_check))
    assert substitute_row(first_check, 1) == multiply({(1, 0)}, cyclotomic)
    assert substitute_row(first_check, 2) == multiply({(2, 0)}, cyclotomic)

    alpha = 2
    beta = field_power(alpha, -1)
    gamma = field_power(alpha, 4)
    assert 1 ^ field_power(alpha, 2) ^ field_power(alpha, 3) == 0
    assert alpha != 1 and field_power(alpha, 7) == 1
    assert beta != gamma and beta ^ gamma == 1
    assert field_multiply(beta, gamma) == field_power(alpha, 3)
    assert field_power(beta, 7) == field_power(gamma, 7) == 1
    assert field_multiply(gamma, field_power(beta, -1)) == field_power(alpha, 5)
    assert field_power(alpha, 3) ^ beta ^ field_power(beta, 2) == 0
    assert field_power(alpha, 3) ^ gamma ^ field_power(gamma, 2) == 0

    fixed_generator = {(-3, -2), (2, 1)}
    fixed_boundary = multiply(second_check, fixed_generator), multiply(first_check, fixed_generator)
    assert detector(fixed_boundary, -1, -2, 7) == set()
    assert detector(fixed_boundary, 4, 3, 7) == set()

    verify_witness(7, stripe_seven, 4, first_check, second_check)
    assert multiply(second_check, stripe_seven) == multiply(cyclotomic, {(0, 0), (7, 0)})

    primitive_twenty_one = multiply(stripe_seven, {(0, 0), (7, 0), (14, 0)})
    ordinary_twenty_one = multiply({(0, 0), (1, 0)},
                                  {(0, 0), (3, 0), (6, 0), (9, 0), (12, 0), (15, 0), (18, 0)})
    period_twenty_one = {(0, 0), (21, 0)}
    assert multiply(primitive, primitive_twenty_one) == period_twenty_one
    assert multiply(cyclotomic, ordinary_twenty_one) == period_twenty_one
    assert len(primitive_twenty_one) == 12
    assert len(ordinary_twenty_one) == 14
    mixed_twenty_one = add(primitive_twenty_one, ordinary_twenty_one)
    verify_witness(21, mixed_twenty_one, 10, first_check, second_check)
    assert multiply(second_check, mixed_twenty_one) == multiply(add(primitive, cyclotomic), period_twenty_one)

    peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    peak_bytes = peak if sys.platform == "darwin" else peak * 1024
    assert peak_bytes <= 768 * 1024 ** 2
    print({"status": "PASS", "prescribed_stripe_weights": {"ell7": 4, "ell21": 10},
           "nonzero_detectors_verified": True,
           "wall_seconds": round(time.monotonic() - started, 6), "peak_rss_bytes": peak_bytes,
           "scope": "fixed identities and explicit witnesses; all-size minima are proved analytically"})


if __name__ == "__main__":
    main()
