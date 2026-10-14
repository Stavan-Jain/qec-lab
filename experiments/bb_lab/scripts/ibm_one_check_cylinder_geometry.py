"""Historical partial physical projection certificates, not a full-proof checker.

The full floor is proved in notes/IBM_ONE_CHECK_CYLINDER_TANGENT.md;
the exact modified family is in notes/IBM_MODIFIED_EXACT_DISTANCE_THEOREM.md.
The projection limitations checked here remain valid after that closure.
"""

import resource
import signal
import sys
import time


def add(*polynomials):
    result = set()
    for polynomial in polynomials:
        result.symmetric_difference_update(polynomial)
    return result


def multiply(left, right, period=None):
    result = set()
    for first_column, first_row in left:
        for second_column, second_row in right:
            column = first_column + second_column
            if period is not None:
                column %= period
            result.symmetric_difference_update({(column, first_row + second_row)})
    return result


def power(polynomial, exponent, period=None):
    result = {(0, 0)}
    for index in range(exponent):
        result = multiply(result, polynomial, period)
    return result


def diagonal(polynomial, period=None):
    result = set()
    for column, row in polynomial:
        exponent = column + row
        if period is not None:
            exponent %= period
        result.symmetric_difference_update({(exponent, 0)})
    return result


def detector(chain, period):
    return add(multiply({(1, 0)}, diagonal(chain[0], period), period),
               diagonal(chain[1], period))


def field_multiply(left, right):
    result = 0
    while right:
        if right & 1:
            result ^= left
        right >>= 1
        left <<= 1
        if left & 4:
            left ^= 7
    return result


def dual_add(*elements):
    constant = 0
    tangent = 0
    for first, second in elements:
        constant ^= first
        tangent ^= second
    return constant, tangent


def dual_multiply(left, right):
    return (field_multiply(left[0], right[0]),
            field_multiply(left[0], right[1]) ^ field_multiply(left[1], right[0]))


def dual_power(element, exponent):
    result = (1, 0)
    for index in range(exponent):
        result = dual_multiply(result, element)
    return result


def dual_evaluate(polynomial):
    result = (0, 0)
    for column, row in polynomial:
        term = dual_multiply(dual_power((2, 1), column % 6),
                             dual_power((3, 3), row % 6))
        result = dual_add(result, term)
    return result


def verify_geometry():
    identity = {(0, 0)}
    horizontal = {(1, 0)}
    chi = {(0, 0), (1, 0), (2, 0)}
    upper = {(0, 1), (1, 0)}
    lower = add(upper, identity)
    first_check = {(3, 0), (0, 1), (0, 2)}
    second_check = multiply(chi, lower)
    assert second_check == {(0, 0), (3, 0), (0, 1), (1, 1), (2, 1)}
    assert add(multiply(chi, first_check), multiply(upper, second_check)) == \
        multiply(horizontal, power(chi, 2))
    assert diagonal(first_check) == multiply(horizontal, chi)
    assert diagonal(second_check) == chi
    assert dual_evaluate(chi) == (0, 1)
    assert dual_evaluate(lower) == (0, 2)
    assert dual_evaluate(upper) == (1, 2)
    assert dual_evaluate(multiply(upper, lower)) == (0, 2)
    assert dual_evaluate(first_check) == dual_evaluate(second_check) == (0, 0)
    assert dual_power((2, 1), 6) == (1, 0)
    assert (3, 3) not in {dual_power((2, 1), residue) for residue in range(6)}

    period = 12
    transverse_period = {(0, 0), (period, 0)}
    repeats = {(0, 0), (6, 0)}
    thin = multiply({(0, 0), (2, 0)}, repeats)
    stripe = multiply(chi, thin)
    assert multiply(chi, stripe) == transverse_period
    assert len(thin) == 4 and len(stripe) == 8
    base = (multiply({(-1, 0)}, stripe),
            multiply(multiply({(-1, 0)}, upper), thin))
    nilpotent = set(), multiply(upper, stripe)
    doubled = multiply(upper, base[0]), multiply(upper, base[1])
    records = (
        ((set(), stripe), lower, 8, stripe),
        (base, identity, 16, stripe),
        (nilpotent, multiply(upper, lower), 16, set()),
        (doubled, upper, 24, set()),
    )
    for chain, period_coefficient, expected_weight, expected_detector in records:
        syndrome = add(multiply(first_check, chain[0]), multiply(second_check, chain[1]))
        assert syndrome == multiply(transverse_period, period_coefficient)
        reduced = multiply(identity, chain[0], period), multiply(identity, chain[1], period)
        assert sum(map(len, reduced)) == expected_weight
        assert not add(multiply(first_check, reduced[0], period),
                       multiply(second_check, reduced[1], period))
        assert detector(reduced, period) == expected_detector
        assert expected_detector == multiply(stripe, diagonal(period_coefficient), period)
    print("PASS: ell=12 stripe and three prescribed physical witnesses, weights 8,16,16,24")


def verify_amplification():
    exponent = 6
    period = 48
    input_row = {(0, exponent), (3, exponent)}
    transformed = multiply({(0, 0), (3, 0)}, power({(0, 0), (1, 0)}, exponent), period)
    explicit = multiply({(0, 0), (3, 0)}, {(0, 0), (2, 0), (4, 0), (6, 0)}, period)
    assert transformed == explicit
    assert len(input_row) == 2 and len(transformed) == 8
    assert max(column for column, row in transformed) == 9 < period
    assert multiply({(0, 0), (1, 0)}, {(0, 0), (1, 0), (2, 0)}) == {(0, 0), (3, 0)}
    print("PASS: prescribed noncontracting y=x+1 substitution, input weight 2, output weight 8")


def main():
    signal.alarm(55)
    resource.setrlimit(resource.RLIMIT_CPU, (50, 50))
    started = time.monotonic()
    verify_geometry()
    verify_amplification()
    print("Historical partial scope: 48 detected classes and small-row cycles; "
          "the projection still misses 15 doubled classes")
    print("Full floor now proved in notes/IBM_ONE_CHECK_CYLINDER_TANGENT.md; "
          "full modified-family theorem: notes/IBM_MODIFIED_EXACT_DISTANCE_THEOREM.md")
    peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    peak_bytes = peak if sys.platform == "darwin" else peak * 1024
    assert peak_bytes <= 768 * 1024 ** 2
    print("Elapsed", round(time.monotonic() - started, 6), "seconds; peak MiB",
          round(peak_bytes / 1024 ** 2, 1))


if __name__ == "__main__":
    main()
