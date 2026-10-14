"""Fixed polynomial certificate of mixedness for the dyadic square witnesses."""

import resource
import signal
import sys
import time

from a40_parity_tor_check import add, multiply
from a42_lib import F2k


ALPHA = frozenset({(0, 0), (0, -1), (-3, 1)})
BETA = frozenset({(0, 0), (-1, 0), (1, 3)})


def polynomial_power(polynomial, exponent):
    result = frozenset({(0, 0)})
    while exponent:
        if exponent & 1:
            result = multiply(result, polynomial)
        polynomial = multiply(polynomial, polynomial)
        exponent >>= 1
    return result


def evaluate(polynomial, field, first_value, second_value):
    result = 0
    for horizontal, vertical in polynomial:
        result ^= field.mul(field.pow(first_value, horizontal), field.pow(second_value, vertical))
    return result


def verify():
    first_flux = frozenset({(-3, 2), (-1, 4), (0, 7)})
    second_flux = frozenset({(1, 0), (0, 1)})
    field = F2k(7)
    first_value, second_value = 2, 3
    assert evaluate(ALPHA, field, first_value, second_value) == 0
    assert evaluate(BETA, field, first_value, second_value) == 0
    assert evaluate(first_flux, field, first_value, second_value) == 3
    assert evaluate(second_flux, field, first_value, second_value) == 1
    results = []
    for height, radius in enumerate((1, 2, 4, 8)):
        norm = frozenset({(0, 0), (2 * radius, 2 * radius), (4 * radius, 4 * radius)})
        first = multiply(multiply(norm, frozenset({(radius, 0)})), polynomial_power(BETA, radius - 1))
        second = multiply(multiply(norm, frozenset({(0, radius)})), polynomial_power(ALPHA, radius - 1))
        syndrome = add(multiply(BETA, first), multiply(ALPHA, second))
        expected = add(multiply(frozenset({(0, 0), (6 * radius, 0)}), polynomial_power(first_flux, radius)),
                       multiply(frozenset({(0, 0), (0, 6 * radius)}), polynomial_power(second_flux, radius)))
        assert syndrome == expected
        first_coordinate = evaluate(polynomial_power(first_flux, radius), field, first_value, second_value)
        second_coordinate = evaluate(polynomial_power(second_flux, radius), field, first_value, second_value)
        assert first_coordinate in (2, 3) and second_coordinate == 1
        reduced_blocks = []
        for block in (first, second):
            reduced = set()
            for horizontal, vertical in block:
                reduced.symmetric_difference_update({(horizontal % (6 * radius), vertical % (6 * radius))})
            reduced_blocks.append(reduced)
        assert sum(map(len, reduced_blocks)) == 6 * 3 ** height
        results.append({"radius": radius, "period_flux_evaluations_F4": [first_coordinate, second_coordinate],
                        "witness_weight": sum(map(len, reduced_blocks))})
    return results


def main():
    signal.alarm(55)
    started = time.monotonic()
    results = verify()
    peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    peak_bytes = peak if sys.platform == "darwin" else 1024 * peak
    assert peak_bytes < 768 * 1024 ** 2
    print({"status": "PASS", "fixed_identity_replays": results,
           "elapsed_seconds": time.monotonic() - started, "peak_rss_bytes": peak_bytes,
           "scope": "fixed identity and field evaluation; universal proof uses Frobenius"})


if __name__ == "__main__":
    main()
