"""Fixed algebraic anchors for the all-size cylinder parity proof; no census."""

import json
import resource
import signal
import sys
import time

from a42_lib import F2k, pdeg, pdivmod, pgcd, pmod, pmul, sylvester_resultant_y


def poly_power(polynomial, exponent, modulus=None):
    result = 1
    while exponent:
        if exponent & 1:
            result = pmul(result, polynomial)
            if modulus is not None:
                result = pmod(result, modulus)
        polynomial = pmul(polynomial, polynomial)
        if modulus is not None:
            polynomial = pmod(polynomial, modulus)
        exponent >>= 1
    return result


def add(*polynomials):
    result = set()
    for polynomial in polynomials:
        result.symmetric_difference_update(polynomial)
    return frozenset(result)


def multiply(first, second):
    result = set()
    for first_x, first_y in first:
        for second_x, second_y in second:
            term = (first_x + second_x, first_y + second_y)
            if term in result:
                result.remove(term)
            else:
                result.add(term)
    return frozenset(result)


def geometric_sum(exponent, axis):
    powers = range(exponent) if exponent >= 0 else range(exponent, 0)
    return frozenset((power, 0) if axis == 0 else (0, power) for power in powers)


def verify():
    scalar = 1 ^ (1 << 4) ^ (1 << 5)
    linear = (1 << 1) ^ (1 << 2) ^ (1 << 4) ^ (1 << 5)
    primary = 7
    exceptional = 213
    assert pdivmod(scalar, primary) == (11, 0)
    assert pdivmod(linear, primary) == (10, 0)
    assert pgcd(scalar, linear) == primary
    determinant = pmul(scalar, scalar) ^ pmul(scalar, linear) ^ pmul(8, pmul(linear, linear))
    assert determinant == pmul(poly_power(primary, 3), exceptional)
    assert determinant == sylvester_resultant_y()
    smith_second, remainder = pdivmod(determinant, primary)
    assert remainder == 0
    assert smith_second == pmul(poly_power(primary, 2), exceptional)
    assert pdeg(primary) + pdeg(smith_second) == 13
    assert poly_power(2, 128, exceptional) == 2
    assert pgcd(exceptional, (1 << 2) ^ (1 << 1)) == 1
    assert pgcd(exceptional, primary) == 1
    assert pdivmod((1 << 6) ^ 1, poly_power(primary, 2))[1] == 0

    field = F2k(exceptional)
    x_value = 2
    y_value = field.mul(field.inv(pmod(linear, exceptional)), pmod(scalar, exceptional))
    assert y_value == field.pow(x_value, 91)
    assert field.pow(x_value, 3) ^ y_value ^ field.mul(y_value, y_value) == 0
    assert field.mul(x_value ^ field.mul(x_value, x_value), field.pow(y_value, 3)) == 1
    assert x_value not in (0, 1) and y_value not in (0, 1)
    assert field.pow(x_value, 127) == field.pow(y_value, 127) == 1

    unit = frozenset({(0, 0)})
    first_period = frozenset({(0, 0), (1, 0)})
    second_period = frozenset({(0, 0), (0, 1)})
    winding_checks = 0
    for first_winding in range(-3, 4):
        for second_winding in range(-3, 4):
            coefficient_first = multiply(frozenset({(0, second_winding)}), geometric_sum(first_winding, 0))
            coefficient_second = geometric_sum(second_winding, 1)
            target = add(unit, frozenset({(first_winding, second_winding)}))
            assert add(multiply(coefficient_first, first_period), multiply(coefficient_second, second_period)) == target
            assert len(coefficient_first) % 2 == first_winding % 2
            assert len(coefficient_second) % 2 == second_winding % 2
            winding_checks += 1

    cover_second = multiply(second_period, second_period)
    pullback_degree_one = (second_period, unit)
    projection_degree_one = (unit, second_period)
    assert tuple(multiply(value, diagonal) for value, diagonal in zip((first_period, cover_second), pullback_degree_one)) == tuple(multiply(second_period, value) for value in (first_period, second_period))
    assert tuple(multiply(diagonal, value) for diagonal, value in zip(pullback_degree_one, (second_period, first_period))) == (cover_second, first_period)
    assert tuple(multiply(value, diagonal) for value, diagonal in zip((first_period, second_period), projection_degree_one)) == (first_period, cover_second)
    assert tuple(multiply(diagonal, value) for diagonal, value in zip(projection_degree_one, (cover_second, first_period))) == tuple(multiply(second_period, value) for value in (second_period, first_period))

    counts = {}
    for dimension in (6, 13):
        nonzero = (1 << dimension) - 1
        counts[str(dimension)] = {"cylinder_classes": 3 * nonzero, "mixed_classes": nonzero * (nonzero - 1)}
        assert sum(counts[str(dimension)].values()) == (1 << (2 * dimension)) - 1
    return {"status": "PASS", "plane_module_dimension": 13, "ordinary_dimension": 6,
            "exceptional_dimension": 7, "winding_identity_checks": winding_checks,
            "cover_chain_map_identities": 4, "counts": counts,
            "scope": "fixed algebraic inputs; universal theorem is in the proof note"}


def timeout_handler(signum, frame):
    raise TimeoutError("55-second diagnostic limit")


def main():
    signal.signal(signal.SIGALRM, timeout_handler)
    signal.alarm(55)
    started = time.monotonic()
    report = verify()
    report["elapsed_seconds"] = time.monotonic() - started
    peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    peak_bytes = peak if sys.platform == "darwin" else 1024 * peak
    if peak_bytes > 768 * 1024 * 1024:
        raise MemoryError("768 MB diagnostic limit")
    report["peak_rss_bytes"] = peak_bytes
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
