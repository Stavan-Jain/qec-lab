"""Fixed q=4 identity checks; not a distance search or a tower verifier."""

import json
import resource
import signal
import sys
import time

sys.dont_write_bytecode = True

from a42_s6_codex_identities import add, field_mul, field_pow, monomial, multiply, pairs, power


def truncated_add(*polynomials):
    result = [0] * 4
    for polynomial in polynomials:
        for exponent, coefficient in enumerate(polynomial):
            result[exponent] ^= coefficient
    return result


def truncated_multiply(first, second):
    result = [0] * 4
    for first_exponent, first_coefficient in enumerate(first):
        for second_exponent, second_coefficient in enumerate(second):
            if first_exponent + second_exponent < 4:
                result[first_exponent + second_exponent] ^= field_mul(
                    first_coefficient, second_coefficient
                )
    return result


def u_coefficients(polynomial, period):
    result = [{} for _ in range(period)]
    for (column, row), coefficient in polynomial.items():
        for exponent in range(row + 1):
            if exponent & row == exponent:
                value = result[exponent].get(column, 0) ^ coefficient
                if value:
                    result[exponent][column] = value
                else:
                    result[exponent].pop(column, None)
    return result


def verify():
    scalar = [0, 3, 3, 0]
    linear = [0, 1, 0, 0]
    inverse_cube = [1, 1, 0, 0]
    determinant = truncated_add(
        truncated_multiply(scalar, scalar),
        truncated_multiply(scalar, linear),
        truncated_multiply(truncated_multiply(linear, linear), inverse_cube),
    )
    assert determinant == [0, 0, 0, 2]

    period = 4
    _, second, first = pairs(period)
    identity = monomial(0, 0, period)
    parameter = add(identity, monomial(0, 1, period))
    repeated = {(0, row): 1 for row in range(period)}
    canonical = (
        repeated,
        add(repeated, multiply(monomial(1, 0, period), repeated, period)),
    )
    local_generator = monomial(0, period - 2, period)
    upper_generator = power(parameter, period - 2, period)
    lower_generator = add(local_generator, upper_generator)
    local_boundary = (
        multiply(second, local_generator, period),
        multiply(first, local_generator, period),
    )
    vector = tuple(add(block, boundary) for block, boundary in zip(canonical, local_boundary))
    normalized = (
        add(vector[0], multiply(second, lower_generator, period)),
        add(vector[1], multiply(first, lower_generator, period)),
    )
    expected_normalized = (
        add(canonical[0], multiply(second, upper_generator, period)),
        add(canonical[1], multiply(first, upper_generator, period)),
    )
    assert normalized == expected_normalized
    assert sum(map(len, canonical)) == 12
    assert sum(map(len, vector)) == 14
    assert sum(map(len, normalized)) == 16
    assert not add(multiply(first, vector[0], period), multiply(second, vector[1], period))
    assert not add(multiply(first, normalized[0], period), multiply(second, normalized[1], period))
    assert u_coefficients(local_generator, period)[period - 2] == {0: 1}
    assert not any(u_coefficients(lower_generator, period)[period - 2:])
    assert not any(u_coefficients(normalized[0], period)[:period - 2])
    assert not any(u_coefficients(normalized[1], period)[:period - 2])

    expected_base_first = {(0, 1): 1, (1, 1): 1, (2, 1): 1}
    expected_base_second = {(3, 0): 1, (1, 0): 1, (0, 0): 2, (1, 1): 1, (0, 1): 3}
    assert normalized[0] == multiply(upper_generator, expected_base_first, period)
    assert normalized[1] == multiply(upper_generator, expected_base_second, period)

    tangent_base_first = {(1, 0): 1, (0, 0): 2}
    tangent_base_second = {(0, 0): 1, (2, 1): 1, (1, 1): 3, (0, 1): 3}
    tangent = (
        multiply(upper_generator, tangent_base_first, period),
        multiply(upper_generator, tangent_base_second, period),
    )
    tangent_vector = tuple(add(block, boundary) for block, boundary in zip(tangent, local_boundary))
    tangent_normalized = (
        add(tangent_vector[0], multiply(second, lower_generator, period)),
        add(tangent_vector[1], multiply(first, lower_generator, period)),
    )
    assert sum(map(len, tangent)) == 12
    assert sum(map(len, tangent_vector)) == 15
    assert sum(map(len, tangent_normalized)) == 18
    assert not add(multiply(first, tangent[0], period), multiply(second, tangent[1], period))
    assert not add(multiply(first, tangent_vector[0], period), multiply(second, tangent_vector[1], period))
    assert not add(multiply(first, tangent_normalized[0], period), multiply(second, tangent_normalized[1], period))
    expected_tangent_first = {(1, 0): 1, (0, 0): 3, (1, 1): 1, (2, 1): 1}
    expected_tangent_second = {(3, 0): 1, (0, 0): 2, (2, 1): 1, (1, 1): 3, (0, 1): 1}
    assert tangent_normalized[0] == multiply(upper_generator, expected_tangent_first, period)
    assert tangent_normalized[1] == multiply(upper_generator, expected_tangent_second, period)
    tangent_leading_first = u_coefficients(tangent[0], period)[period - 2]
    assert tangent_leading_first == {1: 1, 0: 2}
    root_value = 0
    for column, coefficient in tangent_leading_first.items():
        root_value ^= field_mul(coefficient, field_pow(3, column))
    assert root_value == 1
    assert 1 ^ 3 ^ field_mul(3, 3) == 0

    toy_first = power(parameter, period - 2, period)
    toy_logical = power(parameter, period - 1, period)
    assert multiply(parameter, toy_first, period) == toy_logical
    assert not multiply(parameter, toy_logical, period)
    assert len(toy_first) == period // 2
    assert len(toy_logical) == period
    return {
        "status": "PASS",
        "fixed_period": period,
        "plane_determinant_leading_coefficient": determinant,
        "actual_vector_weight": sum(map(len, vector)),
        "normalized_vector_weight": sum(map(len, normalized)),
        "maximal_layer_vector_weight": sum(map(len, tangent_vector)),
        "maximal_layer_normalized_weight": sum(map(len, tangent_normalized)),
        "scope": "Fixed symbolic conventions only; universal proofs and nontriviality are in the note.",
    }


def timeout_handler(signum, frame):
    raise TimeoutError("55-second wall limit")


def main():
    signal.signal(signal.SIGALRM, timeout_handler)
    signal.alarm(55)
    resource.setrlimit(resource.RLIMIT_CPU, (50, 50))
    started = time.perf_counter()
    report = verify()
    report["elapsed_seconds"] = time.perf_counter() - started
    peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    peak_bytes = peak if sys.platform == "darwin" else 1024 * peak
    if peak_bytes > 768 * 1024 * 1024:
        raise MemoryError("768-MiB RSS limit")
    report["peak_rss_bytes"] = peak_bytes
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
