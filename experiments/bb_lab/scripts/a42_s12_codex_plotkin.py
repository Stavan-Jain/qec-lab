"""Fixed symbolic Plotkin identities and matched-parity obstruction; no search."""

import json
import resource
import signal
import sys
import time

sys.dont_write_bytecode = True

from a42_s6_codex_identities import add, evaluate_column, monomial, multiply, pairs


def split(polynomial):
    even = {}
    odd = {}
    for (column, row), coefficient in polynomial.items():
        target = odd if row % 2 else even
        target[(column, row // 2)] = coefficient
    return even, odd


def join(even, odd):
    result = {(column, 2 * row): coefficient for (column, row), coefficient in even.items()}
    result.update({(column, 2 * row + 1): coefficient for (column, row), coefficient in odd.items()})
    return result


def pair_weight(vector):
    return sum(map(len, vector))


def boundary(generator, period):
    _, second, first = pairs(period)
    return multiply(second, generator, period), multiply(first, generator, period)


def syndrome(vector, period):
    _, second, first = pairs(period)
    return add(multiply(first, vector[0], period), multiply(second, vector[1], period))


def parameters(half_period):
    identity = monomial(0, 0, half_period)
    shift = monomial(0, 1, half_period)
    omega = monomial(0, 0, half_period, 2)
    even_entry = add(monomial(3, 0, half_period), monomial(0, 1, half_period, 3))
    odd_entry = add(monomial(1, 1, half_period), monomial(2, 1, half_period))
    gamma = add(even_entry, multiply(omega, multiply(shift, odd_entry, half_period), half_period))
    determinant = add(identity, multiply(shift, multiply(odd_entry, odd_entry, half_period), half_period))
    coupling = add(omega, multiply(even_entry, odd_entry, half_period))
    return identity, shift, omega, even_entry, odd_entry, gamma, determinant, coupling


def reduced_coordinates(vector, half_period):
    _, _, omega, _, odd_entry, _, _, _ = parameters(half_period)
    even_first, odd_first = split(vector[0])
    _, odd_second = split(vector[1])
    reduced_first = add(odd_first, multiply(odd_entry, even_first, half_period))
    reduced_second = add(odd_second, multiply(omega, even_first, half_period))
    return even_first, reduced_first, reduced_second


def verify_recursion(vector, half_period):
    _, shift, omega, even_entry, odd_entry, gamma, determinant, _ = parameters(half_period)
    even_first, odd_first = split(vector[0])
    even_second, odd_second = split(vector[1])
    expected_even = add(multiply(even_entry, even_first, half_period),
                        multiply(omega, multiply(shift, odd_first, half_period), half_period),
                        even_second,
                        multiply(shift, multiply(odd_entry, odd_second, half_period), half_period))
    expected_odd = add(multiply(omega, even_first, half_period),
                       multiply(even_entry, odd_first, half_period),
                       multiply(odd_entry, even_second, half_period), odd_second)
    assert syndrome(vector, 2 * half_period) == join(expected_even, expected_odd)
    if expected_even or expected_odd:
        return
    gauge, reduced_first, reduced_second = reduced_coordinates(vector, half_period)
    assert not add(multiply(gamma, reduced_first, half_period),
                   multiply(determinant, reduced_second, half_period))
    recovered_odd_first = add(reduced_first, multiply(odd_entry, gauge, half_period))
    recovered_even_second = add(multiply(even_entry, gauge, half_period),
                                multiply(omega, multiply(shift, reduced_first, half_period), half_period),
                                multiply(shift, multiply(odd_entry, reduced_second, half_period), half_period))
    recovered_odd_second = add(reduced_second, multiply(omega, gauge, half_period))
    assert vector == (join(gauge, recovered_odd_first),
                      join(recovered_even_second, recovered_odd_second))
    cost = sum(map(len, (gauge, recovered_odd_first, recovered_even_second, recovered_odd_second)))
    assert cost == pair_weight(vector)


def verify():
    period = 2
    identity = monomial(0, 0, period)
    shift = monomial(0, 1, period)
    horizontal = add(monomial(1, 0, period), monomial(2, 0, period))
    cyclotomic = add(identity, horizontal)
    base_even = add(monomial(3, 0, period), monomial(0, 0, period, 3))
    residue_factor = add(monomial(2, 0, period), monomial(1, 0, period, 3),
                         monomial(0, 0, period, 3))
    other_factor = add(monomial(2, 0, period), monomial(1, 0, period),
                       monomial(0, 0, period, 3))
    gamma = multiply(add(monomial(1, 0, period), monomial(0, 0, period, 3)), cyclotomic, period)
    coupling = add(monomial(0, 0, period, 2), multiply(base_even, horizontal, period))
    assert coupling == multiply(gamma, residue_factor, period)
    residue = add(monomial(5, 0, period), monomial(0, 0, period, 3))
    assert evaluate_column(residue, 3, period) == identity
    assert not evaluate_column(coupling, 3, period)
    odd_second = multiply(add(monomial(1, 0, period), monomial(0, 0, period, 2)), other_factor, period)
    odd_first = multiply(add(monomial(1, 0, period), monomial(0, 0, period, 2)), odd_second, period)
    assert residue == multiply(residue_factor, odd_second, period)
    assert multiply(multiply(cyclotomic, cyclotomic, period), residue, period) == multiply(coupling, odd_first, period)
    source = (multiply(shift, odd_first, period), add(residue, multiply(shift, odd_second, period)))
    assert not syndrome(source, period)
    assert pair_weight(source) == 10

    multiplier = multiply(monomial(1, 0, period, 2), residue, period)
    alternate_residue = multiply(multiplier, residue, period)
    assert alternate_residue == add(monomial(11, 0, period, 2), monomial(1, 0, period, 3))
    second_source = tuple(multiply(shift, multiply(multiplier, block, period), period) for block in source)
    assert not syndrome(second_source, period)
    square_linear = add(monomial(2, 0, period), monomial(0, 0, period, 2))
    correction_even = multiply(monomial(2, 0, period, 2), multiply(other_factor, square_linear, period), period)
    assert add(alternate_residue, multiply(horizontal, residue, period)) == multiply(coupling, correction_even, period)
    correction_odd = add(multiply(horizontal, correction_even, period), odd_first)
    correction = add(correction_even, multiply(shift, correction_odd, period))
    assert tuple(add(first, second) for first, second in zip(source, second_source)) == boundary(correction, period)
    source_even = tuple(split(block)[0] for block in source)
    second_source_odd = tuple(split(block)[1] for block in second_source)
    assert pair_weight(source_even) == pair_weight(second_source_odd) == 2
    incompatible = tuple(join(even, odd) for even, odd in zip(source_even, second_source_odd))
    assert pair_weight(incompatible) == 4
    assert syndrome(incompatible, period)

    cover_period = 4
    norm = add(monomial(0, 0, cover_period), monomial(0, 2, cover_period))
    first_cover = tuple(multiply(norm, block, cover_period) for block in source)
    second_cover = tuple(multiply(norm, block, cover_period) for block in second_source)
    assert not syndrome(first_cover, cover_period)
    assert not syndrome(second_cover, cover_period)
    assert tuple(add(first, second) for first, second in zip(first_cover, second_cover)) == boundary(multiply(norm, correction, cover_period), cover_period)
    assert sum(len(split(block)[0]) for block in first_cover) == 4
    assert sum(len(split(block)[1]) for block in second_cover) == 4
    verify_recursion(first_cover, 2)
    arbitrary = (add(monomial(-1, 0, cover_period, 2), monomial(2, 3, cover_period)),
                 add(monomial(0, 1, cover_period, 3), monomial(4, 2, cover_period)))
    verify_recursion(arbitrary, 2)
    generator = add(monomial(0, 0, cover_period), monomial(1, 3, cover_period, 2))
    shifted_cover = tuple(add(block, change) for block, change in zip(first_cover, boundary(generator, cover_period)))
    verify_recursion(shifted_cover, 2)
    old_gauge, old_first, old_second = reduced_coordinates(first_cover, 2)
    new_gauge, new_first, new_second = reduced_coordinates(shifted_cover, 2)
    _, half_shift, _, _, odd_entry, half_gamma, half_determinant, _ = parameters(2)
    generator_even, generator_odd = split(generator)
    assert add(old_gauge, new_gauge) == add(generator_even, multiply(half_shift, multiply(odd_entry, generator_odd, 2), 2))
    assert add(old_first, new_first) == multiply(half_determinant, generator_odd, 2)
    assert add(old_second, new_second) == multiply(half_gamma, generator_odd, 2)
    return {
        "status": "PASS",
        "seed_period": period,
        "seed_source_weights": [pair_weight(source), pair_weight(second_source)],
        "same_class_selected_parity_weights": [pair_weight(source_even), pair_weight(second_source_odd)],
        "incompatible_splice_weight": pair_weight(incompatible),
        "incompatible_splice_syndrome_weight": len(syndrome(incompatible, period)),
        "fixed_recursion_check_period": cover_period,
        "scope": "Exact symbolic identities and nonzero syndrome; no distance or support search.",
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
