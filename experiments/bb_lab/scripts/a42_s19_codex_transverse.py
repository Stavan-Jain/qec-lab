"""Fixed period identities for an image-class counterexample, not a distance search."""

import resource
import signal
import sys
import time

sys.dont_write_bytecode = True

from a42_s6_codex_identities import (
    add,
    evaluate_column,
    field_mul,
    monomial,
    multiply,
    pairs,
    power,
    square,
)


SENTINEL_PERIOD = 16


def scale(polynomial, coefficient, period=SENTINEL_PERIOD):
    return multiply(monomial(0, 0, period, coefficient), polynomial, period)


def fold_columns(polynomial):
    result = {}
    for (column, row), coefficient in polynomial.items():
        position = column % 6, row
        result[position] = result.get(position, 0) ^ coefficient
        if not result[position]:
            del result[position]
    return result


def to_local_parameter(polynomial):
    one_plus_parameter = {(0, 0): 1, (0, 1): 1}
    return add(*(scale(power(one_plus_parameter, row, SENTINEL_PERIOD), coefficient)
                 for (column, row), coefficient in polynomial.items()))


def vertical_detector(first, second, period):
    _, second_check, first_check = pairs(SENTINEL_PERIOD)
    first_evaluation = evaluate_column(first, 2, SENTINEL_PERIOD)
    second_evaluation = evaluate_column(second, 2, SENTINEL_PERIOD)
    syndrome = add(
        multiply(evaluate_column(first_check, 2, SENTINEL_PERIOD),
                 first_evaluation, SENTINEL_PERIOD),
        multiply(evaluate_column(second_check, 2, SENTINEL_PERIOD),
                 second_evaluation, SENTINEL_PERIOD),
    )
    local = to_local_parameter(syndrome)
    assert all(row >= period for column, row in local)
    return local.get((0, period), 0)


def local_module_check():
    identity = {(0, 0): 1}
    parameter = {(0, 1): 1}
    local_t = add(identity, parameter)
    inverse_t_cubed_jet = local_t

    def jet(polynomial):
        return {position: coefficient for position, coefficient in polynomial.items()
                if position[1] < 4}

    assert jet(multiply(power(local_t, 3, SENTINEL_PERIOD),
                        inverse_t_cubed_jet, SENTINEL_PERIOD)) == identity
    first_entry = jet(add(inverse_t_cubed_jet, scale(local_t, 2),
                          scale(square(local_t, SENTINEL_PERIOD), 3)))
    second_entry = add(identity, inverse_t_cubed_jet)
    assert first_entry == {(0, 1): 3, (0, 2): 3}
    assert second_entry == parameter
    determinant = jet(add(
        square(first_entry, SENTINEL_PERIOD),
        multiply(first_entry, second_entry, SENTINEL_PERIOD),
        multiply(square(second_entry, SENTINEL_PERIOD),
                 inverse_t_cubed_jet, SENTINEL_PERIOD),
    ))
    assert determinant == {(0, 3): 2}


def plane_period_identity():
    _, second_check, first_check = pairs(SENTINEL_PERIOD)
    local_parameter = {(0, 0): 1, (0, 1): 1}
    vertical_square = square(local_parameter, SENTINEL_PERIOD)
    vertical_fourth = square(vertical_square, SENTINEL_PERIOD)
    transverse_period = {(0, 0): 1, (6, 0): 1}
    assert add(square(first_check, SENTINEL_PERIOD),
               scale(vertical_square, 3), scale(vertical_fourth, 2)) == transverse_period
    first_lift = multiply(vertical_square, first_check, SENTINEL_PERIOD)
    vertical_coefficient = add({(0, 0): 3}, scale(vertical_square, 2))
    transverse_coefficient = vertical_square
    assert multiply(first_check, first_lift, SENTINEL_PERIOD) == add(
        multiply(vertical_fourth, vertical_coefficient, SENTINEL_PERIOD),
        multiply(transverse_period, transverse_coefficient, SENTINEL_PERIOD),
    )
    assert add(*({(0, 0): coefficient}
                 for coefficient in evaluate_column(vertical_coefficient, 2, SENTINEL_PERIOD).values())) == {(0, 0): 3}
    for check in (first_check, second_check):
        assert add(*({(0, 0): coefficient}
                     for coefficient in evaluate_column(check, 2, SENTINEL_PERIOD).values())) == {}


def height_four_witness():
    period = 4
    _, second_check, first_check = pairs(period)
    transfer = {(0, 0): 1, (0, 2): 1}
    transverse_period = {(0, 0): 1, (6, 0): 1}
    first = multiply(transfer, first_check, period)
    assert first == {
        (3, 0): 1, (3, 2): 1,
        (0, 0): 3, (0, 2): 3,
        (0, 1): 2, (0, 3): 2,
    }
    unfolded_syndrome = multiply(first_check, first, period)
    assert unfolded_syndrome == multiply(transfer, transverse_period, period)
    assert unfolded_syndrome and not fold_columns(unfolded_syndrome)
    assert len(first) == 6 < 3 * period
    assert vertical_detector(first, {}, period) == 3
    generator = {(5, 3): 2}
    first_boundary = fold_columns(multiply(second_check, generator, period))
    second_boundary = fold_columns(multiply(first_check, generator, period))
    assert vertical_detector(first_boundary, second_boundary, period) == 0
    print({"fixed_height_four": "weight 6 < 12", "vertical_detector": "omega^2",
           "transverse_seam_coefficient": "u^2, zero in M"})


def height_two_axis_warning():
    period = 2
    _, _, first_check = pairs(period)
    transverse_period = {(0, 0): 1, (6, 0): 1}
    assert multiply(first_check, first_check, period) == transverse_period
    assert not fold_columns(transverse_period)
    assert len(first_check) == 3
    assert vertical_detector(first_check, {}, period) == 3
    print({"fixed_height_two": "weight 3, but NOT in the cylinder image",
           "period_coordinates": ["omega^2", "1"]})


def main():
    signal.alarm(55)
    resource.setrlimit(resource.RLIMIT_CPU, (50, 50))
    started = time.monotonic()
    local_module_check()
    plane_period_identity()
    height_four_witness()
    height_two_axis_warning()
    peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    peak_bytes = peak if sys.platform == "darwin" else 1024 * peak
    assert peak_bytes <= 768 * 1024 ** 2
    print({"status": "PASS", "wall_seconds": round(time.monotonic() - started, 6),
           "peak_rss_bytes": peak_bytes})
    print("Scope: folded image-class floor obstruction; original cylinder floor remains open.")


if __name__ == "__main__":
    main()
