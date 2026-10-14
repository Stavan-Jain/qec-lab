"""Fixed Laurent identities for the S3 seam obstruction; no search."""

import resource
import signal
import sys
import time


ALPHA = {(0, 0), (1, -1), (-4, 1)}
BETA = {(0, 0), (-1, 0), (-2, 3)}
SEAM = ({(-4, 0), (-2, 0), (-1, 1), (-2, 1), (-3, 1)}, {(-2, 0)})
CHARGE = {
    (-3, 0), (-4, 0), (-5, 0),
    (-3, 2), (-4, 2), (-5, 2),
    (-4, 1), (-6, 1), (-1, -1),
}


def multiply(left, right):
    result = set()
    for horizontal, vertical in left:
        for shift_x, shift_y in right:
            result.symmetric_difference_update({(horizontal + shift_x, vertical + shift_y)})
    return result


def translate(polynomial, position):
    return multiply(polynomial, {position})


def reduce_poly(polynomial, width, height):
    result = set()
    for horizontal, vertical in polynomial:
        result.symmetric_difference_update({(horizontal % width, vertical % height)})
    return result


def syndrome(pair):
    return multiply(BETA, pair[0]) ^ multiply(ALPHA, pair[1])


def evaluate_ordinary_point(polynomial):
    omega_powers = (1, 2, 3)
    result = 0
    for horizontal, vertical in polynomial:
        result ^= omega_powers[horizontal % 3]
    return result


def check_laurent_constants():
    divisor = {(0, 0), (0, 2)}
    assert syndrome(SEAM) == multiply(divisor, CHARGE)
    assert evaluate_ordinary_point(ALPHA) == 0
    assert evaluate_ordinary_point(BETA) == 0
    assert evaluate_ordinary_point(CHARGE) == 1
    quadratic = {(2, 0), (1, 0), (0, 0)}
    cubic = {(3, 0), (1, 0), (0, 0)}
    assert multiply(quadratic, cubic) == {(5, 0), (4, 0), (0, 0)}
    coefficient_a = {(1, 0), (-4, 0)}
    coefficient_b = {(0, 0), (-1, 0)}
    coefficient_c = {(-2, 0)}
    cubic_quadratic = multiply(multiply(quadratic, quadratic), quadratic)
    assert coefficient_c ^ multiply(coefficient_a, coefficient_b) == translate(cubic_quadratic, (-5, 0))
    folded_first = translate(cubic, (-4, 0))
    first_choice = translate(quadratic, (-3, 0))
    second_choice = first_choice ^ folded_first
    assert first_choice == {(-1, 0), (-2, 0), (-3, 0)}
    assert second_choice == {(-2, 0), (-4, 0)}
    assert (-2, 0) not in folded_first
    assert not first_choice <= folded_first
    assert not second_choice <= folded_first
    differences = [left[0] - right[0] for left in cubic for right in cubic if left != right]
    assert sorted(differences) == [-3, -2, -1, 1, 2, 3]
    print("Plane transporter, ordinary charge, and cylinder lower-bound constants PASS")


def check_one_original_member():
    radius = 4
    period = 6 * radius
    repetitions = {(0, 2 * index) for index in range(3 * radius)}
    lifted = tuple(reduce_poly(multiply(repetitions, block), period, period) for block in SEAM)
    assert tuple(map(len, lifted)) == (15 * radius, 3 * radius)
    assert not reduce_poly(syndrome(lifted), period, period)
    assert all(reduce_poly(translate(block, (0, 2)), period, period) == block for block in lifted)
    physical = tuple({((horizontal + vertical) % period, vertical) for horizontal, vertical in block}
                     for block in lifted)
    original_alpha = {(0, 0), (0, -1), (-3, 1)}
    original_beta = {(0, 0), (-1, 0), (1, 3)}
    physical_syndrome = multiply(original_beta, physical[0]) ^ multiply(original_alpha, physical[1])
    assert not reduce_poly(physical_syndrome, period, period)
    assert sum(map(len, physical)) == 72
    assert all(reduce_poly(translate(block, (2, 2)), period, period) == block for block in physical)
    print("Single original (24,24) seam: weight 72, cycle, diagonal-period-two invariance PASS")


def main():
    signal.alarm(55)
    resource.setrlimit(resource.RLIMIT_CPU, (50, 50))
    started = time.monotonic()
    check_laurent_constants()
    check_one_original_member()
    peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    peak_bytes = peak if sys.platform == "darwin" else 1024 * peak
    assert peak_bytes <= 768 * 1024 ** 2
    print("PASS", round(time.monotonic() - started, 4), "seconds; peak MiB",
          round(peak_bytes / 1024 ** 2, 1))


if __name__ == "__main__":
    main()
