"""Fixed symbolic boundary identities; no minimization or coset enumeration."""

import resource
import signal
import sys
import time


ALPHA = {(0, 0), (1, -1), (-4, 1)}
BETA = {(0, 0), (-1, 0), (-2, 3)}


def multiply(left, right, frame):
    result = set()
    for horizontal, vertical in left:
        for shift_x, shift_y in right:
            result.symmetric_difference_update({((horizontal + shift_x) % frame[0],
                                                (vertical + shift_y) % frame[1])})
    return result


def power(poly, exponent, frame):
    result = {(0, 0)}
    for index in range(exponent):
        result = multiply(result, poly, frame)
    return result


def shift(poly, exponent, frame):
    return multiply(poly, {exponent}, frame)


def syndrome(pair, frame):
    return multiply(BETA, pair[0], frame) ^ multiply(ALPHA, pair[1], frame)


def motif(radius):
    frame = (6 * radius, 2 * radius)
    return (shift(power(BETA, radius - 1, frame), (radius, 0), frame),
            shift(power(ALPHA, radius - 1, frame), (-radius, radius), frame))


def collective_rewrite(radius):
    frame = (6 * radius, 2 * radius)
    original = motif(radius)
    first_scaled = shift(BETA, (1, 0), frame)
    second_scaled = shift(ALPHA, (-1, 1), frame)
    total = first_scaled ^ second_scaled
    multiplier = set()
    for index in range(radius - 1):
        multiplier ^= multiply(power(first_scaled, radius - 2 - index, frame),
                               power(second_scaled, index, frame), frame)
    multiplier = shift(multiplier, (0, 1), frame)
    rewritten = (original[0] ^ multiply(ALPHA, multiplier, frame),
                 original[1] ^ multiply(BETA, multiplier, frame))
    expected = (shift(power(total, radius - 1, frame), (1, 0), frame),
                shift(power(total, radius - 1, frame), (-1, 1), frame))
    assert rewritten == expected
    assert not syndrome(rewritten, frame)
    print("collective", radius, "old", tuple(map(len, original)),
          "new", tuple(map(len, rewritten)), "boundary terms", len(multiplier))


def partial_polarization(radius):
    frame = (6 * radius, 2 * radius)
    half = radius // 2
    original = motif(radius)
    first_scaled = shift(BETA, (1, 0), frame)
    second_scaled = shift(ALPHA, (-1, 1), frame)
    total = first_scaled ^ second_scaled
    multiplier = shift(power(multiply(first_scaled, second_scaled, frame), half - 1, frame),
                       (0, 1), frame)
    rewritten = (original[0] ^ multiply(ALPHA, multiplier, frame),
                 original[1] ^ multiply(BETA, multiplier, frame))
    lower = (shift(power(first_scaled, half - 1, frame), (1, 0), frame),
             shift(power(second_scaled, half - 1, frame), (-1, 1), frame))
    expected = tuple(multiply(power(total, half, frame), block, frame) for block in lower)
    assert rewritten == expected
    assert not syndrome(rewritten, frame)
    assert sum(map(len, rewritten)) == 8 * 3 ** (radius.bit_length() - 2)
    print("partial polarization", radius, "blocks", tuple(map(len, rewritten)),
          "boundary terms", len(multiplier))


def verify_floor_constants():
    frame = (12, 4)
    alpha = multiply(ALPHA, {(0, 0)}, frame)
    beta = multiply(BETA, {(0, 0)}, frame)
    for polynomial in (alpha, beta):
        differences = {((left[0] - right[0]) % 12, (left[1] - right[1]) % 4)
                       for left in polynomial for right in polynomial if left != right}
        assert len(differences) == 6
        assert {horizontal % 2 for horizontal, vertical in polynomial} == {0, 1}
        assert {horizontal % 2 for horizontal, vertical in power(polynomial, 2, frame)} == {0}
    cross_differences = {((left[0] - right[0]) % 12, (left[1] - right[1]) % 4)
                         for left in alpha for right in beta}
    assert len(cross_differences) == 9
    alpha_patterns = []
    beta_patterns = []
    for polynomial, directions, patterns in (
        (alpha, ((1, -1), (-4, 1), (-5, 2)), alpha_patterns),
        (beta, ((-1, 0), (-2, 3), (-1, 3)), beta_patterns),
    ):
        for direction in directions:
            paired = polynomial ^ shift(polynomial, direction, frame)
            assert len(paired) == 4
            patterns.append(sorted(horizontal for horizontal, vertical in paired))
    assert alpha_patterns == [[0, 2, 8, 9], [0, 1, 4, 9], [0, 1, 3, 7]]
    assert beta_patterns == [[0, 9, 10, 10], [0, 8, 9, 11], [0, 9, 11, 11]]
    print("Q2 local proof constants PASS; no representative enumeration")
    print("alpha horizontal patterns", alpha_patterns)
    print("beta horizontal patterns", beta_patterns)
    print("Q2 weight-six motif", [sorted(block) for block in motif(2)])


def verify_inherited_identity():
    for radius in (1, 2):
        old_frame = (6 * radius, 2 * radius)
        new_frame = (12 * radius, 4 * radius)
        original = motif(radius)
        correction = multiply({(0, 0), (1, 1), (3, 0)}, {(0, 0)}, old_frame)
        corrected = (original[0] ^ multiply(ALPHA, correction, old_frame),
                     original[1] ^ multiply(BETA, correction, old_frame))
        squares = tuple({(2 * horizontal, 2 * vertical) for horizontal, vertical in block}
                        for block in corrected)
        lifted = multiply(BETA, squares[0], new_frame), multiply(ALPHA, squares[1], new_frame)
        assert sum(map(len, lifted)) == 3 * sum(map(len, corrected))
        multiplier = multiply(multiply(ALPHA, BETA, new_frame),
                              {(2 * horizontal, 2 * vertical) for horizontal, vertical in correction},
                              new_frame)
        raw = motif(2 * radius)
        assert lifted == (raw[0] ^ multiply(ALPHA, multiplier, new_frame),
                          raw[1] ^ multiply(BETA, multiplier, new_frame))
    print("fixed inherited-boundary identity checks PASS")


def main():
    signal.alarm(55)
    resource.setrlimit(resource.RLIMIT_CPU, (55, 55))
    started = time.monotonic()
    for radius in (2, 4, 8):
        collective_rewrite(radius)
        partial_polarization(radius)
    verify_floor_constants()
    verify_inherited_identity()
    peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    peak_bytes = peak if sys.platform == "darwin" else 1024 * peak
    assert peak_bytes <= 768 * 1024 ** 2
    print("PASS", round(time.monotonic() - started, 4), "seconds; peak MiB",
          round(peak_bytes / 1024 ** 2, 1))


if __name__ == "__main__":
    main()
