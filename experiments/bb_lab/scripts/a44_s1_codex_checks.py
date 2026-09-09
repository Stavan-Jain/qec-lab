"""Bounded sparse identities for A44-S1; no solver, census, or output files."""

import json
import resource
import signal
import sys
import time
from pathlib import Path


LAB = Path(__file__).resolve().parent.parent
ALPHA = {(0, 0), (0, -1), (-3, 1)}
BETA = {(0, 0), (-1, 0), (1, 3)}


def reduce_poly(poly, frame=(None, None)):
    result = set()
    for horizontal, vertical in poly:
        exponent = (
            horizontal if frame[0] is None else horizontal % frame[0],
            vertical if frame[1] is None else vertical % frame[1],
        )
        result.symmetric_difference_update({exponent})
    return result


def multiply(left, right, frame=(None, None)):
    return reduce_poly(
        ((horizontal + shift_x, vertical + shift_y)
         for horizontal, vertical in left for shift_x, shift_y in right),
        frame,
    )


def shifted(poly, shift, frame=(None, None)):
    return multiply(poly, {shift}, frame)


def syndrome(pair, frame=(None, None)):
    return multiply(BETA, pair[0], frame) ^ multiply(ALPHA, pair[1], frame)


def split_cells(cells):
    return tuple({(horizontal, vertical) for block, horizontal, vertical in cells
                  if block == wanted} for wanted in (0, 1))


def reflect(poly):
    return {(-horizontal, -vertical) for horizontal, vertical in poly}


def verify_banked_upper():
    bands = json.loads((LAB / "data/a40/s2_ub_bands.json").read_text())
    dual_data = json.loads((LAB / "data/a40/s2_ub_dual3.json").read_text())
    carrier = split_cells(bands["L12"])
    dual = split_cells(dual_data["dual"]["triples"])
    assert not syndrome(carrier, (None, 6))
    assert not (multiply(reflect(ALPHA), dual[0], (6, None))
                ^ multiply(reflect(BETA), dual[1], (6, None)))
    print("L12 plane residual:", sorted(syndrome(carrier)))
    print("dual plane residual:", sorted(multiply(reflect(ALPHA), dual[0])
                                         ^ multiply(reflect(BETA), dual[1])))
    for radius, offset in ((1, 0), (1, 1), (2, 0), (2, 1), (3, 0), (3, 1)):
        frame = (6 * (radius + offset), 6 * radius)
        witness = tuple(reduce_poly(
            ((horizontal, vertical + 6 * index)
             for horizontal, vertical in block for index in range(radius)), frame)
            for block in carrier)
        certificate = tuple(reduce_poly(
            ((horizontal + 6 * index, vertical)
             for horizontal, vertical in block for index in range(radius + offset)), frame)
            for block in dual)
        pairing = sum(len(witness[index] & certificate[index]) for index in (0, 1)) % 2
        assert not syndrome(witness, frame)
        assert not (multiply(reflect(ALPHA), certificate[0], frame)
                    ^ multiply(reflect(BETA), certificate[1], frame))
        assert pairing == 1
        assert sum(map(len, witness)) == 12 * radius
        print("stack", frame, "weight", sum(map(len, witness)), "pairing", pairing)


def frobenius_lift(pair, frame):
    squares = tuple({(2 * horizontal, 2 * vertical) for horizontal, vertical in block}
                    for block in pair)
    return multiply(BETA, squares[0], frame), multiply(ALPHA, squares[1], frame)


def improve_locally(pair, frame):
    result = tuple(set(block) for block in pair)
    changes = []
    for iteration in range(20):
        candidates = set()
        for block, generator in zip(result, (ALPHA, BETA)):
            for horizontal, vertical in block:
                for shift_x, shift_y in generator:
                    candidates.add(((horizontal - shift_x) % frame[0],
                                    (vertical - shift_y) % frame[1]))
        best = None
        for shift in sorted(candidates):
            generator = tuple(shifted(poly, shift, frame) for poly in (ALPHA, BETA))
            overlap = sum(len(block & check) for block, check in zip(result, generator))
            if overlap > 3 and (best is None or overlap > best[0]):
                best = overlap, shift, generator
        if best is None:
            break
        result = tuple(block ^ check for block, check in zip(result, best[2]))
        changes.append(best[1])
    assert not syndrome(result, frame)
    return result, changes


def verify_lift():
    seed = ({(1, 0), (3, 2), (5, 4)}, {(0, 1), (2, 3), (4, 5)})
    assert not syndrome(seed, (6, 6))
    for radius in (2, 4):
        frame = (6 * radius, 6 * radius)
        lifted = frobenius_lift(seed, frame)
        assert not syndrome(lifted, frame)
        reduced, changes = improve_locally(lifted, frame)
        print("Frobenius", radius, "raw", sum(map(len, lifted)),
              "reduced", sum(map(len, reduced)), "shifts", changes,
              "nontrivial", not is_stabilizer(reduced, frame))
        seed = reduced


def vector_bits(pair, frame):
    return sum(1 << (block * frame[0] * frame[1] + horizontal * frame[1] + vertical)
               for block, poly in enumerate(pair) for horizontal, vertical in poly)


def is_stabilizer(pair, frame):
    basis = {}
    for horizontal in range(frame[0]):
        for vertical in range(frame[1]):
            row = vector_bits(tuple(shifted(poly, (horizontal, vertical), frame)
                                    for poly in (ALPHA, BETA)), frame)
            while row:
                pivot = row.bit_length() - 1
                if pivot not in basis:
                    basis[pivot] = row
                    break
                row ^= basis[pivot]
    row = vector_bits(pair, frame)
    while row:
        pivot = row.bit_length() - 1
        if pivot not in basis:
            return False
        row ^= basis[pivot]
    return True


def verify_obstructions():
    frame = (18, 18)
    initial = {(horizontal + shift, 0) for horizontal in (0, 2) for shift in (0, 6, 12)}
    row = initial
    band = set()
    for index in range(6):
        band ^= shifted(row, (0, 1 - 3 * index), frame)
        row = multiply(row, {(-1, 0), (-2, 0)}, frame)
    assert row == initial
    assert len(band) == 54
    assert not syndrome((band, set()), frame)
    dual_data = json.loads((LAB / "data/a40/s2_ub_dual3.json").read_text())
    dual_first = split_cells(dual_data["dual"]["triples"])[0]
    dual_first = reduce_poly(((horizontal + 6 * index, vertical)
                             for horizontal, vertical in dual_first for index in range(3)), frame)
    assert len(band & dual_first) % 2 == 1
    print("r=3 pure-band carrier: weight 54, certified odd dual pairing")
    for radius in (1, 2, 3):
        frame = (6 * radius, 6 * radius)
        second = {(1 + 2 * index, 2 * index) for index in range(3 * radius)}
        second = reduce_poly(second, frame)
        pair = shifted(second, (-3, 1), frame), second
        residual = syndrome(pair, frame)
        expected = multiply(multiply({(0, 0), (0, -1)}, {(0, 0), (-6, 0)}, frame),
                            second, frame)
        assert residual == expected
        assert len(residual) == (0 if radius == 1 else 12 * radius)
        print("single helix", radius, "syndrome weight", len(residual))


def verify_quotient_recurrence():
    frame = (24, 8)
    alpha = {(0, 0), (1, -1), (-4, 1)}
    beta = {(0, 0), (-1, 0), (-2, 3)}
    first = multiply(multiply(beta, beta, frame), beta, frame)
    first = shifted(first, (4, 0), frame)
    second = multiply(multiply(alpha, alpha, frame), alpha, frame)
    second = shifted(second, (-4, 4), frame)
    pair = first, second
    assert len(first) == len(second) == 9
    assert not (multiply(beta, first, frame) ^ multiply(alpha, second, frame))
    checks = []
    for horizontal in range(frame[0]):
        for vertical in range(frame[1]):
            checks.append(tuple(shifted(poly, (horizontal, vertical), frame)
                                for poly in (alpha, beta)))
    overlaps = [sum(len(block & check) for block, check in zip(pair, generator))
                for generator in checks]
    pair_overlaps = [sum(len(left & right) for left, right in zip(checks[0], generator))
                     for generator in checks[1:]]
    assert max(overlaps) == 2
    assert max(pair_overlaps) == 1
    print("r=4 odd-degree quotient: 384 qubits, raw motif weight 18, target 14")
    print("motif f(z,t):", sorted(first))
    print("motif g(z,t):", sorted(second))
    print("max motif/check overlap 2; max distinct check/check overlap 1")


def verify_formal_pairing():
    def field_product(left, right):
        result = 0
        while right:
            if right & 1:
                result ^= left
            right >>= 1
            left <<= 1
            if left & 4:
                left ^= 7
        return result

    def series_sum(left, right):
        return [first ^ second for first, second in zip(left, right)]

    def series_product(left, right):
        result = [0] * len(left)
        for left_degree, left_coefficient in enumerate(left):
            for right_degree, right_coefficient in enumerate(right[:len(left) - left_degree]):
                result[left_degree + right_degree] ^= field_product(left_coefficient, right_coefficient)
        return result

    def series_power(series, exponent):
        result = [1] + [0] * (len(series) - 1)
        for index in range(exponent):
            result = series_product(result, series)
        return result

    def series_inverse(series):
        inverse_constant = field_product(series[0], series[0])
        result = [inverse_constant] + [0] * (len(series) - 1)
        for degree in range(1, len(series)):
            coefficient = 0
            for index in range(1, degree + 1):
                coefficient ^= field_product(series[index], result[degree - index])
            result[degree] = field_product(inverse_constant, coefficient)
        return result

    for radius in (1, 2, 4, 8):
        order = 2 * radius
        one = [1] + [0] * (order - 1)
        diagonal = [1, 1] + [0] * (order - 2)
        horizontal = [2] + [0] * (order - 1)
        for iteration in range(5):
            horizontal = series_sum(diagonal, series_product(series_power(diagonal, 2),
                                    series_power(series_inverse(horizontal), 4)))
        inverse_horizontal = series_inverse(horizontal)
        alpha = series_sum(one, series_sum(
            series_product(horizontal, series_inverse(diagonal)),
            series_product(diagonal, series_power(inverse_horizontal, 4))))
        beta = series_sum(one, series_sum(inverse_horizontal,
            series_product(series_power(diagonal, 3), series_power(inverse_horizontal, 2))))
        assert not any(alpha)
        assert beta[:min(3, order)] == ([0, 0] if order == 2 else [0, 0, 1])
        value = series_product(series_power(horizontal, radius), series_power(beta, radius - 1))
        scale = (1, 2, 3)[(1 - radius) % 3]
        coefficient = field_product(scale, value[2 * radius - 2])
        assert coefficient ^ field_product(coefficient, coefficient) == 1
        assert series_power(horizontal, 6 * radius) == one
        assert series_power(diagonal, 2 * radius) == one
    print("formal F4 pairing: PASS at r=1,2,4,8; universal argument is in the note")


def main():
    signal.alarm(55)
    resource.setrlimit(resource.RLIMIT_CPU, (55, 55))
    started = time.monotonic()
    verify_banked_upper()
    verify_lift()
    verify_obstructions()
    verify_quotient_recurrence()
    verify_formal_pairing()
    peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    peak_bytes = peak if sys.platform == "darwin" else 1024 * peak
    assert peak_bytes < 1024 ** 3
    print("PASS", round(time.monotonic() - started, 3), "seconds; peak MiB",
          round(peak_bytes / 1024 ** 2, 1))


if __name__ == "__main__":
    main()
