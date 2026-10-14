"""Fixed cylinder identities and metric counterexamples; no distance search."""

import resource
import signal
import sys
import time


def add(*polynomials):
    result = set()
    for polynomial in polynomials:
        result.symmetric_difference_update(polynomial)
    return result


def multiply(left, right, circumference):
    result = set()
    for first_column, first_row in left:
        for second_column, second_row in right:
            result.symmetric_difference_update({((first_column + second_column) % circumference,
                                                  first_row + second_row)})
    return result


def vector_add(left, right):
    return tuple(add(first, second) for first, second in zip(left, right))


def weight(vector):
    return sum(map(len, vector))


def overlap(left, right):
    return sum(len(first & second) for first, second in zip(left, right))


def row_mask(polynomial, row):
    return {(column, 0) for column, actual_row in polynomial if actual_row == row}


def shear_evaluate(polynomial, circumference):
    result = set()
    for column, row in polynomial:
        result.symmetric_difference_update({((column + row) % circumference, 0)})
    return result


def main():
    signal.alarm(55)
    resource.setrlimit(resource.RLIMIT_CPU, (50, 50))
    started = time.monotonic()
    circumference = 60
    identity = {(0, 0)}
    column = {(1, 0)}
    row = {(0, 1)}
    column_cubed = {(3, 0)}
    cyclotomic = {(0, 0), (1, 0), (2, 0)}
    first_factor = add(column, row)
    second_factor = add(first_factor, identity)
    first_check = {(3, 0), (0, 1), (0, 2)}
    second_check = {(3, 0), (0, 0), (0, 1), (1, 1), (2, 1)}

    def product(left, right):
        return multiply(left, right, circumference)

    def boundary(generator):
        return product(second_check, generator), product(first_check, generator)

    def syndrome(vector):
        return add(product(first_check, vector[0]), product(second_check, vector[1]))

    def label(vector):
        return add(product(column, shear_evaluate(vector[0], circumference)),
                   shear_evaluate(vector[1], circumference))

    assert second_check == product(cyclotomic, second_factor)
    assert add(product(first_check, cyclotomic), product(second_check, first_factor)) == product(
        column, product(cyclotomic, cyclotomic))
    assert shear_evaluate(first_check, circumference) == product(column, cyclotomic)
    assert shear_evaluate(second_check, circumference) == cyclotomic

    stripe = {(start + offset, 0) for start in range(0, circumference, 3) for offset in (0, 1)}
    target = 2 * circumference // 3
    assert len(stripe) == target == 40
    assert not product(cyclotomic, stripe)
    assert product(column_cubed, stripe) == stripe
    original = set(), product({(0, 2)}, stripe)
    normal = set(), product(add(identity, row), stripe)
    assert not syndrome(original) and not syndrome(normal)
    assert vector_add(original, normal) == boundary(stripe)
    assert weight(original) == target and weight(normal) == 2 * target
    assert label(original) == label(normal) == product({(2, 0)}, stripe)
    assert len(label(original)) == target

    top_mask = {(0, 0), (1, 0)}
    lower_mask = {(column_index, 0) for column_index in (0, 1, 2, 3, 4)}
    top_generator = product(row, top_mask)
    lower_boundary = boundary(lower_mask)
    top_boundary = boundary(top_generator)
    original_boundary = boundary(add(lower_mask, top_generator))
    assert weight(lower_boundary) == 26
    assert weight(top_boundary) == 12
    assert overlap(lower_boundary, top_boundary) == 7
    assert weight(original_boundary) == 24
    assert vector_add(original_boundary, top_boundary) == lower_boundary
    assert row_mask(original_boundary[1], 3) == top_mask
    assert max(actual_row for _, actual_row in original_boundary[1]) == 3
    assert row_mask(original_boundary[0], 2) == product(cyclotomic, top_mask)
    assert max(actual_row for _, actual_row in original_boundary[0]) == 2
    first_update = product(add(identity, column_cubed), top_mask)
    first_overlap = len(row_mask(original_boundary[0], 1) & first_update)
    second_overlap = len(row_mask(original_boundary[1], 2) & top_mask)
    third_overlap = len(row_mask(original_boundary[1], 1) & product(column_cubed, top_mask))
    intersection_sum = first_overlap + second_overlap + third_overlap
    assert intersection_sum == 1
    predicted_change = len(first_update) + len(top_mask) - len(product(cyclotomic, top_mask)) - 2 * intersection_sum
    assert predicted_change == weight(lower_boundary) - weight(original_boundary) == 2
    assert not label(original_boundary) and not syndrome(original_boundary)

    padding = product({(0, -5)}, add(identity, column_cubed))
    padding_boundary = boundary(padding)
    assert weight(padding_boundary) == 14
    assert overlap(padding_boundary, original_boundary) == 0
    assert overlap(padding_boundary, lower_boundary) == 0
    below_threshold = boundary(add(lower_mask, top_generator, padding))
    peeled = vector_add(below_threshold, top_boundary)
    assert weight(below_threshold) == 38 < target
    assert weight(peeled) == target
    assert peeled == boundary(add(lower_mask, padding))
    assert not syndrome(below_threshold) and not syndrome(peeled)
    assert not label(below_threshold) and not label(peeled)

    peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    peak_bytes = peak if sys.platform == "darwin" else 1024 * peak
    assert peak_bytes <= 768 * 1024 ** 2
    print({"status": "PASS", "fixed_circumference": circumference,
           "nontrivial_normalization_weights": [weight(original), weight(normal)],
           "local_peeling_weights": [weight(original_boundary), weight(lower_boundary)],
           "threshold_peeling_weights": [weight(below_threshold), weight(peeled)],
           "wall_seconds": round(time.monotonic() - started, 6), "peak_rss_bytes": peak_bytes,
           "scope": "fixed identities only; the full floor is proved separately by shear and seam"})


if __name__ == "__main__":
    main()
