"""Fixed shear and logical-vector identities; no distance or support search."""

import resource
import signal
import sys
import time
from collections import Counter

sys.dont_write_bytecode = True

from ibm_one_check_modification import add, evaluate_dual, multiply


def shear(polynomial, block_shift=0):
    return {(column + row + block_shift, row) for column, row in polynomial}


def inverse_shear(polynomial, block_shift=0):
    return {(column - row - block_shift, row) for column, row in polynomial}


def vertical_augmentation(polynomial):
    return add(*({(column, 0)} for column, row in polynomial))


def reduce_columns(polynomial, width):
    return add(*({(column % width, row)} for column, row in polynomial))


def fixed_shear_checks():
    first_check = {(3, 0), (0, 1), (0, 2)}
    second_check = {(3, 0), (0, 0), (0, 1), (1, 1), (2, 1)}
    cyclotomic = {(0, 0), (1, 0), (2, 0)}
    transformed_first = shear(first_check, -1)
    transformed_second = shear(second_check)
    assert transformed_first == {(2, 0), (0, 1), (1, 2)}
    assert transformed_second == {(0, 0), (3, 0), (1, 1), (2, 1), (3, 1)}
    assert vertical_augmentation(transformed_first) == cyclotomic
    assert vertical_augmentation(transformed_second) == cyclotomic
    assert {column for column, row in transformed_first} == {0, 1, 2}
    assert {column for column, row in transformed_second} == {0, 1, 2, 3}
    assert (0, 1) in transformed_first and (0, 0) in transformed_second
    first = {(-3, 2), (1, -2), (0, 0)}
    second = {(2, -1), (-1, 3)}
    transformed_vector = shear(first, 1), shear(second)
    assert inverse_shear(transformed_vector[0], 1) == first
    assert inverse_shear(transformed_vector[1]) == second
    assert tuple(map(len, transformed_vector)) == (len(first), len(second))
    original_syndrome = add(multiply(first_check, first), multiply(second_check, second))
    transformed_syndrome = add(multiply(transformed_first, transformed_vector[0]),
                               multiply(transformed_second, transformed_vector[1]))
    assert transformed_syndrome == shear(original_syndrome)
    projected_vector = add(*(vertical_augmentation(block) for block in transformed_vector))
    assert vertical_augmentation(transformed_syndrome) == multiply(cyclotomic, projected_vector)
    generator = {(-2, 1), (0, -1)}
    original_boundary = multiply(second_check, generator), multiply(first_check, generator)
    transformed_boundary = shear(original_boundary[0], 1), shear(original_boundary[1])
    transformed_generator = shear(generator, 1)
    assert transformed_boundary == (multiply(transformed_second, transformed_generator),
                                     multiply(transformed_first, transformed_generator))
    assert add(*(vertical_augmentation(block) for block in transformed_boundary)) == set()
    column_counts = Counter(column for block in transformed_boundary for column, row in block)
    assert all(count % 2 == 0 for count in column_counts.values())
    assert 2 * len(column_counts) <= sum(map(len, transformed_boundary))
    return first_check, second_check, cyclotomic


def fixed_classical_identities(cyclotomic):
    period_three = {(0, 0), (3, 0)}
    prescribed_words = (
        ({(0, 0), (1, 0)}, {(0, 0)}),
        ({(1, 0), (2, 0)}, {(1, 0)}),
        ({(0, 0), (2, 0)}, {(0, 0), (1, 0)}),
    )
    for word, multiplier in prescribed_words:
        assert multiply(cyclotomic, word) == multiply(period_three, multiplier)
        assert len(word) == 2


def fixed_logical_vectors(first_check, second_check, cyclotomic):
    period_six = {(0, 0), (6, 0)}
    first_factor = {(0, 1), (1, 0)}
    second_factor = add(first_factor, {(0, 0)})
    stripe = {(0, 0), (1, 0), (3, 0), (4, 0)}
    assert multiply(cyclotomic, stripe) == period_six
    stripe_syndrome = multiply(second_check, stripe)
    assert stripe_syndrome == multiply(second_factor, period_six)
    assert not reduce_columns(stripe_syndrome, 6)
    assert len(reduce_columns(stripe, 6)) == 4
    assert evaluate_dual(second_factor, (2, 0), (2, 0)) == (1, 0)
    quotient = {(0, 0), (2, 0)}
    assert multiply(multiply(cyclotomic, cyclotomic), quotient) == period_six
    tangent_first = multiply(cyclotomic, quotient)
    tangent_second = multiply({(0, 1), (3, 0), (2, 0)}, quotient)
    tangent_coefficient = multiply({(1, 0)}, first_factor)
    tangent_syndrome = add(multiply(first_check, tangent_first),
                           multiply(second_check, tangent_second))
    assert tangent_syndrome == multiply(period_six, tangent_coefficient)
    assert not reduce_columns(tangent_syndrome, 6)
    assert (len(reduce_columns(tangent_first, 6)),
            len(reduce_columns(tangent_second, 6))) == (4, 6)
    assert evaluate_dual(tangent_coefficient, (2, 0), (2, 0)) == (0, 0)
    assert evaluate_dual(tangent_coefficient, (2, 1), (3, 3)) == (2, 2)
    assert add(vertical_augmentation(shear(tangent_first, 1)),
               vertical_augmentation(shear(tangent_second))) == set()
    print({"fixed_width": 6, "stripe_weight": 4, "pure_double_vector_weight": 10,
           "pure_double_period_coefficient": "x(y+x), a unit on the doubled component"})


def main():
    signal.alarm(55)
    resource.setrlimit(resource.RLIMIT_CPU, (50, 50))
    started = time.monotonic()
    first_check, second_check, cyclotomic = fixed_shear_checks()
    fixed_classical_identities(cyclotomic)
    fixed_logical_vectors(first_check, second_check, cyclotomic)
    peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    peak_bytes = peak if sys.platform == "darwin" else 1024 * peak
    assert peak_bytes <= 768 * 1024 ** 2
    print({"status": "PASS", "wall_seconds": round(time.monotonic() - started, 6),
           "peak_rss_bytes": peak_bytes})
    print("Fixed identities only; the note proves the all-class cylinder distance by an empty seam.")


if __name__ == "__main__":
    main()
