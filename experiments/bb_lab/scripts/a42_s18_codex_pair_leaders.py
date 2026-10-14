"""Finite identities for pair leaders; the unrestricted route remains open."""

import resource
import signal
import sys
import time
from collections import Counter, defaultdict


def add(*polynomials):
    result = set()
    for polynomial in polynomials:
        result.symmetric_difference_update(polynomial)
    return result


def multiply(left, right, height=None, width=None):
    result = set()
    for first_column, first_row in left:
        for second_column, second_row in right:
            column = first_column + second_column
            row = first_row + second_row
            if height is not None:
                row %= height
            if width is not None:
                column %= width
            result.symmetric_difference_update({(column, row)})
    return result


def boundary(generator, height=None, width=None):
    first_check = {(3, 0), (0, 1), (0, 2)}
    second_check = {(0, 0), (1, 3), (2, 3)}
    return (multiply(second_check, generator, height, width),
            multiply(first_check, generator, height, width))


def syndrome(chain, height=None, width=None):
    first_check = {(3, 0), (0, 1), (0, 2)}
    second_check = {(0, 0), (1, 3), (2, 3)}
    return add(multiply(first_check, chain[0], height, width),
               multiply(second_check, chain[1], height, width))


def formal_columns(multiplier):
    generator = {(-3, 0, "a"), (-2, 0, "b"), (-1, 0, "c")}
    output = set()
    for column, row, label in generator:
        for column_shift, row_shift in multiplier:
            output.symmetric_difference_update({(column + column_shift, row + row_shift, label)})
    columns = defaultdict(set)
    for column, row, label in output:
        columns[column].add((label, row))
    return dict(columns)


def verify_strip_identities():
    assert formal_columns({(0, 0), (1, 3), (2, 3)}) == {
        -3: {("a", 0)},
        -2: {("b", 0), ("a", 3)},
        -1: {("c", 0), ("a", 3), ("b", 3)},
        0: {("b", 3), ("c", 3)},
        1: {("c", 3)},
    }
    assert formal_columns({(3, 0), (0, 1), (0, 2)}) == {
        -3: {("a", 1), ("a", 2)},
        -2: {("b", 1), ("b", 2)},
        -1: {("c", 1), ("c", 2)},
        0: {("a", 0)},
        1: {("b", 0)},
        2: {("c", 0)},
    }
    total = Counter(alpha=2, beta=1, gamma=2, R1=1, R2=1, Z=1, U=1)
    intersections = {
        "AB": Counter(alpha=1, Z=1),
        "AC": Counter(beta=1, Z=1),
        "BC": Counter(alpha=1, beta=1),
    }
    decompositions = {
        "AB": [Counter(beta=1, gamma=1, Z=-1),
               Counter(gamma=1, R1=1, R2=1, U=1)],
        "BC": [Counter(Z=1, gamma=1, beta=-1),
               Counter(gamma=1, R1=1, R2=1, U=1)],
        "AC": [Counter(R1=1, alpha=1, beta=-1),
               Counter(R2=1, gamma=1, V=-1),
               Counter(V=1, alpha=1, beta=-1),
               Counter(beta=1, gamma=1, Z=-1), Counter(U=1)],
    }
    for pair, intersection in intersections.items():
        energy = total.copy()
        energy.subtract({name: 2 * coefficient for name, coefficient in intersection.items()})
        certificate = Counter()
        for term in decompositions[pair]:
            certificate.update(term)
        assert energy == certificate


def verify_prescribed_fold_identity():
    height = 8
    norm = {(0, row) for row in range(height)}
    inverse = {(0, row) for row in (1, 2, 4, 5, 7)}
    phi = {(0, 0), (0, 1), (0, 2)}
    chi = {(0, 0), (1, 0), (2, 0)}
    adjacent = {(1, 0), (2, 0)}
    shifted_norm = multiply({(1, 0)}, norm, height)
    assert multiply(phi, inverse, height) == {(0, 0)}
    representative = (add(inverse, multiply(add(norm, inverse), adjacent, height)), chi)
    pair_ab = norm, norm
    pair_ac = norm, shifted_norm
    pair_bc = set(), add(norm, shifted_norm)
    assert sum(map(len, representative)) == 14 < 2 * height
    expected_syndrome = multiply(adjacent, norm, height, 3)
    assert syndrome(representative, height, 3) == expected_syndrome
    for pair, generator in (
        (pair_ab, add(norm, multiply(chi, inverse, height, 3))),
        (pair_ac, add(shifted_norm, multiply(chi, inverse, height, 3))),
    ):
        assert syndrome(pair, height, 3) == expected_syndrome
        assert boundary(generator, height, 3) == (
            add(pair[0], representative[0]), add(pair[1], representative[1]))
    assert boundary(add(norm, shifted_norm), height, 3) == pair_bc
    common = add(multiply(chi, inverse, height),
                 multiply(adjacent, norm, height), {(1, 3)})
    for pair, residual_factor in ((pair_ac, common), (pair_ab, add(common, norm))):
        actual_residual = add(syndrome(representative, height), syndrome(pair, height))
        expected_residual = multiply({(0, 0), (3, 0)}, residual_factor, height)
        assert actual_residual == expected_residual
        assert actual_residual
        assert {row for column, row in actual_residual if column == 5} == {0, 3, 6}
    dressing = multiply({(-3, 0), (-2, 0)}, norm, height)
    assert boundary(dressing, height) == (
        add(norm, multiply({(-3, 0)}, norm, height)), add(norm, shifted_norm))


def main():
    signal.alarm(55)
    resource.setrlimit(resource.RLIMIT_CPU, (50, 50))
    started = time.monotonic()
    verify_strip_identities()
    verify_prescribed_fold_identity()
    print("PASS: formal three-column expansions and three triangle certificates")
    print("PASS: prescribed q=8 folded representatives, exact boundaries, nonzero original residuals")
    print("Scope: unrestricted original pair leaders remain unproved; no original counterexample")
    peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    peak_bytes = peak if sys.platform == "darwin" else peak * 1024
    assert peak_bytes <= 768 * 1024 ** 2
    print("Elapsed", round(time.monotonic() - started, 6), "seconds; peak MiB",
          round(peak_bytes / 1024 ** 2, 1))


if __name__ == "__main__":
    main()
