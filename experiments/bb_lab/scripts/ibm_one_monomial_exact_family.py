"""Fixed identities for the one-monomial theorem, not a distance search."""

import resource
import signal
import sys
import time


def add(*polynomials):
    result = set()
    for polynomial in polynomials:
        result.symmetric_difference_update(polynomial)
    return result


def multiply(left, right):
    result = set()
    for left_column, left_row in left:
        for right_column, right_row in right:
            result.symmetric_difference_update(
                {(left_column + right_column, left_row + right_row)}
            )
    return result


def shear_evaluate(polynomial, slope):
    result = set()
    for column, row in polynomial:
        result.symmetric_difference_update({(column + slope * row, 0)})
    return result


def reduce_periods(polynomial, width, height):
    result = set()
    for column, row in polynomial:
        result.symmetric_difference_update({(column % width, row % height)})
    return result


def main():
    signal.alarm(55)
    resource.setrlimit(resource.RLIMIT_CPU, (50, 50))
    started = time.monotonic()
    identity = {(0, 0)}
    column = {(1, 0)}
    row = {(0, 1)}
    column_squared = {(2, 0)}
    row_squared = {(0, 2)}
    column_cubed = {(3, 0)}
    row_cubed = {(0, 3)}
    column_recurrence = add(identity, column, column_squared)
    row_recurrence = add(identity, row, row_squared)
    first_original = {(0, 0), (0, 1), (3, -1)}
    second_original = {(0, 0), (1, 0), (-1, -3)}
    second_modified = {(0, 0), (1, 0), (-1, 0)}
    first_cleared = add(column_cubed, row, row_squared)
    assert multiply(row, first_original) == first_cleared
    assert multiply(column, second_modified) == column_recurrence
    assert add(second_original, second_modified) == {(-1, -3), (-1, 0)}
    assert first_cleared == add(
        row_recurrence, multiply(add(identity, column), column_recurrence)
    )
    assert multiply(add(identity, column), column_recurrence) == add(
        identity, column_cubed
    )
    assert multiply(add(identity, row), row_recurrence) == add(
        identity, row_cubed
    )
    assert shear_evaluate(first_cleared, 1) == multiply(column, column_recurrence)
    assert shear_evaluate(first_cleared, 2) == multiply(
        column_squared, column_recurrence
    )
    assert add(column, column_squared, identity) == column_recurrence
    stripe = {(0, 0), (1, 0), (3, 0), (4, 0)}
    horizontal_period = {(0, 0), (6, 0)}
    assert multiply(column_recurrence, stripe) == horizontal_period
    assert not reduce_periods(multiply(column_recurrence, stripe), 6, 6)
    assert len(reduce_periods(stripe, 6, 6)) == 4
    assert len(first_cleared) + len(column_recurrence) == 6
    assert add(
        multiply(column, stripe), multiply(column_squared, stripe), stripe
    ) == horizontal_period
    elapsed = time.monotonic() - started
    peak_rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    peak_bytes = peak_rss if sys.platform == "darwin" else peak_rss * 1024
    assert peak_bytes < 768 * 1024 * 1024
    signal.alarm(0)
    print(
        f"PASS prescribed identities and length-six stripe; "
        f"{elapsed:.6f}s, {peak_bytes / 1024**2:.2f} MiB peak RSS. "
        "No distance enumeration."
    )


if __name__ == "__main__":
    main()
