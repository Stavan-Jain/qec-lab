"""Fixed identities and decoder conventions for the canonical recovery obstruction."""

import json
import resource
import signal
import sys
import time


def shift(mask, amount, period):
    amount %= period
    full = (1 << period) - 1
    return ((mask << amount) | (mask >> (period - amount))) & full


def product(first, second, period):
    result = 0
    while second:
        term = second & -second
        result ^= shift(first, term.bit_length() - 1, period)
        second ^= term
    return result


def times_u(mask, period):
    return shift(mask, 1, period) ^ shift(mask, 2, period)


def inverse_phi(period):
    assert period % 3 != 0
    residues = {0, 2} if period % 3 == 1 else {1, 2}
    inverse = sum(1 << row for row in range(period) if row % 3 in residues)
    assert inverse ^ times_u(inverse, period) == 1
    return inverse


def inside_syndrome(inside, period):
    first, second_zero, second_one = inside
    return (times_u(first, period) ^ second_zero,
            second_one ^ shift(second_zero, 3, period),
            shift(second_zero ^ second_one, 3, period),
            first ^ shift(second_one, 3, period))


def recover(syndrome, period):
    zero, one, two, three = syndrome
    second_zero = product(inverse_phi(period), zero ^ times_u(one ^ three ^ two ^ shift(two, -3, period), period), period)
    second_one = one ^ shift(second_zero, 3, period)
    first = three ^ shift(one, 3, period) ^ shift(second_zero, 6, period)
    return first, second_zero, second_one


def impulse_syndrome(block, column, period):
    result = [0, 0, 0, 0]
    terms = ((3, 0), (0, 1), (0, 2)) if block == 0 else ((0, 0), (1, 3), (2, 3))
    for horizontal, vertical in terms:
        target = column + horizontal
        if 0 <= target <= 3:
            result[target] ^= shift(1, vertical, period)
    return tuple(result)


def check_period(period):
    for block in range(3):
        source = tuple(1 if index == block else 0 for index in range(3))
        assert recover(inside_syndrome(source, period), period) == source
    outside_columns = [(0, column) for column in (-3, -2, -1, 1, 2, 3)]
    outside_columns += [(1, column) for column in (-2, -1, 2, 3)]
    responses = {position: recover(impulse_syndrome(*position, period), period) for position in outside_columns}
    loads = {f"{block}:{column}": sum(mask.bit_count() for mask in response)
             for (block, column), response in responses.items()}
    expected_load = 2 * period + (1 if period % 3 == 1 else -1)
    assert loads["0:-3"] == expected_load
    for generator_column in range(-6, 7):
        first_boundary = {generator_column: 1, generator_column + 1: shift(1, 3, period),
                          generator_column + 2: shift(1, 3, period)}
        second_boundary = {generator_column: times_u(1, period), generator_column + 3: 1}
        true_inside = (first_boundary.get(0, 0), second_boundary.get(0, 0), second_boundary.get(1, 0))
        recovered = [0, 0, 0]
        for (block, column), response in responses.items():
            source = (first_boundary, second_boundary)[block].get(column, 0)
            for output in range(3):
                recovered[output] ^= product(response[output], source, period)
        assert tuple(recovered) == true_inside
    return {"height": period, "outside_block_columns": len(outside_columns),
            "loads": loads, "maximum_load": max(loads.values()),
            "sharp_load3_requires_at_least_block_columns": (period + 2) // 3 - 1}


def main():
    signal.alarm(55)
    resource.setrlimit(resource.RLIMIT_CPU, (50, 50))
    started = time.monotonic()
    mask = 3
    assert inside_syndrome((mask, mask, mask), 3) == (0, 0, 0, 0)
    reports = [check_period(period) for period in (1, 2, 4, 8, 16)]
    peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    peak_bytes = peak if sys.platform == "darwin" else peak * 1024
    assert peak_bytes < 768 * 1024 ** 2
    print(json.dumps({"status": "PASS", "period_three_cycle_weight": 6,
                      "fixed_replays": reports, "elapsed_seconds": time.monotonic() - started,
                      "peak_rss_bytes": peak_bytes}, indent=2))


if __name__ == "__main__":
    main()
