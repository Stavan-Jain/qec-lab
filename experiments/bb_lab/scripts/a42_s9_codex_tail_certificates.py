"""Exact original-cylinder tail certificates; no optimization or width search."""

from __future__ import annotations

import json
import resource
import signal
import time


def shift_rows(mask, amount, height):
    amount %= height
    return ((mask << amount) | (mask >> (height - amount))) & ((1 << height) - 1)


def derivative(mask, height):
    return shift_rows(mask, 1, height) ^ shift_rows(mask, 2, height)


def toggle(columns, position, mask):
    columns[position] = columns.get(position, 0) ^ mask
    if not columns[position]:
        del columns[position]


def boundary(generator, height):
    first = {}
    second = {}
    for position, mask in generator.items():
        shifted = shift_rows(mask, 3, height)
        toggle(first, position, mask)
        toggle(first, position + 1, shifted)
        toggle(first, position + 2, shifted)
        toggle(second, position, derivative(mask, height))
        toggle(second, position + 3, mask)
    return first, second


def statistics(generator, height):
    first, second = boundary(generator, height)
    inside = first.get(0, 0).bit_count()
    inside += second.get(0, 0).bit_count() + second.get(1, 0).bit_count()
    total = sum(mask.bit_count() for mask in first.values())
    total += sum(mask.bit_count() for mask in second.values())
    return {"inside": inside, "outside": total - inside, "energy": 3 * total - 4 * inside}


def left_deletion_charge(generator, position, height):
    assert position == min(generator) and position <= -4
    current = generator[position]
    next_one = generator.get(position + 1, 0)
    next_two = generator.get(position + 2, 0)
    next_three = generator.get(position + 3, 0)
    shifted = shift_rows(current, 3, height)
    intersections = (shifted & next_one).bit_count()
    intersections += (shifted & (next_two ^ shift_rows(next_one, 3, height))).bit_count()
    intersections += (current & derivative(next_three, height)).bit_count()
    return 3 * (2 * intersections - 4 * current.bit_count() - derivative(current, height).bit_count())


def right_deletion_charge(generator, position, height):
    assert position == max(generator) and position >= 3
    current = generator[position]
    previous_one = generator.get(position - 1, 0)
    previous_two = generator.get(position - 2, 0)
    previous_three = generator.get(position - 3, 0)
    intersections = (current & shift_rows(previous_one ^ previous_two, 3, height)).bit_count()
    intersections += (current & previous_one).bit_count()
    intersections += (derivative(current, height) & previous_three).bit_count()
    return 3 * (2 * intersections - 4 * current.bit_count() - derivative(current, height).bit_count())


def main():
    signal.alarm(55)
    resource.setrlimit(resource.RLIMIT_CPU, (50, 50))
    started = time.monotonic()
    height = 2
    constant = 3
    alternating = 1
    left = {-4: constant, -3: constant, -1: alternating, 1: alternating}
    left_before = statistics(left, height)
    assert left_before == {"inside": 5, "outside": 9, "energy": 22}
    clipped = {position: mask for position, mask in left.items() if position >= -3}
    assert statistics(clipped, height) == {"inside": 5, "outside": 13, "energy": 34}
    assert left_deletion_charge(left, -4, height) == 12
    moved_energies = {}
    for destination, expected in zip(range(-3, 3), (24, 30, 28, 34, 28, 46), strict=True):
        moved = dict(clipped)
        toggle(moved, destination, constant)
        moved_energies[destination] = statistics(moved, height)["energy"]
        assert moved_energies[destination] == expected
        assert moved_energies[destination] > left_before["energy"]
    translated = {position + 1: mask for position, mask in left.items()}
    assert statistics(translated, height) == {"inside": 3, "outside": 11, "energy": 30}
    split = dict(clipped)
    toggle(split, -3, constant)
    toggle(split, -1, constant)
    assert statistics(split, height) == {"inside": 3, "outside": 7, "energy": 18}
    trace_repair = {-3: constant, -2: constant, -1: constant ^ alternating}
    assert statistics(trace_repair, height) == {"inside": 5, "outside": 7, "energy": 16}
    original_boundary = boundary(left, height)
    repaired_boundary = boundary(trace_repair, height)
    assert original_boundary[0].get(0, 0) == repaired_boundary[0].get(0, 0)
    assert all(original_boundary[1].get(position, 0) == repaired_boundary[1].get(position, 0)
               for position in (0, 1))
    right = {0: constant, 1: constant, 2: alternating, 3: alternating}
    right_before = statistics(right, height)
    right_clipped = {position: mask for position, mask in right.items() if position <= 2}
    assert right_before == {"inside": 2, "outside": 8, "energy": 22}
    assert statistics(right_clipped, height) == {"inside": 2, "outside": 10, "energy": 28}
    assert right_deletion_charge(right, 3, height) == 6
    output = {
        "status": "PASS",
        "height": height,
        "left_original": left_before,
        "left_clipped": statistics(clipped, height),
        "left_single_column_moves": moved_energies,
        "left_whole_translation": statistics(translated, height),
        "left_two_column_repair": statistics(split, height),
        "left_inside_preserving_repair": statistics(trace_repair, height),
        "right_original": right_before,
        "right_clipped": statistics(right_clipped, height),
    }
    output["wall_seconds"] = round(time.monotonic() - started, 4)
    output["peak_rss_bytes_macos"] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    assert output["peak_rss_bytes_macos"] < 768 * 1024**2
    print(json.dumps(output, indent=2))


if __name__ == "__main__":
    main()
