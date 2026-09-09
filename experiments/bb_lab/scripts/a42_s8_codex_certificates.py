"""Tiny exact certificates for the canonical-HM audit; no search."""

from __future__ import annotations

import json
import resource
import signal
import time


def rotate(mask, amount):
    return ((mask << amount) | (mask >> (6 - amount))) & 63


def boundary(rows, folded):
    first = []
    second = []
    for index, current in enumerate(rows):
        oldest = rows[(index - 3) % len(rows)]
        if folded:
            first.append(current ^ rotate(oldest, 1) ^ rotate(oldest, 2))
            second.append(rotate(current, 3) ^ rows[(index - 1) % len(rows)]
                          ^ rows[(index - 2) % len(rows)])
        else:
            first.append(current ^ (oldest << 1) ^ (oldest << 2))
            second.append((current << 3) ^ rows[(index - 1) % len(rows)]
                          ^ rows[(index - 2) % len(rows)])
    return first, second


def statistics(first, second, origin=0):
    inside = sum((mask >> origin) & 1 for mask in first)
    inside += sum(((mask >> origin) & 3).bit_count() for mask in second)
    total = sum(mask.bit_count() for mask in first + second)
    return {"inside": inside, "outside": total - inside, "energy": 3 * total - 4 * inside}


def solve(columns, target):
    pivots = {}
    for index, column in enumerate(columns):
        witness = 1 << index
        while column:
            pivot = column.bit_length() - 1
            if pivot not in pivots:
                pivots[pivot] = (column, witness)
                break
            reducer, certificate = pivots[pivot]
            column ^= reducer
            witness ^= certificate
    answer = 0
    while target:
        pivot = target.bit_length() - 1
        reducer, certificate = pivots[pivot]
        target ^= reducer
        answer ^= certificate
    return answer


def pack(rows):
    return sum(mask << (6 * index) for index, mask in enumerate(rows))


def main():
    signal.alarm(55)
    resource.setrlimit(resource.RLIMIT_CPU, (50, 50))
    started = time.monotonic()
    cycle = [56, 45, 52, 19, 36, 58, 42, 26, 30, 62, 28, 12, 6, 8,
             41, 20, 49, 52, 48, 36, 59, 35, 27, 25, 8, 18, 23, 49,
             44, 51, 37, 42, 33, 29, 62, 4, 15, 7, 25, 43, 18, 22]
    first, second = boundary(cycle, True)
    assert statistics(first, second)["energy"] == -12
    assert first == first[:7] * 6 and second == second[:7] * 6
    columns = []
    for index in range(42):
        basis_rows = [0] * 7
        basis_rows[index // 6] = 1 << (index % 6)
        columns.append(pack(boundary(basis_rows, True)[1]))
    solution = solve(columns, pack(second[:7]))
    short_cycle = [(solution >> (6 * index)) & 63 for index in range(7)]
    short_first, short_second = boundary(short_cycle, True)
    assert short_first == first[:7] and short_second == second[:7]
    output = {"period_seven_rows": short_cycle,
              "period_seven_boundary": [short_first, short_second],
              "period_seven_statistics": statistics(short_first, short_second)}
    padded_first, padded_second = boundary(short_cycle + [0] * 3, True)
    output["seam_rows"] = [[index, padded_first[index], padded_second[index]]
                           for index in [0, 1, 2, 7, 8, 9]]
    assert statistics(padded_first, padded_second) == {
        "inside": 17, "outside": 31, "energy": 76
    }
    power_rows = short_cycle * 71 + [0] * 15
    assert len(power_rows) == 512
    folded = boundary(power_rows, True)
    actual = boundary(power_rows, False)
    output["q512_folded"] = statistics(*folded)
    output["q512_unfolded"] = statistics(*actual)
    assert output["q512_folded"]["energy"] < 0
    assert output["q512_unfolded"]["energy"] >= 0
    assert short_cycle == [61, 10, 32, 13, 53, 3, 47]
    assert output["q512_folded"] == {"inside": 1207, "outside": 381, "energy": -64}
    output["wall_seconds"] = round(time.monotonic() - started, 4)
    output["peak_rss_bytes_macos"] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    assert output["peak_rss_bytes_macos"] < 768 * 1024**2
    print(json.dumps(output, indent=2))


if __name__ == "__main__":
    main()
