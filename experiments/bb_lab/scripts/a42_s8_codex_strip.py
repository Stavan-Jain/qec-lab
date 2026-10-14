"""Exact fixed six-column strip audit, with no x-periodic identification."""

from __future__ import annotations

import json
import resource
import signal
import time

import numpy as np


def main():
    signal.alarm(55)
    resource.setrlimit(resource.RLIMIT_CPU, (50, 50))
    started = time.monotonic()
    masks = np.arange(64, dtype=np.int16)
    first_cost = np.array([3 * mask.bit_count() - 4 * ((mask >> 3) & 1)
                           for mask in range(256)], dtype=np.int16)
    second_cost = np.array([3 * mask.bit_count() - 4 * ((mask >> 3) & 3).bit_count()
                            for mask in range(512)], dtype=np.int16)
    first_edges = first_cost[masks[:, None] ^ ((masks << 1) ^ (masks << 2))[None, :]]
    second_edges = second_cost[(masks << 3)[:, None, None]
                               ^ masks[None, :, None] ^ masks[None, None, :]]
    potential = np.zeros((64, 64, 64), dtype=np.int16)
    history = []
    expected_changes = (896, 15, 2, 2, 0)
    for iteration in range(1, 6):
        candidate = np.empty_like(potential)
        for current in range(64):
            candidate[current] = second_edges[current] + np.min(
                potential + first_edges[current][None, None, :], axis=2
            )
        candidate = np.minimum(candidate, potential)
        changed = int(np.count_nonzero(candidate != potential))
        history.append({"iteration": iteration, "changed": changed,
                        "minimum": int(candidate.min())})
        assert changed == expected_changes[iteration - 1], f"Unexpected F_{iteration} update"
        if iteration < 5:
            assert changed > 0, f"F_{iteration} must differ from F_{iteration - 1}"
        else:
            assert np.array_equal(candidate, potential), "Certificate requires F_5 == F_4"
        potential = candidate
    values, counts = np.unique(potential, return_counts=True)
    histogram = dict(zip(map(int, values), map(int, counts)))
    assert histogram == {-4: 1, -3: 130, -2: 388, -1: 391, 0: 261234}, "Unexpected F_4 histogram"
    minimum_slack = 32767
    for current in range(64):
        slack = (potential + first_edges[current][None, None, :]
                 + second_edges[current][:, :, None] - potential[current][:, :, None])
        current_minimum = int(slack.min())
        assert current_minimum >= 0, f"Negative edge slack for current row {current}"
        minimum_slack = min(minimum_slack, current_minimum)
    output = {"columns": [-3, -2, -1, 0, 1, 2], "settled": True, "history": history,
              "minimum_slack": minimum_slack, "potential_histogram": histogram}
    output["wall_seconds"] = round(time.monotonic() - started, 4)
    output["peak_rss_bytes_macos"] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    assert output["peak_rss_bytes_macos"] < 768 * 1024**2
    print(json.dumps(output, indent=2))


if __name__ == "__main__":
    main()
