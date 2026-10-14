"""Fixed-width Bellman audit for the original-pair canonical HM support."""

from __future__ import annotations

import json
import resource
import signal
import time

import numpy as np


def rotate(mask, amount):
    return ((mask << amount) | (mask >> (6 - amount))) & 63


def main():
    signal.alarm(55)
    resource.setrlimit(resource.RLIMIT_CPU, (50, 50))
    started = time.monotonic()
    masks = np.arange(64, dtype=np.int16)
    forward = np.array([rotate(mask, 1) ^ rotate(mask, 2) for mask in range(64)], dtype=np.int16)
    half_turn = np.array([rotate(mask, 3) for mask in range(64)], dtype=np.int16)
    first_cost = np.array([3 * mask.bit_count() - 4 * (mask & 1) for mask in range(64)], dtype=np.int16)
    second_cost = np.array([3 * mask.bit_count() - 4 * (mask & 3).bit_count()
                            for mask in range(64)], dtype=np.int16)
    first_edges = first_cost[masks[:, None] ^ forward[None, :]]
    second_edges = second_cost[half_turn[:, None, None] ^ masks[None, :, None] ^ masks[None, None, :]]
    potential = np.zeros((64, 64, 64), dtype=np.int16)
    settled = False
    history = []
    predecessors = []
    for iteration in range(1, 129):
        candidate = np.empty_like(potential)
        predecessor = np.empty(potential.shape, dtype=np.uint8)
        for current in range(64):
            options = potential + first_edges[current][None, None, :]
            predecessor[current] = np.argmin(options, axis=2)
            candidate[current] = second_edges[current] + np.min(options, axis=2)
        predecessor[candidate >= potential] = 64
        predecessors.append(predecessor)
        candidate = np.minimum(candidate, potential)
        changed = int(np.count_nonzero(candidate != potential))
        history.append({"iteration": iteration, "changed": changed,
                        "minimum": int(candidate.min())})
        if not changed:
            settled = True
            break
        potential = candidate
        if time.monotonic() - started > 43:
            break
    output = {"states": 64**3, "edges": 64**4, "settled": settled, "history": history}
    if not settled:
        target = tuple(map(int, np.unravel_index(np.argmin(potential), potential.shape)))
        path = [target]
        for predecessor in reversed(predecessors):
            oldest = int(predecessor[target])
            if oldest != 64:
                target = (target[1], target[2], oldest)
                path.append(target)
        path.reverse()
        costs = [int(first_edges[right[0], left[2]]
                     + second_edges[right[0], left[0], left[1]])
                 for left, right in zip(path, path[1:])]
        output["path_cost"] = sum(costs)
        for start in range(len(path)):
            for stop in range(start + 1, len(path)):
                if path[start] == path[stop] and sum(costs[start:stop]) < 0:
                    output["negative_cycle_rows"] = [state[0] for state in path[start:stop]]
                    output["negative_cycle_cost"] = sum(costs[start:stop])
                    break
            if "negative_cycle_rows" in output:
                break
        assert "negative_cycle_rows" in output
        assert output["negative_cycle_cost"] < 0
    if settled:
        minimum_slack = 32767
        for current in range(64):
            slacks = (potential + first_edges[current][None, None, :]
                      + second_edges[current][:, :, None]
                      - potential[current][:, :, None])
            minimum_slack = min(minimum_slack, int(slacks.min()))
        assert minimum_slack >= 0
        values, counts = np.unique(potential, return_counts=True)
        output["minimum_slack"] = minimum_slack
        output["potential_histogram"] = dict(zip(map(int, values), map(int, counts)))
    output["wall_seconds"] = round(time.monotonic() - started, 4)
    output["peak_rss_bytes_macos"] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    assert output["peak_rss_bytes_macos"] < 768 * 1024**2
    print(json.dumps(output, indent=2))


if __name__ == "__main__":
    main()
