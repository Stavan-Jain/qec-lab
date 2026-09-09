"""Fixed 16-state audit for the single slack-two branch; no support search."""

from __future__ import annotations

import json
import resource
import signal
import time
from collections import deque
from itertools import product


def potential(state):
    return int(state not in {(0, 0, 0), (0, 0, 1), (0, 1, 0), (1, 0, 0)})


def label(state):
    return "".join(map(str, state))


def reachable(start, adjacency):
    visited = {start}
    pending = deque([start])
    while pending:
        current = pending.popleft()
        for target in adjacency[current]:
            if target not in visited:
                visited.add(target)
                pending.append(target)
    return visited


def multiply(left, right):
    result = set()
    for first in left:
        for second in right:
            result.symmetric_difference_update({first + second})
    return result


def main():
    signal.alarm(55)
    resource.setrlimit(resource.RLIMIT_CPU, (50, 50))
    started = time.monotonic()
    fourier_support = {}
    for mask in range(8):
        residues = [residue for residue in range(3) if (mask >> residue) & 1]
        field = 0
        for residue in residues:
            field ^= (1, 2, 3)[residue]
        fourier_support[(len(residues) % 2, field)] = residues
    assert len(fourier_support) == 8
    states = list(product(range(2), range(2), range(4)))
    adjacency = {state: [] for state in states}
    edges = []
    for source in states:
        previous, older, previous_field = source
        for current, field in product(range(2), range(4)):
            target = (current, previous, field)
            first = fourier_support[(current, 1 ^ field ^ previous_field)]
            second = fourier_support[(1 ^ current ^ previous ^ older, field)]
            slack = len(first) + len(second) - 1 + potential(source) - potential(target)
            assert slack >= 0
            edges.append((source, target, first, second, slack))
            if slack == 0:
                adjacency[source].append(target)
    reachability = {state: reachable(state, adjacency) for state in states}
    charged = [edge for edge in edges if edge[4] == 2]
    viable = [edge for edge in charged if edge[0] in reachability[edge[1]]]
    assert len(edges) == 128 and len(charged) == 39 and len(viable) == 18
    expected_viable = {
        "001->000", "001->001", "001->100", "002->002", "003->003",
        "010->101", "012->000", "013->000", "100->110", "102->010",
        "103->010", "110->112", "110->113", "111->111", "112->011",
        "112->112", "113->011", "113->113",
    }
    assert {label(source) + "->" + label(target) for source, target, _, _, _ in viable} == expected_viable
    expected_zero = {
        "000->001", "000->100", "002->003", "003->002", "010->102",
        "010->103", "011->000", "100->011", "100->112", "100->113",
        "101->010", "110->011", "110->110", "112->013", "112->110",
        "113->012", "113->110",
    }
    assert {label(source) + "->" + label(target) for source, targets in adjacency.items()
            for target in targets} == expected_zero
    scalar_check = {0, 1, 2}
    shifted_pair = {1, 2}
    assert multiply(shifted_pair, scalar_check) == {1, 4}
    assert multiply(shifted_pair, {-1, 0, 1}) == {0, 3}
    assert multiply({3}, scalar_check) ^ multiply(shifted_pair, shifted_pair) == {2, 3, 5}
    assert {2, 3, 5} != {0}
    left_coefficient = {3} ^ multiply(shifted_pair, scalar_check)
    right_coefficient = {0} ^ multiply(shifted_pair, {0, 3})
    assert left_coefficient == {1, 3, 4} and right_coefficient == {0, 1, 2, 4, 5}
    assert len(left_coefficient) == 3 and len(right_coefficient) == 5
    assert len({-3, 0, 3}) == 3
    assert {0, 3}.isdisjoint({-3, -9})
    fixed_words = [
        ["001", "100", "011", "000", "001"],
        ["012", "000", "100", "113", "012"],
        ["013", "000", "100", "112", "013"],
        ["100", "110", "011", "000", "100"],
    ]
    edge_lookup = {(label(source), label(target)): (first, second, slack)
                   for source, target, first, second, slack in edges}
    for word in fixed_words:
        word_edges = [edge_lookup[(source, target)] for source, target in zip(word, word[1:])]
        assert [slack for _, _, slack in word_edges] == [2, 0, 0, 0]
        assert sum(len(first) + len(second) for first, second, _ in word_edges) == 6
    output = {
        "status": "PASS",
        "fixed_four_row_words_checked": len(fixed_words),
        "viable_edges": [{"edge": label(source) + "->" + label(target),
                          "first_residues": first, "second_residues": second}
                         for source, target, first, second, _ in viable],
        "zero_adjacency": {label(source): list(map(label, targets))
                           for source, targets in adjacency.items() if targets},
        "wall_seconds": round(time.monotonic() - started, 4),
        "peak_rss_bytes_macos": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
    }
    assert output["peak_rss_bytes_macos"] < 768 * 1024**2
    print(json.dumps(output, indent=2))


if __name__ == "__main__":
    main()
