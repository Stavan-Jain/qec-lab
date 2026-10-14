"""Fixed graph and Laurent-identity checks for the excess-two shell."""

from __future__ import annotations

import json
import resource
import signal
import time
from collections import Counter, deque
from itertools import product


def weight(augmentation, field):
    if field == 0:
        return 3 * augmentation
    return 2 - augmentation


def potential(state):
    return int(state not in {(0, 0, 0), (0, 0, 1), (0, 1, 0), (1, 0, 0)})


def shift(polynomial, amount):
    return {exponent + amount for exponent in polynomial}


def multiply(left, right):
    result = set()
    for first in left:
        for second in right:
            result.symmetric_difference_update({first + second})
    return result


def residue_zero(polynomial):
    return {exponent for exponent in polynomial if exponent % 3 == 0}


def pair_corrections(block, row, pair):
    first = set()
    second = set()
    if block == "f":
        if row == 2:
            first = shift(pair, 3)
        elif row in (0, 1):
            first = set(pair)
        if row == 3:
            second = shift(residue_zero(pair), 3)
        elif row in (1, 2):
            second = residue_zero(pair)
    else:
        if row == 2:
            first = set(pair)
        elif row == -1:
            first = multiply({1, 2}, pair)
        if row == 3:
            second = residue_zero(pair)
        elif row == 0:
            second = residue_zero(multiply({1, 2}, pair))
    return first, second


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


def main():
    signal.alarm(55)
    resource.setrlimit(resource.RLIMIT_CPU, (50, 50))
    started = time.monotonic()
    states = list(product(range(2), range(2), range(4)))
    adjacency = {state: [] for state in states}
    edges = []
    for source in states:
        previous, older, previous_field = source
        for current, field in product(range(2), range(4)):
            target = (current, previous, field)
            cost = weight(current, 1 ^ field ^ previous_field)
            cost += weight(1 ^ current ^ previous ^ older, field)
            slack = cost - 1 + potential(source) - potential(target)
            assert slack >= 0
            edges.append((source, target, slack))
            if slack == 0:
                adjacency[source].append(target)
    assert len(edges) == 128
    assert sum(slack == 0 for _, _, slack in edges) == 17
    main_component = {(0, 0, 0), (1, 0, 0), (0, 1, 1),
                      (1, 1, 0), (1, 1, 2), (1, 1, 3)}
    two_cycle = {(0, 0, 2), (0, 0, 3)}
    reachability = {state: reachable(state, adjacency) for state in states}
    components = {frozenset(target for target in states
                            if target in reachability[state] and state in reachability[target])
                  for state in states}
    recurrent = {component for component in components
                 if len(component) > 1 or any(state in adjacency[state] for state in component)}
    assert recurrent == {frozenset(main_component), frozenset(two_cycle)}
    source_d = (1, 1, 0)
    target_c = (0, 1, 1)
    distances = {target_c: 0}
    pending = deque([target_c])
    while pending:
        current = pending.popleft()
        for target in adjacency[current]:
            if target not in distances:
                distances[target] = distances[current] + 1
                pending.append(target)
    assert distances[source_d] == 4
    adjacency[source_d].remove(target_c)
    reduced_reachability = {state: reachable(state, adjacency) for state in states}
    reduced_components = {frozenset(target for target in states
                                    if target in reduced_reachability[state]
                                    and state in reduced_reachability[target])
                          for state in states}
    reduced_recurrent = {component for component in reduced_components
                         if len(component) > 1 or any(state in adjacency[state] for state in component)}
    assert reduced_recurrent == {
        frozenset({source_d}), frozenset({(0, 0, 0), (1, 0, 0), target_c}), frozenset(two_cycle)
    }
    templates = [
        ("f", 0, {0, 3}),
        ("f", 3, {-6, -3}),
        ("g", 0, {-5, -2}),
        ("g", 0, {-4, -1}),
        ("g", 2, {0, 3}),
        ("g", 3, {-3, 0}),
    ]
    for block, row, pair in templates:
        assert len(pair) == 2 and len({exponent % 3 for exponent in pair}) == 1
        first, second = pair_corrections(block, row, pair)
        assert residue_zero(first) == first
        assert first ^ shift(second, 3) == {0, 3}
    assert {-6, -3} != {0, 3}
    assert multiply({1, 2}, {0, 3}) == {1, 2, 4, 5}
    assert {-3} ^ {-6, -3} == {-6}
    assert len({0} ^ {-6, -3}) == 3
    scalar_pair = {0, 1, 2}
    scalar_remainder = {0, 1, 2, 4, 5}
    assert multiply(scalar_pair, {0, 1}) == {0, 3}
    assert multiply(scalar_remainder, {0, 1}) == {0, 3, 4, 6}
    output = {
        "status": "PASS",
        "scope": "128 fixed graph edges and six fixed spatial correction identities; no support search",
        "slack_histogram": dict(sorted(Counter(slack for _, _, slack in edges).items())),
        "minimum_separation_of_D_to_C_edges": distances[source_d] + 1,
        "verified_pair_templates": len(templates),
        "wall_seconds": round(time.monotonic() - started, 4),
        "peak_rss_bytes_macos": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
    }
    assert output["peak_rss_bytes_macos"] < 768 * 1024**2
    print(json.dumps(output, indent=2))


if __name__ == "__main__":
    main()
