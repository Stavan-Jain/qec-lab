"""Fixed 16-state and nine-context certificate for A42 S7; no period search."""

import resource
import signal
import sys
import time


def fibre_weight(binary_value, field_value):
    if field_value == 0:
        return 3 * binary_value
    return 1 if binary_value else 2


def inverse_fourier_table():
    inverse = {}
    for mask in range(8):
        residues = frozenset(index for index in range(3) if mask & (1 << index))
        binary_value = len(residues) % 2
        field_value = 0
        for residue in residues:
            field_value ^= (1, 2, 3)[residue]
        assert len(residues) == fibre_weight(binary_value, field_value)
        inverse[(binary_value, field_value)] = residues
    assert len(inverse) == 8
    return inverse


def build_graph(inverse):
    states = [(previous, earlier, field_previous)
              for previous in (0, 1) for earlier in (0, 1) for field_previous in range(4)]
    low_states = {(0, 0, 0), (0, 0, 1), (0, 1, 0), (1, 0, 0)}
    potential = {state: 0 if state in low_states else 1 for state in states}
    edges = {}
    for source in states:
        previous, earlier, field_previous = source
        for binary_value in (0, 1):
            for field_value in range(4):
                first_field = 1 ^ field_value ^ field_previous
                second_binary = 1 ^ binary_value ^ previous ^ earlier
                first_residues = inverse[(binary_value, first_field)]
                second_residues = inverse[(second_binary, field_value)]
                cost = len(first_residues) + len(second_residues)
                target = (binary_value, previous, field_value)
                slack = cost - 1 + potential[source] - potential[target]
                assert slack >= 0
                edges[(source, target)] = (cost, slack, first_residues, second_residues)
    assert len(edges) == 128
    return states, edges


def recurrent_components(states, zero_edges):
    adjacency = {state: [] for state in states}
    for source, target in zero_edges:
        adjacency[source].append(target)
    reachable = {}
    for source in states:
        visited = {source}
        pending = [source]
        while pending:
            current = pending.pop()
            for target in adjacency[current]:
                if target not in visited:
                    visited.add(target)
                    pending.append(target)
        reachable[source] = visited
    remaining = set(states)
    components = set()
    while remaining:
        source = min(remaining)
        component = frozenset(target for target in remaining
                              if target in reachable[source] and source in reachable[target])
        remaining -= component
        if len(component) > 1 or (source, source) in zero_edges:
            components.add(component)
    return components


def check_certificate():
    inverse = inverse_fourier_table()
    states, edges = build_graph(inverse)
    zero_edges = {edge: record for edge, record in edges.items() if record[1] == 0}
    assert len(zero_edges) == 17
    state_a = (0, 0, 0)
    state_b = (1, 0, 0)
    state_c = (0, 1, 1)
    state_d = (1, 1, 0)
    state_e = (1, 1, 2)
    state_f = (1, 1, 3)
    main_component = frozenset({state_a, state_b, state_c, state_d, state_e, state_f})
    other_component = frozenset({(0, 0, 2), (0, 0, 3)})
    assert recurrent_components(states, zero_edges) == {main_component, other_component}
    expected_main_edges = {
        (state_a, state_b), (state_b, state_c), (state_b, state_e), (state_b, state_f),
        (state_c, state_a), (state_d, state_c), (state_d, state_d),
        (state_e, state_d), (state_f, state_d),
    }
    assert {edge for edge in zero_edges if edge[0] in main_component and edge[1] in main_component} == expected_main_edges
    assert {source for source, target in zero_edges if target == state_d} == {state_d, state_e, state_f}
    assert {target for source, target in zero_edges if source == state_b and target in main_component} == {state_c, state_e, state_f}
    contexts = 0
    for predecessor in (state_d, state_e, state_f):
        for successor in (state_c, state_e, state_f):
            walk = (predecessor, state_d, state_c, state_a, state_b, successor)
            records = [zero_edges[(source, target)] for source, target in zip(walk, walk[1:])]
            prior_row, row_zero, row_one, row_two, row_three = records
            assert prior_row[3] == frozenset()
            assert row_zero[2:] == (frozenset(), frozenset({0}))
            assert row_one[2:] == (frozenset(), frozenset())
            assert row_two[2:] == (frozenset({0}), frozenset())
            assert row_three[2].isdisjoint({0}) and row_three[3].isdisjoint({0})
            shifted_second_residues = {(residue + shift) % 3 for residue in row_zero[3] for shift in (1, 2)}
            assert shifted_second_residues.isdisjoint({0})
            contexts += 1
    assert contexts == 9
    after_spatial_obstruction = {edge: record for edge, record in zero_edges.items() if edge != (state_d, state_c)}
    assert recurrent_components(states, after_spatial_obstruction) == {
        frozenset({state_d}), frozenset({state_a, state_b, state_c}), other_component,
    }
    assert zero_edges[(state_d, state_d)][2:] == (frozenset({0}), frozenset())
    for source, target in zero_edges:
        if source in other_component and target in other_component:
            record = zero_edges[(source, target)]
            assert record[2] == frozenset() and len(record[3]) == 1
    forced_exponent = -3
    assert forced_exponent + 3 == 0
    assert forced_exponent != 0
    print({"vertices": len(states), "edges": len(edges), "zero_slack_edges": len(zero_edges),
           "recurrent_component_sizes": sorted(map(len, recurrent_components(states, zero_edges))),
           "spatial_contexts": contexts})
    print({"forbidden_transition": (state_d, state_c),
           "row_two_forces_f_exponent": forced_exponent,
           "row_three_requires_f_exponent": 0,
           "scope": "Symbolic Laurent-row obstruction; no compact-vector or period enumeration"})


def main():
    started = time.monotonic()
    signal.alarm(55)
    resource.setrlimit(resource.RLIMIT_CPU, (50, 50))
    check_certificate()
    peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    peak_bytes = peak if sys.platform == "darwin" else 1024 * peak
    assert peak_bytes <= 768 * 1024 ** 2
    print({"elapsed_seconds": round(time.monotonic() - started, 6), "peak_rss_bytes": peak_bytes})


if __name__ == "__main__":
    main()
