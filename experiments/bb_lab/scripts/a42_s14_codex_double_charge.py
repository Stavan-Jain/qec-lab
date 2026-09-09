"""Fixed graph and finite Laurent certificates for the two-charge branch."""

import resource
import signal
import sys
import time
from itertools import product


def graph():
    inverse = {}
    for mask in range(8):
        residues = frozenset(index for index in range(3) if mask & (1 << index))
        field = 0
        for residue in residues:
            field ^= (1, 2, 3)[residue]
        inverse[(len(residues) % 2, field)] = residues
    low = {(0, 0, 0), (0, 0, 1), (0, 1, 0), (1, 0, 0)}
    states = list(product(range(2), range(2), range(4)))
    records = {}
    for source in states:
        previous, earlier, previous_field = source
        for current, field in product(range(2), range(4)):
            target = (current, previous, field)
            first = inverse[(current, 1 ^ field ^ previous_field)]
            second = inverse[(1 ^ current ^ previous ^ earlier, field)]
            cost = len(first) + len(second)
            slack = cost - 1 + int(source not in low) - int(target not in low)
            assert slack >= 0
            records[(source, target)] = (cost, slack, first, second)
    return states, records


def reachability(states, edges):
    adjacency = {state: [] for state in states}
    for source, target in edges:
        adjacency[source].append(target)
    closure = {}
    for source in states:
        pending = [source]
        visited = {source}
        while pending:
            current = pending.pop()
            for target in adjacency[current]:
                if target not in visited:
                    visited.add(target)
                    pending.append(target)
        closure[source] = visited
    return closure


def polynomial_product(left, right):
    result = set()
    for first in left:
        for second in right:
            result.symmetric_difference_update({first + second})
    return result


def path_lengths(states, edges, source, target):
    closure = reachability(states, edges)
    relevant = {state for state in closure[source] if target in closure[state]}
    remaining = {state: {prior for prior, later in edges if later == state and prior in relevant}
                 for state in relevant}
    lengths = {state: set() for state in relevant}
    lengths[source].add(0)
    while remaining:
        ready = [state for state, incoming in remaining.items() if not incoming]
        assert ready, "A cycle was encountered in an allegedly finite path family"
        for current in ready:
            del remaining[current]
            for prior, later in edges:
                if prior == current and later in remaining:
                    lengths[later].update(length + 1 for length in lengths[current])
                    remaining[later].discard(current)
    return lengths[target]


def path_residues(states, edges, source, target, modulus):
    adjacency = {state: [] for state in states}
    for prior, later in edges:
        adjacency[prior].append(later)
    visited = {(source, 0)}
    pending = [(source, 0)]
    while pending:
        current, residue = pending.pop()
        for later in adjacency[current]:
            lifted = later, (residue + 1) % modulus
            if lifted not in visited:
                visited.add(lifted)
                pending.append(lifted)
    return {residue for state, residue in visited if state == target}


def verify_dc_prefixes(records, charged):
    state_d, state_c = (1, 1, 0), (0, 1, 1)
    assert not any(target == state_d for source, target in charged)
    prefixes = {}
    pending = [(state_c, ())]
    for depth in range(1, 4):
        following = []
        for current, path in pending:
            for edge, record in records.items():
                if edge[0] != current:
                    continue
                if record[1] == 1:
                    prefixes[charged.index(edge)] = path + (edge,)
                elif record[1] == 0:
                    following.append((edge[1], path + (edge,)))
        pending = following
    expected = {index: 2 for index in range(4)}
    expected.update({index: 3 for index in range(4, 9)})
    expected.update({index: 1 for index in (13, 14, 15)})
    expected[22] = 3
    assert {index: len(path) for index, path in prefixes.items()} == expected
    assert records[charged[22]][2:] == (frozenset({1, 2}), frozenset())
    return set(expected) - {22}


def verify_pair_partition(states, records, charged, zero, ordered_pairs):
    pair_ids = {tuple(sorted((charged.index(first), charged.index(second))))
                for first, second in ordered_pairs}
    assert len(pair_ids) == 46 and all(first != second for first, second in pair_ids)
    groups = {
        "terminal_forbidden_DC": {(17, 24), (17, 28), (20, 25), (20, 31), (18, 28), (21, 31)},
        "terminal_nonbinary_length": {(11, 17), (11, 20), (17, 25), (17, 31), (20, 24), (20, 28)},
        "terminal_four_row_contradiction": {(17, 29), (20, 26)},
        "O_to_J": {(6, 11), (6, 12), (6, 24), (6, 25), (6, 28), (6, 31)},
        "K_to_E_or_L_to_F": {(2, 28), (3, 31), (7, 28), (8, 31), (13, 28),
                                (13, 31), (14, 28), (15, 31), (22, 28), (22, 31)},
        "I_to_A": {(11, 13), (11, 22)},
        "opposite_residue_bend": {(13, 24), (13, 25)},
    }
    eliminated = set().union(*groups.values())
    assert len(eliminated) == sum(map(len, groups.values())) == 34
    residual = {tuple(sorted((entry, departure))) for entry in (0, 1, 4, 5, 16, 19)
                for departure in (9, 10)}
    assert len(residual) == 12
    assert not eliminated & residual
    assert eliminated | residual == pair_ids
    eligible = verify_dc_prefixes(records, charged)
    without_dc = [edge for edge in zero if edge != ((1, 1, 0), (0, 1, 1))]
    closure = reachability(states, without_dc)
    for first, second in groups["terminal_forbidden_DC"]:
        assert first not in eligible and second not in eligible
        first_edge, second_edge = charged[first], charged[second]
        assert not (second_edge[0] in closure[first_edge[1]]
                    and first_edge[0] in closure[second_edge[1]])
    for first, second in ((11, 17), (11, 20)):
        first_edge, second_edge = charged[first], charged[second]
        forward = path_residues(states, without_dc, first_edge[1], second_edge[0], 3)
        backward = path_residues(states, without_dc, second_edge[1], first_edge[0], 3)
        assert {(2 + left + right) % 3 for left in forward for right in backward} == {0}
    exact_lengths = {(17, 25): 3, (17, 31): 5, (20, 24): 3, (20, 28): 5,
                     (17, 29): 4, (20, 26): 4, (6, 12): 3}
    for (first, second), expected in exact_lengths.items():
        first_edge, second_edge = charged[first], charged[second]
        forward = path_lengths(states, without_dc, first_edge[1], second_edge[0])
        backward = path_lengths(states, without_dc, second_edge[1], first_edge[0])
        assert {2 + left + right for left in forward for right in backward} == {expected}
    source_a, target_o = (0, 0, 0), (0, 0, 1)
    relevant = {state for state in closure[source_a] if target_o in closure[state]}
    assert relevant == {source_a, target_o, (1, 0, 0), (0, 1, 1)}
    print("Fixed pair ledger:", {label: len(pairs) for label, pairs in groups.items()},
          "unresolved_pairs", sorted(residual))


def projected_terms(rows, index, period=None, wanted=0):
    result = set()
    for block, row, shift in ((0, index, 3), (0, index - 1, 0), (0, index - 2, 0),
                              (1, index, 0), (1, index - 3, 1), (1, index - 3, 2)):
        physical_row = row if period is None else row % period
        for residue in rows[physical_row][block]:
            if (residue + shift) % 3 == wanted:
                result.symmetric_difference_update({(block, physical_row, residue, shift)})
    return result


def verify_pinch_contexts(records, charged):
    state_a, state_b, state_c, state_d = (0, 0, 0), (1, 0, 0), (0, 1, 1), (1, 1, 0)
    state_e, state_f = (1, 1, 2), (1, 1, 3)
    state_o, state_i, state_j = (0, 0, 1), (0, 1, 0), (1, 0, 1)
    state_k, state_l = (1, 0, 2), (1, 0, 3)
    abbreviations = {
        "DC": (state_d, state_c), "CA": (state_c, state_a), "AB": (state_a, state_b),
        "AO": (state_a, state_o), "JI": (state_j, state_i), "IK": (state_i, state_k),
        "IL": (state_i, state_l), "ED": (state_e, state_d), "FD": (state_f, state_d),
        "DD": (state_d, state_d), "BC": (state_b, state_c), "BE": (state_b, state_e),
        "BF": (state_b, state_f),
    }
    contexts = [
        (["DC", "CA", 2, 28, "ED"], 3, 3),
        (["DC", "CA", 3, 31, "FD"], 3, 3),
        (["DC", "CA", "AO", 7, 28, "ED"], 4, 4),
        (["DC", "CA", "AO", 8, 31, "FD"], 4, 4),
        ([13, "JI", "IK", 28, "ED"], 3, 3),
        ([13, "JI", "IL", 31, "FD"], 3, 3),
    ]
    for predecessor in ("DD", "ED", "FD"):
        contexts.extend([
            ([predecessor, "DC", 14, 28, "ED"], 3, 3),
            ([predecessor, "DC", 15, 31, "FD"], 3, 3),
            ([predecessor, "DC", "CA", "AB", 22], 3, 3),
        ])
    for successor in ("JI", 24, 25):
        contexts.append((["DC", "CA", "AO", 6, successor], 3, 3))
    for successor in ("BC", "BE", "BF"):
        contexts.extend([
            ([13, "JI", 11, "AB", successor], 3, 3),
            (["AB", 22, 11, "AB", successor], 3, 3),
        ])
    for path, equation, monomial_row in contexts:
        edges = [charged[label] if isinstance(label, int) else abbreviations[label] for label in path]
        assert all(first[1] == second[0] for first, second in zip(edges, edges[1:]))
        rows = {index: records[edge][2:] for index, edge in enumerate(edges)}
        assert projected_terms(rows, equation) == {(0, monomial_row, 0, 3)}
        assert projected_terms(rows, equation + 1) == {(0, monomial_row, 0, 0)}
    cyclic = [abbreviations["AO"], charged[6], abbreviations["JI"], charged[11]]
    rows = {index: records[edge][2:] for index, edge in enumerate(cyclic)}
    assert projected_terms(rows, 1, 4) == {(0, 1, 0, 3)}
    assert projected_terms(rows, 2, 4) == {(0, 1, 0, 0)}
    print("Exact residue-zero monomial-pinch contexts PASS:", len(contexts) + 1)


def verify_spatial_constants():
    assert polynomial_product({1, 2}, {-1}) == {0, 1}
    assert polynomial_product({1, 2}, {-2}) == {-1, 0}
    assert polynomial_product({1, 2}, {0, 1}) == {1, 3}
    assert {5, -1, -2} ^ polynomial_product({1, 2}, {0, 1}) == {-2, -1, 1, 3, 5}
    assert {-2, -1, 1, 3, 5} != {0}
    assert {4} != {1}
    assert {-3} != {0}
    assert {1 + 3} != {1}
    assert polynomial_product({0, 3}, {0}) == {0, 3}
    assert max(4 + 9 + 3 * 2, 5 + 9 + 3 * 2,
               6 + 9 + 3 * 2 + 7, 7 + 9 + 3 * 2 + 7) == 29
    print("Fixed Laurent contradictions and the residual height bound 29 PASS")


def verify_opposite_bends(records, charged):
    state_d, state_c = (1, 1, 0), (0, 1, 1)
    for charge_index, exit_state, residue in ((24, (1, 1, 2), 1), (25, (1, 1, 3), 2)):
        opposite = 3 - residue
        base = [(state_d, state_c), charged[13], charged[charge_index], (exit_state, state_d)]
        for closing in (False, True):
            path = base if closing else base + [(state_d, state_d)]
            rows = {index: records[edge][2:] for index, edge in enumerate(path)}
            period = 4 if closing else None
            assert projected_terms(rows, 3, period, residue) == {
                (0, 2, residue, 0), (1, 0, 0, residue)}
            assert projected_terms(rows, 3, period, opposite) == {
                (0, 3, opposite, 3), (1, 0, 0, opposite)}
            assert projected_terms(rows, 4, period, residue) == {
                (0, 2, residue, 0), (1, 1, 0, residue)}
            assert projected_terms(rows, 4, period, opposite) == {
                (0, 3, opposite, 0), (1, 1, 0, opposite)}
    print("Opposite-residue bend contexts, including the height-four cyclic endpoint, PASS")


def main():
    signal.alarm(55)
    resource.setrlimit(resource.RLIMIT_CPU, (50, 50))
    started = time.monotonic()
    states, records = graph()
    zero = [edge for edge, record in records.items() if record[1] == 0]
    charged = [edge for edge, record in records.items() if record[1] == 1]
    closure = reachability(states, zero)
    pairs = [(first, second) for first in charged for second in charged
             if second[0] in closure[first[1]] and first[0] in closure[second[1]]]
    assert (len(states), len(records), len(zero), len(charged), len(pairs)) == (16, 128, 17, 32, 92)
    verify_pair_partition(states, records, charged, zero, pairs)
    verify_pinch_contexts(records, charged)
    verify_opposite_bends(records, charged)
    verify_spatial_constants()
    print("Scope: fixed graph and Laurent constants only; 12 pair types remain unresolved")
    peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    peak_bytes = peak if sys.platform == "darwin" else peak * 1024
    assert peak_bytes <= 768 * 1024 ** 2
    print("PASS", round(time.monotonic() - started, 5), "seconds; peak MiB", round(peak_bytes / 1024 ** 2, 1))


if __name__ == "__main__":
    main()
