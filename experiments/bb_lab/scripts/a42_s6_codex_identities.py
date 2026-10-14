"""Tiny exact checks for A42_S6_codex_residual.md; no searches or data writes."""

import resource
import signal
import sys
import time


def field_mul(left, right):
    product = 0
    while right:
        if right & 1:
            product ^= left
        right >>= 1
        left <<= 1
        if left & 4:
            left ^= 7
    return product


def field_pow(value, exponent):
    result = 1
    for _ in range(exponent):
        result = field_mul(result, value)
    return result


def add(*polynomials):
    result = {}
    for polynomial in polynomials:
        for position, coefficient in polynomial.items():
            result[position] = result.get(position, 0) ^ coefficient
            if not result[position]:
                del result[position]
    return result


def monomial(column, row, period, coefficient=1):
    return {(column, row % period): coefficient}


def multiply(left, right, period):
    result = {}
    for (left_column, left_row), left_coefficient in left.items():
        for (right_column, right_row), right_coefficient in right.items():
            position = (left_column + right_column, (left_row + right_row) % period)
            coefficient = field_mul(left_coefficient, right_coefficient)
            result[position] = result.get(position, 0) ^ coefficient
            if not result[position]:
                del result[position]
    return result


def square(polynomial, period):
    return add(*(monomial(2 * column, 2 * row, period, field_mul(coefficient, coefficient))
                 for (column, row), coefficient in polynomial.items()))


def power(polynomial, exponent, period):
    result = monomial(0, 0, period)
    while exponent:
        if exponent & 1:
            result = multiply(result, polynomial, period)
        polynomial = square(polynomial, period)
        exponent >>= 1
    return result


def norm(polynomial):
    return {position: 1 for position in polynomial}


def evaluate_column(polynomial, value, period):
    return add(*(monomial(0, row, period, field_mul(coefficient, field_pow(value, column % 3)))
                 for (column, row), coefficient in polynomial.items()))


def pairs(period):
    first = add(monomial(3, 0, period), monomial(0, 1, period), monomial(0, 2, period))
    second = add(monomial(0, 0, period), monomial(1, 3, period), monomial(2, 3, period))
    twisted = add(monomial(3, 0, period), monomial(0, 1, period, 2), monomial(0, 2, period, 3))
    return first, second, twisted


def compressed_checks(first_block, second_block, period):
    cyclotomic = add(monomial(0, 0, period), monomial(0, 1, period), monomial(0, 2, period))
    derivative = add(monomial(0, 0, period), monomial(0, 1, period))
    check_one = add(multiply(cyclotomic, evaluate_column(first_block, 1, period), period),
                    evaluate_column(second_block, 1, period))
    check_omega = add(evaluate_column(first_block, 2, period),
                      multiply(derivative, evaluate_column(second_block, 2, period), period))
    return check_one, check_omega


def polarization_check(first_block, second_block, period):
    terms = [
        multiply(monomial(3, 0, period), first_block, period),
        multiply(monomial(0, 1, period, 2), first_block, period),
        multiply(monomial(0, 2, period, 3), first_block, period),
        second_block,
        multiply(monomial(1, 3, period), second_block, period),
        multiply(monomial(2, 3, period), second_block, period),
    ]
    quadratic = {}
    for left_index, left_term in enumerate(terms):
        for right_term in terms[left_index + 1:]:
            for position in left_term.keys() & right_term.keys():
                product = field_mul(left_term[position], field_mul(right_term[position], right_term[position]))
                trace = product ^ field_mul(product, product)
                if trace:
                    quadratic = add(quadratic, {position: trace})
    assert norm(add(*terms)) == add(*(norm(term) for term in terms), quadratic)
    return quadratic


def tower_check(period):
    first, second, twisted = pairs(period)
    identity = monomial(0, 0, period)
    inverse = multiply(monomial(-3 * period, 0, period), power(first, period - 1, period), period)
    transfer = multiply(second, inverse, period)
    depth = period.bit_length() - 1
    assert multiply(first, inverse, period) == identity
    assert len(inverse) == 3 ** depth
    assert len(transfer) == 3 ** (depth + 1)
    expected_power = add(*(monomial(-multiple * period, 0, period) for multiple in (1, 2, 3)))
    assert power(transfer, period, period) == expected_power
    assert multiply(first, transfer, period) == second
    repeated = {(0, row): 1 for row in range(period)}
    omega_first = repeated
    omega_second = add(repeated, multiply(monomial(1, 0, period), repeated, period))
    assert not add(multiply(twisted, omega_first, period), multiply(second, omega_second, period))
    support_first, support_second = norm(omega_first), norm(omega_second)
    assert add(multiply(first, support_first, period), multiply(second, support_second, period)) == repeated
    assert compressed_checks(support_first, support_second, period) == (repeated, repeated)
    assert polarization_check(omega_first, omega_second, period) == repeated
    correction_chain = multiply(add(monomial(-3, 0, period), monomial(-2, 0, period)), repeated, period)
    correction = (multiply(second, correction_chain, period), multiply(first, correction_chain, period))
    assert sum(len(add(support, block)) for support, block in zip((support_first, support_second), correction)) == period
    test_chains = [identity, inverse, correction_chain,
                   add(monomial(-2, 1, period), monomial(1, 0, period), monomial(3, 2, period))]
    support_size = 3 * period
    for chain in test_chains:
        barren_first = multiply(second, chain, period)
        barren_second = multiply(first, chain, period)
        assert compressed_checks(barren_first, barren_second, period) == ({}, {})
        barren_size = len(barren_first) + len(barren_second)
        overlap = len(support_first.keys() & barren_first.keys()) + len(support_second.keys() & barren_second.keys())
        weight = 2 * support_size + 3 * barren_size - 4 * overlap
        assert weight >= 4 * period
        if overlap == support_size:
            assert weight >= 6 * period
        omega_boundary_first = multiply(second, chain, period)
        omega_boundary_second = multiply(twisted, chain, period)
        polarization_check(add(omega_first, omega_boundary_first), add(omega_second, omega_boundary_second), period)
    return {"q": period, "inverse_weight": len(inverse), "graph_impulse_weight": len(transfer) + 1,
            "proved_mixed_bound": 4 * period}


def mirror_check():
    first = {(0, 0): 1, (0, 1): 1, (3, -1): 1}
    second = {(0, 0): 1, (1, 0): 1, (-1, -3): 1}

    def theta(polynomial):
        return {(row, -column): coefficient for (column, row), coefficient in polynomial.items()}

    antipode_first = {(-column, -row): coefficient for (column, row), coefficient in first.items()}
    assert theta(first) == second
    assert theta(second) == antipode_first
    for period in (3, 6, 12):
        syndrome = add(multiply(first, antipode_first, period), square(second, period))
        assert syndrome[(2, 0)] == 1


def fibre_weight(binary_value, field_value):
    if not field_value:
        return 3 * binary_value
    return 1 if binary_value else 2


def zero_slack_components(states, edges, potential):
    zero_edges = [edge for edge in edges if edge[2] - 1 + potential[edge[0]] - potential[edge[1]] == 0]
    adjacency = {state: [] for state in states}
    for source, target, _, _ in zero_edges:
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
    unassigned = set(states)
    recurrent = []
    while unassigned:
        source = min(unassigned)
        component = {target for target in unassigned
                     if target in reachable[source] and source in reachable[target]}
        unassigned -= component
        internal_edges = [edge for edge in zero_edges if edge[0] in component and edge[1] in component]
        if len(component) > 1 or internal_edges:
            recurrent.append({"states": sorted(component), "edges": internal_edges})
    expected_components = {
        frozenset({(0, 0, 0), (0, 1, 1), (1, 0, 0), (1, 1, 0), (1, 1, 2), (1, 1, 3)}),
        frozenset({(0, 0, 2), (0, 0, 3)}),
    }
    assert len(zero_edges) == 17
    assert {frozenset(component["states"]) for component in recurrent} == expected_components
    eight_step_walk = [(0, 0, 0), (1, 0, 0), (1, 1, 2), (1, 1, 0),
                       (1, 1, 0), (1, 1, 0), (1, 1, 0), (0, 1, 1), (0, 0, 0)]
    edge_costs = {(source, target): cost for source, target, cost, _ in zero_edges}
    assert sum(edge_costs[(source, target)] for source, target in zip(eight_step_walk, eight_step_walk[1:])) == 8
    print({"zero_slack_edges": zero_edges, "recurrent_components": recurrent})
    print({"nonconstant_eight_step_equality_walk": eight_step_walk,
           "scope": "Compressed graph witness only; no actual cylinder error lift asserted"})


def compressed_graph_check():
    states = [(previous, earlier, field_previous)
              for previous in (0, 1) for earlier in (0, 1) for field_previous in range(4)]
    edges = []
    for source in states:
        previous, earlier, field_previous = source
        for binary_value in (0, 1):
            for field_value in range(4):
                first_field = 1 ^ field_value ^ field_previous
                second_binary = 1 ^ binary_value ^ previous ^ earlier
                cost = fibre_weight(binary_value, first_field) + fibre_weight(second_binary, field_value)
                target = (binary_value, previous, field_value)
                edges.append((source, target, cost, (binary_value, field_value)))
    potential = dict.fromkeys(states, 0)
    predecessors = {}
    changed = None
    for _ in range(len(states)):
        changed = None
        for source, target, cost, choice in edges:
            candidate = potential[source] + cost - 1
            if potential[target] > candidate:
                potential[target] = candidate
                predecessors[target] = (source, cost, choice)
                changed = target
        if changed is None:
            low_states = {(0, 0, 0), (0, 0, 1), (0, 1, 0), (1, 0, 0)}
            assert potential == {state: -1 if state in low_states else 0 for state in states}
            assert all(cost - 1 + potential[source] - potential[target] >= 0
                       for source, target, cost, _ in edges)
            assert ((1, 1, 0), (1, 1, 0), 1, (1, 0)) in edges
            print({"compressed_graph": "Bellman certificate", "states": len(states), "edges": len(edges),
                   "potential": [(state, potential[state]) for state in states]})
            zero_slack_components(states, edges, potential)
            return
    repeated_state = changed
    for _ in states:
        repeated_state = predecessors[repeated_state][0]
    negative_cycle = []
    current = repeated_state
    while True:
        source, cost, choice = predecessors[current]
        negative_cycle.append((source, current, cost, choice))
        current = source
        if current == repeated_state:
            break
    negative_cycle.reverse()
    assert sum(edge[2] - 1 for edge in negative_cycle) < 0
    print({"compressed_graph": "negative cycle", "cycle": negative_cycle,
           "total_cost": sum(edge[2] for edge in negative_cycle), "length": len(negative_cycle)})
    raise AssertionError("The banked Bellman certificate must have no negative cycle")


def main():
    started = time.monotonic()
    signal.alarm(55)
    resource.setrlimit(resource.RLIMIT_CPU, (50, 50))
    compressed_graph_check()
    mirror_check()
    for period in (1, 2, 4, 8, 16, 32):
        print(tower_check(period))
    peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    peak_bytes = peak if sys.platform == "darwin" else 1024 * peak
    assert peak_bytes < 1024 ** 3
    print({"elapsed_seconds": round(time.monotonic() - started, 3), "peak_rss_bytes": peak_bytes})
    print("All finite identity checks pass; the infinite claims use the proofs in the note.")


if __name__ == "__main__":
    main()
