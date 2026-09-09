"""Check the fixed graph bookkeeping in the canonical stability proof."""

import resource
import signal
import sys
import time

from a42_s7_codex_spatial_certificate import build_graph, inverse_fourier_table, recurrent_components


def verify():
    states, edges = build_graph(inverse_fourier_table())
    zero = {edge: record for edge, record in edges.items() if record[1] == 0}
    state_a, state_b, state_c = (0, 0, 0), (1, 0, 0), (0, 1, 1)
    state_d, state_e, state_f = (1, 1, 0), (1, 1, 2), (1, 1, 3)
    inside_cap = {}
    difference_cap = {}
    for edge, (cost, slack, first, second) in zero.items():
        inside_cap[edge] = int(0 in first) + len(second & {0, 1})
        difference_cap[edge] = cost + 1
        assert inside_cap[edge] <= 1
        assert difference_cap[edge] <= 3
    assert {source for source, target in zero if target == state_d} == {state_d, state_e, state_f}
    assert {target for source, target in zero if source == state_c} == {state_a}
    assert {target for source, target in zero if source == state_a} == {state_b, (0, 0, 1)}
    assert not any(source == (0, 0, 1) for source, target in zero)
    without_forbidden = {edge: record for edge, record in zero.items() if edge != (state_d, state_c)}
    assert recurrent_components(states, without_forbidden) == {
        frozenset({state_a, state_b, state_c}), frozenset({state_d}),
        frozenset({(0, 0, 2), (0, 0, 3)})}
    triangle = (state_a, state_b, state_c)
    triangle_holes, triangle_difference = 0, 0
    for phase in range(3):
        walk = [triangle[(phase + offset) % 3] for offset in range(9)]
        assert any(walk[start:start + 7] == [state_a, state_b, state_c, state_a, state_b, state_c, state_a]
                   for start in range(3))
        for length in range(8):
            path = list(zip(walk[:length], walk[1:length + 1]))
            triangle_holes = max(triangle_holes, sum(inside_cap[edge] for edge in path))
            triangle_difference = max(triangle_difference, sum(difference_cap[edge] for edge in path))
    assert (triangle_holes, triangle_difference) == (5, 15)
    branches = ((state_b, state_e, state_d), (state_b, state_f, state_d),
                (state_b, state_e, (0, 1, 3)), (state_b, state_f, (0, 1, 2)),
                (state_a, (0, 0, 1)))
    assert max(sum(inside_cap[edge] for edge in zip(walk, walk[1:])) for walk in branches) == 1
    assert max(sum(difference_cap[edge] for edge in zip(walk, walk[1:])) for walk in branches) == 5
    assigned = {(triangle[index], triangle[(index + 1) % 3]) for index in range(3)}
    assigned |= {edge for walk in branches for edge in zip(walk, walk[1:])}
    assigned |= {(state_d, state_d), ((0, 0, 2), (0, 0, 3)), ((0, 0, 3), (0, 0, 2))}
    expected_transient = {((1, 0, 1), (0, 1, 0)), ((0, 1, 0), (1, 0, 2)), ((0, 1, 0), (1, 0, 3))}
    assert set(without_forbidden) - assigned == expected_transient
    clean_holes = 4 + triangle_holes + 1 + 4
    clean_difference = 4 * 3 + triangle_difference + 5 + 16
    hole_constant = clean_holes + 3
    difference_constant = clean_difference + 7
    assert (clean_holes, clean_difference, hole_constant, difference_constant) == (14, 48, 17, 55)
    assert 2 * 17 - 3 == 31
    assert 6 * 16 - (2 * 17 - 3) == 65
    return {"vertices": len(states), "edges": len(edges), "zero_edges": len(zero),
            "triangle_caps": [triangle_holes, triangle_difference],
            "clean_interval_caps": [clean_holes, clean_difference],
            "stability_constants": [hole_constant, difference_constant],
            "scope": "Fixed graph bookkeeping; spatial recurrence and all-height arguments are proved in the note."}


if __name__ == "__main__":
    started = time.monotonic()
    signal.alarm(55)
    resource.setrlimit(resource.RLIMIT_CPU, (50, 50))
    result = verify()
    peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    result["peak_rss_bytes"] = peak if sys.platform == "darwin" else peak * 1024
    assert result["peak_rss_bytes"] < 768 * 1024 ** 2
    result["elapsed_seconds"] = time.monotonic() - started
    print(result)
