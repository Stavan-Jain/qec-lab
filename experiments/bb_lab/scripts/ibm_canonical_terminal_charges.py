"""Fixed identities for four terminal-entry charge pairs, not a search."""

import resource
import signal
import sys
import time

from a42_s14_codex_double_charge import graph, polynomial_product, projected_terms, reachability


def main():
    signal.alarm(55)
    resource.setrlimit(resource.RLIMIT_CPU, (50, 50))
    started = time.monotonic()
    states, records = graph()
    charges = [edge for edge, record in records.items() if record[1] == 1]
    zero = [edge for edge, record in records.items() if record[1] == 0]
    state_a, state_b, state_c, state_d = (0, 0, 0), (1, 0, 0), (0, 1, 1), (1, 1, 0)
    assert not any(target == state_d for source, target in charges)
    expected = {
        16: ((0, 1, 2), (0, 0, 3), frozenset(), frozenset({0, 1})),
        19: ((0, 1, 3), (0, 0, 2), frozenset(), frozenset({0, 2})),
        9: ((0, 0, 2), state_b, frozenset({2}), frozenset()),
        10: ((0, 0, 3), state_b, frozenset({1}), frozenset()),
    }
    for charge_id, (source, target, first, second) in expected.items():
        assert charges[charge_id] == (source, target)
        assert records[charges[charge_id]][2:] == (first, second)
    without_dc = [edge for edge in zero if edge != (state_d, state_c)]
    closure = reachability(states, without_dc)
    for terminal, branch, residue in (((0, 1, 2), (1, 1, 3), 1),
                                      ((0, 1, 3), (1, 1, 2), 2)):
        relevant = {state for state in closure[state_b] if terminal in closure[state]}
        assert {edge for edge in without_dc if edge[0] in relevant and edge[1] in relevant} == {
            (state_b, state_c), (state_c, state_a), (state_a, state_b),
            (state_b, branch), (branch, terminal),
        }
        rows = {
            -6: (frozenset(), frozenset()),
            -5: (frozenset(), frozenset({1, 2})),
            -4: (frozenset(), frozenset()),
            -3: (frozenset({0}), frozenset()),
            -2: (frozenset({residue}), frozenset({3 - residue})),
            -1: (frozenset(), frozenset({residue})),
        }
        assert projected_terms(rows, -3) == {(0, -3, 0, 3)}
        assert projected_terms(rows, -1) == {(0, -3, 0, 0)}
    shifted_pair = {1, 2}
    assert {0} ^ polynomial_product(shifted_pair, {-2}) == {-1}
    assert {0} ^ polynomial_product(shifted_pair, {-1}) == {1}
    first = {1: {-2}, 2: {-5}}
    second = {0: {0, -5}, 2: {-1}, 3: {-2}}
    syndromes = []
    for row in range(4):
        syndromes.append(
            polynomial_product({3}, first.get(row, set()))
            ^ first.get((row - 1) % 4, set())
            ^ first.get((row - 2) % 4, set())
            ^ second.get(row, set())
            ^ polynomial_product(shifted_pair, second.get((row - 3) % 4, set()))
        )
    assert syndromes[:3] == [{0}, {0}, {0}]
    assert syndromes[3] == {-5, -4, -3, 1, 2}
    assert polynomial_product({3}, {-4}) ^ polynomial_product(shifted_pair, {-2}) == {0}
    assert {-2} ^ {exponent for exponent in polynomial_product(shifted_pair, {-1})
                    if exponent % 3 == 1} == {-2, 1}
    peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    peak_bytes = peak if sys.platform == "darwin" else peak * 1024
    assert peak_bytes <= 768 * 1024 ** 2
    print({"status": "PASS", "excluded_pairs": [(16, 9), (16, 10), (19, 9), (19, 10)],
           "wall_seconds": round(time.monotonic() - started, 6), "peak_rss_bytes": peak_bytes})


if __name__ == "__main__":
    main()
