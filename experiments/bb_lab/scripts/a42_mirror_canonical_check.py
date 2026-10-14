"""Fixed graph contexts and polynomial identities for A42_MIRROR_CANONICAL.md."""

import resource
import signal
import sys
import time

from a42_s6_codex_identities import add, monomial, multiply, power
from a42_s7_codex_spatial_certificate import build_graph, inverse_fourier_table


def check_contexts():
    _, edges = build_graph(inverse_fourier_table())
    state_a, state_b, state_c = (0, 0, 0), (1, 0, 0), (0, 1, 1)
    state_d, state_e, state_f = (1, 1, 0), (1, 1, 2), (1, 1, 3)
    contexts = 0
    for predecessor in (state_d, state_e, state_f):
        for successor in (state_c, state_e, state_f):
            following = state_a if successor == state_c else state_d
            walk = (predecessor, state_d, state_c, state_a, state_b, successor, following)
            records = [edges[source, target] for source, target in zip(walk, walk[1:])]
            assert all(record[1] == 0 for record in records)
            prior, first, second, third, fourth, fifth = records
            assert not prior[3]
            assert not first[2]
            assert not second[2] and not second[3]
            assert third[2] == frozenset({0}) and not third[3]
            assert 0 not in fourth[2] and 0 not in fifth[2]
            assert not fifth[3]
            contexts += 1
    assert contexts == 9
    return contexts


def check_period(period):
    first = add(monomial(3, 0, period), monomial(3, 1, period), monomial(0, 2, period))
    second = add(monomial(0, 0, period), monomial(1, 3, period), monomial(2, 3, period))
    twisted = add(monomial(3, 0, period), monomial(3, 1, period, 2), monomial(0, 2, period, 3))
    norm = add(*(monomial(0, row, period) for row in range(period)))
    second_block = multiply(add(monomial(0, 0, period), monomial(1, 0, period)), norm, period)
    target = multiply(monomial(3, 0, period), norm, period)
    twisted_first_block = multiply(monomial(0, 0, period, 2), norm, period)
    assert not add(multiply(twisted, twisted_first_block, period),
                   multiply(second, second_block, period))
    assert add(multiply(first, norm, period), multiply(second, second_block, period)) == target
    assert multiply(first, target, period) == target
    boundary_first = multiply(second, second_block, period)
    boundary_second = multiply(first, second_block, period)
    assert boundary_first == add(norm, target)
    assert boundary_second == second_block
    assert len(boundary_first) + len(boundary_second) == 4 * period
    if period & (period - 1) == 0:
        assert power(first, period, period) == monomial(0, 0, period)
        inverse = power(first, period - 1, period)
        assert multiply(first, inverse, period) == monomial(0, 0, period)
    return period


def main():
    signal.alarm(55)
    started = time.monotonic()
    contexts = check_contexts()
    periods = [check_period(period) for period in (1, 2, 4, 5, 7, 8)]
    peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    peak_bytes = peak if sys.platform == "darwin" else peak * 1024
    assert peak_bytes < 768 * 1024 ** 2
    print({"extended_spatial_contexts": contexts, "identity_check_periods": periods,
           "elapsed_seconds": round(time.monotonic() - started, 4),
           "peak_rss_bytes": peak_bytes,
           "tier": "finite conventions; universal uniqueness proved in the note"})


if __name__ == "__main__":
    main()
