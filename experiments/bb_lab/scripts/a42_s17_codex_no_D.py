"""Fixed graph, endpoint, and wrapped-row identities; no word or support search."""

import resource
import signal
import sys
import time

sys.dont_write_bytecode = True

from a42_s14_codex_double_charge import graph, polynomial_product, projected_terms


STATE_A = (0, 0, 0)
STATE_O = (0, 0, 1)
STATE_B = (1, 0, 0)
STATE_C = (0, 1, 1)
STATE_D = (1, 1, 0)
STATE_H2 = (0, 0, 2)
STATE_H3 = (0, 0, 3)


def label(state):
    return "".join(map(str, state))


def rows_for_word(word, records):
    return {index: records[(source, target)][2:]
            for index, (source, target) in enumerate(zip(word, word[1:]))}


def fixed_graph_check():
    states, records = graph()
    charges = [edge for edge, record in records.items() if record[1] == 1]
    zero = {edge for edge, record in records.items() if record[1] == 0}
    expected_zero = {
        "000->001", "000->100", "002->003", "003->002", "010->102",
        "010->103", "011->000", "100->011", "100->112", "100->113",
        "101->010", "110->011", "110->110", "112->013", "112->110",
        "113->012", "113->110",
    }
    assert (len(states), len(records), len(charges)) == (16, 128, 32)
    assert {label(source) + "->" + label(target) for source, target in zero} == expected_zero
    expected_charges = {
        0: (STATE_A, STATE_H2, frozenset({0, 1}), frozenset({1})),
        1: (STATE_A, STATE_H3, frozenset({0, 2}), frozenset({2})),
        4: (STATE_O, STATE_H2, frozenset({0, 2}), frozenset({1})),
        5: (STATE_O, STATE_H3, frozenset({0, 1}), frozenset({2})),
        9: (STATE_H2, STATE_B, frozenset({2}), frozenset()),
        10: (STATE_H3, STATE_B, frozenset({1}), frozenset()),
    }
    for charge_id, (source, target, first, second) in expected_charges.items():
        assert charges[charge_id] == (source, target)
        assert records[(source, target)][2:] == (first, second)
    assert {target for source, target in zero if source == STATE_D} == {STATE_D, STATE_C}
    return records, charges


def fixed_triangle_pinch(records):
    for predecessor in (STATE_A, STATE_H2, STATE_H3):
        word = [predecessor, STATE_B, STATE_C, STATE_A, STATE_B, STATE_C, STATE_A]
        rows = rows_for_word(word, records)
        assert not rows[0][1]
        assert projected_terms(rows, 3) == {(0, 3, 0, 3)}
        assert projected_terms(rows, 5) == {(0, 3, 0, 0)}


def fixed_endpoint_identities():
    shifted_pair = {1, 2}
    endpoint_data = (
        (1, -2, -1, -4),
        (2, -1, 1, -2),
    )
    for residue, forced_second, first_at_ca, first_at_exit in endpoint_data:
        projected = {exponent for exponent in polynomial_product(shifted_pair, {forced_second})
                     if exponent % 3 == 0}
        assert forced_second % 3 == residue and projected == {0}
        assert {0} ^ polynomial_product(shifted_pair, {forced_second}) == {first_at_ca}
        assert polynomial_product({3}, {first_at_exit}) == {first_at_ca}
        assert first_at_exit != first_at_ca
        defect = polynomial_product({3}, {first_at_ca}) ^ {first_at_ca}
        assert len(defect) == 2


def fixed_wrapped_words(records, charges):
    short_words = (
        (0, 9, [STATE_A, STATE_H2, STATE_B, STATE_C, STATE_A]),
        (1, 10, [STATE_A, STATE_H3, STATE_B, STATE_C, STATE_A]),
    )
    for entry_id, exit_id, word in short_words:
        rows = rows_for_word(word, records)
        assert (word[0], word[1]) == charges[entry_id]
        assert (word[1], word[2]) == charges[exit_id]
        assert [records[(source, target)][1] for source, target in zip(word, word[1:])] == [1, 1, 0, 0]
        assert sum(len(first) + len(second) for first, second in rows.values()) == 6
        assert projected_terms(rows, 0, 4) == {(0, 0, 0, 3)}
        assert projected_terms(rows, 2, 4) == {(0, 0, 0, 0)}
    endpoint_words = (
        (0, 9, 5, [STATE_A, STATE_H2, STATE_H3, STATE_H2, STATE_H3,
                    STATE_H2, STATE_B, STATE_C, STATE_A]),
        (1, 10, 5, [STATE_A, STATE_H3, STATE_H2, STATE_H3, STATE_H2,
                     STATE_H3, STATE_B, STATE_C, STATE_A]),
        (4, 10, 4, [STATE_O, STATE_H2, STATE_H3, STATE_H2, STATE_H3,
                     STATE_B, STATE_C, STATE_A, STATE_O]),
        (5, 9, 4, [STATE_O, STATE_H3, STATE_H2, STATE_H3, STATE_H2,
                    STATE_B, STATE_C, STATE_A, STATE_O]),
    )
    for entry_id, exit_id, exit_row, word in endpoint_words:
        rows = rows_for_word(word, records)
        assert (word[0], word[1]) == charges[entry_id]
        assert (word[exit_row], word[exit_row + 1]) == charges[exit_id]
        slacks = [records[(source, target)][1] for source, target in zip(word, word[1:])]
        assert [index for index, slack in enumerate(slacks) if slack] == [0, exit_row]
        assert slacks[0] == slacks[exit_row] == 1
        assert sum(len(first) + len(second) for first, second in rows.values()) == 10
        assert not rows[exit_row - 1][0] and not rows[exit_row - 2][0]
        assert not rows[exit_row][1]
        earlier_residue, = rows[exit_row - 3][1]
        later_residue, = rows[exit_row - 1][1]
        assert earlier_residue == later_residue
        shift = 3 - earlier_residue
        assert projected_terms(rows, exit_row, 8) == {(1, exit_row - 3, earlier_residue, shift)}
        assert projected_terms(rows, exit_row + 2, 8) == {(1, exit_row - 1, later_residue, shift)}
        if entry_id in (4, 5):
            assert projected_terms(rows, 7, 8) == {(1, 7, 0, 0)}
    assert polynomial_product({3}, {-3}) == {0}
    assert {-3} != {0}
    print({"fixed_triangle_contexts": 3, "fixed_height_four_words": 2,
           "fixed_height_eight_endpoint_words": 4})


def main():
    signal.alarm(55)
    resource.setrlimit(resource.RLIMIT_CPU, (50, 50))
    started = time.monotonic()
    records, charges = fixed_graph_check()
    fixed_triangle_pinch(records)
    fixed_endpoint_identities()
    fixed_wrapped_words(records, charges)
    peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    peak_bytes = peak if sys.platform == "darwin" else 1024 * peak
    assert peak_bytes <= 768 * 1024 ** 2
    print({"status": "PASS", "wall_seconds": round(time.monotonic() - started, 6),
           "peak_rss_bytes": peak_bytes,
           "scope": "fixed identities only; the note proves the all-height no-DC theorem"})


if __name__ == "__main__":
    main()
