"""Exact two-row polynomial certificates excluding residual D excursions."""

import resource
import signal
import sys
import time
from itertools import product


STATES = {
    "A": (0, 0, 0), "B": (1, 0, 0), "C": (0, 1, 1),
    "D": (1, 1, 0), "E": (1, 1, 2), "F": (1, 1, 3),
    "G": (0, 0, 2), "H": (0, 0, 3), "O": (0, 0, 1),
    "I": (0, 1, 0), "J": (1, 0, 1), "K": (1, 0, 2),
    "L": (1, 0, 3), "M": (0, 1, 2), "N": (0, 1, 3),
    "Z": (1, 1, 1),
}

EXPECTED = {
    "AB": ({0}, set(), 0), "BC": (set(), {1, 2}, 0),
    "BE": ({2}, {1}, 0), "BF": ({1}, {2}, 0),
    "ED": ({2}, set(), 0), "FD": ({1}, set(), 0),
    "DD": ({0}, set(), 0), "DC": (set(), {0}, 0),
    "CA": (set(), set(), 0), "AO": (set(), {0}, 0),
    "GH": (set(), {2}, 0), "HG": (set(), {1}, 0),
    "AG": ({0, 1}, {1}, 1), "AH": ({0, 2}, {2}, 1),
    "OG": ({0, 2}, {1}, 1), "OH": ({0, 1}, {2}, 1),
    "GB": ({2}, set(), 1), "HB": ({1}, set(), 1),
}


def graph_records():
    inverse = {}
    for mask in range(8):
        residues = frozenset(index for index in range(3) if mask & (1 << index))
        field = 0
        for residue in residues:
            field ^= (1, 2, 3)[residue]
        inverse[len(residues) % 2, field] = residues
    assert len(inverse) == 8
    low = {STATES[name] for name in ("A", "O", "I", "B")}
    names = {state: name for name, state in STATES.items()}
    records = {}
    for source_name, source in STATES.items():
        previous, earlier, previous_field = source
        for current, field in product(range(2), range(4)):
            target = current, previous, field
            first = inverse[current, 1 ^ field ^ previous_field]
            second = inverse[1 ^ current ^ previous ^ earlier, field]
            slack = len(first) + len(second) - 1 + int(source not in low) - int(target not in low)
            records[source_name + names[target]] = first, second, slack
    assert len(records) == 128
    assert {edge for edge, record in records.items() if record[2] == 0} == set(
        "AB BC BE BF CA DC DD ED FD GH HG AO IK IL JI EN FM".split()
    )
    for edge, expected in EXPECTED.items():
        assert records[edge] == expected
    return records


def projected_equation(rows, current, period=None):
    terms = {("constant", 0, 0, 0)}
    stencil = (("f", current, 3), ("f", current - 1, 0),
               ("f", current - 2, 0), ("g", current, 0),
               ("g", current - 3, 1), ("g", current - 3, 2))
    for block, physical_row, shift in stencil:
        if period is not None:
            physical_row %= period
        first, second, slack = rows[physical_row]
        residues = first if block == "f" else second
        for residue in residues:
            if (residue + shift) % 3 == 0:
                terms.symmetric_difference_update({(block, physical_row, residue, shift)})
    return terms


def shifted(terms, amount):
    return {(block, physical_row, residue, shift + amount)
            for block, physical_row, residue, shift in terms}


def certify(rows, period=None):
    first = projected_equation(rows, 0, period)
    second = projected_equation(rows, 1, period)
    assert first == {("f", 0, 0, 3), ("constant", 0, 0, 0)}
    assert second == {("f", 0, 0, 0), ("constant", 0, 0, 0)}
    assert first ^ shifted(second, 3) == {
        ("constant", 0, 0, 0), ("constant", 0, 0, 3)
    }


def verify_local_contexts(records):
    count = 0
    for entry in ("AG", "AH", "OG", "OH"):
        destination = entry[1]
        following_edges = (destination + ("H" if destination == "G" else "G"),
                           destination + "B")
        prefixes = [(incoming, "DC", "CA") for incoming in ("DD", "ED", "FD")] \
            if entry[0] == "A" else [("DC", "CA", "AO")]
        for prefix, following in product(prefixes, following_edges):
            word = (*prefix, entry, following)
            assert all(first[1] == second[0] for first, second in zip(word, word[1:]))
            rows = {position - 3: records[edge] for position, edge in enumerate(word)}
            certify(rows)
            count += 1
    assert count == 16
    return count


def verify_eight_row_shapes(records):
    shapes = (("A", 0, 2), ("A", 1, 1), ("O", 0, 1), ("O", 1, 0))
    count = 0
    for (source, middle_count, loop_count), destination, branch in product(
            shapes, ("G", "H"), (("BE", "ED"), ("BF", "FD"))):
        word = [source + destination]
        current = destination
        if middle_count:
            following = "H" if current == "G" else "G"
            word.append(current + following)
            current = following
        word.append(current + "B")
        word.extend(branch)
        word.extend(["DD"] * loop_count)
        word.extend(("DC", "CA"))
        if source == "O":
            word.append("AO")
        assert len(word) == 8
        assert all(first[1] == second[0] for first, second in zip(word, word[1:] + word[:1]))
        rows = {position: records[edge] for position, edge in enumerate(word)}
        assert sum(len(first) + len(second) for first, second, slack in rows.values()) == 10
        assert sum(slack for first, second, slack in rows.values()) == 2
        certify(rows, period=8)
        count += 1
    assert count == 16
    return count


def main():
    signal.alarm(55)
    resource.setrlimit(resource.RLIMIT_CPU, (50, 50))
    started = time.monotonic()
    records = graph_records()
    local_count = verify_local_contexts(records)
    cyclic_count = verify_eight_row_shapes(records)
    print(f"PASS: {local_count} exact local contexts and {cyclic_count} prescribed cyclic templates")
    print("Each has the formal contradiction E_0 + x^3 E_1 = 1 + x^3")
    print("Scope: residual A/O-entry DC excursions excluded; DC-free sharp bound not proved")
    peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    peak_bytes = peak if sys.platform == "darwin" else peak * 1024
    assert peak_bytes <= 768 * 1024 ** 2
    print("Elapsed", round(time.monotonic() - started, 6), "seconds; peak MiB",
          round(peak_bytes / 1024 ** 2, 1))


if __name__ == "__main__":
    main()
