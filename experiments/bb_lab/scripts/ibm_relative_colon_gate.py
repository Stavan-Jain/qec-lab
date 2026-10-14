"""Reverify banked gross witnesses and their separate block-coset costs."""

import json
from pathlib import Path
import resource
import signal
import sys
import time

from a43_s4_codex_block_purity import apply_columns, independent_basis, matrices, outside_span, positions


def remainder(vector, pivots):
    for pivot in sorted(pivots, reverse=True):
        if vector & (1 << pivot):
            vector ^= pivots[pivot]
    return vector


def clear_x_cycle(cells, width, height):
    size = width * height
    result = 0
    for block, horizontal, vertical in cells:
        if block == 0:
            target = size + ((-horizontal - 1) % width) * height + (-vertical - 3) % height
        else:
            target = ((-horizontal) % width) * height + (-vertical - 1) % height
        result ^= 1 << target
    assert result.bit_count() == len(cells)
    return result


def coordinates(vector, height):
    return [list(divmod(position, height)) for position in positions(vector)]


def verify():
    lab = Path(__file__).resolve().parent.parent
    sys.path.insert(0, str(lab / "src"))
    from bb_lab.tower import validate_banked

    validate_banked(lab / "data")
    data = json.loads((lab / "data/a43/s1_l3_shape.json").read_text())
    objects = data["frames"]["gross(12,6)"]["objects"]
    width, height = 12, 6
    size = width * height
    mask = (1 << size) - 1
    first, second, differential, boundary = matrices(width, height)
    first_basis = independent_basis(first)
    second_basis = independent_basis(second)
    boundary_basis = independent_basis(boundary)
    fixed_sources = [sum(1 << position for position in support) for support in (
        (0, 61, 68, 75, 77, 81, 134, 135, 138, 141, 142, 143),
        (10, 12, 17, 19, 26, 34, 43, 44, 60, 83, 84, 109))]
    for vector in fixed_sources:
        assert vector.bit_count() == 12
        assert apply_columns(differential, vector) == 0
        assert outside_span(vector, boundary_basis)
    fixed_first = fixed_sources[0] & mask
    fixed_second = fixed_sources[1] >> size
    assert fixed_first.bit_count() == fixed_second.bit_count() == 3
    assert remainder(fixed_first, second_basis) == remainder(fixed_sources[1] & mask, second_basis) != 0
    assert remainder(fixed_second, first_basis) == remainder(fixed_sources[0] >> size, first_basis) != 0
    assert apply_columns(first, fixed_first) != apply_columns(second, fixed_second)
    classes = {}
    principal_count = 0
    for record in objects:
        vector = clear_x_cycle(record["cells"], width, height)
        assert vector.bit_count() == record["w"] == 12
        assert apply_columns(differential, vector) == 0
        assert outside_span(vector, boundary_basis)
        first_vector, second_vector = vector & mask, vector >> size
        first_coset = remainder(first_vector, second_basis)
        second_coset = remainder(second_vector, first_basis)
        assert bool(first_coset) == bool(second_coset)
        if not first_coset:
            principal_count += 1
            continue
        entry = classes.setdefault(first_coset, {"second_coset": second_coset,
                                                 "first": first_vector, "second": second_vector,
                                                 "first_source": vector, "second_source": vector})
        assert entry["second_coset"] == second_coset
        if first_vector.bit_count() < entry["first"].bit_count():
            entry["first"] = first_vector
            entry["first_source"] = vector
        if second_vector.bit_count() < entry["second"].bit_count():
            entry["second"] = second_vector
            entry["second_source"] = vector
    best_label, best = min(classes.items(), key=lambda item: item[1]["first"].bit_count()
                          + item[1]["second"].bit_count())
    first_vector, second_vector = best["first"], best["second"]
    assert first_vector.bit_count() + second_vector.bit_count() == 6
    mismatch = apply_columns(first, first_vector) ^ apply_columns(second, second_vector)
    assert best_label and best["second_coset"]
    histogram = {}
    for entry in classes.values():
        weight = entry["first"].bit_count() + entry["second"].bit_count()
        histogram[str(weight)] = histogram.get(str(weight), 0) + 1
    return {"status": "PASS", "banked_population_size": len(objects),
            "principal_population_size": principal_count,
            "observed_nonzero_quotient_labels": len(classes),
            "separate_weight_sum_upper_bounds": histogram,
            "smallest_observed_separate_sum": first_vector.bit_count() + second_vector.bit_count(),
            "separate_block_weights": [first_vector.bit_count(), second_vector.bit_count()],
            "first_block_support": coordinates(first_vector, height),
            "second_block_support": coordinates(second_vector, height),
            "joint_syndrome_support": coordinates(mismatch, height),
            "first_source_cleared_support": list(positions(best["first_source"])),
            "second_source_cleared_support": list(positions(best["second_source"])),
            "scope": "Existing verified witnesses give upper bounds on the relaxed classical terms, not a distance census."}


if __name__ == "__main__":
    signal.signal(signal.SIGALRM, lambda signum, frame: sys.exit("55-second wall guard"))
    signal.alarm(55)
    resource.setrlimit(resource.RLIMIT_CPU, (50, 50))
    started = time.monotonic()
    result = verify()
    result["elapsed_seconds"] = time.monotonic() - started
    result["peak_rss_bytes"] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    assert result["peak_rss_bytes"] < 768 * 1024 * 1024
    print(json.dumps(result, indent=2))
