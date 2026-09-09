"""One gross rank check and one fixed norm-cover check; no distance search."""

import json
import resource
import signal
import sys
import time


def positions(vector):
    while vector:
        lowest = vector & -vector
        yield lowest.bit_length() - 1
        vector ^= lowest


def apply_columns(columns, vector):
    result = 0
    for position in positions(vector):
        result ^= columns[position]
    return result


def independent_basis(vectors):
    pivots = {}
    for vector in vectors:
        while vector:
            pivot = vector.bit_length() - 1
            if pivot not in pivots:
                pivots[pivot] = vector
                break
            vector ^= pivots[pivot]
    return pivots


def outside_span(vector, pivots):
    while vector:
        pivot = vector.bit_length() - 1
        if pivot not in pivots:
            return True
        vector ^= pivots[pivot]
    return False


def kernel_basis(columns):
    pivots = {}
    kernel = []
    for column, value in enumerate(columns):
        source = 1 << column
        while value:
            pivot = value.bit_length() - 1
            if pivot not in pivots:
                pivots[pivot] = (value, source)
                break
            old_value, old_source = pivots[pivot]
            value ^= old_value
            source ^= old_source
        if not value:
            kernel.append(source)
    assert len(kernel) + len(pivots) == len(columns)
    assert all(apply_columns(columns, vector) == 0 for vector in kernel)
    return kernel


def multiplication_columns(width, height, support):
    columns = []
    for horizontal in range(width):
        for vertical in range(height):
            image = 0
            for horizontal_step, vertical_step in support:
                target = ((horizontal + horizontal_step) % width) * height
                target += (vertical + vertical_step) % height
                image ^= 1 << target
            columns.append(image)
    return columns


def matrices(width, height):
    first = multiplication_columns(width, height, ((3, 0), (0, 1), (0, 2)))
    second = multiplication_columns(width, height, ((0, 0), (1, 3), (2, 3)))
    size = width * height
    boundary = [second[column] | (first[column] << size) for column in range(size)]
    differential = first + second
    assert all(apply_columns(differential, vector) == 0 for vector in boundary)
    return first, second, differential, boundary


def norm_to_y_cover(vector, width, height):
    base_size = width * height
    cover_height = 2 * height
    cover_size = width * cover_height
    result = 0
    for position in positions(vector):
        block, coordinate = divmod(position, base_size)
        horizontal, vertical = divmod(coordinate, height)
        first_position = block * cover_size + horizontal * cover_height + vertical
        second_position = first_position + height
        result ^= (1 << first_position) | (1 << second_position)
    assert result.bit_count() == 2 * vector.bit_count()
    return result


def verify():
    width, height = 12, 6
    size = width * height
    first, second, differential, boundary = matrices(width, height)
    boundary_basis = independent_basis(boundary)
    boundary_rank = len(boundary_basis)
    cycles = kernel_basis(differential)
    ann_first = kernel_basis(first)
    ann_second = kernel_basis(second)
    product = [apply_columns(first, vector) for vector in second]
    ann_product = kernel_basis(product)
    common_ann = kernel_basis([first[column] | (second[column] << size) for column in range(size)])
    first_only = ann_first
    second_only = [vector << size for vector in ann_second]
    hidden = [apply_columns(second, vector) for vector in ann_product]
    first_dimension = len(independent_basis(boundary + first_only)) - boundary_rank
    second_dimension = len(independent_basis(boundary + second_only)) - boundary_rank
    joint_dimension = len(independent_basis(boundary + first_only + second_only)) - boundary_rank
    hidden_dimension = len(independent_basis(boundary + hidden)) - boundary_rank
    hidden_quotient_dimension = len(ann_product) - len(independent_basis(ann_first + ann_second))
    assert boundary_rank == 66
    assert len(cycles) - boundary_rank == 12
    assert len(common_ann) == 6
    assert first_dimension == second_dimension == joint_dimension == hidden_dimension == 6
    assert hidden_quotient_dimension == 6
    assert len(independent_basis(boundary + first_only + second_only + hidden)) - boundary_rank == 6

    witness_rows = ((0, 1, 4, 5), (1, 2, 3, 4), (0, 2, 3, 5), (1, 2, 4, 5),
                    (0, 1, 3, 4), (0, 2, 3, 5), (0, 2), (0, 2))
    witness_parameter = sum(1 << (horizontal * height + vertical)
                            for horizontal, rows in enumerate(witness_rows) for vertical in rows)
    assert apply_columns(product, witness_parameter) == 0
    assert outside_span(apply_columns(second, witness_parameter), boundary_basis)
    transverse_second_coefficient = 0
    for position in positions(witness_parameter):
        horizontal, vertical = divmod(position, height)
        if horizontal & 2 == 2:
            transverse_second_coefficient ^= (1, 2, 3)[(horizontal + vertical) % 3]
    assert transverse_second_coefficient == 1
    witness_first = apply_columns(second, witness_parameter ^ 1)
    witness_second = first[0]
    witness = witness_first | (witness_second << size)
    assert witness_first and witness_second
    assert apply_columns(differential, witness) == 0
    assert outside_span(witness, boundary_basis)

    _, _, cover_differential, cover_boundary = matrices(width, 2 * height)
    cover_boundary_rank = len(independent_basis(cover_boundary))
    norm_cycles = [norm_to_y_cover(vector, width, height) for vector in cycles]
    norm_pure = [norm_to_y_cover(vector, width, height) for vector in first_only + second_only]
    assert all(apply_columns(cover_differential, vector) == 0 for vector in norm_cycles)
    norm_rank = len(independent_basis(cover_boundary + norm_cycles)) - cover_boundary_rank
    norm_pure_rank = len(independent_basis(cover_boundary + norm_pure)) - cover_boundary_rank
    assert cover_boundary_rank == 138
    assert norm_rank == 6
    assert norm_pure_rank == 0
    return {
        "status": "PASS",
        "fixed_base": [width, height],
        "fixed_norm_cover": [width, 2 * height],
        "homology_dimension": len(cycles) - boundary_rank,
        "annihilator_dimensions": {"P": len(ann_first), "Q": len(ann_second),
                                   "PQ": len(ann_product), "common": len(common_ann)},
        "pure_image_dimensions": [first_dimension, second_dimension, joint_dimension],
        "hidden_principal_dimension": hidden_dimension,
        "norm_homology_rank": norm_rank,
        "norm_on_pure_images_rank": norm_pure_rank,
        "both_principal_witness_parameter": [list(divmod(position, height))
                                              for position in positions(witness_parameter)],
        "both_principal_witness_weights": [witness_first.bit_count(), witness_second.bit_count()],
        "transverse_U_squared_coefficient": transverse_second_coefficient,
        "scope": "Fixed ranks and an exact nonboundary witness; no minimum-distance assertion.",
    }


def timeout_handler(signum, frame):
    raise TimeoutError("55-second wall limit")


def main():
    signal.signal(signal.SIGALRM, timeout_handler)
    signal.alarm(55)
    resource.setrlimit(resource.RLIMIT_CPU, (50, 50))
    started = time.perf_counter()
    report = verify()
    report["elapsed_seconds"] = time.perf_counter() - started
    peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    peak_bytes = peak if sys.platform == "darwin" else 1024 * peak
    if peak_bytes > 768 * 1024 * 1024:
        raise MemoryError("768-MiB RSS limit")
    report["peak_rss_bytes"] = peak_bytes
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
