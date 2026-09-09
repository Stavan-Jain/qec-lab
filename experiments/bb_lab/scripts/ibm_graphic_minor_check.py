"""Find and verify a fixed Fano minor; no torus distance computation."""

import itertools
import json
import random
import resource
import signal
import sys
import time


ALPHA = ((0, 0), (0, -1), (-3, 1))
BETA = ((0, 0), (-1, 0), (1, 3))
CONTRACTED = (
    (0, 5, 3), (0, 3, 5), (0, 1, 2), (0, 6, 1), (1, 2, 3),
    (1, 5, 3), (1, 4, 5), (0, 0, 4), (1, 7, 6), (0, 3, 1),
    (0, 6, 4), (1, 3, 1), (1, 3, 3), (0, 6, 2), (1, 6, 3),
    (1, 2, 4), (1, 4, 3), (1, 6, 6), (0, 2, 3), (1, 4, 7),
    (0, 2, 4), (0, 4, 2), (1, 7, 7), (0, 5, 1), (1, 3, 4),
    (0, 1, 3), (1, 2, 2), (1, 6, 2), (0, 5, 0), (1, 7, 5),
    (1, 5, 1), (0, 2, 2), (0, 4, 0), (1, 4, 6), (1, 5, 5),
)
FANO = ((0, 0, 3), (0, 3, 4), (1, 6, 5), (0, 3, 2),
        (1, 3, 2), (0, 0, 5), (1, 4, 2))


def syndrome_column(block, horizontal, vertical):
    polynomial = BETA if block == 0 else ALPHA
    return {(horizontal + shift_x, vertical + shift_y)
            for shift_x, shift_y in polynomial}


def encode(points, row_index):
    return sum(1 << row_index[point] for point in points)


def coordinate_basis(columns):
    pivots = {}
    basis_indices = []
    for index, column in enumerate(columns):
        residual = column
        coefficients = 0
        while residual:
            pivot = residual.bit_length() - 1
            if pivot not in pivots:
                basis_bit = 1 << len(basis_indices)
                pivots[pivot] = residual, coefficients ^ basis_bit
                basis_indices.append(index)
                break
            reduced_column, reduced_coefficients = pivots[pivot]
            residual ^= reduced_column
            coefficients ^= reduced_coefficients
    coordinates = []
    for column in columns:
        residual = column
        coefficients = 0
        while residual:
            pivot = residual.bit_length() - 1
            reduced_column, reduced_coefficients = pivots[pivot]
            residual ^= reduced_column
            coefficients ^= reduced_coefficients
        reconstructed = 0
        for basis_position, basis_index in enumerate(basis_indices):
            if coefficients >> basis_position & 1:
                reconstructed ^= columns[basis_index]
        assert reconstructed == column
        coordinates.append(coefficients)
    return basis_indices, coordinates


def rank(columns):
    pivots = {}
    for column in columns:
        while column:
            pivot = column.bit_length() - 1
            if pivot not in pivots:
                pivots[pivot] = column
                break
            column ^= pivots[pivot]
    return len(pivots)


def main():
    signal.alarm(55)
    started = time.monotonic()
    if "--discover" not in sys.argv:
        cells = CONTRACTED + FANO
        assert len(set(cells)) == 42
        supports = [syndrome_column(*cell) for cell in cells]
        rows = sorted(set().union(*supports))
        row_index = {point: index for index, point in enumerate(rows)}
        columns = [encode(support, row_index) for support in supports]
        contracted_columns = columns[:len(CONTRACTED)]
        selected_columns = columns[len(CONTRACTED):]
        base_rank = rank(contracted_columns)
        assert base_rank == 35
        for mask in range(128):
            subset = [selected_columns[position] for position in range(7)
                      if mask >> position & 1]
            expected = [position + 1 for position in range(7)
                        if mask >> position & 1]
            assert rank(contracted_columns + subset) - base_rank == rank(expected)
        assert all(-3 <= horizontal <= 8 and -1 <= vertical <= 10
                   for horizontal, vertical in rows)
        peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        peak_bytes = peak if sys.platform == "darwin" else peak * 1024
        assert peak_bytes < 1024 ** 3
        print(json.dumps({"status": "fixed F7 certificate: all 128 subset ranks verified",
                          "contracted_rank": base_rank, "union_rank": rank(columns),
                          "qubits": len(cells), "valid_periods": "ell,m >= 12",
                          "elapsed_seconds": round(time.monotonic() - started, 4),
                          "peak_rss_bytes": peak_bytes}))
        return
    cells = [(block, horizontal, vertical)
             for block in (0, 1)
             for horizontal in range(8)
             for vertical in range(8)]
    random.Random(0).shuffle(cells)
    supports = [syndrome_column(*cell) for cell in cells]
    rows = sorted(set().union(*supports))
    row_index = {point: index for index, point in enumerate(rows)}
    columns = [encode(support, row_index) for support in supports]
    basis_indices, coordinates = coordinate_basis(columns)
    for retained in itertools.combinations(range(len(basis_indices)), 3):
        representatives = {}
        for column_index, coefficients in enumerate(coordinates):
            projected = sum(((coefficients >> position) & 1) << offset
                            for offset, position in enumerate(retained))
            if projected:
                representatives.setdefault(projected, column_index)
        if len(representatives) != 7:
            continue
        selected = [representatives[value] for value in range(1, 8)]
        fano_row_span = {sum(((functional & value).bit_count() % 2) << (value - 1)
                             for value in range(1, 8))
                         for functional in range(8)}
        contracted = []
        for position, index in enumerate(basis_indices):
            row_pattern = sum(((coordinates[selected[offset]] >> position) & 1) << offset
                              for offset in range(7))
            if row_pattern not in fano_row_span:
                contracted.append(index)
        contracted_columns = [columns[index] for index in contracted]
        base_rank = rank(contracted_columns)
        assert base_rank == len(contracted)
        for mask in range(128):
            subset = [columns[selected[position]] for position in range(7)
                      if mask >> position & 1]
            expected = [position + 1 for position in range(7)
                        if mask >> position & 1]
            assert rank(contracted_columns + subset) - base_rank == rank(expected)
        peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        peak_bytes = peak if sys.platform == "darwin" else peak * 1024
        assert peak_bytes < 1024 ** 3
        print(json.dumps({
            "status": "F7 minor verified on all 128 subsets",
            "syndrome_row_box": [min(point[0] for point in rows),
                                 max(point[0] for point in rows),
                                 min(point[1] for point in rows),
                                 max(point[1] for point in rows)],
            "columns": len(columns),
            "rank": len(basis_indices),
            "contracted_qubits": [cells[index] for index in contracted],
            "fano_qubits_in_binary_order": [cells[index] for index in selected],
            "elapsed_seconds": round(time.monotonic() - started, 4),
            "peak_rss_bytes": peak_bytes,
        }, indent=2))
        return
    print(json.dumps({"status": "no F7 found by the fixed-basis projection test",
                      "rank": len(basis_indices),
                      "elapsed_seconds": round(time.monotonic() - started, 4)}))
    raise SystemExit(1)


if __name__ == "__main__":
    main()
