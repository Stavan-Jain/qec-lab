"""Small structural checks, not distance searches, for GROSS_SURFACE_PRODUCT.md."""

import itertools
import json
import resource
import signal
import sys
import time
from pathlib import Path

import numpy as np

LAB = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(LAB / "src"))

from bb_lab.tower import TowerCode, validate_banked


def bits(positions):
    return sum(1 << int(position) for position in positions)


def indices(vector):
    while vector:
        lowest = vector & -vector
        yield lowest.bit_length() - 1
        vector ^= lowest


def apply_columns(columns, vector):
    result = 0
    for position in indices(vector):
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


def transposed(columns, rows):
    result = [0] * rows
    for column, support in enumerate(columns):
        for row in indices(support):
            result[row] ^= 1 << column
    return result


def guard():
    peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    peak_bytes = peak if sys.platform == "darwin" else peak * 1024
    assert peak_bytes < 768 * 1024 ** 2
    return peak_bytes


def seed():
    code = TowerCode("gross-product-seed", (12, 6),
                     frozenset(((0, 0), (0, 1), (3, -1))),
                     frozenset(((0, 0), (1, 0), (-1, -3))))
    assert code.n == 144 and code.k == 12
    assert np.all(code.HX.sum(axis=1) == 6)
    assert np.all(code.HZ.sum(axis=1) == 6)
    assert np.all(code.HX.sum(axis=0) == 3)
    assert np.all(code.HZ.sum(axis=0) == 3)
    cells = (
        (0, 1, 0), (0, 1, 2), (0, 2, 1), (0, 2, 3), (0, 3, 1), (0, 3, 5),
        (1, 3, 1), (1, 3, 5), (1, 4, 0), (1, 4, 1), (1, 4, 2), (1, 4, 5),
    )
    witness_x = bits(block * 72 + code.G.index((horizontal, vertical))
                     for block, horizontal, vertical in cells)
    witness_z = bits((1 - block) * 72 + code.G.index((-horizontal % 12, -vertical % 6))
                     for block, horizontal, vertical in cells)
    rows_x = [bits(np.flatnonzero(row)) for row in code.HX]
    rows_z = [bits(np.flatnonzero(row)) for row in code.HZ]
    assert all((row & witness_x).bit_count() % 2 == 0 for row in rows_z)
    assert all((row & witness_z).bit_count() % 2 == 0 for row in rows_x)
    assert outside_span(witness_x, independent_basis(rows_x))
    assert outside_span(witness_z, independent_basis(rows_z))
    boundary = {1: [set(np.flatnonzero(code.HZ[:, column])) for column in range(144)],
                2: [set(np.flatnonzero(row)) for row in code.HX]}
    return boundary, witness_x, witness_z


def product(parameter, base_boundary):
    dimensions = ((72, 144, 72), (parameter - 1, parameter), (parameter, parameter - 1))
    factors = (
        base_boundary,
        {1: [{row for row in (column - 1, column) if 0 <= row < parameter - 1}
             for column in range(parameter)]},
        {1: [{column, column + 1} for column in range(parameter - 1)]},
    )
    cells = {degree: [] for degree in range(5)}
    for degrees in itertools.product(range(3), range(2), range(2)):
        spaces = [range(dimensions[factor][degree]) for factor, degree in enumerate(degrees)]
        for coordinates in itertools.product(*spaces):
            cells[sum(degrees)].append(degrees + coordinates)
    lookup = {degree: {cell: index for index, cell in enumerate(level)}
              for degree, level in cells.items()}
    boundary = {}
    for degree in range(1, 5):
        columns = []
        for cell in cells[degree]:
            support = 0
            for factor in range(3):
                if not cell[factor]:
                    continue
                for target in factors[factor][cell[factor]][cell[factor + 3]]:
                    image = list(cell)
                    image[factor] -= 1
                    image[factor + 3] = int(target)
                    support ^= 1 << lookup[degree - 1][tuple(image)]
            columns.append(support)
        boundary[degree] = columns
    return cells, lookup, boundary


def connected(checks, size):
    parent = list(range(size))

    def root(vertex):
        while vertex != parent[vertex]:
            parent[vertex] = parent[parent[vertex]]
            vertex = parent[vertex]
        return vertex

    for check in checks:
        support = list(indices(check))
        if support:
            representative = root(support[0])
            for vertex in support[1:]:
                parent[root(vertex)] = representative
    return len({root(vertex) for vertex in range(size)}) == 1


def check_member(parameter, base_boundary, base_x, base_z):
    cells, lookup, boundary = product(parameter, base_boundary)
    guard()
    for degree in (2, 3, 4):
        assert all(apply_columns(boundary[degree - 1], column) == 0
                   for column in boundary[degree])
    checks_x = boundary[3]
    checks_z = transposed(boundary[2], len(cells[1]))
    basis_x = independent_basis(checks_x)
    basis_z = independent_basis(checks_z)
    size = len(cells[2])
    assert size == 144 * (3 * parameter ** 2 - 3 * parameter + 1)
    assert size - len(basis_x) - len(basis_z) == 12
    weights = max(map(int.bit_count, checks_x)), max(map(int.bit_count, checks_z))
    assert max(weights) <= 8
    participation_x = transposed(checks_x, size)
    assert max(map(int.bit_count, participation_x)) <= 6
    assert max(map(int.bit_count, boundary[2])) <= 6
    assert connected(checks_x + checks_z, size)
    witness_x = bits(lookup[2][(1, 1, 0, coordinate, row, 0)]
                     for coordinate in indices(base_x) for row in range(parameter))
    witness_z = bits(lookup[2][(1, 1, 0, coordinate, 0, column)]
                     for coordinate in indices(base_z) for column in range(parameter))
    assert witness_x.bit_count() == witness_z.bit_count() == 12 * parameter
    assert apply_columns(boundary[2], witness_x) == 0
    assert all((check & witness_z).bit_count() % 2 == 0 for check in checks_x)
    assert outside_span(witness_x, basis_x)
    assert outside_span(witness_z, basis_z)
    if parameter == 1:
        assert checks_x == [bits(support) for support in base_boundary[2]]
        assert boundary[2] == [bits(support) for support in base_boundary[1]]
    guard()
    return {"t": parameter, "n": size, "k": 12, "x_rank": len(basis_x),
            "z_rank": len(basis_z), "check_weight_maxima": weights,
            "nontrivial_x_z_witness_weights": [12 * parameter, 12 * parameter],
            "connected": True}


def main():
    signal.alarm(55)
    started = time.monotonic()
    validate_banked(LAB / "data")
    base_boundary, base_x, base_z = seed()
    result = [check_member(parameter, base_boundary, base_x, base_z)
              for parameter in (1, 2, 3)]
    print(json.dumps({"members": result, "elapsed_seconds": round(time.monotonic() - started, 3),
                      "peak_rss_bytes": guard(),
                      "tier": "conventions/ranks/witnesses; universal distances proved in note"}, indent=2))


if __name__ == "__main__":
    main()
