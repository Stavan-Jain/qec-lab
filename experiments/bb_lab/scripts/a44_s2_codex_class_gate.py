"""The three index-two pullbacks of the single banked W4 vector."""

import json
import resource
import signal
import sys
import time

from a44_s2_codex_identities import motif


ALPHA = {(0, 0), (0, -1), (-3, 1)}
BETA = {(0, 0), (-1, 0), (1, 3)}


def source_pair():
    return tuple({((horizontal + vertical + 8 * index) % 24,
                   (vertical + 8 * index) % 24)
                  for horizontal, vertical in block for index in range(3)}
                 for block in motif(4))


def normalize(position, frame):
    horizontal, vertical = position
    width, height, shear = frame
    quotient, remainder = divmod(vertical, height)
    return (horizontal - quotient * shear) % width, remainder


def translate(poly, position, frame):
    horizontal, vertical = position
    result = set()
    for shift_x, shift_y in poly:
        result.symmetric_difference_update({normalize((horizontal + shift_x, vertical + shift_y), frame)})
    return result


def syndrome(pair, frame):
    result = set()
    for block, polynomial in zip(pair, (BETA, ALPHA)):
        for position in block:
            result ^= translate(polynomial, position, frame)
    return result


def vector_bits(pair, frame):
    width, height, shear = frame
    return sum(1 << (block * width * height + horizontal * height + vertical)
               for block, polynomial in enumerate(pair) for horizontal, vertical in polynomial)


def verify_cover(label, frame, source):
    width, height, shear = frame
    lifted = tuple({(horizontal, vertical) for horizontal in range(width) for vertical in range(height)
                    if (horizontal % 24, vertical % 24) in block} for block in source)
    assert sum(map(len, lifted)) == 108
    assert not syndrome(lifted, frame)
    rows = []
    basis = {}
    for horizontal in range(width):
        for vertical in range(height):
            row = vector_bits(tuple(translate(polynomial, (horizontal, vertical), frame)
                                    for polynomial in (ALPHA, BETA)), frame)
            rows.append(row)
            multiplier = 1 << (horizontal * height + vertical)
            while row:
                pivot = row.bit_length() - 1
                if pivot not in basis:
                    basis[pivot] = row, multiplier
                    break
                reduced, reduced_multiplier = basis[pivot]
                row ^= reduced
                multiplier ^= reduced_multiplier
    vector = vector_bits(lifted, frame)
    residual = vector
    multiplier = 0
    while residual:
        pivot = residual.bit_length() - 1
        if pivot not in basis:
            break
        reduced, reduced_multiplier = basis[pivot]
        residual ^= reduced
        multiplier ^= reduced_multiplier
    record = {"cover": label, "frame_width_height_shear": frame,
              "qubits": 2 * width * height, "pullback_weight": vector.bit_count(),
              "rank_HX": len(basis), "nontrivial": bool(residual)}
    if residual:
        dual = 1 << (residual.bit_length() - 1)
        for pivot in sorted(basis):
            row = basis[pivot][0]
            if (row & dual).bit_count() % 2:
                dual ^= 1 << pivot
        assert all((row & dual).bit_count() % 2 == 0 for row in rows)
        assert (vector & dual).bit_count() % 2 == 1
        record.update(certificate_kind="odd Z-dual", certificate_hex=hex(dual),
                      certificate_weight=dual.bit_count())
    else:
        reconstructed = 0
        for index, row in enumerate(rows):
            if multiplier >> index & 1:
                reconstructed ^= row
        assert reconstructed == vector
        record.update(certificate_kind="X-stabilizer multiplier", certificate_hex=hex(multiplier),
                      certificate_weight=multiplier.bit_count())
    return record


def main():
    signal.alarm(55)
    resource.setrlimit(resource.RLIMIT_CPU, (55, 55))
    started = time.monotonic()
    source = source_pair()
    assert sum(map(len, source)) == 54
    assert not syndrome(source, (24, 24, 0))
    records = [verify_cover(label, frame, source) for label, frame in (
        ("x-period doubled", (48, 24, 0)),
        ("y-period doubled", (24, 48, 0)),
        ("diagonal cover", (48, 24, 24)),
    )]
    peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    peak_bytes = peak if sys.platform == "darwin" else peak * 1024
    assert peak_bytes <= 768 * 1024 ** 2
    print(json.dumps({"source": "W4 on original (24,24)", "source_weight": 54,
                      "covers": records, "elapsed_seconds": round(time.monotonic() - started, 4),
                      "peak_rss_bytes": peak_bytes}, indent=2))


if __name__ == "__main__":
    main()
