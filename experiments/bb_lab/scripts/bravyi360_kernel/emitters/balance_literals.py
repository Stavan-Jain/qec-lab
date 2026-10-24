"""Render large natural literals as balanced kernel-reducible expressions."""
import pathlib
import re
import sys


def lit(value):
    if value.bit_length() <= 4096:
        return hex(value)
    shift = ((value.bit_length() + 4095) // 4096 // 2) * 4096
    return ('(Nat.lor ' + lit(value & ((1 << shift) - 1)) +
            ' (Nat.shiftLeft ' + lit(value >> shift) + ' ' + str(shift) + '))')


if __name__ == '__main__':
    for filename in sys.argv[1:]:
        source = pathlib.Path(filename)
        text = re.sub(r'0x[0-9a-f]{1025,}',
                      lambda match: lit(int(match[0], 16)), source.read_text())
        target = source.with_name(source.stem + 'Balanced.lean')
        target.write_text(text)
        print(target, len(text))
