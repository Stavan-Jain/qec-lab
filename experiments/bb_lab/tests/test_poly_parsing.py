"""Parser normalisation: Unicode superscripts, dots, implicit products.

The lab's own notes write polynomials as `xy⁶ + xy¹⁰ + x²y¹²`; pasting
them into the UI or CLI must parse identically to the strict
`x*y^6 + x*y^10 + x^2*y^12` form.  Malformed input must still fail
loudly (no silent misparse — the original v0 concern).
"""

import pytest

from bb_lab.group import AbelianGroup
from bb_lab.poly import Poly, normalize_poly_text


G = AbelianGroup((5, 15))


def canon(s: str) -> str:
    return Poly.from_string(s, G).canonical_string()


def test_unicode_and_implicit_forms_match_strict():
    strict = canon("x*y^6 + x*y^10 + x^2*y^12")
    assert canon("xy^6 + xy^10 + x^2y^12") == strict
    assert canon("xy⁶ + xy¹⁰ + x²y¹²") == strict
    assert canon("x·y⁶ + x·y¹⁰ + x²·y¹²") == strict


def test_normalizer_is_identity_on_valid_input():
    for s in ("1 + y + x", "x^3 + y + y^2", "y^4 + x^8*y^2 + x^13"):
        assert normalize_poly_text(s) == s


def test_two_digit_superscripts():
    assert normalize_poly_text("y¹⁰") == "y^10"
    assert canon("y¹⁰") == canon("y^10")


def test_malformed_still_fails_loudly():
    with pytest.raises(ValueError):
        Poly.from_string("ab^2 + x", G)      # unknown variables
    with pytest.raises(ValueError):
        Poly.from_string("x^ + y", G)        # dangling caret
