"""F₂-valued polynomials over an abelian group (i.e. elements of F₂[G]).

A polynomial is represented by its support: a frozenset of group elements.
F₂ coefficients are 0/1, so support determines the polynomial.

The string format matches `pipeline/attempts/*/state.yaml`:

    'x^3 + y + y^2'       # grossA
    'y^3 + x + x^2'       # grossB
    '1 + x^2 + x^7'       # from the [[90,8,10]] instance

Convention: variables are `x, y, z, ...` for the 1st, 2nd, 3rd, ... cyclic
factor. Exponents reduce mod the corresponding order. The constant monomial
is spelled `1` (the empty product).

Input is normalised before parsing (`normalize_poly_text`): Unicode
superscripts become `^n` (`y⁶` → `y^6`), the product dots `·⋅×∗` become
`*`, and implicit products between single-letter factors gain the `*`
(`xy^6` → `x*y^6`, `x²y¹²` → `x^2*y^12`).  Each rule is unambiguous for
single-letter variables and is the identity on already-valid input;
unknown variables and malformed factors still raise, so nothing is
silently misparsed.
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from typing import Iterable

from .group import AbelianGroup


_VAR_NAMES = "xyzwvu"  # extend if the lab ever needs rank > 6
_TERM_RE = re.compile(
    r"\A\s*(?P<var>[a-z])(?:\s*\^\s*(?P<exp>\d+))?\s*\Z",
    re.IGNORECASE,
)
_CONST_RE = re.compile(r"\A\s*1\s*\Z")

_SUPERSCRIPTS = str.maketrans("⁰¹²³⁴⁵⁶⁷⁸⁹", "0123456789")
_SUP_RUN_RE = re.compile(r"[⁰¹²³⁴⁵⁶⁷⁸⁹]+")
_IMPLICIT_RE = re.compile(r"([a-zA-Z](?:\s*\^\s*\d+)?)\s*(?=[a-zA-Z])")


def normalize_poly_text(s: str) -> str:
    """Widen the accepted syntax without ambiguity (identity on valid input).

    Unicode superscripts → `^n`; product dots → `*`; implicit products
    between single-letter factors → explicit `*`.  Unknown variables and
    malformed factors are left for the strict parser to reject loudly.
    """
    s = _SUP_RUN_RE.sub(lambda m: "^" + m.group().translate(_SUPERSCRIPTS), s)
    for dot in "·⋅×∗":
        s = s.replace(dot, "*")
    return _IMPLICIT_RE.sub(r"\1*", s)


def _parse_term(term: str, group: AbelianGroup) -> tuple[int, ...]:
    """Parse a monomial: a `*`-separated product of single-variable powers.

    Examples:  '1', 'x', 'y^3', 'x*y', 'x^2*y^3'.
    The constant 1 is the empty product. Implicit products like 'xy'
    (without a `*`) are *not* accepted — raise to avoid silently
    misparsing `state.yaml` inputs.
    """
    term = term.strip()
    if _CONST_RE.match(term):
        return tuple(0 for _ in group.orders)
    factors = [f.strip() for f in term.split("*") if f.strip()]
    out = [0] * group.rank
    for f in factors:
        if _CONST_RE.match(f):
            continue
        m = _TERM_RE.match(f)
        if not m:
            raise ValueError(
                f"poly parser: cannot parse factor {f!r} (in term {term!r}). "
                "Factors must be single-variable powers like 'x', 'y^3', or '1'. "
                "Implicit products like 'xy' (no '*') are unsupported."
            )
        var = m.group("var").lower()
        if var not in _VAR_NAMES[: group.rank]:
            raise ValueError(
                f"poly parser: variable {var!r} out of range for group of rank "
                f"{group.rank} (allowed: {list(_VAR_NAMES[: group.rank])})"
            )
        axis = _VAR_NAMES.index(var)
        exp = int(m.group("exp") or "1")
        out[axis] = (out[axis] + exp) % group.orders[axis]
    return tuple(out)


@dataclass(frozen=True, slots=True)
class Poly:
    """An element of F₂[G] represented by its support set."""

    support: frozenset[tuple[int, ...]]
    group: AbelianGroup

    @classmethod
    def from_string(cls, s: str, group: AbelianGroup) -> "Poly":
        """Parse a polynomial like ``'x^3 + y + y^2'``.

        Empty terms ('+ +', leading/trailing '+') are tolerated and produce
        the zero polynomial.  Input is normalised first, so `xy^6`,
        `x²y¹²`, and `x·y` parse as their explicit-`*` forms.
        """
        terms = [t.strip() for t in normalize_poly_text(s).split("+")]
        # Two monomials at the same group element cancel (we're over F₂);
        # accumulate as a symmetric multiset, then take the parity.
        counts: dict[tuple[int, ...], int] = {}
        for t in terms:
            if not t:
                continue
            g = _parse_term(t, group)
            counts[g] = counts.get(g, 0) + 1
        support = frozenset(g for g, c in counts.items() if c % 2 == 1)
        return cls(support=support, group=group)

    @classmethod
    def from_support(
        cls, support: Iterable[tuple[int, ...]], group: AbelianGroup
    ) -> "Poly":
        return cls(
            support=frozenset(group.reduce(g) for g in support),
            group=group,
        )

    @classmethod
    def zero(cls, group: AbelianGroup) -> "Poly":
        return cls(support=frozenset(), group=group)

    def weight(self) -> int:
        return len(self.support)

    def coef(self, g: tuple[int, ...]) -> int:
        return 1 if self.group.reduce(g) in self.support else 0

    def __call__(self, g: tuple[int, ...]) -> int:
        return self.coef(g)

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, Poly):
            return NotImplemented
        return self.group == other.group and self.support == other.support

    def __hash__(self) -> int:
        return hash((self.group, self.support))

    def canonical_string(self) -> str:
        """Deterministic string form, ordered by axis then exponent.

        Two polys with the same support always produce the same string.
        Used as the canonical key for the corpus.
        """
        if not self.support:
            return "0"
        terms: list[str] = []
        for g in sorted(self.support):
            terms.append(_monomial_to_string(g))
        return " + ".join(terms)


def _monomial_to_string(g: tuple[int, ...]) -> str:
    parts: list[str] = []
    for axis, exp in enumerate(g):
        if exp == 0:
            continue
        var = _VAR_NAMES[axis]
        parts.append(var if exp == 1 else f"{var}^{exp}")
    if not parts:
        return "1"
    return "*".join(parts) if len(parts) > 1 else parts[0]
