"""Everything about a BB code that is cheap enough to compute on keystroke.

The web UI splits a code's report in two: this module holds the half
that costs milliseconds — `[[n, k, ·]]`, check and qubit weights, the
CSS commutation guard, and the solver *premises* — while `solver.py`
holds the half that costs seconds to hours (the distance).

Nothing here is new mathematics; it is the same `checks` / `codeparams`
/ `linalg` layer the sweep scripts use, packaged as one JSON-able
record so the front end never has to know the shapes.
"""

from __future__ import annotations

import os
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any

import numpy as np

from ..checks import bb_check_matrices, CheckMatrices
from ..codeparams import code_params
from ..group import AbelianGroup
from ..linalg import nullspace_f2, quotient_complement_basis
from ..poly import Poly

LAB_ROOT = Path(__file__).resolve().parents[3]  # experiments/bb_lab

# The report is recomputed as the user types, so it has to stay quick.
# Measured (Apple Silicon): |G| = 144 → 0.02 s, 576 → 0.6 s, 1152 → 4.8 s,
# 2000 → 27 s, dominated by `checks.circulant`'s |G|² Python loop. 1200 keeps
# the worst case around five seconds and still sits far above anything this
# program studies (the corpus tops out at n = 288). The cap is about
# responsiveness, not correctness.
MAX_GROUP_ORDER = 1200


class CodeInputError(ValueError):
    """A malformed group/polynomial spec, phrased for the end user."""


@dataclass(frozen=True, slots=True)
class WeightProfile:
    """Row and column weights of one check matrix."""

    row_min: int
    row_max: int
    col_min: int
    col_max: int

    @property
    def rows_uniform(self) -> bool:
        return self.row_min == self.row_max

    @property
    def cols_uniform(self) -> bool:
        return self.col_min == self.col_max


def _profile(H: np.ndarray) -> WeightProfile:
    rows = H.sum(axis=1)
    cols = H.sum(axis=0)
    return WeightProfile(
        row_min=int(rows.min()), row_max=int(rows.max()),
        col_min=int(cols.min()), col_max=int(cols.max()),
    )


def parse_orders(spec: str | list[int]) -> AbelianGroup:
    """Parse a cyclic-factor spec: '12x6', '12, 6', '[12, 6]', '6x6x3'.

    The BB literature writes ℓ × m; the lab's group layer is rank-generic,
    so the UI accepts any number of factors up to the six variable names
    `poly.py` knows (x, y, z, w, v, u).
    """
    if isinstance(spec, (list, tuple)):
        raw = [str(s) for s in spec]
    else:
        cleaned = str(spec).strip().strip("[]()")
        for sep in ("×", "*", "x", "X", ";"):
            cleaned = cleaned.replace(sep, ",")
        raw = [t for t in (p.strip() for p in cleaned.split(",")) if t]
    if not raw:
        raise CodeInputError("Give at least one cyclic factor, e.g. '12x6'.")
    if len(raw) > 6:
        raise CodeInputError(
            f"{len(raw)} cyclic factors, but the polynomial parser only "
            "names six variables (x, y, z, w, v, u)."
        )
    orders: list[int] = []
    for tok in raw:
        try:
            n = int(tok)
        except ValueError:
            raise CodeInputError(
                f"{tok!r} is not a whole number — a group spec looks like '12x6'."
            ) from None
        if n <= 0:
            raise CodeInputError(f"Cyclic order must be positive, got {n}.")
        orders.append(n)
    size = 1
    for n in orders:
        size *= n
    if size > MAX_GROUP_ORDER:
        raise CodeInputError(
            f"|G| = {size} exceeds the UI's cap of {MAX_GROUP_ORDER} "
            f"(n would be {2 * size} qubits). Drive bb_lab directly for that."
        )
    return AbelianGroup(tuple(orders))


def parse_poly(text: str, group: AbelianGroup, label: str) -> Poly:
    """Parse one polynomial, re-phrasing parser errors for the UI."""
    if not text or not text.strip():
        raise CodeInputError(f"Polynomial {label} is empty.")
    try:
        p = Poly.from_string(text, group)
    except ValueError as e:
        raise CodeInputError(f"{label}: {e}") from None
    if not p.support:
        raise CodeInputError(
            f"{label} reduces to 0 over F₂ (terms cancelled in pairs). "
            "A BB code needs both polynomials non-zero."
        )
    return p


@dataclass
class CodeReport:
    """The instant half of a code's report."""

    orders: list[int]
    group_label: str
    group_order: int
    A: str
    B: str
    A_weight: int
    B_weight: int
    n: int
    k: int
    rank_HX: int
    rank_HZ: int
    num_checks: int
    check_weight: int | None          # None when the blocks disagree
    qubit_degree: int | None
    x_profile: dict[str, int]
    z_profile: dict[str, int]
    css_commutes: bool
    rate: float
    premises: dict[str, Any] = field(default_factory=dict)
    known: dict[str, Any] | None = None
    warnings: list[str] = field(default_factory=list)


def analyse(
    orders_spec: str | list[int],
    A_text: str,
    B_text: str,
    *,
    lookup_corpus: bool = True,
) -> tuple[CodeReport, CheckMatrices]:
    """Full instant report plus the check matrices (for the solver)."""
    G = parse_orders(orders_spec)
    A = parse_poly(A_text, G, "A")
    B = parse_poly(B_text, G, "B")

    checks = bb_check_matrices(A, B)
    params = code_params(checks)

    xp, zp = _profile(checks.H_X), _profile(checks.H_Z)
    # For a BB code every check touches wt(A) + wt(B) qubits and every
    # qubit sits in that many checks; report the measured values and only
    # collapse to a single number when the matrices really are uniform.
    uniform = (
        xp.rows_uniform and zp.rows_uniform
        and xp.row_min == zp.row_min
    )
    deg_uniform = (
        xp.cols_uniform and zp.cols_uniform
        and xp.col_min == zp.col_min
    )

    prod = (checks.H_X @ checks.H_Z.T) % 2
    css_ok = not bool(prod.any())

    warnings: list[str] = []
    if params.k == 0:
        warnings.append(
            "k = 0: this code encodes nothing, so distance is undefined."
        )
    if not css_ok:
        warnings.append(
            "H_X · H_Zᵀ ≠ 0 — the CSS commutation guard failed. "
            "That should be impossible for an abelian group; treat any "
            "downstream number as suspect."
        )
    expected = A.weight() + B.weight()
    if uniform and xp.row_min != expected:
        warnings.append(
            f"check weight {xp.row_min} ≠ wt(A) + wt(B) = {expected}; "
            "the supports must overlap after the circulant build."
        )

    report = CodeReport(
        orders=list(G.orders),
        group_label=G.label(),
        group_order=G.cardinality,
        A=A.canonical_string(),
        B=B.canonical_string(),
        A_weight=A.weight(),
        B_weight=B.weight(),
        n=params.n,
        k=params.k,
        rank_HX=params.rank_HX,
        rank_HZ=params.rank_HZ,
        num_checks=2 * G.cardinality,
        check_weight=xp.row_min if uniform else None,
        qubit_degree=xp.col_min if deg_uniform else None,
        x_profile=asdict(xp),
        z_profile=asdict(zp),
        css_commutes=css_ok,
        rate=round(params.k / params.n, 4) if params.n else 0.0,
        premises=premises(checks),
        warnings=warnings,
    )
    if lookup_corpus:
        report.known = corpus_lookup(G, A, B)
    return report, checks


# ---------------------------------------------------------------- premises


def premises(checks: CheckMatrices) -> dict[str, Any]:
    """Machine-checked facts a solver flag may be allowed to assume.

    Keys here are referenced by `solver.FLAG_NOTES[...].requires`, so a
    new Tandem flag with a caller obligation is wired up by naming the
    premise it needs — no changes to the request path.

    `coset_parity_even` is the `-cost-step=2` premise: every H_X row has
    even weight (so weight parity is constant on cosets of rowspan H_X)
    **and** every logical-class representative is even (so the constant
    is 0 in every class). Together: all feasible costs are even. This is
    the same test `ladder_sweep.py` and `tandem_verify.py` run before
    passing the flag.
    """
    hx_even = not any(int(r.sum()) % 2 for r in checks.H_X)
    V = quotient_complement_basis(checks.H_X, nullspace_f2(checks.H_Z))
    classes_even = not any(int(v.sum()) % 2 for v in V)
    return {
        "coset_parity_even": {
            "holds": bool(hx_even and classes_even),
            "label": "All feasible logical weights are even",
            "detail": (
                f"H_X rows all even: {hx_even}; "
                f"logical-class reps all even: {classes_even}"
            ),
        }
    }


# ------------------------------------------------------------ corpus lookup


def corpus_db() -> Path | None:
    """Locate the (gitignored) BB corpus, if this checkout has one."""
    env = os.environ.get("BB_LAB_DB")
    if env:
        p = Path(env).expanduser()
        return p if p.exists() else None
    p = LAB_ROOT / "data" / "bb_instances.duckdb"
    return p if p.exists() else None


def corpus_lookup(G: AbelianGroup, A: Poly, B: Poly) -> dict[str, Any] | None:
    """Has this code already been settled? Returns the stored row, if any.

    Matches on the canonical orbit representative under Aut(G) ⋉ G plus
    block swap — the same key `cell_hunt_ingest.py` writes — so a code
    entered in any equivalent presentation still finds its row. Any
    failure (no DB, no duckdb, unreadable file) degrades to `None`; the
    lookup is a convenience, never a dependency.
    """
    db = corpus_db()
    if db is None:
        return None
    try:
        import duckdb

        from ..canonical import canonical_pair
        from ..store import canonical_hash

        canon = canonical_pair(A.support, B.support, G)
        cA = Poly(support=frozenset(canon.A_support), group=G).canonical_string()
        cB = Poly(support=frozenset(canon.B_support), group=G).canonical_string()
        iid = canonical_hash(G.label(), cA, cB)
        con = duckdb.connect(str(db), read_only=True)
        try:
            row = con.execute(
                "SELECT code_id, d_exact, d_lb, d_ub, d_method "
                "FROM bb_instances WHERE instance_id = ?",
                [iid],
            ).fetchone()
        finally:
            con.close()
        if row is None:
            return {"instance_id": iid, "found": False,
                    "orbit_size": canon.orbit_size}
        return {
            "instance_id": iid,
            "found": True,
            "orbit_size": canon.orbit_size,
            "code_id": row[0],
            "d_exact": row[1],
            "d_lb": row[2],
            "d_ub": row[3],
            "d_method": row[4],
        }
    except Exception:  # corpus is optional; never fail the report over it
        return None
