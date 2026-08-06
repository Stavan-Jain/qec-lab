"""General (ε,δ) fiber-sweep engine for BB codes — the A22/A23 technique
as a portable library (A28).

Ancestry: the f2a6f17e safe-floor deletion (A22 census fibering, A23 site
sweep), the per-code ports a19_fibering.py / a20_fibering.py, and the A27
feasibility probe. This module replaces all of those hand-built frames with
one code-agnostic engine. See ``notes/A28_general_fibering.md`` for the
theory write-up and the validation record.

Mathematical content (pure F₂ throughout — no field arithmetic):

* **Fiber frame.** For any z ∈ G of odd order q ≥ 3, G partitions into
  |G|/q fibers (cosets of ⟨z⟩). Fixing a section, every 1-chain block
  splits per site into a fiber vector t ∈ F₂^q.

* **Residue coordinates.** D(t)_j := t_j + t_{q−1} (j < q−1) has kernel
  exactly {0, N} (N = all-ones), and (parity, D) : F₂^q ≅ F₂ × F₂^{q−1}.
  The two lifts of a residue r weigh |r| and q − |r| (choose by parity):
  the *exact* weight table of A22 in closed form, valid for every odd q —
  the CRT/GF-field structure is scaffolding, not engine input.
  Consequences: w0(r) := min(|r|, q−|r|) ≥ 1 for r ≠ 0, and the parity
  gap is |q − 2|r|| ≥ 1 (q odd).

* **Cost lower bounds.** |u| + |v| = Σ_sites (fiber weights). Unpaired:
  every δ-active position (site, block) costs ≥ w0 ≥ 1. Paired: when the
  ε-quotient polynomials satisfy B_ε = x̄^m · A_ε (site-monomial link),
  the fiber parities of A⋆f at s and B⋆f at s+m agree up to a fixed
  offset twist τ(s); a parity-matched active site pair costs ≥ 2
  (τ = 0) because the even lift of a nonzero residue weighs ≥ 2 and the
  odd lift of the zero residue weighs q. This is the general form of
  A23's V3 table; the ≥2-per-site bound is a *theorem exactly in the
  linked case* — unlinked codes get the weaker unpaired bound (this
  corrects A27 §3.3 P5, whose ≤9-active estimate presumed the link).

* **The sweep.** A weight-< t element of an affine family
  (u, v) = (u₀ + A⋆f, v₀ + B⋆f) has ≤ maxact active units (sites when
  paired, positions when unpaired). For every maximal admissible mask,
  "off-mask fibers are ε-only" is an affine F₂ system in f; enumerate
  its solutions modulo the cost-invariant kernel K₀ = {f : A⋆f and B⋆f
  both ε-only} and check costs. Completeness: any violator's active set
  is contained in some maximal mask.

* **ε-exact costs.** For a fixed δ-config, cost(ε-choice) is *linear* in
  the ε-bits (per position: W(r,0) + (W(r,1)−W(r,0))·ε), and the
  reachable ε-patterns form base + {(A_ε ḡ, B_ε ḡ)}: the exact minimum
  is a minimum of a linear functional over an affine F₂-space — computed
  by enumeration when the ε-quotient is small. Relaxing ε (Φ ≥ 0) is
  always sound for floors.

* **Safe classes.** Under (R), safe classes are Δ(ker ∂₂ ∖ 0), coset
  minima are constant on G-orbits, and δ₂[ζ] = [seamC ζ] with seamC the
  cover-window carry chain (Prop A14.1). The safe-floor query feeds
  seamC orbit representatives as offsets (u₀|v₀).

Conventions match bb_lab.checks / the Lean side: conv A f g =
Σ_h A(h) f(g−h); M_A[g,h] = A(g−h); ∂₂ f = (A⋆f | B⋆f) with the A-block
first. All chains are numpy uint8; masks and linear algebra run on
int-packed bitrows for speed.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from math import comb, gcd
from typing import Iterable, Optional

import numpy as np

from .checks import circulant
from .group import AbelianGroup
from .poly import Poly

__all__ = [
    "FiberFrame",
    "FiberScore",
    "SweepResult",
    "SafeFloorReport",
    "enumerate_fiber_generators",
    "score_fiber",
    "best_frames",
    "kernel_basis",
    "kernel_orbit_reps",
    "seam_offsets",
    "safe_floor_certify",
]


# --------------------------------------------------------------------------
# int-bitrow F₂ linear algebra (fast path for the mask sweep)
# --------------------------------------------------------------------------


def _rref_int(rows: list[int]) -> tuple[list[int], dict[int, int]]:
    """RREF of int-bitmask rows. Returns (reduced rows, {pivot col: row idx})."""
    red: list[int] = []
    pivots: dict[int, int] = {}
    for r in rows:
        cur = r
        for c, ri in pivots.items():
            if (cur >> c) & 1:
                cur ^= red[ri]
        if cur:
            c = (cur & -cur).bit_length() - 1
            for i, rr in enumerate(red):
                if (rr >> c) & 1:
                    red[i] = rr ^ cur
            pivots[c] = len(red)
            red.append(cur)
    return red, pivots


def _solve_affine_int(
    rows: list[int], consts: list[int], ncols: int
) -> Optional[tuple[int, list[int]]]:
    """Solve M x = c over F₂ (rows = int bitmasks over ncols, consts 0/1).

    Returns None if inconsistent, else (particular solution bitmask,
    kernel basis bitmasks).
    """
    aug = [r | (c << ncols) for r, c in zip(rows, consts)]
    return _solve_affine_aug(aug, ncols)


def _solve_affine_aug(
    aug: list[int], ncols: int
) -> Optional[tuple[int, list[int]]]:
    """Like _solve_affine_int but rows carry the constant at bit ncols."""
    red, pivots = _rref_int(aug)
    if ncols in pivots:
        return None  # row 0 = 1
    part = 0
    piv_cols = set(pivots)
    for c, ri in pivots.items():
        if (red[ri] >> ncols) & 1:
            part |= 1 << c
    kernel: list[int] = []
    for free in range(ncols):
        if free in piv_cols:
            continue
        v = 1 << free
        for c, ri in pivots.items():
            if (red[ri] >> free) & 1:
                v |= 1 << c
        kernel.append(v)
    return part, kernel


def _int_to_vec(x: int, n: int) -> np.ndarray:
    return np.frombuffer(
        np.binary_repr(x, width=n)[::-1].encode(), dtype=np.uint8
    ) - ord("0")


def _vec_to_int(v: np.ndarray) -> int:
    out = 0
    for i in np.flatnonzero(v):
        out |= 1 << int(i)
    return out


# --------------------------------------------------------------------------
# fiber frame
# --------------------------------------------------------------------------


@dataclass
class FiberFrame:
    """The (ε,δ) fiber decomposition of F₂[G] for one fiber generator z."""

    A: Poly
    B: Poly
    z: tuple[int, ...]

    def __post_init__(self) -> None:  # noqa: PLR0915
        G = self.A.group
        if self.B.group != G:
            raise ValueError("A and B must share a group")
        self.group = G
        n = G.cardinality
        self.n = n
        # fiber subgroup
        q = 1
        g = G.reduce(self.z)
        zero = tuple(0 for _ in G.orders)
        cur = g
        while cur != zero:
            q += 1
            cur = G.add(cur, g)
        if q % 2 == 0 or q < 3:
            raise ValueError(f"fiber order must be odd ≥ 3, got {q}")
        self.q = q
        elems = list(G)
        self.elems = elems
        eidx = {e: i for i, e in enumerate(elems)}
        # sites = cosets of <z>; canonical rep = min element in iteration order
        site_of = [-1] * n
        fibers: list[list[int]] = []
        site_rep: list[tuple[int, ...]] = []
        for i, e in enumerate(elems):
            if site_of[i] >= 0:
                continue
            s = len(fibers)
            cells = []
            cur = e
            for _ in range(q):
                ci = eidx[cur]
                site_of[ci] = s
                cells.append(ci)
                cur = G.add(cur, g)
            fibers.append(cells)  # ordered by fiber coordinate k: rep + k z
            site_rep.append(e)
        self.S = len(fibers)
        self.FIBERS = np.array(fibers)  # (S, q) cell indices
        self.site_of = np.array(site_of)
        self.site_rep = site_rep
        # quotient group as Cayley table on site indices
        rep_idx = {e: s for s, e in enumerate(site_rep)}
        self.site_add = np.array(
            [
                [site_of[eidx[G.add(site_rep[a], site_rep[b])]] for b in range(self.S)]
                for a in range(self.S)
            ]
        )
        self.site_neg = np.array(
            [site_of[eidx[G.neg(site_rep[a])]] for a in range(self.S)]
        )
        del rep_idx
        # operators
        self.MA = circulant(self.A)
        self.MB = circulant(self.B)
        # ε-quotient polynomials over the site set (multiset parity)
        self.A_eps = self._eps_poly(self.A)
        self.B_eps = self._eps_poly(self.B)
        # δ-condition rows: for (block, site): (q−1) rows over f-columns,
        # row j = M[cells[j]] ^ M[cells[q−1]]  (residue of that fiber of M f)
        self._cond_rows: dict[tuple[int, int], list[int]] = {}
        MA_int = [_vec_to_int(r) for r in self.MA]
        MB_int = [_vec_to_int(r) for r in self.MB]
        for blk, Mint in ((0, MA_int), (1, MB_int)):
            for s in range(self.S):
                cells = fibers[s]
                last = Mint[cells[q - 1]]
                self._cond_rows[(blk, s)] = [
                    Mint[cells[j]] ^ last for j in range(q - 1)
                ]
        # cost-invariant kernel K₀ = {f : all residues of Af and Bf vanish}
        all_rows: list[int] = []
        for key in self._cond_rows:
            all_rows.extend(self._cond_rows[key])
        sol = _solve_affine_int(all_rows, [0] * len(all_rows), n)
        assert sol is not None
        _, self.K0 = sol
        self.r_delta = n - len(self.K0)  # rank of the δ-data map
        # ε-link detection: B_eps == x̄^m · A_eps for some site-monomial m?
        self.link_shift: Optional[int] = None
        for m in range(self.S):
            if self._site_shift(self.B_eps) == self._site_shift(
                self._mono_mul(self.A_eps, m)
            ):
                self.link_shift = m
                break

    # -- small helpers -----------------------------------------------------

    def _eps_poly(self, P: Poly) -> frozenset[int]:
        acc: set[int] = set()
        eidx = {e: i for i, e in enumerate(self.elems)}
        for mono in P.support:
            s = int(self.site_of[eidx[self.group.reduce(mono)]])
            if s in acc:
                acc.remove(s)
            else:
                acc.add(s)
        return frozenset(acc)

    def _mono_mul(self, P: frozenset[int], m: int) -> frozenset[int]:
        return frozenset(int(self.site_add[p, m]) for p in P)

    @staticmethod
    def _site_shift(P: frozenset[int]) -> frozenset[int]:
        return P

    # -- chain-level maps --------------------------------------------------

    def fiber_weights(self, v: np.ndarray) -> np.ndarray:
        return v[self.FIBERS].sum(axis=1)

    def eps_of(self, v: np.ndarray) -> np.ndarray:
        return self.fiber_weights(v) % 2

    def residues(self, v: np.ndarray) -> np.ndarray:
        """(S, q−1) residue bits of each fiber of v."""
        fib = v[self.FIBERS]
        return (fib[:, : self.q - 1] ^ fib[:, self.q - 1 :][:, [0] * (self.q - 1)])

    def residue_pop(self, v: np.ndarray) -> np.ndarray:
        return self.residues(v).sum(axis=1)

    def w0_of_pop(self, pop: np.ndarray) -> np.ndarray:
        return np.minimum(pop, self.q - pop)

    def W(self, pop: np.ndarray, eps: np.ndarray) -> np.ndarray:
        """Exact fiber weight from residue popcount + parity bit."""
        even_is_pop = (pop % 2) == (eps % 2)
        return np.where(even_is_pop, pop, self.q - pop)

    def chain_weight_check(self, v: np.ndarray) -> None:
        """Assert Σ W(residue, ε) == |v| — the exact weight formula."""
        pop = self.residue_pop(v)
        eps = self.eps_of(v)
        assert int(self.W(pop, eps).sum()) == int(v.sum())

    # -- ε-quotient operators as site matrices -----------------------------

    def eps_op(self, P: frozenset[int]) -> np.ndarray:
        """S×S matrix of convolution by the ε-image P on F₂[sites]."""
        M = np.zeros((self.S, self.S), dtype=np.uint8)
        for s in range(self.S):
            for p in P:
                M[self.site_add[p, s], s] ^= 1
        return M


# --------------------------------------------------------------------------
# feasibility scoring
# --------------------------------------------------------------------------


@dataclass
class FiberScore:
    z: tuple[int, ...]
    q: int
    S: int
    budget: int
    mode: str  # "paired" | "unpaired" | "vacuous"
    maxact: int
    n_masks: int
    r_delta: int
    conds_at_max: int
    margin: int
    link_shift: Optional[int]
    feasible: bool
    note: str = ""

    def sort_key(self) -> tuple:
        return (not self.feasible, self.n_masks)


MASK_CAP = 2_000_000


def enumerate_fiber_generators(G: AbelianGroup) -> list[tuple[int, ...]]:
    """One generator per cyclic subgroup of odd order ≥ 3."""
    seen: set[frozenset[tuple[int, ...]]] = set()
    out = []
    zero = tuple(0 for _ in G.orders)
    for g in G:
        if g == zero:
            continue
        # order of g
        q = 1
        cur = g
        while cur != zero:
            cur = G.add(cur, g)
            q += 1
        if q % 2 == 0 or q < 3:
            continue
        sub = set()
        cur = zero
        for _ in range(q):
            sub.add(cur)
            cur = G.add(cur, g)
        key = frozenset(sub)
        if key in seen:
            continue
        seen.add(key)
        out.append(g)
    return out


def score_fiber(
    A: Poly, B: Poly, z: tuple[int, ...], target: int, *, parity_even: bool = True
) -> FiberScore:
    """Cheap feasibility triage for one (code, fiber, floor-target)."""
    fr = FiberFrame(A, B, z)
    budget = target - 2 if parity_even else target - 1
    if fr.link_shift is not None:
        maxact = budget // 2
        if maxact >= fr.S:
            return FiberScore(
                z, fr.q, fr.S, budget, "vacuous", maxact, 0, fr.r_delta, 0, 0,
                fr.link_shift, False, "paired maxact >= S",
            )
        n_masks = comb(fr.S, maxact)
        conds = 2 * (fr.S - maxact) * (fr.q - 1)
        margin = conds - fr.r_delta
        return FiberScore(
            z, fr.q, fr.S, budget, "paired", maxact, n_masks, fr.r_delta,
            conds, margin, fr.link_shift,
            n_masks <= MASK_CAP and margin > -24,
            "",
        )
    # unpaired: positions = 2S, per-position ≥ 1; the sweep is a
    # consistency-pruned DFS, so the raw mask count is an upper bound,
    # not a feasibility verdict — margin is the real gate.
    maxact = budget
    if maxact >= 2 * fr.S:
        return FiberScore(
            z, fr.q, fr.S, budget, "vacuous", maxact, 0, fr.r_delta, 0, 0,
            None, False, "unpaired maxact >= 2S",
        )
    n_masks = comb(2 * fr.S, maxact)
    conds = (2 * fr.S - maxact) * (fr.q - 1)
    margin = conds - fr.r_delta
    return FiberScore(
        z, fr.q, fr.S, budget, "unpaired", maxact, n_masks, fr.r_delta,
        conds, margin, None,
        margin > -24,
        "dfs-pruned" if n_masks > MASK_CAP else "",
    )


def best_frames(
    A: Poly, B: Poly, target: int, *, parity_even: bool = True
) -> list[FiberScore]:
    scores = [
        score_fiber(A, B, z, target, parity_even=parity_even)
        for z in enumerate_fiber_generators(A.group)
    ]
    return sorted(scores, key=FiberScore.sort_key)


# --------------------------------------------------------------------------
# the sweep
# --------------------------------------------------------------------------


@dataclass
class SweepResult:
    mode: str
    target: int
    budget: int
    maxact: int
    n_masks: int
    n_incons: int
    n_cons: int
    kdims: dict[int, int]
    n_reps: int
    min_relaxed: int
    n_suspects: int
    min_exact: Optional[int]
    violations: list[dict]
    tight: int
    certified: bool
    aborted: str = ""
    undecided: int = 0

    def summary(self) -> str:
        s = (
            f"[{self.mode}] target {self.target} budget {self.budget} "
            f"maxact {self.maxact}: {self.n_cons} cons / {self.n_incons} incons "
            f"of {self.n_masks} masks; kdims {self.kdims}; reps {self.n_reps}; "
            f"min_relaxed {self.min_relaxed}"
        )
        if self.n_suspects:
            s += f"; suspects {self.n_suspects} min_exact {self.min_exact}"
        s += (
            f"; violations {len(self.violations)}; tight {self.tight}; "
            f"{'CERTIFIED' if self.certified else 'NOT-CERTIFIED'}"
        )
        if self.aborted:
            s += f" [ABORTED: {self.aborted}]"
        return s


KEXTRA_CAP = 16
EPS_ENUM_CAP = 22
REPS_CAP = 2_000_000
NODE_CAP = 20_000_000


class FiberSweep:
    """Floor / census sweeps for one frame and one affine family."""

    def __init__(self, frame: FiberFrame):
        self.fr = frame
        fr = frame
        # ε-image spaces for the exact ε-minimization: rows of the joint
        # map ḡ ↦ (A_ε ḡ | B_ε ḡ) on F₂[sites] (2S bits per generator)
        MAe = fr.eps_op(fr.A_eps)
        MBe = fr.eps_op(fr.B_eps)
        gens = []
        for j in range(fr.S):
            vec = np.concatenate([MAe[:, j], MBe[:, j]])
            gens.append(_vec_to_int(vec))
        red, piv = _rref_int(gens)
        self.eps_image_basis = red  # independent rows (ints over 2S bits)

    # -- exact ε-min for one δ-config --------------------------------------

    FLIPS_DIM_CAP = 20

    def _flips_matrix(self) -> Optional[np.ndarray]:
        """(2^dim, 2S) bool matrix of all ε-image flip patterns (cached)."""
        if hasattr(self, "_FLIPS"):
            return self._FLIPS
        basis = self.eps_image_basis
        dim = len(basis)
        npos = 2 * self.fr.S
        if dim > self.FLIPS_DIM_CAP:
            self._FLIPS = None
            return None
        bas = np.zeros((max(dim, 1), npos), dtype=bool)
        for i, b in enumerate(basis):
            bas[i] = _int_to_vec(b, npos).astype(bool)
        F = np.zeros((1 << dim, npos), dtype=bool)
        for x in range(1, 1 << dim):
            lsb = (x & -x).bit_length() - 1
            F[x] = F[x ^ (1 << lsb)] ^ bas[lsb]
        self._FLIPS = F
        return F

    def _config_cost_vector(
        self, u: np.ndarray, v: np.ndarray
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """(costs over all flips, Wat0, eps0) for one δ-config."""
        fr = self.fr
        pops = np.concatenate([fr.residue_pop(u), fr.residue_pop(v)])
        eps0 = np.concatenate([fr.eps_of(u), fr.eps_of(v)])
        Wat0 = fr.W(pops, eps0).astype(np.int64)
        Wat1 = fr.W(pops, eps0 ^ 1).astype(np.int64)
        F = self._flips_matrix()
        assert F is not None
        costs = int(Wat0.sum()) + F @ (Wat1 - Wat0)
        return costs, pops, eps0

    def exact_cost(
        self, u: np.ndarray, v: np.ndarray, *, exclude_zero: bool = False
    ) -> tuple[int, bool, Optional[np.ndarray]]:
        """min |u'| + |v'| over the ε-fiber of (u, v)'s δ-config.

        Returns (value, is_exact, best_flip_pattern). When the ε-image is
        too large to enumerate, falls back to the relaxed (ε-free) bound
        with is_exact=False — sound for floors, unusable for refutations.
        With exclude_zero, the zero chain (only reachable when the
        δ-config vanishes) is skipped — census/boundary-floor semantics.
        """
        fr = self.fr
        F = self._flips_matrix()
        if F is None:
            pops = np.concatenate([fr.residue_pop(u), fr.residue_pop(v)])
            return int(fr.w0_of_pop(pops).sum()), False, None
        costs, pops, eps0 = self._config_cost_vector(u, v)
        if exclude_zero and int(pops.sum()) == 0:
            # the zero chain is the flip pattern equal to eps0 (if present)
            zero_rows = np.flatnonzero((F == eps0.astype(bool)).all(axis=1))
            if zero_rows.size:
                costs = costs.copy()
                costs[zero_rows] = 10 ** 9
        idx = int(np.argmin(costs))
        best = int(costs[idx])
        if best >= 10 ** 9:
            return 10 ** 9, True, None
        return best, True, F[idx].astype(np.uint8)

    def dress_witness(
        self, f: np.ndarray, offsets: tuple[np.ndarray, np.ndarray],
        flip: np.ndarray,
    ) -> Optional[tuple[np.ndarray, np.ndarray, np.ndarray]]:
        """Realize an ε-flip pattern as a family member.

        Solves (A_ε ḡ | B_ε ḡ) = flip for a site function ḡ, lifts it to
        κ = Σ ḡ_s N_s ∈ K₀, and returns (f + κ, u', v'). None if the flip
        is not in the ε-image (cannot happen for exact_cost outputs).
        """
        fr = self.fr
        if not hasattr(self, "_joint_rows"):
            MAe = fr.eps_op(fr.A_eps)
            MBe = fr.eps_op(fr.B_eps)
            joint = np.vstack([MAe, MBe])  # (2S, S)
            self._joint_rows = [
                _vec_to_int(joint[i]) for i in range(2 * fr.S)
            ]
        sol = _solve_affine_int(
            self._joint_rows, [int(b) for b in flip], fr.S
        )
        if sol is None:
            return None
        gbar = _int_to_vec(sol[0], fr.S)
        kappa = np.zeros(fr.n, dtype=np.uint8)
        for s in np.flatnonzero(gbar):
            kappa[fr.FIBERS[int(s)]] ^= 1
        f2 = (f + kappa) % 2
        u0, v0 = offsets
        u2 = (fr.MA @ f2 + u0) % 2
        v2 = (fr.MB @ f2 + v0) % 2
        return f2, u2, v2

    # -- mask machinery ----------------------------------------------------

    def _mask_universe(self, mode: str) -> list[tuple[int, int]]:
        """Positions: (block, site) pairs; paired mode groups them."""
        fr = self.fr
        if mode == "paired":
            m = fr.link_shift
            assert m is not None
            # site-unit s covers u-position (0, s) and v-position (1, s+m)
            return [(s, int(fr.site_add[s, m])) for s in range(fr.S)]
        return [(blk, s) for blk in (0, 1) for s in range(fr.S)]

    def _unit_block_rows(
        self, mode: str, unit, offs: tuple[np.ndarray, np.ndarray]
    ) -> dict[str, tuple[list[int], list[int]]]:
        """Residue-vanishing conditions per block of one unit.

        Keys: "u"/"v" (paired); unpaired units expose their single block
        as "u"."""
        fr = self.fr
        u0, v0 = offs

        def block(blk: int, s: int, off: np.ndarray):
            cells = fr.FIBERS[s]
            last = int(off[cells[fr.q - 1]])
            rows = list(fr._cond_rows[(blk, s)])
            consts = [int(off[cells[j]]) ^ last for j in range(fr.q - 1)]
            return rows, consts

        if mode == "paired":
            su, sv = unit
            return {"u": block(0, su, u0), "v": block(1, sv, v0)}
        blk, s = unit
        return {"u": block(blk, s, u0 if blk == 0 else v0)}

    # -- the main sweep ----------------------------------------------------

    def floor_sweep(  # noqa: PLR0912, PLR0915
        self,
        offsets: tuple[np.ndarray, np.ndarray],
        target: int,
        *,
        mode: Optional[str] = None,
        collect_survivors: bool = False,
        progress: bool = False,
        exclude_zero: Optional[bool] = None,
        exact_pass: bool = True,
        reps_cap: int = REPS_CAP,
    ) -> SweepResult:
        """Certify min_{f} |u₀+Af| + |v₀+Bf| ≥ target, or find violations.

        offsets = (u₀, v₀) as length-n uint8 vectors (zero for census use).
        exclude_zero defaults to True exactly when both offsets vanish
        (nonzero-boundary floor semantics).
        """
        fr = self.fr
        n = fr.n
        u0, v0 = offsets
        if exclude_zero is None:
            exclude_zero = not (u0.any() or v0.any())
        # family parity: |u|+|v| ≡ |u₀|+|v₀| (mod 2) since |A|,|B| odd ⇒
        # |Af|+|Bf| ≡ 0. Guard: require odd polynomial weights.
        assert len(fr.A.support) % 2 == 1 and len(fr.B.support) % 2 == 1
        par = (int(u0.sum()) + int(v0.sum())) % 2
        budget = target - 1
        if (budget - par) % 2 == 1:
            budget -= 1  # parity: weights ≡ par (mod 2)
        if mode is None:
            mode = "paired" if fr.link_shift is not None else "unpaired"
        if mode == "paired" and fr.link_shift is None:
            raise ValueError("paired mode requires the ε-monomial link")
        units = self._mask_universe(mode)
        # Per-unit state lattice (state = (min cost, which blocks' residue
        # conditions apply)). Paired τ=0 sites: OFF(0, both) < ON(2, none).
        # Paired τ=1 sites (offset ε-twist): OFF(q, both) vs
        # U0(1, u-side) / V0(1, v-side) < ON(3, none) — τ1 site costs are
        # odd, and a cost-1 site has one block residue-free (the general
        # correct accounting; see A28 note). Unpaired positions:
        # OFF(0, self) < ON(1, none).
        tau: list[int] = [0] * len(units)
        if mode == "paired":
            eu, ev = fr.eps_of(u0), fr.eps_of(v0)
            for i, (su_, sv_) in enumerate(units):
                tau[i] = int(eu[su_]) ^ int(ev[sv_])
        # states per unit: list of (cost, cond_key), condition-strength
        # DEcreasing (last = weakest = "active"). τ1-OFF (cost q, both
        # conds) is dominated by U0/V0 (cost 1, fewer conds) and dropped.
        q = fr.q
        unit_states: list[list[tuple[int, str]]] = []
        for i in range(len(units)):
            if mode == "unpaired":
                unit_states.append([(0, "u"), (1, "")])
            elif tau[i] == 0:
                unit_states.append([(0, "uv"), (2, "")])
            else:
                unit_states.append([(1, "u"), (1, "v"), (3, "")])
        min_unit_cost = [min(c for c, _ in st) for st in unit_states]
        total_min = sum(min_unit_cost)
        maxact = budget  # informational
        if total_min > budget:
            # every family member exceeds the budget: certified trivially
            return SweepResult(
                mode, target, budget, maxact, 0, 0, 0, {}, 0, 10 ** 9, 0,
                None, [], 0, True,
            )
        suffix_min = [0] * (len(units) + 1)
        for i in range(len(units) - 1, -1, -1):
            suffix_min[i] = suffix_min[i + 1] + min_unit_cost[i]
        # "active" cost per unit (the weakest, condition-free state)
        active_cost = [st[-1][0] for st in unit_states]

        # precompute per-unit per-block rows/consts (augmented ints: bit n
        # is the affine constant)
        unit_rows = [
            self._unit_block_rows(mode, units[i], offsets)
            for i in range(len(units))
        ]
        unit_aug: list[dict[str, list[int]]] = []
        for i in range(len(units)):
            d = {}
            for ch, (r, c) in unit_rows[i].items():
                d[ch] = [rr | (cc << n) for rr, cc in zip(r, c)]
            unit_aug.append(d)
        # K₀ RREF for quotienting
        k0_red, k0_piv = _rref_int(list(fr.K0))
        MAu = fr.MA
        MBu = fr.MB
        n_incons = n_cons = n_reps = 0
        n_masks = 0
        n_nodes = 0
        kdims: dict[int, int] = {}
        min_relaxed = 10 ** 9
        suspects: list[tuple[np.ndarray, np.ndarray, int]] = []
        survivors: list[dict] = []
        violations: list[dict] = []
        tight = 0
        abort_msg: list[str] = []

        # incremental forward elimination for subtree inconsistency pruning:
        # elim maps lowest-set-bit pivots to reduced augmented rows; a row
        # reducing to the bare constant bit (1 << n) is a 0 = 1 witness.
        elim: dict[int, int] = {}
        path_rows: list[int] = []  # raw augmented rows along the DFS path

        def push(rows_aug: list[int]) -> tuple[bool, list[int]]:
            added: list[int] = []
            for row in rows_aug:
                cur = row
                while cur:
                    c = (cur & -cur).bit_length() - 1
                    if c == n:
                        for a in added:
                            del elim[a]
                        return False, []
                    if c in elim:
                        cur ^= elim[c]
                    else:
                        elim[c] = cur
                        added.append(c)
                        break
            return True, added

        state_idx = [0] * len(units)

        def process_leaf() -> bool:
            nonlocal n_cons, n_reps, min_relaxed
            sol = _solve_affine_aug(list(path_rows), n)
            assert sol is not None  # elim guarantees consistency
            n_cons += 1
            part, ker = sol
            # quotient kernel by K₀: stack K₀ then kernel rows
            red = list(k0_red)
            piv = dict(k0_piv)
            extras: list[int] = []
            for kv in ker:
                cur = kv
                for c, ri in piv.items():
                    if (cur >> c) & 1:
                        cur ^= red[ri]
                if cur:
                    c = (cur & -cur).bit_length() - 1
                    for i2, rr in enumerate(red):
                        if (rr >> c) & 1:
                            red[i2] = rr ^ cur
                    piv[c] = len(red)
                    red.append(cur)
                    extras.append(cur)
            k_extra = len(extras)
            kdims[k_extra] = kdims.get(k_extra, 0) + 1
            if k_extra > KEXTRA_CAP:
                abort_msg.append(f"k_extra {k_extra} > cap")
                return False
            if n_reps + (1 << k_extra) > reps_cap:
                abort_msg.append(f"reps cap {reps_cap} exceeded")
                return False
            base_f = _int_to_vec(part, n)
            exvecs = [_int_to_vec(e, n) for e in extras]
            for t in range(1 << k_extra):
                f = base_f.copy()
                for j in range(k_extra):
                    if (t >> j) & 1:
                        f ^= exvecs[j]
                u = (MAu @ f + u0) % 2
                v = (MBu @ f + v0) % 2
                n_reps += 1
                pops = np.concatenate(
                    [fr.residue_pop(u), fr.residue_pop(v)]
                )
                relaxed = int(fr.w0_of_pop(pops).sum())
                if relaxed < min_relaxed:
                    min_relaxed = relaxed
                if relaxed < target:
                    suspects.append((f.copy(), u, v, relaxed))
                if collect_survivors and relaxed <= budget:
                    survivors.append(
                        {"f": f.copy(), "u": u, "v": v, "relaxed": relaxed}
                    )
            return True

        def dfs(i: int, spent: int) -> bool:
            nonlocal n_nodes, n_incons, n_masks
            if n_nodes > NODE_CAP or n_masks > MASK_CAP:
                abort_msg.append(
                    f"node/leaf cap ({n_nodes} nodes, {n_masks} leaves)"
                )
                return False
            if i == len(units):
                spent_here = spent
                for jj in range(len(units)):
                    cj, kj = unit_states[jj][state_idx[jj]]
                    if kj and spent_here - cj + active_cost[jj] <= budget:
                        return True  # dominated leaf ⇒ skip
                n_masks += 1
                return process_leaf()
            for si, (c, key) in enumerate(unit_states[i]):
                if spent + c + suffix_min[i + 1] > budget:
                    continue
                n_nodes += 1
                rows_here: list[int] = []
                for ch in key:
                    rows_here.extend(unit_aug[i][ch])
                ok, added = push(rows_here)
                if not ok:
                    n_incons += 1  # whole subtree pruned
                    continue
                state_idx[i] = si
                path_rows.extend(rows_here)
                cont = dfs(i + 1, spent + c)
                for a in added:
                    del elim[a]
                del path_rows[len(path_rows) - len(rows_here):]
                if not cont:
                    return False
            return True

        completed = dfs(0, 0)
        if not completed:
            return SweepResult(
                mode, target, budget, maxact, n_masks, n_incons, n_cons,
                kdims, n_reps, min_relaxed, len(suspects), None, [], tight,
                False, aborted=abort_msg[0] if abort_msg else "dfs abort",
            )
        # exact pass on suspects
        min_exact: Optional[int] = None
        undecided = 0
        if not exact_pass:
            undecided = len(suspects)
            suspects = []
        for f_s, u, v, _rel in suspects:
            ex, is_exact, flip = self.exact_cost(u, v, exclude_zero=exclude_zero)
            if min_exact is None or ex < min_exact:
                min_exact = ex
            if ex < target:
                if is_exact:
                    dressed = self.dress_witness(f_s, offsets, flip)
                    assert dressed is not None
                    f2, u2, v2 = dressed
                    assert int(u2.sum() + v2.sum()) == ex, (
                        "dressed witness weight mismatch"
                    )
                    violations.append(
                        {"f": f2, "u": u2, "v": v2, "exact": ex}
                    )
                else:
                    undecided += 1
            elif ex == target:
                tight += 1
        certified = not violations and undecided == 0
        res = SweepResult(
            mode, target, budget, maxact, n_masks, n_incons, n_cons, kdims,
            n_reps, min_relaxed, len(suspects), min_exact, violations, tight,
            certified, undecided=undecided,
        )
        if collect_survivors:
            res.survivors = survivors  # type: ignore[attr-defined]
        return res


    # -- census (LightClassification analog) --------------------------------

    def census(self, budget: int, *, mode: Optional[str] = None) -> dict:
        """All G-translation classes of nonzero boundaries of weight ≤ budget.

        Runs the floor sweep at target budget+1 with zero offsets and
        survivor collection, then expands each surviving δ-config through
        its ε-fiber (weights are linear in the flip bits) and
        canonicalizes chains under diagonal G-translation.

        Returns {"classes": [(weight, u_support, v_support)], "hist": {w: count},
        "sweep": SweepResult}.
        """
        fr = self.fr
        z0 = np.zeros(fr.n, dtype=np.uint8)
        res = self.floor_sweep(
            (z0, z0), budget + 1, mode=mode, collect_survivors=True,
            exclude_zero=True, exact_pass=False, reps_cap=12_000_000,
        )
        if res.aborted:
            return {"classes": [], "hist": {}, "sweep": res}
        survivors = getattr(res, "survivors", [])
        G = fr.group
        elems = list(G)
        # precompute translation index permutations
        perms = []
        for g in elems:
            perm = np.array(
                [G.index(G.add(e, g)) for e in elems]
            )
            perms.append(perm)

        def canonical(u: np.ndarray, v: np.ndarray) -> bytes:
            best = None
            for perm in perms:
                tu = np.zeros_like(u)
                tv = np.zeros_like(v)
                tu[perm] = u
                tv[perm] = v
                b = tu.tobytes() + tv.tobytes()
                if best is None or b < best:
                    best = b
            return best

        seen: dict[bytes, tuple[int, np.ndarray, np.ndarray]] = {}
        F = self._flips_matrix()
        if F is None:
            return {"classes": [], "hist": {},
                    "sweep": res, "aborted": "eps dim > cap"}
        seen_delta: set[bytes] = set()
        for sv in survivors:
            u, v = sv["u"], sv["v"]
            dkey = (
                fr.residues(u.astype(np.uint8)).tobytes()
                + fr.residues(v.astype(np.uint8)).tobytes()
            )
            if dkey in seen_delta:
                continue
            seen_delta.add(dkey)
            costs, _pops, _eps0 = self._config_cost_vector(u, v)
            f_s = sv["f"]
            for x in np.flatnonzero(costs <= budget):
                dressed = self.dress_witness(
                    f_s, (z0, z0), F[int(x)].astype(np.uint8)
                )
                assert dressed is not None
                _f2, u2, v2 = dressed
                w = int(u2.sum() + v2.sum())
                assert w == int(costs[int(x)]), (w, int(costs[int(x)]))
                if w == 0:
                    continue
                key = canonical(u2, v2)
                if key not in seen or seen[key][0] > w:
                    seen[key] = (w, u2.copy(), v2.copy())
        hist: dict[int, int] = {}
        for w, _u, _v in seen.values():
            hist[w] = hist.get(w, 0) + 1
        classes = [
            (w, sorted(int(i) for i in np.flatnonzero(u)),
             sorted(int(i) for i in np.flatnonzero(v)))
            for w, u, v in seen.values()
        ]
        classes.sort(key=lambda t: (t[0], t[1], t[2]))
        return {"classes": classes, "hist": dict(sorted(hist.items())),
                "sweep": res}


# --------------------------------------------------------------------------
# safe classes: kernel, orbits, seamC offsets
# --------------------------------------------------------------------------


def kernel_basis(A: Poly, B: Poly) -> np.ndarray:
    """Basis of K = Ann(A) ∩ Ann(B) = ker ∂₂ as rows (uint8, length n)."""
    from .linalg import nullspace_f2

    M = np.vstack([circulant(A), circulant(B)])
    return nullspace_f2(M)


def _translate_vec(G: AbelianGroup, g: tuple[int, ...], v: np.ndarray) -> np.ndarray:
    out = np.zeros_like(v)
    elems = list(G)
    for i, e in enumerate(elems):
        out[G.index(G.add(e, g))] = v[i]
    return out


def kernel_orbit_reps(A: Poly, B: Poly) -> list[np.ndarray]:
    """G-translation orbit representatives of ker ∂₂ ∖ 0."""
    G = A.group
    K = kernel_basis(A, B)
    dim = K.shape[0]
    if dim == 0:
        return []
    all_elts = []
    for mask in range(1, 1 << dim):
        v = np.zeros(K.shape[1], dtype=np.uint8)
        for j in range(dim):
            if (mask >> j) & 1:
                v ^= K[j]
        all_elts.append(v)
    seen: set[bytes] = set()
    reps = []
    for v in all_elts:
        if v.tobytes() in seen:
            continue
        reps.append(v)
        for g in G:
            seen.add(_translate_vec(G, g, v).tobytes())
    return reps


def seam_offsets(
    A: Poly, B: Poly, axis: int
) -> list[tuple[np.ndarray, np.ndarray, np.ndarray]]:
    """[(ζ, seamC_u, seamC_v)] for kernel orbit reps, doubling along `axis`.

    seamC per Prop A14.1(3): lift ζ by the canonical section (axis
    coordinate < ℓ), apply the cover ∂₂, read the carry window
    (coordinate in [ℓ, 2ℓ)).
    """
    G = A.group
    orders = list(G.orders)
    ell = orders[axis]
    cov_orders = list(orders)
    cov_orders[axis] = 2 * ell
    Gc = AbelianGroup(tuple(cov_orders))
    Ac = Poly.from_support(A.support, Gc)
    Bc = Poly.from_support(B.support, Gc)
    MAc = circulant(Ac)
    MBc = circulant(Bc)
    out = []
    for zeta in kernel_orbit_reps(A, B):
        lift = np.zeros(Gc.cardinality, dtype=np.uint8)
        for i, e in enumerate(G):
            lift[Gc.index(e)] = zeta[i]  # section: coords unchanged (< ℓ)
        tu = (MAc @ lift) % 2
        tv = (MBc @ lift) % 2
        su = np.zeros(G.cardinality, dtype=np.uint8)
        sv = np.zeros(G.cardinality, dtype=np.uint8)
        for i, e in enumerate(G):
            ec = list(e)
            ec[axis] += ell
            j = Gc.index(tuple(ec))
            su[i] = tu[j]
            sv[i] = tv[j]
        out.append((zeta, su, sv))
    return out


@dataclass
class SafeFloorReport:
    code_label: str
    axis: int
    fiber: tuple[int, ...]
    target: int
    per_class: list[SweepResult]
    certified: bool
    refuted: bool
    notes: list[str] = field(default_factory=list)

    def summary(self) -> str:
        head = (
            f"{self.code_label} axis={self.axis} fiber z={self.fiber} "
            f"target {self.target}: "
        )
        if self.certified:
            head += "SAFE-FLOOR CERTIFIED"
        elif self.refuted:
            head += "SAFE-FLOOR REFUTED"
        else:
            head += "UNDECIDED"
        return head + "".join(f"\n  {r.summary()}" for r in self.per_class)


def safe_floor_certify(
    A: Poly,
    B: Poly,
    axis: int,
    target: int,
    z: Optional[tuple[int, ...]] = None,
    *,
    mode: Optional[str] = None,
    progress: bool = False,
) -> SafeFloorReport:
    """End-to-end safe-floor decision for one (code, axis) at `target`.

    Sweeps every G-orbit class offset seamC(ζ). Per class, feasible fiber
    frames are tried in score order until one decides (floors are
    frame-independent, so classes may mix frames). Certified iff every
    class certifies on some frame; refuted iff some class has a verified
    non-boundary violation witness.
    """
    if z is not None:
        cand = [z]
    else:
        feas = [s for s in best_frames(A, B, target) if s.feasible]
        # linked frames first (cheap paired sweeps), then at most two
        # unpaired fallbacks — further unlinked frames share the failure
        # shape and only burn node budget.
        linked = [s.z for s in feas if s.link_shift is not None]
        unlinked = [s.z for s in feas if s.link_shift is None]
        cand = linked + unlinked[:2]
        if not cand:
            raise ValueError("no feasible fiber found; inspect best_frames()")
    frames = [FiberFrame(A, B, zz) for zz in cand]
    sweeps = [FiberSweep(fr) for fr in frames]
    label = f"{A.group.label()}|A={A.canonical_string()}|B={B.canonical_string()}"
    reports = []
    notes = []
    refuted = False
    # boundary space RREF for non-boundary verification of witnesses:
    # im ∂₂ is spanned by the 2n-length columns (A δ_g | B δ_g), g ∈ G
    bnd_rows = [
        _vec_to_int(np.concatenate([frames[0].MA[:, j], frames[0].MB[:, j]]))
        for j in range(frames[0].n)
    ]
    bnd_red, bnd_piv = _rref_int(bnd_rows)

    def _in_boundaries(w: np.ndarray) -> bool:
        cur = _vec_to_int(w)
        for c, ri in bnd_piv.items():
            if (cur >> c) & 1:
                cur ^= bnd_red[ri]
        return cur == 0

    used_frames: list[tuple[int, ...]] = []
    for zeta, su, sv in seam_offsets(A, B, axis):
        best_res: Optional[SweepResult] = None
        best_z: Optional[tuple[int, ...]] = None
        for fr, sw in zip(frames, sweeps):
            res = sw.floor_sweep((su, sv), target, mode=mode, progress=progress)
            if best_res is None:
                best_res, best_z = res, fr.z
            decided = res.certified or (res.violations and not res.aborted)
            if decided:
                best_res, best_z = res, fr.z
                break
        res = best_res
        assert res is not None and best_z is not None
        used_frames.append(best_z)
        # verify any violations end-to-end: weight + membership + non-boundary
        real_viols = []
        for viol in res.violations:
            u, v = viol["u"], viol["v"]
            w = int(u.sum() + v.sum())
            viol["weight"] = w
            assert w == viol["exact"]
            chain = np.concatenate([u, v])
            assert not _in_boundaries(chain), (
                "witness is a boundary — safe-class setup violated"
            )
            viol["verified"] = True
            real_viols.append(viol)
        if real_viols:
            refuted = True
        if res.undecided:
            notes.append(f"{res.undecided} suspects undecided (ε-enum cap)")
        reports.append(res)
    if len(set(used_frames)) > 1:
        notes.append(f"classes decided on mixed frames: {used_frames}")
    certified = bool(reports) and all(r.certified for r in reports)
    return SafeFloorReport(
        label, axis, used_frames[0] if used_frames else cand[0], target,
        reports, certified and not refuted, refuted, notes,
    )
