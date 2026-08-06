"""A28 — the abelian shift-bound game engine (van Lint–Wilson transported).

Single-block game (Theorem S of the A28 note).  For a nonzero word c in
F2[G] (G odd here) with EXACT zero set Z(c), and any O subseteq Z(c):
a set I built by the moves

    close:  requires I subseteq O and b notin Z(c);  I -> I u {b}
    shift:  I -> mu*I  (any character mu; additive on the (s,t) grid)

satisfies rank{chi|_supp(c) : chi in I} = |I|, hence |c| >= |I|.
(Proof: close adds a vector outside span(V_I) subseteq c-perp; shift is a
Hadamard multiplication by a unit vector.)

In census mode the oracle O = Z_hyp u Z(A) for a hypothesized exact
Z(f) = Z_hyp, and closers b split into FREE (hypothesized nonzero, i.e.
in Z_out) and BRANCHING (unknown — using one forks the dichotomy tree).

The search is a beam over canonical states; any history found is a valid
certificate (search incompleteness only weakens bounds, never soundness).
"""

from __future__ import annotations

import sys
from itertools import combinations
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from a28_lsc_lib import GroupAlg, rref
from a28_spectral import Spectral


class ShiftGame:
    def __init__(self, sp: Spectral):
        self.sp = sp
        self.N = sp.G.N
        # translate_table[mu] = tuple perm of indices
        self.tt = [tuple(sp.translate_chi(ci, mu) for ci in range(self.N))
                   for mu in range(self.N)]

    def translate_mask(self, mask: int, mu: int) -> int:
        t = self.tt[mu]
        out = 0
        while mask:
            low = mask & -mask
            out |= 1 << t[low.bit_length() - 1]
            mask ^= low
        return out

    def fingerprint(self, mask: int) -> tuple:
        """Translation-invariant state key: sorted difference multiset."""
        idxs = []
        m = mask
        while m:
            low = m & -m
            idxs.append(low.bit_length() - 1)
            m ^= low
        G = self.sp.G
        diffs = []
        for i in idxs:
            s1, t1 = divmod(i, G.m)
            for j in idxs:
                s2, t2 = divmod(j, G.m)
                diffs.append(((s1 - s2) % G.l) * G.m + (t1 - t2) % G.m)
        return tuple(sorted(diffs))

    def best_bound(self, oracle: int, free: int | None = None,
                   beam: int = 600, max_level: int = 24,
                   want_history: bool = False):
        """Max |I| found.  oracle = mask of certified-zero characters;
        closers allowed only from `free` (mask) if given, else any
        non-oracle character.  Returns (bound, history|None).
        History entries: (mu, b) — shift by mu then close b."""
        full = (1 << self.N) - 1
        closer_pool = (free if free is not None else full) & ~oracle
        if closer_pool == 0:
            return 0, []
        # states: dict fingerprint -> (mask, history)
        states = {(): (0, [])}
        best, best_hist = 0, []
        for level in range(max_level):
            nxt = {}
            for fp, (mask, hist) in states.items():
                for mu in range(self.N):
                    tm = self.translate_mask(mask, mu) if mask else 0
                    if tm & ~oracle:
                        continue
                    pool = closer_pool & ~tm
                    while pool:
                        low = pool & -pool
                        b = low.bit_length() - 1
                        pool ^= low
                        nm = tm | low
                        nfp = self.fingerprint(nm)
                        if nfp not in nxt:
                            nxt[nfp] = (nm, hist + [(mu, b)] if want_history else hist)
                            if len(nxt) >= beam * 4:
                                break
                    if len(nxt) >= beam * 4:
                        break
            if not nxt:
                break
            # keep a diverse beam (dict order is insertion; sample evenly)
            items = list(nxt.items())
            if len(items) > beam:
                step = len(items) / beam
                items = [items[int(i * step)] for i in range(beam)]
            states = dict(items)
            best = level + 1
            if want_history:
                best_hist = next(iter(states.values()))[1]
        return best, best_hist

    def verify_history(self, oracle: int, hist: list) -> int:
        """Replay legality: each step shifts current I into the oracle and
        closes a non-oracle character.  Returns |I| (the certified bound)."""
        mask = 0
        for mu, b in hist:
            tm = self.translate_mask(mask, mu) if mask else 0
            assert tm & ~oracle == 0, "close with I not inside oracle"
            assert not (oracle >> b) & 1, "closer inside oracle"
            assert not (tm >> b) & 1
            mask = tm | (1 << b)
        return bin(mask).count("1")


# ----------------------------------------------------------------------
# exact minima for validation: ideal(O) = {u : Z(u) >= O}
# ----------------------------------------------------------------------

def ideal_basis(sp: Spectral, O_mask: int) -> list[int]:
    """F2 basis of {u : uhat vanishes on O} via kernel of the eval map."""
    N = sp.G.N
    # build rows: for each u-basis monomial i, the GF16 evals at O as 4 bits each
    ocols = [ci for ci in range(N) if (O_mask >> ci) & 1]
    mrows = []
    for i in range(N):
        row = 0
        for j, ci in enumerate(ocols):
            row |= sp.chi[ci][i] << (4 * j)
        mrows.append(row)
    # kernel of the map F2^N -> (F2^4)^|O| given by mrows (rows = images of e_i)
    # augment: [image | e_i], reduce on image columns; kernel rows have image 0
    W = 4 * len(ocols)
    aug = [mrows[i] | (1 << (W + i)) for i in range(N)]
    # eliminate on image columns; surviving rows with zero image = kernel
    kernel = []
    red = aug[:]
    for c in range(W):
        bit = 1 << c
        pivrow = None
        for r in red:
            if r & bit:
                pivrow = r
                break
        if pivrow is None:
            continue
        red = [(r ^ pivrow) if (r is not pivrow and r & bit) else r for r in red]
        red.remove(pivrow)
    for r in red:
        assert r & ((1 << W) - 1) == 0
        kernel.append(r >> W)
    return kernel


def min_weight_bz(rows: list[int], ncols: int, cap: int = 40) -> int:
    """Exact min weight of the span via BZ escalation (small dims only)."""
    basis, piv = rref(rows, ncols)
    kappa = len(basis)
    if kappa == 0:
        return 0
    best = min(bin(b).count("1") for b in basis)
    # one systematic set: after enumerating all combos of size <= r, every
    # unseen codeword has > r ones in the pivot set, so min >= min(best, r+1);
    # escalate until r + 1 >= best, then best is exact.
    r = 0
    while r + 1 < best and r < min(kappa, cap):
        r += 1
        for comb in combinations(range(kappa), r):
            v = 0
            for i in comb:
                v ^= basis[i]
            w = bin(v).count("1")
            if 0 < w < best:
                best = w
    return best
