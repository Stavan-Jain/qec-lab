"""A28 — the two-block SYZYGY shift game (Theorem J of the note).

For b = (u, v) = (fA, fB), work in Fbar2^(supp u ⊔ supp v).  For a
character chi and a projective line lambda = [l1 : l2] define

    V_{chi,lambda} = (l1 * chi|_supp(u) ; l2 * chi|_supp(v)).

Pairing: <V_{chi,lambda}, b> = l1*uhat(chi) + l2*vhat(chi)
                             = fhat(chi) * (l1*Ahat(chi) + l2*Bhat(chi)).

So V_{chi,lambda} is orthogonal to b iff
    chi in Z(f)  (hypothesis-dependent), or
    chi in Z_A cap Z_B (the kernel — always), or
    lambda = ann(chi) := [Bhat(chi) : Ahat(chi)]  (the SYZYGY line — always,
    with no hypothesis on f: this is the two-block enrichment; its slope
    function rho = Bhat/Ahat is the abstract home of gross's "C-ratios").

Moves: shift (chi,lambda) -> (mu*chi, lambda) [Hadamard by a unit vector];
close b* = (chi,lambda) with pairing provably nonzero: chi in P_out
(hypothesized fhat != 0) and lambda != ann(chi), requiring the whole
current I orthogonal.  Per-chi at most 2 lines (lambda -> V is linear in
lambda, so a third line is dependent).  Invariant rank V_I = |I|, hence
|u| + |v| >= |I|.

Lines coded 0..16: code c in 0..15 = [c : 1], code 16 = [1 : 0].
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from a28_spectral import Spectral, gf_mul, _EXP, _LOG

LINF = 16  # the line [1 : 0]


def line_code(a: int, b: int) -> int:
    """Projective line through (a, b) != 0."""
    if b == 0:
        return LINF
    if a == 0:
        return 0
    return _EXP[(_LOG[a] - _LOG[b]) % 15]


class JointGame:
    def __init__(self, sp: Spectral, Ahat: list[int], Bhat: list[int]):
        self.sp = sp
        self.N = sp.G.N
        self.kernel = 0
        self.ann = [None] * self.N
        for ci in range(self.N):
            if Ahat[ci] == 0 and Bhat[ci] == 0:
                self.kernel |= 1 << ci
            else:
                self.ann[ci] = line_code(Bhat[ci], Ahat[ci])
        self.tt = [tuple(sp.translate_chi(ci, mu) for ci in range(self.N))
                   for mu in range(self.N)]

    def _orthogonal(self, pairs, P_in: int) -> bool:
        zero = self.kernel | P_in
        for chi, lam in pairs:
            if (zero >> chi) & 1:
                continue
            if lam != self.ann[chi]:
                return False
        return True

    def _canon(self, pairs) -> tuple:
        best = None
        for mu in range(self.N):
            t = self.tt[mu]
            s = tuple(sorted((t[chi], lam) for chi, lam in pairs))
            if best is None or s < best:
                best = s
        return best

    def best_bound(self, P_in: int, P_out: int, target: int,
                   beam: int = 500, lam_menu=None, want_history=False):
        """Max |I| found (early exit at `target`).  Free closers only:
        chi in P_out, lambda != ann(chi).  Returns (bound, history).
        History: (mu, chi, lam) = shift all by mu, then close (chi,lam)."""
        if P_out == 0:
            return 0, []
        if lam_menu is None:
            lam_menu = sorted({self.ann[c] for c in range(self.N)
                               if self.ann[c] is not None} | {0, LINF})
        out_chis = [c for c in range(self.N) if (P_out >> c) & 1]
        states = {(): ((), [])}
        best, best_hist = 0, []
        for level in range(target):
            nxt = {}
            for _, (pairs, hist) in states.items():
                # translations keeping I orthogonal
                for mu in range(self.N):
                    t = self.tt[mu]
                    tp = tuple((t[chi], lam) for chi, lam in pairs)
                    if not self._orthogonal(tp, P_in):
                        continue
                    occ = {}
                    for chi, lam in tp:
                        occ.setdefault(chi, set()).add(lam)
                    for chi in out_chis:
                        have = occ.get(chi, set())
                        if len(have) >= 2:
                            continue
                        for lam in lam_menu:
                            if lam == self.ann[chi] or lam in have:
                                continue
                            np_ = tuple(sorted(tp + ((chi, lam),)))
                            key = self._canon(np_)
                            if key not in nxt:
                                nxt[key] = (np_, hist + [(mu, chi, lam)]
                                            if want_history else hist)
                            if len(nxt) >= beam * 3:
                                break
                        if len(nxt) >= beam * 3:
                            break
                    if len(nxt) >= beam * 3:
                        break
                    if not pairs:
                        break  # empty state: one mu suffices
            if not nxt:
                break
            items = list(nxt.items())
            if len(items) > beam:
                step = len(items) / beam
                items = [items[int(i * step)] for i in range(beam)]
            states = dict(items)
            best = level + 1
            if want_history:
                best_hist = next(iter(states.values()))[1]
            if best >= target:
                break
        return best, best_hist

    def verify_history(self, P_in: int, P_out: int, hist) -> int:
        """Replay legality; returns the certified |I|."""
        pairs = ()
        for mu, chi, lam in hist:
            t = self.tt[mu]
            tp = tuple((t[c], l) for c, l in pairs)
            assert self._orthogonal(tp, P_in), "close with non-orthogonal I"
            assert (P_out >> chi) & 1, "closer chi not hypothesized nonzero"
            assert not (self.kernel >> chi) & 1 and not (P_in >> chi) & 1
            assert lam != self.ann[chi], "closer on the syzygy line"
            have = {l for c, l in tp if c == chi}
            assert lam not in have and len(have) < 2
            pairs = tuple(sorted(tp + ((chi, lam),)))
        return len(pairs)
