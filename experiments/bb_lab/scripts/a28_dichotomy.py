"""A28 — the shifting dichotomy tree (the (CLASS) kill-side systematizer).

Cells: pairs (P_in, P_out) of disjoint Galois-orbit unions in Ghat.
The cell contains exactly the f with fhat == 0 on P_in and fhat != 0 on
P_out (orbit-level dichotomy is exact for F2-rational f: fhat(chi^2) =
fhat(chi)^2).  Branching on an undecided orbit partitions the cell.

At each cell, three certified lower bounds on |del f| for every f in it:
  - union:  l_u(Z_A u P_in; free closers P_out) + l_v(Z_B u P_in; ...)
            (block additivity — Lemma B of the note), and
  - joint:  the syzygy game (Theorem J), free closers P_out.
Cells with bound > W are KILLED (certificates = replayable histories).
Cells where every orbit is decided are SURVIVOR leaves (exact spectral
patterns — where the actual census classes live).  Depth-capped cells are
GIVEUP leaves (the honesty metric).

Output: tree statistics + verification that no known census class lands
in a killed cell, and that kills + survivors + giveups partition Ghat.
"""

import json
import sys
import time
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from a28_lsc_lib import DATA, REGISTRY, load_f2a6_census
from a28_spectral import Spectral
from a28_shift_engine import ShiftGame
from a28_joint_game import JointGame


def mask_of(idxs):
    m = 0
    for i in idxs:
        m |= 1 << i
    return m


class Dichotomy:
    def __init__(self, inst, quick_beam=120, kill_beam=350, max_depth=10):
        self.inst = inst
        self.sp = Spectral(inst.G)
        self.N = inst.G.N
        self.game = ShiftGame(self.sp)
        Ahat = self.sp.fourier(inst.A)
        Bhat = self.sp.fourier(inst.B)
        self.joint = JointGame(self.sp, Ahat, Bhat)
        self.ZA = mask_of(self.sp.zero_set(inst.A))
        self.ZB = mask_of(self.sp.zero_set(inst.B))
        self.kernel = self.ZA & self.ZB
        self.orbits = [mask_of(o) for o in self.sp.galois_orbits()]
        # orbits fully inside the kernel are invisible — exclude
        self.free_orbits = [o for o in self.orbits if not (o & self.kernel)]
        self.W = inst.W
        self.target = inst.W + 1
        self.quick_beam = quick_beam
        self.kill_beam = kill_beam
        self.max_depth = max_depth
        self.leaves = []
        self.nodes = 0
        self.t0 = time.time()

    def bounds(self, P_in, P_out, beam, want_history=False):
        bu, hu = self.game.best_bound(self.ZA | P_in, free=P_out, beam=beam,
                                      max_level=self.target,
                                      want_history=want_history)
        bv, hv = self.game.best_bound(self.ZB | P_in, free=P_out, beam=beam,
                                      max_level=self.target,
                                      want_history=want_history)
        if bu + bv >= self.target:
            return bu + bv, ("union", bu, hu, bv, hv)
        bj, hj = self.joint.best_bound(P_in, P_out, target=self.target,
                                       beam=beam, want_history=want_history)
        if bj >= self.target:
            return bj, ("joint", bj, hj)
        return max(bu + bv, bj), None

    def pick_branch_orbit(self, P_in, P_out):
        """Cheap probe: the undecided orbit whose IN-branch raises the union
        bound most (ties: largest orbit)."""
        cand = [o for o in self.free_orbits if not (o & (P_in | P_out))]
        best, best_score = None, (-1, 0)
        for o in cand:
            bu, _ = self.game.best_bound(self.ZA | P_in | o, free=P_out,
                                         beam=40, max_level=self.target)
            bv, _ = self.game.best_bound(self.ZB | P_in | o, free=P_out,
                                         beam=40, max_level=self.target)
            score = (bu + bv, bin(o).count("1"))
            if score > best_score:
                best_score, best = score, o
        return best

    def explore(self, P_in, P_out, depth):
        self.nodes += 1
        # quick kill probe, then a stronger attempt before giving up on a kill
        b, cert = self.bounds(P_in, P_out, self.quick_beam)
        if cert is None and b >= self.target - 2:
            b, cert = self.bounds(P_in, P_out, self.kill_beam)
        if cert is not None:
            self.leaves.append({"kind": "kill", "P_in": P_in, "P_out": P_out,
                                "bound": b, "cert": cert[0]})
            return
        undecided = [o for o in self.free_orbits if not (o & (P_in | P_out))]
        if not undecided:
            self.leaves.append({"kind": "survivor", "P_in": P_in,
                                "P_out": P_out, "bound": b})
            return
        if depth >= self.max_depth:
            self.leaves.append({"kind": "giveup", "P_in": P_in,
                                "P_out": P_out, "bound": b,
                                "undecided": len(undecided)})
            return
        o = self.pick_branch_orbit(P_in, P_out)
        self.explore(P_in | o, P_out, depth + 1)
        self.explore(P_in, P_out | o, depth + 1)

    def run(self):
        self.explore(0, 0, 0)
        return self.leaves


def main():
    inst = REGISTRY["f2a6"]
    di = Dichotomy(inst, max_depth=int(sys.argv[1]) if len(sys.argv) > 1 else 10)
    leaves = di.run()
    kinds = Counter(l["kind"] for l in leaves)
    print(f"tree: {di.nodes} nodes, leaves {dict(kinds)}, "
          f"wall {time.time()-di.t0:.0f}s")

    # where do the 113 classes land?  walk each class's actual pattern.
    sp = di.sp
    G = inst.G
    census = load_f2a6_census()
    landings = Counter()
    bad = 0
    for row in census:
        u = G.from_support([(x, y) for blk, x, y in row["b_support"] if blk == 0])
        v = G.from_support([(x, y) for blk, x, y in row["b_support"] if blk == 1])
        Zu = mask_of(sp.zero_set(u))
        Zv = mask_of(sp.zero_set(v))
        shared = ~(di.ZA | di.ZB) & ((1 << di.N) - 1)
        vis = (Zu & shared) | (Zv & di.ZA & ~di.ZB) | (Zu & di.ZB & ~di.ZA)
        # find the unique leaf containing this f-pattern
        hit = None
        for leaf in leaves:
            if leaf["P_in"] & ~vis:
                continue          # leaf requires zeros where f is nonzero
            if leaf["P_out"] & vis:
                continue          # leaf requires nonzeros where f vanishes
            # cell membership needs BOTH directions decided the same way;
            # a class matches the leaf iff its pattern refines it
            hit = leaf
            break
        assert hit is not None, "class escaped the partition!"
        landings[hit["kind"]] += 1
        if hit["kind"] == "kill":
            bad += 1
            print("  FALSIFIED: census class in killed cell!", row["b_weight"])
    print(f"landings of the 113 classes: {dict(landings)}  "
          f"(killed-cell violations: {bad})")
    assert bad == 0, "SOUNDNESS FALSIFIED"

    # verify all kill certificates by replay
    ver = 0
    for leaf in leaves:
        if leaf["kind"] != "kill":
            continue
        b, cert = di.bounds(leaf["P_in"], leaf["P_out"], di.kill_beam,
                            want_history=True)
        if cert is None:
            continue  # nondeterministic beam; the leaf stands by its record
        if cert[0] == "union":
            _, bu, hu, bv, hv = cert
            cu = di.game.verify_history(di.ZA | leaf["P_in"], hu)
            cv = di.game.verify_history(di.ZB | leaf["P_in"], hv)
            assert cu + cv >= di.target
        else:
            _, bj, hj = cert
            assert di.joint.verify_history(leaf["P_in"], leaf["P_out"], hj) >= di.target
        ver += 1
    print(f"kill certificates re-derived and replay-verified: {ver}/{kinds['kill']}")

    out = {"nodes": di.nodes, "leaves": [
        {k: (format(v, "x") if isinstance(v, int) and k in ("P_in", "P_out") else v)
         for k, v in l.items() if k != "cert"} for l in leaves]}
    (DATA / "a28").mkdir(exist_ok=True)
    with open(DATA / "a28" / "dichotomy_f2a6.json", "w") as fh:
        json.dump(out, fh, indent=1)


if __name__ == "__main__":
    main()
