#!/usr/bin/env python3
"""Force the pivot's colour onto a *proper subset* of an odd magic circle.

A ball grown from the unit steps is invariant under the polygon's own
rotation, because the step set ±zeta_n^k is closed under multiplying by
zeta_n.  In an invariant graph no proper subset of the target orbit is ever
forced -- checked exhaustively at order 3, where every pair separates and only
all three force.  So the core comes out equal to the orbit, 3 at order 3 and
11 at order 11, while the capacities are 2 and 5.  The orbit is always one
size too many.

Adding an *off-centre* ball fixes that.  No rotation about the pivot preserves
it, it only adds constraints so forcing can only get easier, and at order 3 it
took the core from 3 straight down to 1.

So: symmetric ball for the targets, off-centre balls for the asymmetry, then
shrink the core greedily and compare it against the capacity.
"""
import os
import sys
import time
from fractions import Fraction

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from hn.cyclograph import CycloPoint, build_cyclo_graph, cyclo_generated
from hn.cyclotomic import CycloField, moser_rotation, unit_polygon
from hn.spindle import SeparationTest
from hn.transversal import trapping_bound

K = int(os.environ.get("HN_K", "4"))
ORDER = int(os.environ.get("HN_ORDER", "11"))
RADIUS = float(os.environ.get("HN_RADIUS", "2.2"))
OFF = float(os.environ.get("HN_OFF", "1.6"))
DEPTH = int(os.environ.get("HN_DEPTH", "8"))
CAP = int(os.environ.get("HN_CAP", "9000"))
SHIFTS = int(os.environ.get("HN_SHIFTS", "2"))


def ball(F, gens, centre, radius, cap, extra=()):
    seeds = [CycloPoint(F, centre)] + list(extra)
    cz = F.to_complex(centre)
    r2 = radius * radius
    return [p for p in cyclo_generated(F, seeds, gens, DEPTH, cap=cap)
            if abs(p.z - cz) ** 2 <= r2]


def shrink(st, core):
    """Drop targets while the disjunction stays forced."""
    core = list(core)
    changed = True
    while changed and len(core) > 1:
        changed = False
        for t in list(core):
            trial = [x for x in core if x != t]
            sep, _ = st.run(subset=trial)
            if not sep:
                core = trial
                changed = True
                break
    return core


def main() -> None:
    F = CycloField(33)
    poly = unit_polygon(F, ORDER)
    rho = moser_rotation(F)
    # A low-rank step set, because rank decides whether there is an interior.
    # The 66 roots of unity generate Z[zeta_33], rank 20, and a ball over it is
    # all boundary: depth 3 gave an empty 8-core. omega and rho generate
    # Q(omega, sqrt(-11)), degree 4, with denominators bounded by 6 -- the
    # classical setting, where depth pays off and interior points keep their
    # full complement of neighbours.
    w = F.zeta(F.n // 3)
    base = [F.one(), w, F.mul(w, w)]
    steps = []
    for b in base:
        for mult in (F.one(), rho, F.conj(rho)):
            z = F.mul(b, mult)
            steps += [z, F.neg(z)]
    gens = [("add", z) for z in steps]
    cap = trapping_bound(ORDER)
    print(f"order {ORDER}, k={K}: capacity {cap}, polygon radius "
          f"{abs(F.to_complex(poly[0])):.4f}", flush=True)

    pts = {p.c: p for p in ball(F, gens, F.zero(), RADIUS, CAP,
                                [CycloPoint(F, v) for v in poly])}
    print(f"  symmetric ball: {len(pts)}", flush=True)
    steps = F.unit_steps()
    for s in range(SHIFTS):
        centre = tuple(Fraction(x) for x in steps[1 + s * 3])
        for p in ball(F, gens, centre, OFF, CAP):
            pts.setdefault(p.c, p)
        print(f"  + off-centre ball {s + 1}: {len(pts)}", flush=True)

    t0 = time.time()
    g = build_cyclo_graph(list(pts.values()))
    idx = {p.c: i for i, p in enumerate(g.vertices)}
    pivot = idx[F.zero()]
    targets = [idx[v] for v in poly if v in idx]
    print(f"  {g}  {len(targets)} targets  [{time.time() - t0:.0f}s]", flush=True)

    st = SeparationTest(g, K, pivot, targets)
    try:
        t0 = time.time()
        sep, core = st.run(subset=targets)
        if sep:
            print(f"  separable -- no forcing at all  [{time.time() - t0:.0f}s]",
                  flush=True)
            return
        print(f"  FORCED, raw core {len(core)}  [{time.time() - t0:.0f}s]",
              flush=True)
        t0 = time.time()
        small = shrink(st, core)
        verdict = "WITHIN CAPACITY" if len(small) <= cap else "still too many"
        print(f"  minimal core {len(small)} against capacity {cap}: {verdict}"
              f"  [{time.time() - t0:.0f}s]", flush=True)
        if len(small) <= cap:
            print(f"  *** targets {small} ***", flush=True)
    finally:
        st.close()


if __name__ == "__main__":
    main()
