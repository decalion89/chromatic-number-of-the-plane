#!/usr/bin/env python3
"""The order-15 magic circle: reachable, and worth 4 instead of 2.

Two obstructions decide which magic circles a unit-step walk can ever touch.

*Niven* rules out every rational radius but d^2 = 1/3, whose cycle is the
triangle and whose capacity is 2 -- the classical wall.

*Ramification* rules out the prime-power orders.  A unit step has u conj(u) =
1, so at a prime complex conjugation fixes its valuation is 0, and every point
of a connected component inherits v >= 0.  The order-n magic radius squared is
1/((1 - zeta_n)(1 - zeta_n^-1)), which for n a prime power is the ramified
prime above p and has valuation -2.  Measured: 0 points on the order-5 circle
and 0 on the order-11 one, with 36481 points and 11-gon vertices sitting in a
different component from the pivot.

But 15 is not a prime power, so 1 - zeta_15 is a unit, the radius is an
algebraic integer, and nothing forbids it.  Measured: 30 points on it at depth
3.  Capacity 4.

So: find the circle's points, group them into orbits of the order-15 rotation,
and ask whether the pivot's colour is forced onto at most four of one orbit.
"""
import os
import sys
import time
from fractions import Fraction

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from hn.cyclograph import CycloPoint, UnitDistanceGraph
from hn.cyclotomic import CycloField
from hn.quadext import QuadExtField
from hn.spindle import SeparationTest
from hn.transversal import trapping_bound

K = int(os.environ.get("HN_K", "4"))
ORDER = int(os.environ.get("HN_ORDER", "15"))
DEPTH = int(os.environ.get("HN_DEPTH", "4"))
RADIUS = float(os.environ.get("HN_RADIUS", "2.8"))
CAP = int(os.environ.get("HN_CAP", "60000"))
OFF = int(os.environ.get("HN_OFF", "2"))
OFFDIST = int(os.environ.get("HN_OFFDIST", "3"))
OFFSCALE = float(os.environ.get("HN_OFFSCALE", "0.8"))
MULT = int(os.environ.get("HN_MULT", "0"))


def unit_step_set(F, mult: int = 0):
    """Roots of unity, the two split-prime quotients, and their products.

    rho = (5 + sqrt(-11))/6 is the split-prime quotient at 3 -- the Moser
    rotation -- and u = (-1 + 3 sqrt(-11))/10 the one at 5.  Multiplying the
    roots of unity by their powers gives more steps, and breaks the closure
    under zeta_15 that would otherwise make every ball invariant.
    """
    rho = F.moser_rotation()
    u5 = F.scale(F.add(F.rational(-1), F.scale(F.radical(), 3)), Fraction(1, 10))
    roots = [F.root_of_unity(15, k) for k in range(15)]
    steps = list(roots) + [rho, F.conj(rho), u5, F.conj(u5)]
    for _ in range(mult):
        grown = list(steps)
        for a in (rho, F.conj(rho), u5, F.conj(u5)):
            grown += [F.mul(r, a) for r in roots]
        steps = grown
    steps += [F.neg(s) for s in list(steps)]
    return [s for s in {t: None for t in steps} if F.norm2(s) == F.one()]


def grow(F, steps, seeds, depth, radius, cap, centre=0j):
    seen = {c: None for c in seeds}
    frontier = list(seen)
    for _ in range(depth):
        nxt = []
        for c in frontier:
            for s in steps:
                q = F.add(c, s)
                if q in seen or abs(F.to_complex(q) - centre) > radius:
                    continue
                seen[q] = None
                nxt.append(q)
                if len(seen) >= cap:
                    return list(seen)
        frontier = nxt
        if not frontier:
            break
    return list(seen)


def wire(F, coeffs, steps):
    index = {c: i for i, c in enumerate(coeffs)}
    adj = [set() for _ in coeffs]
    for s in steps:
        for i, c in enumerate(coeffs):
            j = index.get(F.add(c, s))
            if j is not None and j != i:
                adj[i].add(j)
                adj[j].add(i)
    return UnitDistanceGraph([CycloPoint(F, c) for c in coeffs], adj), index


def shrink(st, core):
    core = list(core)
    changed = True
    while changed and len(core) > 1:
        changed = False
        for t in list(core):
            trial = [x for x in core if x != t]
            if not st.run(subset=trial)[0]:
                core, changed = trial, True
                break
    return core


def main() -> None:
    F = QuadExtField(CycloField(15), -11)
    steps = unit_step_set(F, MULT)
    cap_r = trapping_bound(ORDER)
    print(f"order {ORDER}, k={K}: capacity {cap_r}, {len(steps)} unit steps",
          flush=True)

    t0 = time.time()
    coeffs = grow(F, steps, [F.zero()], DEPTH, RADIUS, CAP)
    print(f"  ball: {len(coeffs)} points  [{time.time() - t0:.0f}s]", flush=True)

    # break the rotation symmetry: an invariant ball never forces a proper
    # subset of an orbit, which is what all of today's cores equal to the
    # orbit size were saying.
    # Centred on a unit step the extra ball lands almost entirely inside the
    # first one -- 41 new points out of 7871 -- yet that alone took the core
    # from the full orbit of 15 down to 5.  Pushing the centres out to where
    # they actually add points is the obvious next turn of the same screw.
    for s in range(OFF):
        c = F.zero()
        for _ in range(OFFDIST):
            c = F.add(c, steps[1 + s * 3])
        extra = grow(F, steps, [c], DEPTH - 1, RADIUS * OFFSCALE, CAP,
                     centre=F.to_complex(c))
        before = len(coeffs)
        coeffs = list({**{c2: None for c2 in coeffs}, **{c2: None for c2 in extra}})
        print(f"  + off-centre ball {s + 1}: {before} -> {len(coeffs)}", flush=True)

    t0 = time.time()
    g, index = wire(F, coeffs, steps)
    print(f"  {g}  deg~{2 * g.m / max(g.n, 1):.2f}  [{time.time() - t0:.0f}s]",
          flush=True)

    z = F.root_of_unity(ORDER)
    witness = F.sub(F.rational(2), F.add(z, F.conj(z)))
    one = F.one()
    circle = [c for c in coeffs if F.mul(F.norm2(c), witness) == one]
    print(f"  {len(circle)} points on the magic circle", flush=True)
    if not circle:
        return

    orbits, used = [], set()
    for c in circle:
        if c in used:
            continue
        orb, x = [], c
        for _ in range(ORDER):
            orb.append(x)
            used.add(x)
            x = F.mul(x, z)
        orbits.append(orb)
    print(f"  {len(orbits)} orbits of the order-{ORDER} rotation", flush=True)

    pivot = index[F.zero()]
    for n, orb in enumerate(orbits):
        targets = [index[c] for c in orb if c in index]
        if len(targets) < 2:
            continue
        st = SeparationTest(g, K, pivot, targets)
        try:
            t0 = time.time()
            sep, core = st.run(subset=targets)
            if sep:
                print(f"  orbit {n}: {len(targets)} targets  separable"
                      f"  [{time.time() - t0:.0f}s]", flush=True)
                continue
            small = shrink(st, core)
            ok = "WITHIN CAPACITY" if len(small) <= cap_r else "too many"
            print(f"  orbit {n}: {len(targets)} targets  FORCED, minimal core "
                  f"{len(small)} vs capacity {cap_r}: {ok}"
                  f"  [{time.time() - t0:.0f}s]", flush=True)
            if len(small) <= cap_r:
                print(f"  *** {small} ***", flush=True)
                return
        finally:
            st.close()


if __name__ == "__main__":
    main()
