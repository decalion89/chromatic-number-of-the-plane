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
from itertools import combinations

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from hn.cyclograph import scaled_graph
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

    # Several centres, because a ball a rotation preserves never forces a
    # proper subset of a target orbit -- which is what every core equal to the
    # orbit size was saying. Off-centre growth took the core from 15 to 5.
    seeds, centres = [F.zero()], [0j]
    for s in range(OFF):
        c = F.zero()
        for _ in range(OFFDIST):
            c = F.add(c, steps[1 + s * 3])
        seeds.append(c)
        centres.append(F.to_complex(c))

    t0 = time.time()
    g = scaled_graph(F, steps, DEPTH, RADIUS, cap=CAP, seeds=seeds,
                     centres=centres)
    index = {p.c: i for i, p in enumerate(g.vertices)}
    print(f"  {g}  deg~{2 * g.m / max(g.n, 1):.2f}  [{time.time() - t0:.0f}s]",
          flush=True)

    z = F.root_of_unity(ORDER)
    witness = F.sub(F.rational(2), F.add(z, F.conj(z)))
    one = F.one()
    circle = [p.c for p in g.vertices if F.mul(F.norm2(p.c), witness) == one]
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
        if all(o in index for o in orb):
            orbits.append(orb)
    print(f"  {len(orbits)} complete orbits of the order-{ORDER} rotation",
          flush=True)

    blocking = blocking_subsets(ORDER, cap_r)
    print(f"  {len(blocking)} blocking {cap_r}-subsets of C_{ORDER}", flush=True)

    pivot = index[F.zero()]
    for n, orb in enumerate(orbits):
        idx = [index[c] for c in orb]
        st = SeparationTest(g, K, pivot, idx)
        try:
            t0 = time.time()
            sep, core = st.run(subset=idx)
            if sep:
                print(f"  orbit {n}: separable  [{time.time() - t0:.0f}s]",
                      flush=True)
                continue
            small = shrink(st, core)
            pos = {t: k for k, t in enumerate(idx)}
            arc = sorted(pos[t] for t in small)
            print(f"  orbit {n}: FORCED, minimal core {len(small)} at {arc}"
                  f" vs capacity {cap_r}  [{time.time() - t0:.0f}s]", flush=True)
            # a minimal core is not unique, so ask the blocking sets directly
            for T in blocking:
                if not st.run(subset=[idx[t] for t in T])[0]:
                    print(f"  *** FORCED on blocking subset {T} ***", flush=True)
                    return
            print(f"  none of the {len(blocking)} blocking subsets is forced",
                  flush=True)
        finally:
            st.close()


def blocking_subsets(n: int, r: int):
    """The r-subsets of C_n that no independent set can meet in every shift.

    The copies are the n rotations, so their image sets are the n shifts of T,
    and a colouring escapes through an independent S with S - T = Z_n.  For
    C_15 and r = 4 exactly 45 of the 1365 subsets survive, in three classes up
    to rotation.  A minimal core is not unique, so these get asked directly
    rather than hoped for from shrinking.
    """
    indep = [()]
    for size in range(1, n // 2 + 1):
        for S in combinations(range(n), size):
            if all((a - b) % n not in (1, n - 1) for a in S for b in S if a != b):
                indep.append(S)
    full = set(range(n))
    return [T for T in combinations(range(n), r)
            if not any({(a - b) % n for a in S for b in T} == full
                       for S in indep)]


if __name__ == "__main__":
    main()
