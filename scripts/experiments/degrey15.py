#!/usr/bin/env python3
"""de Grey's rotations and the order-15 circle, in one field.

`hn.quadext.degrey_field` is Q(zeta_15, sqrt(-7), sqrt(-11)), degree 32. Read
as complex numbers de Grey's three rotations need only sqrt(-15), sqrt(-7) and
sqrt(-11), and sqrt(-15) = sqrt(-3) sqrt(5) is already cyclotomic -- so his
whole 2018 construction fits there, beside zeta_15.

That matters because the two halves of the argument had never been able to
share a field. His construction forces at k = 4, which nothing built here
does. The order-15 magic circle has capacity 4, where every multiquadratic
field is capped at 2 by n | 24. Both, now, in degree 32.

So: grow with his rotations as unit steps, find the order-15 circle, and ask
whether the pivot's colour is forced onto at most four of one orbit.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import os
import sys
import time
from fractions import Fraction
from itertools import combinations

sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]

from hn.cyclograph import scaled_graph
from hn.quadext import degrey_field, degrey_rotations
from hn.spindle import SeparationTest
from hn.transversal import trapping_bound

K = int(os.environ.get("HN_K", "4"))
ORDER = int(os.environ.get("HN_ORDER", "15"))
T_STEP = int(os.environ.get("HN_T", "1"))
DEPTH = int(os.environ.get("HN_DEPTH", "3"))
RADIUS = float(os.environ.get("HN_RADIUS", "3.0"))
CAP = int(os.environ.get("HN_CAP", "40000"))
OFF = int(os.environ.get("HN_OFF", "2"))
MULT = int(os.environ.get("HN_MULT", "1"))
WHICH = [w for w in os.environ.get("HN_ROTS", "").split(",") if w] or None


def roots(L):
    base = L.base.base
    return [L.embed(L.base.embed(base.zeta(k))) for k in range(base.n)]


def unit_steps(L, mult: int, which=None):
    """Roots of unity times short products of de Grey's rotations.

    Folding in all four at once gives 270 steps, and then a cap truncates the
    breadth-first frontier inside the first layer: 25000 points came out at
    average degree 6.43 with *no* circle points at all, because the circle sits
    two or three steps out. Fewer rotations, more depth.
    """
    allr = degrey_rotations(L)
    rots = [allr[k] for k in (which or list(allr))]
    steps = list(roots(L))
    cur = list(steps)
    for _ in range(mult):
        grown = list(cur)
        for r in rots:
            grown += [L.mul(s, r) for s in cur]
            grown += [L.mul(s, L.conj(r)) for s in cur]
        cur = grown
    steps = cur + [L.neg(s) for s in cur]
    return [s for s in {t: None for t in steps} if L.norm2(s) == L.one()]


def blocking_subsets(n: int, r: int):
    indep = [()]
    for size in range(1, n // 2 + 1):
        for S in combinations(range(n), size):
            if all((a - b) % n not in (1, n - 1) for a in S for b in S if a != b):
                indep.append(S)
    full = set(range(n))
    return [T for T in combinations(range(n), r)
            if not any({(a - b) % n for a in S for b in T} == full for S in indep)]


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
    L = degrey_field()
    steps = unit_steps(L, MULT, WHICH)
    print(f"degree {L.degree}, k={K}, order {ORDER}, t={T_STEP}: "
          f"{len(steps)} unit steps", flush=True)

    seeds, centres = [L.zero()], [0j]
    for s in range(OFF):
        c = L.zero()
        for _ in range(3):
            c = L.add(c, steps[1 + s * 3])
        seeds.append(c)
        centres.append(L.to_complex(c))

    t0 = time.time()
    g = scaled_graph(L, steps, DEPTH, RADIUS, cap=CAP, seeds=seeds, centres=centres)
    index = {p.c: i for i, p in enumerate(g.vertices)}
    print(f"  {g}  deg~{2 * g.m / max(g.n, 1):.2f}  [{time.time() - t0:.0f}s]",
          flush=True)

    base = L.base.base
    z = L.embed(L.base.embed(base.zeta(base.n // ORDER)))
    zt = L.embed(L.base.embed(base.zeta((base.n // ORDER) * T_STEP)))
    witness = L.sub(L.rational(2), L.add(zt, L.conj(zt)))
    one = L.one()
    circle = [p.c for p in g.vertices if L.mul(L.norm2(p.c), witness) == one]
    print(f"  {len(circle)} points on the order-{ORDER} circle (t={T_STEP})",
          flush=True)
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
            x = L.mul(x, z)
        if all(o in index for o in orb):
            orbits.append(orb)
    from math import gcd
    cyc = ORDER // gcd(ORDER, T_STEP)
    cap_r = trapping_bound(cyc)
    blocking = [tuple((p * T_STEP) % ORDER for p in T)
                for T in blocking_subsets(cyc, cap_r)]
    print(f"  {len(orbits)} complete orbits | C_{cyc}, capacity {cap_r}, "
          f"{len(blocking)} blocking subsets", flush=True)

    pivot = index[L.zero()]
    for n, orb in enumerate(orbits):
        idx = [index[c] for c in orb]
        st = SeparationTest(g, K, pivot, idx)
        try:
            t0 = time.time()
            sep, core = st.run(subset=idx)
            if sep:
                print(f"  orbit {n}: separable  [{time.time() - t0:.0f}s]", flush=True)
                continue
            small = shrink(st, core)
            pos = {t: k for k, t in enumerate(idx)}
            print(f"  orbit {n}: FORCED, minimal core {len(small)} at "
                  f"{sorted(pos[t] for t in small)} vs capacity {cap_r}"
                  f"  [{time.time() - t0:.0f}s]", flush=True)
            for T in blocking:
                if not st.run(subset=[idx[t] for t in T])[0]:
                    print(f"  *** FORCED on blocking subset {T} ***", flush=True)
                    return
            print(f"  none of the {len(blocking)} blocking subsets is forced",
                  flush=True)
        finally:
            st.close()


if __name__ == "__main__":
    main()
