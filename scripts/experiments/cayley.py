"""Balls in the Cayley graph of a unit-vector set closed under differences.

The sumset sweep climbed to 6.59 edges per vertex by adding generators, which
beat the whole de Grey family -- but generic generators give a hypercube, and a
hypercube is bipartite.  What raises the chromatic number is RELATIONS between
the generators, and there is a natural source of them nobody used here.

If two unit vectors meet at sixty degrees their DIFFERENCE is also a unit
vector, since the triangle is equilateral.  So the set of unit vectors is not
arbitrary: close it under "u - v whenever |u - v| = 1" and it grows, each new
vector arriving with a relation attached.

Then take the ball of radius t in the Cayley graph of the additive group those
vectors generate: start at the origin and step by any of them, repeatedly.
Every point has an edge to every point one step away, so the degree is |U|
before coincidences and the density is |U|/2 per vertex by construction --
which makes density a dial rather than a discovery.

Reported per ball: the generator count, the points, the density, and whether
it colours with five.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
from fractions import Fraction as Fr
from collections import deque
from hn.field import Field
from hn.geometry import Point, rotation_joining, _rot60, Rotation
from hn.graph import build_graph
from pysat.solvers import Solver

k = 5
t0 = time.time()
K = Field((3, 5, 7, 11))
ONE = Point(K.one(), K.zero())
ZERO = Point(K.zero(), K.zero())


def closure(seed, cap=64):
    """Unit vectors, closed under differences that are themselves unit."""
    U = list(seed)
    seen = set(U)
    changed = True
    while changed and len(U) < cap:
        changed = False
        for a in list(U):
            for b in list(U):
                if a == b:
                    continue
                d = Point(a.x - b.x, a.y - b.y)
                if d.norm2() == 1 and d not in seen:
                    seen.add(d)
                    U.append(d)
                    changed = True
                    if len(U) >= cap:
                        break
            if len(U) >= cap:
                break
    return U


def ball(U, cap):
    """Points reachable from the origin by steps in +-U, up to the cap."""
    seen = {ZERO}
    out = [ZERO]
    q = deque([ZERO])
    steps = U + [Point(-u.x, -u.y) for u in U]
    while q and len(out) < cap:
        p = q.popleft()
        for u in steps:
            w = Point(p.x + u.x, p.y + u.y)
            if w not in seen:
                seen.add(w)
                out.append(w)
                q.append(w)
                if len(out) >= cap:
                    break
    return out


def solve(P):
    g = build_graph(P)
    E = sorted(set((min(a, b), max(a, b)) for a, b in g.edges()))
    cls = [[1 + v * k + c for c in range(k)] for v in range(g.n)]
    for a, b in E:
        for c in range(k):
            cls.append([-(1 + a * k + c), -(1 + b * k + c)])
    sv = Solver(name="cd15", bootstrap_with=cls)
    ok = sv.solve()
    sv.delete()
    return g.n, len(E), ok


w = _rot60(K)(ONE)
SEEDS = {
    "1, w": [ONE, w],
    "1, w, r4": [ONE, w, rotation_joining(Fr(4), K)(ONE)],
    "1, w, r3, r4": [ONE, w, rotation_joining(Fr(3), K)(ONE),
                     rotation_joining(Fr(4), K)(ONE)],
    "1, w, r3, r4, r7": [ONE, w, rotation_joining(Fr(3), K)(ONE),
                         rotation_joining(Fr(4), K)(ONE),
                         rotation_joining(Fr(7), K)(ONE)],
}
for label, seed in SEEDS.items():
    U = closure(seed)
    print(f"\nseed {label}: closes to {len(U)} unit vectors"
          f"  [{time.time()-t0:.0f}s]", flush=True)
    for cap in (400, 1500, 4000):
        P = ball(U, cap)
        n, e, ok = solve(P)
        print(f"   ball of {n}: {e} edges, {e/n:.2f} per vertex -> "
              f"{'5-colourable' if ok else '*** NOT 5-COLOURABLE ***'}"
              f"  [{time.time()-t0:.0f}s]", flush=True)
        if not ok:
            print("*** chi >= 6 ***", flush=True)
            sys.exit(0)
print("\nDONE", flush=True)
