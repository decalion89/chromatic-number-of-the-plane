"""Bite the bitten graph about the origin, where the answer can only be 3 or 0.

Sa u Sb keeps three of Sa's ten patterns and they form one orbit under the
sixty-degree rotation, so any operation commuting with that rotation keeps all
three or none.  Rotations about the origin commute with it.  So a second bite
about the origin has exactly two possible outcomes, and one of them is that the
union does not 4-colour -- a 5-chromatic unit-distance graph in about 1586
points by a route that is not de Grey's spindle.

Small and complete: every rational ring the field can join, both directions,
about one centre.  The colourability check is the whole test, and it is the
satisfiable answer that is fast, so a negative sweep is cheap and a positive
one announces itself.
"""
import sys, time
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
import numpy as np
from fractions import Fraction as Fr
from hn.degrey import build_Sa, build_Sb
from hn.geometry import DEGREY_FIELD as K, Point, Rotation, rotation_joining
from hn.fast import IntBasis, fast_edges_complete
from hn.homcol import closable_distance
from pysat.solvers import Solver

k = 4
t0 = time.time()
base, seen = [], set()
for p in build_Sa(K) + build_Sb(K):
    if p not in seen:
        seen.add(p)
        base.append(p)
ZERO = Point(K.zero(), K.zero())
print(f"Sa u Sb: {len(base)} points  [{time.time()-t0:.0f}s]", flush=True)

Ds = sorted({Fr(x, y) for y in range(1, 17) for x in range(1, 17 * y + 1)}
            - {Fr(1)})
Ds = [D for D in Ds if closable_distance(D)]
print(f"{len(Ds)} rings the field can join about the origin"
      f"  [{time.time()-t0:.0f}s]", flush=True)

tried = 0
for D in Ds:
    rot = rotation_joining(D, K)
    for sgn in (+1, -1):
        f = (rot if sgn > 0 else Rotation(rot.cos, -rot.sin)).about(ZERO)
        s2 = set(base)
        U = list(base) + [q for q in (f(p) for p in base) if q not in s2]
        b = IntBasis.covering(U)
        r = b.rows(U)
        if b.overflow_headroom(r) >= 1.0:
            continue
        E = sorted(set((min(a, c), max(a, c))
                       for a, c in fast_edges_complete(b, r)))
        n = len(U)
        cls = [[1 + v * k + c for c in range(k)] for v in range(n)]
        for a, c in E:
            for col in range(k):
                cls.append([-(1 + a * k + col), -(1 + c * k + col)])
        sv = Solver(name="cd15", bootstrap_with=cls)
        ok = sv.solve()
        sv.delete()
        tried += 1
        if not ok:
            print(f"*** ring D={D}, dir {sgn}: {n} points, {len(E)} edges, "
                  f"NOT 4-COLOURABLE -- 5-chromatic ***"
                  f"  [{time.time()-t0:.0f}s]", flush=True)
            sys.exit()
        if tried % 20 == 0:
            print(f"   {tried} second bites, latest D={D} -> {n} points, "
                  f"colours  [{time.time()-t0:.0f}s]", flush=True)
print(f"\n{tried} second bites about the origin, all 4-colourable"
      f"  [{time.time()-t0:.0f}s]", flush=True)
print("DONE", flush=True)
