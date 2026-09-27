"""Bite from somewhere other than the origin, and watch the orbit break.

The three surviving patterns of Sa u rho(Sa) are one orbit under sixty-degree
rotation about the ORIGIN, which is why the census cannot read one or two.
Every operation used so far -- the bite, the stacks, the high-contact unions --
rotates about that same point and so preserves the orbit.  A single added
point does not break it either: twenty-five candidates with up to twelve
neighbours each left the census at three.

Rotating about a DIFFERENT centre breaks the symmetry wholesale, and the union
is still a unit-distance graph.  If the census then reads one, the centre and
its three antipodal pairs are pinned up to renaming colours -- strictly more
forced structure than de Grey's chain uses.  If it reads zero, the union is
5-chromatic, and its size is worth knowing against G's 1581.

Centres are taken from the graph itself, rotations from the rings the field
can join, and the three patterns are the whole test.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time, math
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
import numpy as np
from fractions import Fraction as Fr
from collections import defaultdict
from hn.degrey import build_Sa, build_Sb
from hn.geometry import DEGREY_FIELD as K, Point, Rotation, rotation_joining
from hn.fast import IntBasis, fast_edges_complete
from hn.homcol import closable_distance
from pysat.solvers import Solver

k = 4
t0 = time.time()
THREE = [
    [[0, 1, 4], [2, 3, 5, 6]],
    [[1, 2, 4, 5], [0, 3, 6]],
    [[0, 2, 5], [1, 3, 4, 6]],
]
ZERO = Point(K.zero(), K.zero())
U, seen = [], set()
for p in build_Sa(K) + build_Sb(K):
    if p not in seen:
        seen.add(p)
        U.append(p)
b0 = IntBasis.covering(U)
r0 = b0.rows(U)
E0 = set((min(a, c), max(a, c)) for a, c in fast_edges_complete(b0, r0))
deg = defaultdict(int)
for a, c in E0:
    deg[a] += 1
    deg[c] += 1
print(f"Sa u Sb: {len(U)} pts, {len(E0)} edges  [{time.time()-t0:.0f}s]",
      flush=True)

Ds = [D for D in (Fr(x, y) for y in range(1, 13) for x in range(1, 13 * y + 1))
      if D != 1 and closable_distance(D)]
Ds = sorted(set(Ds))
centres = sorted(range(len(U)), key=lambda i: -deg[i])[:6]
centres = [c for c in centres if U[c] != ZERO][:4]
print(f"{len(Ds)} rings, centres {centres} with degrees "
      f"{[deg[c] for c in centres]}  [{time.time()-t0:.0f}s]", flush=True)


def seven(P):
    b = IntBasis.covering(P)
    r = b.rows(P)
    if b.overflow_headroom(r) >= 1.0:
        return None
    dm, d2 = b.dim, b.D * b.D
    E = sorted(set((min(a, c), max(a, c))
                   for a, c in fast_edges_complete(b, r)))
    ci = P.index(ZERO)
    key = {tuple(r[i]): i for i in range(len(P))}
    d = r - r[ci]
    sq = b._field_square(d[:, :dm]) + b._field_square(d[:, dm:])
    ok = np.ones(len(sq), dtype=bool)
    for j in range(1, dm):
        ok &= sq[:, j] == 0
    ring = [int(o) for o in np.nonzero(ok)[0] if Fr(int(sq[o, 0]), d2) == 4]
    ms = set(ring)
    pairs = [(j, key[tuple(2 * r[ci] - r[j])]) for j in ring
             if tuple(2 * r[ci] - r[j]) in key
             and key[tuple(2 * r[ci] - r[j])] in ms
             and key[tuple(2 * r[ci] - r[j])] > j]
    if len(pairs) < 3:
        return None
    six = [x for pr in pairs[:3] for x in pr]
    six.sort(key=lambda i: math.atan2(float(P[i].y), float(P[i].x)))
    return len(P), E, [ci] + six


tried = 0
for ci in centres:
    C = U[ci]
    for D in Ds:
        rot = rotation_joining(D, K)
        for sgn in (+1, -1):
            f = (rot if sgn > 0 else Rotation(rot.cos, -rot.sin)).about(C)
            s2 = set(U)
            V = list(U) + [q for q in (f(p) for p in U) if q not in s2]
            st = seven(V)
            if st is None:
                continue
            n, E, W = st
            cls = [[1 + v * k + c for c in range(k)] for v in range(n)]
            for a, c in E:
                for col in range(k):
                    cls.append([-(1 + a * k + col), -(1 + c * k + col)])
            sv = Solver(name="cd15", bootstrap_with=cls)
            tried += 1
            if not sv.solve():
                print(f"*** centre {ci}, ring D={D}, dir {sgn}: {n} points "
                      f"NOT 4-COLOURABLE ***  [{time.time()-t0:.0f}s]",
                      flush=True)
                sv.delete()
                sys.exit()
            keep = sum(1 for part in THREE
                       if sv.solve(assumptions=[1 + W[e] * k + bi
                                                for bi, blk in enumerate(part)
                                                for e in blk]))
            sv.delete()
            if keep < 3:
                print(f"*** centre {ci}, ring D={D}, dir {sgn}: {n} pts, "
                      f"census {keep} of 3 -- ORBIT BROKEN ***"
                      f"  [{time.time()-t0:.0f}s]", flush=True)
            if tried % 8 == 0:
                print(f"   {tried} tried, latest centre {ci} D={D}: {n} pts, "
                      f"census {keep}  [{time.time()-t0:.0f}s]", flush=True)
print(f"\n{tried} off-centre bites  [{time.time()-t0:.0f}s]", flush=True)
print("DONE", flush=True)
