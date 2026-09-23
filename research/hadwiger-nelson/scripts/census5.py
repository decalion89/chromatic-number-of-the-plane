"""Every rotation the class of five offers, on the 1021-point carrier.

The hub and its four forced neighbours are one colour class in every
4-colouring, so all ten pairs inside it are forced.  With the hub as pivot all
four radii are equal and the free angle collapses to the classical spindle,
but with a NEIGHBOUR as pivot the distances to the hub and to the other
neighbours differ, and there the free angle is a genuinely different rotation.

This enumerates, for every ordered (pivot, x, y) inside the class: whether the
circles meet, what radical the free angle needs, and whether the classical
spindle at |pivot - x| is available too.  The point is to find the largest set
of independent rotations the SMALL carrier supports, since its one-pivot object
is the only one whose five-colour test can be asked.
"""
import sys, time, random
from fractions import Fraction as Fr
from collections import defaultdict
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
from hn.degrey import build_Sa
from hn.geometry import Point, rotation_joining, _rot60
from hn.field import Field
from hn.graph import build_graph
from pysat.solvers import Solver

F = Field((3, 11, 247)); rot60 = _rot60(F); Sa = build_Sa(F)
def orb(p):
    out, q = [], p
    for _ in range(6):
        if q not in out: out.append(q)
        q = rot60(q)
    for q in list(out):
        r = Point(q.x, -q.y)
        if r not in out: out.append(r)
    return out
t0 = time.time()
seen, U = set(Sa), list(Sa)
for w in orb(Sa[25]):
    rot = rotation_joining(Fr(1), F).about(w)
    for p in Sa:
        q = rot(p)
        if q not in seen: seen.add(q); U.append(q)
g = build_graph(U); n = g.n
print(f"carrier n={n}   [{time.time()-t0:.0f}s]", flush=True)
K = 4
X = lambda v, c: 1 + v * K + c
cnf = [[X(v, c) for c in range(K)] for v in range(n)]
for u, v in g.edges():
    for c in range(K):
        cnf.append([-X(u, c), -X(v, c)])
HUB = 416
s = Solver(name="m22", bootstrap_with=cnf)
cls = [HUB]
for v in range(n):
    if v == HUB: continue
    d = (U[HUB] - U[v]).norm2()
    if not d.is_rational() or Fr(d.c[0]) != Fr(64, 3): continue
    if not s.solve(assumptions=[X(HUB, 0), X(v, 1)]):
        cls.append(v)
s.delete()
print(f"  class of {len(cls)}: {cls}   [{time.time()-t0:.0f}s]", flush=True)
def sqfree(x):
    r, d = x, 2
    while d * d <= r:
        while r % (d * d) == 0: r //= d * d
        d += 1
    return r
FREE = {1, 3, 11, 33, 247, 741, 2717, 8151}
rots = []
for i, piv in enumerate(cls):
    ds = {}
    for v in cls:
        if v == piv: continue
        dd = (U[piv] - U[v]).norm2()
        ds[v] = Fr(dd.c[0]) if dd.is_rational() else None
    # the classical spindle, one per distinct radius
    for v, q in ds.items():
        if q is None or q <= Fr(1, 4): continue
        c = Fr(1) - Fr(1, 2) / q
        s2 = 1 - c * c
        rad = 1 if s2 == 0 else sqfree(s2.numerator * s2.denominator)
        rots.append(("classical", piv, v, None, q, rad))
    # the free angle, one ordered pair at a time
    for v, q1 in ds.items():
        for w, q2 in ds.items():
            if v >= w or q1 is None or q2 is None or q1 == q2: continue
            lo = abs(float(q1) ** .5 - float(q2) ** .5); hi = float(q1) ** .5 + float(q2) ** .5
            if not (lo <= 1 <= hi): continue
            R = (q1 + q2 - 1) / 2
            disc = q1 * q2 - R * R
            if disc <= 0: continue
            rad = sqfree(disc.numerator * disc.denominator)
            rots.append(("free", piv, v, w, (q1, q2), rad))
byrad = defaultdict(int)
for r in rots: byrad[r[5]] += 1
print(f"  {len(rots)} rotations available; radicals {dict(byrad)}", flush=True)
freeones = [r for r in rots if r[5] in FREE]
print(f"  of which {len(freeones)} need no new radical", flush=True)
for r in rots[:14]:
    kind, piv, v, w, q, rad = r
    tag = "FREE" if rad in FREE else f"needs sqrt{rad}"
    print(f"    {kind:<9} pivot v{piv:<5} {v}"
          + (f",{w}" if w is not None else "")
          + f"  r={q}  radical {rad}  {tag}", flush=True)
