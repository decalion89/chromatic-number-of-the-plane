"""What are the five vertices that force all four colours?

Unioning rotated copies of Sa drops rho from 7 to 5 and holds it there, one
above the floor rho >= k.  Five vertices out of 619 such that EVERY proper
4-colouring puts all four colours on them is extreme rigidity, and the shape
of that set is the thing to copy at five colours.

So: what do they induce, how far apart are they, and is one of them a pivot
with the others on its circle?
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
from fractions import Fraction
from hn import degrey
from hn.geometry import Rotation
from hn.graph import build_graph
from hn.forced import ColourRelations, shrink_forcing_set
from pysat.solvers import Solver

FLD = degrey.DEGREY_FIELD
ROT = Rotation(FLD.sqrt(33) * FLD.rational(Fraction(1, 6)),
               FLD.sqrt(3) * FLD.rational(Fraction(1, 6)))

pts = list(degrey.build_Sa())
seen = set(pts)
for p in list(pts):
    q = ROT(p)
    if q not in seen:
        seen.add(q)
        pts.append(q)
g = build_graph(pts)
print(f"Sa u rot(Sa): {g.n} vertices, {g.m} edges", flush=True)

k = 4
NV = g.n * k


def x(v, c):
    return 1 + v * k + c


def sel(v):
    return NV + 1 + v


cls = [[x(v, c) for c in range(k)] for v in range(g.n)]
for a, b in g.edges():
    for c in range(k):
        cls.append([-x(a, c), -x(b, c)])
for v in range(g.n):
    cls.append([-sel(v), -x(v, 0)])

S = list(range(g.n))
for _ in range(14):
    with Solver(name="cd19", bootstrap_with=cls) as s:
        if s.solve(assumptions=[sel(v) for v in S]):
            break
        core = sorted({abs(l) - NV - 1 for l in s.get_core()})
    if len(core) >= len(S):
        break
    S = core
S = shrink_forcing_set(ColourRelations(g, k), S)
print(f"  minimal forcing set: {len(S)} vertices {sorted(S)}", flush=True)

deg = g.degrees()
print(f"  degrees: {[deg[v] for v in S]}", flush=True)
sub = g.induced(sorted(S))
print(f"  induced: {sub.n} vertices, {sub.m} edges", flush=True)
import itertools
print("  pairwise squared distances:")
for a, b in itertools.combinations(sorted(S), 2):
    d = g.vertices[a].dist2(g.vertices[b])
    print(f"    {a:4} {b:4}  {float(d):8.4f}"
          + ("   (unit)" if d == FLD.one() else ""), flush=True)
# is any of them a pivot with the rest on its circle?
for v in sorted(S):
    on = [u for u in S if u != v and u in g.adj[v]]
    if len(on) >= 2:
        print(f"  vertex {v} has {len(on)} of the others on its unit circle",
              flush=True)
