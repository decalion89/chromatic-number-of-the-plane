"""Independent check of rho = 5 on Sa u rot(Sa).

Five vertices out of 619 forcing all four colours is the headline the strategy
now rests on, and it came from the same code that computed it.  So: rebuild
the CNF from scratch, and for each of the four colours ask directly whether a
proper 4-colouring exists leaving that colour off the set.  Then confirm
minimality by putting each vertex back.

Also checks the thing that makes rho meaningful at all -- that the union is
still 4-colourable, so the question is not vacuous.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
import sys, time
sys.path.insert(0, HN_DIR)
from fractions import Fraction
from hn import degrey
from hn.geometry import Rotation
from hn.graph import build_graph
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
S = [107, 208, 502, 568, 618]
k = 4


def x(v, c):
    return 1 + v * k + c


base = [[x(v, c) for c in range(k)] for v in range(g.n)]
for a, b in g.edges():
    for c in range(k):
        base.append([-x(a, c), -x(b, c)])

print(f"Sa u rot(Sa): {g.n} vertices, {g.m} edges", flush=True)
with Solver(name="cd19", bootstrap_with=base) as s:
    print(f"  4-colourable at all: {s.solve()}  (else rho is vacuous)",
          flush=True)

ok = True
for c in range(k):
    with Solver(name="cd19",
                bootstrap_with=base + [[-x(v, c)] for v in S]) as s:
        r = s.solve()
    print(f"  a 4-colouring leaving colour {c} off S: {r}", flush=True)
    ok &= not r
print(f"  => S is {'FORCING' if ok else 'NOT forcing'}", flush=True)

broken = 0
for drop in S:
    T = [v for v in S if v != drop]
    if any(Solver(name="cd19",
                  bootstrap_with=base + [[-x(v, c)] for v in T]).solve()
           for c in range(k)):
        broken += 1
print(f"  minimal: {broken} of {len(S)} single deletions break it "
      f"({'all' if broken == len(S) else 'NOT all'})", flush=True)

sub = g.induced(sorted(S))
for kk in (1, 2, 3, 4):
    cls = [[1 + i * kk + c for c in range(kk)] for i in range(sub.n)]
    for a, b in sub.edges():
        for c in range(kk):
            cls.append([-(1 + a * kk + c), -(1 + b * kk + c)])
    with Solver(name="cd19", bootstrap_with=cls) as s:
        if s.solve():
            print(f"  the five induce {sub.m} edges and are {kk}-chromatic "
                  f"-- the forcing is ambient", flush=True)
            break
