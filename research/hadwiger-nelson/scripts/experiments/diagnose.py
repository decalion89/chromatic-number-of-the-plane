"""Why is the minimum hitting set 2?  Look at the classes being generated.

A lower bound of 2 on rho is below the floor rho >= k = 5, so either the
family is degenerate or something in the loop is wrong.  The obvious suspect
is the SIZE of the colour-0 classes the solver returns: asked only to avoid a
63-vertex set, it is free to make that class as small as it likes, and tiny
classes are trivially hit.

Earlier sampling, which asked for classes disjoint from a growing union, found
sizes 267 to 326.  This asks what the escaping-class loop actually produces.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time, collections
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
from pysat.card import CardEnc, EncType
from pysat.formula import IDPool
from pysat.solvers import Solver
from hn import degrey

g = degrey.build_G()
n, k, BUD = g.n, 5, 63


def x(v, c):
    return 1 + v * k + c


cls = [[x(v, c) for c in range(k)] for v in range(n)]
for a, b in g.edges():
    for c in range(k):
        cls.append([-x(a, c), -x(b, c)])
colour = Solver(name="cd19", bootstrap_with=cls)
pool = IDPool(start_from=n + 1)
card = CardEnc.atmost(lits=list(range(1, n + 1)), bound=BUD, vpool=pool,
                      encoding=EncType.seqcounter)
hit = Solver(name="cd19", bootstrap_with=card.clauses)

sizes, common = [], None
for rnd in range(400):
    if not hit.solve():
        print("no hitting set")
        break
    m = hit.get_model()
    S = [v for v in range(n) if m[v] > 0]
    if not colour.solve(assumptions=[-x(v, 0) for v in S]):
        print(f"forcing at round {rnd}")
        break
    mm = set(colour.get_model())
    c0 = frozenset(v for v in range(n) if x(v, 0) in mm)
    sizes.append(len(c0))
    common = c0 if common is None else (common & c0)
    hit.add_clause([v + 1 for v in c0])
colour.delete()
hit.delete()

print(f"{len(sizes)} escaping classes generated")
print(f"  sizes: min {min(sizes)}, max {max(sizes)}, "
      f"median {sorted(sizes)[len(sizes)//2]}")
print(f"  vertices common to ALL of them: {len(common)}  {sorted(common)[:10]}")
print("\nA colour class of a proper 5-colouring of a 5-chromatic graph is")
print("never empty, but nothing stops it being SMALL -- and a family of")
print("classes sharing a couple of vertices is hit by a couple of vertices.")
