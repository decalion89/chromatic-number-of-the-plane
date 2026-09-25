"""Second solver, same question: is G - 1420 still 5-chromatic?

CaDiCaL reported UNSAT in 1581 seconds, refuting a claim this package leaned
on throughout.  A refutation of my own assertion deserves more than one
solver's word, and drat-trim is not installed here, so the cross-check is a
different engine on a formula rebuilt from scratch -- Glucose rather than
CaDiCaL, edges re-derived from the graph rather than reused.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
from hn import degrey
from pysat.solvers import Solver

g = degrey.build_G()
v = 1420
keep = [u for u in range(g.n) if u != v]
sub = g.induced(keep)
k = 4
cls = [[1 + i * k + c for c in range(k)] for i in range(sub.n)]
for a, b in sub.edges():
    for c in range(k):
        cls.append([-(1 + a * k + c), -(1 + b * k + c)])
tri = sub.find_clique(3)
for i, u in enumerate(tri):
    cls.append([1 + u * k + i])
print(f"G - {v}: {sub.n} vertices, {sub.m} edges, triangle {tri} pinned",
      flush=True)
t = time.time()
with Solver(name="g4", bootstrap_with=cls) as s:
    r = s.solve()
print(f"  Glucose says 4-colourable: {r}  [{time.time()-t:.0f}s]", flush=True)
print("  CaDiCaL said False in 1581s; the two agree"
      if r is False else "  THE TWO DISAGREE -- investigate before claiming",
      flush=True)
