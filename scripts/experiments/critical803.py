"""Is the 803-point graph vertex-critical?  The cheap direction answers it.

Batch removal stopped dead: dropping the single lowest-degree vertex already
makes the graph 4-colourable.  That is the signature of a vertex-critical
graph -- one where every vertex is essential and no deletion can shrink it.

Checking it in full looks like 803 UNSAT proofs at 47 s each, which is ten
hours.  It is not.  Asking "does H - v still refuse four?" is asking the solver
to REFUTE H - v, and the expensive direction is exactly the one that does not
happen if the graph is critical: if v is essential, H - v is 4-COLOURABLE, and
finding a colouring is fast.  Only a removable vertex costs a slow refutation,
and by hypothesis there are none.  So the full check is cheap precisely when
the answer is yes.

A triangle is pinned in every call.  Deciding colourability is invariant under
permuting colours, so fixing one clique to 0,1,2 is free, and it deletes the
4! symmetric copies of every refutation -- the difference between 47 s and
more than five minutes on this graph.  The triangle is recomputed for each
H - v, since v may be in it.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time, json
from fractions import Fraction as Fr
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
from hn.field import Field
from hn.geometry import Point
from hn.graph import build_graph
from pysat.solvers import Solver

ROOT = HN_DIR
t0 = time.time()
F = Field((3, 11, 247))
d = json.load(open(f"{ROOT}/data/five_247_c.json"))
PTS = [Point(F.element([Fr(a, b) for a, b in x]),
             F.element([Fr(a, b) for a, b in y])) for x, y in d["points"]]
K = 4

def four_colourable(sub):
    g = build_graph([PTS[i] for i in sub]); n = g.n
    X = lambda v, c: 1 + v * K + c
    cl = [[X(v, c) for c in range(K)] for v in range(n)]
    for v in range(n):
        for a in range(K):
            for b in range(a + 1, K):
                cl.append([-X(v, a), -X(v, b)])
    for x, y in g.edges():
        for c in range(K):
            cl.append([-X(x, c), -X(y, c)])
    for i, x in enumerate(g.find_clique(3) or []):
        cl.append([X(x, i)])
        for c in range(K):
            if c != i:
                cl.append([-X(x, c)])
    s = Solver(name="cd19", bootstrap_with=cl); r = s.solve(); s.delete()
    return r

full = list(range(len(PTS)))
g = build_graph(PTS)
print(f"n={g.n} m={sum(len(a) for a in g.adj)//2}", flush=True)
removable, slow = [], []
for v in full:
    t1 = time.time()
    ok = four_colourable([u for u in full if u != v])
    dt = time.time() - t1
    if not ok:
        removable.append(v)
        print(f"  *** v{v} is REMOVABLE (H-v still refuses four) "
              f"[{dt:.0f}s] ***", flush=True)
    if dt > 5:
        slow.append((v, round(dt, 1)))
    if v % 50 == 0:
        print(f"    ..{v}/{g.n}  {len(removable)} removable"
              f"   [{time.time()-t0:.0f}s]", flush=True)
print(f"  {len(removable)} removable vertices of {g.n}   "
      f"[{time.time()-t0:.0f}s]", flush=True)
if not removable:
    print("  *** VERTEX-CRITICAL: every vertex is essential ***", flush=True)
    json.dump({"vertex_critical": True, "n": g.n,
               "m": sum(len(a) for a in g.adj) // 2,
               "slowest_calls": sorted(slow, key=lambda t: -t[1])[:10]},
              open(f"{ROOT}/data/five_247_c_critical.json", "w"))
else:
    keep = [u for u in full if u not in set(removable)]
    gg = build_graph([PTS[i] for i in keep])
    mm = sum(len(a) for a in gg.adj) // 2
    print(f"  after one pass: n={gg.n} m={mm}", flush=True)
    json.dump({"field_generators": list(F.gens), "n": gg.n, "m": mm,
               "points": [[[[c.numerator, c.denominator] for c in p.x.c],
                           [[c.numerator, c.denominator] for c in p.y.c]]
                          for p in gg.vertices]},
              open(f"{ROOT}/data/five_247_min.json", "w"))
