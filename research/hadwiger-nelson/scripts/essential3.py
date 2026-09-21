"""How much of G can be peeled at all?  Exhaustively, at full scale.

Thirty vertices were probed once before: nineteen essential -- G - v is
4-COLOURABLE, so v lies in every 5-critical subgraph and cannot be removed --
none dispensable, eleven undecided.  That sampled two per cent.

Exhibiting a colouring of G - v proves v essential outright; failing inside a
budget proves nothing either way.  So a budgeted pass is sound for every
vertex it decides, and it is the satisfiable direction, which is the fast one.

Two encodings were tried and discarded before this one.  Rebuilding the
instance per vertex costs eighteen seconds in clause construction alone.
Guarding every vertex with a selector AND forbidding an absent vertex to carry
a colour slows the solver badly -- five of fifty decided at 25000 conflicts.

The right encoding is the lighter half of the second.  Guard only the
"some colour" clause: an absent vertex is then free to carry no colour at all,
and a colouring of G - v extends to exactly that, while any solution restricts
to a proper colouring of G - v.  The relaxation IS G - v, no extra clauses
needed, and one solver keeps everything it learns across all 1581 calls.
"""
import sys, time, json
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
import numpy as np
from hn.degrey import build_G
from hn.geometry import DEGREY_FIELD as K
from hn.fast import IntBasis, fast_edges_complete
from pysat.solvers import Solver

k = 4
BUDGET = int(sys.argv[1]) if len(sys.argv) > 1 else 1000000
t0 = time.time()
P = build_G(K, as_graph=False)
n = len(P)
b = IntBasis.covering(P)
r = b.rows(P)
assert b.overflow_headroom(r) < 1.0
E = sorted(set((min(a, c), max(a, c)) for a, c in fast_edges_complete(b, r)))
print(f"G: {n} pts, {len(E)} edges; budget {BUDGET} conflicts"
      f"  [{time.time()-t0:.0f}s]", flush=True)

nv = n * k
sel = [nv + 1 + v for v in range(n)]
cls = [[-sel[v]] + [1 + v * k + c for c in range(k)] for v in range(n)]
for a, c in E:
    for col in range(k):
        cls.append([-(1 + a * k + col), -(1 + c * k + col)])
sv = Solver(name="cd15", bootstrap_with=cls)
print(f"solver built, {len(cls)} clauses  [{time.time()-t0:.0f}s]", flush=True)
assume_all = [sel[u] for u in range(n)]
# G itself is not 4-colourable -- that is de Grey's theorem, and the
# repository verifies it in its own slow test.  Re-proving it here would spend
# the entire budget before the census starts, so it is taken as given and the
# per-vertex answers stand on their own: each is a COLOURING, exhibited and
# checked edge by edge below.

essential, undecided = [], []
for v in range(n):
    a = list(assume_all)
    a[v] = -sel[v]
    sv.conf_budget(BUDGET)
    res = sv.solve_limited(assumptions=a)
    if res is True:
        mdl = set(x for x in sv.get_model() if x > 0)
        bad = [(x, y) for x, y in E if x != v and y != v
               and any((1 + x * k + col) in mdl and (1 + y * k + col) in mdl
                       for col in range(k))]
        assert not bad, f"claimed colouring of G-{v} has a monochromatic edge"
        assert all(any((1 + u * k + col) in mdl for col in range(k))
                   for u in range(n) if u != v), "a vertex went uncoloured"
        essential.append(v)
    else:
        undecided.append(v)
    if (v + 1) % 20 == 0:
        print(f"   {v+1}/{n}: {len(essential)} essential, {len(undecided)} "
              f"undecided  [{time.time()-t0:.0f}s]", flush=True)
        json.dump({"budget": BUDGET, "done": v + 1,
                   "essential": essential, "undecided": undecided},
                  open("essential.json", "w"))
print(f"\n{len(essential)} of {n} PROVED essential "
      f"({100*len(essential)/n:.1f} per cent), {len(undecided)} undecided"
      f"  [{time.time()-t0:.0f}s]", flush=True)
print("DONE", flush=True)
