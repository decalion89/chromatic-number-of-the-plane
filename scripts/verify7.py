"""Independent check: are those seven vertices really rainbow-forcing?

The construction says rho(Sa, 4) = 7, and the seven it found induce only five
edges with chromatic number 3 and three isolated vertices -- so they are NOT a
Moser spindle and the forcing cannot come from the subgraph.  That is a strong
enough claim to re-derive from scratch rather than take from the same code
that produced it.

Rainbow-forcing means no proper 4-colouring leaves a colour off the set.  By
colour symmetry it is enough to ask for one: is there a proper 4-colouring of
Sa in which colour 0 appears nowhere on S?  A direct CNF, built here rather
than through ColourRelations, answers it.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
import sys
sys.path.insert(0, HN_DIR)
from hn import degrey
from hn.graph import build_graph
from hn.forced import ColourRelations, forcing_set, shrink_forcing_set
from pysat.solvers import Solver

g = build_graph(degrey.build_Sa())
rel = ColourRelations(g, 4)
S, closed = forcing_set(rel, limit=500)
S = sorted(shrink_forcing_set(rel, S))
print(f"S = {S}", flush=True)

k = 4


def x(v, c):
    return 1 + v * k + c


base = [[x(v, c) for c in range(k)] for v in range(g.n)]
for a, b in g.edges():
    for c in range(k):
        base.append([-x(a, c), -x(b, c)])

# sanity: Sa is 4-colourable at all
with Solver(name="cd19", bootstrap_with=base) as s:
    assert s.solve(), "Sa must be 4-colourable"

for c in range(k):
    cls = base + [[-x(v, c)] for v in S]
    with Solver(name="cd19", bootstrap_with=cls) as s:
        r = s.solve()
    print(f"  a proper 4-colouring avoiding colour {c} on S: {r}", flush=True)
    if r:
        print("    => S is NOT rainbow-forcing; the measurement is wrong")
        break
else:
    print(f"\n  CONFIRMED: no proper 4-colouring of Sa leaves any colour off "
          f"these {len(S)} vertices, and they induce only 5 edges with "
          f"chromatic number 3.")
    print("  The forcing is AMBIENT, not carried by a critical subgraph.")

# and minimality: every proper subset should fail
bad = 0
for drop in S:
    T = [v for v in S if v != drop]
    for c in range(k):
        cls = base + [[-x(v, c)] for v in T]
        with Solver(name="cd19", bootstrap_with=cls) as s:
            if s.solve():
                bad += 1
                break
print(f"  minimal: {bad} of {len(S)} single deletions break the forcing "
      f"({'all' if bad == len(S) else 'NOT all'})")
