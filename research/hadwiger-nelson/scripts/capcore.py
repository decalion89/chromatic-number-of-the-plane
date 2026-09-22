"""The atom of de Grey's lemma: the smallest anti-rainbow gadget.

The lemma is that the six points at distance 2 from the origin show at most 2
of 4 colours.  Those six points are pairwise non-adjacent (their distances are
2, 2sqrt3 and 4), so nothing local forbids them from being rainbow.  The whole
force comes from the other 390 points of Sa.

"Shows at most 2 colours" is the same as "no three of them are pairwise
differently coloured".  So the lemma is 20 separate UNSAT statements, one per
triple, and each has a core.  The smallest of those cores is the real object:
a unit-distance graph H together with three pairwise non-adjacent vertices
that no proper 4-colouring can make rainbow.

Call that an anti-rainbow gadget of type (k=4, t=3).  Knowing its size tells
us what to look for at (k=5, t=4) -- which is exactly what a cap at five
colours would be, and therefore the first rung of the ladder to chi >= 6.
"""
import sys, time, itertools, pickle
from fractions import Fraction as Fr
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
from hn.degrey import build_Sa
from hn.geometry import Point
from hn.field import Field
from hn.graph import build_graph
from pysat.solvers import Solver

F = Field((3, 5, 7, 11))
ORI = Point(F.zero(), F.zero())
K = 4
pts = build_Sa(F)
g = build_graph(pts)
n = g.n
R = [i for i, p in enumerate(pts) if (p - ORI).norm2() == 4]
print(f"Sa: n={n}, m={sum(len(a) for a in g.adj)//2};  circle r^2=4: {R}", flush=True)
for i, j in itertools.combinations(R, 2):
    assert j not in g.adj[i], "circle points must be non-adjacent"
print("  the six are pairwise non-adjacent, as expected", flush=True)

X = lambda v, c: 1 + v * K + c
A = lambda v: 1 + n * K + v
base = [[-A(v)] + [X(v, c) for c in range(K)] for v in range(n)]
for u, v in g.edges():
    for c in range(K):
        base.append([-X(u, c), -X(v, c)])

def unsat_core(triple, subset):
    """Is 'these three are rainbow' unsatisfiable inside `subset`?  If so
    return a sub-subset that already refuses."""
    s = Solver(name="cd15", bootstrap_with=base)
    for idx, v in enumerate(triple):
        s.add_clause([X(v, idx)])
    assum = [A(v) for v in sorted(subset)]
    ok = s.solve(assumptions=assum)
    if ok:
        s.delete(); return None
    core = s.get_core() or assum
    s.delete()
    return set(l - 1 - n * K for l in core)

allv = set(range(n))
best = None
t0 = time.time()
for triple in itertools.combinations(R, 3):
    c = unsat_core(list(triple), allv)
    if c is None:
        print(f"  triple {triple}: RAINBOW IS POSSIBLE -- lemma false?", flush=True)
        continue
    # iterate the core
    for _ in range(40):
        c2 = unsat_core(list(triple), c | set(triple))
        if c2 is None or len(c2) >= len(c):
            break
        c = c2
    if best is None or len(c) < len(best[1]):
        best = (triple, c)
    print(f"  triple {triple}: core {len(c)}  [{time.time()-t0:.0f}s]", flush=True)

triple, core = best
core = core | set(triple)
print(f"\nsmallest core: triple {triple}, {len(core)} vertices.  Greedy shrink...", flush=True)
changed = True
while changed:
    changed = False
    for v in sorted(core, key=lambda v: len(g.adj[v])):
        if v in triple or v not in core:
            continue
        trial = core - {v}
        c2 = unsat_core(list(triple), trial)
        if c2 is not None:
            core = (c2 | set(triple)) if len(c2 | set(triple)) < len(trial) else trial
            changed = True
print(f"\nANTI-RAINBOW GADGET (k=4, t=3): {len(core)} vertices", flush=True)
sub = g.induced(sorted(core))
print(f"  induced: n={sub.n}, m={sum(len(a) for a in sub.adj)//2}", flush=True)
pickle.dump({"triple": list(triple), "core": sorted(core),
             "pts": [pts[i].approx() for i in sorted(core)]},
            open("/tmp/claude-0/-home-user-darwin-50/aceaa9ec-f432-5848-a506-39c59179b415/scratchpad/antirainbow4.pkl", "wb"))
