"""How many colour patterns the interface can carry -- the four/five gap, as a number.

The glue forces because agreeing on the shared set pins the rest of the
colouring.  So the quantity that decides it is not slack and not density: it is
how many DISTINCT patterns a small set of vertices can carry over all proper
k-colourings.  A carrier whose interface can only be coloured a handful of ways
leaves the second copy almost no freedom, and the forcing follows; one whose
interface can be coloured a thousand ways leaves the sigma-argument intact.

Measured as a saturation curve: generate colourings one at a time, each blocked
away from the last on a random sample so it is genuinely different, and count
the distinct patterns on a fixed set T up to permutation of the colours.  A
carrier that saturates has a small interface space; one whose count keeps
climbing with every colouring does not.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, json, time, random
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
from fractions import Fraction as Fr
from hn.degrey import build_Sa, build_G
from hn.geometry import Point, rotation_joining
from hn.field import Field
from hn.graph import build_graph
from pysat.solvers import Solver

WHICH, K, T, M = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4])

def load(which):
    if which in ("Sa", "H"):
        F = Field((3, 11, 247))
        pts = build_Sa(F)
        if which == "H":
            rot = rotation_joining(Fr(1), F).about(pts[25])
            seen, out = set(), []
            for p in pts:
                for q in (p, rot(p)):
                    if q not in seen: seen.add(q); out.append(q)
            return out
        return pts
    if which == "G":
        return build_G(Field((3, 5, 7, 11)), as_graph=False)
    d = json.load(open(f"{HN_DIR}/data/{which}"))
    F = Field(tuple(d["field_generators"]))
    return [Point(F.element([Fr(a, b) for a, b in x]),
                  F.element([Fr(a, b) for a, b in y])) for x, y in d["points"]]

pts = load(WHICH)
g = build_graph(pts); n = g.n
X = lambda v, c: 1 + v * K + c
cnf = [[X(v, c) for c in range(K)] for v in range(n)]
for u, v in g.edges():
    for c in range(K):
        cnf.append([-X(u, c), -X(v, c)])
# At-most-one is normally omitted -- the edge clauses already make adjacent
# colour SETS disjoint, so picking any colour from each set is proper.  But
# here the colour is READ OFF the model, and a vertex carrying three true
# colour variables reads as whichever is smallest, which is a degenerate
# function of the model rather than of the colouring.  That is what made an
# eight-vertex interface report a single pattern over four hundred colourings.
for v in range(n):
    for c in range(K):
        for dd in range(c + 1, K):
            cnf.append([-X(v, c), -X(v, dd)])
rng = random.Random(31)
# T vertices of high degree, which is where the glue actually bites
Tset = sorted(range(n), key=lambda v: -len(g.adj[v]))[:T]
print(f"{WHICH}: n={n} k={K}, interface of {T} high-degree vertices", flush=True)

def canon(pat):
    seen, out = {}, []
    for c in pat:
        if c not in seen: seen[c] = len(seen)
        out.append(seen[c])
    return tuple(out)

s = Solver(name="m22", bootstrap_with=cnf)
if not s.solve(): sys.exit("not colourable")
pos = set(l for l in s.get_model() if l > 0)
col = [next(c for c in range(K) if X(v, c) in pos) for v in range(n)]
seenpat = {canon([col[v] for v in Tset])}
t0 = time.time()
marks = {25, 50, 100, 200, 400, 800}
for i in range(1, M + 1):
    # Blocking alone is too weak for a readout: the solver satisfies the clause
    # by changing one vertex of the forty and leaves the rest where they were,
    # so a fixed interface can report one pattern over four hundred colourings
    # purely out of laziness.  Randomised polarity moves the whole assignment;
    # blocking guarantees it moves at all.
    samp = rng.sample(range(n), 40)
    s.add_clause([-X(v, col[v]) for v in samp])
    s.set_phases([(1 if rng.random() < 1.0 / K else -1) * X(v, c)
                  for v in range(n) for c in range(K)])
    if not s.solve(): break
    pos = set(l for l in s.get_model() if l > 0)
    col = [next(c for c in range(K) if X(v, c) in pos) for v in range(n)]
    seenpat.add(canon([col[v] for v in Tset]))
    if i in marks:
        print(f"  after {i:>4} colourings: {len(seenpat):>5} distinct patterns"
              f"   [{time.time()-t0:.0f}s]", flush=True)
s.delete()
print(f"  final: {len(seenpat)} distinct patterns from {i} colourings "
      f"({100.0*len(seenpat)/max(i,1):.0f}% new)   [{time.time()-t0:.0f}s]", flush=True)
