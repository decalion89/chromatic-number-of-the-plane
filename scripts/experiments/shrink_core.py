"""Extract a small 5-chromatic subgraph by UNSAT core, not by deletion.

Everything downstream is priced in copies of a 5-chromatic graph: a forced
pair's certificate is the whole near-critical graph, the 1/sqrt3 closure needs
three copies of it, the square needs two.  So its ORDER is the one quantity that
makes every other attack cheaper, and this project's best is 803 against a
published 509.

Greedy deletion is how the 803 was found, and it costs one refutation per
vertex -- the slow direction, minutes each on a thousand vertices.  A core costs
one.  Give every vertex a selector s_v and make "v has a colour" conditional on
it; an inactive vertex simply takes no colour and its edges go vacuous, so
satisfiability is unchanged, but an UNSAT answer under all the selectors comes
back with the subset that caused it -- a 5-chromatic subgraph, in one solve.
Re-assuming just the core shrinks it again, and a greedy pass over what survives
makes it vertex-critical.

Colour symmetry is 120 refutations of every refutation, so a triangle is pinned
to colours 0, 1, 2 and forced active.  That is sound for deciding colourability
and took the 803-graph's own refutation from over five minutes to 47 seconds.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time, json, random
from fractions import Fraction as Fr
from collections import defaultdict
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
from hn.field import Field
from hn.geometry import Point
from hn.graph import build_graph
from pysat.solvers import Solver

ROOT = HN_DIR
# FOUR, not five.  A 5-chromatic graph IS 5-colourable -- that is what chi = 5
# means -- so the property to certify is that it refuses FOUR, and the core to
# extract is the one that makes the 4-colouring instance unsatisfiable.  A first
# run asked for five, got SAT, and reported nothing to extract.
t0 = time.time(); K = 4
NAME = sys.argv[1] if len(sys.argv) > 1 else "five_247.json"
SEED = int(sys.argv[2]) if len(sys.argv) > 2 else 0
d = json.load(open(f"{ROOT}/data/{NAME}"))
F = Field(tuple(d["field_generators"]))
P = [Point(F.element([Fr(a, b) for a, b in x]),
           F.element([Fr(a, b) for a, b in y])) for x, y in d["points"]]
g = build_graph(P); n = g.n
E = list(g.edges())
adj = defaultdict(set)
for x, y in E:
    adj[x].add(y); adj[y].add(x)
tri = None
for u in range(n):
    for v in sorted(adj[u]):
        w = (adj[u] & adj[v])
        if w:
            tri = (u, v, min(w)); break
    if tri: break
print(f"{NAME} n={n} edges={len(E)}; pinned triangle {tri}   "
      f"[{time.time()-t0:.0f}s]", flush=True)
X = lambda v, c: 1 + v * K + c
S = lambda v: n * K + 1 + v
cnf = [[-S(v)] + [X(v, c) for c in range(K)] for v in range(n)]
for x, y in E:
    for c in range(K):
        cnf.append([-X(x, c), -X(y, c)])
for i, v in enumerate(tri):
    cnf.append([X(v, i)])
s = Solver(name="cd19", bootstrap_with=cnf)
ass = [S(v) for v in range(n)]
ok = s.solve(assumptions=ass)
print(f"  whole graph 4-colourable: {ok}   [{time.time()-t0:.0f}s]", flush=True)
if ok:
    print("  4-colourable -- not 5-chromatic, nothing to extract", flush=True)
    sys.exit(0)
core = set(s.get_core() or [])
sizes = [len(core)]
print(f"  first core: {len(core)} vertices   [{time.time()-t0:.0f}s]", flush=True)
for _ in range(40):
    a2 = sorted(core)
    if s.solve(assumptions=a2): break
    c2 = set(s.get_core() or [])
    if len(c2) >= len(core): break
    core = c2; sizes.append(len(core))
    print(f"    core -> {len(core)}   [{time.time()-t0:.0f}s]", flush=True)
cur = sorted(core)
random.seed(SEED); random.shuffle(cur)
i = 0; passes = 0
while i < len(cur):
    trial = cur[:i] + cur[i+1:]
    keep = set(trial)
    if all(S(v) in keep for v in tri) and trial and not s.solve(assumptions=trial):
        c2 = set(s.get_core() or trial)
        cur = [a for a in trial if a in c2] or trial
        i = 0; passes += 1
        if passes % 25 == 0:
            print(f"    greedy: {len(cur)} left   [{time.time()-t0:.0f}s]", flush=True)
    else:
        i += 1
verts = sorted(a - (n * K + 1) for a in cur)
sub = [(x, y) for x, y in E if x in set(verts) and y in set(verts)]
print(f"\n  SUBGRAPH THAT REFUSES FOUR: {len(verts)} vertices, {len(sub)} edges   "
      f"[{time.time()-t0:.0f}s]", flush=True)
json.dump({"source": NAME, "seed": SEED, "n": len(verts),
           "field_generators": list(F.gens),
           "points": [[[[t.numerator, t.denominator] for t in g.vertices[v].x.c],
                       [[t.numerator, t.denominator] for t in g.vertices[v].y.c]]
                      for v in verts]},
          open(f"{ROOT}/data/shrunk_{len(verts)}.json", "w"))
print(f"  written data/shrunk_{len(verts)}.json", flush=True)
