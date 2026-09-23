"""How many RINGS about one vertex does it take to pin its colour?

A scattered disjunction cannot be consumed: copies of a graph compose only when
their conclusions name the same point.  The consumable shape is

    every 5-colouring gives  c(h) = c(v)  for some v at one of r distances
    d_1..d_r from h

because rotations about h fix h, so every rotated copy says this about the SAME
h, and all the partners it can name live on r circles about h.  Two points on
those circles are a unit apart at a computable angle, which is exactly the
multispindle: with enough copies at well-chosen angles, every choice of partner
has to clash with another copy's.

So: greedily forbid h from sharing with the richest rings about it until the
graph dies, then minimise.  The number r of rings is the width of the
disjunction and the price of consuming it.

Cheap in the right direction.  Forbidding EVERY pair (h,v) asks whether G - h
is 4-colourable, a refutation that costs a minute on a thousand vertices; a
partial ring set is usually satisfiable, and satisfiable is the fast direction,
so the greedy walk is quick until the step that actually wins.

Pigeonhole-free without a check: a unit-distance graph in the plane has clique
number 3, and every extra edge here is incident to h, so no K6 can appear.
"""
import sys, time, json, random
from fractions import Fraction as Fr
from collections import defaultdict
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
from hn.field import Field
from hn.geometry import Point
from hn.graph import build_graph
from pysat.solvers import Solver

ROOT = "/home/user/darwin-50/research/hadwiger-nelson"
t0 = time.time(); K = 5
NAME = sys.argv[1] if len(sys.argv) > 1 else "five_247.json"
HUBS = int(sys.argv[2]) if len(sys.argv) > 2 else 12
MAXR = int(sys.argv[3]) if len(sys.argv) > 3 else 8
d = json.load(open(f"{ROOT}/data/{NAME}"))
F = Field(tuple(d["field_generators"]))
P = [Point(F.element([Fr(a, b) for a, b in x]),
           F.element([Fr(a, b) for a, b in y])) for x, y in d["points"]]
g = build_graph(P); n = g.n
E = list(g.edges())
adj = defaultdict(set)
for x, y in E:
    adj[x].add(y); adj[y].add(x)
hx = [q.fx for q in g.vertices]; hy = [q.fy for q in g.vertices]
print(f"{NAME} n={n} unit edges={len(E)}   [{time.time()-t0:.0f}s]", flush=True)

X = lambda v, c: 1 + v * K + c
base = [[X(v, c) for c in range(K)] for v in range(n)]
for x, y in E:
    for c in range(K):
        base.append([-X(x, c), -X(y, c)])

def dead(h, vs):
    cnf = list(base)
    for v in vs:
        for c in range(K):
            cnf.append([-X(h, c), -X(v, c)])
    s = Solver(name="cd19", bootstrap_with=cnf)
    r = s.solve(); s.delete()
    return not r

random.seed(11)
order = sorted(range(n), key=lambda v: -len(adj[v]))
order = order[:HUBS // 2] + random.sample(range(n), HUBS)
seen = set(); results = []
for h in order:
    if h in seen: continue
    seen.add(h)
    rings = defaultdict(list)
    for v in range(n):
        if v == h or v in adj[h]: continue
        rings[round((hx[h]-hx[v])**2 + (hy[h]-hy[v])**2, 7)].append(v)
    ranked = sorted(rings.items(), key=lambda kv: -len(kv[1]))[:MAXR]
    # Forbidding more rings only adds clauses, so "dead" is monotone in the
    # prefix length and a binary search finds the smallest one in about nine
    # solves instead of four hundred.  Ask the whole list first: if even that
    # is alive then h's colour class can be {h} alone, G - h is 4-colourable,
    # and no prefix can win.
    # The filter has to forbid EVERY non-neighbour, not just the ranked prefix:
    # leaving any partner free lets h share with it, so a SAT answer there says
    # nothing.  With all of them forbidden, SAT means exactly that h's colour
    # class can be {h} alone, i.e. G - h is 4-colourable, and then no subset
    # can win either.
    allv = [v for v in range(n) if v != h and v not in adj[h]]
    if not dead(h, allv):
        print(f"  h={h} deg={len(adj[h])}: alive with ALL {len(allv)} partners "
              f"forbidden -- G-h is 4-colourable, no hub here   "
              f"[{time.time()-t0:.0f}s]", flush=True)
        continue
    print(f"  h={h} deg={len(adj[h])}: HUB (G-h is 5-chromatic); "
          f"binary searching rings   [{time.time()-t0:.0f}s]", flush=True)
    lo, hi = 1, len(ranked)
    while lo < hi:
        mid = (lo + hi) // 2
        if dead(h, [v for _, part in ranked[:mid] for v in part]):
            hi = mid
        else:
            lo = mid + 1
    chosen = [d2 for d2, _ in ranked[:hi]]
    print(f"  h={h} deg={len(adj[h])}: smallest prefix {hi} rings   "
          f"[{time.time()-t0:.0f}s]", flush=True)
    # minimise the ring set
    minimal = list(chosen)
    for d2 in list(chosen):
        trial = [x for x in minimal if x != d2]
        if trial and dead(h, [v for x in trial for v in rings[x]]):
            minimal = trial
    tot = sum(len(rings[x]) for x in minimal)
    print(f"  *** h={h} deg={len(adj[h])}: {len(minimal)} rings suffice "
          f"({tot} partners) d^2 = {[round(x,6) for x in minimal]}   "
          f"[{time.time()-t0:.0f}s]", flush=True)
    results.append({"hub": h, "rings": minimal, "partners": tot,
                    "by_ring": {str(x): len(rings[x]) for x in minimal}})
    json.dump({"graph": NAME, "results": results},
              open(f"{ROOT}/data/hubrings_{NAME}", "w"))
print(f"\n  {len(results)} hubs pinned   [{time.time()-t0:.0f}s]", flush=True)
