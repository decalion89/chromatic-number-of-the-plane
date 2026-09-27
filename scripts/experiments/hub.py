"""A hub disjunction: in every 5-colouring, h shares its colour with one of these.

delta <= 10 on the 803-graph says ten specific non-unit pairs, forbidden, kill
five colours -- a disjunction of width ten, which is de Grey's shape.  But those
ten pairs are scattered, and a scattered disjunction cannot be consumed: copies
compose only when the conclusions share a point.

So ask for the disjunction at a HUB.  Forbid every pair (h, v) at once: that
says h's colour class is {h} alone, which is possible exactly when G - h is
4-colourable.  On a vertex-critical graph it always is, so no hub disjunction
exists there -- and that is why criticality, which made the ceiling free, is
the wrong property here.  On a graph that is 5-chromatic but NOT critical,
G - h is still 5-chromatic for most h, the instance dies, and the UNSAT core is

    P(h) = a set of partners such that every 5-colouring has c(h) = c(p)
           for some p in P(h)

which is consumable: rotations about h fix h, so copies of G rotated about h
each contribute their own conclusion about the SAME point, and the conclusions
clash when two of the rotated partners land a unit apart.  That is exactly the
multispindle, and the width of P(h) is what it costs.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time, json
from fractions import Fraction as Fr
from collections import defaultdict, Counter
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
from hn.field import Field
from hn.geometry import Point
from hn.graph import build_graph
from pysat.solvers import Solver

ROOT = HN_DIR
t0 = time.time(); K = 5
NAME = sys.argv[1] if len(sys.argv) > 1 else "five_247.json"
HOWMANY = int(sys.argv[2]) if len(sys.argv) > 2 else 40
d = json.load(open(f"{ROOT}/data/{NAME}"))
F = Field(tuple(d["field_generators"]))
P = [Point(F.element([Fr(a, b) for a, b in x]),
           F.element([Fr(a, b) for a, b in y])) for x, y in d["points"]]
g = build_graph(P); n = g.n
E = list(g.edges())
adj = defaultdict(set)
for x, y in E:
    adj[x].add(y); adj[y].add(x)
print(f"{NAME} n={n} unit edges={len(E)}   [{time.time()-t0:.0f}s]", flush=True)

X = lambda v, c: 1 + v * K + c
S = lambda v: n * K + 1 + v          # selector for "h may not share with v"
cnf = [[X(v, c) for c in range(K)] for v in range(n)]
for x, y in E:
    for c in range(K):
        cnf.append([-X(x, c), -X(y, c)])

best = []
order = sorted(range(n), key=lambda v: -len(adj[v]))[:HOWMANY]
for h in order:
    part = [v for v in range(n) if v != h and v not in adj[h]]
    local = list(cnf)
    for v in part:
        for c in range(K):
            local.append([-S(v), -X(h, c), -X(v, c)])
    s = Solver(name="cd19", bootstrap_with=local)
    ass = [S(v) for v in part]
    if s.solve(assumptions=ass):
        s.delete(); continue                      # G - h is 4-colourable
    core = set(s.get_core() or [])
    for _ in range(25):
        a2 = sorted(core)
        if s.solve(assumptions=a2): break
        c2 = set(s.get_core() or [])
        if len(c2) >= len(core): break
        core = c2
    cur = sorted(core)
    i = 0
    while i < len(cur):
        trial = cur[:i] + cur[i+1:]
        if trial and not s.solve(assumptions=trial):
            c2 = set(s.get_core() or trial)
            cur = [a for a in trial if a in c2] or trial
            i = 0
        else:
            i += 1
    parts = [a - (n * K + 1) for a in cur]
    ds = Counter(round(float((g.vertices[h] - g.vertices[v]).norm2()), 6) for v in parts)
    print(f"  hub {h} (deg {len(adj[h])}): width {len(parts)}  distances^2 "
          f"{dict(sorted(ds.items()))}   [{time.time()-t0:.0f}s]", flush=True)
    best.append((len(parts), h, parts))
    s.delete()
best.sort()
if best:
    w, h, parts = best[0]
    print(f"\n  NARROWEST HUB: {h}, width {w}   [{time.time()-t0:.0f}s]", flush=True)
    json.dump({"graph": NAME, "hub": h, "partners": parts,
               "width": w}, open(f"{ROOT}/data/hub_{NAME}", "w"))
else:
    print(f"\n  no hub disjunction: every vertex can stand alone in its colour "
          f"(the graph is vertex-critical here)   [{time.time()-t0:.0f}s]", flush=True)
