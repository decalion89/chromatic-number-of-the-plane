"""chi of N(h) united with the whole 1/sqrt3 ring -- can the local floor pass 3?

Two ceilings have been proved here and both are exact.  On the unit circle the
induced graph is a union of paths and hexagons, so chi = 2.  On the 1/sqrt3
ring it is a disjoint union of triangles, so chi = 3.  Niven says there is no
third radius with an odd cycle.

But the two rings TOGETHER are a different graph.  A point u at radius 1 and a
point t at radius 1/sqrt3 are a unit apart when their angles differ by the
angle with cos = sqrt3/6, about 73.22 degrees, so cross edges exist -- and the
union is no longer a disjoint anything.

With ONE triangle the union is still 3-colourable, always, and the argument is
list colouring: two triangle vertices are 120 degrees apart while u sits at
+-73.22 from any neighbour of its own, so u is adjacent to AT MOST ONE vertex
of the triangle and keeps a list of size >= 2 out of the triangle's three
colours; N(h) is a union of paths and even cycles, and those are 2-choosable.
So mu(N(h) u T) <= 3 is a theorem, which is exactly what the scan measured --
3 at all 436 pairs of the 803-graph and all 450 of the 1139.

With SEVERAL triangles the lists can empty.  u may be adjacent to one vertex of
each of up to six triangles, and if those carry three distinct colours its list
is gone.  So chi of the union is not bounded by 3 any more, and the question is
what it actually is.  That is a purely local computation -- fifty vertices, no
colouring of the ambient graph needed -- and it is the first place in this
project where the local floor has any room at all.
"""
import sys, time, json
from fractions import Fraction as Fr
from collections import defaultdict, Counter
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
from hn.field import Field
from hn.geometry import Point
from hn.graph import build_graph
from pysat.solvers import Solver

ROOT = "/home/user/darwin-50/research/hadwiger-nelson"
t0 = time.time()
NAMES = sys.argv[1:] or ["five_247_c.json"]

def chi_of(verts, edges, cap=6):
    for k in range(1, cap + 1):
        idx = {v: i for i, v in enumerate(verts)}
        X = lambda v, c: 1 + idx[v] * k + c
        cnf = [[X(v, c) for c in range(k)] for v in verts]
        for a, b in edges:
            for c in range(k):
                cnf.append([-X(a, c), -X(b, c)])
        s = Solver(name="cd19", bootstrap_with=cnf)
        ok = s.solve(); s.delete()
        if ok: return k
    return cap + 1

for NAME in NAMES:
    d = json.load(open(f"{ROOT}/data/{NAME}"))
    F = Field(tuple(d["field_generators"]))
    P = [Point(F.element([Fr(a, b) for a, b in x]),
               F.element([Fr(a, b) for a, b in y])) for x, y in d["points"]]
    g = build_graph(P); n = g.n
    adj = defaultdict(set)
    for x, y in g.edges():
        adj[x].add(y); adj[y].add(x)
    hx = [q.fx for q in g.vertices]; hy = [q.fy for q in g.vertices]
    third = F.rational(Fr(1, 3))
    ring = defaultdict(list)
    for i in range(n):
        for j in range(i + 1, n):
            dd = (hx[i]-hx[j])**2 + (hy[i]-hy[j])**2
            if abs(dd - 1/3) < 1e-9 and (g.vertices[i]-g.vertices[j]).norm2() == third:
                ring[i].append(j); ring[j].append(i)
    tally = Counter(); best = None
    order = sorted(ring, key=lambda h: -(len(adj[h]) + len(ring[h])))
    for h in order[:400]:
        V = sorted(set(adj[h]) | set(ring[h]))
        if len(V) < 4: continue
        S = set(V)
        Eloc = [(a, b) for a in V for b in adj[a] if b in S and b > a]
        k = chi_of(V, Eloc)
        tally[k] += 1
        if best is None or k > best[0]:
            best = (k, h, len(adj[h]), len(ring[h]), len(V), len(Eloc))
            print(f"    chi = {k} at h={h}: |N|={len(adj[h])} |ring|={len(ring[h])} "
                  f"|X|={len(V)} edges={len(Eloc)}   [{time.time()-t0:.0f}s]",
                  flush=True)
    print(f"\n{NAME} n={n}: chi(N(h) u ring) over {sum(tally.values())} hubs: "
          f"{dict(sorted(tally.items()))}   [{time.time()-t0:.0f}s]", flush=True)
