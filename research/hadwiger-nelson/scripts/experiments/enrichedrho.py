"""Use monotonicity: any 5-chromatic graph with rho <= 63 will do, not just G.

rho(W) <= rho(G) whenever G sits inside W, and a core of three needs only SOME
5-chromatic graph with rho <= deg(p) + 3.  So the chord-enriched graph is a
strictly better place to look than G itself: it contains G, so it is
5-chromatic and its rho is no larger, and it carries 1680 auxiliaries around
one pivot where G carries 150.

That local structure is exactly what ambient forcing is made of, and the
pivot's degree there can be pushed to 134, which relaxes the target from 63 to
137.

Same decision loop, driven by a greedy cover so the escaping colourings stay
diverse.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time, itertools
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
from pysat.solvers import Solver
from hn import degrey
from hn.geometry import Point
from hn.graph import build_graph


def greedy_cover(fam, n):
    if not fam:
        return []
    sets = [set(cl) for cl in fam]
    live = set(range(len(fam)))
    count = [0] * n
    for i in live:
        for v in sets[i]:
            count[v] += 1
    chosen = []
    while live:
        v = max(range(n), key=lambda u: count[u])
        chosen.append(v)
        for i in [i for i in live if v in sets[i]]:
            live.discard(i)
            for u in sets[i]:
                count[u] -= 1
    return chosen


g = degrey.build_G()
pts = list(g.vertices)
deg = g.degrees()
p_idx = max(range(g.n), key=lambda v: deg[v])
p = pts[p_idx]
circle = sorted(g.adj[p_idx])
have, extra = set(pts), []
for i, j in itertools.combinations(circle, 2):
    q = Point(pts[i].x + pts[j].x - p.x, pts[i].y + pts[j].y - p.y)
    if q != p and q not in have:
        have.add(q)
        extra.append(q)
w = build_graph(pts + extra)
pi = w.index_of(p)
BUD = len(w.adj[pi]) + 3
n, k = w.n, 5
print(f"chord-enriched G: {n} vertices, {w.m} edges; pivot degree "
      f"{len(w.adj[pi])}, so a core of three needs rho <= {BUD}", flush=True)


def x(v, c):
    return 1 + v * k + c


cls = [[x(v, c) for c in range(k)] for v in range(n)]
for a, b in w.edges():
    for c in range(k):
        cls.append([-x(a, c), -x(b, c)])
colour = Solver(name="cd19", bootstrap_with=cls)
fam, t0, best = [], time.time(), 0
try:
    for rnd in range(200000):
        S = greedy_cover(fam, n)
        if len(S) > BUD:
            print(f"  round {rnd}: greedy needs {len(S)} > {BUD}  "
                  f"[{time.time()-t0:.0f}s]", flush=True)
            break
        if not colour.solve(assumptions=[-x(v, 0) for v in S]):
            print(f"  round {rnd}: S of {len(S)} is FORCING -- "
                  f"rho <= {len(S)}   *** AND {len(S)} <= {BUD} ***  "
                  f"[{time.time()-t0:.0f}s]", flush=True)
            break
        mm = set(colour.get_model())
        fam.append([v for v in range(n) if x(v, 0) in mm])
        if len(S) > best:
            best = len(S)
            print(f"  round {rnd}: |F| = {len(fam)}, greedy cover {len(S)}  "
                  f"[{time.time()-t0:.0f}s]", flush=True)
finally:
    colour.delete()
print(f"largest greedy cover: {best} against the target {BUD}")
