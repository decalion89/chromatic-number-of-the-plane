"""Shrink Z to a small subgraph that still refuses four colours.

find_uncolorable_core turns off symmetry breaking, which multiplies the search
by 4! on every UNSAT proof and is why it made no progress on 1139 vertices in
an hour.  Pinning a triangle is sound here as long as the pinning is made
conditional on the vertex selectors and the triangle is kept in every subset
tried -- then a vertex that is switched off carries no constraint at all.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time, json
from fractions import Fraction as Fr
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
from hn.degrey import build_Sa
from hn.geometry import rotation_joining
from hn.field import Field
from hn.graph import build_graph
from pysat.solvers import Solver

K1 = Field((3, 11, 247))
pts = build_Sa(K1)
r1 = rotation_joining(Fr(1), K1).about(pts[25])
seen, H = set(), []
for p in pts:
    for q in (p, r1(p)):
        if q not in seen: seen.add(q); H.append(q)
r2 = rotation_joining(Fr(64, 9), K1).about(H[157])
seen, Z = set(), []
for p in H:
    for q in (p, r2(p)):
        if q not in seen: seen.add(q); Z.append(q)
g = build_graph(Z).k_core(4)
print(f"Z's 4-core: n={g.n}, m={sum(len(a) for a in g.adj)//2}", flush=True)
n, K = g.n, 4
X = lambda v, c: 1 + v * K + c
A = lambda v: 1 + n * K + v
cnf = [[-A(v)] + [X(v, c) for c in range(K)] for v in range(n)]
for u, v in g.edges():
    for c in range(K):
        cnf.append([-X(u, c), -X(v, c)])
tri = g.find_clique(3)
print(f"  pinning triangle {tri} conditionally on its selectors", flush=True)
for i, v in enumerate(tri):
    cnf.append([-A(v), X(v, i)])
    for c in range(K):
        if c != i:
            cnf.append([-A(v), -X(v, c)])

def refuses(sub):
    s = Solver(name="m22", bootstrap_with=cnf)
    ok = s.solve(assumptions=[A(v) for v in sorted(sub)])
    if ok:
        s.delete(); return None
    core = s.get_core()
    s.delete()
    return set(l - 1 - n * K for l in core) if core else set(sub)

t0 = time.time()
sub = set(range(n))
assert refuses(sub) is not None, "Z's 4-core is 4-colourable?!"
for r in range(60):
    c = refuses(sub | set(tri))
    if c is None or len(c) >= len(sub): break
    sub = c | set(tri)
    print(f"  core round {r+1}: {len(sub)}   [{time.time()-t0:.0f}s]", flush=True)
print(f"  after cores: {len(sub)}   [{time.time()-t0:.0f}s]", flush=True)
order = sorted(sub, key=lambda v: len(g.adj[v]))
for v in order:
    if v in tri or v not in sub: continue
    trial = sub - {v}
    c = refuses(trial)
    if c is not None:
        sub = (c | set(tri)) if len(c | set(tri)) < len(trial) else trial
print(f"\n4-critical-ish core of Z: {len(sub)} vertices   [{time.time()-t0:.0f}s]",
      flush=True)
core = build_graph([Z[v] for v in sorted(sub)]).k_core(4).largest_component()
print(f"  cleaned: n={core.n}, m={sum(len(a) for a in core.adj)//2}", flush=True)
s = Solver(name="m22")
gg = core
XX = lambda v, c: 1 + v * K + c
cc = [[XX(v, c) for c in range(K)] for v in range(gg.n)]
for u, v in gg.edges():
    for c in range(K):
        cc.append([-XX(u, c), -XX(v, c)])
s = Solver(name="m22", bootstrap_with=cc)
print(f"  still refuses four colours: {not s.solve()}", flush=True)
s.delete()
def dump(pt):
    return [[[c.numerator, c.denominator] for c in pt.x.c],
            [[c.numerator, c.denominator] for c in pt.y.c]]
json.dump({"field_generators": list(K1.gens), "n": core.n,
           "parent": "data/five_247.json",
           "points": [dump(p) for p in core.vertices]},
          open(HN_DIR + "/data/five_247_core.json", "w"))
print("  written data/five_247_core.json", flush=True)
