"""Tightness at k: the fraction of vertices whose closed neighbourhood does
NOT already use all k colours.  Sa sits at 0.0% at four -- every vertex sees
every colour -- and that is the carrier whose glue manufactured forcing.  If
tightness is what the glue needs, the number to watch at five is this one."""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, random
from fractions import Fraction as Fr
from collections import Counter
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
from hn.degrey import build_Sa, build_Y, build_G
from hn.geometry import rotation_joining
from hn.field import Field
from hn.graph import build_graph
from pysat.solvers import Solver

def Zpts():
    F = Field((3, 11, 247))
    base = build_Sa(F)
    r1 = rotation_joining(Fr(1), F).about(base[25])
    seen, H = set(), []
    for p in base:
        for q in (p, r1(p)):
            if q not in seen: seen.add(q); H.append(q)
    r2 = rotation_joining(Fr(64, 9), F).about(H[157])
    seen, Z = set(), []
    for p in H:
        for q in (p, r2(p)):
            if q not in seen: seen.add(q); Z.append(q)
    return Z, H

F4 = Field((3, 5, 7, 11))
Z, H = Zpts()
Zg = build_graph(Z)
ZH = build_graph(H)
F2 = Field((3, 11, 247))
Zglue = None

carriers = [
    ("Sa", build_graph(build_Sa(F4)), 4),
    ("Y",  build_graph(build_Y(F4)),  4),
    ("H (Sa glued at a vertex)", ZH, 4),
    ("G",  build_G(F4),               5),
    ("Z (the new one)", Zg,           5),
]
# and Z glued once more, to see whether iterating tightens anything
r = rotation_joining(Fr(1), F2).about(Z[0])
seen, Z2 = set(), []
for p in Z:
    for q in (p, r(p)):
        if q not in seen: seen.add(q); Z2.append(q)
carriers.append(("Z glued at a vertex", build_graph(Z2), 5))

print(f"{'carrier':<28} {'n':>6} {'m':>7} {'deg':>6} {'k':>3} {'free':>8}")
for name, g, K in carriers:
    n = g.n
    m = sum(len(a) for a in g.adj) // 2
    X = lambda v, c: 1 + v * K + c
    cnf = [[X(v, c) for c in range(K)] for v in range(n)]
    for u, v in g.edges():
        for c in range(K):
            cnf.append([-X(u, c), -X(v, c)])
    s = Solver(name="m22", bootstrap_with=cnf)
    best = None
    rng = random.Random(5)
    for _ in range(12):
        s.set_phases([(1 if rng.random() < 1.0 / K else -1) * X(v, c)
                      for v in range(n) for c in range(K)])
        if not s.solve(): continue
        pos = set(l for l in s.get_model() if l > 0)
        col = [next(c for c in range(K) if X(v, c) in pos) for v in range(n)]
        free = sum(1 for v in range(n)
                   if len({col[u] for u in g.adj[v]} | {col[v]}) < K)
        if best is None or free < best: best = free
    s.delete()
    print(f"{name:<28} {n:>6} {m:>7} {2.0*m/n:>6.2f} {K:>3} "
          f"{100.0*best/n:>7.2f}%", flush=True)
