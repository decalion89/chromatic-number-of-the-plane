"""Are the dense symmetric carriers still 4-chromatic, and how rigid?

Adding orbits of glue centres takes the mean degree from 9.94 to 18.32 while
keeping the full group and, so far, perfect rigidity at four colours.  Each
extra orbit is a fresh set of constraints on the same colouring, so at some
degree the carrier must stop being 4-colourable -- and where it stops is
exactly the question, because a carrier that refuses four at high degree is a
much better thing to spindle than Sa was.
"""
import sys, time, random
from fractions import Fraction as Fr
from collections import defaultdict
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
from hn.degrey import build_Sa
from hn.geometry import Point, rotation_joining, _rot60
from hn.field import Field
from hn.graph import build_graph
from pysat.solvers import Solver

F = Field((3, 11, 247)); rot60 = _rot60(F); Sa = build_Sa(F)
g60 = rotation_joining(Fr(1), F)
def orbit(p):
    out, q = [], p
    for _ in range(6):
        if q not in out: out.append(q)
        q = rot60(q)
    for q in list(out):
        r = Point(q.x, -q.y)
        if r not in out: out.append(r)
    return out
BEST = [265, 253, 211, 43, 139, 67, 241, 79, 277, 307]
LO = int(sys.argv[1]); HI = int(sys.argv[2])
t0 = time.time()
for k in range(LO, HI + 1):
    cs = []
    for b in BEST[:k]:
        for w in orbit(Sa[b]):
            if w not in cs: cs.append(w)
    seen, U = set(Sa), list(Sa)
    for w in cs:
        rot = g60.about(w)
        for p in Sa:
            q = rot(p)
            if q not in seen: seen.add(q); U.append(q)
    g = build_graph(U); n = g.n; m = sum(len(a) for a in g.adj) // 2
    tri = g.find_clique(3)
    verdict = None
    for K in (4, 5, 6):
        X = lambda v, c: 1 + v * K + c
        cnf = [[X(v, c) for c in range(K)] for v in range(n)]
        for u, v in g.edges():
            for c in range(K):
                cnf.append([-X(u, c), -X(v, c)])
        for i, v in enumerate(tri):
            cnf.append([X(v, i)])
            for c in range(K):
                if c != i: cnf.append([-X(v, c)])
        t1 = time.time()
        s = Solver(name="cd19", bootstrap_with=cnf)
        ok = s.solve()
        model = s.get_model() if ok else None
        s.delete()
        print(f"  {k} orbits: n={n} m={m} deg={2.0*m/n:.2f}  "
              f"{K}-colourable={ok}   [{time.time()-t1:.0f}s]", flush=True)
        if ok:
            verdict = K
            pos = set(l for l in model if l > 0)
            col = [next(c for c in range(K) if X(v, c) in pos) for v in range(n)]
            free = sum(1 for v in range(n)
                       if len({col[u] for u in g.adj[v]} | {col[v]}) < K)
            print(f"      chi = {K}, free@{K} = {100.0*free/n:.2f}% "
                  f"(one colouring)   [{time.time()-t0:.0f}s]", flush=True)
            break
        print(f"      *** it refuses {K} colours ***", flush=True)
