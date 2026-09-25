"""Glue symmetrically, so the union is invariant by construction.

Sa carries a group of order 12; the graphs built from it carry a group of
order 1, because the glue rotation is about ONE vertex and that breaks every
symmetry Sa had.  Taking the D6 orbit of the finished graph afterwards
imposes Sa's group on an object that is not shaped for it, and it loosens:
free@5 went from 11.46% to 13.14%.

The alternative is not to repair the symmetry but never to break it.  Glue at
every vertex of an orbit AT ONCE:

    U  =  Sa  u  U_{w in orbit(v)} rot_w(Sa)

For a rotation g of the group, g rot60_w g^-1 = rot60_{g(w)} and g(Sa) = Sa,
so g permutes the pieces and U is invariant.  The union is then a carrier with
a real group of its own -- one that came from the construction rather than
from Sa -- and any forced pair it has arrives in whole orbits.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time, random, json
from fractions import Fraction as Fr
from collections import defaultdict, Counter
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
from hn.degrey import build_Sa
from hn.geometry import Point, rotation_joining, _rot60
from hn.field import Field
from hn.graph import build_graph
from pysat.solvers import Solver

K = 4
NCOL = 16
F = Field((3, 11, 247))
Sa = build_Sa(F)
S = set(Sa)
rot60 = _rot60(F)
GLUE = {"rot60": rotation_joining(Fr(1), F),
        "rot120": rotation_joining(Fr(1, 3), F)}
WHICH = sys.argv[1] if len(sys.argv) > 1 else "rot60"
VERTS = [int(x) for x in sys.argv[2].split(",")] if len(sys.argv) > 2 else [25]

def orbit(p, reflect=False):
    out, q = [], p
    for _ in range(6):
        if q not in out: out.append(q)
        q = rot60(q)
    if reflect:
        for q in list(out):
            r = Point(q.x, -q.y)
            if r not in out: out.append(r)
    return out

def isometry_order(P):
    n = len(P)
    sx, sy = P[0].x, P[0].y
    for p in P[1:]: sx, sy = sx + p.x, sy + p.y
    inv = F.rational(Fr(1, n))
    C = Point(sx * inv, sy * inv)
    rel = [p - C for p in P]
    a = max(rel, key=lambda r: float(r.x * r.x + r.y * r.y))
    na = a.x * a.x + a.y * a.y
    sc = {Point(r.x * na, r.y * na) for r in rel}
    tot = 0
    for b in rel:
        if b.x * b.x + b.y * b.y != na: continue
        for sign in (1, -1):
            if sign > 0:
                co = a.x * b.x + a.y * b.y; si = a.x * b.y - a.y * b.x
            else:
                co = a.x * b.x - a.y * b.y; si = a.x * b.y + a.y * b.x
            ok = True
            for r in rel:
                q = (Point(co * r.x - si * r.y, si * r.x + co * r.y) if sign > 0
                     else Point(co * r.x + si * r.y, si * r.x - co * r.y))
                if q not in sc: ok = False; break
            if ok: tot += 1
    return tot

def analyse(P, tag):
    g = build_graph(P); n = g.n; m = sum(len(a) for a in g.adj) // 2
    X = lambda v, c: 1 + v * K + c
    cnf = [[X(v, c) for c in range(K)] for v in range(n)]
    for u, v in g.edges():
        for c in range(K):
            cnf.append([-X(u, c), -X(v, c)])
    s = Solver(name="m22", bootstrap_with=cnf)
    if not s.solve():
        s.delete()
        print(f"  {tag}: n={n} m={m} *** NOT {K}-COLOURABLE ***", flush=True)
        return
    pos = set(l for l in s.get_model() if l > 0)
    cols = [[next(c for c in range(K) if X(v, c) in pos) for v in range(n)]]
    rng = random.Random(7); tries = 0
    while len(cols) < NCOL and tries < NCOL * 4:
        tries += 1
        s.set_phases([(1 if rng.random() < 1.0 / K else -1) * X(v, c)
                      for v in range(n) for c in range(K)])
        if not s.solve(): continue
        p2 = set(l for l in s.get_model() if l > 0)
        cols.append([next(c for c in range(K) if X(v, c) in p2) for v in range(n)])
    s.delete()
    free = min(sum(1 for v in range(n)
                   if len({col[u] for u in g.adj[v]} | {col[v]}) < K) for col in cols)
    buck = defaultdict(list)
    for v in range(n):
        buck[tuple(c[v] for c in cols)].append(v)
    cand = [(a, b) for vs in buck.values() if len(vs) > 1
            for i, a in enumerate(vs) for b in vs[i+1:]]
    conf = []
    for a, b in cand[:4000]:
        s2 = Solver(name="m22", bootstrap_with=cnf)
        d = s2.solve(assumptions=[X(a, 0), X(b, 1)]); s2.delete()
        if not d:
            dd = (P[a] - P[b]).norm2()
            conf.append(str(dd) if dd.is_rational() else "irrational")
    print(f"  {tag}: n={n} m={m} deg={2.0*m/n:.2f} free@4={100.0*free/n:.2f}% "
          f"cand={len(cand)} FORCED={len(conf)} at {dict(Counter(conf).most_common(3))}"
          f"  group order {isometry_order(g.vertices)}", flush=True)

t0 = time.time()
print(f"Sa: 397 points, group order {isometry_order(Sa)}", flush=True)
for vi in VERTS:
    for reflect in (False, True):
        O = orbit(Sa[vi], reflect)
        seen, U = set(Sa), list(Sa)
        for w in O:
            rot = GLUE[WHICH].about(w)
            for p in Sa:
                q = rot(p)
                if q not in seen: seen.add(q); U.append(q)
        tag = (f"v{vi} {WHICH} over its {'D6' if reflect else 'C6'} orbit "
               f"({len(O)} centres)")
        analyse(U, tag)
        print(f"      [{time.time()-t0:.0f}s]", flush=True)
