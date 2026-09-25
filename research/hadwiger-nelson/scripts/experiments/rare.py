"""Hunt carriers with a RARE forced distance.

Forced orbits are not equal.  The 1021-point carrier has 144 pairs at
d^2 = 64/9 and 9 at 64/3; spent alone, each gives a 5-chromatic graph of the
same size, and cadical needs 322 seconds for the common one and 1060 for the
rare one.  Rarity is worth more than count, and density destroys it -- the
densest carrier has 102 forced orbits and not one of them is rare.

So scan for rarity, not for count.  The distance of a candidate pair is known
before any solver call, so the scan confirms only the representatives whose
distance is not the common one; the bulk is skipped and the rare ones cost a
handful of calls each.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time, random
from fractions import Fraction as Fr
from collections import defaultdict, Counter
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
from hn.degrey import build_Sa
from hn.geometry import Point, rotation_joining, _rot60
from hn.field import Field
from hn.graph import build_graph
from pysat.solvers import Solver

F = Field((3, 11, 247)); rot60 = _rot60(F); Sa = build_Sa(F)
g60 = rotation_joining(Fr(1), F)
COMMON = Fr(64, 9)
LO, HI = int(sys.argv[1]), int(sys.argv[2])
def orb(p):
    out, q = [], p
    for _ in range(6):
        if q not in out: out.append(q)
        q = rot60(q)
    for q in list(out):
        r = Point(q.x, -q.y)
        if r not in out: out.append(r)
    return out
def sqfree(r):
    d = 2
    while d * d <= r:
        while r % (d * d) == 0: r //= d * d
        d += 1
    return r
def radical(q):
    if q <= Fr(1, 4): return None
    c = Fr(1) - Fr(1, 2) / q
    s2 = 1 - c * c
    if s2 == 0: return 1
    return sqfree(s2.numerator * s2.denominator)
t0 = time.time()
seenorb = set()
for b in range(LO, min(HI, len(Sa))):
    O = orb(Sa[b])
    key = frozenset(O)
    if key in seenorb: continue
    seenorb.add(key)
    seen, U = set(Sa), list(Sa)
    for w in O:
        rot = g60.about(w)
        for p in Sa:
            q = rot(p)
            if q not in seen: seen.add(q); U.append(q)
    if len(U) == len(Sa): continue
    g = build_graph(U); n = g.n
    pos = {p: i for i, p in enumerate(g.vertices)}
    Sset = set(g.vertices)
    maps = []
    for refl in (False, True):
        base = [Point(p.x, -p.y) for p in g.vertices] if refl else list(g.vertices)
        if any(p not in Sset for p in base): continue
        cur = [pos[p] for p in base]
        for _ in range(6):
            maps.append(cur)
            cur = [pos[rot60(g.vertices[v])] for v in cur]
    K = 4
    X = lambda v, c: 1 + v * K + c
    cnf = [[X(v, c) for c in range(K)] for v in range(n)]
    for u, v in g.edges():
        for c in range(K):
            cnf.append([-X(u, c), -X(v, c)])
    for v in range(n):
        for c in range(K):
            for e in range(c + 1, K):
                cnf.append([-X(v, c), -X(v, e)])
    s = Solver(name="m22", bootstrap_with=cnf)
    if not s.solve():
        print(f"  v{b}: n={n} *** NOT 4-COLOURABLE ***", flush=True); s.delete(); continue
    def ro():
        p = set(l for l in s.get_model() if l > 0)
        return [next(c for c in range(K) if X(v, c) in p) for v in range(n)]
    raw = [ro()]; rng = random.Random(5)
    while len(raw) < 8:
        s.add_clause([-X(v, raw[-1][v]) for v in rng.sample(range(n), 30)])
        s.set_phases([(1 if rng.random() < 0.25 else -1) * X(v, c)
                      for v in range(n) for c in range(K)])
        if not s.solve(): break
        raw.append(ro())
    s.delete()
    eff = [[c[pm[v]] for v in range(n)] for c in raw for pm in maps]
    buck = defaultdict(list)
    for v in range(n): buck[tuple(c[v] for c in eff)].append(v)
    cand = [(a, b2) for vs in buck.values() if len(vs) > 1
            for i, a in enumerate(vs) for b2 in vs[i+1:]]
    seenp, reps = set(), []
    for a, b2 in cand:
        if (min(a, b2), max(a, b2)) in seenp: continue
        for pm in maps:
            x, y = pm[a], pm[b2]
            seenp.add((min(x, y), max(x, y)))
        reps.append((a, b2))
    # the distance is free; only the unusual ones are worth a solver call
    odd = []
    for a, b2 in reps:
        dd = (U[a] - U[b2]).norm2()
        if dd.is_rational() and Fr(dd.c[0]) != COMMON:
            odd.append((a, b2, Fr(dd.c[0])))
    found = defaultdict(int)
    s2 = Solver(name="m22", bootstrap_with=cnf)
    for a, b2, q in odd:
        if not s2.solve(assumptions=[X(a, 0), X(b2, 1)]):
            found[q] += 1
    s2.delete()
    m = sum(len(a) for a in g.adj) // 2
    tag = ("  <<<< " + ", ".join(f"{q} (radical {radical(q)})" for q in found)) if found else ""
    print(f"  v{b}: n={n} deg={2.0*m/n:.2f} reps={len(reps)} "
          f"unusual-distance reps={len(odd)} RARE FORCED={dict(found)}{tag}"
          f"   [{time.time()-t0:.0f}s]", flush=True)
