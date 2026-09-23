"""A group that is not Sa's: D12, order 24, free in the field already used.

Everything symmetric so far has used C6 or D6 -- Sa's own group.  The field
contains cos 30 = sqrt3/2 with no adjunction at all, so the rotation of order
twelve is free, and D12 has order 24.  It was tried once, by taking the D12
orbit of a FINISHED graph, and it loosened the census; but that is exactly the
mistake the symmetric construction identified -- imposing a group on an object
not shaped for it.  Building with it from the start has not been tried.

Two reasons to want it.  The carrier is a genuinely different object, not a
denser version of the same one.  And the filter gets a factor of 24 rather
than 12 for free, which is what decides whether the forced-pair test at five
is affordable at all.

    Sa12 = Sa u rot30(Sa)       C12-invariant, and D12 since Sa is
                                reflection-closed
"""
import sys, time, random
from fractions import Fraction as Fr
from collections import defaultdict
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
from hn.degrey import build_Sa
from hn.geometry import Point, Rotation, rotation_joining, _rot60
from hn.field import Field
from hn.graph import build_graph
from pysat.solvers import Solver

F = Field((3, 11, 247))
half = F.rational(Fr(1, 2))
rot30 = Rotation(F.sqrt(3) * half, half)          # cos 30 = sqrt3/2, sin 30 = 1/2
rot60 = _rot60(F)
Sa = build_Sa(F)
t0 = time.time()
seen, S12 = set(), []
for p in Sa:
    for q in (p, rot30(p)):
        if q not in seen: seen.add(q); S12.append(q)
g0 = build_graph(S12); m0 = sum(len(a) for a in g0.adj) // 2
SS = set(g0.vertices)
print(f"Sa12 = Sa u rot30(Sa): n={g0.n} m={m0} deg={2.0*m0/g0.n:.2f}", flush=True)
print(f"  rot30-invariant: {all(rot30(p) in SS for p in g0.vertices)}   "
      f"reflection-invariant: {all(Point(p.x, -p.y) in SS for p in g0.vertices)}"
      f"   [{time.time()-t0:.0f}s]", flush=True)
K = 4
def analyse(P, tag, nraw=8):
    g = build_graph(P); n = g.n; m = sum(len(a) for a in g.adj) // 2
    pos = {p: i for i, p in enumerate(g.vertices)}
    Sset = set(g.vertices)
    maps = []
    for refl in (False, True):
        base = [Point(p.x, -p.y) for p in g.vertices] if refl else list(g.vertices)
        if any(p not in Sset for p in base): continue
        cur = [pos[p] for p in base]
        for _ in range(12):
            maps.append(cur)
            cur = [pos[rot30(g.vertices[v])] for v in cur]
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
        print(f"  {tag}: n={n} m={m} *** NOT 4-COLOURABLE ***", flush=True)
        s.delete(); return None
    def ro():
        p = set(l for l in s.get_model() if l > 0)
        return [next(c for c in range(K) if X(v, c) in p) for v in range(n)]
    raw = [ro()]; rng = random.Random(5)
    while len(raw) < nraw:
        s.add_clause([-X(v, raw[-1][v]) for v in rng.sample(range(n), 30)])
        s.set_phases([(1 if rng.random() < 0.25 else -1) * X(v, c)
                      for v in range(n) for c in range(K)])
        if not s.solve(): break
        raw.append(ro())
    s.delete()
    eff = [[c[pm[v]] for v in range(n)] for c in raw for pm in maps]
    free = min(sum(1 for v in range(n)
                   if len({col[u] for u in g.adj[v]} | {col[v]}) < K) for col in raw)
    buck = defaultdict(list)
    for v in range(n): buck[tuple(c[v] for c in eff)].append(v)
    cand = [(a, b) for vs in buck.values() if len(vs) > 1
            for i, a in enumerate(vs) for b in vs[i+1:]]
    seenp, reps = set(), []
    for a, b in cand:
        if (min(a, b), max(a, b)) in seenp: continue
        for pm in maps:
            x, y = pm[a], pm[b]
            seenp.add((min(x, y), max(x, y)))
        reps.append((a, b))
    byd = defaultdict(int); tot = 0
    s2 = Solver(name="m22", bootstrap_with=cnf)
    for a, b in reps:
        if not s2.solve(assumptions=[X(a, 0), X(b, 1)]):
            tot += 1
            dd = (P[a] - P[b]).norm2()
            byd[str(dd) if dd.is_rational() else "irr"] += 1
    s2.delete()
    print(f"  {tag}: n={n} m={m} deg={2.0*m/n:.2f} group {len(maps)} "
          f"free@4={100.0*free/n:.2f}% cand={len(cand)} reps={len(reps)} "
          f"FORCED={tot} {dict(sorted(byd.items(), key=lambda kv: -kv[1])[:4])}"
          f"   [{time.time()-t0:.0f}s]", flush=True)
    return g
analyse(S12, "Sa12 itself")
# glue over a D12 orbit of centres
def orb12(p):
    out, q = [], p
    for _ in range(12):
        if q not in out: out.append(q)
        q = rot30(q)
    for q in list(out):
        r = Point(q.x, -q.y)
        if r not in out: out.append(r)
    return out
g60 = rotation_joining(Fr(1), F)
for vi in (25, 265):
    cs = orb12(S12[vi])
    seen2, U = set(S12), list(S12)
    for w in cs:
        rot = g60.about(w)
        for p in S12:
            q = rot(p)
            if q not in seen2: seen2.add(q); U.append(q)
    analyse(U, f"Sa12 glued over the D12 orbit of v{vi} ({len(cs)} centres)")
