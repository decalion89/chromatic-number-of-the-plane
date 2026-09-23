"""Does a SMALL carrier have a forced class of three?

The class on Sa[199] came from two RARE distances sharing a pivot, and that
carrier is 3025 points, so its one-pivot object is 12097 and its minimisation
is out of reach.  But a class of three needs only two forced pairs sharing a
vertex, and nothing says the two distances have to be rare.  A small carrier
with a class of three would give a small one-pivot object -- and those are the
only ones whose five-colour test can actually be asked.

So: for each carrier, confirm the forced pairs and look at the graph they
form.  Any vertex of degree two or more in it is the pivot of a class.
"""
import sys, time, random
from fractions import Fraction as Fr
from collections import defaultdict, Counter
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
from hn.degrey import build_Sa
from hn.geometry import Point, rotation_joining, _rot60
from hn.field import Field
from hn.graph import build_graph
from pysat.solvers import Solver

F = Field((3, 11, 247)); rot60 = _rot60(F); Sa = build_Sa(F)
g60 = rotation_joining(Fr(1), F)
BASES = [int(x) for x in sys.argv[1].split(",")]
def orb(p):
    out, q = [], p
    for _ in range(6):
        if q not in out: out.append(q)
        q = rot60(q)
    for q in list(out):
        r = Point(q.x, -q.y)
        if r not in out: out.append(r)
    return out
t0 = time.time()
for b in BASES:
    seen, U = set(Sa), list(Sa)
    for w in orb(Sa[b]):
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
    if not s.solve(): s.delete(); continue
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
    cand = [(x, y) for vs in buck.values() if len(vs) > 1
            for i, x in enumerate(vs) for y in vs[i+1:]]
    seenp, reps = set(), []
    for x, y in cand:
        if (min(x, y), max(x, y)) in seenp: continue
        for pm in maps:
            p1, p2 = pm[x], pm[y]
            seenp.add((min(p1, p2), max(p1, p2)))
        reps.append((x, y))
    # confirm, then look at the graph the forced pairs form
    edges = []
    s2 = Solver(name="m22", bootstrap_with=cnf)
    for x, y in reps:
        if not s2.solve(assumptions=[X(x, 0), X(y, 1)]):
            for pm in maps:
                edges.append((pm[x], pm[y]))
    s2.delete()
    deg = Counter()
    for x, y in edges: deg[x] += 1; deg[y] += 1
    hubs = [(v, d) for v, d in deg.items() if d >= 2]
    m = sum(len(a) for a in g.adj) // 2
    print(f"  v{b}: n={n} deg={2.0*m/n:.2f}  {len(reps)} forced reps, "
          f"{len(edges)} forced pairs, {len(hubs)} vertices in two or more"
          + ("   <<<< CLASS OF THREE" if hubs else "")
          + f"   [{time.time()-t0:.0f}s]", flush=True)
    for v, d in sorted(hubs, key=lambda t: -t[1])[:4]:
        nb = [y for x, y in edges if x == v] + [x for x, y in edges if y == v]
        ds = [str((U[v] - U[w]).norm2()) for w in nb[:4]]
        print(f"      v{v} is in {d} forced pairs, distances^2 {ds}", flush=True)
