"""The precondition for six is tightness at five, and it is cheap to measure.

Forcing at k colours needs chi >= k: if chi(H) < k there is a spare colour, so
the pivot can always be recoloured with it and nothing is forced.  That is
necessary, not sufficient.  What actually makes a pair forced is RIGIDITY: a
vertex whose neighbourhood already shows the other k-1 colours cannot move.

Sa reaches free@4 = 0.00 % -- in a 4-colouring not one of its 397 vertices has
a spare colour -- and that is why gluing it manufactures forcing.  Every
5-chromatic graph built here sits at free@5 = 11-13 %: at five colours they are
loose, and a loose graph cannot force anything.  So the target is not size, not
density and not rarity.  It is free@5 -> 0.

The measurement is O(n + m) on a single colouring, no solver call per vertex,
so hundreds of candidate unions can be ranked in the time one forced-pair
filter would take.

How Sa got tight is the recipe to copy: Sa is not one diamond, it is the union
of the 39-point diamond over a group of rotations about a common centre, so
every vertex lies in several copies at once.  The same move one rung up takes
the 803-point 5-chromatic graph as the diamond and unions it over a group --
which is exactly the group this project needs and does not yet have.  C12 costs
nothing here: cos 30 = sqrt3/2 and sin 30 = 1/2 are already in the field.
"""
import sys, time, json, random
from fractions import Fraction as Fr
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
from hn.field import Field
from hn.geometry import Point, Rotation, rotation_joining, _rot60
from hn.graph import build_graph
from pysat.solvers import Solver

ROOT = "/home/user/darwin-50/research/hadwiger-nelson"
t0 = time.time()
F = Field((3, 11, 247))
half = F.rational(Fr(1, 2))
d = json.load(open(ROOT + "/data/five_247_c.json"))
Z = [Point(F.element([Fr(a, b) for a, b in x]),
           F.element([Fr(a, b) for a, b in y])) for x, y in d["points"]]
gZ = build_graph(Z)
print(f"diamond n={gZ.n} m={sum(len(a) for a in gZ.adj)//2}", flush=True)

r60 = _rot60(F)
r30 = Rotation(F.sqrt(3) * half, half)
r64_9 = rotation_joining(Fr(64, 9), F)          # the spindle letter it was built with
def powers(r, k):
    out, cur = [], Rotation(F.rational(Fr(1)), F.rational(Fr(0)))
    for _ in range(k):
        out.append(cur)
        cur = Rotation(cur.cos * r.cos - cur.sin * r.sin,
                       cur.cos * r.sin + cur.sin * r.cos)
    return out
FAMILIES = [("C6", powers(r60, 6)), ("C12", powers(r30, 12)),
            ("C3", powers(r60, 6)[::2]), ("spindle64/9 x6", powers(r64_9, 6))]

def free_at(g, col, K):
    n = g.n; fr = 0
    for v in range(n):
        seen = set()
        for u in g.adj[v]:
            seen.add(col[u])
        if len(seen) < K - 1:
            fr += 1
    return fr / n

def colour(g, K, budget=8_000_000):
    n = g.n
    X = lambda v, c: 1 + v * K + c
    cl = [[X(v, c) for c in range(K)] for v in range(n)]
    for v in range(n):
        for a in range(K):
            for b in range(a + 1, K):
                cl.append([-X(v, a), -X(v, b)])
    for x, y in g.edges():
        for c in range(K):
            cl.append([-X(x, c), -X(y, c)])
    s = Solver(name="cd19", bootstrap_with=cl)
    s.conf_budget(budget)
    r = s.solve_limited()
    if r is not True:
        s.delete(); return None, r
    pos = set(l for l in s.get_model() if l > 0); s.delete()
    return [next(c for c in range(K) if X(v, c) in pos) for v in range(n)], True

deg = sorted(range(gZ.n), key=lambda v: -len(gZ.adj[v]))
ORIGIN = Point(F.rational(Fr(0)), F.rational(Fr(0)))
centres = [("origin", ORIGIN)] + [(f"v{v}", Z[v]) for v in deg[:3]] \
        + [(f"lowv{v}", Z[v]) for v in deg[-2:]]

rows = []
for cname, c in centres:
    for fname, fam in FAMILIES:
        seen, V = set(), []
        for r in fam:
            rot = r.about(c)
            for p in Z:
                q = rot(p)
                if q not in seen:
                    seen.add(q); V.append(q)
        if len(V) > 11000:
            print(f"  {cname:8s} {fname:14s} n={len(V)}  (skipped, too big)",
                  flush=True)
            continue
        g = build_graph(V); n = g.n; m = sum(len(a) for a in g.adj) // 2
        col, st = colour(g, 5)
        if col is None:
            print(f"  {cname:8s} {fname:14s} n={n} m={m} deg={2.0*m/n:.2f}"
                  f"  *** 5-colour {'UNSAT' if st is False else 'budget out'} ***"
                  f"   [{time.time()-t0:.0f}s]", flush=True)
            if st is False:
                json.dump({"field_generators": list(F.gens), "n": n, "m": m,
                           "centre": cname, "family": fname,
                           "points": [[[[q.numerator, q.denominator] for q in p.x.c],
                                       [[q.numerator, q.denominator] for q in p.y.c]]
                                      for p in g.vertices]},
                          open(ROOT + "/data/six_candidate.json", "w"))
                print("  *** written data/six_candidate.json ***", flush=True)
            rows.append((1.0, cname, fname, n, m, "hard"))
            continue
        f5 = free_at(g, col, 5)
        rows.append((f5, cname, fname, n, m, ""))
        print(f"  {cname:8s} {fname:14s} n={n} m={m} deg={2.0*m/n:.2f}"
              f"  free@5={100*f5:.2f}%   [{time.time()-t0:.0f}s]", flush=True)

rows.sort()
print("\n  tightest first:", flush=True)
for f5, cn, fn, n, m, tag in rows[:8]:
    print(f"    free@5={100*f5:6.2f}%  {cn:8s} {fn:14s} n={n} deg={2.0*m/n:.2f} {tag}",
          flush=True)
