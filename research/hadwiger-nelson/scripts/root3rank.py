"""Rank the sqrt3 pairs by the target they present, and ask the best directly.

Forcing u and v apart at five colours is exactly the contracted graph refusing
five, and the contracted vertex is adjacent to N(u) u N(v).  So the size of
that UNION is the natural measure of how good a candidate the pair is: a pair
whose two unit circles between them touch thirty points is a far bigger target
than one touching fifteen.

That replaces the filter.  Sweeping 46 680 pairs with repeated colourings of a
6925-vertex graph costs a colouring per round and eliminates five per cent a
time -- hours to reach a testable set.  Ranking costs one pass over the
adjacency lists and no solving at all, and the pairs it puts on top are the
ones worth a call.

Each call is the same one-line question: can u and v share colour 0?  UNSAT is
a pair forced apart at five, at distance sqrt3 -- the whole first step, and the
thing twenty-five of thirty candidate points in the 803-graph asked for.
"""
import sys, time, json
from fractions import Fraction as Fr
from collections import defaultdict
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
from hn.field import Field
from hn.geometry import Point
from hn.graph import build_graph
from pysat.solvers import Solver

ROOT = "/home/user/darwin-50/research/hadwiger-nelson"
t0 = time.time()
K = 5
for name in ("five_dense_2.json", "five_twotune_small.json"):
    d = json.load(open(f"{ROOT}/data/{name}"))
    F = Field(tuple(d["field_generators"]))
    P = [Point(F.element([Fr(a, b) for a, b in x]),
               F.element([Fr(a, b) for a, b in y])) for x, y in d["points"]]
    g = build_graph(P); n = g.n
    print(f"\n  {name}  n={n} deg="
          f"{2.0*sum(len(a) for a in g.adj)/2/n:.2f}   [{time.time()-t0:.0f}s]",
          flush=True)
    three = F.rational(Fr(3))
    cells = defaultdict(list)
    for i, p in enumerate(g.vertices):
        cells[(int(p.fx // 2), int(p.fy // 2))].append(i)
    gx = [q.fx for q in g.vertices]; gy = [q.fy for q in g.vertices]
    pairs = []
    for i in range(n):
        cx, cy = int(gx[i] // 2), int(gy[i] // 2)
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                for j in cells.get((cx + dx, cy + dy), ()):
                    if j <= i:
                        continue
                    ex, ey = gx[j] - gx[i], gy[j] - gy[i]
                    if abs(ex * ex + ey * ey - 3.0) < 1e-9 and \
                       (g.vertices[i] - g.vertices[j]).norm2() == three:
                        pairs.append((i, j))
    ranked = sorted(pairs, key=lambda ab: -len(set(g.adj[ab[0]]) |
                                               set(g.adj[ab[1]])))
    if not ranked:
        continue
    top = len(set(g.adj[ranked[0][0]]) | set(g.adj[ranked[0][1]]))
    print(f"    {len(pairs)} pairs at sqrt3; largest |N(u) u N(v)| = {top}"
          f"   [{time.time()-t0:.0f}s]", flush=True)
    X = lambda v, c: 1 + v * K + c
    cnf = [[X(v, c) for c in range(K)] for v in range(n)]
    for x, y in g.edges():
        for c in range(K):
            cnf.append([-X(x, c), -X(y, c)])
    s = Solver(name="cd19", bootstrap_with=cnf)
    t1 = time.time()
    if not s.solve():
        print("    refuses five", flush=True); s.delete(); continue
    print(f"    base colouring in {time.time()-t1:.0f}s", flush=True)
    forced, done, out = [], 0, 0
    for (a, b) in ranked[:400]:
        done += 1
        s.conf_budget(3_000_000)
        r = s.solve_limited(assumptions=[X(a, 0), X(b, 0)])
        if r is False:
            u = len(set(g.adj[a]) | set(g.adj[b]))
            forced.append((a, b, u))
            print(f"    *** FORCED APART AT FIVE: {a},{b}, "
                  f"|N(u) u N(v)| = {u} ***", flush=True)
            json.dump({"graph": name, "forced": forced},
                      open(f"{ROOT}/data/root3_forced.json", "w"))
        elif r is None:
            out += 1
        if done % 50 == 0:
            print(f"      ..{done}/400  {len(forced)} forced, {out} over budget"
                  f"   [{time.time()-t0:.0f}s]", flush=True)
    s.delete()
    print(f"    {done} of the richest pairs tested, {len(forced)} forced apart"
          f"   [{time.time()-t0:.0f}s]", flush=True)
