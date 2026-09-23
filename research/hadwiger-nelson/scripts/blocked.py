"""The blocked point: one SAT call per candidate, and equivalent to the answer.

A forced pair at five is SUFFICIENT for chi(R^2) >= 6 -- one rotation closes it
-- but nobody knows one exists, and nothing in this project has produced one.
A BLOCKED POINT is equivalent, in both directions:

  (<=)  if G is a 5-colourable unit-distance graph and p a point whose unit
        circle meets G in a set N(p) that uses all five colours in EVERY proper
        5-colouring of G, then G + p has no 5-colouring.
  (=>)  if chi(R^2) >= 6, take a 6-chromatic unit-distance graph, make it
        vertex-critical and delete any vertex p.  What remains is 5-colourable
        and N(p) must use all five colours in every 5-colouring, or p could be
        put back.

So this is not a shortcut to the answer; it IS the answer, asked directly.  Its
virtue is cost.  Colour symmetry collapses the test to a single call: if some
colour can be missing from N(p), then colour 0 can, so

        solve(assumptions = [-X(u, 0) for u in N(p)])

UNSAT means blocked.  No refutation of a whole graph, no filter over millions
of pairs -- one assumption call per candidate point, on a warm solver.

Candidates have to lie in the field, so they are taken as images of the graph's
own vertices under rotations the field already holds -- rot60 about each of a
sample of vertices, and the tuning rotations.  Each candidate is scored by how
many graph vertices sit at distance exactly 1 from it, and only the richest are
worth a call: a blocked point needs at least five neighbours, and since a
neighbourhood is bipartite it needs many more than five to have any chance of
carrying five colours.

Reported graded rather than yes/no: the MINIMUM number of colours the solver
can leave on N(p).  A blocked point is that minimum reaching 5.  On G the
earlier pass found every candidate sitting at 2.
"""
import sys, json, time
from fractions import Fraction as Fr
from collections import defaultdict
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
from hn.field import Field
from hn.geometry import Point, Rotation, _rot60
from hn.graph import build_graph
from pysat.solvers import Solver

ROOT = "/home/user/darwin-50/research/hadwiger-nelson"
t0 = time.time()
K = 5
FILES = ["five_tuned_1_1.json", "five_twotune_small.json", "five_247_c.json"]
for name in FILES:
    try:
        d = json.load(open(f"{ROOT}/data/{name}"))
    except FileNotFoundError:
        continue
    F = Field(tuple(d["field_generators"]))
    P = [Point(F.element([Fr(a, b) for a, b in x]),
               F.element([Fr(a, b) for a, b in y])) for x, y in d["points"]]
    g = build_graph(P); n = g.n; m = sum(len(a) for a in g.adj) // 2
    print(f"\n  {name}  n={n} m={m} deg={2.0*m/n:.2f}   [{time.time()-t0:.0f}s]",
          flush=True)
    S = set(g.vertices)
    idx = {p: i for i, p in enumerate(g.vertices)}
    one = F.rational(Fr(1))
    r60 = _rot60(F)
    r30 = Rotation(F.sqrt(3) * F.rational(Fr(1, 2)), F.rational(Fr(1, 2)))
    # candidates: images of the vertex set under rotations about a spread of
    # centres, which keeps every coordinate inside the field
    deg_order = sorted(range(n), key=lambda v: -len(g.adj[v]))
    cand = defaultdict(int)
    for c in deg_order[:40]:
        for rot in (r60.about(g.vertices[c]), r30.about(g.vertices[c])):
            for p in g.vertices:
                z = rot(p)
                if z not in S:
                    cand[z] += 0          # register it
    print(f"    {len(cand)} candidate points off the graph"
          f"   [{time.time()-t0:.0f}s]", flush=True)
    # score each by how many graph vertices sit at distance exactly 1
    cells = defaultdict(list)
    for i, p in enumerate(g.vertices):
        cells[(int(p.fx // 1), int(p.fy // 1))].append(i)
    # Float first, exact second.  Scoring 126 090 candidates by exact field
    # arithmetic alone means millions of multiplications in an 8-dimensional
    # field and takes longer than every solver call put together; a float
    # distance discards all but a handful before any exact work happens, and
    # the exact check still decides every survivor, so nothing is assumed.
    gx = [q.fx for q in g.vertices]
    gy = [q.fy for q in g.vertices]
    scored = []
    for z in cand:
        zx, zy = z.fx, z.fy
        cx, cy = int(zx // 1), int(zy // 1)
        nb = []
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                for i in cells.get((cx + dx, cy + dy), ()):
                    ex, ey = gx[i] - zx, gy[i] - zy
                    if abs(ex * ex + ey * ey - 1.0) < 1e-9 and \
                       (g.vertices[i] - z).norm2() == one:
                        nb.append(i)
        if len(nb) >= 5:
            scored.append((len(nb), z, tuple(sorted(nb))))
    scored.sort(key=lambda t: -t[0])
    seen = set()
    uniq = []
    for k, z, nb in scored:
        if nb in seen:
            continue
        seen.add(nb); uniq.append((k, z, nb))
    print(f"    {len(uniq)} distinct neighbourhoods of size >= 5, "
          f"largest {uniq[0][0] if uniq else 0}   [{time.time()-t0:.0f}s]",
          flush=True)
    if not uniq:
        continue
    X = lambda v, c: 1 + v * K + c
    cnf = [[X(v, c) for c in range(K)] for v in range(n)]
    for v in range(n):
        for a in range(K):
            for b in range(a + 1, K):
                cnf.append([-X(v, a), -X(v, b)])
    for x, y in g.edges():
        for c in range(K):
            cnf.append([-X(x, c), -X(y, c)])
    s = Solver(name="cd19", bootstrap_with=cnf)
    t1 = time.time()
    assert s.solve(), "must be 5-colourable"
    print(f"    base colouring in {time.time()-t1:.0f}s", flush=True)
    # Each assumption call restarts cadical, so they cost tens of seconds
    # rather than nothing, and four hundred of them is hours.  The question
    # only lives at the top of the list anyway: a blocked point needs five
    # colours on a set that is bipartite on its own, so it needs MANY
    # neighbours, and a candidate with six cannot possibly carry five.  Take
    # the richest fifty, bound every call, and report the grade rather than a
    # yes/no -- how many colours the solver could not avoid putting on N(p).
    blocked, worst = [], 0
    grades = []
    for k, z, nb in uniq[:50]:
        s.conf_budget(2_000_000)
        r = s.solve_limited(assumptions=[-X(u, 0) for u in nb])
        if r is None:
            print(f"    |N|={k}: budget out", flush=True)
            continue
        if r is False:
            blocked.append((k, [ [[c.numerator, c.denominator] for c in z.x.c],
                                 [[c.numerator, c.denominator] for c in z.y.c] ]))
            print(f"    *** BLOCKED POINT with {k} neighbours -- "
                  f"chi(R^2) >= 6 ***", flush=True)
        else:
            pos = set(l for l in s.get_model() if l > 0)
            used = len({next(c for c in range(K) if X(u, c) in pos) for u in nb})
            worst = max(worst, used)
            grades.append((k, used))
            print(f"    |N|={k}: a colouring leaves {used} colours on N(p)"
                  f"   [{time.time()-t0:.0f}s]", flush=True)
    s.delete()
    print(f"    {len(grades)} graded, {len(blocked)} blocked; the most a "
          f"colouring was forced to put on N(p) is {worst} of 5"
          f"   [{time.time()-t0:.0f}s]", flush=True)
    if blocked:
        json.dump({"graph": name, "blocked": blocked},
                  open(f"{ROOT}/data/blocked_{name}", "w"))
        print(f"    *** written data/blocked_{name} ***", flush=True)
