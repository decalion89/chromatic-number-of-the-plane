"""Hunt a 3.  The quantity IS the problem, and it has been measured once.

For a unit-distance graph G and a point p, let

        mu(G, p) = min over proper 5-colourings c of  |c(N(p))|

where N(p) is the set of graph points at distance exactly 1 from p.  Then

        mu = 5 somewhere   <=>   chi(R^2) >= 6
        mu <= 2 always     <=>   chi(R^2) = 5

both by the same argument: G + p has no 5-colouring exactly when N(p) is
forced to carry all five colours, and a vertex-critical 6-chromatic graph
minus a vertex is precisely such a G and p.

So mu is not a proxy for the problem, it is the problem, restated as a number
between 1 and 5.  It has been computed once in this project, on G, where it
came back 2 for all four hundred candidates examined; and now on one tuned
object, where it came back 2 again.  Nowhere else.

mu <= 2 is also what the geometry expects.  N(p) lies on the unit circle about
p, two of its points are adjacent only at 60 degrees, so the induced graph has
maximum degree 2 and even cycles only: chi(N(p)) = 2 on its own, always, in
every planar unit-distance graph.  Anything above 2 has to be imposed by the
rest of G, which is a five-colour forcing statement -- and every search for
one of those has come back empty.

A single 3 would therefore be the first crack in that: the first evidence that
the ambient graph can constrain a unit circle at all at five colours.  This
sweeps every graph in data/ rather than one, takes the richest candidate
neighbourhoods of each, and reports the distribution of mu.
"""
import sys, json, time, os
from fractions import Fraction as Fr
from collections import defaultdict, Counter
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
from hn.field import Field
from hn.geometry import Point, Rotation, _rot60
from hn.graph import build_graph
from pysat.solvers import Solver

ROOT = "/home/user/darwin-50/research/hadwiger-nelson"
t0 = time.time()
K = 5
FILES = ["five_247_b.json", "five_247.json", "five_tuned_4_1.json",
         "five_tuned_1_3.json", "five_twotune_small.json", "five_23_v7.json"]
best_overall = 0
for name in FILES:
    path = f"{ROOT}/data/{name}"
    if not os.path.exists(path):
        continue
    d = json.load(open(path))
    F = Field(tuple(d["field_generators"]))
    P = [Point(F.element([Fr(a, b) for a, b in x]),
               F.element([Fr(a, b) for a, b in y])) for x, y in d["points"]]
    g = build_graph(P); n = g.n
    if n > 5000:
        print(f"\n  {name}: n={n}, skipped (one colouring costs too much)",
              flush=True)
        continue
    m = sum(len(a) for a in g.adj) // 2
    print(f"\n  {name}  n={n} m={m} deg={2.0*m/n:.2f}   [{time.time()-t0:.0f}s]",
          flush=True)
    S = set(g.vertices); one = F.rational(Fr(1))
    r60 = _rot60(F)
    r30 = Rotation(F.sqrt(3) * F.rational(Fr(1, 2)), F.rational(Fr(1, 2)))
    deg_order = sorted(range(n), key=lambda v: -len(g.adj[v]))
    cand = set()
    for c in deg_order[:30]:
        for rot in (r60.about(g.vertices[c]), r30.about(g.vertices[c])):
            for p in g.vertices:
                z = rot(p)
                if z not in S:
                    cand.add(z)
    cells = defaultdict(list)
    for i, p in enumerate(g.vertices):
        cells[(int(p.fx // 1), int(p.fy // 1))].append(i)
    gx = [q.fx for q in g.vertices]; gy = [q.fy for q in g.vertices]
    scored, seen = [], set()
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
        if len(nb) >= 7:
            t = tuple(sorted(nb))
            if t not in seen:
                seen.add(t); scored.append((len(nb), t))
    scored.sort(key=lambda u: -u[0])
    print(f"    {len(scored)} neighbourhoods of size >= 7, largest "
          f"{scored[0][0] if scored else 0}   [{time.time()-t0:.0f}s]",
          flush=True)
    if not scored:
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
    if not s.solve():
        print(f"    *** refuses five ***", flush=True); s.delete(); continue
    print(f"    base colouring in {time.time()-t1:.0f}s", flush=True)
    # mu <= 2 is decided by ONE call: can three colours be absent at once?
    # A yes settles the candidate as mu <= 2 with no further work; only a NO
    # is interesting, and then the chain pins the value.
    hist = Counter()
    for k, nb in scored[:25]:
        s.conf_budget(3_000_000)
        r = s.solve_limited(assumptions=[-X(u, c) for u in nb for c in (0, 1, 2)])
        if r is True:
            hist[2] += 1
            continue
        if r is None:
            hist[-1] += 1
            continue
        # not 2 -- find the real value
        mu = 5
        for j in (2, 1):
            s.conf_budget(3_000_000)
            if s.solve_limited(assumptions=[-X(u, c) for u in nb
                                            for c in range(j)]) is True:
                mu = K - j
                break
        hist[mu] += 1
        best_overall = max(best_overall, mu)
        print(f"    *** |N|={k}: mu = {mu}, ABOVE TWO ***", flush=True)
        if mu == K:
            print("    *** BLOCKED POINT -- chi(R^2) >= 6 ***", flush=True)
            json.dump({"graph": name, "neighbourhood": list(nb)},
                      open(f"{ROOT}/data/blocked_point.json", "w"))
    s.delete()
    print(f"    mu distribution over the richest 25: {dict(hist)}"
          f"   [{time.time()-t0:.0f}s]", flush=True)
print(f"\n  best mu found anywhere: {max(best_overall, 2)} of 5"
      f"   [{time.time()-t0:.0f}s]", flush=True)
