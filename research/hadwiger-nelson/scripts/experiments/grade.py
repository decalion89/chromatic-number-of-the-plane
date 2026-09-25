"""How many colours can a 5-colouring be FORCED to put on N(p)?  The real grade.

The previous run asked "can colour 0 be absent from N(p)?" and, when the answer
was yes, reported how many colours the resulting colouring happened to use.
That number can never exceed 4, because a colour had just been forbidden, so
"4 of 5" was a tautology dressed as progress.  It says the point is PLACEABLE,
which is the opposite of what is wanted.

The quantity that matters is the minimum over all proper 5-colourings of
|c(N(p))|.  A blocked point -- equivalently, chi(R^2) >= 6 -- is that minimum
reaching 5.  Colour symmetry makes it cheap to compute exactly: if j colours
can be simultaneously absent from N(p) in some colouring, then colours
0..j-1 can, so the test is a nested chain of assumption calls

        forbid {0}          SAT -> minimum <= 4
        forbid {0,1}        SAT -> minimum <= 3
        forbid {0,1,2}      SAT -> minimum <= 2
        forbid {0,1,2,3}    SAT -> minimum <= 1

and the first UNSAT in the chain pins the minimum exactly.  Four calls per
candidate instead of one, and the answer is a number between 1 and 5 rather
than a yes/no that reads backwards.

The earlier pass measured this on G and got 2 everywhere -- a blocked point is
not one colour away there, it is three.  This asks the same question of the
densest objects built since.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, json, time
from fractions import Fraction as Fr
from collections import defaultdict
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
from hn.field import Field
from hn.geometry import Point, Rotation, _rot60
from hn.graph import build_graph
from pysat.solvers import Solver

ROOT = HN_DIR
t0 = time.time()
K = 5
for name in ("five_tuned_1_1.json", "five_247_c.json"):
    d = json.load(open(f"{ROOT}/data/{name}"))
    F = Field(tuple(d["field_generators"]))
    P = [Point(F.element([Fr(a, b) for a, b in x]),
               F.element([Fr(a, b) for a, b in y])) for x, y in d["points"]]
    g = build_graph(P); n = g.n; m = sum(len(a) for a in g.adj) // 2
    print(f"\n  {name}  n={n} m={m} deg={2.0*m/n:.2f}   [{time.time()-t0:.0f}s]",
          flush=True)
    S = set(g.vertices); one = F.rational(Fr(1))
    r60 = _rot60(F)
    r30 = Rotation(F.sqrt(3) * F.rational(Fr(1, 2)), F.rational(Fr(1, 2)))
    deg_order = sorted(range(n), key=lambda v: -len(g.adj[v]))
    cand = set()
    for c in deg_order[:40]:
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
        if len(nb) >= 6:
            t = tuple(sorted(nb))
            if t not in seen:
                seen.add(t); scored.append((len(nb), t))
    scored.sort(key=lambda u: -u[0])
    print(f"    {len(scored)} distinct neighbourhoods of size >= 6, "
          f"largest {scored[0][0] if scored else 0}   [{time.time()-t0:.0f}s]",
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
    assert s.solve()
    print(f"    base colouring in {time.time()-t1:.0f}s", flush=True)
    best = 0
    for k, nb in scored[:12]:
        forbidden = 0
        for j in range(1, K):
            ass = [-X(u, c) for u in nb for c in range(j)]
            s.conf_budget(3_000_000)
            r = s.solve_limited(assumptions=ass)
            if r is True:
                forbidden = j
            else:
                break
        mn = K - forbidden
        best = max(best, mn)
        tag = "  *** BLOCKED -- chi(R^2) >= 6 ***" if mn == K else ""
        print(f"    |N|={k}: minimum colours on N(p) is {mn} of 5{tag}"
              f"   [{time.time()-t0:.0f}s]", flush=True)
    print(f"    best over the richest 12: {best} of 5 -- a blocked point needs 5"
          f"   [{time.time()-t0:.0f}s]", flush=True)
