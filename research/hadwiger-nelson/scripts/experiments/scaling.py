"""Does coupling SCALE?  The question that decides whether this reaches six.

Six copies of the vertex-critical 803, coupled on one hexagon, need 207 370
conflicts to 5-colour; the same six copies placed apart need 205.  Same edges,
fewer vertices, a thousand times the difficulty -- so the coupling is the lever.

If that lever is multiplicative in the number of coupled copies, then enough of
them drive the number of 5-colourings to one -- every non-adjacent pair forced,
and a single spindle-able forced pair closes the problem -- or straight to zero,
which is a 6-chromatic unit-distance graph with nothing further needed.  So
measure the curve.

Copies sit on a patch of the triangular lattice, each with its densest vertex on
one lattice point and turned by (index x 60 degrees) about it, the turn that was
fifteen times tighter than translation:

     1   the centre alone
     7   centre + hexagon
    13   + the six points at distance sqrt3
    19   + the six points at distance 2

and each is compared with the same number of copies placed far apart, which is
what independence alone would give.  A conflict budget keeps any one solve from
eating the machine; running out of it is reported as such, never as a result.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time, json, random, statistics
from fractions import Fraction as Fr
from collections import defaultdict
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
from hn.field import Field
from hn.geometry import Point
from hn.graph import build_graph
from pysat.solvers import Solver

ROOT = HN_DIR
t0 = time.time(); K = 5
BUDGET = int(sys.argv[1]) if len(sys.argv) > 1 else 30_000_000
SIZES = [int(x) for x in sys.argv[2:]] or [1, 7, 13, 19]
d = json.load(open(f"{ROOT}/data/five_247_c.json"))
F = Field(tuple(d["field_generators"]))
P = [Point(F.element([Fr(a, b) for a, b in x]),
           F.element([Fr(a, b) for a, b in y])) for x, y in d["points"]]
g0 = build_graph(P)
adj0 = defaultdict(set)
for x, y in g0.edges(): adj0[x].add(y); adj0[y].add(x)
v0 = max(range(g0.n), key=lambda v: len(adj0[v]))
r3 = F.sqrt(3); half = F.rational(Fr(1, 2)); c60, s60 = half, r3 * half
zero, one = F.rational(0), F.rational(1)
def rot(x, y, k):
    for _ in range(k % 6): x, y = x * c60 - y * s60, x * s60 + y * c60
    return x, y
sites = [(zero, zero)]
for k in range(6): sites.append(rot(one, zero, k))
a3x, a3y = F.rational(Fr(3, 2)), r3 * half                  # distance sqrt3
for k in range(6): sites.append(rot(a3x, a3y, k))
for k in range(6): sites.append(rot(F.rational(2), zero, k))  # distance 2

def solve_stats(pts):
    G = build_graph(pts); n = G.n; E = list(G.edges())
    X = lambda v, c: 1 + v * K + c
    cnf = [[X(v, c) for c in range(K)] for v in range(n)]
    for a, b in E:
        for c in range(K): cnf.append([-X(a, c), -X(b, c)])
    s = Solver(name="cd19", bootstrap_with=cnf)
    rng = random.Random(n)
    v = rng.randrange(n); c = rng.randrange(K)
    s.conf_budget(BUDGET)
    r = s.solve_limited(assumptions=[X(v, c)])
    st = s.accum_stats(); s.delete()
    return n, len(E), r, st.get("conflicts", 0)

bx, by = g0.vertices[v0].x, g0.vertices[v0].y
for m in SIZES:
    coupled = {}
    for idx, (sx, sy) in enumerate(sites[:m]):
        cc, ss = rot(one, zero, idx)
        for p in P:
            dx, dy = p.x - bx, p.y - by
            q = Point(sx + dx * cc - dy * ss, sy + dx * ss + dy * cc)
            coupled[(round(q.fx, 9), round(q.fy, 9))] = q
    apart = []
    for idx in range(m):
        off = F.rational(100 * idx)
        apart += [Point(p.x + off, p.y) for p in P]
    n1, e1, r1, c1 = solve_stats(list(coupled.values()))
    n2, e2, r2, c2 = solve_stats(apart)
    tag = {True: "5-colourable", False: "*** NOT 5-COLOURABLE ***",
           None: "budget out"}[r1]
    print(f"  {m:>2d} copies  coupled: n={n1:<6d} e={e1:<7d} conflicts {c1:>10d} "
          f"{tag:<24s} | apart: n={n2:<6d} conflicts {c2:>6d}  "
          f"ratio {c1/max(1,c2):>9.1f}   [{time.time()-t0:.0f}s]", flush=True)
    if r1 is False:
        json.dump({"copies": m, "field_generators": list(F.gens),
                   "points": [[[[t.numerator, t.denominator] for t in q.x.c],
                               [[t.numerator, t.denominator] for t in q.y.c]]
                              for q in coupled.values()]},
                  open(f"{ROOT}/data/scaling_hit_{m}.json", "w"))
        print("  written -- VERIFY INDEPENDENTLY BEFORE BELIEVING", flush=True)
        break
