"""de Grey's bite, applied for the first time to a graph that is already tight.

Coupling copies on a lattice saturates, because each copy is tied to the rest
only through a neighbourhood and brings far more freedom than that removes.  A
coupling that touches a whole ring does better: rotate about c by
theta = 2 arcsin(1/(2r)) and EVERY vertex at distance r from c lands exactly a
unit from its own image.  In G u rho_c(G) the entire ring is coupled at once.
That is de Grey's bite, the step that capped his radius-2 ring at two colours
and produced the forced pair he spindled.  On his own family at five colours it
fails -- the ring and its images form three disjoint ladders -- but it has never
been applied to a graph that is already tight.

The rotation lies in the field exactly when sqrt(4r^2 - 1) does:
cos = 1 - 1/(2r^2), sin = sqrt(4r^2 - 1)/(2r^2).  In Q(sqrt3, sqrt11, sqrt247)
that admits, among others,

    r^2 = 3      cos 5/6,   sin sqrt11/6      the Moser angle
    r^2 = 7      cos 13/14, sin 3 sqrt3/14
    r^2 = 13/4   cos 11/13, sin 4 sqrt3/13
    r^2 = 7/3    cos 11/14, sin 5 sqrt3/14

and the ring to bite is whichever is most POPULATED about some vertex, since the
bite couples every point on it.
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
NAME = sys.argv[1] if len(sys.argv) > 1 else "tight_hexagon_4159.json"
HOWMANY = int(sys.argv[2]) if len(sys.argv) > 2 else 3
d = json.load(open(f"{ROOT}/data/{NAME}"))
F = Field(tuple(d["field_generators"]))
P = [Point(F.element([Fr(a, b) for a, b in x]),
           F.element([Fr(a, b) for a, b in y])) for x, y in d["points"]]
G0 = build_graph(P); n0 = G0.n
hx = [q.fx for q in G0.vertices]; hy = [q.fy for q in G0.vertices]
r3, r11 = F.sqrt(3), F.sqrt(11)
Q = lambda a, b: F.rational(Fr(a, b))
BITES = {3.0: (Q(5, 6), r11 * Q(1, 6)),
         7.0: (Q(13, 14), r3 * Q(3, 14)),
         13/4: (Q(11, 13), r3 * Q(4, 13)),
         7/3: (Q(11, 14), r3 * Q(5, 14))}
for r2, (c, s) in BITES.items():
    assert c * c + s * s == F.rational(1)
cells = defaultdict(list)
for i in range(n0): cells[(int(hx[i] // 3), int(hy[i] // 3))].append(i)
best = []
for c0 in range(n0):
    cx, cy = int(hx[c0] // 3), int(hy[c0] // 3)
    cnt = defaultdict(int)
    for dx in (-1, 0, 1):
        for dy in (-1, 0, 1):
            for j in cells.get((cx + dx, cy + dy), ()):
                dd = (hx[j]-hx[c0])**2 + (hy[j]-hy[c0])**2
                for r2 in BITES:
                    if abs(dd - r2) < 1e-9: cnt[r2] += 1
    for r2, k in cnt.items(): best.append((k, c0, r2))
best.sort(reverse=True)
print(f"{NAME}: richest bitable rings {[(k, r2) for k, _, r2 in best[:8]]}   "
      f"[{time.time()-t0:.0f}s]", flush=True)

def tight(G):
    n = G.n; E = list(G.edges())
    X = lambda v, c: 1 + v * K + c
    cnf = [[X(v, c) for c in range(K)] for v in range(n)]
    for a, b in E:
        for c in range(K): cnf.append([-X(a, c), -X(b, c)])
    s = Solver(name="cd19", bootstrap_with=cnf)
    v = random.Random(n).randrange(n)
    s.conf_budget(40_000_000)
    r = s.solve_limited(assumptions=[X(v, 0)])
    st = s.accum_stats(); s.delete()
    return r, st.get("conflicts", 0), len(E)

done = set()
for k, c0, r2 in best:
    if (c0, r2) in done or len(done) >= HOWMANY: continue
    done.add((c0, r2))
    cc, ss = BITES[r2]
    cx, cy = G0.vertices[c0].x, G0.vertices[c0].y
    pts = {(round(q.fx, 9), round(q.fy, 9)): q for q in P}
    for p in P:
        dx, dy = p.x - cx, p.y - cy
        q = Point(cx + dx * cc - dy * ss, cy + dx * ss + dy * cc)
        pts[(round(q.fx, 9), round(q.fy, 9))] = q
    G = build_graph(list(pts.values()))
    r, conf, m = tight(G)
    tag = {True: "5-colourable", False: "*** NOT 5-COLOURABLE ***",
           None: "budget out"}[r]
    print(f"  bite ring r^2={r2} ({k} points) about {c0}: n={G.n} edges={m} "
          f"conflicts {conf}  {tag}   [{time.time()-t0:.0f}s]", flush=True)
    if r is False:
        json.dump({"source": NAME, "centre": c0, "r2": r2,
                   "field_generators": list(F.gens),
                   "points": [[[[t.numerator, t.denominator] for t in q.x.c],
                               [[t.numerator, t.denominator] for t in q.y.c]]
                              for q in G.vertices]},
                  open(f"{ROOT}/data/bite_hit.json", "w"))
        print("  written -- VERIFY INDEPENDENTLY BEFORE BELIEVING", flush=True)
        break
