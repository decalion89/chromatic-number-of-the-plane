"""The whole dihedral group: reflections as well as rotations, two copies a site.

The tightness lever is overlap WITH a twist: 60-degree turns keep the copies on
the triangular lattice, so they share points, and put every shared point in a
different role in each copy, so the constraints are new.  Translation overlaps
without twisting and is twelve times looser; a generic angle twists without
overlapping and is as loose as copies that never touch.

The dihedral group of the lattice is exactly the set of moves that keep both, and
the reflections have not been used.  So put TWO copies on every hexagon vertex --
one turned by i * 60 degrees, one reflected across the line through the centre
and then turned -- and compare with the six-copy 60-degree graph, whose median is
2 415 157 conflicts.  A reflection is (x, y) -> (x, -y) about the axis, rational,
so the field never moves.

Also the control that matters for any growth: the same twelve copies but with
only rotations, two different rotations per site, so that any gain can be put
on the reflections rather than on having twelve copies instead of six.
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
MODES = sys.argv[1:] or ["dihedral", "rot2"]
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
hexa = [rot(one, zero, k) for k in range(6)]
bx, by = g0.vertices[v0].x, g0.vertices[v0].y

def build(mode):
    pts = {}
    for i, (hxv, hyv) in enumerate(hexa):
        variants = [(i, False)]
        if mode == "dihedral": variants.append((i, True))
        elif mode == "rot2": variants.append((i + 3, False))
        for k, refl in variants:
            cc, ss = rot(one, zero, k)
            for p in P:
                dx, dy = p.x - bx, p.y - by
                if refl: dy = -dy
                q = Point(hxv + dx * cc - dy * ss, hyv + dx * ss + dy * cc)
                pts[(round(q.fx, 9), round(q.fy, 9))] = q
    hk = {(round(float(a), 9), round(float(b), 9)) for a, b in hexa}
    keep = [p for k, p in pts.items()
            if not (abs(p.fx*p.fx + p.fy*p.fy - 1.0) < 1e-9 and k not in hk)]
    return build_graph(keep)

for mode in MODES:
    G = build(mode); n = G.n; E = list(G.edges())
    X = lambda v, c: 1 + v * K + c
    cnf = [[X(v, c) for c in range(K)] for v in range(n)]
    for a, b in E:
        for c in range(K): cnf.append([-X(a, c), -X(b, c)])
    rng = random.Random(0); out = []
    for _ in range(3):
        s = Solver(name="cd19", bootstrap_with=cnf)
        v = rng.randrange(n); c = rng.randrange(K)
        s.conf_budget(60_000_000)
        r = s.solve_limited(assumptions=[X(v, c)])
        st = s.accum_stats(); s.delete()
        out.append((st.get("conflicts", 0), r))
    conf = [o[0] for o in out]
    print(f"  {mode:<9s} 12 copies: n={n} edges={len(E)} shared "
          f"{12*803 - n} conflicts {conf} median {statistics.median(conf):.0f} "
          f"colourable {[o[1] for o in out]}   [{time.time()-t0:.0f}s]", flush=True)
    if any(o[1] is False for o in out):
        json.dump({"mode": mode, "field_generators": list(F.gens),
                   "points": [[[[t.numerator, t.denominator] for t in q.x.c],
                               [[t.numerator, t.denominator] for t in q.y.c]]
                              for q in G.vertices]},
                  open(f"{ROOT}/data/dihedral_hit.json", "w"))
        print("  written -- VERIFY INDEPENDENTLY BEFORE BELIEVING", flush=True)
