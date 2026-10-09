#!/usr/bin/env python3
"""Recount the vertices/edges quoted for the union of the nine rotated copies of the 76-vertex graph
(data file q11.json of the repository, read as data only).  Exact arithmetic in Q(sqrt 11)(i)."""
import json
import os
import sys
from fractions import Fraction as Fr

PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "..", "..", "quadratic_planes", "q11.json")
g = json.load(open(PATH))
d, D = g["d"], g["D"]
assert d == 11
FAILS = []


def check(name, cond, detail=""):
    print(("PASS  " if cond else "FAIL  ") + name + (("   [" + detail + "]") if detail else ""))
    if not cond:
        FAILS.append(name)


pts = [tuple(Fr(x, D) for x in p) for p in g["points"]]  # (x0, x1, y0, y1): x = x0 + x1 s, y = y0 + y1 s


def is_unit(p, q):
    a0, a1, b0, b1 = (q[0] - p[0], q[1] - p[1], q[2] - p[2], q[3] - p[3])
    return a0 * a0 + d * a1 * a1 + b0 * b0 + d * b1 * b1 == 1 and a0 * a1 + b0 * b1 == 0


n0 = len(pts)
E0 = {tuple(sorted(e)) for e in g["edges"]}
allpairs0 = {(i, j) for i in range(n0) for j in range(i + 1, n0) if is_unit(pts[i], pts[j])}
check("76 vertices; the stored edges are exactly the unit-distance pairs", n0 == 76 and E0 == allpairs0,
      "%d edges" % len(E0))


def rot(gam, p):
    g0, g1 = gam
    x0, x1, y0, y1 = p
    return (g0 * x0 - g1 * y0, g0 * x1 - g1 * y1, g0 * y0 + g1 * x0, g0 * y1 + g1 * x1)


def cmul(a, b):
    return (a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0])


rho, rhoi = (Fr(3, 5), Fr(4, 5)), (Fr(3, 5), Fr(-4, 5))
sig, sigi = (Fr(5, 13), Fr(12, 13)), (Fr(5, 13), Fr(-12, 13))
ONE = (Fr(1), Fr(0))
rots = []
for j in (-1, 0, 1):
    for l in (-1, 0, 1):
        r = ONE
        r = cmul(r, rho if j == 1 else rhoi if j == -1 else ONE)
        r = cmul(r, sig if l == 1 else sigi if l == -1 else ONE)
        rots.append(r)
idx = {}
copy_edges = set()
for r in rots:
    loc = []
    for p in pts:
        q = rot(r, p)
        if q not in idx:
            idx[q] = len(idx)
        loc.append(idx[q])
    for (a, b) in E0:
        x, y = loc[a], loc[b]
        copy_edges.add((min(x, y), max(x, y)))
P = list(idx)
n = len(P)
allpairs = {(i, j) for i in range(n) for j in range(i + 1, n) if is_unit(P[i], P[j])}
check("union of the nine copies has 628 vertices", n == 628, "%d" % n)
check("the nine copies have 1494 edges", len(copy_edges) == 1494, "%d" % len(copy_edges))
check("1506 edges when all pairs at distance 1 are joined", len(allpairs) == 1506, "%d" % len(allpairs))
print()
print("FAILED:", FAILS if FAILS else "none")
sys.exit(1 if FAILS else 0)
