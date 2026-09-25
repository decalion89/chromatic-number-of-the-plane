"""Thicken the ring by iterating the bite, and ask the only question.

The centre scan settled why the template stalls: G's richest doubly-usable
ring, at ANY centre, is twelve points -- exactly Sa's.  Four times the size
bought no extra coupling, so the bite at five colours has the same twelve
cross edges that sufficed at four, against a harder problem.

The bite itself is the lever.  rho = rotation_joining(4) has |1 - rho| = 1/2,
so a point at radius 2 and its image are at distance exactly one; and so are
its image and ITS image.  Iterating rho threads a PATH along the circle,
twelve of them, one per point of the D6 orbit, and every step adds a unit
edge that was not there before.  The angle arccos(7/8) is not a rational part
of a turn, so the paths never close and the ring thickens without limit.

At each step there is exactly one question worth asking and it takes one SAT
call: is the union still 5-colourable.  If it ever is not, that is chi >= 6
outright, with no spindle and no forced pair needed.  Failing that, test the
ring's antipodal pairs exactly -- they are the ones the template forces.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
from fractions import Fraction as Fr
from collections import Counter
from hn.degrey import build_G
from hn.geometry import DEGREY_FIELD as F, Point, rotation_joining
from hn.graph import build_graph
from pysat.solvers import Solver

k, D = 5, Fr(4)
t0 = time.time()
O = Point(F.zero(), F.zero())
rho = rotation_joining(D, F)
G = build_G(F, as_graph=False)
print(f"G: {len(G)} pts; rho = {rho.cos} + i {rho.sin}, "
      f"|1-rho| = 1/2 so the ring at radius 2 threads  "
      f"[{time.time()-t0:.0f}s]", flush=True)

U, seen = list(G), set(G)
cur = list(G)
for m in range(1, 13):
    cur = [rho(p) for p in cur]
    added = 0
    for p in cur:
        if p not in seen:
            seen.add(p)
            U.append(p)
            added += 1
    ring = sum(1 for p in U
               if (lambda d: all(x == 0 for x in d.c[1:]) and d.c[0] == D)
               (p.norm2()))
    g = build_graph(U)
    E = sorted(set((min(a, b), max(a, b)) for a, b in g.edges()))
    cls = [[1 + v * k + c for c in range(k)] for v in range(g.n)]
    for a, b in E:
        for c in range(k):
            cls.append([-(1 + a * k + c), -(1 + b * k + c)])
    sv = Solver(name="cd15", bootstrap_with=cls)
    ok = sv.solve()
    print(f"m={m}: +{added} -> {g.n} pts, {len(E)} edges, "
          f"{ring} on the D=4 ring, "
          f"{'5-colourable' if ok else '*** NOT 5-COLOURABLE ***'}"
          f"  [{time.time()-t0:.0f}s]", flush=True)
    if not ok:
        print("*** chi >= 6 ***", flush=True)
        break
    idx = {p: i for i, p in enumerate(U)}
    anti, hits = [], []
    for p in U:
        d = p.norm2()
        if all(x == 0 for x in d.c[1:]) and d.c[0] == D:
            q = Point(-p.x, -p.y)
            if q in idx and idx[p] < idx[q]:
                anti.append((idx[p], idx[q]))
    for i, j in anti:
        if not sv.solve(assumptions=[1 + i * k, -(1 + j * k)]):
            hits.append((i, j))
    print(f"    {len(anti)} antipodal pairs at squared distance {4*D}, "
          f"{len(hits)} FORCED  [{time.time()-t0:.0f}s]", flush=True)
    sv.delete()
    if hits:
        print(f"*** FORCED: {hits[:5]} -- spindle gives chi >= 6 ***",
              flush=True)
        break
print("DONE", flush=True)
