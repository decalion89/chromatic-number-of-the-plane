"""Thicken the ring, at the centre that actually has one.

The first attempt centred at the origin and the ring grew by one point per
step: G has a single vertex at radius 2 from the origin, so there was one
thread and nothing to thicken.  The centre scan had already said where the
ring is -- G[0], with twelve points on D = 4, the same twelve-point D6 orbit
Sa carried -- and G[0] is not the origin.

So centre there, and thread BOTH ways: rho and rho inverse.  Twelve threads
growing in two directions, each new point at distance exactly one from the
one it came from, and the D6 orbit is antipodally closed so the template's
pair is present from the first step.

One SAT call per step asks the only question that ends this outright: is the
union still 5-colourable.  Then the antipodal pairs, exactly.
"""
import sys, time
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
from fractions import Fraction as Fr
from hn.degrey import build_G
from hn.geometry import DEGREY_FIELD as F, Point, rotation_joining, Rotation
from hn.graph import build_graph
from pysat.solvers import Solver

k, D = 5, Fr(4)
t0 = time.time()
G = build_G(F, as_graph=False)
C = G[0]
print(f"G: {len(G)} pts; centre G[0] = ({C.fx:+.6f}, {C.fy:+.6f})",
      flush=True)
base = rotation_joining(D, F)
fwd = base.about(C)
inv = Rotation(base.cos, -base.sin).about(C)


def ringsize(pts):
    return sum(1 for p in pts
               if (lambda d: all(x == 0 for x in d.c[1:]) and d.c[0] == D)
               (p.dist2(C)))


print(f"  ring at D=4 about it: {ringsize(G)} points  "
      f"[{time.time()-t0:.0f}s]", flush=True)

U, seen = list(G), set(G)
hi, lo = list(G), list(G)
for m in range(1, 13):
    hi = [fwd(p) for p in hi]
    lo = [inv(p) for p in lo]
    added = 0
    for p in hi + lo:
        if p not in seen:
            seen.add(p)
            U.append(p)
            added += 1
    g = build_graph(U)
    E = sorted(set((min(a, b), max(a, b)) for a, b in g.edges()))
    cls = [[1 + v * k + c for c in range(k)] for v in range(g.n)]
    for a, b in E:
        for c in range(k):
            cls.append([-(1 + a * k + c), -(1 + b * k + c)])
    sv = Solver(name="cd15", bootstrap_with=cls)
    ok = sv.solve()
    print(f"m={m}: +{added} -> {g.n} pts, {len(E)} edges, "
          f"{ringsize(U)} on the ring, "
          f"{'5-colourable' if ok else '*** NOT 5-COLOURABLE ***'}"
          f"  [{time.time()-t0:.0f}s]", flush=True)
    if not ok:
        print("*** chi >= 6 ***", flush=True)
        break
    idx = {p: i for i, p in enumerate(U)}
    anti, hits = [], []
    for p in U:
        d = p.dist2(C)
        if all(x == 0 for x in d.c[1:]) and d.c[0] == D:
            q = Point(C.x + (C.x - p.x), C.y + (C.y - p.y))
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
