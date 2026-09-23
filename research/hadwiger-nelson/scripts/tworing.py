"""The universal two-ring configuration about a point, and its chromatic number.

Fix h and keep only the points at distance 1 and at distance 1/sqrt3 from it.
All the unit-distance edges among them are decided by angles alone:

    radius 1     to radius 1      a unit apart at   60 degrees
    radius 1/r3  to radius 1/r3   a unit apart at  120 degrees
    radius 1     to radius 1/r3   a unit apart at  theta0 = arccos(sqrt3/6)

So every such configuration is an induced subgraph of ONE infinite graph, and
that graph has coordinates.  Every reachable angle is 60i + theta0*j, the ring
being radius 1 for j even and 1/sqrt3 for j odd, and the two coordinates never
collide: 2*cos(theta0) = 1/sqrt3 is not an algebraic integer (its minimal
polynomial 3x^2 - 1 is not monic over Z), so by Niven-Lehmer theta0 is not a
rational multiple of pi.  Since 60*6 = 360, i lives in Z/6 and j in Z:

    level j even   a hexagon      i ~ i +- 1 (mod 6)
    level j odd    two triangles  i ~ i +- 2 (mod 6)
    between levels a perfect matching  (i,j) ~ (i,j+-1)

which is the whole content of "the unit circle is bipartite" and "the 1/sqrt3
ring is a union of triangles", now stacked with the cross edges that join them.
Its chromatic number is the ceiling on everything a two-ring argument can ever
reach, in any unit-distance graph whatsoever -- so it is worth knowing exactly,
and it is a small computation.
"""
import sys, time
from pysat.solvers import Solver

t0 = time.time()
def build(levels):
    V = [(i, j) for j in range(levels) for i in range(6)]
    idx = {v: k for k, v in enumerate(V)}
    E = set()
    for j in range(levels):
        step = 1 if j % 2 == 0 else 2
        for i in range(6):
            for s in (step, -step):
                E.add(tuple(sorted((idx[(i, j)], idx[((i + s) % 6, j)]))))
        if j + 1 < levels:
            for i in range(6):
                E.add(tuple(sorted((idx[(i, j)], idx[(i, j + 1)]))))
    return len(V), sorted(E)

def chi(n, E, cap=6):
    for k in range(1, cap + 1):
        X = lambda v, c: 1 + v * k + c
        cnf = [[X(v, c) for c in range(k)] for v in range(n)]
        for a, b in E:
            for c in range(k):
                cnf.append([-X(a, c), -X(b, c)])
        s = Solver(name="cd19", bootstrap_with=cnf)
        ok = s.solve()
        mod = s.get_model() if ok else None
        s.delete()
        if ok:
            col = [next(c for c in range(k) if mod[X(v, c) - 1] > 0) for v in range(n)]
            return k, col
    return cap + 1, None

for L in (2, 3, 4, 5, 6, 8, 10, 14, 20, 30):
    n, E = build(L)
    k, col = chi(n, E)
    print(f"  levels={L:<3d} n={n:<4d} edges={len(E):<5d}  chi = {k}   "
          f"[{time.time()-t0:.0f}s]", flush=True)
    if k >= 4:
        print("   *** the two-ring configuration passes three ***", flush=True)
        break
n, E = build(6)
k, col = chi(n, E)
print(f"\n  a 3-colouring of six levels, by level:", flush=True)
for j in range(6):
    print(f"    j={j} ({'hexagon' if j%2==0 else 'triangles'}): "
          f"{[col[j*6+i] for i in range(6)]}", flush=True)
