"""Colour the Erdos lattice graphs that defeat every coset colouring.

A patch of Z^2 scaled by 1/sqrt(N), with edges between points differing by a
vector of norm N -- so the edges are exactly unit distances, and the graph is
a genuine unit-distance graph in the plane, from a family with no ancestor in
common with the Moser spindle.

The twelve N found by the sublattice search are the ones where no sublattice
of index at most five avoids the connection set, so no coset colouring with
five colours exists.  That is necessary and not sufficient: the solver
decides whether some other 5-colouring does.
"""
import sys, time
from itertools import product
from pysat.solvers import Solver

N = int(sys.argv[1]) if len(sys.argv) > 1 else 46800
SIDE = int(sys.argv[2]) if len(sys.argv) > 2 else 40
KC = int(sys.argv[3]) if len(sys.argv) > 3 else 5
t0 = time.time()
conn = set()
a = 0
while a * a <= N:
    b2 = N - a * a
    b = int(b2 ** 0.5)
    for bb in (b - 1, b, b + 1):
        if bb >= 0 and a * a + bb * bb == N:
            for sa, sb in product((1, -1), repeat=2):
                conn.add((sa * a, sb * bb))
                conn.add((sb * bb, sa * a))
    a += 1
conn.discard((0, 0))
print(f"N = {N}: connection set has {len(conn)} vectors (degree {len(conn)})"
      f"  [{time.time()-t0:.0f}s]", flush=True)
pts = [(x, y) for x in range(SIDE) for y in range(SIDE)]
idx = {p: i for i, p in enumerate(pts)}
n = len(pts)
E = []
for (x, y) in pts:
    for (dx, dy) in conn:
        q = (x + dx, y + dy)
        if q in idx and idx[q] > idx[(x, y)]:
            E.append((idx[(x, y)], idx[q]))
print(f"patch {SIDE}x{SIDE} = {n} points, {len(E)} unit edges, mean degree "
      f"{2*len(E)/n:.2f}  [{time.time()-t0:.0f}s]", flush=True)
if not E:
    print("no edges -- the patch is smaller than the connection vectors",
          flush=True)
    sys.exit()
cls = [[1 + v * KC + c for c in range(KC)] for v in range(n)]
for a2, c2 in E:
    for col in range(KC):
        cls.append([-(1 + a2 * KC + col), -(1 + c2 * KC + col)])
s = Solver(name="cd15", bootstrap_with=cls)
ok = s.solve()
print(f"{KC}-colourable: {ok}  [{time.time()-t0:.0f}s]", flush=True)
if not ok:
    print(f"*** the patch is NOT {KC}-colourable ***", flush=True)
else:
    for k in range(KC - 1, 2, -1):
        cls2 = [[1 + v * k + c for c in range(k)] for v in range(n)]
        for a2, c2 in E:
            for col in range(k):
                cls2.append([-(1 + a2 * k + col), -(1 + c2 * k + col)])
        s2 = Solver(name="cd15", bootstrap_with=cls2)
        ok2 = s2.solve()
        s2.delete()
        print(f"   {k}-colourable: {ok2}  [{time.time()-t0:.0f}s]", flush=True)
        if not ok2:
            print(f"   chi = {k+1}", flush=True)
            break
print("DONE", flush=True)
