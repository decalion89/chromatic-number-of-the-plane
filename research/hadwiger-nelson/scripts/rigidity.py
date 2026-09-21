"""Is the base's colouring space rigid?  That is the quantity, not the ring.

Y forced with SIX antipodal pairs.  The unions here now carry 192 and force
nothing, so the count of antipodal pairs is not what the mechanism turns on.

What Sa has and G does not may be rigidity.  G's colour relation at five
colours was already shown to be exactly its edge set: no two non-adjacent
vertices are forced together, and none are forced apart either -- its
5-colourings are generic.  If Sa's 4-colourings are NOT generic, if
non-adjacent pairs are forced to differ, then the bite is acting on a rigid
space and that is why twelve cross edges sufficed.

Both directions are one exact SAT call, by the same colour symmetry:
  forced SAME      c(i)=0, c(j)!=0  UNSAT
  forced DIFFERENT c(i)=0, c(j)=0   UNSAT
so the whole relation of a 397-point graph is affordable outright.

Also here, because it is one union and cheap: the D=1 bite about G[0].  Its
ring is 60 points and the bite fixes it setwise -- none of the sixty moves --
so it adds no cross edge on the ring, which is why doubly_usable_ring rules
D=1 out.  But it does add 792 points elsewhere, and the ring's sixty antipodal
pairs at squared distance 4 are already present.  Free to ask whether the new
points force any of them.
"""
import sys, time
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
from fractions import Fraction as Fr
from hn.degrey import build_G, build_Y, build_Sa
from hn.geometry import DEGREY_FIELD as F, Point, rotation_joining
from hn.graph import build_graph
from hn.homcol import forced_same
from pysat.solvers import Solver

t0 = time.time()


def solver_for(P, k):
    g = build_graph(P)
    E = set((min(a, b), max(a, b)) for a, b in g.edges())
    cls = [[1 + v * k + c for c in range(k)] for v in range(g.n)]
    for a, b in sorted(E):
        for c in range(k):
            cls.append([-(1 + a * k + c), -(1 + b * k + c)])
    sv = Solver(name="cd15", bootstrap_with=cls)
    return sv, E, g.n


def relation(name, P, k, cap=None):
    sv, E, n = solver_for(P, k)
    if not sv.solve():
        print(f"{name}: not {k}-colourable", flush=True)
        return
    same = diff = tested = 0
    for i in range(n - 1):
        for j in range(i + 1, n):
            if (i, j) in E:
                continue
            tested += 1
            if not sv.solve(assumptions=[1 + i * k, 1 + j * k]):
                diff += 1
                if diff <= 6:
                    print(f"    forced DIFFERENT: {i},{j}", flush=True)
            elif forced_same(sv, i, j, k):
                same += 1
                if same <= 6:
                    print(f"    forced SAME: {i},{j}", flush=True)
        if cap and time.time() - t0 > cap:
            print(f"  ({name}: stopped at row {i})", flush=True)
            break
    print(f"{name} at {k} colours: {n} pts, {len(E)} edges, {tested} "
          f"non-edge pairs tested, {diff} forced DIFFERENT, {same} forced "
          f"SAME  [{time.time()-t0:.0f}s]", flush=True)
    sv.delete()


relation("Sa", build_Sa(F), 4)
relation("Y", build_Y(F), 4, cap=1500)
print("\n=== the D=1 bite about G[0] ===", flush=True)
G = build_G(F, as_graph=False)
C = G[0]
rot = rotation_joining(Fr(1), F).about(C)
seen, U = set(), []
for p in list(G) + [rot(p) for p in G]:
    if p not in seen:
        seen.add(p)
        U.append(p)
sv, E, n = solver_for(U, 5)
print(f"union: {n} pts, {len(E)} edges, "
      f"{'5-colourable' if sv.solve() else 'NOT 5-COLOURABLE'}"
      f"  [{time.time()-t0:.0f}s]", flush=True)
idx = {p: i for i, p in enumerate(U)}
anti = []
for p in U:
    if p.dist2(C) == 1:
        q = Point(C.x + (C.x - p.x), C.y + (C.y - p.y))
        if q in idx and idx[p] < idx[q]:
            anti.append((idx[p], idx[q]))
hits = [(i, j) for i, j in anti if forced_same(sv, i, j, 5)]
print(f"  {len(anti)} antipodal pairs at squared distance 4, {len(hits)} "
      f"FORCED  [{time.time()-t0:.0f}s]", flush=True)
print("DONE", flush=True)
