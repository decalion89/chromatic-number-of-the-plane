"""The ring sqrt(17) unlocks: D = 5/3, forty-eight points, free.

The template's bottleneck was never scale.  De Grey's field closes twenty-two
doubly-usable rings and G populates exactly one of them -- D = 4, twelve
points, Sa's own ring, inherited and never improved on -- so the bite at five
colours ran with the same twelve cross edges that sufficed at four.

Adjoining a radical costs no vertices.  The points of G do not move; only the
biting rotation needs sqrt(4D - 1) and the spindle sqrt(16D - 1).  Censusing
every rational ring at every vertex centre and asking which single adjunction
turns the biggest one usable:

    adjoin sqrt(17)  ->  D = 5/3 becomes doubly usable, 48 points
    adjoin sqrt(13)  ->  D = 1/3 becomes doubly usable, 36 points
    de Grey's field  ->  D = 4,                         12 points

and the arithmetic is clean: 4D - 1 = 17/3, squarefree part 51 = 3 * 17, so
the BITE is what wants the new radical; 16D - 1 = 77/3, squarefree part
231 = 3 * 7 * 11, already there, so the SPINDLE is free.  Four times de
Grey's coupling for one generator and not one extra point.

Centre is G[0], which is where the ring is.  Both the forward and inverse
bites, because threading both ways is what made the ring grow before.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
from fractions import Fraction as Fr
from hn.degrey import build_G
from hn.field import Field, embed
from hn.geometry import DEGREY_FIELD as K, Point, rotation_joining, Rotation
from hn.graph import build_graph
from hn.homcol import forced_same, closable_distance, doubly_usable_ring
from pysat.solvers import Solver

k, D = 5, Fr(5, 3)
t0 = time.time()
K17 = Field((3, 5, 7, 11, 17))
print(f"D = {D}: closable {closable_distance(D, (3,5,7,11,17))}, "
      f"doubly usable {doubly_usable_ring(D, (3,5,7,11,17))}; "
      f"in de Grey's own field doubly usable "
      f"{doubly_usable_ring(D, (3,5,7,11))}", flush=True)

G0 = build_G(K, as_graph=False)
G = [Point(embed(p.x, K17), embed(p.y, K17)) for p in G0]
C = G[0]
print(f"G re-embedded: {len(G)} pts; centre ({C.fx:+.6f}, {C.fy:+.6f})"
      f"  [{time.time()-t0:.0f}s]", flush=True)


def ring(pts):
    out = []
    for i, p in enumerate(pts):
        d = p.dist2(C)
        if all(x == 0 for x in d.c[1:]) and d.c[0] == D:
            out.append(i)
    return out


print(f"  ring at D={D} about it: {len(ring(G))} points"
      f"  [{time.time()-t0:.0f}s]", flush=True)
base = rotation_joining(D, K17)
print(f"  bite: cos {base.cos}, sin {base.sin}", flush=True)
fwd = base.about(C)
inv = Rotation(base.cos, -base.sin).about(C)
p = G[ring(G)[0]]
print(f"  check: a ring point and its image are {p.dist2(fwd(p))} apart",
      flush=True)

U, seen = list(G), set(G)
hi, lo = list(G), list(G)
for m in range(1, 9):
    hi = [fwd(q) for q in hi]
    lo = [inv(q) for q in lo]
    for q in hi + lo:
        if q not in seen:
            seen.add(q)
            U.append(q)
    g = build_graph(U)
    E = sorted(set((min(a, b), max(a, b)) for a, b in g.edges()))
    cls = [[1 + v * k + c for c in range(k)] for v in range(g.n)]
    for a, b in E:
        for c in range(k):
            cls.append([-(1 + a * k + c), -(1 + b * k + c)])
    sv = Solver(name="cd15", bootstrap_with=cls)
    ok = sv.solve()
    R = ring(U)
    print(f"m={m}: {g.n} pts, {len(E)} edges, {len(R)} on the ring, "
          f"{'5-colourable' if ok else '*** NOT 5-COLOURABLE ***'}"
          f"  [{time.time()-t0:.0f}s]", flush=True)
    if not ok:
        print("*** chi >= 6 ***", flush=True)
        break
    idx = {q: i for i, q in enumerate(U)}
    anti, hits = [], []
    for i in R:
        q = U[i]
        a = Point(C.x + (C.x - q.x), C.y + (C.y - q.y))
        if a in idx and i < idx[a]:
            anti.append((i, idx[a]))
    for i, j in anti:
        if forced_same(sv, i, j, k):
            hits.append((i, j))
    print(f"    {len(anti)} antipodal pairs at squared distance {4*D}, "
          f"{len(hits)} FORCED  [{time.time()-t0:.0f}s]", flush=True)
    sv.delete()
    if hits:
        print(f"*** FORCED: {hits[:5]} -- spindle gives chi >= 6 ***",
              flush=True)
        break
print("DONE", flush=True)
