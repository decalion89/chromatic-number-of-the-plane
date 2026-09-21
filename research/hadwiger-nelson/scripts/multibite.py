"""Bite every unlocked ring at once, about the one centre that has them all.

sqrt(17) turned D = 5/3 into a doubly-usable ring of forty-eight points, four
times de Grey's twelve, for one generator and no extra vertices.  sqrt(13)
does the same for D = 1/3, thirty-six points.  And all three rings -- 1/3,
5/3, 4 -- sit about the SAME centre, G[0], which is where G's ring structure
actually lives.

There is no reason to bite one at a time.  Each rho_D contributes its ring's
worth of cross edges and they are different edges, so biting all three about
one centre compounds the coupling in a single union: 36 + 48 + 12 = 96 ring
points against de Grey's 12, and every one of them gains a unit edge it did
not have.

The unit ring D = 1 is left out on purpose.  It closes both radicals, but its
rotation is the 60-degree turn, which fixes the Eisenstein core -- it bites
nothing, which is why doubly_usable_ring excludes it.
"""
import sys, time
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
from fractions import Fraction as Fr
from hn.degrey import build_G
from hn.field import Field, embed
from hn.geometry import DEGREY_FIELD as K, Point, rotation_joining, Rotation
from hn.graph import build_graph
from hn.homcol import forced_same, doubly_usable_ring
from pysat.solvers import Solver

k = 5
RINGS = [Fr(5, 3), Fr(1, 3), Fr(4)]
RAD = (3, 5, 7, 11, 13, 17)
t0 = time.time()
KB = Field(RAD)
for D in RINGS:
    print(f"  D={D}: doubly usable over {RAD} -> "
          f"{doubly_usable_ring(D, RAD)}", flush=True)

G = [Point(embed(p.x, KB), embed(p.y, KB))
     for p in build_G(K, as_graph=False)]
C = G[0]
print(f"G: {len(G)} pts; centre ({C.fx:+.6f}, {C.fy:+.6f})"
      f"  [{time.time()-t0:.0f}s]", flush=True)


def ring(pts, D):
    return [i for i, p in enumerate(pts)
            if (lambda d: all(x == 0 for x in d.c[1:]) and d.c[0] == D)
            (p.dist2(C))]


turns = []
for D in RINGS:
    b = rotation_joining(D, KB)
    print(f"  D={D}: {len(ring(G, D))} points on the ring, bite "
          f"cos {b.cos} sin {b.sin}", flush=True)
    turns.append((D, b.about(C), Rotation(b.cos, -b.sin).about(C)))

U, seen = list(G), set(G)
front = {D: list(G) for D in RINGS}
for m in range(1, 7):
    for D, f, iv in turns:
        nxt = []
        for q in front[D]:
            for r in (f, iv):
                w = r(q)
                nxt.append(w)
                if w not in seen:
                    seen.add(w)
                    U.append(w)
        front[D] = nxt[:len(G)]
    g = build_graph(U)
    E = sorted(set((min(a, b), max(a, b)) for a, b in g.edges()))
    cls = [[1 + v * k + c for c in range(k)] for v in range(g.n)]
    for a, b in E:
        for c in range(k):
            cls.append([-(1 + a * k + c), -(1 + b * k + c)])
    sv = Solver(name="cd15", bootstrap_with=cls)
    ok = sv.solve()
    sizes = {str(D): len(ring(U, D)) for D in RINGS}
    print(f"m={m}: {g.n} pts, {len(E)} edges, rings {sizes}, "
          f"{'5-colourable' if ok else '*** NOT 5-COLOURABLE ***'}"
          f"  [{time.time()-t0:.0f}s]", flush=True)
    if not ok:
        print("*** chi >= 6 ***", flush=True)
        break
    idx = {q: i for i, q in enumerate(U)}
    tot = hit = 0
    for D in RINGS:
        anti = []
        for i in ring(U, D):
            q = U[i]
            a = Point(C.x + (C.x - q.x), C.y + (C.y - q.y))
            if a in idx and i < idx[a]:
                anti.append((i, idx[a]))
        h = [(i, j) for i, j in anti if forced_same(sv, i, j, k)]
        tot += len(anti)
        hit += len(h)
        if h:
            print(f"  *** FORCED at D={D}, spindle distance {4*D}: "
                  f"{h[:5]} ***", flush=True)
    print(f"    {tot} antipodal pairs over the three rings, {hit} FORCED"
          f"  [{time.time()-t0:.0f}s]", flush=True)
    sv.delete()
    if hit:
        break
print("DONE", flush=True)
