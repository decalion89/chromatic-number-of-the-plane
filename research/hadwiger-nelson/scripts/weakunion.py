"""The weak property on the unions, and on Sa's own ancestors.

G fails the gateway everywhere: 0 of 1664 (centre, ring) candidates carry
"some antipodal pair of this ring is monochromatic in every 5-colouring", in
fifty-five seconds.  That is why every strong test failed -- the strong
property implies the weak one, so a graph without the weak one cannot have it.

Two things to ask now.

The unions.  Sa has the weak property and Y sharpens it; maybe a union has it
where G does not.  Their construction centre is the one that matters, so scan
every ring about it rather than every centre -- instant, and it is the centre
the whole construction is organised around.

The ancestors.  What is the SMALLEST graph with the weak property at four?
Sa is 397 points and has it on two rings.  S is 39.  If the property appears
well below Sa, then what it needs is not size, and the calibration for five
colours changes.
"""
import sys, time
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
from fractions import Fraction as Fr
from collections import defaultdict
from hn.degrey import build_G, build_Y, build_Sa, build_S, build_Sb
from hn.field import Field, embed
from hn.geometry import DEGREY_FIELD as K, Point, rotation_joining, Rotation
from hn.graph import build_graph
from hn.homcol import closable_distance
from pysat.solvers import Solver

t0 = time.time()


def weak(name, P, C, k, quiet_ok=True):
    g = build_graph(P)
    E = set((min(a, b), max(a, b)) for a, b in g.edges())
    idx = {p: i for i, p in enumerate(P)}
    by = defaultdict(list)
    for i, p in enumerate(P):
        d = p.dist2(C)
        if not all(x == 0 for x in d.c[1:]) or d.c[0] == 0:
            continue
        q = Point(C.x + (C.x - p.x), C.y + (C.y - p.y))
        j = idx.get(q)
        if j is not None and i < j and (i, j) not in E:
            by[d.c[0]].append((i, j))
    cls = [[1 + v * k + c for c in range(k)] for v in range(g.n)]
    for a, b in sorted(E):
        for c in range(k):
            cls.append([-(1 + a * k + c), -(1 + b * k + c)])
    sel0 = g.n * k + 1
    order = sorted(by, key=lambda d: -len(by[d]))
    for t, D in enumerate(order):
        for i, j in by[D]:
            for c in range(k):
                cls.append([-(1 + i * k + c), -(1 + j * k + c), sel0 + t])
    sv = Solver(name="cd15", bootstrap_with=cls)
    if not sv.solve():
        print(f"{name}: NOT {k}-COLOURABLE ({g.n} pts)", flush=True)
        return
    hits = []
    for t, D in enumerate(order):
        if not sv.solve(assumptions=[-(sel0 + t)]):
            hits.append((D, len(by[D]), closable_distance(4 * D)))
    sv.delete()
    print(f"{name} at {k}: {g.n} pts, {len(E)} edges, {len(order)} rings "
          f"about the centre, {len(hits)} carry the weak property"
          f"  [{time.time()-t0:.0f}s]", flush=True)
    for D, m, spin in hits:
        print(f"   *** D={D}: {m} antipodal pairs at squared distance "
              f"{4*D}, spindle closable {spin} ***", flush=True)


O = Point(K.zero(), K.zero())
print("=== the ancestors, at four colours ===", flush=True)
weak("S", build_S(K), O, 4)
weak("Sa", build_Sa(K), O, 4)
weak("Sb", build_Sb(K), O, 4)
weak("Y", build_Y(K), O, 4)

print("\n=== G and its unions, at five ===", flush=True)
G = build_G(K, as_graph=False)
weak("G about G[0]", G, G[0], 5)
weak("G about the origin", G, O, 5)

K17 = Field((3, 5, 7, 11, 17))
GB = [Point(embed(p.x, K17), embed(p.y, K17)) for p in G]
C = GB[0]
for D, label in ((Fr(5, 3), "5/3"), (Fr(4), "4"), (Fr(1), "1")):
    base = rotation_joining(D, K17)
    f, iv = base.about(C), Rotation(base.cos, -base.sin).about(C)
    seen, U = set(), []
    for p in list(GB) + [f(p) for p in GB] + [iv(p) for p in GB]:
        if p not in seen:
            seen.add(p)
            U.append(p)
    weak(f"G u rho_{label}(G) u rho_{label}^-1(G) about G[0]", U, C, 5)
print("\nDONE", flush=True)
