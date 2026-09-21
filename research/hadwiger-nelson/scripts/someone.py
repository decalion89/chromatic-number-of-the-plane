"""Not "this pair is forced" but "SOME pair of the ring is" -- one SAT call.

The hypothesis tested so far has been the strong one: a named pair, the same
colour in every proper colouring.  It is what the spindling lemma consumes,
and every union built here fails it -- 336 antipodal pairs across three rings
and not one forced.

The weaker hypothesis is the interesting one and it has never been measured.
"In every proper k-colouring, at least one antipodal pair of this ring is
monochromatic" is exactly "the graph with ALL those antipodal edges added is
not k-colourable", because an added edge forbids precisely that pair from
being monochromatic.  One call for the whole ring, not one per pair.

It is strictly weaker than forcing a named pair, so it can hold where the
strong one fails, and it is the natural quantity: it says the ring as a whole
is constrained even when no single pair is.  A ring that satisfies it is what
a multispindle consumes -- m copies rather than two -- so it is worth knowing
which rings are close.

Measured on Sa and Y at four, where the mechanism is known to work, and on G
and the unions at five, where it has not.
"""
import sys, time, pickle, os
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
from fractions import Fraction as Fr
from collections import defaultdict
from hn.degrey import build_G, build_Y, build_Sa
from hn.geometry import DEGREY_FIELD as F, Point, rotation_joining, Rotation
from hn.graph import build_graph
from hn.homcol import closable_distance
from pysat.solvers import Solver

t0 = time.time()


def antipodal_by_ring(P, C):
    """Antipodal pairs about C, grouped by the ring they sit on."""
    idx = {p: i for i, p in enumerate(P)}
    out = defaultdict(list)
    for p in P:
        d = p.dist2(C)
        if not all(x == 0 for x in d.c[1:]):
            continue
        q = Point(C.x + (C.x - p.x), C.y + (C.y - p.y))
        if q in idx and idx[p] < idx[q]:
            out[d.c[0]].append((idx[p], idx[q]))
    return out


def someone(name, P, C, k):
    g = build_graph(P)
    E = set((min(a, b), max(a, b)) for a, b in g.edges())
    base = [[1 + v * k + c for c in range(k)] for v in range(g.n)]
    for a, b in sorted(E):
        for c in range(k):
            base.append([-(1 + a * k + c), -(1 + b * k + c)])
    sv0 = Solver(name="cd15", bootstrap_with=base)
    ok = sv0.solve()
    print(f"\n{name} at {k}: {g.n} pts, {len(E)} edges, "
          f"{'colourable' if ok else 'NOT COLOURABLE'}"
          f"  [{time.time()-t0:.0f}s]", flush=True)
    sv0.delete()
    if not ok:
        return
    rings = antipodal_by_ring(P, C)
    for D in sorted(rings, key=lambda d: -len(rings[d])):
        pr = [(i, j) for i, j in rings[D] if (i, j) not in E]
        if not pr:
            continue
        cls = list(base)
        for i, j in pr:
            for c in range(k):
                cls.append([-(1 + i * k + c), -(1 + j * k + c)])
        sv = Solver(name="cd15", bootstrap_with=cls)
        good = sv.solve()
        sv.delete()
        spin = closable_distance(4 * D)
        verdict = ("still colourable" if good else
                   "*** UNCOLOURABLE: SOME PAIR IS ALWAYS MONOCHROMATIC ***")
        print(f"   ring D={D}: {len(pr)} antipodal pairs at squared "
              f"distance {4*D}, all forbidden -> {verdict}"
              f"  (spindle closable: {spin})  [{time.time()-t0:.0f}s]",
              flush=True)


O = Point(F.zero(), F.zero())
someone("Sa", build_Sa(F), O, 4)
someone("Y", build_Y(F), O, 4)
G = build_G(F, as_graph=False)
someone("G", G, G[0], 5)
someone("G", G, O, 5)

# the union that has the most antipodal pairs: the sqrt(17) ring, one bite
from hn.field import Field, embed
K17 = Field((3, 5, 7, 11, 17))
GB = [Point(embed(p.x, K17), embed(p.y, K17)) for p in G]
C = GB[0]
rot = rotation_joining(Fr(5, 3), K17).about(C)
seen, U = set(), []
for p in list(GB) + [rot(p) for p in GB]:
    if p not in seen:
        seen.add(p)
        U.append(p)
someone("G u rho_{5/3}(G)", U, C, 5)
print("\nDONE", flush=True)
