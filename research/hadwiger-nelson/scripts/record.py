"""Push a mu_4 = 3 point to 4, under a budget of 110 extra vertices.

Sa has candidate points p with thirteen graph vertices on their unit circle and
mu_4(Sa, p) = 3: in every 4-colouring of Sa, N(p) already carries three
colours.  One more and p is BLOCKED -- Sa + p would have no 4-colouring at all,
a 5-chromatic unit-distance graph on 398 vertices.

The record for a 5-chromatic unit-distance graph is the 509-vertex Parts graph,
so there are 110 vertices of headroom: any A with |A| <= 110 for which

        mu_4(Sa u A, p) = 4

gives Sa u A u {p} at 509 or fewer.  That is a far smaller target than six, it
is the same instrument, and nothing about it needs a new idea -- only the right
A.

Which points can help?  Adding a vertex can only constrain N(p) if it sees it,
so the pool is the in-field points at distance exactly 1 from at least two
members of N(p): each such point forbids its own colour to two of them at once.
Points touching only one member add a neighbour but no relation between the
neighbours, and relations between the neighbours are the whole question, since
N(p) is bipartite on its own and will take two colours unless something outside
stops it.

Greedy, in batches, re-measuring mu after each: the measurement is a few
assumption calls on a warm solver and the graph is four hundred points, so the
loop is cheap and the budget, not the clock, is the binding constraint.
"""
import sys, time, json
from fractions import Fraction as Fr
from collections import defaultdict
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
from hn.field import Field
from hn.degrey import build_Sa
from hn.geometry import Point, Rotation, _rot60
from hn.graph import build_graph
from hn.blocked import MuSolver, neighbourhood

ROOT = "/home/user/darwin-50/research/hadwiger-nelson"
t0 = time.time()
F = Field((3, 5, 7, 11))
one = F.rational(Fr(1))
half = F.rational(Fr(1, 2))
rot60 = _rot60(F)
r30 = Rotation(F.sqrt(3) * half, half)
Sa = build_Sa(F)
g = build_graph(Sa); n = g.n
S = set(g.vertices)
print(f"Sa n={n} m={sum(len(a) for a in g.adj)//2}   [{time.time()-t0:.0f}s]",
      flush=True)

# in-field candidate points, scored by how much of Sa sits on their unit circle
deg_order = sorted(range(n), key=lambda v: -len(g.adj[v]))
cand = set()
for c in deg_order:
    for rot in (rot60.about(g.vertices[c]), r30.about(g.vertices[c])):
        for p in g.vertices:
            z = rot(p)
            if z not in S:
                cand.add(z)
print(f"  {len(cand)} candidates   [{time.time()-t0:.0f}s]", flush=True)
cells = defaultdict(list)
for i, p in enumerate(g.vertices):
    cells[(int(p.fx // 1), int(p.fy // 1))].append(i)
gx = [q.fx for q in g.vertices]; gy = [q.fy for q in g.vertices]
def nbr(z, verts, cls, fx, fy):
    zx, zy = z.fx, z.fy
    cx, cy = int(zx // 1), int(zy // 1)
    out = []
    for dx in (-1, 0, 1):
        for dy in (-1, 0, 1):
            for i in cls.get((cx + dx, cy + dy), ()):
                ex, ey = fx[i] - zx, fy[i] - zy
                if abs(ex * ex + ey * ey - 1.0) < 1e-9 and \
                   (verts[i] - z).norm2() == one:
                    out.append(i)
    return out
scored = []
for z in cand:
    nb = nbr(z, g.vertices, cells, gx, gy)
    if len(nb) >= 11:
        scored.append((len(nb), z, tuple(sorted(nb))))
scored.sort(key=lambda u: -u[0])
print(f"  {len(scored)} candidates with >= 11 neighbours, largest "
      f"{scored[0][0] if scored else 0}   [{time.time()-t0:.0f}s]", flush=True)

target = None
with MuSolver(g, k=4, budget=6_000_000) as ms:
    for k, z, nb in scored[:20]:
        v = ms.mu(list(nb))
        print(f"    |N|={k}: mu_4 = {v}   [{time.time()-t0:.0f}s]", flush=True)
        if v == 4:
            print("    *** already blocked -- Sa + p is 5-chromatic on 398 ***",
                  flush=True)
            target = (z, nb); break
        if v == 3 and target is None:
            target = (z, nb)
if target is None:
    print("  no mu_4 = 3 candidate found", flush=True); sys.exit()
P, NB = target
print(f"  target p with |N|={len(NB)}, mu_4 = 3   [{time.time()-t0:.0f}s]",
      flush=True)

# The pool has to be built FROM N(p), and without needing a square root the
# field does not have.
#
# The analytic route -- intersect the unit circles of a pair u, v -- needs
# sqrt(1/d2 - 1/4), and 59 of the 78 pairs of a thirteen-point neighbourhood
# have an IRRATIONAL d2, so that root is not a rational one and this field
# offers no general square root on an element.  (The pairs whose d2 is
# rational need only sqrt3 and sqrt11, which are already here; they give nine
# points, and nine is not enough to move anything.)
#
# The synthetic route avoids it completely.  Every point at distance 1 from a
# member u of N(p) is u + e for a UNIT VECTOR e of the field, and unit vectors
# are free: each spindle letter (cos, sin) = (1 - 1/(2r), sqrt(4r-1)/(2r)) is
# one, so is every power of the 60 and 30 degree rotations, and so is every
# product of them.  Sweeping u over N(p) and e over a few hundred unit vectors
# gives an exactly-represented pool, and the ones worth adding are those that
# land a unit from TWO or more members at once -- an added vertex can only
# speak about the RELATION between two neighbours if it sees both, and
# relations between them are the whole question, N(p) being bipartite alone.
NBset = set(NB)
NBpts = [g.vertices[i] for i in NB]
units = []
seenu = set()
base = [rot60, r30]
for d2n, d2d in ((1, 1), (3, 1), (4, 1), (16, 1), (64, 9), (64, 3), (256, 9),
                 (1, 3), (4, 3), (9, 1), (25, 1), (49, 1), (2, 1), (12, 1)):
    try:
        base.append(rotation_joining(Fr(d2n, d2d), F))
    except Exception:
        pass
cur_rot = Rotation(one, F.zero())
for r1 in base:
    acc = Rotation(one, F.zero())
    for _ in range(12):
        acc = Rotation(acc.cos * r1.cos - acc.sin * r1.sin,
                       acc.cos * r1.sin + acc.sin * r1.cos)
        key = (acc.cos, acc.sin)
        if key not in seenu:
            seenu.add(key)
            units.append(Point(acc.cos, acc.sin))
            units.append(Point(acc.cos, -acc.sin))
print(f"  {len(units)} unit vectors in the field   [{time.time()-t0:.0f}s]",
      flush=True)
pool = {}
for u in NBpts:
    for e in units:
        z = Point(u.x + e.x, u.y + e.y)
        if z in S or z == P:
            continue
        if z in pool:
            continue
        hits = [i for i in nbr(z, g.vertices, cells, gx, gy) if i in NBset]
        if len(hits) >= 2:
            pool[z] = len(hits)
print(f"  pool of {len(pool)} points seeing two or more of N(p)"
      f"   [{time.time()-t0:.0f}s]", flush=True)
order = sorted(pool, key=lambda z: -pool[z])

BUDGET = 110
added, cur = [], list(g.vertices)
for batch in (20, 20, 20, 20, 15, 15):
    if len(added) + batch > BUDGET:
        batch = BUDGET - len(added)
    if batch <= 0:
        break
    take = [z for z in order if z not in set(added)][:batch]
    if not take:
        break
    added += take
    gg = build_graph(cur + added)
    nb2 = neighbourhood(gg, P)
    with MuSolver(gg, k=4, budget=6_000_000) as ms:
        if not ms.colourable:
            print(f"  *** Sa + {len(added)} points already refuses four "
                  f"({gg.n} vertices) ***", flush=True)
            json.dump({"field_generators": list(F.gens), "n": gg.n,
                       "points": [[[[c.numerator, c.denominator] for c in q.x.c],
                                   [[c.numerator, c.denominator] for c in q.y.c]]
                                  for q in gg.vertices]},
                      open(f"{ROOT}/data/five_small_record.json", "w"))
            break
        v = ms.mu(nb2)
    print(f"  +{len(added)} points (n={gg.n}), |N(p)|={len(nb2)}: mu_4 = {v}"
          f"   [{time.time()-t0:.0f}s]", flush=True)
    if v is not None and v >= 4:
        total = gg.n + 1
        print(f"  *** BLOCKED: Sa + {len(added)} points + p is 5-chromatic on "
              f"{total} vertices (record is 509) ***", flush=True)
        json.dump({"field_generators": list(F.gens), "n": total,
                   "points": [[[[c.numerator, c.denominator] for c in q.x.c],
                               [[c.numerator, c.denominator] for c in q.y.c]]
                              for q in list(gg.vertices) + [P]]},
                  open(f"{ROOT}/data/five_small_record.json", "w"))
        break
