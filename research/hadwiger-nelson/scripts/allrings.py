"""All the rings at once, and how much slack Sa has to spare.

The gateway test asks one ring at a time: forbid that ring's antipodal pairs
and see whether the graph still colours.  Forbidding SEVERAL rings at once is
strictly stronger and costs the same single call, and it is still a usable
property -- "in every k-colouring some antipodal pair about this centre is
monochromatic" is what a multispindle over the rings consumes.

G about its hub has 228 antipodal pairs spread over eleven rings and not one
ring carries the property alone.  Whether all eleven together do is one call
and has never been asked.

The second question is calibration.  Sa's D=1 ring carries the property with
fifteen pairs and its D=4 ring with three.  How much of that is slack?  Drop
pairs one at a time and see how few still suffice -- the smallest surviving
subset is what the property actually costs, and it is the number to compare
against when asking what a five-colour analogue would need.
"""
import sys, time
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
from fractions import Fraction as Fr
from collections import defaultdict
from hn.degrey import build_G, build_Sa, build_Y
from hn.geometry import DEGREY_FIELD as F, Point, _rot60
from hn.graph import build_graph
from pysat.solvers import Solver

t0 = time.time()


def setup(P, C, k):
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
    base = [[1 + v * k + c for c in range(k)] for v in range(g.n)]
    for a, b in sorted(E):
        for c in range(k):
            base.append([-(1 + a * k + c), -(1 + b * k + c)])
    return g.n, E, by, base


def forbid(base, n, k, pairs):
    cls = list(base)
    for i, j in pairs:
        for c in range(k):
            cls.append([-(1 + i * k + c), -(1 + j * k + c)])
    sv = Solver(name="cd15", bootstrap_with=cls)
    ok = sv.solve()
    sv.delete()
    return ok


def allrings(name, P, C, k):
    n, E, by, base = setup(P, C, k)
    every = [pr for D in by for pr in by[D]]
    ok = forbid(base, n, k, every)
    print(f"{name} at {k}: {n} pts, {len(by)} rings, {len(every)} antipodal "
          f"pairs in all -> "
          f"{'still colourable' if ok else '*** UNCOLOURABLE ***'}"
          f"  [{time.time()-t0:.0f}s]", flush=True)
    return not ok


O = Point(F.zero(), F.zero())
G = build_G(F, as_graph=False)
allrings("G about G[0]", G, G[0], 5)
allrings("G about the origin", G, O, 5)

rot = _rot60(F).about(G[0])
seen, H = set(), []
for refl in (False, True):
    for j in range(6):
        for p in G:
            q = p
            for _ in range(j):
                q = rot(q)
            if refl:
                q = Point(q.x, G[0].y + (G[0].y - q.y))
            if q not in seen:
                seen.add(q)
                H.append(q)
allrings("G closed about its hub", H, G[0], 5)
allrings("Sa", build_Sa(F), O, 4)

# --- how few pairs does Sa actually need? ----------------------------------
print("\n=== how much slack Sa has ===", flush=True)
Sa = build_Sa(F)
n, E, by, base = setup(Sa, O, 4)
for D in (Fr(4), Fr(1)):
    pr = list(by[D])
    keep = list(pr)
    for x in list(pr):
        trial = [y for y in keep if y != x]
        if trial and not forbid(base, n, 4, trial):
            keep = trial
    print(f"  Sa D={D}: {len(pr)} antipodal pairs, and {len(keep)} of them "
          f"already suffice: {keep}  [{time.time()-t0:.0f}s]", flush=True)
print("\nDONE", flush=True)
