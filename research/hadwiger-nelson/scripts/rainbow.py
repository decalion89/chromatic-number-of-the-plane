"""Look for the local structure the explanation says is missing: a 5-rainbow.

Everything measured here says the obstruction is the plane's, not the graph's.
A unit-distance graph has no clique bigger than a triangle, so against four
colours a triangle spends three and leaves one spare, and against five it
leaves two.  What would close that gap is a set of five points pairwise forced
onto different colours -- a 5-rainbow -- which cannot be made of edges and so
must use pairs forced apart at distances other than 1.

The search in hn/forced.py asks the solver about each pair, which is why it
has only ever been run small.  Sampling makes the question cheap: a vertex can
only complete a triangle if, in EVERY sampled colouring, it avoids all three
of the triangle's colours.  A random vertex manages that with probability
(2/5) per sample, so twenty samples leave 1e-8 of them by chance and every
survivor is a real candidate.  Vectorised, it is one numpy reduction per
triangle.

A triangle plus two such vertices, pairwise forced apart, is a 5-rainbow, and
a 5-rainbow in a 5-colourable graph is a gadget that constrains five colours
the way a triangle constrains four.
"""
import sys, time, pickle, random
from itertools import combinations
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
import numpy as np
from hn.degrey import build_G, build_Sa, build_Y
from hn.geometry import DEGREY_FIELD as F
from hn.graph import build_graph
from pysat.solvers import Solver

SC = ("/tmp/claude-0/-home-user-darwin-50/aceaa9ec-f432-5848-a506-"
      "39c59179b415/scratchpad/")
t0 = time.time()


def hunt(name, pts, k, samples=24, seed=1729, verify=True):
    g = build_graph(pts)
    n = g.n
    E = sorted((min(a, b), max(a, b)) for a, b in g.edges())
    cls = [[1 + v * k + c for c in range(k)] for v in range(n)]
    for a, b in E:
        for c in range(k):
            cls.append([-(1 + a * k + c), -(1 + b * k + c)])
    sv = Solver(name="g4", bootstrap_with=cls)
    if not sv.solve():
        print(f"  {name}: NOT {k}-colourable", flush=True)
        sv.delete()
        return
    rng, cols = random.Random(seed), []
    for s in range(samples):
        sv.set_phases([(1 if rng.random() < .5 else -1) * (1 + w)
                       for w in range(n * k)])
        sv.solve()
        m = sv.get_model()
        cols.append([next(c for c in range(k) if m[w * k + c] > 0)
                     for w in range(n)])
    C = np.array(cols, dtype=np.int8)
    tri = []
    for a, b in E:
        for c in g.adj[a] & g.adj[b]:
            if c > b:
                tri.append((a, b, c))
    print(f"  {name}: {n} points, {len(E)} edges, {len(tri)} triangles, "
          f"{samples} samples  [{time.time()-t0:.0f}s]", flush=True)

    best, quads = 0, []
    for a, b, c in tri:
        avoid = ((C != C[:, a:a + 1]) & (C != C[:, b:b + 1])
                 & (C != C[:, c:c + 1])).all(axis=0)
        idx = np.nonzero(avoid)[0]
        idx = [int(v) for v in idx if v not in (a, b, c)]
        if len(idx) > best:
            best = len(idx)
        if idx:
            quads.append((a, b, c, idx))
    print(f"    {len(quads)} triangles have a vertex avoiding all three "
          f"colours in every sample; best has {best}  "
          f"[{time.time()-t0:.0f}s]", flush=True)

    # Verification is cheap at five colours -- G's whole colouring costs
    # 6410 conflicts -- and ruinous at four, where a single forced-different
    # proof runs to minutes on these graphs.  So the four-colour controls
    # report the sampling only, which is the discovery step, and the
    # five-colour cases are verified in full.
    found = []
    if not verify:
        print(f"    (four colours: sampling only, verification skipped)",
              flush=True)
        sv.delete()
        return
    for a, b, c, idx in quads[:200]:
        good = []
        for v in idx:
            if all(not sv.solve(assumptions=[1 + v * k, 1 + u * k])
                   for u in (a, b, c)):
                good.append(v)
        if len(good) >= 2:
            for u, w in combinations(good, 2):
                if not sv.solve(assumptions=[1 + u * k, 1 + w * k]):
                    found.append((a, b, c, u, w))
                    print(f"    *** 5-RAINBOW: {(a, b, c, u, w)} ***",
                          flush=True)
                    break
        elif good:
            found.append((a, b, c, good[0]))
    sv.delete()
    four = [f for f in found if len(f) == 4]
    five = [f for f in found if len(f) == 5]
    print(f"    verified: {len(four)} 4-rainbows, {len(five)} 5-RAINBOWS  "
          f"[{time.time()-t0:.0f}s]", flush=True)
    if found:
        with open(SC + f"rainbow_{name.replace(' ', '_')}.pkl", "wb") as fh:
            pickle.dump((name, k, [(str(p.x), str(p.y)) for p in pts], found),
                        fh)


# The control first: at four colours a triangle leaves one spare, so a
# 4-rainbow is a triangle plus one vertex forced off all three -- the same
# question one level down, where the answer should be easy.
hunt("Sa at four", build_Sa(F), 4, verify=False)
hunt("Y at four", build_Y(F), 4, verify=False)
hunt("G at five", build_G(F, as_graph=False), 5)
hunt("Y at five", build_Y(F), 5)
from hn.geometry import Point, _rot60
r60 = _rot60(F)
PIV = Point(F.rational(-2), F.zero())
rot = r60.about(PIV)
G = build_G(F, as_graph=False)
Gs, seen = [], set()
for refl in (False, True):
    for j in range(6):
        for q0 in G:
            q = q0
            for _ in range(j):
                q = rot(q)
            if refl:
                q = Point(q.x, -q.y)
            if q not in seen:
                seen.add(q)
                Gs.append(q)
hunt("Gstar at five", Gs, 5)
