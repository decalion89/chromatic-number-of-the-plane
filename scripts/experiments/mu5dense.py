"""mu_5 on the dense 5-chromatic analogue.  The decisive asymmetry test.

At FOUR colours enrichment saturates mu:

    Sa                397 pts, deg  9.94   mu_4 = 3 for  42 of 600  (7 %)
    Sa + 1 orbit     2689 pts, deg 15.92   mu_4 = 3 for 150 of 150  (100 %)

Every rich point of the dense carrier is one colour short of blocking.  None
reaches 4, but the population moves all the way to k-1.

At FIVE colours nothing has ever exceeded 2, which is k-3 -- but the objects
measured there were sparse, 803 to 2041 points at degree 10 to 13.4.  The
comparison is only fair against the dense analogue, and the tuned chain built
exactly that: five_dense_2, 6925 points at mean degree 16.81, 5-chromatic, the
same regime as the 2689-point carrier that saturated at 3.

So: does enrichment move mu_5 from 2 to 3, as it moves mu_4 from 2 to 3?

  yes -> the mechanism scales, and six is a question of size
  no  -> the asymmetry is real, and five colours is not four colours with one
         more colour

The base colouring of this graph cost cadical 8921 s once already, and it is
the only expensive part; every candidate afterwards is assumption calls on the
warm solver.  So the run is one long wait and then a fast sweep.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time, json
from fractions import Fraction as Fr
from collections import Counter, defaultdict
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
from hn.field import Field
from hn.geometry import Point, Rotation, _rot60
from hn.graph import build_graph
from hn.blocked import MuSolver

ROOT = HN_DIR
t0 = time.time()
d = json.load(open(f"{ROOT}/data/five_dense_2.json"))
F = Field(tuple(d["field_generators"]))
P = [Point(F.element([Fr(a, b) for a, b in x]),
           F.element([Fr(a, b) for a, b in y])) for x, y in d["points"]]
g = build_graph(P); n = g.n; m = sum(len(a) for a in g.adj) // 2
print(f"  five_dense_2  n={n} m={m} deg={2.0*m/n:.2f}   [{time.time()-t0:.0f}s]",
      flush=True)
one = F.rational(Fr(1)); half = F.rational(Fr(1, 2))
rot60 = _rot60(F); r30 = Rotation(F.sqrt(3) * half, half)
S = set(g.vertices)
deg_order = sorted(range(n), key=lambda v: -len(g.adj[v]))
cand = set()
for c in deg_order[:25]:
    for rot in (rot60.about(g.vertices[c]), r30.about(g.vertices[c])):
        for p in g.vertices:
            z = rot(p)
            if z not in S:
                cand.add(z)
cells = defaultdict(list)
for i, p in enumerate(g.vertices):
    cells[(int(p.fx // 1), int(p.fy // 1))].append(i)
gx = [q.fx for q in g.vertices]; gy = [q.fy for q in g.vertices]
scored, seenn = [], set()
for z in cand:
    zx, zy = z.fx, z.fy
    cx, cy = int(zx // 1), int(zy // 1)
    nb = []
    for dx in (-1, 0, 1):
        for dy in (-1, 0, 1):
            for i in cells.get((cx + dx, cy + dy), ()):
                ex, ey = gx[i] - zx, gy[i] - zy
                if abs(ex * ex + ey * ey - 1.0) < 1e-9 and \
                   (g.vertices[i] - z).norm2() == one:
                    nb.append(i)
    if len(nb) >= 9:
        t = tuple(sorted(nb))
        if t not in seenn:
            seenn.add(t); scored.append((len(nb), t))
scored.sort(key=lambda u: -u[0])
print(f"    {len(scored)} neighbourhoods of size >= 9, largest "
      f"{scored[0][0] if scored else 0}   [{time.time()-t0:.0f}s]", flush=True)
hist = Counter()
t1 = time.time()
with MuSolver(g, k=5, budget=6_000_000) as ms:
    print(f"    base colouring in {time.time()-t1:.0f}s "
          f"(colourable={ms.colourable})", flush=True)
    if ms.colourable:
        # Only "is mu above 2?" matters, and that is ONE call -- forbid three
        # colours on N(p) and see whether a colouring survives.  The full
        # chain starts by asking whether N(p) can be monochromatic, which on a
        # seventeen-point neighbourhood of a 6925-vertex graph is both
        # expensive and beside the point.
        for sz, nb in scored[:60]:
            cheap = ms.at_most_two(list(nb))
            if cheap is True:
                hist[2] += 1
                if sum(hist.values()) % 5 == 0:
                    print("      ..%d: %s   [%.0fs]"
                          % (sum(hist.values()), dict(hist), time.time() - t0),
                          flush=True)
                continue
            if cheap is None:
                hist[-1] += 1
                continue
            v = ms.mu(list(nb))
            hist[v] += 1
            if v is not None and v >= 3:
                print(f"    *** mu_5 = {v} with |N|={sz} -- enrichment MOVED "
                      f"it ***", flush=True)
            if sum(hist.values()) % 10 == 0:
                print(f"      ..{sum(hist.values())}: {dict(hist)}"
                      f"   [{time.time()-t0:.0f}s]", flush=True)
print(f"    mu_5 over the richest {sum(hist.values())}: {dict(hist)}"
      f"   [{time.time()-t0:.0f}s]", flush=True)
