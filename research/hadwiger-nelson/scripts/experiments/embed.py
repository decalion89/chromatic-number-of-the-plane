"""De Grey's gadget EMBEDDED, exhausted in positions.

The isolated exhaustion said the bite imposes a floor on the arc span, not a
ceiling -- the opposite of what a spindle needs.  But de Grey's lemma is not
about the gadget alone.  It holds for the gadget sitting inside Sa, and the
ambient graph is what does the work: seven points on their own cannot be
capped at two colours either, since a hexagon plus a centre joined to it has
perfectly good 3-colourings.

So put the question to the carrier.  Find the origin and the six points at
distance 2 from it -- (+-2,0) and (+-1,+-sqrt 3), de Grey's own ring -- and
ask the solver, over ALL homomorphisms of the whole graph to K(9/2), how
wide an arc those seven points span and how close the antipodal pairs are
forced.  That is the exhaustion he did in colours, done in positions, with
the carrier included this time.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time, pickle
from fractions import Fraction as Fr
from collections import defaultdict, Counter
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
from hn.degrey import build_G, build_Sa, build_Y
from hn.geometry import DEGREY_FIELD as K
from hn.fast import IntBasis, fast_edges_complete
from pysat.solvers import Solver

SC = ("/tmp/hn/")
CAR = sys.argv[1] if len(sys.argv) > 1 else "G"
R = Fr(*map(int, (sys.argv[2] if len(sys.argv) > 2 else "9/2").split("/")))
t0 = time.time()
P = ({"G": lambda k: build_G(k, as_graph=False), "Sa": build_Sa,
      "Y": build_Y}[CAR](K) if CAR in ("G", "Sa", "Y")
     else pickle.load(open(SC + CAR, "rb")))
b = IntBasis.covering(P)
rows = b.rows(P)
E = sorted(set((min(a, c), max(a, c)) for a, c in fast_edges_complete(b, rows)))
n, p, q = len(P), R.numerator, R.denominator
xy = [(float(t.x), float(t.y)) for t in P]
# the anchor is build_G's first point; the ring is what sits at distance 2
cx, cy = xy[0]
ring = []
for v in range(n):
    d2 = (xy[v][0] - cx) ** 2 + (xy[v][1] - cy) ** 2
    if abs(d2 - 4.0) < 1e-9:
        ring.append(v)
print(f"{CAR}: {n} pts, {len(E)} edges; centre v0 at ({cx:.3f},{cy:.3f}), "
      f"{len(ring)} points at distance 2  [{time.time()-t0:.0f}s]", flush=True)
import math
ring.sort(key=lambda v: math.atan2(xy[v][1] - cy, xy[v][0] - cx))
S = [0] + ring
cls = [[1 + v * p + j for j in range(p)] for v in range(n)]
for a, c in E:
    for j in range(p):
        for d in range(-(q - 1), q):
            cls.append([-(1 + a * p + j), -(1 + c * p + (j + d) % p)])
cls += [[-(1 + j)] for j in range(1, p)] + [[1]]   # pin the centre at 0
s = Solver(name="cd15", bootstrap_with=cls)
if not s.solve():
    print("carrier refuses the ratio outright", flush=True)
    sys.exit()
print(f"centre pinned at position 0; sampling ring patterns"
      f"  [{time.time()-t0:.0f}s]", flush=True)
pats = Counter()
spans = Counter()
seen = set()
for it in range(4000):
    if not s.solve():
        break
    m = s.get_model()
    pos = {}
    for v in S:
        for j in range(p):
            if m[v * p + j] > 0:
                pos[v] = j
                break
    key = tuple(pos[v] for v in S)
    if key in seen:
        break
    seen.add(key)
    used = sorted(set(key))
    gap = max((used[(i + 1) % len(used)] - used[i]) % p
              for i in range(len(used))) if len(used) > 1 else p
    spans[p - gap] += 1
    pats[len(used)] += 1
    s.add_clause([-(1 + v * p + pos[v]) for v in S])
print(f"\n{len(seen)} distinct patterns of centre+ring found"
      f"  [{time.time()-t0:.0f}s]", flush=True)
print(f"   arc span: {dict(sorted(spans.items()))}", flush=True)
print(f"   distinct positions used: {dict(sorted(pats.items()))}", flush=True)
print(f"   narrowest arc seen: {min(spans) if spans else '-'} of {p}; "
      f"a cap below {2*q} would force independence", flush=True)
print("DONE", flush=True)
