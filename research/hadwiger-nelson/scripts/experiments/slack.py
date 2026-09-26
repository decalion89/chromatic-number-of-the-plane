"""How much slack is there, measured instead of inferred.

The sharpest obstruction found today is that G has NO non-adjacent pair
forced to different colours at five -- so for every such pair, some proper
5-colouring gives them the same colour.  The colourings are free.  How free
has never been measured, and it needs no solver: take a colouring and count,
for each vertex, how many of the k colours its neighbours leave available.
A vertex with two or more free colours can be recoloured on the spot.

Run at the level where forcing is known to exist (Sa and Y at four) and at
the level where it does not (G at five), the two numbers say how much room
separates a problem that yields from one that does not.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time, pickle, random
from collections import defaultdict, Counter
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
from hn.degrey import build_G, build_Sa, build_Y
from hn.geometry import DEGREY_FIELD as K
from hn.fast import IntBasis, fast_edges_complete
from pysat.solvers import Solver

SC = ("/tmp/hn/")
t0 = time.time()
jobs = [("Sa", 4), ("Y", 4), ("tightmin.pkl", 4), ("Sa", 5), ("G", 5),
        ("deepgrow.pkl", 5), ("union3.pkl", 5)]
print(f"{'carrier':16s} {'k':>2s} {'points':>7s} {'deg':>6s} "
      f"{'free>=1':>8s} {'free>=2':>8s} {'mean free':>10s}", flush=True)
for name, KC in jobs:
    try:
        P = ({"G": lambda k: build_G(k, as_graph=False), "Sa": build_Sa,
              "Y": build_Y}[name](K) if name in ("G", "Sa", "Y")
             else pickle.load(open(SC + name, "rb")))
    except Exception:
        continue
    b = IntBasis.covering(P)
    rows = b.rows(P)
    if b.overflow_headroom(rows) >= 1.0:
        continue
    E = sorted(set((min(a, c), max(a, c))
                   for a, c in fast_edges_complete(b, rows)))
    n = len(P)
    if not E or n > 45000:
        continue
    adj = defaultdict(set)
    for a, c in E:
        adj[a].add(c)
        adj[c].add(a)
    cls = [[1 + v * KC + c for c in range(KC)] for v in range(n)]
    for a, c in E:
        for col in range(KC):
            cls.append([-(1 + a * KC + col), -(1 + c * KC + col)])
    s = Solver(name="cd15", bootstrap_with=cls)
    if not s.solve():
        print(f"{name:16s} {KC:2d} {n:7d}   not {KC}-colourable", flush=True)
        s.delete()
        continue
    m = s.get_model()
    s.delete()
    col = {}
    for v in range(n):
        for c in range(KC):
            if m[v * KC + c] > 0:
                col[v] = c
                break
    freecnt = Counter()
    tot = 0
    for v in range(n):
        used = {col[u] for u in adj[v]}
        free = KC - len(used) - (0 if col[v] in used else 1)
        free = max(0, KC - len(used | {col[v]}))
        freecnt[free] += 1
        tot += free
    ge1 = sum(c for f, c in freecnt.items() if f >= 1)
    ge2 = sum(c for f, c in freecnt.items() if f >= 2)
    print(f"{name:16s} {KC:2d} {n:7d} {2*len(E)/n:6.2f} "
          f"{100*ge1/n:7.1f}% {100*ge2/n:7.1f}% {tot/n:10.3f}"
          f"   [{time.time()-t0:.0f}s]", flush=True)
print("\n'free' counts colours a vertex could switch to without breaking the "
      "colouring; a graph with forcing has vertices at zero", flush=True)
print("DONE", flush=True)
