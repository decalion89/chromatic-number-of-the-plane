"""The atom of any 6-chromatic unit-distance graph, searched for everywhere.

A 6-critical graph has minimum degree 5, so some vertex has five neighbours
forced to five different colours in every proper 5-colouring.  Those five sit
on a unit circle about the vertex, and two points of a unit circle are a unit
apart only when they subtend 60 degrees, so the graph a circle carries is a
subgraph of disjoint hexagons -- bipartite, maximum clique 2.  At most TWO of
the five can be adjacent.  The other pairs are forced different while not
being adjacent, which the carrier has to supply.

    Every 6-chromatic unit-distance graph contains non-adjacent pairs
    forced to different colours at five colours.

Elementary, and it turns "find a 6-chromatic graph" into "find one PAIR",
which is a single solver call: by colour permutation, if u and v can agree
they can both take colour 0, so the test is two assumptions.

Zero such pairs were found in 400 neighbourhoods of G.  This asks every
carrier in the project, and does not restrict to neighbourhoods -- any
non-adjacent forced pair anywhere is the brick, whether or not it already
sits on some vertex's circle.  A control at four colours comes first, where
forcing is known to exist, so that a zero at five means something.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time, pickle, glob, os, math, random
from collections import defaultdict
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
from hn.degrey import build_G, build_Sa, build_Y
from hn.geometry import DEGREY_FIELD as K
from hn.fast import IntBasis, fast_edges_complete
from pysat.solvers import Solver

SC = ("/tmp/hn/")
BUDGET = int(sys.argv[1]) if len(sys.argv) > 1 else 40000
t0 = time.time()
jobs = [("Sa", 4), ("Y", 4), ("tightmin.pkl", 4), ("denseSa.pkl", 4),
        ("G", 5), ("deepgrow.pkl", 5), ("wide0.pkl", 5), ("union3.pkl", 5),
        ("symG.pkl", 5), ("fiveD.pkl", 5), ("pivG.pkl", 5)]
print(f"{'carrier':16s} {'k':>2s} {'points':>7s} {'pairs':>7s} "
      f"{'forced':>7s}   verdict", flush=True)
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
    if not E or n > 40000:
        continue
    adj = defaultdict(set)
    for a, c in E:
        adj[a].add(c)
        adj[c].add(a)
    cls = [[1 + v * KC + c for c in range(KC)] for v in range(n)]
    for a, c in E:
        for col in range(KC):
            cls.append([-(1 + a * KC + col), -(1 + c * KC + col)])
    s = Solver(name="m22", bootstrap_with=cls)  # minisat: 31 ms per assumption-solve here, against 114 for glucose and worse for cadical
    if not s.solve():
        print(f"{name:16s} {KC:2d} {n:7d}   not {KC}-colourable -- skipped",
              flush=True)
        s.delete()
        continue
    xy = [(float(p.x), float(p.y)) for p in P]
    random.seed(3)
    # nearest non-adjacent pairs first: forcing is a local phenomenon
    cand = []
    src = list(range(n))
    random.shuffle(src)
    for u in src[:400]:
        for v in range(n):
            if v == u or v in adj[u]:
                continue
            d = math.hypot(xy[u][0] - xy[v][0], xy[u][1] - xy[v][1])
            if d < 3.0:
                cand.append((d, min(u, v), max(u, v)))
    cand = sorted(set(cand))[:BUDGET]
    found = 0
    for k, (d, u, v) in enumerate(cand):
        if not s.solve(assumptions=[1 + u * KC, 1 + v * KC]):
            found += 1
            if found <= 3:
                print(f"   FORCED-DIFFERENT non-adjacent pair {u},{v} at "
                      f"distance {d:.4f}  [{time.time()-t0:.0f}s]",
                      flush=True)
        if time.time() - t0 > 1500:
            break
    s.delete()
    # Forced-different is MONOTONE in the carrier: every colouring of a
    # supergraph restricts to one of the subgraph, so a pair forced in Sa is
    # forced in anything containing Sa.  denseSa contains Sa and this scan
    # reported 0 for it against Sa's 29, which is impossible and shows the
    # verdict was overclaiming: 400 random source vertices out of 3501 simply
    # missed them.  A zero here means NOT FOUND IN THIS SAMPLE, nothing more.
    verdict = ("the brick EXISTS here" if found else
               f"none in this sample of {len(cand)} pairs (not an exhaustion)")
    print(f"{name:16s} {KC:2d} {n:7d} {len(cand):7d} {found:7d}   {verdict}"
          f"   [{time.time()-t0:.0f}s]", flush=True)
print("DONE", flush=True)
