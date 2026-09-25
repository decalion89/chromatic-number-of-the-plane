"""Count the 5-cycles, which is the screen the analysis actually asks for.

An orientation stays at or below 9/2 only if every 5-cycle gets at least two
edges on its minority side -- pattern 3-2, never 4-1.  A lone 5-cycle can
always manage that.  Overlapping ones sharing edges cannot always, and the
question of chi_c > 9/2 is exactly whether enough of them overlap to make
every acyclic orientation fail somewhere.

So the useful measure of a candidate graph is not its independence ratio --
which the corpus screen showed is a threshold nothing crosses -- but how
many unit-distance 5-cycles it carries and how heavily they share edges.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time, pickle, glob, os
from collections import defaultdict
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
from hn.degrey import build_G, build_Sa, build_Y
from hn.geometry import DEGREY_FIELD as K
from hn.fast import IntBasis, fast_edges_complete

SC = ("/tmp/claude-0/-home-user-darwin-50/"
      "aceaa9ec-f432-5848-a506-39c59179b415/scratchpad/")
t0 = time.time()
cands = [("G", build_G(K, as_graph=False)), ("Sa", build_Sa(K)),
         ("Y", build_Y(K))]
for f in sorted(glob.glob(SC + "*.pkl")):
    try:
        P = pickle.load(open(f, "rb"))
    except Exception:
        continue
    if isinstance(P, list) and 100 <= len(P) <= 20000 and hasattr(P[0], "x"):
        cands.append((os.path.basename(f), P))
print(f"{'graph':20s} {'pts':>6s} {'edges':>7s} {'triangles':>10s} "
      f"{'5-cycles':>10s} {'per edge':>9s} {'shared':>7s}", flush=True)
for name, P in cands:
    b = IntBasis.covering(P)
    r = b.rows(P)
    if b.overflow_headroom(r) >= 1.0:
        continue
    E = sorted(set((min(a, c), max(a, c)) for a, c in fast_edges_complete(b, r)))
    n = len(P)
    if not E or n > 14000:
        continue
    adj = [set() for _ in range(n)]
    for a, c in E:
        adj[a].add(c)
        adj[c].add(a)
    tri = sum(len(adj[a] & adj[c]) for a, c in E) // 3
    # 5-cycles: for each edge (a,c), count paths a - x - y - c of length 3
    # with all five vertices distinct, then divide by the 5 edges each cycle
    # can be opened at.
    cnt = 0
    edge_hits = defaultdict(int)
    for a, c in E:
        for x in adj[a]:
            if x == c:
                continue
            for y in adj[x]:
                if y == a or y == c or y == x:
                    continue
                if c in adj[y]:
                    cnt += 1
                    for e in ((min(a, x), max(a, x)), (min(x, y), max(x, y)),
                              (min(y, c), max(y, c)), (min(a, c), max(a, c))):
                        edge_hits[e] += 1
    fives = cnt // 5
    hot = sum(1 for e in E if edge_hits.get(e, 0) >= 5)
    print(f"{name:20s} {n:6d} {len(E):7d} {tri:10d} {fives:10d} "
          f"{fives/max(1,len(E)):9.3f} {hot/max(1,len(E)):7.3f}"
          f"   [{time.time()-t0:.0f}s]", flush=True)
print("\n'per edge' is 5-cycles divided by edges; 'shared' is the fraction "
      "of edges lying on five or more of them", flush=True)
print("DONE", flush=True)
