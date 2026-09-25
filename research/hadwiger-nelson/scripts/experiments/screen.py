"""Screen every graph here by independence ratio, which is the real filter.

A map to K(9/2) needs nine independent sets covering each point twice, so
their sizes average (2/9)n = 0.2222n.  A graph whose independent sets are
all comfortably larger than that has plenty of room and will map; one whose
independence ratio is near 2/9 has none.  G's ratio is at least 0.3011,
which is why it maps with room to spare, and the best unit-distance graphs
known reach about 0.2565.

Greedy gives a LOWER bound on alpha, so a high greedy ratio proves the graph
is a poor candidate while a low one proves nothing -- a one-sided filter,
and one-sided in the useful direction, since it eliminates cheaply.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time, pickle, glob, os, random
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
    if isinstance(P, list) and 50 <= len(P) <= 60000 and hasattr(P[0], "x"):
        cands.append((os.path.basename(f), P))
print(f"{'graph':22s} {'points':>7s} {'edges':>8s} {'deg':>6s} "
      f"{'alpha>=':>8s} {'ratio':>7s}   verdict", flush=True)
rows = []
for name, P in cands:
    b = IntBasis.covering(P)
    r = b.rows(P)
    if b.overflow_headroom(r) >= 1.0:
        continue
    E = sorted(set((min(a, c), max(a, c)) for a, c in fast_edges_complete(b, r)))
    n = len(P)
    if not E:
        continue
    adj = [set() for _ in range(n)]
    for a, c in E:
        adj[a].add(c)
        adj[c].add(a)
    best = 0
    random.seed(1)
    for trial in range(30):
        order = sorted(range(n), key=lambda v: (len(adj[v]), random.random()))
        banned, cnt = set(), 0
        for v in order:
            if v not in banned:
                cnt += 1
                banned |= adj[v] | {v}
        best = max(best, cnt)
    ratio = best / n
    verdict = ("room to spare, will map" if ratio > 0.26 else
               "worth a solver run" if ratio > 2 / 9 else
               "REFUSES K(9/2) BY COUNTING")
    print(f"{name:22s} {n:7d} {len(E):8d} {2*len(E)/n:6.2f} "
          f"{best:8d} {ratio:7.4f}   {verdict}", flush=True)
    rows.append((ratio, name))
print(f"\nthreshold for K(9/2) is 2/9 = {2/9:.4f}; "
      f"best known unit-distance graphs reach about 0.2565"
      f"  [{time.time()-t0:.0f}s]", flush=True)
print("tightest first:", ", ".join(f"{nm} {rt:.4f}"
                                   for rt, nm in sorted(rows)[:5]), flush=True)
print("DONE", flush=True)
