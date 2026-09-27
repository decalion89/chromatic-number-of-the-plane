"""Hunt for a unit-distance graph with independence ratio below 1/4.

By averaging over random isometries, any finite unit-distance graph H gives
m_1 <= alpha(H)/n(H), where m_1 is the largest density of a measurable set
avoiding distance 1.  Four measurable colour classes cover the plane only if
4*m_1 >= 1, so a graph with ratio below 0.25 proves the measurable chromatic
number is at least 5 -- by density alone, with no finite obstruction and no
spindle.  That is a third route to a theorem currently known two ways
(Falconer 1981 for measurable colourings, de Grey 2018 in general).

The record is m_1 <= 0.2544.  The threshold is 0.25.  The gap is 1.7 per
cent, and this project's best graph sits at 0.2531 -- above it.

Greedy bounds alpha from BELOW, so a greedy ratio below 0.25 proves nothing
by itself but is the only thing worth spending an exact computation on.  The
lever not yet pulled: high k-cores.  Peeling low-degree vertices raises the
mean degree and should lower the ratio, and a core is still a unit-distance
graph, so its ratio bounds m_1 just as well.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time, pickle, glob, os, random
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
from hn.degrey import build_G, build_Sa, build_Y
from hn.geometry import DEGREY_FIELD as K
from hn.fast import IntBasis, fast_edges_complete

SC = ("/tmp/hn/")
t0 = time.time()
cands = [("G", build_G(K, as_graph=False)), ("Sa", build_Sa(K))]
for f in sorted(glob.glob(SC + "*.pkl")):
    try:
        P = pickle.load(open(f, "rb"))
    except Exception:
        continue
    if isinstance(P, list) and 200 <= len(P) <= 42000 and hasattr(P[0], "x"):
        cands.append((os.path.basename(f), P))
print(f"{'graph':18s} {'k':>3s} {'points':>7s} {'deg':>6s} {'alpha>=':>8s} "
      f"{'ratio>=':>8s}   below 0.25?", flush=True)
best = (1.0, None)
for name, P in cands:
    b = IntBasis.covering(P)
    r = b.rows(P)
    if b.overflow_headroom(r) >= 1.0:
        continue
    E = sorted(set((min(a, c), max(a, c)) for a, c in fast_edges_complete(b, r)))
    n0 = len(P)
    if not E:
        continue
    adj0 = [set() for _ in range(n0)]
    for a, c in E:
        adj0[a].add(c)
        adj0[c].add(a)
    for k in (0, 6, 8, 10, 12, 14, 16):
        keep = set(range(n0))
        if k:
            while True:
                drop = {v for v in keep if len(adj0[v] & keep) < k}
                if not drop:
                    break
                keep -= drop
        if len(keep) < 100:
            continue
        ks = sorted(keep)
        idx = {v: i for i, v in enumerate(ks)}
        adj = [set() for _ in range(len(ks))]
        m = 0
        for v in ks:
            for u in adj0[v] & keep:
                adj[idx[v]].add(idx[u])
                m += 1
        m //= 2
        n = len(ks)
        random.seed(2)
        bestset = 0
        for trial in range(25):
            order = sorted(range(n),
                           key=lambda v: (len(adj[v]), random.random()))
            banned, cnt = set(), 0
            for v in order:
                if v not in banned:
                    cnt += 1
                    banned |= adj[v] | {v}
            bestset = max(bestset, cnt)
        ratio = bestset / n
        mark = "  <<< BELOW" if ratio < 0.25 else ""
        print(f"{name:18s} {k:3d} {n:7d} {2*m/n:6.2f} {bestset:8d} "
              f"{ratio:8.4f}{mark}   [{time.time()-t0:.0f}s]", flush=True)
        if ratio < best[0]:
            best = (ratio, f"{name} k={k}, {n} points")
print(f"\nlowest greedy ratio: {best[0]:.4f} on {best[1]}", flush=True)
print(f"threshold for a new proof of measurable chi >= 5 is 0.25; "
      f"the published record is 0.2544", flush=True)
print("DONE", flush=True)
