"""The 5-cycle screen, counted correctly this time.

The first version walked a path a-x-y-c and closed it with the edge (a,c).
That is four vertices, not five: it counted QUADRILATERALS, and divided them
by 5.  Checked against an independent enumeration on Sa, which finds 714
triangles, 2712 4-cycles and 15948 5-cycles: the old script's "2169
five-cycles" is exactly 2712 * 4 / 5.  Every figure in that table was a
4-cycle count in disguise.

Counted properly here: walk four edges from a vertex and close the fifth,
each cycle seen once by requiring the start to be the smallest vertex and
the second vertex smaller than the last, which kills the rotations and the
reflection.
"""
import sys, time, pickle, glob, os
from collections import defaultdict
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
from hn.degrey import build_G, build_Sa, build_Y
from hn.geometry import DEGREY_FIELD as K
from hn.fast import IntBasis, fast_edges_complete

SC = ("/tmp/claude-0/-home-user-darwin-50/"
      "aceaa9ec-f432-5848-a506-39c59179b415/scratchpad/")
LIMIT = int(sys.argv[1]) if len(sys.argv) > 1 else 12000
t0 = time.time()
cands = [("G", build_G(K, as_graph=False)), ("Sa", build_Sa(K)),
         ("Y", build_Y(K))]
for f in sorted(glob.glob(SC + "*.pkl")):
    try:
        P = pickle.load(open(f, "rb"))
    except Exception:
        continue
    if isinstance(P, list) and 100 <= len(P) <= LIMIT and hasattr(P[0], "x"):
        cands.append((os.path.basename(f), P))
print(f"{'graph':20s} {'pts':>6s} {'edges':>7s} {'tri':>8s} {'4-cyc':>9s} "
      f"{'5-cyc':>9s} {'5/edge':>8s}", flush=True)
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
    tri = sum(len(adj[a] & adj[c]) for a, c in E) // 3
    c4 = c5 = 0
    for a in range(n):
        na = [u for u in adj[a] if u > a]
        for x in na:
            for y in adj[x]:
                if y <= a or y == x:
                    continue
                for z in adj[y]:
                    if z <= a or z in (x, y):
                        continue
                    if a in adj[z] and x < z:
                        c4 += 1
                    for w in adj[z]:
                        if w <= a or w in (x, y, z):
                            continue
                        if a in adj[w] and x < w:
                            c5 += 1
    print(f"{name:20s} {n:6d} {len(E):7d} {tri:8d} {c4:9d} {c5:9d} "
          f"{c5/len(E):8.3f}   [{time.time()-t0:.0f}s]", flush=True)
    rows.append((c5 / len(E), name, n))
print("\nrichest in 5-cycles per edge: " +
      ", ".join(f"{nm} {rt:.3f} ({pt} pts)"
                for rt, nm, pt in sorted(rows, reverse=True)[:4]), flush=True)
print("DONE", flush=True)
