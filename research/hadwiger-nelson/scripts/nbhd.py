"""The radius-one ring, which every scan here skipped as trivial.

A ring at distance exactly 1 from its centre is capped for free: the centre
touches all of it, so the ring cannot use the centre's colour and its palette
is at most k-1 with no geometry involved.  That is why every scan in this work
skipped it.

But "at most k-1" is the free part.  A radius-one ring whose palette falls
BELOW k-1 is capped by the graph rather than by its centre, and that is not
trivial at all -- it is a cap on a vertex's own neighbourhood, and the
neighbourhood is the one ring every vertex is guaranteed to have.  Nothing
here has ever asked.

The bound still travels upward, so a neighbourhood capped below k-1 anywhere
is capped in every supergraph.  And the neighbourhood is a set a bite can use,
since the turn that stitches a radius-one ring to itself is the sixty-degree
rotation, which is in the field for nothing.

Asked of every vertex, on the whole graph.
"""
import sys, time
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
import pickle
from collections import defaultdict
from hn.degrey import build_G, build_Sa, build_Y
from hn.geometry import DEGREY_FIELD as K
from hn.fast import IntBasis, fast_edges_complete
from hn.homcol import ring_palette_bound

SEED = sys.argv[1] if len(sys.argv) > 1 else "G"
k = int(sys.argv[2]) if len(sys.argv) > 2 else 5
t0 = time.time()
SC = ("/tmp/claude-0/-home-user-darwin-50/"
      "aceaa9ec-f432-5848-a506-39c59179b415/scratchpad/")
if SEED.endswith(".pkl"):
    P = pickle.load(open(SC + SEED, "rb"))
else:
    P = {"Sa": build_Sa, "Y": build_Y,
         "G": lambda f: build_G(f, as_graph=False)}[SEED](K)
b = IntBasis.covering(P)
r = b.rows(P)
assert b.overflow_headroom(r) < 1.0
E = sorted(set((min(a, c), max(a, c)) for a, c in fast_edges_complete(b, r)))
n = len(P)
adj = defaultdict(list)
for a, c in E:
    adj[a].append(c)
    adj[c].append(a)
print(f"{SEED}: {n} points, {len(E)} edges, k={k}; the free bound on a "
      f"neighbourhood is {k-1}  [{time.time()-t0:.0f}s]", flush=True)
cls = [[1 + v * k + c for c in range(k)] for v in range(n)]
for a, c in E:
    for col in range(k):
        cls.append([-(1 + a * k + col), -(1 + c * k + col)])
hist = defaultdict(int)
below = []
order = sorted(range(n), key=lambda v: -len(adj[v]))
for i, v in enumerate(order):
    nb = adj[v]
    if len(nb) < k - 1:
        continue
    p = ring_palette_bound(cls, n * k, nb, k)
    hist[p] += 1
    if p < k - 1:
        below.append((p, v, len(nb)))
        print(f"*** neighbourhood BELOW the free bound: vertex {v}, "
              f"degree {len(nb)}, palette {p} (free bound {k-1})"
              f"  [{time.time()-t0:.0f}s]", flush=True)
    if i % 200 == 0:
        print(f"   .. {i}/{n}, palettes {dict(sorted(hist.items()))}"
              f"  [{time.time()-t0:.0f}s]", flush=True)
print(f"\n{SEED} at k={k}: neighbourhood palettes "
      f"{dict(sorted(hist.items()))}  [{time.time()-t0:.0f}s]", flush=True)
print(f"below the free bound of {k-1}: {len(below)}", flush=True)
print("DONE", flush=True)
