"""The cycles that hold G at four and a half.

In its best orientation G scores exactly 4.5, so some cycle realises that
ratio: |C| / min(|C+|, |C-|) = 9/2, which means a 9-cycle with two edges the
minority way, or an 18-cycle with four, or any multiple.  Those cycles are
what pins the graph to the rung.  Pushing chi_c above 4.5 means making the
ratio unachievable, so knowing what the witnesses look like -- how long they
are, whether they are geometric objects or sprawling combinatorial ones --
says whether there is anything there to attack.
"""
import sys, time, collections
import numpy as np
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
from hn.degrey import build_G
from hn.geometry import DEGREY_FIELD as K
from hn.fast import IntBasis, fast_edges_complete
from pysat.solvers import Solver

t0 = time.time()
P = build_G(K, as_graph=False)
b = IntBasis.covering(P)
rows = b.rows(P)
E = sorted(set((min(a, c), max(a, c)) for a, c in fast_edges_complete(b, rows)))
n, m, p, q = len(P), len(E), 9, 2
xy = np.array([[float(t.x), float(t.y)] for t in P])
EA = np.array([a for a, _ in E])
EB = np.array([c for _, c in E])
cls = [[1 + v * p + j for j in range(p)] for v in range(n)]
for a, c in E:
    for j in range(p):
        for d in range(-(q - 1), q):
            cls.append([-(1 + a * p + j), -(1 + c * p + (j + d) % p)])
s = Solver(name="cd15", bootstrap_with=cls)
assert s.solve()
mod = s.get_model()
pos = np.zeros(n, dtype=np.int64)
for v in range(n):
    for j in range(p):
        if mod[v * p + j] > 0:
            pos[v] = j
            break
orient = pos[EA] < pos[EB]
print(f"G: {n} pts, {m} edges; orientation from the K(9/2) map"
      f"  [{time.time()-t0:.0f}s]", flush=True)

# A cycle realises the ratio exactly when, with r = 9/2, its total weight is
# zero: forward edges r-1, backward -1.  Tight cycles therefore live in the
# subgraph of arcs that are tight for the shortest-path potential at r.
r = 4.5
src = np.where(orient, EA, EB)
dst = np.where(orient, EB, EA)
u = np.concatenate([src, dst])
v = np.concatenate([dst, src])
w = np.concatenate([np.full(m, r - 1.0), np.full(m, -1.0)])
d = np.zeros(n)
for it in range(n + 2):
    nd = d.copy()
    np.minimum.at(nd, v, d[u] + w)
    if np.array_equal(nd, d):
        break
    d = nd
print(f"potential settled after {it} rounds, range "
      f"[{d.min():.2f}, {d.max():.2f}]  [{time.time()-t0:.0f}s]", flush=True)
tight = np.nonzero(d[u] + w <= d[v] + 1e-9)[0]
print(f"{len(tight)} tight arcs of {2*m}  [{time.time()-t0:.0f}s]", flush=True)
nxt = collections.defaultdict(list)
for i in tight:
    nxt[int(u[i])].append((int(v[i]), int(i)))
found, seen = [], set()
for s0 in list(nxt)[:400]:
    # shortest tight cycle through s0, by breadth-first search
    par, dq = {s0: None}, collections.deque([s0])
    hit = None
    while dq and hit is None:
        x = dq.popleft()
        for y, ai in nxt.get(x, []):
            if y == s0:
                hit = (x, ai)
                break
            if y not in par:
                par[y] = (x, ai)
                dq.append(y)
    if hit is None:
        continue
    cyc, x = [hit[1]], hit[0]
    while par[x] is not None:
        px, ai = par[x]
        cyc.append(ai)
        x = px
    verts = tuple(sorted({int(u[a]) for a in cyc}))
    if verts in seen:
        continue
    seen.add(verts)
    fwd = sum(1 for a in cyc if a < m)
    found.append((len(cyc), fwd, verts))
found.sort()
print(f"\n{len(found)} distinct tight cycles found; shortest ten:",
      flush=True)
for L, fwd, verts in found[:10]:
    pts = xy[list(verts)]
    lo = min(fwd, L - fwd)
    print(f"   length {L:3d}, minority side {lo:2d}, ratio "
          f"{L/lo if lo else float('inf'):.3f}; diameter "
          f"{np.max(np.linalg.norm(pts[:,None]-pts[None],axis=2)):.3f}",
          flush=True)
hist = collections.Counter(L for L, _, _ in found)
print(f"\nlengths: {dict(sorted(hist.items()))}", flush=True)
print("DONE", flush=True)
