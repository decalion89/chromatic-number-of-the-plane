"""Local repair on orientations: look for one that scores below a threshold.

An orientation whose every cycle has ratio at most r is a CERTIFICATE that
chi_c <= r -- no search has to be exhaustive for the bound to hold.  So
finding one below 22/5 would answer the pending question from the other
side, without waiting on an UNSAT.

The earlier descent rescored by binary search after every trial flip, about
twenty-one Bellman-Ford passes each, which is why it managed only 200 moves.
Here the threshold is fixed: find a cycle that violates it, flip one of that
cycle's edges, repeat.  One pass per move instead of twenty-one, and only
edges on a violating cycle are ever touched, since no other flip can help.
The threshold is then stepped down, restarting from the last orientation
that met it.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time, random, collections
import numpy as np
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
from fractions import Fraction as Fr
from hn.degrey import build_G
from hn.geometry import DEGREY_FIELD as K
from hn.fast import IntBasis, fast_edges_complete
from pysat.solvers import Solver

SC = ("/tmp/claude-0/-home-user-darwin-50/"
      "aceaa9ec-f432-5848-a506-39c59179b415/scratchpad/")
TARGET = Fr(*map(int, (sys.argv[1] if len(sys.argv) > 1 else "22/5").split("/")))
BUDGET = int(sys.argv[2]) if len(sys.argv) > 2 else 20000
t0 = time.time()
P = build_G(K, as_graph=False)
b = IntBasis.covering(P)
rows = b.rows(P)
E = sorted(set((min(a, c), max(a, c)) for a, c in fast_edges_complete(b, rows)))
n, m = len(P), len(E)
EA = np.array([a for a, _ in E])
EB = np.array([c for _, c in E])
p, q = 9, 2
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
print(f"G: {n} pts, {m} edges; seeded from the K(9/2) map, target "
      f"{TARGET} = {float(TARGET):.4f}  [{time.time()-t0:.0f}s]", flush=True)


def violating(orient, r):
    """A cycle with |C| / |C+| > r, as edge indices, or None."""
    src = np.where(orient, EA, EB)
    dst = np.where(orient, EB, EA)
    u = np.concatenate([src, dst])
    v = np.concatenate([dst, src])
    w = np.concatenate([np.full(m, r - 1.0), np.full(m, -1.0)])
    d = np.zeros(n)
    changed = False
    for _ in range(n + 2):
        nd = d.copy()
        np.minimum.at(nd, v, d[u] + w)
        changed = not np.array_equal(nd, d)
        d = nd
        if not changed:
            return None
    ok = np.nonzero(d[u] + w <= d[v] + 1e-9)[0]
    nxt = collections.defaultdict(list)
    for i in ok:
        nxt[int(u[i])].append((int(v[i]), int(i)))
    colour, path = {}, []
    for s0 in list(nxt):
        if colour.get(s0):
            continue
        stack = [(s0, iter(nxt.get(s0, [])))]
        colour[s0] = 1
        path = []
        while stack:
            x, it = stack[-1]
            nx = next(it, None)
            if nx is None:
                colour[x] = 2
                stack.pop()
                if path:
                    path.pop()
                continue
            y, ai = nx
            path.append(ai)
            if colour.get(y) == 1:
                cyc, k = [path[-1]], len(path) - 1
                while k >= 0 and int(u[cyc[-1]]) != y:
                    k -= 1
                    if k >= 0:
                        cyc.append(path[k])
                return [a % m for a in cyc]
            if colour.get(y) is None:
                colour[y] = 1
                stack.append((y, iter(nxt.get(y, []))))
            else:
                path.pop()
    return []


random.seed(13)
r = 4.5
best = None
while r > float(TARGET) - 1e-6:
    r -= 0.005
    tries = 0
    while tries < BUDGET:
        cyc = violating(orient, r)
        if cyc is None:
            best = r
            print(f"   orientation found at r = {r:.4f} after {tries} "
                  f"flips  [{time.time()-t0:.0f}s]", flush=True)
            break
        if not cyc:
            tries += 1
            continue
        orient[random.choice(cyc)] ^= True
        tries += 1
    else:
        print(f"   stuck at r = {r:.4f} after {BUDGET} flips"
              f"  [{time.time()-t0:.0f}s]", flush=True)
        break
if best is not None and best <= float(TARGET) + 1e-9:
    print(f"\nCERTIFICATE: chi_c(G) <= {best:.4f} <= {float(TARGET):.4f}, "
          f"so G MAPS to K({TARGET})", flush=True)
else:
    print(f"\nno orientation found below {float(TARGET):.4f}; best reached "
          f"{best if best else 4.5:.4f}  [{time.time()-t0:.0f}s]", flush=True)
print("DONE", flush=True)
