"""Hunt for the five-colour analogue of the fan.

At four colours ambient forcing is carried by a fan: a centre joined to a path
of four, three triangles overlapping, chromatic number 3 and no 4-chromatic
subgraph inside.  Found twice, on two different graphs.

At five colours the analogue would be a small structured set that forces all
five.  Fans are cheap to enumerate -- centre plus a path among its neighbours
-- and cheap to reject, since a set that fails produces an escaping colouring
in one fast SAT call.  Longer paths give bigger fans, and a fan on a path of
p has p + 1 vertices, so this sweeps p upward.

Any hit would be a forcing set far below the 63 a core of three needs.
"""
import sys, time, itertools
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
from pysat.solvers import Solver
from hn import degrey

g = degrey.build_G()
n, k = g.n, 5


def x(v, c):
    return 1 + v * k + c


cls = [[x(v, c) for c in range(k)] for v in range(n)]
for a, b in g.edges():
    for c in range(k):
        cls.append([-x(a, c), -x(b, c)])
solver = Solver(name="cd19", bootstrap_with=cls)


def forces(S):
    """By colour symmetry one colour is enough to test."""
    return not solver.solve(assumptions=[-x(v, 0) for v in S])


def paths_in(adj_sub, nodes, length):
    """Simple paths of `length` vertices inside a neighbourhood."""
    out = []
    for start in nodes:
        stack = [(start, [start])]
        while stack:
            u, path = stack.pop()
            if len(path) == length:
                out.append(tuple(path))
                continue
            for w in adj_sub[u]:
                if w not in path:
                    stack.append((w, path + [w]))
    return out


deg = g.degrees()
centres = sorted(range(n), key=lambda v: -deg[v])
print(f"G: {n} vertices, {g.m} edges; hunting fans at five colours",
      flush=True)
t0, tested, best = time.time(), 0, 0
try:
    for plen in (4, 5, 6, 7, 8):
        for c in centres[:40]:
            nb = sorted(g.adj[c])
            sub = {u: sorted(g.adj[u] & g.adj[c]) for u in nb}
            for path in paths_in(sub, nb, plen)[:300]:
                S = [c] + list(path)
                tested += 1
                if forces(S):
                    print(f"  *** FAN OF {len(S)} FORCES: centre {c}, "
                          f"path {path} ***  [{time.time()-t0:.0f}s]",
                          flush=True)
                    raise SystemExit
                best = max(best, len(S))
        print(f"  path length {plen}: {tested} fans tested, none forcing  "
              f"[{time.time()-t0:.0f}s]", flush=True)
finally:
    solver.delete()
print(f"largest fan tried: {best} vertices")
