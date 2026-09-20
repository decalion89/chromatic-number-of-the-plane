"""Every vertex as a fan centre, and pairs of fans joined together.

The first hunt used the forty highest-degree vertices as centres, but the
witness on Sa is centred at 1/sqrt(3) from the hub -- an ordinary vertex, not
a pivot.  So every vertex gets a turn.

And a structural cap is now clear: the neighbourhood of any vertex is a set of
points on a unit circle, whose own unit-distance graph has maximum degree 2,
being paths and 6-cycles.  So a fan tops out at seven vertices, centre
included, and at five colours seven is probably too few whatever the geometry.

Hence the second family: DOUBLE fans, two centres with their paths, joined
where they meet.  Those reach thirteen vertices and are still small enough to
reject in one SAT call each.
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
    return not solver.solve(assumptions=[-x(v, 0) for v in S])


def longest_paths(c, cap=40):
    nb = sorted(g.adj[c])
    sub = {u: sorted(g.adj[u] & g.adj[c]) for u in nb}
    out, seen = [], 0
    for start in nb:
        stack = [(start, [start])]
        while stack and seen < cap:
            u, path = stack.pop()
            ext = [w for w in sub[u] if w not in path]
            if not ext:
                out.append(tuple(path))
                seen += 1
                continue
            for w in ext:
                stack.append((w, path + [w]))
    return out


t0, tested = time.time(), 0
fans = []
print(f"G: {n} vertices; every vertex as a fan centre", flush=True)
try:
    for c in range(n):
        for path in longest_paths(c, cap=6):
            S = [c] + list(path)
            if len(S) < 4:
                continue
            tested += 1
            fans.append(S)
            if forces(S):
                print(f"  *** FAN OF {len(S)} FORCES: {S} ***", flush=True)
                raise SystemExit
    print(f"  {tested} fans over all centres, none forcing  "
          f"[{time.time()-t0:.0f}s]", flush=True)

    print(f"  now double fans, from {len(fans)} singles", flush=True)
    fans.sort(key=len, reverse=True)
    dbl = 0
    for a, b in itertools.combinations(fans[:600], 2):
        S = sorted(set(a) | set(b))
        if not (8 <= len(S) <= 14):
            continue
        dbl += 1
        if forces(S):
            print(f"  *** DOUBLE FAN OF {len(S)} FORCES: {S} ***", flush=True)
            raise SystemExit
        if dbl >= 40000:
            break
    print(f"  {dbl} double fans, none forcing  [{time.time()-t0:.0f}s]",
          flush=True)
finally:
    solver.delete()
