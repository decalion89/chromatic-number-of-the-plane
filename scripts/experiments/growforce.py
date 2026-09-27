"""How few vertices of G force all five colours?  Built by adding hubs.

Each test is a single SAT call and comes back instantly, because avoiding a
colour on a small set is easy -- so the space can be swept rather than
sampled.  One hub's closed neighbourhood (61 vertices) does not force; nor do
adjacent pairs.  This adds closed neighbourhoods, highest degree first, until
the union forces all five colours, then shrinks the result greedily.

The number that comes out is an upper bound on rho(G,5), and the target is 63
-- what a core of three at a degree-60 pivot would need.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
from hn import degrey
from pysat.solvers import Solver

g = degrey.build_G()
k = 5


def x(v, c):
    return 1 + v * k + c


base = [[x(v, c) for c in range(k)] for v in range(g.n)]
for a, b in g.edges():
    for c in range(k):
        base.append([-x(a, c), -x(b, c)])


def forcing(S):
    for c in range(k):
        with Solver(name="cd19",
                    bootstrap_with=base + [[-x(v, c)] for v in S]) as s:
            if s.solve():
                return False
    return True


deg = g.degrees()
order = sorted(range(g.n), key=lambda v: -deg[v])
print(f"G: {g.n} vertices, degrees {min(deg)}..{max(deg)}", flush=True)

S, added, t0 = set(), 0, time.time()
for p in order:
    S |= {p} | g.adj[p]
    added += 1
    if added % 5 == 0 or added < 4:
        f = forcing(sorted(S))
        print(f"  {added} neighbourhoods, {len(S)} vertices: forcing={f}"
              f"  [{time.time()-t0:.0f}s]", flush=True)
        if f:
            break
    if len(S) > 1200:
        print(f"  {added} neighbourhoods, {len(S)} vertices and still not "
              f"forcing -- stopping", flush=True)
        break

if forcing(sorted(S)):
    cur = sorted(S)
    print(f"\nshrinking from {len(cur)}:", flush=True)
    changed = True
    while changed:
        changed = False
        for v in list(cur):
            trial = [u for u in cur if u != v]
            if forcing(trial):
                cur = trial
                changed = True
        print(f"  down to {len(cur)}  [{time.time()-t0:.0f}s]", flush=True)
    print(f"\nrho(G,5) <= {len(cur)}, against the 63 a core of three needs")
else:
    print(f"\nno forcing set found this way; rho(G,5) is not reached by "
          f"unions of {added} closed neighbourhoods")
