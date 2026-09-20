"""An UNSAT core answers both questions at once, and faster than deletion.

Give every vertex a selector s_v and weaken its "at least one colour" clause
to (not s_v or x_v0 or ... or x_v3).  Assuming every s_v recovers the original
formula; the solver's UNSAT CORE is then a subset of vertices that is already
uncolourable -- which is to say a 5-CHROMATIC SUBGRAPH of G, found in one
solve instead of 1581.

Two things fall out together:

  If the core is smaller than G, then G minus a vertex outside it is still
  5-chromatic, so G is NOT 5-vertex-critical and the corollary this package
  leans on does not apply to it.

  And rho(G,5) <= |core|, since a 5-chromatic subgraph uses all five colours
  in every 5-colouring.  A core of three at a degree-60 pivot needs rho <= 63,
  so the core's size is the distance left, measured rather than assumed.
"""
import sys, time
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
from hn import degrey
from pysat.solvers import Solver

g = degrey.build_G()
k = 4
NV = g.n * k


def x(v, c):
    return 1 + v * k + c


def sel(v):
    return NV + 1 + v


cls = []
for v in range(g.n):
    cls.append([-sel(v)] + [x(v, c) for c in range(k)])
for a, b in g.edges():
    for c in range(k):
        cls.append([-x(a, c), -x(b, c)])
tri = g.find_clique(3)
for i, u in enumerate(tri):
    cls.append([x(u, i)])
print(f"G: {g.n} vertices, {g.m} edges; pinned triangle {tri}", flush=True)

t = time.time()
with Solver(name="g4", bootstrap_with=cls) as s:
    r = s.solve(assumptions=[sel(v) for v in range(g.n)])
    print(f"  4-colourable with every vertex on: {r}  "
          f"[{time.time()-t:.0f}s]", flush=True)
    if r:
        print("  the graph is 4-colourable -- something is wrong upstream")
        sys.exit()
    core = sorted(-l - NV - 1 if l < 0 else l - NV - 1 for l in s.get_core())
core = sorted({v if v >= 0 else -v for v in core})
print(f"  UNSAT core: {len(core)} of {g.n} vertices  "
      f"[{time.time()-t:.0f}s]", flush=True)
if len(core) < g.n:
    print(f"  *** a proper subgraph of {len(core)} vertices is already "
          f"5-chromatic, so G is NOT 5-vertex-critical ***", flush=True)
    print(f"  and rho(G,5) <= {len(core)}, against the 63 a core of three "
          f"at a degree-60 pivot needs", flush=True)
else:
    print("  the core is the whole graph: consistent with vertex-criticality")

# tighten it: re-solve on the core alone, repeatedly, until it stops shrinking
cur = core
for rnd in range(6):
    sub = g.induced(cur)
    cl2 = []
    for v in range(sub.n):
        cl2.append([-(sub.n * k + 1 + v)] + [1 + v * k + c for c in range(k)])
    for a, b in sub.edges():
        for c in range(k):
            cl2.append([-(1 + a * k + c), -(1 + b * k + c)])
    t2 = time.time()
    with Solver(name="g4", bootstrap_with=cl2) as s:
        if s.solve(assumptions=[sub.n * k + 1 + v for v in range(sub.n)]):
            print(f"    round {rnd+1}: became 4-colourable, stop", flush=True)
            break
        c2 = sorted({abs(l) - sub.n * k - 1 for l in s.get_core()})
    nxt = [cur[i] for i in c2 if 0 <= i < len(cur)]
    print(f"    round {rnd+1}: {len(cur)} -> {len(nxt)}  "
          f"[{time.time()-t2:.0f}s]", flush=True)
    if len(nxt) >= len(cur):
        break
    cur = nxt
print(f"\nsmallest 5-chromatic subgraph found: {len(cur)} vertices")
