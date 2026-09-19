"""Vertex-criticality, vertex by vertex, with the ambiguity stated up front.

chi(G - v) < 5 for every v is what 5-vertex-criticality asserts, and this
package leans on it to explain why nothing is ever forced at five colours.
The evidence offered for it was separability, which is a CONSEQUENCE of
criticality rather than a proof.

Testing it means 1581 four-colourability calls.  A YES is cheap -- the solver
exhibits a colouring -- while a NO is the same hard UNSAT that took kissat 522
seconds on the whole graph.  So each call gets a conflict budget and an
exhausted budget is reported as ambiguous rather than counted either way; a
budget that runs out is weak evidence for NO, since YES is normally quick.

The sample is adversarial rather than random: the highest-degree hubs, most
likely to be essential, and the lowest-degree vertices, most likely to be
removable.
"""
import sys, time
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
from hn import degrey
from pysat.solvers import Solver

g = degrey.build_G()
deg = g.degrees()
order = sorted(range(g.n), key=lambda v: -deg[v])
sample = order[:10] + order[-10:] + order[g.n // 2 - 5:g.n // 2 + 5]

print(f"G: {g.n} vertices.  chi(G - v) <= 4 on {len(sample)} vertices "
      f"(hubs, leaves, middle), budget 3e6 conflicts each.", flush=True)
ok = amb = red = 0
for v in sample:
    sub = g.induced([u for u in range(g.n) if u != v])
    k = 4
    cls = [[1 + i * k + c for c in range(k)] for i in range(sub.n)]
    for a, b in sub.edges():
        for c in range(k):
            cls.append([-(1 + a * k + c), -(1 + b * k + c)])
    t = time.time()
    with Solver(name="g4", bootstrap_with=cls) as s:
        s.conf_budget(3000000)
        r = s.solve_limited()
    el = time.time() - t
    if r is True:
        ok += 1
    elif r is False:
        red += 1
        print(f"  v={v} (deg {deg[v]}): G - v is STILL 5-chromatic -- G is "
              f"NOT vertex-critical  [{el:.0f}s]", flush=True)
    else:
        amb += 1
        print(f"  v={v} (deg {deg[v]}): budget exhausted, ambiguous "
              f"[{el:.0f}s]", flush=True)
print(f"\n{ok} essential (G - v became 4-colourable), {red} redundant, "
      f"{amb} ambiguous")
if red:
    print("=> the vertex-criticality claim is FALSE and must be retracted")
elif amb:
    print("=> not settled: the ambiguous ones are the interesting ones")
else:
    print("=> consistent with criticality on this sample, not a proof for all")
