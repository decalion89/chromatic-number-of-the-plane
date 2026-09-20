"""Make G critical before rotating it.  That is what every negative points at.

Y is nearly 4-critical: 787 of its 789 non-pair vertices are individually
indispensable to its forcing, so its 4-colourings are scarce and six cross
edges are enough to pin seventy-five pairs.  G is 5-chromatic with 1581
vertices where the literature gets the same property from about five hundred,
so it is redundant by a factor of three, its 5-colourings are plentiful, and
nothing is pinned -- 44 distinct biting rotations about arbitrary pivots, up
to 812 cross edges, and not one pair agrees across fourteen colourings.

Criticality is the missing ingredient, and an unsatisfiable core extracts it
in one solve rather than 1581.  Gate each vertex behind a selector: the clause
"v gets a colour" becomes "s_v implies v gets a colour", so not asserting s_v
lets v go uncoloured, which is deletion -- the edge clauses stay satisfied by
an uncoloured endpoint.  Solve with every selector assumed; the proof of
unsatisfiability names a subset of them that already suffices, and that subset
is a smaller graph with the same chromatic number.

Then iterate.  The first solve is de Grey's own theorem and is slow; each one
after it is on a smaller graph.
"""
import sys, time, pickle
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
from hn.degrey import build_G
from hn.geometry import DEGREY_FIELD as K
from hn.graph import build_graph
from pysat.solvers import Solver

SC = ("/tmp/claude-0/-home-user-darwin-50/aceaa9ec-f432-5848-a506-"
      "39c59179b415/scratchpad/")
t0 = time.time()
k = 4
pts = build_G(K, as_graph=False)
g = build_graph(pts)
keep = list(range(g.n))
print(f"G: {g.n} points, {g.m} edges  [{time.time()-t0:.0f}s]", flush=True)

for rnd in range(12):
    sub = g.induced(keep)
    sub = sub.k_core(k)
    idx = {p: i for i, p in enumerate(sub.vertices)}
    n = sub.n
    E = sorted((min(a, b), max(a, b)) for a, b in sub.edges())
    cls, sel = [], []
    for v in range(n):
        s = 1 + n * k + v
        sel.append(s)
        cls.append([-s] + [1 + v * k + c for c in range(k)])
    for a, b in E:
        for c in range(k):
            cls.append([-(1 + a * k + c), -(1 + b * k + c)])
    sv = Solver(name="cd19", bootstrap_with=cls)
    ok = sv.solve(assumptions=sel)
    if ok:
        print(f"  round {rnd}: {n} points is 4-COLOURABLE -- the previous "
              f"round was already minimal for this method", flush=True)
        sv.delete()
        break
    core = set(sv.get_core() or sel)
    sv.delete()
    kept = [v for v in range(n) if sel[v] in core]
    print(f"  round {rnd}: {n} points, {len(E)} edges -> core of "
          f"{len(kept)}  [{time.time()-t0:.0f}s]", flush=True)
    pts2 = [sub.vertices[v] for v in kept]
    g = build_graph(pts2)
    keep = list(range(g.n))
    with open(SC + "shrink.pkl", "wb") as fh:
        pickle.dump([(str(p.x), str(p.y)) for p in g.vertices], fh)
    if len(kept) == n:
        print("  no further shrinkage from the core", flush=True)
        break

final = build_graph([p for p in g.vertices]).k_core(k)
print(f"final: {final.n} points, {final.m} edges  [{time.time()-t0:.0f}s]",
      flush=True)
with open(SC + "shrink.pkl", "wb") as fh:
    pickle.dump([(str(p.x), str(p.y)) for p in final.vertices], fh)
