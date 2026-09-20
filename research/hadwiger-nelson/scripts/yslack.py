"""Resolve the inconclusive vertices: is there a smaller forcer inside Y?

The triage left 22 of 61 sampled vertices undecided at a 60000-conflict
budget.  That budget is ambiguous by construction -- Y's own unforced pairs
cost up to 35365 conflicts, so exceeding 60000 could mean either "the pair is
still forced and this needs a full refutation" (the vertex is SLACK, and a
smaller forcer exists inside Y) or "the pair separates but the search is hard"
(the vertex is essential).

Run a few of them to the end.  A refutation settles it one way, a model the
other, and either answer is worth the minutes.
"""
import sys, time
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
from hn.degrey import build_Y
from hn.geometry import DEGREY_FIELD as F, Point
from hn.graph import build_graph
from pysat.solvers import Solver

t0 = time.time()
Y = build_Y()
g = build_graph(Y)
E = list(g.edges())
ia = next(i for i, p in enumerate(Y) if p == Point(F.rational(2), F.zero()))
ib = next(i for i, p in enumerate(Y) if p == Point(F.rational(-2), F.zero()))
deg = {}
for a, b in E:
    deg[a] = deg.get(a, 0) + 1
    deg[b] = deg.get(b, 0) + 1
allv = sorted((v for v in range(len(Y)) if v not in (ia, ib)),
              key=lambda v: deg.get(v, 0))
order = allv[::max(1, len(allv) // 60)]
print(f"Y: {len(Y)} vertices; sampling the same 61  [{time.time()-t0:.0f}s]",
      flush=True)


def query(drop, budget):
    keep = [i for i in range(len(Y)) if i != drop]
    idx = {v: i for i, v in enumerate(keep)}
    cls = [[1 + v * 4 + c for c in range(4)] for v in range(len(keep))]
    for a, b in E:
        if a in idx and b in idx:
            for c in range(4):
                cls.append([-(1 + idx[a] * 4 + c), -(1 + idx[b] * 4 + c)])
    with Solver(name="cd19", bootstrap_with=cls) as sv:
        sv.solve()
        if budget is None:
            return sv.solve(assumptions=[1 + idx[ia] * 4,
                                         -(1 + idx[ib] * 4)])
        sv.conf_budget(budget)
        return sv.solve_limited(assumptions=[1 + idx[ia] * 4,
                                             -(1 + idx[ib] * 4)])


undecided = [v for v in order if query(v, 60000) is None]
print(f"  {len(undecided)} undecided at 60000 conflicts: "
      f"{undecided[:12]}  [{time.time()-t0:.0f}s]", flush=True)
for v in undecided[:4]:
    t1 = time.time()
    r = query(v, None)
    print(f"  vertex {v} (degree {deg.get(v,0)}): Y - v "
          f"{'SEPARATES -- essential' if r else 'STILL FORCES -- SLACK'}  "
          f"[{time.time()-t1:.0f}s for this one]", flush=True)
