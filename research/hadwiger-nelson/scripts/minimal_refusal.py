"""Name the half of de Grey's argument the odd cycle does not cover.

Squeezing the pivot's circle into two colours fixes a colouring on it, up to
which two colours and an orientation bit per hexagon. In sixteen of the
thirty-two orientations the confined points -- those seeing two
differently-coloured circle points, hence barred from both -- form a BIPARTITE
graph, so they fit in the two remaining colours with room to spare and the odd
cycle explains nothing there. Something else refuses those orientations,
because the pressure is 3 for all of them.

Fix such an orientation as assumptions, confirm the graph has no 4-colouring
extending it, then delete every vertex that can go while it stays impossible.
What survives is that half of the machine, written out rather than guessed.
"""
import sys, time, collections
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
from pysat.solvers import Solver
from hn.certify import load_certificate
from hn.graph import build_graph

R = "/home/user/darwin-50/research/hadwiger-nelson"
pts, doc = load_certificate(f"{R}/certificates/pressure3_witness_47.json")
g = build_graph(pts)
K, piv = 4, 0
circle = set(g.adj[piv])
extras = [v for v in range(g.n) if v != piv and v not in circle]

seen, comps = set(), []
for s in sorted(circle):
    if s in seen:
        continue
    comp, st = [], [s]
    seen.add(s)
    while st:
        x = st.pop()
        comp.append(x)
        for y in g.adj[x] & circle:
            if y not in seen:
                seen.add(y)
                st.append(y)
    comps.append(set(comp))
par, comp_of = {}, {}
for ci, comp in enumerate(comps):
    s = min(comp)
    par[s], front = 0, [s]
    for v in comp:
        comp_of[v] = ci
    while front:
        x = front.pop()
        for y in g.adj[x] & comp:
            if y not in par:
                par[y] = 1 - par[x]
                front.append(y)
print(f"{g}  circle {len(circle)} in {len(comps)} hexagons, "
      f"{len(extras)} extras", flush=True)


def colouring(bits):
    return {v: par[v] ^ ((bits >> comp_of[v]) & 1) for v in circle}


def impossible(keep, fixed):
    """No K-colouring of the induced subgraph extends this circle colouring."""
    idx = {v: i for i, v in enumerate(sorted(keep))}
    sub = build_graph([g.vertices[v] for v in sorted(keep)])

    def x(i, c):
        return 1 + i * K + c

    cls = [[x(i, c) for c in range(K)] for i in range(sub.n)]
    for a, b in sub.edges():
        for c in range(K):
            cls.append([-x(a, c), -x(b, c)])
    ass = [x(idx[v], fixed[v]) for v in fixed if v in keep]
    with Solver(name="cd19", bootstrap_with=cls) as s:
        return not s.solve(assumptions=ass)


t0 = time.time()
for bits in (0b00000, 0b00010, 0b00100):
    fixed = colouring(bits)
    barred = [v for v in extras
              if len({fixed[c] for c in g.adj[v] & circle}) == 2]
    full = set(range(g.n))
    if not impossible(full, fixed):
        print(f"  orientation {bits:05b}: the squeeze SUCCEEDS -- "
              f"pressure would be 2", flush=True)
        continue
    cur = set(full)
    for v in list(extras) + sorted(circle):
        if v == piv:
            continue
        trial = cur - {v}
        if impossible(trial, {k: c for k, c in fixed.items() if k in trial}):
            cur = trial
    kept_extra = sorted(set(cur) & set(extras))
    kept_circle = sorted(set(cur) & circle)
    sub = build_graph([g.vertices[v] for v in sorted(cur)])
    pv = g.vertices[piv]
    d2 = collections.Counter(str(pv.dist2(g.vertices[v])) for v in kept_extra)
    print(f"  orientation {bits:05b}: {len(barred)} barred (bipartite); "
          f"minimal refusal uses {len(cur)} vertices = {len(kept_circle)} "
          f"circle + {len(kept_extra)} extras{' + pivot' if piv in cur else ''}"
          f"  [{time.time()-t0:.0f}s]", flush=True)
    print(f"      extras kept at d^2 {dict(d2)}", flush=True)
    print(f"      of those, barred: "
          f"{sorted(set(kept_extra) & set(barred))}, free: "
          f"{sorted(set(kept_extra) - set(barred))}", flush=True)
    print(f"      hexagons touched: "
          f"{sorted({comp_of[v] for v in kept_circle})}  subgraph {sub}",
          flush=True)
