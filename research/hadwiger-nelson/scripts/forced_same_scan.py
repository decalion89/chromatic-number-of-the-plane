"""Where does forced-same live, one level down?

A graph is uniquely k-colourable exactly when the forced-same relation has k
classes, so forced-same pairs are the raw material of rigidity. de Grey's G
has none at five colours and neither does the tightest union available -- 8549
queries, zero. Before concluding that rigidity is unreachable, find out where
it DOES occur.

Sa and Y are 4-colourable, so the question is not vacuous at four, and the
triangular lattice is the known rigid example at three, where sqrt(3) forces
agreement. If Sa has forced-same pairs, the gadget behind them can be
extracted the same way the pressure gadget was, and that is the machine to
lift.
"""
import sys, time, random, collections
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
from hn.degrey import build_S, build_Sa, build_Y
from hn.forced import ColourRelations
from hn.graph import build_graph

random.seed(29)
for name, builder, k in (("Sa", build_Sa, 4), ("Y", build_Y, 4)):
    g = builder()
    g = g if hasattr(g, "vertices") else build_graph(g)
    rel = ColourRelations(g, k)
    print(f"{name} {g}  {k}-colourable {rel.colourable}", flush=True)
    if not rel.colourable:
        rel.close()
        continue
    t0, forced, checked = time.time(), [], 0
    order = sorted(range(g.n), key=lambda v: -len(g.adj[v]))
    for n, u in enumerate(order):
        cand = set(range(g.n)) - g.adj[u] - {u}
        for _ in range(5):
            if not cand or not rel._s.solve(assumptions=[rel._x(u, 0)]):
                break
            m = set(rel._s.get_model())
            cand -= {v for v in cand if rel._x(v, 0) not in m}
            rel._s.set_phases([random.choice([1, -1])
                               * rel._x(v, random.randrange(k))
                               for v in random.sample(range(g.n),
                                                      min(80, g.n))])
        for v in sorted(cand):
            checked += 1
            if rel.same(u, v):
                forced.append((u, v))
                if len(forced) <= 6:
                    print(f"  *** FORCED SAME {u},{v}  d^2 = "
                          f"{g.vertices[u].dist2(g.vertices[v])}",
                          flush=True)
        if n % 60 == 59:
            print(f"    ... {n+1}/{g.n}, {checked} queried, {len(forced)} "
                  f"forced same  [{time.time()-t0:.0f}s]", flush=True)
    d2 = collections.Counter(str(g.vertices[u].dist2(g.vertices[v]))
                             for u, v in forced)
    print(f"  {name}: {len(forced)} forced-same pairs over {checked} queries; "
          f"distances {dict(sorted(d2.items(), key=lambda kv: -kv[1])[:6])}  "
          f"[{time.time()-t0:.0f}s]", flush=True)
    rel.close()
