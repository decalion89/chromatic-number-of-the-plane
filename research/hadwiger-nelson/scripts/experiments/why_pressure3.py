"""Reverse-engineer the only pressure 3 in the package.

Sa reaches pressure 3 at four colours and nothing here reaches it at five.
Knowing WHY is worth more than another stack of copies: if the mechanism is
nameable it can be built to order instead of hunted.

Take the pivot, keep its neighbourhood, and delete everything else that can be
deleted while the neighbourhood still refuses to be squeezed into two colours.
What survives is the machine.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time, collections
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
from hn.degrey import build_Sa
from hn.forced import ColourRelations, min_colours_on
from hn.graph import build_graph

g = build_Sa()
g = g if hasattr(g, "vertices") else build_graph(g)
bp = 0
circle = sorted(g.adj[bp])
print(f"Sa {g}  pivot {bp}, circle {len(circle)}", flush=True)

keep = set(circle) | {bp}
rest = [v for v in range(g.n) if v not in keep]
t0 = time.time()


def pressure_of(idx):
    pts = [g.vertices[v] for v in sorted(idx)]
    h = build_graph(pts)
    j = h.vertices.index(g.vertices[bp])
    rel = ColourRelations(h, 4)
    if not rel.colourable:
        rel.close()
        return None
    t = min_colours_on(rel, sorted(h.adj[j]))
    rel.close()
    return t


cur = set(range(g.n))
print(f"  full graph pressure: {pressure_of(cur)}", flush=True)
dropped = 0
for v in rest:
    trial = cur - {v}
    if pressure_of(trial) == 3:
        cur = trial
        dropped += 1
    if dropped and dropped % 25 == 0:
        print(f"    ... dropped {dropped}, {len(cur)} left  "
              f"[{time.time()-t0:.0f}s]", flush=True)
print(f"  minimal witness: {len(cur)} vertices keep pressure 3 "
      f"(circle {len(circle)} + pivot + {len(cur)-len(circle)-1} others)  "
      f"[{time.time()-t0:.0f}s]", flush=True)
sub = build_graph([g.vertices[v] for v in sorted(cur)])
print(f"  witness graph: {sub}", flush=True)
pv = g.vertices[bp]
extra = [v for v in sorted(cur) if v != bp and v not in set(circle)]
d2 = collections.Counter(str(pv.dist2(g.vertices[v])) for v in extra)
print(f"  the {len(extra)} non-circle vertices sit at d^2 "
      f"{dict(sorted(d2.items(), key=lambda kv: -kv[1])[:8])}", flush=True)
import pickle
pickle.dump(sorted(cur), open(f"{SC if False else '/tmp/hn'}/witness.pkl", "wb"))
