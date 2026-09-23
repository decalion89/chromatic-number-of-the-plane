"""How few extra distances does it take to break five colours?

No single distance is forced: the 803-graph stays 5-colourable with both 1 and
d forbidden, for every one of its richest distance classes, 3500 extra edges
included.  So ask the next question up, which is still a usable disjunction:

    find a SET D of distances with  chi(H; {1} u D) > 5

and then every 5-colouring of H has a monochromatic pair at some distance in D.
Forbidding everything works trivially, so the content is HOW FEW suffice -- a
narrow D is a narrow disjunction, and a narrow disjunction is what de Grey's
argument consumed at four colours.

Greedy up, then minimise down: add classes richest-first until the instance
dies, then try removing each one and keep the removal when it stays dead.  What
comes back is a minimal forcing set, and its size is the number worth knowing.
"""
import sys, time, json
from fractions import Fraction as Fr
from collections import defaultdict
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
from hn.field import Field
from hn.geometry import Point
from hn.graph import build_graph
from pysat.solvers import Solver

ROOT = "/home/user/darwin-50/research/hadwiger-nelson"
t0 = time.time(); K = 5
NAME = sys.argv[1] if len(sys.argv) > 1 else "five_247_c.json"
POOL = int(sys.argv[2]) if len(sys.argv) > 2 else 40
d = json.load(open(f"{ROOT}/data/{NAME}"))
F = Field(tuple(d["field_generators"]))
P = [Point(F.element([Fr(a, b) for a, b in x]),
           F.element([Fr(a, b) for a, b in y])) for x, y in d["points"]]
g = build_graph(P); n = g.n
E = list(g.edges())
hx = [q.fx for q in g.vertices]; hy = [q.fy for q in g.vertices]
print(f"{NAME} n={n} unit edges={len(E)}   [{time.time()-t0:.0f}s]", flush=True)

buckets = defaultdict(list)
for i in range(n):
    for j in range(i + 1, n):
        dd = (hx[i]-hx[j])**2 + (hy[i]-hy[j])**2
        if dd > 16.0: continue
        k = round(dd, 7)
        if abs(k - 1.0) < 1e-9: continue
        buckets[k].append((i, j))
classes = sorted(buckets.items(), key=lambda kv: -len(kv[1]))[:POOL]
print(f"  pool of {len(classes)} distance classes, "
      f"{classes[0][1].__len__()} down to {classes[-1][1].__len__()} pairs   "
      f"[{time.time()-t0:.0f}s]", flush=True)

X = lambda v, c: 1 + v * K + c
base = [[X(v, c) for c in range(K)] for v in range(n)]
for x, y in E:
    for c in range(K):
        base.append([-X(x, c), -X(y, c)])

def colourable(idxs):
    cnf = list(base)
    for t in idxs:
        for i, j in classes[t][1]:
            for c in range(K):
                cnf.append([-X(i, c), -X(j, c)])
    s = Solver(name="cd19", bootstrap_with=cnf)
    r = s.solve(); s.delete()
    return r

chosen = []
for t in range(len(classes)):
    chosen.append(t)
    ok = colourable(chosen)
    print(f"  |D| = {len(chosen)}  (added d^2={classes[t][0]:.6f}, "
          f"{len(classes[t][1])} pairs)  5-colourable: {ok}   "
          f"[{time.time()-t0:.0f}s]", flush=True)
    if not ok:
        break
else:
    print(f"\n  still 5-colourable with all {len(classes)} classes forbidden   "
          f"[{time.time()-t0:.0f}s]", flush=True)
    sys.exit(0)

print(f"\n  minimising a forcing set of {len(chosen)}   [{time.time()-t0:.0f}s]",
      flush=True)
minimal = list(chosen)
for t in list(chosen):
    trial = [u for u in minimal if u != t]
    if trial and not colourable(trial):
        minimal = trial
        print(f"    dropped d^2={classes[t][0]:.6f}; {len(minimal)} left   "
              f"[{time.time()-t0:.0f}s]", flush=True)
print(f"\n  MINIMAL FORCING SET: {len(minimal)} distances", flush=True)
for t in minimal:
    print(f"      d^2 = {classes[t][0]:.6f}   ({len(classes[t][1])} pairs)", flush=True)
json.dump({"graph": NAME,
           "minimal_forcing_set": [{"d2": classes[t][0],
                                    "pairs": len(classes[t][1])} for t in minimal]},
          open(f"{ROOT}/data/forcing_set_{NAME}", "w"))
print(f"  [{time.time()-t0:.0f}s]", flush=True)
