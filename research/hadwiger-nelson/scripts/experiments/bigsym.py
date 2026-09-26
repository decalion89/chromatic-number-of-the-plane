"""Escalate the symmetric glue: the whole group, and several orbits of centres.

Gluing Sa at the six vertices of ONE C6 orbit gave 1021 points, mean degree
13.34, and 153 forced pairs where the sequential chain reached 52.  Two ways to
push it, both of which keep the invariance for the same reason -- conjugation
sends a glue about w to a glue about g(w):

  * use the full group.  For a reflection g, g rot60_w g^-1 = rot(-60)_{g(w)},
    so gluing with BOTH senses of the rotation at all twelve centres of a D6
    orbit leaves the union invariant under all twelve isometries, not six.
  * use several orbits of centres at once.  The union of orbits is an orbit.

Both cost points, and the question each time is whether the forcing grows
faster than the bulk.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time, random, json
from fractions import Fraction as Fr
from collections import defaultdict
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
from hn.degrey import build_Sa
from hn.geometry import Point, Rotation, rotation_joining, _rot60
from hn.field import Field
from hn.graph import build_graph
from pysat.solvers import Solver

F = Field((3, 11, 247))
rot60 = _rot60(F)
BASES = [int(x) for x in sys.argv[1].split(",")] if len(sys.argv) > 1 else [25]
FULL = (len(sys.argv) > 2 and sys.argv[2] == "D6")
Sa = build_Sa(F)

def orbit(p, reflect):
    out, q = [], p
    for _ in range(6):
        if q not in out: out.append(q)
        q = rot60(q)
    if reflect:
        for q in list(out):
            r = Point(q.x, -q.y)
            if r not in out: out.append(r)
    return out

g60 = rotation_joining(Fr(1), F)
g60i = Rotation(g60.cos, -g60.sin)
rots = [g60, g60i] if FULL else [g60]
centres = []
for b in BASES:
    for w in orbit(Sa[b], FULL):
        if w not in centres: centres.append(w)
t0 = time.time()
seen, U = set(Sa), list(Sa)
for w in centres:
    for r0 in rots:
        rot = r0.about(w)
        for p in Sa:
            q = rot(p)
            if q not in seen: seen.add(q); U.append(q)
g = build_graph(U); n = g.n; m = sum(len(a) for a in g.adj) // 2
tag = ("D6" if FULL else "C6") + " orbits of " + ",".join(f"v{b}" for b in BASES)
print(f"{tag}: {len(centres)} centres x {len(rots)} senses -> "
      f"n={n} m={m} deg={2.0*m/n:.2f}   [{time.time()-t0:.0f}s]", flush=True)
S = set(U)
print(f"  rot60-invariant: {all(rot60(p) in S for p in U)}   "
      f"reflection-invariant: {all(Point(p.x, -p.y) in S for p in U)}", flush=True)

K = 4
X = lambda v, c: 1 + v * K + c
cnf = [[X(v, c) for c in range(K)] for v in range(n)]
for u, v in g.edges():
    for c in range(K):
        cnf.append([-X(u, c), -X(v, c)])
s = Solver(name="m22", bootstrap_with=cnf)
if not s.solve():
    print("  *** the carrier already refuses four colours ***", flush=True)
    s.delete(); sys.exit()
pos = set(l for l in s.get_model() if l > 0)
cols = [[next(c for c in range(K) if X(v, c) in pos) for v in range(n)]]
rng = random.Random(7)
while len(cols) < 16:
    last = cols[-1]
    s.add_clause([-X(v, last[v]) for v in rng.sample(range(n), 40)])
    s.set_phases([(1 if rng.random() < 0.25 else -1) * X(v, c)
                  for v in range(n) for c in range(K)])
    if not s.solve(): break
    p2 = set(l for l in s.get_model() if l > 0)
    cols.append([next(c for c in range(K) if X(v, c) in p2) for v in range(n)])
s.delete()
free = min(sum(1 for v in range(n)
               if len({col[u] for u in g.adj[v]} | {col[v]}) < K) for col in cols)
buck = defaultdict(list)
for v in range(n): buck[tuple(c[v] for c in cols)].append(v)
cand = [(a, b) for vs in buck.values() if len(vs) > 1
        for i, a in enumerate(vs) for b in vs[i+1:]]
print(f"  free@4={100.0*free/n:.2f}%  candidates {len(cand)}"
      f"   [{time.time()-t0:.0f}s]", flush=True)
byd = defaultdict(int)
tot = 0
for a, b in cand[:8000]:
    s2 = Solver(name="m22", bootstrap_with=cnf)
    d = s2.solve(assumptions=[X(a, 0), X(b, 1)]); s2.delete()
    if not d:
        tot += 1
        dd = (U[a] - U[b]).norm2()
        byd[str(dd) if dd.is_rational() else "irrational"] += 1
print(f"  FORCED {tot} pairs: {dict(sorted(byd.items(), key=lambda kv: -kv[1])[:5])}"
      f"   [{time.time()-t0:.0f}s]", flush=True)
json.dump({"n": n, "m": m, "centres": len(centres), "senses": len(rots),
           "free4": 100.0 * free / n, "forced": tot, "by_distance": dict(byd)},
          open(f"/tmp/hn/bigsym_{'-'.join(map(str,BASES))}_{'D6' if FULL else 'C6'}.json", "w"))
