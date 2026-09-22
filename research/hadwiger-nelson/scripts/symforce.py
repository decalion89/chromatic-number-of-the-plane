"""Forced pairs on a symmetric carrier, with the group doing the work.

The scan of the densest carrier -- 2689 points, mean degree 15.92 -- left 1163
candidates after forty colourings and then spent an hour confirming them one
by one.  Both halves of that are wasteful on an invariant graph.

If c is a proper colouring and g is in the group, c . g is proper too, and it
agrees on (u,v) exactly when c agrees on (g(u), g(v)).  So expanding a handful
of colourings by the group is worth a multiple of them for free, and a pair
survives only if one colouring identifies its whole orbit.  And forcing is
equivariant, so the survivors need one solver call per orbit rather than one
per pair.

The carrier is invariant under the full group of order 12 -- the glue centres
form a reflection-closed orbit -- so both factors are twelve, not six.
"""
import sys, time, random, json
from fractions import Fraction as Fr
from collections import defaultdict, Counter
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
from hn.degrey import build_Sa
from hn.geometry import Point, rotation_joining, _rot60
from hn.field import Field
from hn.graph import build_graph
from pysat.solvers import Solver

F = Field((3, 11, 247)); rot60 = _rot60(F); Sa = build_Sa(F)
g60 = rotation_joining(Fr(1), F)
BASES = [int(x) for x in sys.argv[1].split(",")]
NRAW = int(sys.argv[2]) if len(sys.argv) > 2 else 8
def orbit(p):
    out, q = [], p
    for _ in range(6):
        if q not in out: out.append(q)
        q = rot60(q)
    for q in list(out):
        r = Point(q.x, -q.y)
        if r not in out: out.append(r)
    return out
t0 = time.time()
cs = []
for b in BASES:
    for w in orbit(Sa[b]):
        if w not in cs: cs.append(w)
seen, U = set(Sa), list(Sa)
for w in cs:
    rot = g60.about(w)
    for p in Sa:
        q = rot(p)
        if q not in seen: seen.add(q); U.append(q)
g = build_graph(U); n = g.n; m = sum(len(a) for a in g.adj) // 2
pos = {p: i for i, p in enumerate(g.vertices)}
# the group as vertex permutations: rotations, then those composed with the
# reflection, keeping only the ones that really are symmetries of the set
S = set(g.vertices)
maps = []
for refl in (False, True):
    q = list(range(n))
    base = [Point(p.x, -p.y) for p in g.vertices] if refl else list(g.vertices)
    if any(p not in S for p in base):
        continue
    cur = [pos[p] for p in base]
    for _ in range(6):
        maps.append(cur)
        cur = [pos[rot60(g.vertices[v])] for v in cur]
print(f"carrier {BASES}: n={n} m={m} deg={2.0*m/n:.2f}, group of order {len(maps)}"
      f"   [{time.time()-t0:.0f}s]", flush=True)
K = 4
X = lambda v, c: 1 + v * K + c
cnf = [[X(v, c) for c in range(K)] for v in range(n)]
for u, v in g.edges():
    for c in range(K):
        cnf.append([-X(u, c), -X(v, c)])
for v in range(n):
    for c in range(K):
        for e in range(c + 1, K):
            cnf.append([-X(v, c), -X(v, e)])
# Which engine diversifies depends on the instance, and the candidate count is
# the diagnostic.  On a tight dense carrier at its own chromatic number,
# cadical with blocking gave 61575 candidates where minisat with randomised
# polarity gives a small fraction of that; on a large loose graph the verdict
# reverses, and cadical ignores set_phases entirely.
s = Solver(name="m22", bootstrap_with=cnf)
if not s.solve():
    print("  *** the carrier refuses four colours ***", flush=True); sys.exit()
def readout():
    p = set(l for l in s.get_model() if l > 0)
    return [next(c for c in range(K) if X(v, c) in p) for v in range(n)]
raw = [readout()]
rng = random.Random(11)
while len(raw) < NRAW:
    s.add_clause([-X(v, raw[-1][v]) for v in rng.sample(range(n), 30)])
    s.set_phases([(1 if rng.random() < 1.0 / K else -1) * X(v, c)
                  for v in range(n) for c in range(K)])
    if not s.solve(): break
    raw.append(readout())
eff = [[c[pm[v]] for v in range(n)] for c in raw for pm in maps]
buck = defaultdict(list)
for v in range(n): buck[tuple(c[v] for c in eff)].append(v)
cand = [(a, b) for vs in buck.values() if len(vs) > 1
        for i, a in enumerate(vs) for b in vs[i+1:]]
print(f"  {len(raw)} raw -> {len(eff)} effective colourings, {len(cand)} candidates"
      f"   [{time.time()-t0:.0f}s]", flush=True)
seenp, reps = set(), []
for a, b in cand:
    if (min(a, b), max(a, b)) in seenp: continue
    for pm in maps:
        x, y = pm[a], pm[b]
        seenp.add((min(x, y), max(x, y)))
    reps.append((a, b))
print(f"  -> {len(reps)} orbit representatives   [{time.time()-t0:.0f}s]", flush=True)
byd = defaultdict(int)
tot = 0
s2 = Solver(name="m22", bootstrap_with=cnf)
for a, b in reps:
    if not s2.solve(assumptions=[X(a, 0), X(b, 1)]):
        tot += 1
        dd = (U[a] - U[b]).norm2()
        byd[str(dd) if dd.is_rational() else "irrational"] += 1
s2.delete(); s.delete()
print(f"  FORCED {tot} orbits: {dict(sorted(byd.items(), key=lambda kv: -kv[1])[:6])}"
      f"   [{time.time()-t0:.0f}s]", flush=True)
