"""The full symmetric chain on the best carrier found.

Density turned out not to be the lever: the glue takes the mean degree from
9.94 to 18.32 and every carrier on the way is still 4-colourable in under a
second, perfectly rigid.  What does grow is the forcing -- 8 pairs, then 52 on
the sequential chain, then 153 on the first symmetric one -- and the forcing is
what the spindle spends.  So measure it properly on the densest carrier, then
spindle the largest forced orbit over its whole orbit of pivots and ask the
only question that matters.

Forty colourings rather than sixteen, because the filter's cost is linear in
the number of colourings and quadratic in what it fails to eliminate.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time, random, json
from fractions import Fraction as Fr
from collections import defaultdict, Counter
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
from hn.degrey import build_Sa
from hn.geometry import Point, rotation_joining, _rot60
from hn.field import Field
from hn.graph import build_graph
from pysat.solvers import Solver

F = Field((3, 11, 247)); rot60 = _rot60(F); Sa = build_Sa(F)
g60 = rotation_joining(Fr(1), F)
BASES = [int(x) for x in sys.argv[1].split(",")]
NPIV = int(sys.argv[2]) if len(sys.argv) > 2 else 6
def orbit(p, full=True):
    out, q = [], p
    for _ in range(6):
        if q not in out: out.append(q)
        q = rot60(q)
    if full:
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
print(f"carrier from {BASES}: n={n} m={m} deg={2.0*m/n:.2f}   [{time.time()-t0:.0f}s]",
      flush=True)
K = 4
X = lambda v, c: 1 + v * K + c
cnf = [[X(v, c) for c in range(K)] for v in range(n)]
for u, v in g.edges():
    for c in range(K):
        cnf.append([-X(u, c), -X(v, c)])
s = Solver(name="m22", bootstrap_with=cnf)
if not s.solve():
    print("  *** the carrier refuses four colours ***", flush=True); sys.exit()
pos = set(l for l in s.get_model() if l > 0)
cols = [[next(c for c in range(K) if X(v, c) in pos) for v in range(n)]]
rng = random.Random(7)
while len(cols) < 40:
    last = cols[-1]
    s.add_clause([-X(v, last[v]) for v in rng.sample(range(n), 40)])
    s.set_phases([(1 if rng.random() < 0.25 else -1) * X(v, c)
                  for v in range(n) for c in range(K)])
    if not s.solve(): break
    p2 = set(l for l in s.get_model() if l > 0)
    cols.append([next(c for c in range(K) if X(v, c) in p2) for v in range(n)])
s.delete()
buck = defaultdict(list)
for v in range(n): buck[tuple(c[v] for c in cols)].append(v)
cand = [(a, b) for vs in buck.values() if len(vs) > 1
        for i, a in enumerate(vs) for b in vs[i+1:]]
print(f"  {len(cols)} colourings, {len(cand)} candidates   [{time.time()-t0:.0f}s]",
      flush=True)
byd = defaultdict(list)
for a, b in cand:
    s2 = Solver(name="m22", bootstrap_with=cnf)
    d = s2.solve(assumptions=[X(a, 0), X(b, 1)]); s2.delete()
    if not d:
        dd = (U[a] - U[b]).norm2()
        if dd.is_rational(): byd[Fr(dd.c[0])].append((a, b))
print(f"  FORCED by distance: { {str(k): len(v) for k, v in byd.items()} }"
      f"   [{time.time()-t0:.0f}s]", flush=True)
if not byd: sys.exit()
d2, pairs = max(byd.items(), key=lambda kv: len(kv[1]))
spin = rotation_joining(d2, F)
piv = orbit(U[pairs[0][0]], full=False)[:NPIV]
seen2, V = set(U), list(U)
for w in piv:
    rot = spin.about(w)
    for p in U:
        z = rot(p)
        if z not in seen2: seen2.add(z); V.append(z)
gv = build_graph(V); nv = gv.n; mv = sum(len(a) for a in gv.adj) // 2
print(f"  spindled at d^2={d2} over {len(piv)} pivots: n={nv} m={mv} "
      f"deg={2.0*mv/nv:.2f}   [{time.time()-t0:.0f}s]", flush=True)
tri = gv.find_clique(3)
for KK in (4, 5, 6):
    XX = lambda v, c: 1 + v * KK + c
    cc = [[XX(v, c) for c in range(KK)] for v in range(nv)]
    for u, v in gv.edges():
        for c in range(KK):
            cc.append([-XX(u, c), -XX(v, c)])
    for i, v in enumerate(tri):
        cc.append([XX(v, i)])
        for c in range(KK):
            if c != i: cc.append([-XX(v, c)])
    t1 = time.time()
    sv = Solver(name="cd19", bootstrap_with=cc); ok = sv.solve(); sv.delete()
    print(f"  {KK}-colourable: {ok}   [{time.time()-t1:.0f}s]", flush=True)
    if ok: break
    print(f"  *** IT REFUSES {KK} COLOURS ***", flush=True)
    if KK == 5:
        json.dump({"field_generators": list(F.gens), "n": nv, "m": mv,
                   "points": [[[[c.numerator, c.denominator] for c in p.x.c],
                               [[c.numerator, c.denominator] for c in p.y.c]]
                              for p in gv.vertices]},
                  open(HN_DIR + "/data/six_candidate.json", "w"))
        print("  written data/six_candidate.json", flush=True)
