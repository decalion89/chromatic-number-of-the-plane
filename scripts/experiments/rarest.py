"""Spend the rarest forced orbit found, whatever field it costs.

Rarity has been the one reliable predictor: the orbit at d^2 = 64/3 occurs
nine times in 153 and its spindle is three times harder to colour than the one
at 64/9, which occurs 144 times, at identical size.  The scan over the orbits
of glue centres turned up rarer ones still --

    Sa[7]    1477 points   1 unusual orbit of 50,  at 64/3   (radical 759)
    Sa[151]  2965 points   1 of 43,                at 64/3
    Sa[139]  3541 points   1 of 101,               at 256/9  (radical 1015)

and 1015 = 5 * 7 * 29, so that last spindle needs sqrt5, sqrt7 and sqrt29 on
top of what Sa already uses: Q(sqrt3, sqrt5, sqrt7, sqrt11, sqrt29), a field of
degree 32 and a fifth radical this project has never gone near.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time, json, random, os
from fractions import Fraction as Fr
from collections import defaultdict
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
from hn.degrey import build_Sa
from hn.geometry import Point, rotation_joining, _rot60
from hn.field import Field
from hn.graph import build_graph
from pysat.solvers import Solver

BASE = int(sys.argv[1]); D2 = Fr(sys.argv[2])
GENS = tuple(int(x) for x in sys.argv[3].split(","))
F = Field(GENS)
print(f"field {F.gens}, dimension {F.dim}; base v{BASE}, spindle d^2={D2}",
      flush=True)
rot60 = _rot60(F); Sa = build_Sa(F)
def orb(p, refl=False):
    out, q = [], p
    for _ in range(6):
        if q not in out: out.append(q)
        q = rot60(q)
    if refl:
        for q in list(out):
            r = Point(q.x, -q.y)
            if r not in out: out.append(r)
    return out
t0 = time.time()
seen, U = set(Sa), list(Sa)
for w in orb(Sa[BASE], True):
    rot = rotation_joining(Fr(1), F).about(w)
    for p in Sa:
        q = rot(p)
        if q not in seen: seen.add(q); U.append(q)
g = build_graph(U); n = g.n
print(f"carrier n={n} m={sum(len(a) for a in g.adj)//2}   [{time.time()-t0:.0f}s]",
      flush=True)
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
s = Solver(name="m22", bootstrap_with=cnf); s.solve()
p0 = set(l for l in s.get_model() if l > 0)
cols = [[next(c for c in range(K) if X(v, c) in p0) for v in range(n)]]
rng = random.Random(5)
while len(cols) < 8:
    s.add_clause([-X(v, cols[-1][v]) for v in rng.sample(range(n), 30)])
    s.set_phases([(1 if rng.random() < 0.25 else -1) * X(v, c)
                  for v in range(n) for c in range(K)])
    if not s.solve(): break
    p2 = set(l for l in s.get_model() if l > 0)
    cols.append([next(c for c in range(K) if X(v, c) in p2) for v in range(n)])
s.delete()
pos = {p: i for i, p in enumerate(g.vertices)}
maps = []
for refl in (False, True):
    base = [Point(p.x, -p.y) for p in g.vertices] if refl else list(g.vertices)
    if any(p not in set(g.vertices) for p in base): continue
    cur = [pos[p] for p in base]
    for _ in range(6):
        maps.append(cur)
        cur = [pos[rot60(g.vertices[v])] for v in cur]
eff = [[c[pm[v]] for v in range(n)] for c in cols for pm in maps]
buck = defaultdict(list)
for v in range(n): buck[tuple(c[v] for c in eff)].append(v)
pair = None
for vs in buck.values():
    for i, a in enumerate(vs):
        for b in vs[i+1:]:
            dd = (U[a] - U[b]).norm2()
            if not dd.is_rational() or Fr(dd.c[0]) != D2: continue
            s2 = Solver(name="m22", bootstrap_with=cnf)
            d = s2.solve(assumptions=[X(a, 0), X(b, 1)]); s2.delete()
            if not d: pair = (a, b); break
        if pair: break
    if pair: break
if pair is None:
    print("  the rare pair did not confirm here", flush=True); sys.exit()
a, b = pair
print(f"  forced pair v{a},v{b} at d^2={D2}   [{time.time()-t0:.0f}s]", flush=True)
spin = rotation_joining(D2, F)
seen2, V = set(U), list(U)
for w in orb(U[a]):
    rot = spin.about(w)
    for p in U:
        z = rot(p)
        if z not in seen2: seen2.add(z); V.append(z)
gv = build_graph(V); nv = gv.n; mv = sum(len(x) for x in gv.adj) // 2
print(f"  spindled over 6 pivots: n={nv} m={mv} deg={2.0*mv/nv:.2f}"
      f"   [{time.time()-t0:.0f}s]", flush=True)
tri = gv.find_clique(3)
for KK in (4, 5):
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
    print(f"  *** refuses {KK} ***", flush=True)
    if KK == 4:
        json.dump({"field_generators": list(F.gens), "n": nv, "m": mv,
                   "base": BASE, "d2": str(D2),
                   "points": [[[[c.numerator, c.denominator] for c in p.x.c],
                               [[c.numerator, c.denominator] for c in p.y.c]]
                              for p in gv.vertices]},
                  open(f"{HN_DIR}/data/rarest_v{BASE}.json", "w"))
        print(f"  written data/rarest_v{BASE}.json", flush=True)
