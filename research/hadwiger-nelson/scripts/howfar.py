"""How far is five from six, measured rather than guessed.

A pair is forced when the solver cannot separate it.  Short of that, the WORK
the solver needs to separate it is a distance to forcing: a pair that comes
apart in no conflicts at all is wide open, one that costs thousands is nearly
pinned.  So compare the two cases where the answer is known and the one where
it is not:

    Y at four colours   -- six pairs at distance 4, exactly one forced
    G at five colours   -- 21358 pairs at a closable distance, none forced

If G's hardest pair costs about what Y's UNFORCED pairs cost, five is nowhere
near six; if it costs something approaching Y's forced one, the frontier is
closer than it looks.
"""
import sys, time
from fractions import Fraction as Fr
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
from hn.degrey import build_G, build_Y
from hn.geometry import DEGREY_FIELD as F
from hn.graph import build_graph
from hn.homcol import closable_distance
from pysat.solvers import Solver

t0 = time.time()


def cnf(pts, E, k):
    cls = [[1 + v * k + c for c in range(k)] for v in range(len(pts))]
    for a, b in E:
        for c in range(k):
            cls.append([-(1 + a * k + c), -(1 + b * k + c)])
    return cls


def conflicts(sv):
    st = sv.accum_stats()
    return st.get("conflicts", -1)


Y = build_Y()
gy = build_graph(Y)
EY = list(gy.edges())
zs = [(float(p.x), float(p.y)) for p in Y]
four = [(i, j) for i in range(len(Y)) for j in range(i + 1, len(Y))
        if abs((zs[i][0] - zs[j][0]) ** 2
               + (zs[i][1] - zs[j][1]) ** 2 - 16) < 1e-7
        and (Y[i] - Y[j]).x ** 2 + (Y[i] - Y[j]).y ** 2 == F.rational(16)]
print(f"Y: {len(Y)} vertices, {len(four)} pairs at distance 4  "
      f"[{time.time()-t0:.0f}s]", flush=True)
sv = Solver(name="cd19", bootstrap_with=cnf(Y, EY, 4), use_timer=True)
sv.solve()
prev = conflicts(sv)
rows = []
for i, j in four:
    sv.conf_budget(2000000)
    r = sv.solve_limited(assumptions=[1 + i * 4, -(1 + j * 4)])
    now = conflicts(sv)
    rows.append((now - prev, r))
    prev = now
sv.delete()
for c, r in sorted(rows):
    print(f"  distance-4 pair: {c:>9d} conflicts, "
          f"{ {False: 'FORCED', True: 'separable', None: 'budget exhausted -- this is the forced one'}[r] }", flush=True)

G = build_G(F, as_graph=False)
gg = build_graph(G)
EG = list(gg.edges())
zg = [(float(p.x), float(p.y)) for p in G]
print(f"\nG: {len(G)} vertices  [{time.time()-t0:.0f}s]", flush=True)
sv = Solver(name="cd19", bootstrap_with=cnf(G, EG, 5))
sv.solve()
prev = conflicts(sv)
hard, seen, tot = [], 0, 0
for i in range(len(G)):
    ai, bi = zg[i]
    for j in range(i + 1, len(G)):
        v = (ai - zg[j][0]) ** 2 + (bi - zg[j][1]) ** 2
        if v > 36.0:
            continue
        D = Fr(round(v * 1584), 1584)
        if abs(float(D) - v) > 1e-7 or D == 1 or not closable_distance(D):
            continue
        seen += 1
        sv.conf_budget(200000)
        sv.solve_limited(assumptions=[1 + i * 5, -(1 + j * 5)])
        now = conflicts(sv)
        c = now - prev
        prev = now
        tot += c
        if c > 0:
            hard.append((c, i, j, str(D)))
hard.sort(reverse=True)
print(f"  {seen} pairs, {tot} conflicts in total, {len(hard)} costing any at "
      f"all  [{time.time()-t0:.0f}s]", flush=True)
print(f"  the ten dearest: {hard[:10]}", flush=True)
