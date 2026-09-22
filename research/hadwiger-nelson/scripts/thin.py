"""Is rigidity the gate?

Two data points suggest it: Sa is perfectly rigid at four colours (0.00% of its
vertices have a spare colour) and its glue manufactures eight forced pairs; G
is loose at five (17%) and its glue manufactures none.  If rigidity really is
what the glue needs, then the 10% slack every five-colour carrier here shows is
the whole obstruction and the programme is to remove it.  If forcing survives
in a loose carrier, the slack is a red herring and the obstruction is elsewhere.

The test thins Sa -- deleting vertices raises its slack at four colours in a
controlled way -- and at each thinning glues and counts the forced pairs.  The
deletion also weakens the graph directly, so the control is the curve of forced
pairs against slack, not either number alone.
"""
import sys, time, random
from fractions import Fraction as Fr
from collections import defaultdict
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
from hn.degrey import build_Sa
from hn.geometry import rotation_joining
from hn.field import Field
from hn.graph import build_graph
from pysat.solvers import Solver

F = Field((3, 11, 247))
K = 4
full = build_Sa(F)
rng = random.Random(97)
import os
MODE = os.environ.get("MODE", "random")

def measure(pts, tag):
    g = build_graph(pts)
    n = g.n
    X = lambda v, c: 1 + v * K + c
    cnf = [[X(v, c) for c in range(K)] for v in range(n)]
    for u, v in g.edges():
        for c in range(K):
            cnf.append([-X(u, c), -X(v, c)])
    s = Solver(name="m22", bootstrap_with=cnf)
    if not s.solve():
        s.delete(); return None, None, n, g
    pos = set(l for l in s.get_model() if l > 0)
    cols = [[next(c for c in range(K) if X(v, c) in pos) for v in range(n)]]
    r2 = random.Random(5)
    for _ in range(60):
        if len(cols) >= 24: break
        s.set_phases([(1 if r2.random() < 1.0 / K else -1) * X(v, c)
                      for v in range(n) for c in range(K)])
        if not s.solve(): continue
        p2 = set(l for l in s.get_model() if l > 0)
        cols.append([next(c for c in range(K) if X(v, c) in p2) for v in range(n)])
    s.delete()
    free = min(sum(1 for v in range(n)
                   if len({col[u] for u in g.adj[v]} | {col[v]}) < K) for col in cols)
    buck = defaultdict(list)
    for v in range(n):
        buck[tuple(c[v] for c in cols)].append(v)
    nf = 0
    for vs in buck.values():
        for i, a in enumerate(vs):
            for b in vs[i+1:]:
                s2 = Solver(name="m22", bootstrap_with=cnf)
                d = s2.solve(assumptions=[X(a, 0), X(b, 1)]); s2.delete()
                if not d: nf += 1
    return 100.0 * free / n, nf, n, g

# the glue that works on the full carrier: the 60-degree rotation about Sa[25]
GLUE = rotation_joining(Fr(1), F).about(full[25])
print(f"{'deleted':>8} {'n':>5} {'deg':>6} {'free@4':>8} {'forced in the glue':>20}",
      flush=True)
t0 = time.time()
for drop in (1, 3, 4, 5, 6, 7, 8, 9, 10, 12):
    keep = list(range(len(full)))
    if MODE == "random":
        rng2 = random.Random(1000 + drop)
        for _ in range(drop):
            keep.remove(rng2.choice(keep))
    else:
        g_ = build_graph(full)
        keep = sorted(keep, key=lambda v: len(g_.adj[v]))[drop:]
    pts = [full[i] for i in keep]
    g0 = build_graph(pts)
    m0 = sum(len(a) for a in g0.adj) // 2
    fr, _, _, _ = measure(pts, "base")
    if fr is None:
        print(f"{drop:>8} {g0.n:>5} -- base is not 4-colourable", flush=True); continue
    # The glue must be the SAME map at every thinning, or the comparison is
    # meaningless.  A first version re-chose the best-overlap glue for each
    # thinned carrier and picked near-symmetries -- overlap 393 of 395, union
    # 397 -- so the zeros it reported measured that choice and nothing else.
    rot = GLUE
    ov = sum(1 for p in pts if rot(p) in set(pts))
    seen, U = set(), []
    for p in pts:
        for q in (p, rot(p)):
            if q not in seen: seen.add(q); U.append(q)
    _, nf, nu, _ = measure(U, "glued")
    print(f"{drop:>8} {g0.n:>5} {2.0*m0/g0.n:>6.2f} {fr:>7.2f}% "
          f"{nf if nf is not None else '--':>20}   (glue overlap {ov}, union {nu})"
          f"   [{time.time()-t0:.0f}s]", flush=True)
