"""Spend all the forcing at once instead of one pair at a time.

H has eight pairs forced to share a colour at four colours, all at squared
distance 64/9.  Spindling any ONE of them gives a 5-chromatic graph.  The
other seven are simply discarded, which is odd: each one is an independent
reason the carrier is rigid, and each spindle rotation fixes a different
pivot, so the copies land differently.

    U_j  =  H  u  rho_1(H)  u  ...  u  rho_j(H)

where rho_i is the spindle at the i-th forced pair.  Every U_j refuses four
colours for j >= 1.  The question is whether piling them up reaches five.

This is not the redundancy experiment that failed earlier: those were
TRANSLATED copies of a finished graph, sharing almost nothing.  Here every
copy shares the whole of H.
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

K1 = Field((3, 11, 247))
MAXJ = int(sys.argv[1]) if len(sys.argv) > 1 else 8
pts = build_Sa(K1)
r1 = rotation_joining(Fr(1), K1).about(pts[25])
seen, H = set(), []
for p in pts:
    for q in (p, r1(p)):
        if q not in seen: seen.add(q); H.append(q)
g = build_graph(H)
n, K = g.n, 4
X = lambda v, c: 1 + v * K + c
cnf = [[X(v, c) for c in range(K)] for v in range(n)]
for u, v in g.edges():
    for c in range(K):
        cnf.append([-X(u, c), -X(v, c)])
s = Solver(name="m22", bootstrap_with=cnf)
if not s.solve(): sys.exit("H is not 4-colourable")
pos = set(l for l in s.get_model() if l > 0)
cols = [[next(c for c in range(K) if X(v, c) in pos) for v in range(n)]]
rng = random.Random(3)
for _ in range(60):
    if len(cols) >= 24: break
    s.set_phases([(1 if rng.random() < 0.25 else -1) * X(v, c)
                  for v in range(n) for c in range(K)])
    if not s.solve(): continue
    p2 = set(l for l in s.get_model() if l > 0)
    cols.append([next(c for c in range(K) if X(v, c) in p2) for v in range(n)])
s.delete()
buck = defaultdict(list)
for v in range(n):
    buck[tuple(c[v] for c in cols)].append(v)
pairs = []
for vs in buck.values():
    for i, a in enumerate(vs):
        for b in vs[i+1:]:
            s2 = Solver(name="m22", bootstrap_with=cnf)
            d = s2.solve(assumptions=[X(a, 0), X(b, 1)]); s2.delete()
            if not d: pairs.append((a, b))
print(f"H: n={n}; {len(pairs)} forced pairs at four colours", flush=True)

t0 = time.time()
seen, U = set(), list(H)
seen.update(H)
for j, (a, b) in enumerate(pairs[:MAXJ], start=1):
    d2 = (H[a] - H[b]).norm2()
    rot = rotation_joining(Fr(d2.c[0]), K1).about(H[a])
    for p in H:
        q = rot(p)
        if q not in seen: seen.add(q); U.append(q)
    gu = build_graph(U)
    mu = sum(len(x) for x in gu.adj) // 2
    # The 4-colour verdict is known by construction -- every U_j contains a
    # spindled forced pair -- and proving it again costs a full UNSAT search on
    # thousands of vertices.  Only the five-colour question is open.
    XX = lambda v, c: 1 + v * 5 + c
    four = "by construction"
    for KK in (5,):
        XX = lambda v, c, KK=KK: 1 + v * KK + c
        cc = [[XX(v, c) for c in range(KK)] for v in range(gu.n)]
        for u, v in gu.edges():
            for c in range(KK):
                cc.append([-XX(u, c), -XX(v, c)])
        sv = Solver(name="m22", bootstrap_with=cc)
        ok = sv.solve()
        model = sv.get_model() if ok else None
        sv.delete()
        if True:
            five = ok
            if ok:
                pz = set(l for l in model if l > 0)
                col = [next(c for c in range(5) if XX(v, c) in pz) for v in range(gu.n)]
                fr = sum(1 for v in range(gu.n)
                         if len({col[u] for u in gu.adj[v]} | {col[v]}) < 5)
    print(f"  j={j}: n={gu.n} m={mu} deg={2.0*mu/gu.n:.2f} "
          f"4-col={four} 5-col={five} "
          f"free@5={'--' if not five else f'{100.0*fr/gu.n:.2f}%'}   "
          f"[{time.time()-t0:.0f}s]", flush=True)
    if not five:
        print("  *** SIX COLOURS ***", flush=True); break
