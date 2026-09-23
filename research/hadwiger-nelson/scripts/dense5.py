"""A 5-chromatic graph at the degree the law says rigidity needs.

free@k ~ (k-1)(1 - 1/(k-1))^deg.  At five colours that needs mean degree about
35 to fall below 1/n at n ~ 5000 -- unreachable -- UNLESS the graph is
structured, and Sa shows structure is worth a factor of two: at four it reaches
free@4 = 0.00 % at degree 9.94 where the law predicts 5.2 %, which a random
graph would need degree 18 to match.  Apply the same factor at five and the
target becomes degree about 19.

Every 5-chromatic graph here has degree 10 to 13.8.  But the 4-chromatic
carriers reach 18.32, by stacking orbits of glue centres -- and the tuned chain
now turns ANY carrier with a forced pair into a 5-chromatic graph without
changing its density: compose the forcing with a rotated copy of the carrier,
tune the composite distance to exactly 1, and the pair becomes an edge.  No
spindle, no extra rotation, and the union is two copies of the carrier, so its
mean degree is the carrier's.

With D^2 = 64/9 the tuning angle is cos phi = -119/128, sin phi = 3 sqrt247/128
-- already the carrier's own field, so the densification costs no new radical
either.

What this measures: free@5 of a 5-chromatic graph at degree 17-18, against the
law.  If structure buys the same factor at five that it buys at four, free@5
should land near 0.1 %, and forcing at five becomes worth hunting there.  If it
tracks the law at ~2 %, the local-rigidity route is closed for good.
"""
import sys, time, json, random
from fractions import Fraction as Fr
from collections import defaultdict
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
from hn.field import Field
from hn.degrey import build_Sa
from hn.geometry import Point, Rotation, rotation_joining, _rot60
from hn.graph import build_graph
from pysat.solvers import Solver

ROOT = "/home/user/darwin-50/research/hadwiger-nelson"
t0 = time.time()
F = Field((3, 11, 247))
rot60 = _rot60(F); Sa = build_Sa(F); g60 = rotation_joining(Fr(1), F)
def orbit(p):
    out, q = [], p
    for _ in range(6):
        if q not in out: out.append(q)
        q = rot60(q)
    for q in list(out):
        r = Point(q.x, -q.y)
        if r not in out: out.append(r)
    return out
BEST = [265, 253, 211, 43, 139, 67, 241, 79, 277, 307]
D2 = Fr(64, 9)
COS = F.rational(Fr(1) / (2 * D2) - 1)
SIN = F.sqrt(247) * F.rational(Fr(3, 128))
assert COS * COS + SIN * SIN == F.rational(1)
R = Rotation(COS, SIN)
def law(deg, K):
    return (K - 1) * ((K - 2) / (K - 1)) ** deg

def cnf_for(gr, K, pin=True):
    N = gr.n
    X = lambda v, c: 1 + v * K + c
    cl = [[X(v, c) for c in range(K)] for v in range(N)]
    for v in range(N):
        for a in range(K):
            for b in range(a + 1, K):
                cl.append([-X(v, a), -X(v, b)])
    for x, y in gr.edges():
        for c in range(K):
            cl.append([-X(x, c), -X(y, c)])
    if pin:
        for i, x in enumerate(gr.find_clique(3) or []):
            cl.append([X(x, i)])
            for c in range(K):
                if c != i:
                    cl.append([-X(x, c)])
    return cl, X

for k in (2, 3, 4, 5):
    cs = []
    for b in BEST[:k]:
        for w in orbit(Sa[b]):
            if w not in cs: cs.append(w)
    seen, U = set(Sa), list(Sa)
    for w in cs:
        rot = g60.about(w)
        for p in Sa:
            q = rot(p)
            if q not in seen: seen.add(q); U.append(q)
    g = build_graph(U); n = g.n; m = sum(len(a) for a in g.adj) // 2
    print(f"\n  {k} orbits: carrier n={n} m={m} deg={2.0*m/n:.2f}"
          f"   [{time.time()-t0:.0f}s]", flush=True)
    # a forced pair at 64/9, found by the filter and confirmed by the solver
    cl, X = cnf_for(g, 4, pin=False)
    s = Solver(name="m22", bootstrap_with=cl)
    if not s.solve():
        print("    (carrier already refuses four -- skip)", flush=True)
        s.delete(); continue
    rng = random.Random(11)
    def read():
        p = set(l for l in s.get_model() if l > 0)
        return [next(c for c in range(4) if X(v, c) in p) for v in range(n)]
    cols = [read()]
    while len(cols) < 8:
        s.add_clause([-X(v, cols[-1][v]) for v in rng.sample(range(n), 30)])
        s.set_phases([(1 if rng.random() < 0.25 else -1) * X(v, c)
                      for v in range(n) for c in range(4)])
        if not s.solve(): break
        cols.append(read())
    buck = defaultdict(list)
    for v in range(n):
        buck[tuple(c[v] for c in cols)].append(v)
    DD = F.rational(D2); pair = None
    for vs in buck.values():
        for i, x in enumerate(vs):
            for y in vs[i+1:]:
                if (g.vertices[x] - g.vertices[y]).norm2() == DD and \
                   not s.solve(assumptions=[X(x, 0), X(y, 1)]):
                    pair = (x, y); break
            if pair: break
        if pair: break
    s.delete()
    if not pair:
        print("    (no forced pair at 64/9 found)", flush=True); continue
    v, q = pair
    V, Q = g.vertices[v], g.vertices[q]
    def tau(p):
        dd = Point(p.x - V.x, p.y - V.y); r = R(dd)
        return Point(Q.x + r.x, Q.y + r.y)
    tq = tau(Q)
    assert (V - tq).norm2() == F.rational(1), "tuning missed distance 1"
    seen2, W = set(g.vertices), list(g.vertices)
    for p in g.vertices:
        z = tau(p)
        if z not in seen2: seen2.add(z); W.append(z)
    g2 = build_graph(W); n2 = g2.n; m2 = sum(len(a) for a in g2.adj) // 2
    deg2 = 2.0 * m2 / n2
    print(f"    chained n={n2} m={m2} deg={deg2:.2f}   [{time.time()-t0:.0f}s]",
          flush=True)
    cl4, _ = cnf_for(g2, 4)
    t1 = time.time()
    s4 = Solver(name="cd19", bootstrap_with=cl4); ok4 = s4.solve(); s4.delete()
    print(f"    4-colourable: {ok4}   [{time.time()-t1:.0f}s]", flush=True)
    if ok4:
        print("    (chain did not close -- skip)", flush=True); continue
    cl5, X5 = cnf_for(g2, 5)
    t1 = time.time()
    s5 = Solver(name="cd19", bootstrap_with=cl5)
    s5.conf_budget(40_000_000)
    r5 = s5.solve_limited()
    if r5 is not True:
        print(f"    *** 5-colour {'UNSAT' if r5 is False else 'budget out'} "
              f"*** [{time.time()-t1:.0f}s]", flush=True)
        if r5 is False:
            json.dump({"field_generators": list(F.gens), "n": n2, "m": m2,
                       "points": [[[[c.numerator, c.denominator] for c in p.x.c],
                                   [[c.numerator, c.denominator] for c in p.y.c]]
                                  for p in g2.vertices]},
                      open(f"{ROOT}/data/six_candidate.json", "w"))
            print("    *** written data/six_candidate.json ***", flush=True)
        s5.delete(); continue
    pos5 = set(l for l in s5.get_model() if l > 0); s5.delete()
    col = [next(c for c in range(5) if X5(u, c) in pos5) for u in range(n2)]
    free = sum(1 for u in range(n2)
               if len(set(col[w] for w in g2.adj[u])) < 4) / n2
    print(f"    5-chromatic, free@5 = {100*free:.3f} %  "
          f"(law {100*law(deg2,5):.3f} %, {law(deg2,5)/max(free,1e-9):.1f}x "
          f"better)   [{time.time()-t1:.0f}s]", flush=True)
    json.dump({"field_generators": list(F.gens), "n": n2, "m": m2,
               "orbits": k, "free_at_5": free,
               "points": [[[[c.numerator, c.denominator] for c in p.x.c],
                           [[c.numerator, c.denominator] for c in p.y.c]]
                          for p in g2.vertices]},
              open(f"{ROOT}/data/five_dense_{k}.json", "w"))
    print(f"    written data/five_dense_{k}.json", flush=True)
