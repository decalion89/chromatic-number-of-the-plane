"""The smallest union whose glue still forces.

The 5-chromatic graph this construction gives is twice the size of the union
that forces, so that union is the thing to minimise -- not the carrier.  And
the union is 2n - overlap, which peeling and glue choice trade against each
other: Sa itself gives 570, while Sa peeled to 340 points with a glue of
overlap 186 gives 494 and still forces eight pairs.

So search the two together.  Peel by lowest degree (the thinning experiment
showed that preserves the forcing where random deletion does not), try every
glue, and only pay for the forcing test on unions smaller than the best so far.
"""
import sys, time, json, random
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
gfull = build_graph(full)
order = sorted(range(len(full)), key=lambda v: len(gfull.adj[v]))

def forced_pairs(pts):
    g = build_graph(pts); n = g.n
    X = lambda v, c: 1 + v * K + c
    cnf = [[X(v, c) for c in range(K)] for v in range(n)]
    for u, v in g.edges():
        for c in range(K):
            cnf.append([-X(u, c), -X(v, c)])
    s = Solver(name="m22", bootstrap_with=cnf)
    if not s.solve(): s.delete(); return None
    pos = set(l for l in s.get_model() if l > 0)
    cols = [[next(c for c in range(K) if X(v, c) in pos) for v in range(n)]]
    r = random.Random(5)
    for _ in range(72):
        if len(cols) >= 24: break
        s.set_phases([(1 if r.random() < 0.25 else -1) * X(v, c)
                      for v in range(n) for c in range(K)])
        if not s.solve(): continue
        p2 = set(l for l in s.get_model() if l > 0)
        cols.append([next(c for c in range(K) if X(v, c) in p2) for v in range(n)])
    s.delete()
    buck = defaultdict(list)
    for v in range(n):
        buck[tuple(c[v] for c in cols)].append(v)
    out = []
    for vs in buck.values():
        for i, a in enumerate(vs):
            for b in vs[i+1:]:
                s2 = Solver(name="m22", bootstrap_with=cnf)
                d = s2.solve(assumptions=[X(a, 0), X(b, 1)]); s2.delete()
                if not d: out.append((a, b))
    return out

ROTS = [("rot60", rotation_joining(Fr(1), F)),
        ("rot120", rotation_joining(Fr(1, 3), F))]
best = 570
t0 = time.time()
for drop in (0, 20, 40, 57, 70, 85, 100, 115, 130, 150):
    keep = sorted(order[drop:])
    pts = [full[i] for i in keep]
    S = set(pts); n0 = len(pts)
    cands = []
    for rname, r0 in ROTS:
        seenov = set()
        for c in pts:
            rot = r0.about(c)
            ov = sum(1 for p in pts if rot(p) in S)
            if ov <= 1 or ov == n0 or ov in seenov: continue
            seenov.add(ov)
            cands.append((2 * n0 - ov, rname, rot, ov))
    cands.sort()
    tried = 0
    for usz, rname, rot, ov in cands:
        if usz >= best: break
        seen, U = set(), []
        for p in pts:
            for q in (p, rot(p)):
                if q not in seen: seen.add(q); U.append(q)
        fp = forced_pairs(U)
        tried += 1
        if fp:
            best = len(U)
            a, b = fp[0]
            d2 = (U[a] - U[b]).norm2()
            print(f"  peel {drop} (n={n0}) {rname} overlap {ov}: union {len(U)} "
                  f"FORCES {len(fp)} pairs, first at d^2={d2}   "
                  f"[{time.time()-t0:.0f}s]  <<<< new best", flush=True)
            json.dump({"field_generators": list(F.gens), "n": len(U),
                       "peel": drop, "glue": rname, "overlap": ov,
                       "forced_pairs": [[int(x), int(y)] for x, y in fp],
                       "points": [[[[c.numerator, c.denominator] for c in p.x.c],
                                   [[c.numerator, c.denominator] for c in p.y.c]]
                                  for p in U]},
                      open("/tmp/claude-0/-home-user-darwin-50/aceaa9ec-f432-5848-a506-39c59179b415/scratchpad/bestunion.json", "w"))
            break
    print(f"  peel {drop}: n={n0}, {tried} glues tried, best union so far {best}"
          f"   [{time.time()-t0:.0f}s]", flush=True)
