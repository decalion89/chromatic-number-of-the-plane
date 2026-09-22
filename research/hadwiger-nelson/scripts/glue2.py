"""Glue a carrier to a rotated copy at a legal circle and hunt forced pairs.

Established on de Grey's own construction:
  * Sa has no forced-equal pair at 4 colours at all (certified: 40 colourings
    separate all 78606 pairs).
  * Gluing Sa to a rotated copy creates them.
  * Spindling one of them is exactly how G gets its fifth colour.
  * And the glue circle does NOT have to be capped -- uncapped circles
    manufacture forcing just as well, often more of it.

So the route to six colours is: glue G to a rotated copy at some legal circle,
find a forced-equal pair at FIVE colours in the union, and spindle it.  Every
candidate pair here is confirmed by the solver, never left as a filter
survivor.
"""
import sys, time, random, pickle
from fractions import Fraction as Fr
from collections import defaultdict
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
from hn.degrey import build_Sa, build_Y, build_G
from hn.geometry import Point, rotation_joining
from hn.field import Field
from hn.graph import build_graph
from pysat.solvers import Solver

F = Field((3, 5, 7, 11))
ORI = Point(F.zero(), F.zero())
WHICH = sys.argv[1]
K = int(sys.argv[2])
MINR = int(sys.argv[3]) if len(sys.argv) > 3 else 4
LIMIT = int(sys.argv[4]) if len(sys.argv) > 4 else 10 ** 9
SKIP = int(sys.argv[5]) if len(sys.argv) > 5 else 0
NCOL = 40
IN_FIELD = {1, 3, 5, 7, 11, 15, 21, 33, 35, 55, 77, 105, 165, 231, 385, 1155}

def sqfree(r):
    d = 2
    while d * d <= r:
        while r % (d * d) == 0: r //= d * d
        d += 1
    return r

def radical(q):
    if q <= Fr(1, 4): return None
    c = Fr(1) - Fr(1, 2) / q
    s2 = 1 - c * c
    if s2 == 0: return 1
    return sqfree(s2.numerator * s2.denominator)

pts = {"Sa": build_Sa, "Y": build_Y,
       "G": lambda f: build_G(f, as_graph=False)}[WHICH](F)
n0 = len(pts)
print(f"carrier {WHICH}: {n0} points, k={K}", flush=True)

cands = {}
for name, c in [("origin", ORI)] + [(f"v{i}", p) for i, p in enumerate(pts)]:
    b = defaultdict(list)
    for i, p in enumerate(pts):
        d2 = (p - c).norm2()
        if d2.is_rational() and d2 != 0:
            b[Fr(d2.c[0])].append(i)
    for d2, R in b.items():
        if len(R) < MINR: continue
        rad = radical(d2)
        if rad is None or rad not in IN_FIELD: continue
        # a rotation that is already a symmetry of the carrier buys nothing
        cc = Fr(1) - Fr(1, 2) / d2
        if cc in (Fr(1, 2), Fr(-1, 2), Fr(1), Fr(-1)) and name == "origin":
            continue
        key = (d2, tuple(sorted(R)))
        if key not in cands: cands[key] = (name, c)
order = sorted(cands.items(), key=lambda kv: (-len(kv[0][1]), kv[0][0]))
print(f"  {len(order)} distinct legal circles with >= {MINR} points; "
      f"running [{SKIP}:{SKIP+LIMIT}]", flush=True)

t0 = time.time()
hits = []
for (d2, R), (name, c) in order[SKIP:SKIP + LIMIT]:
    rot = rotation_joining(d2, F).about(c)
    seen, pts2 = set(), []
    for p in pts:
        for q in (p, rot(p)):
            if q not in seen:
                seen.add(q); pts2.append(q)
    g = build_graph(pts2)
    n = g.n
    X = lambda v, cc: 1 + v * K + cc
    cnf = [[X(v, cc) for cc in range(K)] for v in range(n)]
    for u, v in g.edges():
        for cc in range(K):
            cnf.append([-X(u, cc), -X(v, cc)])
    s = Solver(name="m22", bootstrap_with=cnf)
    rng = random.Random(7)
    cols, tries = [], 0
    while len(cols) < NCOL and tries < NCOL * 4:
        tries += 1
        s.set_phases([(1 if rng.random() < 1.0 / K else -1) * X(v, cc)
                      for v in range(n) for cc in range(K)])
        if not s.solve(): continue
        pos = set(l for l in s.get_model() if l > 0)
        cols.append([next(cc for cc in range(K) if X(v, cc) in pos) for v in range(n)])
    s.delete()
    if not cols:
        print(f"  r^2={str(d2):<8} |R|={len(R):<3} centre {name:<6} n={n:<5} "
              f"NOT {K}-COLOURABLE  <<<<<<<<<<", flush=True)
        hits.append((str(d2), name, n, "uncolourable"))
        continue
    buck = defaultdict(list)
    for v in range(n):
        buck[tuple(col[v] for col in cols)].append(v)
    cand = [(a, b2) for vs in buck.values() if len(vs) > 1
            for i, a in enumerate(vs) for b2 in vs[i + 1:]]
    conf = []
    for a, b2 in cand:
        s2 = Solver(name="m22", bootstrap_with=cnf)
        differ = s2.solve(assumptions=[X(a, 0), X(b2, 1)])
        s2.delete()
        if not differ:
            dd = (pts2[a] - pts2[b2]).norm2()
            rr = radical(Fr(dd.c[0])) if dd.is_rational() else None
            conf.append((a, b2, str(dd), rr))
    tag = ""
    spin = [x for x in conf if x[3] in IN_FIELD]
    if spin:
        tag = f"   *** {len(spin)} SPINDLEABLE ***"
        hits.append((str(d2), name, n, spin))
    print(f"  r^2={str(d2):<8} |R|={len(R):<3} centre {name:<6} n={n:<5} "
          f"cand={len(cand):<4} FORCED={len(conf):<4}{tag}   [{time.time()-t0:.0f}s]",
          flush=True)
    for x in conf[:6]:
        print(f"        v{x[0]} v{x[1]}  d^2={x[2]}  radical {x[3]}", flush=True)

print(f"\n{len(hits)} circles produced a spindleable forced pair", flush=True)
pickle.dump(hits, open(f"/tmp/claude-0/-home-user-darwin-50/aceaa9ec-f432-5848-a506-39c59179b415/scratchpad/glue2_{WHICH}_{K}_{SKIP}.pkl", "wb"))
