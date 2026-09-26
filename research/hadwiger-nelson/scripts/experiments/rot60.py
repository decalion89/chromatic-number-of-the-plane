"""The recipe that produced Z, applied at k colours to any carrier.

  for each vertex v:
      glue the carrier to its image under the 60-degree rotation about v
      (the glue circle is N(v), the unit circle about v -- never capped,
       and the rotation costs nothing: cos 1/2, sin sqrt3/2)
      look for pairs forced to share a colour in the union
      report the ones whose distance admits a spindle

At k=4 on Sa this produces eight forced pairs at distance 8/3 and, after the
spindle, a 1139-vertex 5-chromatic graph.  At k=5 the same step would give a
6-chromatic one.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time, random, pickle
from fractions import Fraction as Fr
from collections import defaultdict
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
from hn.degrey import build_Sa, build_Y, build_G
from hn.geometry import Point, rotation_joining
from hn.field import Field
from hn.graph import build_graph
from pysat.solvers import Solver

WHICH, K = sys.argv[1], int(sys.argv[2])
LO = int(sys.argv[3]); HI = int(sys.argv[4])
D2 = Fr(sys.argv[5]) if len(sys.argv) > 5 else Fr(1)
NCOL = 40

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

if WHICH == "Z":
    F = Field((3, 11, 247))
    base = build_Sa(F)
    r1 = rotation_joining(Fr(1), F).about(base[25])
    seen, H = set(), []
    for p in base:
        for q in (p, r1(p)):
            if q not in seen: seen.add(q); H.append(q)
    r2 = rotation_joining(Fr(64, 9), F).about(H[157])
    seen, pts = set(), []
    for p in H:
        for q in (p, r2(p)):
            if q not in seen: seen.add(q); pts.append(q)
else:
    F = Field((3, 5, 7, 11))
    pts = {"Sa": build_Sa, "Y": build_Y,
           "G": lambda f: build_G(f, as_graph=False)}[WHICH](F)
n0 = len(pts)
print(f"carrier {WHICH}: {n0} points, k={K}, glue radius^2={D2}, "
      f"vertices [{LO},{HI})", flush=True)

rot0 = rotation_joining(D2, F)
t0 = time.time()
for vi in range(LO, min(HI, n0)):
    rot = rot0.about(pts[vi])
    seen, pts2 = set(), []
    for p in pts:
        for q in (p, rot(p)):
            if q not in seen: seen.add(q); pts2.append(q)
    if len(pts2) == n0:
        continue                       # the rotation is a symmetry here
    g = build_graph(pts2)
    n = g.n
    X = lambda v, c: 1 + v * K + c
    cnf = [[X(v, c) for c in range(K)] for v in range(n)]
    for u, v in g.edges():
        for c in range(K):
            cnf.append([-X(u, c), -X(v, c)])
    s = Solver(name="m22", bootstrap_with=cnf)
    rng = random.Random(11)
    cols, tries = [], 0
    while len(cols) < NCOL and tries < NCOL * 4:
        tries += 1
        s.set_phases([(1 if rng.random() < 1.0 / K else -1) * X(v, c)
                      for v in range(n) for c in range(K)])
        if not s.solve(): continue
        pos = set(l for l in s.get_model() if l > 0)
        cols.append([next(c for c in range(K) if X(v, c) in pos) for v in range(n)])
    s.delete()
    if not cols:
        print(f"  v{vi}: n={n}  NOT {K}-COLOURABLE  <<<<<<<<<<<<<<<<", flush=True)
        continue
    buck = defaultdict(list)
    for v in range(n):
        buck[tuple(col[v] for col in cols)].append(v)
    cand = [(a, b) for vs in buck.values() if len(vs) > 1
            for i, a in enumerate(vs) for b in vs[i + 1:]]
    if not cand:
        if vi % 25 == 0:
            print(f"  v{vi}: n={n} shared={2*n0-n} cand=0   [{time.time()-t0:.0f}s]",
                  flush=True)
        continue
    conf = []
    for a, b in cand:
        s2 = Solver(name="m22", bootstrap_with=cnf)
        differ = s2.solve(assumptions=[X(a, 0), X(b, 1)])
        s2.delete()
        if not differ:
            dd = (pts2[a] - pts2[b]).norm2()
            conf.append((a, b, str(dd),
                         radical(Fr(dd.c[0])) if dd.is_rational() else None))
    if conf:
        print(f"  v{vi}: n={n} shared={2*n0-n} cand={len(cand)} "
              f"FORCED={len(conf)}  <<<<   [{time.time()-t0:.0f}s]", flush=True)
        for x in conf[:8]:
            print(f"      v{x[0]} v{x[1]}  d^2={x[2]}  spindle radical {x[3]}", flush=True)
        pickle.dump((vi, conf), open(
            f"/tmp/hn/rot60_{WHICH}_{K}_{vi}.pkl", "wb"))
    elif vi % 25 == 0:
        print(f"  v{vi}: n={n} shared={2*n0-n} cand={len(cand)} FORCED=0"
              f"   [{time.time()-t0:.0f}s]", flush=True)
print(f"done [{time.time()-t0:.0f}s]", flush=True)
