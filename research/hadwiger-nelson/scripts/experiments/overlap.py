"""Forcing needs overlap, so rank the glues by it.

Measured: gluing Sa to a rotated copy that shares ONE point produces no forced
pair at all, over every rational rotation tried; gluing it to the 60-degree
image about a vertex, which shares 224 of 397, produces eight.  That is the
sigma-argument showing through -- with no overlap the second copy is free to
take sigma . c and nothing is forced -- so the search should be ordered by
overlap, not by hinge size or by radical.

The highest-overlap glues are cheap to find:

  * a TRANSLATION by t overlaps in exactly the multiplicity of t as a
    difference of the carrier, so the most frequent differences ARE the best
    translations -- and translations cost no radical at all;
  * the 60-, 120- and 180-degree rotations about a vertex, which inherit the
    carrier's local hexagonal structure.

Each candidate is scored by overlap, the best are built, and the forced-pair
filter is run on each.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time, random, pickle
from fractions import Fraction as Fr
from collections import defaultdict, Counter
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
from hn.degrey import build_Sa, build_Y, build_G
from hn.geometry import Point, rotation_joining
from hn.field import Field
from hn.graph import build_graph
from pysat.solvers import Solver

WHICH, K, TOP = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
NCOL = 20

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

if WHICH in ("Z", "H"):
    F = Field((3, 11, 247))
    base = build_Sa(F)
    r1 = rotation_joining(Fr(1), F).about(base[25])
    seen, H = set(), []
    for p in base:
        for q in (p, r1(p)):
            if q not in seen: seen.add(q); H.append(q)
    if WHICH == "H":
        pts = H
    else:
        r2 = rotation_joining(Fr(64, 9), F).about(H[157])
        seen, pts = set(), []
        for p in H:
            for q in (p, r2(p)):
                if q not in seen: seen.add(q); pts.append(q)
    REP = {1, 3, 11, 33, 247, 741, 2717, 8151}
else:
    F = Field((3, 5, 7, 11))
    pts = {"Sa": build_Sa, "Y": build_Y,
           "G": lambda f: build_G(f, as_graph=False)}[WHICH](F)
    REP = {1, 3, 5, 7, 11, 15, 21, 33, 35, 55, 77, 105, 165, 231, 385, 1155}
n0 = len(pts)
S = set(pts)
print(f"carrier {WHICH}: {n0} points, k={K}", flush=True)

t0 = time.time()
cands = []
# translations: overlap = multiplicity of t as a difference
diff = Counter()
for i, p in enumerate(pts):
    for q in pts:
        if p is not q:
            diff[(p - q)] += 1
print(f"  {len(diff)} distinct difference vectors [{time.time()-t0:.0f}s]", flush=True)
for t, mult in diff.most_common(4 * TOP):
    cands.append((mult, "translate", t, None))
# rotations about a vertex
for rname, d2 in (("rot60", Fr(1)), ("rot120", Fr(1, 3)), ("rot180", Fr(1, 4)),
                  ("moser", Fr(3))):
    if radical(d2) not in REP: continue
    r0 = rotation_joining(d2, F)
    for ci in range(n0):
        rot = r0.about(pts[ci])
        ov = sum(1 for p in pts if rot(p) in S)
        if ov < n0:
            cands.append((ov, rname, pts[ci], rot))
cands.sort(key=lambda c: -c[0])
print(f"  {len(cands)} candidate glues; best overlaps "
      f"{[c[0] for c in cands[:8]]}   [{time.time()-t0:.0f}s]", flush=True)

for rank, (ov, kind, arg, rot) in enumerate(cands[:TOP]):
    if kind == "translate":
        images = [p + arg for p in pts]
    else:
        images = [rot(p) for p in pts]
    seen, pts2 = set(), []
    for p, q in zip(pts, images):
        for z in (p, q):
            if z not in seen: seen.add(z); pts2.append(z)
    g = build_graph(pts2)
    n = g.n
    X = lambda v, c: 1 + v * K + c
    cnf = [[X(v, c) for c in range(K)] for v in range(n)]
    for u, v in g.edges():
        for c in range(K):
            cnf.append([-X(u, c), -X(v, c)])
    s = Solver(name="m22", bootstrap_with=cnf)
    if not s.solve():
        s.delete()
        print(f"  !! rank{rank} {kind} overlap={ov}: n={n} NOT {K}-COLOURABLE "
              f"<<<<<<<<<<<<", flush=True)
        pickle.dump((WHICH, K, kind, rank), open(
            f"/tmp/hn/OHIT_{WHICH}_{K}_{rank}.pkl", "wb"))
        continue
    pos0 = set(l for l in s.get_model() if l > 0)
    cols = [[next(c for c in range(K) if X(v, c) in pos0) for v in range(n)]]
    rng = random.Random(rank + 1)
    tries = 0
    while len(cols) < NCOL and tries < NCOL * 3:
        tries += 1
        s.set_phases([(1 if rng.random() < 1.0 / K else -1) * X(v, c)
                      for v in range(n) for c in range(K)])
        if not s.solve(): continue
        pos = set(l for l in s.get_model() if l > 0)
        cols.append([next(c for c in range(K) if X(v, c) in pos) for v in range(n)])
    s.delete()
    buck = defaultdict(list)
    for v in range(n):
        buck[tuple(col[v] for col in cols)].append(v)
    cand = [(a, b) for vs in buck.values() if len(vs) > 1
            for i, a in enumerate(vs) for b in vs[i+1:]]
    conf = []
    for a, b in cand:
        s2 = Solver(name="m22", bootstrap_with=cnf)
        differ = s2.solve(assumptions=[X(a, 0), X(b, 1)]); s2.delete()
        if not differ:
            dd = (pts2[a] - pts2[b]).norm2()
            conf.append((a, b, str(dd),
                         radical(Fr(dd.c[0])) if dd.is_rational() else None))
    good = [x for x in conf if x[3] in REP]
    mark = "  *** SPINDLEABLE ***" if good else ""
    print(f"  rank{rank:<3} {kind:<10} overlap={ov:<5} n={n:<5} "
          f"cand={len(cand):<4} FORCED={len(conf):<4}{mark}   [{time.time()-t0:.0f}s]",
          flush=True)
    for x in conf[:6]:
        print(f"        v{x[0]} v{x[1]} d^2={x[2]} radical {x[3]}", flush=True)
    if good:
        pickle.dump((WHICH, K, kind, rank, conf), open(
            f"/tmp/hn/OFEQ_{WHICH}_{K}_{rank}.pkl", "wb"))
