"""Spindle both forced orbits at once, keeping the symmetry.

The symmetric carrier -- Sa glued at all six vertices of one C6 orbit, 1021
points, mean degree 13.34, perfectly rigid at four colours -- has 153 forced
pairs in TWO orbits:

    144 at d^2 = 64/9,  whose spindle needs sqrt(247) = sqrt(13*19)
      9 at d^2 = 64/3,  whose spindle needs sqrt(759) = sqrt3 sqrt11 sqrt23

Spindling the first orbit alone gives a 7141-vertex graph that is not
4-colourable and whose 5-colouring cadical needs 322 seconds to find, against
4 seconds for a 6607-vertex graph of the same shape -- eighty times harder at
comparable size, which is the first quantitative sign of the boundary.

Both orbits can be spent at once.  Conjugation carries each spindle to another
of its own orbit, so the union of all twelve stays invariant, and the two
rotations are independent constraints on the same carrier rather than one
applied twice.  The price is a third radical: the field becomes
Q(sqrt3, sqrt11, sqrt23, sqrt247).
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

F = Field((3, 11, 23, 247))
print(f"field {F.gens}, dimension {F.dim}", flush=True)
rot60 = _rot60(F)
Sa = build_Sa(F)

def orbit(p):
    out, q = [], p
    for _ in range(6):
        if q not in out: out.append(q)
        q = rot60(q)
    return out

t0 = time.time()
glue = rotation_joining(Fr(1), F)
seen, U = set(Sa), list(Sa)
for w in orbit(Sa[25]):
    rot = glue.about(w)
    for p in Sa:
        q = rot(p)
        if q not in seen: seen.add(q); U.append(q)
gU = build_graph(U)
print(f"symmetric carrier: n={gU.n} m={sum(len(a) for a in gU.adj)//2}"
      f"   [{time.time()-t0:.0f}s]", flush=True)

CACHE = "/tmp/hn/forced1021.json"
import os
K = 4
X = lambda v, c: 1 + v * K + c
cnf = [[X(v, c) for c in range(K)] for v in range(gU.n)]
for u, v in gU.edges():
    for c in range(K):
        cnf.append([-X(u, c), -X(v, c)])
if os.path.exists(CACHE):
    byd = {Fr(k): [tuple(p) for p in v] for k, v in json.load(open(CACHE)).items()}
    print(f"  forced pairs read from cache: "
          f"{ {str(k): len(v) for k, v in byd.items()} }   [{time.time()-t0:.0f}s]",
          flush=True)
else:
  s = Solver(name="m22", bootstrap_with=cnf); s.solve()
  pos = set(l for l in s.get_model() if l > 0)
  cols = [[next(c for c in range(K) if X(v, c) in pos) for v in range(gU.n)]]
  rng = random.Random(7)
  while len(cols) < 16:
      last = cols[-1]
      s.add_clause([-X(v, last[v]) for v in rng.sample(range(gU.n), 40)])
      s.set_phases([(1 if rng.random() < 0.25 else -1) * X(v, c)
                    for v in range(gU.n) for c in range(K)])
      if not s.solve(): break
      p2 = set(l for l in s.get_model() if l > 0)
      cols.append([next(c for c in range(K) if X(v, c) in p2) for v in range(gU.n)])
  s.delete()
  buck = defaultdict(list)
  for v in range(gU.n): buck[tuple(c[v] for c in cols)].append(v)
  byd = defaultdict(list)
  for vs in buck.values():
      for i, a in enumerate(vs):
          for b in vs[i+1:]:
              s2 = Solver(name="m22", bootstrap_with=cnf)
              d = s2.solve(assumptions=[X(a, 0), X(b, 1)]); s2.delete()
              if not d:
                  dd = (U[a] - U[b]).norm2()
                  if dd.is_rational(): byd[Fr(dd.c[0])].append((a, b))
  json.dump({str(k): [list(p) for p in v] for k, v in byd.items()},
            open(CACHE, "w"))
  print(f"  forced pairs by distance: "
      f"{ {str(k): len(v) for k, v in byd.items()} }   [{time.time()-t0:.0f}s]",
      flush=True)

seen2, V = set(U), list(U)
for d2, pairs in sorted(byd.items()):
    spin = rotation_joining(d2, F)
    a, b = pairs[0]
    for w in orbit(U[a]):
        rot = spin.about(w)
        for p in U:
            z = rot(p)
            if z not in seen2: seen2.add(z); V.append(z)
    print(f"  after spindling the d^2={d2} orbit: {len(V)} points"
          f"   [{time.time()-t0:.0f}s]", flush=True)
g = build_graph(V); n = g.n; m = sum(len(a) for a in g.adj) // 2
print(f"double symmetric spindle: n={n} m={m} deg={2.0*m/n:.2f}"
      f"   [{time.time()-t0:.0f}s]", flush=True)
tri = g.find_clique(3)
for KK in (5, 6):
    XX = lambda v, c: 1 + v * KK + c
    cc = [[XX(v, c) for c in range(KK)] for v in range(n)]
    for u, v in g.edges():
        for c in range(KK):
            cc.append([-XX(u, c), -XX(v, c)])
    for i, v in enumerate(tri):
        cc.append([XX(v, i)])
        for c in range(KK):
            if c != i: cc.append([-XX(v, c)])
    t1 = time.time()
    sv = Solver(name="cd19", bootstrap_with=cc); ok = sv.solve(); sv.delete()
    print(f"  {KK}-colourable: {ok}   [{time.time()-t1:.0f}s]", flush=True)
    if not ok:
        print(f"  *** IT REFUSES {KK} COLOURS ***", flush=True)
        json.dump({"field_generators": list(F.gens), "n": n, "m": m,
                   "points": [[[[c.numerator, c.denominator] for c in p.x.c],
                               [[c.numerator, c.denominator] for c in p.y.c]]
                              for p in g.vertices]},
                  open(HN_DIR + "/data/six_candidate.json", "w"))
        print("  written data/six_candidate.json", flush=True)
    else:
        break
