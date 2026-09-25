"""How small can a carrier be and still have its glue force a pair?

The size of the 5-chromatic graph this construction produces is twice the size
of the carrier whose glue forces.  Sa needs 397 and gives 1139.  A carrier of
120 would give under 250, which is less than half the 509-vertex record.  So
it is worth asking the question directly, from the bottom rather than by
shrinking from the top.

Tried here: the 24-point graph that is tight at four colours (data/tight_four
.json), and Sa peeled to a range of sizes.  For each carrier, every glue whose
overlap is neither 1 nor the whole carrier is built and the forced pairs of the
union counted.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time, json, random
from fractions import Fraction as Fr
from collections import defaultdict, Counter
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
from hn.degrey import build_Sa
from hn.geometry import Point, rotation_joining
from hn.field import Field
from hn.graph import build_graph
from pysat.solvers import Solver

K = 4
def forced_count(pts):
    g = build_graph(pts); n = g.n
    X = lambda v, c: 1 + v * K + c
    cnf = [[X(v, c) for c in range(K)] for v in range(n)]
    for u, v in g.edges():
        for c in range(K):
            cnf.append([-X(u, c), -X(v, c)])
    s = Solver(name="m22", bootstrap_with=cnf)
    if not s.solve(): s.delete(); return None, n
    pos = set(l for l in s.get_model() if l > 0)
    cols = [[next(c for c in range(K) if X(v, c) in pos) for v in range(n)]]
    r = random.Random(5)
    for _ in range(60):
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
    nf = 0
    for vs in buck.values():
        for i, a in enumerate(vs):
            for b in vs[i+1:]:
                s2 = Solver(name="m22", bootstrap_with=cnf)
                d = s2.solve(assumptions=[X(a, 0), X(b, 1)]); s2.delete()
                if not d: nf += 1
    return nf, n

def scan(name, pts, F, rots):
    S = set(pts); n0 = len(pts)
    seenkeys = set(); best = 0
    for rname, r0 in rots:
        for c in pts:
            rot = r0.about(c)
            ov = sum(1 for p in pts if rot(p) in S)
            if ov <= 1 or ov == n0: continue
            key = (rname, ov)
            if key in seenkeys: continue
            seenkeys.add(key)
            seen, U = set(), []
            for p in pts:
                for q in (p, rot(p)):
                    if q not in seen: seen.add(q); U.append(q)
            nf, n = forced_count(U)
            if nf:
                print(f"  {name}: {rname} overlap {ov} -> union {n}: "
                      f"FORCED {nf}  <<<<", flush=True)
                best = max(best, nf)
    print(f"  {name} (n={n0}): best forcing over "
          f"{len(seenkeys)} distinct glues = {best}", flush=True)
    return best

# the 24-point graph that is tight at four
d = json.load(open(HN_DIR + "/data/tight_four.json"))
F24 = Field(tuple(d["field"]["generators"]))
# each entry is 16 rationals of x followed by 16 of y, written as strings
P24 = [Point(F24.element([Fr(t) for t in row[:F24.dim]]),
             F24.element([Fr(t) for t in row[F24.dim:]])) for row in d["points"]]
rots24 = [(f"r^2={r}", rotation_joining(Fr(r), F24))
          for r in ("1", "1/3", "3", "4", "1/4")
          if True]
t0 = time.time()
print(f"the tight 24-point graph: {len(P24)} points", flush=True)
try:
    scan("tight24", P24, F24, rots24)
except Exception as e:
    print(f"  tight24 failed: {e}", flush=True)

# Sa peeled to a range of sizes, lowest degree first (which the thinning
# experiment showed preserves the forcing best)
F = Field((3, 11, 247))
full = build_Sa(F)
gfull = build_graph(full)
order = sorted(range(len(full)), key=lambda v: len(gfull.adj[v]))
rots = [("rot60", rotation_joining(Fr(1), F)),
        ("rot120", rotation_joining(Fr(1, 3), F))]
for keepn in (120, 180, 240, 300, 340, 370):
    pts = [full[i] for i in sorted(order[len(full) - keepn:])]
    g = build_graph(pts)
    print(f"\nSa peeled to {g.n} points, m={sum(len(a) for a in g.adj)//2}"
          f"   [{time.time()-t0:.0f}s]", flush=True)
    scan(f"Sa[{keepn}]", pts, F, rots)
