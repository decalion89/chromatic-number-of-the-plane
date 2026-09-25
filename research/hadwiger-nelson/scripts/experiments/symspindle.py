"""A 5-chromatic graph with a group of its own.

The chain so far breaks symmetry twice and repairs it never: Sa has isometry
group order 12, the glue rotation is about one vertex, and the spindle is about
one pivot, so the finished graph has group order 1.  Every attempt to glue it
again then tops out at 39% overlap, where Sa reaches 56-85%, and overlap is
what defeats the sigma-argument.

Both breaks can be avoided, because both operations respect conjugation:

    g rot_w g^-1 = rot_{g(w)}        for g in the rotation group

so gluing at EVERY vertex of an orbit at once, and spindling at EVERY pair of a
forced orbit at once, leaves the result invariant.  The forcing is transitive
and equivariant, so a forced pair in an invariant carrier drags its whole orbit
with it and all six spindles are legitimate at the same time.

The result is a 5-chromatic unit-distance graph carrying a real group -- the
first in this project -- and therefore the first one whose own glues can reach
the overlaps that make forcing appear.
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

F = Field((3, 11, 247))
rot60 = _rot60(F)
VI = int(sys.argv[1]) if len(sys.argv) > 1 else 25
GL = sys.argv[2] if len(sys.argv) > 2 else "1"          # glue radius^2
Sa = build_Sa(F)

def orbit(p):
    out, q = [], p
    for _ in range(6):
        if q not in out: out.append(q)
        q = rot60(q)
    return out

def cnf_of(g, K):
    X = lambda v, c: 1 + v * K + c
    out = [[X(v, c) for c in range(K)] for v in range(g.n)]
    for u, v in g.edges():
        for c in range(K):
            out.append([-X(u, c), -X(v, c)])
    return out, X

def colourings(g, K, ncol, seed=7):
    """Minisat with randomised polarity AND blocking.

    Which diversifier works depends on the instance, and the candidate count is
    the diagnostic.  On this carrier -- 1021 vertices, mean degree 13.34, tight
    at four colours -- sixteen colourings from cadical with blocking alone left
    48491 surviving pairs, while minisat with random phases and the same
    blocking left 168, both in well under a second.  On a large loose instance
    the verdict reverses.  So check the count rather than trusting either.
    """
    cnf, X = cnf_of(g, K)
    s = Solver(name="m22", bootstrap_with=cnf)
    if not s.solve(): s.delete(); return [], cnf, X
    pos = set(l for l in s.get_model() if l > 0)
    cols = [[next(c for c in range(K) if X(v, c) in pos) for v in range(g.n)]]
    rng = random.Random(seed); t = 0
    while len(cols) < ncol and t < ncol * 4:
        t += 1
        last = cols[-1]
        s.add_clause([-X(v, last[v]) for v in rng.sample(range(g.n), min(40, g.n))])
        s.set_phases([(1 if rng.random() < 1.0 / K else -1) * X(v, c)
                      for v in range(g.n) for c in range(K)])
        if not s.solve(): break
        p2 = set(l for l in s.get_model() if l > 0)
        cols.append([next(c for c in range(K) if X(v, c) in p2) for v in range(g.n)])
    s.delete()
    return cols, cnf, X


def forced(g, P, cols, cnf, X, cap=6000):
    buck = defaultdict(list)
    for v in range(g.n): buck[tuple(c[v] for c in cols)].append(v)
    cand = [(a, b) for vs in buck.values() if len(vs) > 1
            for i, a in enumerate(vs) for b in vs[i+1:]]
    out = []
    for a, b in cand[:cap]:
        s = Solver(name="m22", bootstrap_with=cnf)
        d = s.solve(assumptions=[X(a, 0), X(b, 1)]); s.delete()
        if not d: out.append((a, b, (P[a] - P[b]).norm2()))
    return len(cand), out

t0 = time.time()
# 1. the symmetric carrier
glue = rotation_joining(Fr(GL), F)
seen, U = set(Sa), list(Sa)
for w in orbit(Sa[VI]):
    rot = glue.about(w)
    for p in Sa:
        q = rot(p)
        if q not in seen: seen.add(q); U.append(q)
gU = build_graph(U); mU = sum(len(a) for a in gU.adj) // 2
print(f"symmetric carrier: n={gU.n} m={mU} deg={2.0*mU/gU.n:.2f}"
      f"   [{time.time()-t0:.0f}s]", flush=True)
cols, cnf, X = colourings(gU, 4, 16)
if not cols:
    print("  *** the carrier already refuses four colours ***", flush=True); sys.exit()
free = min(sum(1 for v in range(gU.n)
               if len({cols[0][u] for u in gU.adj[v]} | {cols[0][v]}) < 4)
           for _ in [0])
ncand, fp = forced(gU, U, cols, cnf, X)
print(f"  free@4={100.0*free/gU.n:.2f}%  candidates {ncand}  FORCED {len(fp)}"
      f"  distances {dict(Counter(str(f[2]) for f in fp).most_common(3))}"
      f"   [{time.time()-t0:.0f}s]", flush=True)
if not fp:
    print("  no forced pair -- nothing to spindle", flush=True); sys.exit()

# 2. spindle the whole orbit of one forced pair at once
a, b, d2 = next((x for x in fp if x[2].is_rational()), (None, None, None))
if a is None:
    print("  every forced pair is at an irrational distance", flush=True); sys.exit()
q = Fr(d2.c[0])
spin = rotation_joining(q, F)
pivots = orbit(U[a])
seen2, V = set(U), list(U)
for w in pivots:
    rot = spin.about(w)
    for p in U:
        z = rot(p)
        if z not in seen2: seen2.add(z); V.append(z)
gV = build_graph(V); mV = sum(len(x) for x in gV.adj) // 2
print(f"  spindled at d^2={q} over {len(pivots)} pivots: n={gV.n} m={mV} "
      f"deg={2.0*mV/gV.n:.2f}   [{time.time()-t0:.0f}s]", flush=True)
c4, _, _ = colourings(gV, 4, 1)
print(f"  4-colourable: {bool(c4)}", flush=True)
c5, cnf5, X5 = colourings(gV, 5, 20)
if not c5:
    print("  *** IT REFUSES FIVE COLOURS ***", flush=True)
    json.dump({"n": gV.n}, open(f"{'/tmp/claude-0/-home-user-darwin-50/aceaa9ec-f432-5848-a506-39c59179b415/scratchpad'}/SIXSYM.json", "w"))
    sys.exit()
free5 = min(sum(1 for v in range(gV.n)
                if len({col[u] for u in gV.adj[v]} | {col[v]}) < 5) for col in c5)
n5, f5 = forced(gV, V, c5, cnf5, X5)
print(f"  free@5={100.0*free5/gV.n:.2f}%  candidates {n5}  FORCED at five: {len(f5)}"
      + ("   <<<<<<<<<<<<<<<<" if f5 else ""), flush=True)
for x in f5[:8]:
    print(f"      v{x[0]} v{x[1]} d^2={x[2]}", flush=True)
json.dump({"field_generators": list(F.gens), "n": gV.n, "from_vertex": VI,
           "glue_r2": GL, "spindle_d2": str(q),
           "points": [[[[c.numerator, c.denominator] for c in p.x.c],
                       [[c.numerator, c.denominator] for c in p.y.c]]
                      for p in gV.vertices]},
          open("/tmp/claude-0/-home-user-darwin-50/aceaa9ec-f432-5848-a506-39c59179b415/scratchpad/symspindle.json", "w"))
print(f"  saved   [{time.time()-t0:.0f}s]", flush=True)
