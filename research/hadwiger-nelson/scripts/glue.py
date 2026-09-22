"""Is the cap necessary for a glue to manufacture forcing?

The ladder isolated on de Grey's own construction:

    Sa has a capped circle (6 points, <= 2 of 4 colours) and NO forced-equal
    pair at all.  Glue Sa to a rotated copy, hinged on that capped circle, and
    the union Y suddenly has forced-equal pairs.  Spindle one of those and the
    chromatic number goes up.

The open question is whether the cap is doing the work, or whether any hinge
would do.  Sa offers both kinds of circle: the capped one at r^2=4, and others
at r^2 = 5/9, 7/3, 3, ... that cost no new radical.  So the experiment is to
glue at every circle, test the cap and the forcing separately, and see whether
they move together.

If forcing appears only where a cap is, the cap is the bottleneck and chi >= 6
needs one at five colours.  If forcing appears without a cap, the bottleneck is
somewhere else entirely and the route is open.
"""
import sys, time, random, itertools
from fractions import Fraction as Fr
from collections import defaultdict
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
from hn.degrey import build_Sa, build_Y
from hn.geometry import Point, rotation_joining
from hn.field import Field
from hn.graph import build_graph
from pysat.solvers import Solver

F = Field((3, 5, 7, 11))
ORI = Point(F.zero(), F.zero())
K = int(sys.argv[1]) if len(sys.argv) > 1 else 4
MINR = int(sys.argv[2]) if len(sys.argv) > 2 else 4
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

pts = build_Sa(F)
n0 = len(pts)
g0 = build_graph(pts)
X0 = lambda v, c: 1 + v * K + c
cnf0 = [[X0(v, c) for c in range(K)] for v in range(n0)]
for u, v in g0.edges():
    for c in range(K):
        cnf0.append([-X0(u, c), -X0(v, c)])
NV0 = n0 * K

def capped(R):
    """Can R show K-1 or more colours?  UNSAT means a cap at K-2 or below."""
    cnf = list(cnf0)
    Ycol = lambda c: NV0 + 1 + c
    for c in range(K):
        cnf.append([-Ycol(c)] + [X0(v, c) for v in R])
    # at least K-1 of the K colour-indicators
    for pair in itertools.combinations(range(K), 2):
        cnf.append([Ycol(c) for c in range(K) if c not in pair])
    s = Solver(name="m22", bootstrap_with=cnf)
    ok = s.solve(); s.delete()
    return not ok          # True == capped below K-1

def forcing(pts2):
    g = build_graph(pts2)
    n = g.n
    X = lambda v, c: 1 + v * K + c
    cnf = [[X(v, c) for c in range(K)] for v in range(n)]
    for u, v in g.edges():
        for c in range(K):
            cnf.append([-X(u, c), -X(v, c)])
    s = Solver(name="m22", bootstrap_with=cnf)
    rng = random.Random(7)
    cols, tries = [], 0
    while len(cols) < NCOL and tries < NCOL * 4:
        tries += 1
        s.set_phases([(1 if rng.random() < 1.0 / K else -1) * X(v, c)
                      for v in range(n) for c in range(K)])
        if not s.solve(): continue
        pos = set(l for l in s.get_model() if l > 0)
        cols.append([next(c for c in range(K) if X(v, c) in pos) for v in range(n)])
    s.delete()
    if not cols: return None, n, g
    buck = defaultdict(list)
    for v in range(n):
        buck[tuple(col[v] for col in cols)].append(v)
    return sum(len(b) * (len(b) - 1) // 2 for b in buck.values()), n, g

# every circle of Sa whose spindle exists in this field
cands = {}
centres = [("origin", ORI)] + [(f"v{i}", p) for i, p in enumerate(pts)]
for name, c in centres:
    b = defaultdict(list)
    for i, p in enumerate(pts):
        d2 = (p - c).norm2()
        if d2.is_rational() and d2 != 0:
            b[Fr(d2.c[0])].append(i)
    for d2, R in b.items():
        if len(R) < MINR: continue
        rad = radical(d2)
        if rad is None or rad not in IN_FIELD: continue
        key = (d2, tuple(sorted(R)))
        if key not in cands: cands[key] = (name, c)
print(f"Sa, k={K}: {len(cands)} distinct circles with >= {MINR} points and a legal spindle",
      flush=True)

rows, t0 = [], time.time()
for (d2, R), (name, c) in sorted(cands.items(), key=lambda kv: (-len(kv[0][1]), kv[0][0])):
    cp = capped(list(R))
    rot = rotation_joining(d2, F).about(c)
    seen, pts2 = set(), []
    for p in pts:
        for q in (p, rot(p)):
            if q not in seen:
                seen.add(q); pts2.append(q)
    fp, n, g = forcing(pts2)
    shared = 2 * n0 - n
    rows.append((len(R), str(d2), radical(d2), name, cp, shared, n, fp))
    print(f"  r^2={str(d2):<8} |R|={len(R):<3} rad {radical(d2):<3} centre {name:<6} "
          f"capped={str(cp):<5} shared={shared:<3} n={n:<5} forced-pair candidates={fp}"
          f"   [{time.time()-t0:.0f}s]", flush=True)

print("\n--- does forcing follow the cap? ---", flush=True)
yes = [r for r in rows if r[4]]
no = [r for r in rows if not r[4]]
for tag, grp in (("CAPPED circles", yes), ("uncapped circles", no)):
    if not grp: print(f"  {tag}: none"); continue
    withf = [r for r in grp if r[7]]
    print(f"  {tag}: {len(grp)};  glues that manufacture forcing: {len(withf)}")
    for r in withf[:12]:
        print(f"      r^2={r[1]} |R|={r[0]} centre={r[3]} -> {r[7]} candidate pairs")
