"""rotunion.py d STATE.json OUT_PREFIX t0 t1 [t0 t1 ...]: the union of a unit-distance graph G over Q(sqrt d) (a growth
state: points [a,b,c,e]/D, unit vectors U; it may well be 3-colourable) with its copies u G rotated about the origin,
one for each t = t0 + t1 sqrt d given (u = (t + i)/(t - i), |u| = 1; t0, t1 fractions). Edges: those of G and of each
copy (all of length 1); the copies meet G in the points they share. If the union is not 3-colourable (CaDiCaL), it is
shrunk over RUNS random deletion orders (default 100) to a vertex-critical graph; all its unit pairs are listed and
OUT_PREFIX.json is written in the certify_q.py format."""
import sys, json, time, os, random
from fractions import Fraction as Fr
from math import lcm
from pysat.solvers import Solver
d = int(sys.argv[1]); st = json.load(open(sys.argv[2])); OUT = sys.argv[3]
ts = [(Fr(sys.argv[i]), Fr(sys.argv[i + 1])) for i in range(4, len(sys.argv), 2)]
D = st["D"]; U = [tuple(u) for u in st["U"]]; P = [tuple(p) for p in st["P"]]
T0 = time.time()
def log(*a): print(f"[{time.time() - T0:6.0f}s]", *a, flush=True)
def kmul(a, b): return (a[0] * b[0] + d * a[1] * b[1], a[0] * b[1] + a[1] * b[0])
def kadd(a, b): return (a[0] + b[0], a[1] + b[1])
def ksub(a, b): return (a[0] - b[0], a[1] - b[1])
def kinv(a):
    n = a[0] * a[0] - d * a[1] * a[1]; return (a[0] / n, -a[1] / n)
one = (Fr(1), Fr(0))
rots = []
for t in ts:
    t2 = kmul(t, t); den = kinv(kadd(t2, one))
    al = kmul(ksub(t2, one), den); be = kmul((2 * t[0], 2 * t[1]), den)
    assert kadd(kmul(al, al), kmul(be, be)) == one
    rots.append((al, be))
G = [((Fr(p[0], D), Fr(p[1], D)), (Fr(p[2], D), Fr(p[3], D))) for p in P]
copies = [G] + [[(ksub(kmul(al, x), kmul(be, y)), kadd(kmul(al, y), kmul(be, x))) for x, y in G] for al, be in rots]
D2 = D
for C in copies[1:]:
    for x, y in C:
        D2 = lcm(D2, x[0].denominator, x[1].denominator, y[0].denominator, y[1].denominator)
def to_int(pt):
    (x0, x1), (y0, y1) = pt
    return tuple(int(w * D2) for w in (x0, x1, y0, y1))
pts, idx, maps = [], {}, []
for C in copies:
    mp = []
    for p in C:
        q = to_int(p)
        if q not in idx:
            idx[q] = len(pts); pts.append(q)
        mp.append(idx[q])
    maps.append(mp)
Pidx = {p: i for i, p in enumerate(P)}
EG = []
for i, p in enumerate(P):
    for uu in U:
        j = Pidx.get((p[0] + uu[0], p[1] + uu[1], p[2] + uu[2], p[3] + uu[3]))
        if j is not None and i < j: EG.append((i, j))
E = set()
for mp in maps:
    for i, j in EG:
        a, b = mp[i], mp[j]; E.add((min(a, b), max(a, b)))
E = sorted(E); n = len(pts)
shared = [len(set(mp) & set(maps[0])) for mp in maps[1:]]
log(f"G: {len(P)} points; {len(rots)} rotated copies sharing {shared} points with G; union {n} vertices, {len(E)} edges; D' = {D2}")
def sql(i, j):
    a, b, c, e = (pts[j][k] - pts[i][k] for k in range(4))
    return (a * a + d * b * b + c * c + d * e * e, 2 * (a * b + c * e))
assert all(sql(i, j) == (D2 * D2, 0) for i, j in random.Random(1).sample(E, min(2000, len(E))))
adj = [set() for _ in range(n)]
for i, j in E: adj[i].add(j); adj[j].add(i)
x = lambda v, k: 3 * v + k + 1
sel = lambda v: 3 * n + v + 1
s = Solver(name="cadical153")
for v in range(n): s.add_clause([-sel(v), x(v, 0), x(v, 1), x(v, 2)])
for i, j in E:
    for k in range(3): s.add_clause([-x(i, k), -x(j, k)])
def core3(S):
    S = set(S); stack = [v for v in S if len(adj[v] & S) < 3]
    while stack:
        v = stack.pop()
        if v not in S: continue
        S.discard(v); stack += [w for w in adj[v] if w in S and len(adj[w] & S) < 3]
    return S
def refute(T):
    if s.solve(assumptions=[sel(v) for v in sorted(T)]): return None
    return {l - 3 * n - 1 for l in s.get_core() if l > 3 * n}
S = core3(range(n))
C = refute(S)
if C is None:
    log("the union is 3-colourable"); sys.exit(0)
S0 = core3(C); log(f"NOT 3-colourable; first core {len(S0)} vertices")
RUNS = int(os.environ.get("RUNS", "100")); best = None
for run in range(RUNS):
    S = set(S0); order = sorted(S, key=lambda v: (len(adj[v] & S), v))
    if run: random.Random(run).shuffle(order)
    for v in order:
        if v not in S: continue
        C = refute(core3(S - {v}))
        if C is not None: S = core3(C)
    if best is None or len(S) < len(best): best = S
log(f"vertex-critical: {len(best)} vertices (best of {RUNS} orders)")
L = sorted(best)
allE = [[i, j] for i in range(len(L)) for j in range(i + 1, len(L)) if sql(L[i], L[j]) == (D2 * D2, 0)]
json.dump({"D": D2, "U": [], "points": [list(pts[v]) for v in L], "edges": allE,
           "rotations": [[str(t[0]), str(t[1])] for t in ts]}, open(OUT + ".json", "w"))
log(f"wrote {OUT}.json: {len(L)} points, {len(allE)} unit pairs")
