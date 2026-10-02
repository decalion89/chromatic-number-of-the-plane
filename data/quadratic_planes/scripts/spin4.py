"""spin4.py d STATE.json OUT_PREFIX: spin.py at K colours (env KCOL, default 4): the spindle of a forced pair.
STATE is a growforce.py state with status 'forced': a unit-distance graph G over Q(sqrt d) (points [a,b,c,e]/D, unit
vectors U) in which the origin and m have the same colour in every K-colouring, and 4|m|^2 - 1 = t^2 with t in
Q(sqrt d). With u = (t + i)/(t - i) (|u| = 1), the copy uG rotated about the origin contains u m, and
|m - u m| = |m| |1 - u| = 1. So G + uG + the edge {m, u m} has no K-colouring. This script builds it exactly, checks
the new edge and the rotated edges, confirms with kissat, shrinks it to a vertex-critical graph with one incremental
solver (cores of failed assumptions, as min3inc.py), lists all unit pairs of the result and writes OUT_PREFIX.json in
the certify_q.py format. If every point of the result lies in the rotated copy, it also writes OUT_PREFIX_back.json:
the same graph rotated back by conj(u), with the denominator D of G.
G must be K-colourable, or the pair says nothing; with PLAIN=1 a graph G with no K-colouring is accepted, and the union
is only a larger graph to shrink (at 3 colours, that is how the graph for d = 11 was found: its 76 points lie in the rotated copy)."""
import sys, json, subprocess, time, os
from fractions import Fraction as Fr
from math import isqrt, lcm
from pysat.solvers import Solver
d = int(sys.argv[1]); st = json.load(open(sys.argv[2])); OUT = sys.argv[3]
K = int(os.environ.get('KCOL', '4'))
assert st["status"] == "forced"
D = st["D"]; U = [tuple(u) for u in st["U"]]; P = [tuple(p) for p in st["P"]]; m = tuple(st["m"])
i0, im = st["i0"], st["im"]
assert P[i0] == (0, 0, 0, 0) and P[im] == m
T0 = time.time()
def log(*a): print(f"[{time.time() - T0:6.0f}s]", *a, flush=True)
# exact arithmetic in K = Q(sqrt d): pairs (x, y) = x + y sqrt d with Fractions
def kmul(a, b): return (a[0] * b[0] + d * a[1] * b[1], a[0] * b[1] + a[1] * b[0])
def kadd(a, b): return (a[0] + b[0], a[1] + b[1])
def ksub(a, b): return (a[0] - b[0], a[1] - b[1])
def kinv(a):
    nrm = a[0] * a[0] - d * a[1] * a[1]
    return (a[0] / nrm, -a[1] / nrm)
def ksqrt_int(P0, Q0):
    N = P0 * P0 - d * Q0 * Q0
    if N < 0: return None
    n = isqrt(N)
    if n * n != N: return None
    for s in (n, -n):
        if (P0 + s) % 2: continue
        p2 = (P0 + s) // 2
        if p2 < 0: continue
        p = isqrt(p2)
        if p * p != p2 or p == 0: continue
        if Q0 % (2 * p) == 0:
            q = Q0 // (2 * p)
            if p * p + d * q * q == P0: return (p, q)
    if Q0 == 0 and P0 % d == 0:
        q = isqrt(P0 // d)
        if d * q * q == P0: return (0, q)
    return None
a, b, c, e = m
X = a * a + d * b * b + c * c + d * e * e; Y = 2 * (a * b + c * e)
r = ksqrt_int(4 * X - D * D, 4 * Y)
assert r is not None, "m is not at a spindle distance"
t = (Fr(r[0], D), Fr(r[1], D))
t2 = kmul(t, t); one = (Fr(1), Fr(0))
den = kinv(kadd(t2, one))
alpha = kmul(ksub(t2, one), den); beta = kmul((2 * t[0], 2 * t[1]), den)      # u = alpha + i beta
assert kadd(kmul(alpha, alpha), kmul(beta, beta)) == one
log(f"|m|^2 = ({X} + {Y} sqrt{d})/{D}^2, t = {t[0]} + {t[1]} sqrt{d}; u = ({alpha[0]} + {alpha[1]} s) + i ({beta[0]} + {beta[1]} s)")
def as_K(p):  # point -> (x, y) in K
    return ((Fr(p[0], D), Fr(p[1], D)), (Fr(p[2], D), Fr(p[3], D)))
def rot(p):
    x, y = as_K(p)
    return (ksub(kmul(alpha, x), kmul(beta, y)), kadd(kmul(alpha, y), kmul(beta, x)))
G = [as_K(p) for p in P]
R = [rot(p) for p in P]
D2 = D
for x, y in R:
    for z in (x, y):
        D2 = lcm(D2, z[0].denominator, z[1].denominator)
def to_int(pt):
    (x0, x1), (y0, y1) = pt
    v = (x0 * D2, x1 * D2, y0 * D2, y1 * D2)
    assert all(w.denominator == 1 for w in v)
    return tuple(int(w) for w in v)
pts = [to_int(p) for p in G]
idx = {p: i for i, p in enumerate(pts)}
rmap = []
for p in R:
    q = to_int(p)
    if q not in idx:
        idx[q] = len(pts); pts.append(q)
    rmap.append(idx[q])
log(f"common denominator {D2}; {len(P)} points, rotated copy shares {sum(1 for j in rmap if j < len(P))}; union {len(pts)}")
Pidx = {p: i for i, p in enumerate(P)}
EG = set()
for i, p in enumerate(P):
    for uu in U:
        j = Pidx.get((p[0] + uu[0], p[1] + uu[1], p[2] + uu[2], p[3] + uu[3]))
        if j is not None and i < j:
            EG.add((i, j))
def sql(i, j):
    a_, b_, c_, e_ = (pts[j][k] - pts[i][k] for k in range(4))
    return (a_ * a_ + d * b_ * b_ + c_ * c_ + d * e_ * e_, 2 * (a_ * b_ + c_ * e_))
unit = (D2 * D2, 0)
E = set(EG)
for i, j in EG:
    a_, b_ = rmap[i], rmap[j]
    E.add((min(a_, b_), max(a_, b_)))
spin = (min(im, rmap[im]), max(im, rmap[im]))
assert sql(*spin) == unit, "the spindle edge is not of length 1"
E.add(spin)
assert all(sql(i, j) == unit for i, j in E)
E = sorted(E); n = len(pts)
log(f"{n} vertices, {len(E)} edges (all of length 1, checked exactly); the spindle edge {spin}")
# shrink with one incremental solver: selector per vertex
adj = [set() for _ in range(n)]
for i, j in E:
    adj[i].add(j); adj[j].add(i)
x = lambda v, k: K * v + k + 1
sel = lambda v: K * n + v + 1
s = Solver(name="cadical153")
for v in range(n):
    s.add_clause([-sel(v)] + [x(v, k) for k in range(K)])
for i, j in E:
    for k in range(K):
        s.add_clause([-x(i, k), -x(j, k)])
def core3(S):
    S = set(S); stack = [v for v in S if len(adj[v] & S) < K]
    while stack:
        v = stack.pop()
        if v not in S: continue
        S.discard(v)
        stack += [w for w in adj[v] if w in S and len(adj[w] & S) < K]
    return S
def refute(T):
    if s.solve(assumptions=[sel(v) for v in sorted(T)]): return None
    return {l - K * n - 1 for l in s.get_core() if l > K * n}
def colourable_G():      # G must be 3-colourable, or the pair says nothing (a separate solver: the shrink is unchanged)
    with Solver(name="cadical153") as s1:
        for v in range(len(P)):
            s1.add_clause([x(v, k) for k in range(K)])
        for i, j in EG:
            for k in range(K):
                s1.add_clause([-x(i, k), -x(j, k)])
        return s1.solve()
if not colourable_G():
    if os.environ.get("PLAIN") != "1":
        log(f"G itself has no {K}-colouring: the pair is not forced, and the spindle is not needed"); sys.exit(1)
    log(f"G itself has no {K}-colouring (PLAIN=1: shrinking the union anyway)")
S = core3(range(n))
C = refute(S)
assert C is not None, f"the spindle is {K}-colourable: the pair was not forced"
S0 = core3(C)
log(f"not {K}-colourable (CaDiCaL); first core {len(S0)} vertices")
import random
RUNS = int(os.environ.get("RUNS", "1"))
best = None
for run in range(RUNS):
    S = set(S0)
    order = sorted(S, key=lambda v: (len(adj[v] & S), v))
    if run:
        random.Random(run).shuffle(order)
    for v in order:
        if v not in S: continue
        T = core3(S - {v})
        C = refute(T)
        if C is not None:
            S = core3(C)
    if best is None or len(S) < len(best):
        best = S
    if RUNS > 1:
        log(f"run {run}: {len(S)} vertices (best {len(best)})")
S = best
log(f"vertex-critical: {len(S)} vertices")
L = sorted(S); pos = {v: i for i, v in enumerate(L)}
crit_pts = [list(pts[v]) for v in L]
allE = [[i, j] for i in range(len(L)) for j in range(i + 1, len(L)) if sql(L[i], L[j]) == unit]
json.dump({"D": D2, "U": [], "points": crit_pts, "edges": allE, "spindle": {"m": list(m), "D_m": D, "t": [str(t[0]), str(t[1])]}},
          open(OUT + ".json", "w"))
log(f"wrote {OUT}.json: {len(L)} points, {len(allE)} unit pairs (all of them)")
pre = {j: i for i, j in enumerate(rmap)}
if all(v in pre for v in L):       # all in uG: rotate back by conj(u) into G (a rotation keeps every unit pair)
    json.dump({"D": D, "U": [], "points": [list(P[pre[v]]) for v in L], "edges": allE}, open(OUT + "_back.json", "w"))
    log(f"all {len(L)} points lie in the rotated copy; rotated back (denominator {D}): wrote {OUT}_back.json")
