"""How far from blocked, as a number -- and the cheapest target in the plane.

mu_5 = 2 everywhere is a coarse reading: it says a point is placeable, not how
nearly it is not.  The fine reading is the escape count.  Fix p and enumerate
the colourings that 5-colourings of G induce on N(p) while leaving colour 0
unused there -- by colour symmetry exactly "p is still placeable":

    E(p) = # colourings of N(p) induced by 5-colourings of G, avoiding colour 0
    M(p) = # proper colourings of N(p) in 4 colours, the ceiling on E

E(p) = 0 is precisely mu_5 = 5, precisely a blocked point, precisely
chi(R^2) >= 6.  M(p) - E(p) is how many patterns the graph already kills, and
M(p) is the price of admission: no graph can block p without killing all M.

The grading makes M computable with no solving.  N(p) is a disjoint union of
paths of at most six points and hexagons, so

    M = product over components of  4 * 3^(m-1)   for a path of m
                                    3^6 + 3 = 732 for a hexagon

and the floor is |N| = 5, since five colours need five points.  Among the
|N| = 5 shapes:

    one path of five      4 * 3^4            =  324      <- the cheapest
    path of four + point  4*3^3 * 4          =  432
    path of three + pair  4*3^2 * 4*3        =  432
    five isolated points  4^5                = 1024

So the cheapest blocked point in the entire problem is a point whose five
neighbours form a single 60-degree arc, and it costs 324 refutations.  That is
the whole of chi(R^2) >= 6, stated as a finite number.

Enumeration is by solve-and-block, guarded by a per-point selector so the
blocking clauses of one point do not silently constrain the graph for the next,
and bounded by the best seen so far, since only the minimum is wanted.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time, json, math
from fractions import Fraction as Fr
from collections import defaultdict, Counter, deque
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
from hn.field import Field
from hn.geometry import Point
from hn.graph import build_graph
from pysat.solvers import Solver

ROOT = HN_DIR
t0 = time.time(); K = 5
NAME = sys.argv[1] if len(sys.argv) > 1 else "five_247_c.json"
LO, HI = (int(sys.argv[2]), int(sys.argv[3])) if len(sys.argv) > 3 else (5, 6)
CAP = int(sys.argv[4]) if len(sys.argv) > 4 else 1100

d = json.load(open(f"{ROOT}/data/{NAME}"))
F = Field(tuple(d["field_generators"]))
P = [Point(F.element([Fr(a, b) for a, b in x]),
           F.element([Fr(a, b) for a, b in y])) for x, y in d["points"]]
g = build_graph(P); n = g.n
hx = [q.fx for q in g.vertices]; hy = [q.fy for q in g.vertices]
X = lambda v, c: 1 + v * K + c
cnf = [[X(v, c) for c in range(K)] for v in range(n)]
for v in range(n):
    for a in range(K):
        for b in range(a + 1, K):
            cnf.append([-X(v, a), -X(v, b)])
for x, y in g.edges():
    for c in range(K):
        cnf.append([-X(x, c), -X(y, c)])
SEL = n * K + 1
s = Solver(name="cd19", bootstrap_with=cnf)
print(f"{NAME} n={n}  base colouring: {s.solve()}   [{time.time()-t0:.0f}s]", flush=True)

cells = defaultdict(list)
for i in range(n): cells[(int(hx[i] // 2), int(hy[i] // 2))].append(i)
cand = {}
for i in range(n):
    cx, cy = int(hx[i] // 2), int(hy[i] // 2)
    for dx in (-1, 0, 1):
        for dy in (-1, 0, 1):
            for j in cells.get((cx + dx, cy + dy), ()):
                if j <= i: continue
                ex, ey = hx[j] - hx[i], hy[j] - hy[i]
                dd = ex * ex + ey * ey
                if dd <= 1e-12 or dd >= 4.0: continue
                h = math.sqrt(max(0.0, 1.0 - dd / 4.0)); r = math.sqrt(dd)
                mx, my = (hx[i] + hx[j]) / 2, (hy[i] + hy[j]) / 2
                ux, uy = -ey / r, ex / r
                for sg in (+1, -1):
                    cand[(round(mx + sg*h*ux, 9), round(my + sg*h*uy, 9))] = 1

def shape(nb):
    """component sizes, and whether each is a hexagon, with no solving"""
    nadj = {u: [v for v in nb if v != u and
                abs((hx[u]-hx[v])**2 + (hy[u]-hy[v])**2 - 1.0) < 1e-9] for u in nb}
    seen = set(); out = []
    for u in nb:
        if u in seen: continue
        seen.add(u); q = deque([u]); comp = [u]; deg = 0
        while q:
            w = q.popleft(); deg += len(nadj[w])
            for z in nadj[w]:
                if z not in seen: seen.add(z); comp.append(z); q.append(z)
        out.append((len(comp), deg // 2 >= len(comp)))
    return out

def ceiling(sh):
    m = 1
    for size, cyc in sh:
        m *= 732 if cyc else 4 * 3 ** (size - 1)
    return m

targets = []
for (kx, ky) in cand:
    nb = [i for i in range(n) if abs((hx[i]-kx)**2 + (hy[i]-ky)**2 - 1.0) < 1e-9]
    if LO <= len(nb) <= HI:
        sh = shape(nb)
        targets.append((ceiling(sh), kx, ky, nb, sh))
targets.sort(key=lambda t: t[0])
print(f"  {len(targets)} candidate points with {LO} <= |N| <= {HI}; "
      f"cheapest ceilings {[t[0] for t in targets[:8]]}   [{time.time()-t0:.0f}s]",
      flush=True)

def escapes(nb, sel, lim):
    ass = [-X(u, 0) for u in nb] + [-sel]
    cnt = 0
    while cnt < lim:
        if not s.solve(assumptions=ass):
            break
        mod = s.get_model()
        pat = []
        for u in nb:
            for c in range(K):
                if mod[X(u, c) - 1] > 0: pat.append((u, c)); break
        s.add_clause([sel] + [-X(u, c) for u, c in pat])
        cnt += 1
    return cnt

bound = CAP; best = None; hist = Counter()
for t, (M, kx, ky, nb, sh) in enumerate(targets):
    e = escapes(nb, SEL + t, min(bound, M))
    cut = e >= min(bound, M)
    hist[("cut" if cut else e)] += 1
    if not cut and (best is None or e < best[0]):
        best = (e, M, kx, ky, len(nb), sh); bound = max(1, e)
        print(f"    E = {e} of M = {M} at ({kx:.6f},{ky:.6f}) |N|={len(nb)} "
              f"shape={sh}   [{time.time()-t0:.0f}s]", flush=True)
    if e == 0:
        print("    *** BLOCKED POINT ***", flush=True); break
print(f"\n  best: {best}   [{time.time()-t0:.0f}s]", flush=True)
print(f"  histogram of E: {dict(sorted(hist.items(), key=str))}", flush=True)
