"""Calibrate the decay curve against a graph that actually forces.

Y forces a monochromatic pair in every proper 4-colouring.  That is not a
measurement, it is the reason Y exists, so Y's survivor curve CANNOT reach
zero: the forced pair agrees in every colouring there is.  Whatever floor Y's
curve settles at is the shape of forcing, and every other curve is read
against it.

G alone at five colours, and the control -- two copies of G a thousand apart,
so no cross edge is geometrically possible -- give the other end: graphs with
no forcing anywhere near them.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time, random
from fractions import Fraction as Fr
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
from hn.degrey import build_G, build_Y
from hn.geometry import DEGREY_FIELD as F
from hn.graph import build_graph
from hn.homcol import closable_distance
from pysat.solvers import Solver

MARKS = (1, 2, 3, 5, 10, 20, 40)
random.seed(8191)
t0 = time.time()


def curve_of(name, n, edges, k, cand, cap=60, stall=25):
    cls = [[1 + w * k + c for c in range(k)] for w in range(n)]
    for a, b in edges:
        for c in range(k):
            cls.append([-(1 + a * k + c), -(1 + b * k + c)])
    sv = Solver(name="cd19", bootstrap_with=cls)
    if not sv.solve():
        print(f"  {name}: NOT {k}-colourable", flush=True)
        sv.delete()
        return
    surv, out, since = cand, [], 0
    for s in range(cap):
        for _t in range(6):
            vs = random.sample(range(n), 8)
            if sv.solve(assumptions=[1 + v * k + random.randrange(k)
                                     for v in vs]):
                break
        else:
            sv.solve()
        m = sv.get_model()
        col = [0] * n
        for w in range(n):
            for c in range(k):
                if m[w * k + c] > 0:
                    col[w] = c
                    break
        before = len(surv)
        surv = [(i, j) for i, j in surv if col[i] == col[j]]
        since = 0 if len(surv) < before else since + 1
        if s + 1 in MARKS:
            out.append((s + 1, len(surv)))
        if not surv or since >= stall:
            break
    out.append((s + 1, len(surv)))
    hits = [(i, j) for i, j in surv
            if not sv.solve(assumptions=[1 + i * k, -(1 + j * k)])]
    sv.delete()
    print(f"  {name}: {n} pts, {len(edges)} edges, {k} colours, "
          f"{len(cand)} pairs, decay {out}, {len(hits)} FORCED"
          f"  [{time.time()-t0:.0f}s]", flush=True)


def pairs_of(zf, n):
    out = []
    for i in range(n):
        ai, bi = zf[i]
        for j in range(i + 1, n):
            v = (ai - zf[j][0]) ** 2 + (bi - zf[j][1]) ** 2
            if v > 36.0:
                continue
            D = Fr(round(v * 1584), 1584)
            if abs(float(D) - v) > 1e-7 or D == 1 or not closable_distance(D):
                continue
            out.append((i, j))
    return out


ypts = build_Y(F)
yz = [(float(p.x), float(p.y)) for p in ypts]
ye = [(min(a, b), max(a, b)) for a, b in build_graph(ypts).edges()]
print(f"Y: {len(yz)} vertices, {len(ye)} edges  [{time.time()-t0:.0f}s]",
      flush=True)
curve_of("Y at four colours (forces)", len(yz), ye, 4, pairs_of(yz, len(yz)))
curve_of("Y at five colours", len(yz), ye, 5, pairs_of(yz, len(yz)))

pts = build_G(F, as_graph=False)
gz = [(float(p.x), float(p.y)) for p in pts]
ge = [(min(a, b), max(a, b)) for a, b in build_graph(pts).edges()]
n = len(gz)
print(f"G: {n} points, {len(ge)} edges  [{time.time()-t0:.0f}s]", flush=True)
curve_of("G alone", n, ge, 5, pairs_of(gz, n))

cz = gz + [(x + 1000.0, y2) for x, y2 in gz]
ce = ge + [(a + n, b + n) for a, b in ge]
curve_of("control: two copies, 1000 apart", 2 * n, ce, 5, pairs_of(cz, 2 * n))
