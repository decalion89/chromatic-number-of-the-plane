"""How close is the five-colour step, measured the way the four-colour one is.

Z = Sa u rho(Sa) at four colours keeps all six antipodal pairs of its ring
through twenty samples, which is what forcing looks like.  The question is
what Z5 = G* u rho(G*) does at five, and the answer is a number rather than a
verdict: how many pairs survive, and how fast the rest die.

Sampling a 27673-point union is the expensive part -- each sample is a proper
5-colouring of 138365 variables -- so the samples are few and the bucketing
does the rest for nothing.
"""
import sys, time
from fractions import Fraction as Fr
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
from hn.degrey import build_G, build_Sa
from hn.geometry import DEGREY_FIELD as F, Point, _rot60, rotation_joining
from hn.graph import build_graph
from hn.homcol import agreeing_pairs, closable_distance
from pysat.solvers import Solver
import random

t0 = time.time()


def ring_anti(pts, pivot, D):
    idx = {p: i for i, p in enumerate(pts)}
    out = []
    for i, p in enumerate(pts):
        if (p.x - pivot.x) ** 2 + (p.y - pivot.y) ** 2 != F.rational(D):
            continue
        q = Point(pivot.x + pivot.x - p.x, pivot.y + pivot.y - p.y)
        j = idx.get(q)
        if j is not None and i < j:
            out.append((i, j))
    return out


def curve(name, pts, pivot, D, k, samples=14):
    g = build_graph(pts)
    E = sorted((min(a, b), max(a, b)) for a, b in g.edges())
    n = len(pts)
    cls = [[1 + w * k + c for c in range(k)] for w in range(n)]
    for a, b in E:
        for c in range(k):
            cls.append([-(1 + a * k + c), -(1 + b * k + c)])
    sv = Solver(name="g4", bootstrap_with=cls)
    if not sv.solve():
        print(f"  {name}: NOT {k}-colourable", flush=True)
        sv.delete()
        return
    anti = ring_anti(pts, pivot, D)
    rng = random.Random(1618)
    cols, surv, hist = [], list(anti), []
    for s in range(samples):
        sv.set_phases([(1 if rng.random() < .5 else -1) * (1 + v)
                       for v in range(n * k)])
        sv.solve()
        m = sv.get_model()
        col = [next(c for c in range(k) if m[w * k + c] > 0)
               for w in range(n)]
        cols.append(col)
        surv = [(i, j) for i, j in surv if col[i] == col[j]]
        hist.append(len(surv))
        print(f"    {name} sample {s+1}: {len(surv)} of {len(anti)} "
              f"antipodal still agree  [{time.time()-t0:.0f}s]", flush=True)
        if not surv and s >= 3:
            break
    allp = agreeing_pairs(cols, colours=k)
    close = [(i, j) for i, j in allp
             if closable_distance(
                 Fr(round(float((pts[i].x - pts[j].x) ** 2
                                + (pts[i].y - pts[j].y) ** 2) * 1584), 1584))]
    sv.delete()
    print(f"  {name}: {n} points, {len(E)} edges, {len(anti)} antipodal, "
          f"curve {hist}, {len(allp)} pairs agree everywhere "
          f"({len(close)} at a closable distance)  [{time.time()-t0:.0f}s]",
          flush=True)


ORIG = Point(F.zero(), F.zero())
Sa = build_Sa(F)
Z4 = list(Sa)
s4 = set(Sa)
rho4 = rotation_joining(4, F)
for p in Sa:
    q = rho4(p)
    if q not in s4:
        s4.add(q)
        Z4.append(q)
curve("Z4 = Sa u rho(Sa) at four", Z4, ORIG, Fr(4), 4)

PIV = Point(F.rational(-2), F.zero())
rot = _rot60(F).about(PIV)
G = build_G(F, as_graph=False)
Gs, seen = [], set()
for refl in (False, True):
    for j in range(6):
        for p in G:
            q = p
            for _ in range(j):
                q = rot(q)
            if refl:
                q = Point(q.x, -q.y)
            if q not in seen:
                seen.add(q)
                Gs.append(q)
for D in (Fr(4), Fr(17, 2)):
    rho = rotation_joining(D, F).about(PIV)
    Z5 = list(Gs)
    s5 = set(Gs)
    for p in Gs:
        q = rho(p)
        if q not in s5:
            s5.add(q)
            Z5.append(q)
    curve(f"Z5 = G* u rho_{D}(G*) at five", Z5, PIV, D, 5)
