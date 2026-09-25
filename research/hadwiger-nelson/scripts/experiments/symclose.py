"""Symmetrise first, then ask.  The order is the whole point.

Sa is the 12-fold dihedral orbit of S about the origin, and it is exactly
invariant -- 397 points in, 397 out.  That invariance is what makes rotating
it useful: every point stays on its own ring, so a rotation by any angle whose
chord at that radius is 1 produces cross edges in bulk.  Rotating a lopsided
graph instead produces nothing, which is what the earlier scan of G found:
987 rotations of de Grey's field bite G, and the best of them contributes four
cross edges.

So close the graphs under the symmetry BEFORE rotating, and about the centre
each one was actually built around: the origin for Y, whose core Sa lives
there, and the pivot (-2,0) for G, which is the point its two copies of Y were
turned about.  About the origin G's orbit is 6 x 1581 points and 6 x 7877
edges, to the last unit -- six copies that never touch.  About its own pivot
it is 13873 points, five thousand short of twelve disjoint copies.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time, random
from fractions import Fraction as Fr
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
from hn.degrey import build_G, build_Y
from hn.geometry import DEGREY_FIELD as F, Point, _rot60
from hn.graph import build_graph
from hn.homcol import closable_distance
from pysat.solvers import Solver

SC = ("/tmp/claude-0/-home-user-darwin-50/aceaa9ec-f432-5848-a506-"
      "39c59179b415/scratchpad/")
t0 = time.time()
r60 = _rot60(F)
PIV = Point(F.rational(-2), F.zero())


def dihedral(pts, about=None):
    rot = r60.about(about) if about is not None else r60
    out, seen = [], set()
    for refl in (False, True):
        for j in range(6):
            for p in pts:
                q = p
                for _ in range(j):
                    q = rot(q)
                if refl:
                    q = Point(q.x, -q.y)
                if q not in seen:
                    seen.add(q)
                    out.append(q)
    return out


def pairs_of(zf):
    out, n = [], len(zf)
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


def interrogate(name, pts, k):
    g = build_graph(pts)
    E = [(min(a, b), max(a, b)) for a, b in g.edges()]
    n = len(pts)
    cls = [[1 + w * k + c for c in range(k)] for w in range(n)]
    for a, b in E:
        for c in range(k):
            cls.append([-(1 + a * k + c), -(1 + b * k + c)])
    sv = Solver(name="cd19", bootstrap_with=cls)
    ok = sv.solve()
    print(f"  {name}: {n} points, {len(E)} edges, {k}-colourable {ok}"
          f"  [{time.time()-t0:.0f}s]", flush=True)
    if not ok:
        print(f"  *** {name} IS NOT {k}-COLOURABLE ***", flush=True)
        import pickle
        with open(SC + f"symclose_{name.replace(' ', '_')}.pkl", "wb") as fh:
            pickle.dump(([(str(p.x), str(p.y)) for p in pts], E), fh)
        sv.delete()
        return
    sv.delete()
    zf = [(float(p.x), float(p.y)) for p in pts]
    cand = pairs_of(zf)
    sam = Solver(name="g4", bootstrap_with=cls)
    sam.solve()
    rng = random.Random(31337)
    surv, curve, since = cand, [], 0
    for s in range(40):
        sam.set_phases([(1 if rng.random() < .5 else -1) * (1 + i)
                        for i in range(n * k)])
        sam.solve()
        m = sam.get_model()
        col = [next(c for c in range(k) if m[w * k + c] > 0)
               for w in range(n)]
        before = len(surv)
        surv = [(i, j) for i, j in surv if col[i] == col[j]]
        since = 0 if len(surv) < before else since + 1
        if s + 1 in (1, 2, 3, 5, 10, 20, 40):
            curve.append((s + 1, len(surv)))
        if not surv or since >= 15:
            break
    sam.delete()
    chk = Solver(name="cd19", bootstrap_with=cls)
    hits = [(i, j) for i, j in surv
            if not chk.solve(assumptions=[1 + i * k, -(1 + j * k)])]
    chk.delete()
    print(f"    {len(cand)} pairs, decay {curve}, {len(hits)} FORCED"
          f"  [{time.time()-t0:.0f}s]", flush=True)
    if hits:
        import pickle
        with open(SC + f"symclose_forced_{k}.pkl", "wb") as fh:
            pickle.dump(([(str(p.x), str(p.y)) for p in pts], E, hits), fh)
        print(f"  *** FORCED PAIR AT {k} COLOURS in {name} ***", flush=True)


Ystar = dihedral(build_Y(F))
print(f"Y*: {len(Ystar)} points  [{time.time()-t0:.0f}s]", flush=True)
interrogate("Y* at four colours", Ystar, 4)
interrogate("Y* at five colours", Ystar, 5)

Gstar = dihedral(build_G(F, as_graph=False), about=PIV)
print(f"G*: {len(Gstar)} points  [{time.time()-t0:.0f}s]", flush=True)
interrogate("G* at five colours", Gstar, 5)
