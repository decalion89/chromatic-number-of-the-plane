"""The three-fold union, which is the one that stays symmetric.

rho and the 60-degree rotation share a pivot, so they commute and rho(G*) is
still invariant under the rotation part of the group.  Reflection is the part
that breaks: reflecting rho(G*) gives rho^-1(G*), not itself.  So G* u rho(G*)
is only C6-invariant, and the smallest union that keeps the whole dihedral
group is G* u rho(G*) u rho^-1(G*).

That is not an ornament.  Symmetry is the property that makes the NEXT
rotation bite, and de Grey used exactly this shape at the next level up: G is
Y turned by +theta and -theta about a pivot, not Y turned once.  Y itself is
Sa u Sb alone, and its own dihedral closure Sa u Sb u Sb' is the same
three-fold shape -- 1189 points, which is what the closure of Y actually is.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time, pickle, random
from collections import Counter
from fractions import Fraction as Fr
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
from hn.degrey import build_G
from hn.geometry import DEGREY_FIELD as F, Point, _rot60, rotation_joining
from hn.graph import build_graph
from hn.homcol import closable_distance, decay_rates
from pysat.solvers import Solver

SC = ("/tmp/hn/")
t0 = time.time()
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
print(f"G*: {len(Gs)} points  [{time.time()-t0:.0f}s]", flush=True)

rings = Counter()
for p in Gs:
    dx, dy = p.x - PIV.x, p.y - PIV.y
    rings[float(dx * dx + dy * dy)] += 1
cand = []
for r2f, cnt in rings.items():
    D = Fr(round(r2f * 55440), 55440)
    if abs(float(D) - r2f) > 1e-9 or D < Fr(1, 4) or D == 1:
        continue
    if closable_distance(D):
        cand.append((cnt, D))
cand.sort(reverse=True)
print(f"closable rings: {[(c, str(d)) for c, d in cand]}  "
      f"[{time.time()-t0:.0f}s]", flush=True)


def pairs_of(zf, lim=36.0):
    out, n = [], len(zf)
    for i in range(n):
        ai, bi = zf[i]
        for j in range(i + 1, n):
            v = (ai - zf[j][0]) ** 2 + (bi - zf[j][1]) ** 2
            if v > lim:
                continue
            D = Fr(round(v * 1584), 1584)
            if abs(float(D) - v) > 1e-7 or D == 1 or not closable_distance(D):
                continue
            out.append((i, j))
    return out


for cnt, D in cand[:8]:
    base = rotation_joining(D, F)
    rho, inv = base.about(PIV), base.inverse().about(PIV)
    allp, seen2 = list(Gs), set(Gs)
    for p in Gs:
        for q in (rho(p), inv(p)):
            if q not in seen2:
                seen2.add(q)
                allp.append(q)
    n = len(allp)
    g = build_graph(allp)
    E = sorted((min(a, b), max(a, b)) for a, b in g.edges())
    cls = [[1 + w * 5 + c for c in range(5)] for w in range(n)]
    for a, b in E:
        for c in range(5):
            cls.append([-(1 + a * 5 + c), -(1 + b * 5 + c)])
    sv = Solver(name="cd19", bootstrap_with=cls)
    ok = sv.solve()
    cf = sv.accum_stats().get("conflicts", 0)
    sv.delete()
    print(f"  D={D} ({cnt} on the ring): {n} points, {len(E)} edges, "
          f"5-colourable {ok} ({cf} conflicts)  [{time.time()-t0:.0f}s]",
          flush=True)
    if not ok:
        print("  *** NOT 5-COLOURABLE -- SIX COLOURS ***", flush=True)
        with open(SC + "gstar3_six.pkl", "wb") as fh:
            pickle.dump((str(D), [(str(p.x), str(p.y)) for p in allp], E),
                        fh)
        break
