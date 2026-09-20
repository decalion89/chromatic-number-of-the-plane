"""Build Y one level up: G* u rho(G*), by the same rule that built Y from Sa.

de Grey's step from Sa to Y is a rotation about the centre of symmetry through
the angle that brings two points of the ring of radius 2 to distance exactly 1
-- sin(theta/2) = 1/4.  Nothing about the number 2 is special except that the
ring of radius 2 is populous in Sa and that 4 is closable: cos theta = 1 -
1/(2D) and sin theta = sqrt(4D-1)/(2D) both have to live in the field.

G* is the analogue of Sa at five colours.  It is 5-chromatic, because it
contains G, and it is exactly invariant under the 12-element dihedral group
about (-2,0), which is the property that makes a rotation of it bite: every
point stays on its own ring, so the whole ring produces cross edges at once.
So enumerate G*'s rings about its own pivot, keep the closable ones, rank them
by population, and turn G* through each in turn.
"""
import sys, time, pickle
from collections import Counter
from fractions import Fraction as Fr
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
from hn.degrey import build_G
from hn.geometry import DEGREY_FIELD as F, Point, _rot60, rotation_joining
from hn.graph import build_graph
from hn.homcol import closable_distance
from pysat.solvers import Solver

SC = ("/tmp/claude-0/-home-user-darwin-50/aceaa9ec-f432-5848-a506-"
      "39c59179b415/scratchpad/")
t0 = time.time()
PIV = Point(F.rational(-2), F.zero())
rot = _rot60(F).about(PIV)

G = build_G(F, as_graph=False)
out, seen = [], set()
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
                out.append(q)
n0 = len(out)
GE = [(min(a, b), max(a, b)) for a, b in build_graph(out).edges()]
print(f"G*: {n0} points, {len(GE)} edges  [{time.time()-t0:.0f}s]", flush=True)

# Ring populations about the pivot, by exact squared radius.  Only rational
# radii admit rotation_joining; irrational ones would need a root inside F.
rings = Counter()
for p in out:
    dx, dy = p.x - PIV.x, p.y - PIV.y
    r2 = dx * dx + dy * dy
    rings[float(r2)] += 1
cand = []
for r2f, cnt in rings.items():
    D = Fr(round(r2f * 55440), 55440)
    if abs(float(D) - r2f) > 1e-9 or D < Fr(1, 4):
        continue
    if closable_distance(D):
        cand.append((cnt, D))
cand.sort(reverse=True)
print(f"{len(rings)} rings, {len(cand)} of them rational and closable; "
      f"top: {[(c, str(d)) for c, d in cand[:10]]}  [{time.time()-t0:.0f}s]",
      flush=True)

zf = [(float(p.x), float(p.y)) for p in out]
for cnt, D in cand[:12]:
    rho = rotation_joining(D, F).about(PIV)
    allp, allz = list(out), list(zf)
    idx = {p: i for i, p in enumerate(out)}
    img = []
    for p in out:
        q = rho(p)
        if q in idx:
            img.append(idx[q])
        else:
            idx[q] = len(allp)
            img.append(len(allp))
            allp.append(q)
            allz.append((float(q.x), float(q.y)))
    E = set(GE)
    for a, b in GE:
        x, y = img[a], img[b]
        E.add((min(x, y), max(x, y)))
    n = len(allp)
    shared = 2 * n0 - n
    g2 = build_graph(allp)
    E2 = set((min(a, b), max(a, b)) for a, b in g2.edges())
    cross = len(E2 - E)
    E = E2
    cls = [[1 + w * 5 + c for c in range(5)] for w in range(n)]
    for a, b in sorted(E):
        for c in range(5):
            cls.append([-(1 + a * 5 + c), -(1 + b * 5 + c)])
    sv = Solver(name="cd19", bootstrap_with=cls)
    ok = sv.solve()
    cf = sv.accum_stats().get("conflicts", 0)
    sv.delete()
    print(f"  ring D={D} ({cnt} points): union {n} points, {len(E)} edges, "
          f"{shared} shared, {cross} CROSS, 5-colourable {ok} ({cf} "
          f"conflicts)  [{time.time()-t0:.0f}s]", flush=True)
    if not ok:
        print("  *** NOT 5-COLOURABLE -- SIX COLOURS ***", flush=True)
        with open(SC + "gstarY_six.pkl", "wb") as fh:
            pickle.dump((str(D), [(str(p.x), str(p.y)) for p in allp],
                         sorted(E)), fh)
        break
