"""The remaining rings of Z, with the edges carried over instead of rebuilt.

A union C u rho(C) has three kinds of edge: C's own, the images of C's own
under rho, and the cross edges between the halves.  The first two are known
the moment C's edge list is, so only the cross edges need searching -- and
they are the only ones that need exact arithmetic in the big field, because
everything else is a relabelling.  Rebuilding the whole graph at 41618 points
over a 32-dimensional field costs four minutes; carrying the edges over costs
seconds.

Per-sample reporting matters here.  The loop stops early when every antipodal
pair has been separated, so a run that keeps going is a run where pairs are
SURVIVING, and that is precisely the case worth watching rather than waiting
out in silence.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time, pickle, random
from fractions import Fraction as Fr
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
from hn.degrey import build_G
from hn.field import Field, embed
from hn.geometry import DEGREY_FIELD as K, Point, _rot60, rotation_joining
from hn.graph import build_graph
from pysat.solvers import Solver

SC = ("/tmp/claude-0/-home-user-darwin-50/aceaa9ec-f432-5848-a506-"
      "39c59179b415/scratchpad/")
K17 = Field((3, 5, 7, 11, 17))
ONE17 = K17.rational(1)
t0 = time.time()
PIV = Point(K17.rational(-2), K17.zero())


def turn(pts, E, rot):
    """C u rot(C): points, edges, and the cross edges found by grid."""
    out, idx = list(pts), {p: i for i, p in enumerate(pts)}
    img = []
    for p in pts:
        q = rot(p)
        j = idx.get(q)
        if j is None:
            j = len(out)
            idx[q] = j
            out.append(q)
        img.append(j)
    ed = set(E)
    for a, b in E:
        x, y = img[a], img[b]
        ed.add((min(x, y), max(x, y)))
    zf = [(float(p.x), float(p.y)) for p in out]
    cell = {}
    for i, (a, b) in enumerate(zf):
        cell.setdefault((int(a // 1), int(b // 1)), []).append(i)
    cross = 0
    for j in range(len(pts), len(out)):
        a, b = zf[j]
        cx, cy = int(a // 1), int(b // 1)
        for da in (-1, 0, 1):
            for db in (-1, 0, 1):
                for i in cell.get((cx + da, cy + db), ()):
                    if i >= j:
                        continue
                    if abs((a - zf[i][0]) ** 2 + (b - zf[i][1]) ** 2
                           - 1) > 1e-7:
                        continue
                    dx = out[i].x - out[j].x
                    dy = out[i].y - out[j].y
                    if dx * dx + dy * dy == ONE17 and (i, j) not in ed:
                        ed.add((i, j))
                        cross += 1
    return out, sorted(ed), cross, zf


rot60 = _rot60(K).about(Point(K.rational(-2), K.zero()))
G = build_G(K, as_graph=False)
Gs, seen = [], set()
for refl in (False, True):
    for j in range(6):
        for p in G:
            q = p
            for _ in range(j):
                q = rot60(q)
            if refl:
                q = Point(q.x, -q.y)
            if q not in seen:
                seen.add(q)
                Gs.append(q)
GE = sorted((min(a, b), max(a, b)) for a, b in build_graph(Gs).edges())
print(f"G*: {len(Gs)} points, {len(GE)} edges (built over K)  "
      f"[{time.time()-t0:.0f}s]", flush=True)
Gs = [Point(embed(p.x, K17), embed(p.y, K17)) for p in Gs]

Z, ZE, zc, _ = turn(Gs, GE, rotation_joining(16, K17).about(PIV))
print(f"Z: {len(Z)} points, {len(ZE)} edges (+{zc} cross)  "
      f"[{time.time()-t0:.0f}s]", flush=True)

for D in (Fr(4), Fr(17, 2), Fr(16)):
    U, UE, uc, zf = turn(Z, ZE, rotation_joining(D, K17).about(PIV))
    n = len(U)
    idx = {p: i for i, p in enumerate(U)}
    anti = []
    for i, p in enumerate(U):
        if (p.x - PIV.x) ** 2 + (p.y - PIV.y) ** 2 != K17.rational(D):
            continue
        j = idx.get(Point(PIV.x + PIV.x - p.x, PIV.y + PIV.y - p.y))
        if j is not None and i < j:
            anti.append((i, j))
    print(f"  D={D}: {n} points, {len(UE)} edges (+{uc} cross), "
          f"{len(anti)} antipodal pairs  [{time.time()-t0:.0f}s]", flush=True)
    k = 5
    cls = [[1 + w * k + c for c in range(k)] for w in range(n)]
    for a, b in UE:
        for c in range(k):
            cls.append([-(1 + a * k + c), -(1 + b * k + c)])
    sv = Solver(name="g4", bootstrap_with=cls)
    if not sv.solve():
        print(f"  *** D={D}: NOT 5-COLOURABLE ***", flush=True)
        with open(SC + "iterate2_six.pkl", "wb") as fh:
            pickle.dump((str(D), [(str(p.x), str(p.y)) for p in U], UE), fh)
        sv.delete()
        break
    rng, surv = random.Random(271828), list(anti)
    for s in range(14):
        sv.set_phases([(1 if rng.random() < .5 else -1) * (1 + v)
                       for v in range(n * k)])
        sv.solve()
        m = sv.get_model()
        col = {}
        for i, j in surv:
            for w in (i, j):
                if w not in col:
                    col[w] = next(c for c in range(k) if m[w * k + c] > 0)
        surv = [(i, j) for i, j in surv if col[i] == col[j]]
        print(f"      sample {s+1}: {len(surv)} of {len(anti)} still agree "
              f"[{time.time()-t0:.0f}s]", flush=True)
        if not surv:
            break
    sv.delete()
    if surv:
        print(f"  *** D={D}: {len(surv)} antipodal pairs survive all "
              f"samples -- spindle them ***", flush=True)
        with open(SC + f"iterate2_surv_{D.numerator}_{D.denominator}.pkl",
                  "wb") as fh:
            pickle.dump((str(D), [(str(p.x), str(p.y)) for p in U], UE, surv),
                        fh)
