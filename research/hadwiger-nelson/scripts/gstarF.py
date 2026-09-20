"""Does G* u rho(G*) force a pair at five colours?

This is de Grey's own step, one level up and with his own angle.  Sa is the
dihedral core at four colours and Y = Sa u rho(Sa) forces a pair; G* is the
dihedral core at five and the question is whether G* u rho(G*) does the same.
rho here is rotation_joining(4) about the pivot -- cos 7/8, sin sqrt(15)/8 --
which is 2*arcsin(1/4), the very rotation that made Y out of Sa.

Colourability alone is not the interesting answer.  Y is 4-colourable and its
value is entirely in the pair it forces, so the union has to be interrogated
for forced pairs whether or not it colours, and at 27000 points that needs the
sampling filter rather than an enumeration.
"""
import sys, time, pickle
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
sys.path.insert(0, "/tmp/claude-0/-home-user-darwin-50/"
                   "aceaa9ec-f432-5848-a506-39c59179b415/scratchpad")
from fractions import Fraction as Fr
from hn.degrey import build_G
from hn.geometry import DEGREY_FIELD as F, Point, _rot60, rotation_joining
from hn.graph import build_graph
from bigfilter import sample_colourings, survivors, forced_among

SC = ("/tmp/claude-0/-home-user-darwin-50/aceaa9ec-f432-5848-a506-"
      "39c59179b415/scratchpad/")
RINGS = [Fr(4), Fr(31, 3), Fr(17, 2), Fr(7), Fr(20, 3), Fr(13, 3), Fr(7, 3),
         Fr(3, 2), Fr(16)]
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

for D in RINGS:
    base = rotation_joining(D, F)
    for tag, turns in (("two-fold", (base.about(PIV),)),
                       ("three-fold", (base.about(PIV),
                                       base.inverse().about(PIV)))):
        allp, seen2 = list(Gs), set(Gs)
        for p in Gs:
            for tr in turns:
                q = tr(p)
                if q not in seen2:
                    seen2.add(q)
                    allp.append(q)
        n = len(allp)
        E = sorted((min(a, b), max(a, b))
                   for a, b in build_graph(allp).edges())
        k = 5
        cls = [[1 + w * k + c for c in range(k)] for w in range(n)]
        for a, b in E:
            for c in range(k):
                cls.append([-(1 + a * k + c), -(1 + b * k + c)])
        cols = sample_colourings(cls, n, k, samples=14, seed=hash(str(D)) % 997)
        if cols is None:
            print(f"  D={D} {tag}: {n} points, {len(E)} edges -- NOT "
                  f"5-COLOURABLE  [{time.time()-t0:.0f}s]", flush=True)
            with open(SC + f"gstarF_six_{D.numerator}_{D.denominator}.pkl",
                      "wb") as fh:
                pickle.dump(([(str(p.x), str(p.y)) for p in allp], E), fh)
            sys.exit(0)
        zf = [(float(p.x), float(p.y)) for p in allp]
        surv = survivors(cols, zf)
        hits = forced_among(cls, k, surv)
        print(f"  D={D} {tag}: {n} points, {len(E)} edges, "
              f"{len(set(map(tuple, cols)))} distinct colourings, "
              f"{len(surv)} survive 14 samples, {len(hits)} FORCED"
              f"  [{time.time()-t0:.0f}s]", flush=True)
        if hits:
            with open(SC + "gstarF_forced.pkl", "wb") as fh:
                pickle.dump((str(D), tag,
                             [(str(p.x), str(p.y)) for p in allp], E, hits),
                            fh)
            print("  *** FORCED PAIR AT FIVE COLOURS -- SPINDLE IT ***",
                  flush=True)
            sys.exit(0)
