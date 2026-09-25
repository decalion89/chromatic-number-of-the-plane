"""Rebuild the intersection universe with the condition the geometry asks for.

The two points a unit from both A and B sit at the midpoint plus or minus
(1/2)*sqrt((4-dd)/dd) times (-dy, dx).  Only that single square root has to
lie in the field.  Every universe in this work computed it as sqrt(4-dd)
divided by sqrt(dd) and required both separately, which is strictly stronger
and silently dropped pairs whose ratio is a square when neither part is.

Measured on G: the old test accepts 4.50 per cent of in-range pairs, the
correct one accepts 8.86 -- so every universe here was built from about half
the points the field actually provides, and "adding points changes nothing"
was tested on the thin half.

Recognising the square root is done by characters rather than by search: the
conjugates of y are plus or minus the roots of the conjugates of t, with the
sign pattern a character of the Galois group, so there are sixteen patterns
to try in this field.  Each gives one linear solve, and a candidate is
accepted only when squaring it returns t exactly -- the numerics propose, the
field decides.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time, pickle
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
sys.path.insert(0, "/tmp/claude-0/-home-user-darwin-50/"
                   "aceaa9ec-f432-5848-a506-39c59179b415/scratchpad")
from fractions import Fraction as Fr
from hn.degrey import build_G, build_Sa, build_Y
from hn.geometry import DEGREY_FIELD as K, Point
from hn.fast import IntBasis, fast_edges_complete
from sqrtK import Roots
from pairsat import pairs_at
from pysat.solvers import Solver

SEED = sys.argv[1] if len(sys.argv) > 1 else "G"
RAD = float(sys.argv[2]) if len(sys.argv) > 2 else 3.0
CAP = int(sys.argv[3]) if len(sys.argv) > 3 else 60000
k = int(sys.argv[4]) if len(sys.argv) > 4 else 5
t0 = time.time()
SC = ("/tmp/claude-0/-home-user-darwin-50/"
      "aceaa9ec-f432-5848-a506-39c59179b415/scratchpad/")
R = Roots(K)
half = K.rational(Fr(1, 2))
four = K.rational(4)
P = list({"G": lambda: build_G(K, as_graph=False),
          "Sa": lambda: build_Sa(K), "Y": lambda: build_Y(K)}[SEED]())
n0 = len(P)
cx = sum(float(p.x) for p in P) / n0
cy = sum(float(p.y) for p in P) / n0
print(f"{SEED}: {n0} points, centroid ({cx:.3f},{cy:.3f}), radius {RAD}"
      f"  [{time.time()-t0:.0f}s]", flush=True)
have = set(P)
fresh = []
for i in range(n0):
    A = P[i]
    for j in range(i + 1, n0):
        B = P[j]
        D = A.dist2(B)
        if not .05 < float(D) < 3.99:
            continue
        s = R.sqrt((four - D) / D)
        if s is None:
            continue
        mx, my = (A.x + B.x) * half, (A.y + B.y) * half
        nx = -(B.y - A.y) * s * half
        ny = (B.x - A.x) * s * half
        for q in (Point(mx + nx, my + ny), Point(mx - nx, my - ny)):
            if q in have:
                continue
            if (float(q.x) - cx) ** 2 + (float(q.y) - cy) ** 2 > RAD * RAD:
                continue
            have.add(q)
            fresh.append(q)
    if len(fresh) > CAP:
        print(f"   cap reached at seed {i}", flush=True)
        break
P = P + fresh
print(f"universe: {len(P)} points (+{len(fresh)})  [{time.time()-t0:.0f}s]",
      flush=True)
pickle.dump(P, open(SC + f"uni2_{SEED}.pkl", "wb"))
b = IntBasis.covering(P)
r = b.rows(P)
hr = b.overflow_headroom(r)
print(f"overflow headroom {hr:.3f}  [{time.time()-t0:.0f}s]", flush=True)
if hr >= 1.0:
    raise SystemExit("overflow")
E1 = set((min(a, c), max(a, c)) for a, c in fast_edges_complete(b, r))
E2 = set(pairs_at(b, r, int(Fr(3) * b.D * b.D)))
n = len(P)
print(f"{n} points, {len(E1)} unit ({2*len(E1)/n:.2f}/v), {len(E2)} at sqrt3"
      f"  [{time.time()-t0:.0f}s]", flush=True)
for name, ed in (("unit", sorted(E1)), ("with sqrt3", sorted(E1 | E2))):
    for kk in (k, k + 1, k + 2):
        cls = [[1 + v * kk + c for c in range(kk)] for v in range(n)]
        for a, c in ed:
            for col in range(kk):
                cls.append([-(1 + a * kk + col), -(1 + c * kk + col)])
        for col in range(1, kk):
            cls.append([-(1 + col)])
        s = Solver(name="cd15", bootstrap_with=cls)
        ok = s.solve()
        s.delete()
        if ok:
            print(f"   {name}: chi = {kk}"
                  + ("   <<< SIX OR MORE" if kk > k else "")
                  + f"  [{time.time()-t0:.0f}s]", flush=True)
            break
        print(f"   {name}: not {kk}-colourable  [{time.time()-t0:.0f}s]",
              flush=True)
print("DONE", flush=True)
