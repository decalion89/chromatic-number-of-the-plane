"""The strongest local carrier a configuration can have, and its cap.

The lemma search asked G for caps at five colours and got none from 75500
configurations.  G is the WEAKEST carrier available for that question: 1581
points spread over a disc of radius 3, mean degree 10.  The cap is monotone
downward in the carrier -- more points, more constraints, a smaller maximum
-- so the right test builds the densest patch the geometry allows around a
configuration and asks there.

Density comes from a condition this pass derived and never used.  Two points
at squared distance D have unit circles meeting when D < 4, and the meeting
points are

    midpoint  +-  (1/2) sqrt((4-D)/D) * (-dy, dx),

which needs ONE square root of the ratio, not two separate roots of D and
4-D.  Requiring both is strictly stronger and was throwing away about half
the constructible points -- measured at the time: 4.50 per cent of in-range
pairs accepted against 8.86 per cent, a factor of 1.97.

So: seed with the 24-point tight core, close it under the corrected
condition inside a small disc, and run the cap search on the result.  The
control at four colours has to keep working, or the answer at five means
nothing.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time, pickle
from fractions import Fraction as Fr
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
sys.path.insert(0, "/tmp/claude-0/-home-user-darwin-50/"
                   "aceaa9ec-f432-5848-a506-39c59179b415/scratchpad")
from hn.geometry import DEGREY_FIELD as K, Point
from hn.fast import IntBasis, fast_edges_complete
from sqrtK import Roots

SC = ("/tmp/claude-0/-home-user-darwin-50/"
      "aceaa9ec-f432-5848-a506-39c59179b415/scratchpad/")
RAD = float(sys.argv[1]) if len(sys.argv) > 1 else 1.6
CAP = int(sys.argv[2]) if len(sys.argv) > 2 else 2500
ROUNDS = int(sys.argv[3]) if len(sys.argv) > 3 else 3
SEED = sys.argv[4] if len(sys.argv) > 4 else "tightmin.pkl"
OUT = sys.argv[5] if len(sys.argv) > 5 else "densepatch.pkl"
t0 = time.time()
R = Roots(K)
if SEED == "G":
    from hn.degrey import build_G
    base = build_G(K, as_graph=False)
else:
    base = pickle.load(open(SC + SEED, "rb"))
half = K.rational(Fr(1, 2))
four = K.rational(4)
cx = sum(float(p.x) for p in base) / len(base)
cy = sum(float(p.y) for p in base) / len(base)
_root = {}


def ratio_root(D):
    k = tuple(D.c)
    if k not in _root:
        _root[k] = R.sqrt((four - D) / D)
    return _root[k]


def key(p):
    return (tuple(p.x.c), tuple(p.y.c))


pts = {}
for p in base:
    pts[key(p)] = p
print(f"seed {len(pts)} points, disc radius {RAD} about "
      f"({cx:.3f},{cy:.3f})  [{time.time()-t0:.0f}s]", flush=True)
for rnd in range(ROUNDS):
    cur = list(pts.values())
    added = 0
    for i, A in enumerate(cur):
        if len(pts) >= CAP:
            break
        for B in cur[i + 1:]:
            D = A.dist2(B)
            fd = float(D)
            if not (1e-9 < fd < 4.0):
                continue
            s = ratio_root(D)
            if s is None:
                continue
            mx, my = (A.x + B.x) * half, (A.y + B.y) * half
            nx = -(B.y - A.y) * s * half
            ny = (B.x - A.x) * s * half
            for P2 in (Point(mx + nx, my + ny), Point(mx - nx, my - ny)):
                if (float(P2.x) - cx) ** 2 + (float(P2.y) - cy) ** 2 > RAD ** 2:
                    continue
                k = key(P2)
                if k not in pts:
                    pts[k] = P2
                    added += 1
            if len(pts) >= CAP:
                break
    print(f"   round {rnd}: +{added} -> {len(pts)} points, "
          f"{len(_root)} distinct distances  [{time.time()-t0:.0f}s]",
          flush=True)
    if added == 0 or len(pts) >= CAP:
        break
P = list(pts.values())
b = IntBasis.covering(P)
rows = b.rows(P)
hr = b.overflow_headroom(rows)
print(f"overflow headroom {hr:.3f}  [{time.time()-t0:.0f}s]", flush=True)
if hr < 1.0:
    E = sorted(set((min(a, c), max(a, c))
                   for a, c in fast_edges_complete(b, rows)))
    print(f"dense patch: {len(P)} points, {len(E)} edges, mean degree "
          f"{2*len(E)/len(P):.2f}  [{time.time()-t0:.0f}s]", flush=True)
    pickle.dump(P, open(SC + OUT, "wb"))
    print(f"saved to {OUT}", flush=True)
else:
    print("integer basis overflows; patch not usable", flush=True)
print("DONE", flush=True)
