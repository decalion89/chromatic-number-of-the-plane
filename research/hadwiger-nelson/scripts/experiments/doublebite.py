"""Two bites, the way de Grey chains them -- because his first one does not
raise the chromatic number.

Sa is 4-chromatic.  Y = Sa u Sb is 4-chromatic.  G = Ya u Yb is 5-chromatic.
The rise happens on the SECOND bite, not the first, and the two are not
independent: the first turns Sa about the origin by cosine 7/8, the bite on
the ring of radius 2; the second turns Y about (-2, 0) by a relative cosine
31/32, the bite on the ring of radius 4.  And (-2, 0) is not an arbitrary new
centre.  It is a point ON the first bitten ring -- the origin's radius-2 orbit
in Sa contains it -- and the second radius is twice the first.

So the step is: bite a ring of radius rho about c, move the centre to a point
of that ring, bite a ring of radius 2*rho about it.  Applied to Gp, which is
G closed under the order-twelve group about its own pivot and the first
carrier here that is tight and symmetric at once, that is de Grey's whole
construction run one level up.

Both angles have to exist in the field, which for radius 2 and radius 4 they
do: sqrt15/8 and 3*sqrt7/32, his own two numbers.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time, pickle
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
import numpy as np
from fractions import Fraction as Fr
from hn.geometry import DEGREY_FIELD as K, Point, Rotation
from hn.fast import IntBasis, fast_edges_complete
from pysat.solvers import Solver

k = int(sys.argv[1]) if len(sys.argv) > 1 else 5
D1 = Fr(sys.argv[2]) if len(sys.argv) > 2 else Fr(4)
D2 = Fr(sys.argv[3]) if len(sys.argv) > 3 else Fr(16)
t0 = time.time()
CLASSES = [1, 3, 5, 7, 11, 15, 21, 33, 35, 55, 77, 105, 165, 231, 385, 1155]
SC = ("/tmp/hn/")


def turn_for(D):
    ct = 1 - Fr(1, 2) / D
    s2 = 1 - ct * ct
    for rr in CLASSES:
        q = s2 / rr
        a, c = q.numerator, q.denominator
        ra, rc = int(round(a ** .5)), int(round(c ** .5))
        if ra * ra == a and rc * rc == c:
            s = K.rational(Fr(ra, rc))
            if rr != 1:
                s = K.sqrt(rr) * s
            cc = K.rational(ct)
            assert cc * cc + s * s == K.rational(1)
            return cc, s, f"cos={ct}, sin={Fr(ra,rc)}*sqrt{rr}"
    raise SystemExit(f"no bite for D={D} in K")


def grow(P, centre, D, label):
    c, s, txt = turn_for(D)
    on = [p for p in P
          if (p.x - centre.x) * (p.x - centre.x)
          + (p.y - centre.y) * (p.y - centre.y) == K.rational(D)]
    print(f"{label}: ring D={D} about ({float(centre.x):.3f},"
          f"{float(centre.y):.3f}) holds {len(on)} points; {txt}"
          f"  [{time.time()-t0:.0f}s]", flush=True)
    if not on:
        raise SystemExit(f"{label}: that ring is empty, nothing to bite")
    t = Rotation(c, s).about(centre)
    seen, out = set(), []
    for p in P:
        for q in (p, t(p)):
            if q not in seen:
                seen.add(q)
                out.append(q)
    return out, on


def report(P, label, solve=True):
    b = IntBasis.covering(P)
    r = b.rows(P)
    hr = b.overflow_headroom(r)
    if hr >= 1.0:
        print(f"{label}: ABORTED, overflow headroom {hr:.2f}", flush=True)
        raise SystemExit(1)
    E = sorted(set((min(a, c), max(a, c))
                   for a, c in fast_edges_complete(b, r)))
    n = len(P)
    print(f"{label}: {n} points, {len(E)} edges ({2*len(E)/n:.2f}/v), "
          f"headroom {hr:.3f}  [{time.time()-t0:.0f}s]", flush=True)
    if not solve:
        return E
    cls = [[1 + v * k + col for col in range(k)] for v in range(n)]
    for a, c in E:
        for col in range(k):
            cls.append([-(1 + a * k + col), -(1 + c * k + col)])
    for col in range(1, k):
        cls.append([-(1 + col)])
    sv = Solver(name="cd15", bootstrap_with=cls)
    ok = sv.solve()
    sv.delete()
    if ok:
        print(f"{label}: {k}-colourable  [{time.time()-t0:.0f}s]", flush=True)
    else:
        print(f"*** {label}: NOT {k}-COLOURABLE -- chi > {k} ***"
              f"  [{time.time()-t0:.0f}s]", flush=True)
        with open(SC + f"WITNESS_double_{D1}_{D2}.pkl".replace("/", "_"),
                  "wb") as f:
            pickle.dump(P, f)
        print("   witness written", flush=True)
    return E


Gp = pickle.load(open(SC + "pivG.pkl", "rb"))
PIV = Point(K.rational(-2), K.zero())
print(f"Gp: {len(Gp)} points  [{time.time()-t0:.0f}s]", flush=True)
M, ring1 = grow(Gp, PIV, D1, "first bite")
report(M, f"after bite 1 (D={D1})")
# the new centre is a point OF the first bitten ring, as (-2,0) is for Sa
C1 = ring1[0]
N, _ = grow(M, C1, D2, "second bite")
report(N, f"after bite 2 (D={D1} then D={D2})")
print("DONE", flush=True)
