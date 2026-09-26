"""de Grey's own step, one level up, on the first carrier that can take it.

His construction is two operations alternated.  Close under the order-twelve
group, so that every point sits on a complete ring about the origin; then
rotate the whole thing by the one angle that makes a chosen ring touch its own
image, and take the union.  S becomes Sa becomes Y; Y becomes Ya, Yb, G.

The angle is not free.  For a ring of radius rho the image is a unit away when
2*rho*sin(theta/2) = 1, and theta must live in the field or the union has no
exact edges at all.  Working those out for Ga's rings recovers de Grey's own
numbers: radius 2/3 wants cos = -1/8 and sin = 3*sqrt7/8, which is his
rotation about (-2, 0); radius 2 wants cos = 7/8 and sin = sqrt15/8, which is
his Sb.  Two more survive the field that he had no reason to use: radius
sqrt5/3 with cos = 1/10, and radius sqrt3 with cos = 5/6, the Moser hinge.

So there are exactly four bites available on Ga, and this runs all four.  The
expectation is not high -- a bite sharpens a disjunction, and Ga has no capped
ring for it to sharpen -- but the union is a large, tight, symmetric graph
that nothing in this work has built, and the cost of asking is one solve.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time, pickle
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
from fractions import Fraction as Fr
from hn.geometry import DEGREY_FIELD as K, Point, Rotation
from hn.fast import IntBasis, fast_edges_complete
from hn.homcol import ring_palette_bound
from pysat.solvers import Solver

k = int(sys.argv[1]) if len(sys.argv) > 1 else 5
t0 = time.time()
with open("/tmp/hn/symG.pkl",
          "rb") as f:
    Ga = pickle.load(f)
print(f"Ga: {len(Ga)} points  [{time.time()-t0:.0f}s]", flush=True)

# (name, ring D, cos, sin) -- the four angles the field admits.
BITES = [
    ("D=4/9 (de Grey's pivot turn)", Fr(4, 9),
     K.rational(Fr(-1, 8)), K.sqrt(7) * K.rational(Fr(3, 8))),
    ("D=4 (de Grey's Sb turn)", Fr(4),
     K.rational(Fr(7, 8)), K.sqrt(15) * K.rational(Fr(1, 8))),
    ("D=5/9 (unused by de Grey)", Fr(5, 9),
     K.rational(Fr(1, 10)), K.sqrt(11) * K.rational(Fr(3, 10))),
    ("D=3 (the Moser hinge)", Fr(3),
     K.rational(Fr(5, 6)), K.sqrt(11) * K.rational(Fr(1, 6))),
]
for name, D, c, s in BITES:
    chk = c * c + s * s
    assert chk == K.rational(1), (name, chk)
    rot = Rotation(c, s)
    seen, P = set(), []
    for p in Ga:
        for q in (p, rot(p)):
            if q not in seen:
                seen.add(q)
                P.append(q)
    b = IntBasis.covering(P)
    r = b.rows(P)
    hr = b.overflow_headroom(r)
    if hr >= 1.0:
        print(f"{name}: SKIPPED, overflow headroom {hr:.2f}"
              f"  [{time.time()-t0:.0f}s]", flush=True)
        continue
    E = sorted(set((min(a, cc), max(a, cc))
                   for a, cc in fast_edges_complete(b, r)))
    n = len(P)
    cls = [[1 + v * k + col for col in range(k)] for v in range(n)]
    for a, cc in E:
        for col in range(k):
            cls.append([-(1 + a * k + col), -(1 + cc * k + col)])
    sym = list(cls)
    for col in range(1, k):
        sym.append([-(1 + col)])
    sv = Solver(name="cd15", bootstrap_with=sym)
    ok = sv.solve()
    sv.delete()
    print(f"{name}: {n} points, {len(E)} edges ({2*len(E)/n:.2f}/v), "
          f"headroom {hr:.3f} -> "
          f"{'%d-colourable' % k if ok else '*** NOT %d-COLOURABLE ***' % k}"
          f"  [{time.time()-t0:.0f}s]", flush=True)
    if not ok:
        with open(f"/tmp/hn/"
                  f"WITNESS_{str(D).replace('/','_')}.pkl", "wb") as f:
            pickle.dump(P, f)
        print("   witness written", flush=True)
print("DONE", flush=True)
