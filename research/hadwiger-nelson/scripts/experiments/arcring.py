"""The arc bound of de Grey's centre and ring, inside the carrier.

Sampling 400 assignments to see how wide the ring spreads is the expensive
way to ask an extremal question.  What matters is the bound: the largest,
over ALL homomorphisms, of the shortest arc containing the centre and its
ring.  Small means a cap in de Grey's sense and a device for the spindle to
use; large means his lemma has no positional analogue even with the carrier
helping.

One SAT call per candidate width, bisected.  The image fits in an arc of L
exactly when the circle has a run of p-L positions the set misses, so
asserting that NO such run exists asks for an assignment spreading wider
than L; UNSAT at L means every assignment is confined to L.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time, math, pickle
from fractions import Fraction as Fr
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
from hn.degrey import build_G, build_Sa, build_Y
from hn.geometry import DEGREY_FIELD as K
from hn.fast import IntBasis, fast_edges_complete
from pysat.solvers import Solver

SC = ("/tmp/hn/")
CAR = sys.argv[1] if len(sys.argv) > 1 else "G"
R = Fr(*map(int, (sys.argv[2] if len(sys.argv) > 2 else "9/2").split("/")))
RAD2 = float(sys.argv[3]) if len(sys.argv) > 3 else 4.0
t0 = time.time()
P = ({"G": lambda k: build_G(k, as_graph=False), "Sa": build_Sa,
      "Y": build_Y}[CAR](K) if CAR in ("G", "Sa", "Y")
     else pickle.load(open(SC + CAR, "rb")))
b = IntBasis.covering(P)
rows = b.rows(P)
E = sorted(set((min(a, c), max(a, c)) for a, c in fast_edges_complete(b, rows)))
n, p, q = len(P), R.numerator, R.denominator
xy = [(float(t.x), float(t.y)) for t in P]
cx, cy = xy[0]
ring = [v for v in range(n)
        if abs((xy[v][0] - cx) ** 2 + (xy[v][1] - cy) ** 2 - RAD2) < 1e-9]
S = [0] + ring
print(f"{CAR}: {n} pts; centre + {len(ring)} points at distance "
      f"{math.sqrt(RAD2):.3f}, ratio {R}  [{time.time()-t0:.0f}s]", flush=True)
base = [[1 + v * p + j for j in range(p)] for v in range(n)]
for a, c in E:
    for j in range(p):
        for d in range(-(q - 1), q):
            base.append([-(1 + a * p + j), -(1 + c * p + (j + d) % p)])
USED = n * p + 1
anchor = [[-(USED + j)] + [1 + v * p + j for v in S] for j in range(p)]


def spreads_wider_than(L):
    gap = p - L
    if gap <= 0:
        return False
    cls = base + anchor + [[USED + (i + t) % p for t in range(gap)]
                           for i in range(p)]
    s = Solver(name="cd15", bootstrap_with=cls)
    ok = s.solve()
    s.delete()
    return ok


lo, hi = 1, p
while lo < hi:
    mid = (lo + hi) // 2
    wide = spreads_wider_than(mid)
    print(f"   can it spread wider than {mid}? "
          f"{'yes' if wide else 'NO'}  [{time.time()-t0:.0f}s]", flush=True)
    if wide:
        lo = mid + 1
    else:
        hi = mid
print(f"\narc bound of centre+ring in {CAR}: {lo} of {p}"
      f"   (a cap below 2q = {2*q} would force independence, and the spindle "
      f"wants pairs inside {q//2})  [{time.time()-t0:.0f}s]", flush=True)
print("DONE", flush=True)
