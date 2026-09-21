"""chi of the two-distance graph, pushed up the carriers.

A set whose {1, d} graph needs SIX colours forces a monochromatic pair at
distance d in every 5-colouring of the plane's points, with no cap anywhere:
the five classes are independent for distance 1 by construction, so if they
were independent for d as well they would 5-colour a graph that needs six.

On G, with 1581 points, that number is 5.  On Gc, 11047 points, it is still 5.
The carriers left are bigger and denser, and chromatic number feeds on
density, so this asks them.
"""
import sys, time, pickle
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
sys.path.insert(0, "/tmp/claude-0/-home-user-darwin-50/aceaa9ec-f432-5848-a506-39c59179b415/scratchpad/")
import numpy as np
from fractions import Fraction as Fr
from hn.geometry import DEGREY_FIELD as K
from hn.fast import IntBasis, fast_edges_complete
from pysat.solvers import Solver

PKL = sys.argv[1] if len(sys.argv) > 1 else "symG.pkl"
DSQ = Fr(sys.argv[2]) if len(sys.argv) > 2 else Fr(3)
KK = int(sys.argv[3]) if len(sys.argv) > 3 else 5
t0 = time.time()
SC = ("/tmp/claude-0/-home-user-darwin-50/"
      "aceaa9ec-f432-5848-a506-39c59179b415/scratchpad/")
P = pickle.load(open(SC + PKL, "rb"))
b = IntBasis.covering(P)
r = b.rows(P)
hr = b.overflow_headroom(r)
assert hr < 1.0, hr
E1 = set((min(a, c), max(a, c)) for a, c in fast_edges_complete(b, r))
n = len(P)
dm, D2 = b.dim, b.D * b.D
tval = int(DSQ * D2)
from pairsat import pairs_at
s2 = set(pairs_at(b, r, tval))
both = sorted(E1 | s2)
print(f"{PKL}: {n} points, {len(E1)} unit + {len(s2)} at distance^2 {DSQ} "
      f"= {len(both)} edges  [{time.time()-t0:.0f}s]", flush=True)
for kk in range(KK, 9):
    cls = [[1 + v * kk + c for c in range(kk)] for v in range(n)]
    for a, c in both:
        for col in range(kk):
            cls.append([-(1 + a * kk + col), -(1 + c * kk + col)])
    for col in range(1, kk):
        cls.append([-(1 + col)])
    s = Solver(name="cd15", bootstrap_with=cls)
    ok = s.solve()
    s.delete()
    print(f"   {kk}-colourable: {ok}  [{time.time()-t0:.0f}s]", flush=True)
    if ok:
        mark = ("   <<< chi >= 6: a monochromatic pair at that distance is "
                "forced in every 5-colouring" if kk >= 6 else "")
        print(f"   chi = {kk}{mark}", flush=True)
        break
print("DONE", flush=True)
