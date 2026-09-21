"""Rings of every radius, not just the ones that fit inside a unit distance.

Every ring scan in this work filtered to distance squared at most 4, on the
reflex that a ring wider than that cannot have two adjacent points on it.  For
a palette bound that reflex is harmless.  For a BITE it is fatal, because a
bite does not need points of the ring to touch each other -- it needs each
point to touch its own image under the turn, and that works at any radius:
2*rho*sin(theta/2) = 1 has a solution for every rho >= 1/2.

de Grey's own construction proves the point.  Sb turns Sa by 2*arcsin(1/4),
cosine 7/8, which is the bite on the ring of radius 2 -- inside the filter.
But G turns Y by a relative 2*arcsin(1/8), cosine 31/32, and that is the bite
on the ring of radius FOUR, distance squared 16, thrown away by every scan
here.  The radius doubles from one level to the next, and the filter kept only
the first level.

The turn also has to exist in the field: sin theta = sqrt(4*rho^2 - 1)/(2*rho^2),
so 4*rho^2 - 1 must be a square times one of the radicands K carries.  rho = 2
gives 15 and rho = 4 gives 63 = 9*7, which is exactly why de Grey's two angles
are available.  rho = 8, the obvious third level, gives 255 = 3*5*17 and is
not available -- but rho = 3 gives 35 and rho = 5 gives 99 = 9*11, and both
are, and neither has been used.

This lists every ring at every radius, and which of them the field can bite.
"""
import sys, time, pickle
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
import numpy as np
from fractions import Fraction as Fr
from collections import defaultdict
from hn.geometry import DEGREY_FIELD as K, Point
from hn.fast import IntBasis, fast_edges_complete

PKL = sys.argv[1] if len(sys.argv) > 1 else "pivG.pkl"
CX = Fr(sys.argv[2]) if len(sys.argv) > 2 else Fr(-2)
CY = Fr(sys.argv[3]) if len(sys.argv) > 3 else Fr(0)
MINPTS = int(sys.argv[4]) if len(sys.argv) > 4 else 5
t0 = time.time()
CLASSES = [1, 3, 5, 7, 11, 15, 21, 33, 35, 55, 77, 105, 165, 231, 385, 1155]
SC = ("/tmp/claude-0/-home-user-darwin-50/"
      "aceaa9ec-f432-5848-a506-39c59179b415/scratchpad/")
P = pickle.load(open(SC + PKL, "rb"))
b = IntBasis.covering(P)
r = b.rows(P)
assert b.overflow_headroom(r) < 1.0
n = len(P)
print(f"{PKL}: {n} points, centre ({CX},{CY})  [{time.time()-t0:.0f}s]",
      flush=True)
C = Point(K.rational(CX), K.rational(CY))
ci = next((i for i, p in enumerate(P) if p == C), None)
print(f"centre is a vertex: {ci is not None}", flush=True)
rc = b.rows([C])[0]
dm, D2 = b.dim, b.D * b.D
d = r - rc
sq = b._field_square(d[:, :dm]) + b._field_square(d[:, dm:])
ok = np.ones(len(sq), dtype=bool)
for j in range(1, dm):
    ok &= sq[:, j] == 0
grp = defaultdict(list)
for off in np.nonzero(ok)[0]:
    D = Fr(int(sq[off, 0]), D2)
    if D > 0:
        grp[D].append(int(off))


def bite_of(D):
    """cos and sin of the turn that sends a radius-sqrt(D) ring onto itself
    at unit distance, when the field holds it."""
    if D < Fr(1, 4):
        return None
    ct = 1 - Fr(1, 2) / D
    s2 = 1 - ct * ct
    num, den = s2.numerator, s2.denominator
    for rr in CLASSES:
        q = Fr(num, den) / rr
        a, c = q.numerator, q.denominator
        ra, rc2 = int(round(a ** .5)), int(round(c ** .5))
        if ra * ra == a and rc2 * rc2 == c:
            return ct, Fr(ra, rc2), rr
    return None


rows = []
for D, mem in grp.items():
    if len(mem) < MINPTS:
        continue
    bi = bite_of(D)
    rows.append((len(mem), D, bi))
rows.sort(reverse=True, key=lambda t: (t[2] is not None, t[0]))
print(f"\n{len(rows)} rings with at least {MINPTS} points, "
      f"{sum(1 for x in rows if x[2])} of them biteable in K:", flush=True)
for cnt, D, bi in rows[:40]:
    if bi:
        ct, co, rr = bi
        rad = f"{float(D)**.5:.3f}"
        tag = (f"BITEABLE: cos={ct}, sin={co}"
               + (f"*sqrt{rr}" if rr != 1 else ""))
    else:
        rad, tag = f"{float(D)**.5:.3f}", "not in K"
    print(f"   D={D} (rho={rad}), {cnt} pts: {tag}", flush=True)
print("DONE", flush=True)
