"""Palette depth of the WIDE rings -- the ones a bite can actually use.

Every palette measurement in this work was taken on a ring of distance
squared at most 4, and every bite that matters lives outside that range: de
Grey's own second-level turn is the bite on a ring of radius four.  So the
two halves of his argument have never been measured on the same object.  A
cap gives a disjunction; the bite on THAT ring sharpens it into a named pair.
A cap on a narrow ring whose bite does not exist in the field is worth
nothing, and a wide ring that the field can bite has never been asked whether
it is capped.

This asks, for each of the twelve rings about Gp's pivot that K can bite, and
on the whole graph rather than a ball: how many colours can the ring show,
and how many can the ring and its centre show between them.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time, pickle
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
import numpy as np
from fractions import Fraction as Fr
from collections import defaultdict
from hn.geometry import DEGREY_FIELD as K, Point
from hn.fast import IntBasis, fast_edges_complete
from hn.homcol import ring_palette_bound

k = int(sys.argv[1]) if len(sys.argv) > 1 else 5
PKL = sys.argv[2] if len(sys.argv) > 2 else "pivG.pkl"
CX = Fr(sys.argv[3]) if len(sys.argv) > 3 else Fr(-2)
CY = Fr(sys.argv[4]) if len(sys.argv) > 4 else Fr(0)
AUTO = (len(sys.argv) > 5 and sys.argv[5] == "auto")
t0 = time.time()
CLASSES = [1, 3, 5, 7, 11, 15, 21, 33, 35, 55, 77, 105, 165, 231, 385, 1155]
SC = ("/tmp/claude-0/-home-user-darwin-50/"
      "aceaa9ec-f432-5848-a506-39c59179b415/scratchpad/")
P = pickle.load(open(SC + PKL, "rb"))
b = IntBasis.covering(P)
r = b.rows(P)
assert b.overflow_headroom(r) < 1.0
E = sorted(set((min(a, c), max(a, c)) for a, c in fast_edges_complete(b, r)))
n = len(P)
print(f"{PKL}: {n} points, {len(E)} edges, k={k}  [{time.time()-t0:.0f}s]",
      flush=True)
cls = [[1 + v * k + c for c in range(k)] for v in range(n)]
for a, c in E:
    for col in range(k):
        cls.append([-(1 + a * k + col), -(1 + c * k + col)])
deg = defaultdict(int)
for _a, _c in E:
    deg[_a] += 1
    deg[_c] += 1
if str(CX) == "0" and str(CY) == "0" and AUTO:
    ci = max(range(n), key=lambda v: deg[v])
    C = P[ci]
    print(f"centre chosen by degree: vertex {ci}, degree {deg[ci]}",
          flush=True)
else:
    C = Point(K.rational(CX), K.rational(CY))
    ci = next((i for i, p in enumerate(P) if p == C), None)
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


def biteable(D):
    if D < Fr(1, 4):
        return None
    ct = 1 - Fr(1, 2) / D
    s2 = 1 - ct * ct
    for rr in CLASSES:
        q = s2 / rr
        a, c = q.numerator, q.denominator
        ra, rc2 = int(round(a ** .5)), int(round(c ** .5))
        if ra * ra == a and rc2 * rc2 == c:
            return f"cos={ct}, sin={Fr(ra,rc2)}" + (f"*sqrt{rr}" if rr != 1
                                                    else "")
    return None


rows = [(D, mem, biteable(D)) for D, mem in grp.items()
        if len(mem) >= k and biteable(D)]
rows.sort(key=lambda t: -len(t[1]))
print(f"{len(rows)} biteable rings with at least {k} points"
      f"  [{time.time()-t0:.0f}s]", flush=True)
found = []
for D, mem, txt in rows:
    pr = ring_palette_bound(cls, n * k, mem, k)
    pc = (ring_palette_bound(cls, n * k, mem + [ci], k)
          if ci is not None else None)
    note = ""
    if D == 1:
        note = "  (free: the centre touches the ring)"
    elif pr < k:
        note = "   <<< CAPPED AND BITEABLE"
        found.append((pr, pc, D, len(mem), txt))
    print(f"   D={D} (rho={float(D)**.5:.3f}), {len(mem)} pts, {txt}: "
          f"ring {pr}, centre+ring {pc}{note}  [{time.time()-t0:.0f}s]",
          flush=True)
print(f"\ncapped AND biteable rings: {len(found)}", flush=True)
for e in sorted(found):
    print(f"   ring palette {e[0]}, centre+ring {e[1]}, D={e[2]}, "
          f"{e[3]} pts, {e[4]}", flush=True)
print("DONE", flush=True)
