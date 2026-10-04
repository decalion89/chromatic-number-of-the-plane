"""Independent check of the two-prime probe at window (1,1), N = 65 (shares no code with probe2d.py / gdyn.py).

* The rotations: all gamma = (x + iy)/65 with x^2 + y^2 = 65^2 (brute force); these are exactly the rational
  rotations with denominator dividing 65, i.e. G(1,1) = {u rho^j sigma^l : |j|, |l| <= 1}.
* Direct enumeration: every unit cell [a + r, a + 1 - r] x [b + r, b + 1 - r] (a, b mod 65) is cut by the strips of
  every functional c -> Re(conj(c) gamma), Im(conj(c) gamma); polygons in H-representation with exact vertex
  enumeration (hgeom_ref.py, the referee's geometry from probe/indep/, copied unchanged).
* A component is 'main' if a type point (65 h or (65/3)(alpha + beta i)) has the same strip indices.
usage: python3 indep2d.py r"""
import sys, time
from fractions import Fraction as Fr
from hgeom_ref import Poly, ev, fl
r = Fr(sys.argv[1]); lo, hi = r, 1 - r
N = 65
gams = [(Fr(x, N), Fr(y, N)) for x in range(-N, N + 1) for y in range(-N, N + 1) if x * x + y * y == N * N]
fs = []
for (A, B) in gams:
    for f in ((A, B), (B, -A)):          # conj(c) gamma = (A c1 + B c2) + i (B c1 - A c2)
        if f not in fs and (-f[0], -f[1]) not in fs:
            fs.append(f)
print(f"r = {r}: {len(gams)} rotations, {len(fs)} distinct functionals (up to sign)")
t0 = time.time()
comps = []
for a in range(N):
    for b in range(N):
        sq = Poly([(Fr(1), Fr(0), a + hi), (Fr(-1), Fr(0), -(a + lo)), (Fr(0), Fr(1), b + hi), (Fr(0), Fr(-1), -(b + lo))])
        pieces = [sq]
        for f in fs:
            nxt = []
            for P in pieces:
                nxt.extend(P.strip_split(f, lo, hi))
            pieces = nxt
            if not pieces:
                break
        comps.extend(pieces)
print(f"{len(comps)} components ({time.time()-t0:.0f}s)")
types = [('C', (Fr(N, 2), Fr(N, 2)))] + [('Q', (Fr(N * al, 3), Fr(N * be, 3))) for al in (1, 2) for be in (1, 2)]
def idx(p):
    return tuple(fl(ev(f, p) - lo) for f in fs)
tidx = [(nm, idx(TP)) for nm, TP in types]
from collections import Counter
lab = Counter()
extra = []
for P in comps:
    V = P.V
    pt = (sum(v[0] for v in V) / len(V), sum(v[1] for v in V) / len(V))
    ip = idx(pt)
    nm = next((nm for nm, ti in tidx if ti == ip), 'X')
    lab[nm] += 1
    if nm == 'X':
        extra.append(pt)
print("labels:", dict(lab))
def kap(p):
    m = Fr(1, 2)
    for f in fs:
        v = ev(f, p); fr = v - fl(v); m = min(m, fr, 1 - fr)
    return m
for pt in extra[:8]:
    print("   extra component, a point:", pt, " min||.|| at that point =", kap(pt))
