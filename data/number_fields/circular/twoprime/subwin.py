"""Direct exact enumeration (hgeom_ref geometry) of the probe for an arbitrary finite set W of exponent pairs (j, l):
S_W^r = { c : Re(conj(c) u rho^j sigma^l) in [r, 1-r] + Z, u in mu_4, (j, l) in W }, periodic mod 65 Z[i] when
|j|, |l| <= 1.  Prints the number of components, how many contain a type point, and (with an exact LP) the largest
kappa of the others.  usage: python3 subwin.py r 'j,l;j,l;...'"""
import sys, time
from fractions import Fraction as Fr
from hgeom_ref import Poly, ev, fl
from lpexact import lp_max
r = Fr(sys.argv[1]); lo, hi = r, 1 - r
W = [tuple(int(t) for t in p.split(',')) for p in sys.argv[2].split(';')]
N = 65
def gpow(z, e):
    w = (Fr(1), Fr(0)); zz = z if e >= 0 else (z[0], -z[1])
    for _ in range(abs(e)):
        w = (w[0] * zz[0] - w[1] * zz[1], w[0] * zz[1] + w[1] * zz[0])
    return w
RHO = (Fr(3, 5), Fr(4, 5)); SIG = (Fr(5, 13), Fr(12, 13))
fs = []
for (j, l) in W:
    a, b = gpow(RHO, j); c, d = gpow(SIG, l)
    A, B = a * c - b * d, a * d + b * c
    fs += [(A, B), (B, -A)]
if (0, 0) not in W:
    fs = [(Fr(1), Fr(0)), (Fr(0), Fr(1))] + fs     # the cell structure uses gamma = 1 anyway (a weaker set is not periodic otherwise)
t0 = time.time()
comps = []
for a in range(N):
    for b in range(N):
        pieces = [Poly([(Fr(1), Fr(0), a + hi), (Fr(-1), Fr(0), -(a + lo)), (Fr(0), Fr(1), b + hi), (Fr(0), Fr(-1), -(b + lo))])]
        for f in fs:
            nxt = []
            for P in pieces:
                nxt.extend(P.strip_split(f, lo, hi))
            pieces = nxt
            if not pieces: break
        comps.extend(pieces)
types = [(Fr(N, 2), Fr(N, 2))] + [(Fr(N * x, 3), Fr(N * y, 3)) for x in (1, 2) for y in (1, 2)]
def idx(p): return tuple(fl(ev(f, p) - lo) for f in fs)
tid = [idx(T) for T in types]
best = None; nx = 0
for P in comps:
    V = P.V; pt = (sum(v[0] for v in V) / len(V), sum(v[1] for v in V) / len(V))
    if idx(pt) in tid: continue
    nx += 1
    cons = []
    for f in fs:
        v = ev(f, pt); n = fl(v); base = v - n
        cons += [(f[0], f[1], base), (-f[0], -f[1], 1 - base)]
    t = lp_max(cons)[0]
    best = t if best is None or t > best else best
print(f"W={W} r={r}: {len(comps)} components, {len(comps)-nx} main, {nx} other; max kappa of others = {best} ({float(best) if best else 0:.6f}) [{time.time()-t0:.0f}s]")
