"""pentagons.py d D [limit]: unit pentagons over Q(sqrt d): u1 + u2 - v1 = Q and Q = v2 + v3 with all five unit
vectors in Q(sqrt d)^2; u1, u2, v1 range over U_D, and v2, v3 are the two intersections of the unit circles about 0
and Q, which lie in K^2 exactly when |Q|^2 (4 - |Q|^2) is a square in K (with |Q|^2 != 0, 4). Prints the new
directions v2, v3 (exact, as (x, y) in K) and their denominators."""
import sys, json
from fractions import Fraction as Fr
from math import isqrt, lcm
import numpy as np
sys.path.insert(0, __import__("os").path.dirname(__import__("os").path.abspath(__file__)))
from units_fast import units_fast
d, D = int(sys.argv[1]), int(sys.argv[2]); LIMIT = int(sys.argv[3]) if len(sys.argv) > 3 else 20
U = [tuple(u) for u in units_fast(d, D)]
def is_sq_rat(q):
    if q < 0: return None
    a, b = q.numerator, q.denominator; ra, rb = isqrt(a), isqrt(b)
    return Fr(ra, rb) if ra * ra == a and rb * rb == b else None
def sqrt_K(x, y):
    if y == 0:
        r = is_sq_rat(x)
        if r is not None: return (r, Fr(0))
        r = is_sq_rat(x / d)
        return (Fr(0), r) if r is not None else None
    nrm = is_sq_rat(x * x - d * y * y)
    if nrm is None: return None
    for s in (nrm, -nrm):
        p = is_sq_rat((x + s) / 2)
        if p is not None and p != 0:
            q = y / (2 * p)
            if p * p + d * q * q == x and 2 * p * q == y: return (p, q)
    return None
def kmul(a, b): return (a[0] * b[0] + d * a[1] * b[1], a[0] * b[1] + a[1] * b[0])
def kinv(a):
    n = a[0] * a[0] - d * a[1] * a[1]; return (a[0] / n, -a[1] / n)
seen, found = set(), []
Us = sorted(set(U))
# Q = u1 + u2 - v1 as an integer 4-vector over D; dedupe Q
Qs = {}
for i, u1 in enumerate(Us):
    for u2 in Us[i:]:
        s2 = (u1[0] + u2[0], u1[1] + u2[1], u1[2] + u2[2], u1[3] + u2[3])
        for v1 in Us:
            Q = (s2[0] - v1[0], s2[1] - v1[1], s2[2] - v1[2], s2[3] - v1[3])
            if Q not in Qs: Qs[Q] = (u1, u2, v1)
print(f"d={d} D={D}: {len(U)} directions, {len(Qs)} distinct 3-step points Q", flush=True)
for Q, (u1, u2, v1) in Qs.items():
    a, b, c, e = Q
    X = Fr(a * a + d * b * b + c * c + d * e * e, D * D); Y = Fr(2 * (a * b + c * e), D * D)
    if (X, Y) in seen: continue
    if (X, Y) == (0, 0) or (X, Y) == (4, 0) or (X, Y) == (1, 0): continue
    # |Q|^2 (4 - |Q|^2) square in K?
    w = kmul((X, Y), (4 - X, -Y))
    r = sqrt_K(*w)
    seen.add((X, Y))
    if r is None: continue
    # v = Q/2 +- i Q * sqrt(4 - |Q|^2)/(2|Q|) ; with r = |Q| sqrt(4 - |Q|^2): factor = r / (2 |Q|^2)
    fac = kmul(r, kinv((2 * X, 2 * Y)))
    qx, qy = (Fr(a, D), Fr(b, D)), (Fr(c, D), Fr(e, D))
    sols = []
    for sg in (1, -1):
        # i*Q = (-qy, qx); v = Q/2 + sg * fac * (i Q)
        vx = (qx[0] / 2 - sg * kmul(fac, qy)[0], qx[1] / 2 - sg * kmul(fac, qy)[1])
        vy = (qy[0] / 2 + sg * kmul(fac, qx)[0], qy[1] / 2 + sg * kmul(fac, qx)[1])
        n2 = (kmul(vx, vx)[0] + kmul(vy, vy)[0], kmul(vx, vx)[1] + kmul(vy, vy)[1])
        assert n2 == (1, 0), n2
        den = lcm(*(z.denominator for z in (vx[0], vx[1], vy[0], vy[1])))
        sols.append((den, (vx, vy)))
    found.append((min(s[0] for s in sols), (X, Y), Q, (u1, u2, v1), sols))
found.sort(key=lambda f: f[0])
print(f"{len(found)} pentagon classes (|Q|^2 values) found", flush=True)
for den, XY, Q, uuv, sols in found[:LIMIT]:
    print(f"  |Q|^2 = {XY[0]} + {XY[1]} sqrt{d}; new directions with denominators {[s[0] for s in sols]}; Q = {Q}/{D}")
json.dump([{"Q": list(Q), "D": D, "uuv": [list(x) for x in uuv], "dens": [s[0] for s in sols],
            "v": [[[str(c) for c in s[1][0]], [str(c) for c in s[1][1]]] for s in sols]} for den, XY, Q, uuv, sols in found[:200]],
          open(f"pentagons_{d}_{D}.json", "w"))
