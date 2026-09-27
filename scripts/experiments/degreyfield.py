"""de Grey's recipe, over a field that can block.

K = Q(m, sqrt-3, sqrt-11) with m a root of x^3 - 10x^2 + 26x - 11, degree 12.
It carries zeta_6 for the triangles and the Moser rotation (5 + sqrt-11)/6 for
the spindle, and 5 is inert in Q(m) and in Q(sqrt33), so the residue degree in
F = Q(m, sqrt33) is lcm(3,2) = 6 -- far above the bound of 3 that a sixth
colour demands, where every multiquadratic field sits at 2 or less.

The tower is L = Q[T]/(cubic), then a + b.s with s^2 = -11 over L(zeta_3), and
complex conjugation is zeta_3 -> zeta_3^2 together with s -> -s.

What is built on it: the Moser spindle as a seed, closed under the rotations
the field provides, in the manner of Sa -- which is twelve images of a
39-point set.  Measured each round: points, edges, chromatic number, and
whether the direction set blocks.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time, itertools
from fractions import Fraction as Fr
from math import gcd
import cmath
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
from hn.homcol import has_homomorphism
from pysat.solvers import Solver

CUB = (-11, 26, -10)          # x^3 - 10x^2 + 26x - 11, constant first


def lmul(x, y):
    r = [Fr(0)] * 5
    for i, a in enumerate(x):
        if a:
            for j, b in enumerate(y):
                r[i + j] += a * b
    for k in (4, 3):
        c = r[k]
        if c:
            r[k] = Fr(0)
            r[k - 3] -= c * CUB[0]
            r[k - 2] -= c * CUB[1]
            r[k - 1] -= c * CUB[2]
    return tuple(r[:3])


def ladd(x, y):
    return tuple(a + b for a, b in zip(x, y))


def lsub(x, y):
    return tuple(a - b for a, b in zip(x, y))


def lsc(c, x):
    return tuple(Fr(c) * a for a in x)


LZ = (Fr(0),) * 3
L1 = (Fr(1), Fr(0), Fr(0))
LM = (Fr(0), Fr(1), Fr(0))


def linv(x):
    cols = [lmul(x, tuple(Fr(1 if k == j else 0) for k in range(3)))
            for j in range(3)]
    A = [[cols[j][i] for j in range(3)] + [Fr(1 if i == 0 else 0)]
         for i in range(3)]
    for c in range(3):
        p = next(t for t in range(c, 3) if A[t][c])
        A[c], A[p] = A[p], A[c]
        sc = Fr(1) / A[c][c]
        A[c] = [v * sc for v in A[c]]
        for t in range(3):
            if t != c and A[t][c]:
                f = A[t][c]
                A[t] = [u - f * v for u, v in zip(A[t], A[c])]
    return tuple(A[i][3] for i in range(3))


# K element: (a, b, c, d) meaning a + b.w + (c + d.w).s, w = zeta_3, s = sqrt-11
# with w^2 = -1 - w.  Each component is an L element.
def kz():
    return (LZ, LZ, LZ, LZ)


def krat(q):
    return (lsc(q, L1), LZ, LZ, LZ)


def kadd(x, y):
    return tuple(ladd(a, b) for a, b in zip(x, y))


def ksub(x, y):
    return tuple(lsub(a, b) for a, b in zip(x, y))


def _wmul(a0, a1, b0, b1):
    """(a0 + a1 w)(b0 + b1 w) with w^2 = -1 - w."""
    p = lmul(a0, b0)
    q = ladd(lmul(a0, b1), lmul(a1, b0))
    r = lmul(a1, b1)
    return lsub(p, r), lsub(q, r)


def kmul(x, y):
    a0, a1, c0, c1 = x
    b0, b1, d0, d1 = y
    p0, p1 = _wmul(a0, a1, b0, b1)
    q0, q1 = _wmul(c0, c1, d0, d1)
    r0, r1 = _wmul(a0, a1, d0, d1)
    t0, t1 = _wmul(c0, c1, b0, b1)
    return (ladd(p0, lsc(-11, q0)), ladd(p1, lsc(-11, q1)),
            ladd(r0, t0), ladd(r1, t1))


def kconj(x):
    a0, a1, c0, c1 = x
    # w -> w^2 = -1 - w ; s -> -s
    return (lsub(a0, a1), lsc(-1, a1), lsub(lsc(-1, c0), lsc(-1, c1)),
            lsc(1, c1))


K1 = krat(1)
W = (LZ, L1, LZ, LZ)                 # zeta_3
Z6 = ksub(kz(), kmul(W, W))          # -w^2 = zeta_6
S = (LZ, LZ, L1, LZ)                 # sqrt-11
RHO = kmul(krat(Fr(1, 6)), kadd(krat(5), S))


def knorm2(x):
    return kmul(x, kconj(x))


print("zeta_6 unit? ", knorm2(Z6) == K1, flush=True)
print("rho unit?    ", knorm2(RHO) == K1, flush=True)
print("|1-rho|^2=1/3?", knorm2(ksub(K1, RHO)) == krat(Fr(1, 3)), flush=True)


# A rotation of chord c lies in K iff (2-c)^2 - 4 is a square there; with
# K = F(sqrt-3) and the value negative that is 3c(4-c) being a square in
# F = Q(m, sqrt33).  For rational c the value is rational, so it must be a
# rational square or 33 times one -- the hexagonal chord 1 gives 9, and the
# Moser chord 1/3 gives 11/3, which is 33 over a square.
import math


def sqrt_in_F(val):
    """sqrt(val) as a K element, for rational val > 0 that F contains."""
    num, den = val.numerator, val.denominator
    t = num * den
    r = math.isqrt(t)
    if r * r == t:
        return krat(Fr(r, den))
    if t % 33 == 0:
        r = math.isqrt(t // 33)
        if r * r * 33 == t:
            # sqrt33 = sqrt-3 . sqrt-11 / ... : (2w+1) is sqrt-3
            s3 = kadd(kmul(krat(2), W), K1)          # 2 zeta_3 + 1 = sqrt-3
            s33 = kmul(s3, S)                        # sqrt-3 . sqrt-11
            return kmul(krat(Fr(-r, den)), s33) if False else \
                kmul(krat(Fr(r, den)), s33)
    return None


assert knorm2(kadd(kmul(krat(2), W), K1)) == krat(3), "2w+1 must be sqrt-3"
s33 = kmul(kadd(kmul(krat(2), W), K1), S)
assert kmul(s33, s33) == krat(33), "sqrt-3 . sqrt-11 must square to 33"

rots = []
for num in range(1, 40):
    for den in (1, 2, 3, 4, 6, 8, 9, 12):
        c = Fr(num, den)
        if not (0 < c < 4) or any(c == cc for cc, _ in rots):
            continue
        disc = 3 * c * (4 - c)
        r = sqrt_in_F(disc)
        if r is None:
            continue
        # rho = ((2-c) + sqrt(disc)/(-3) . sqrt-3 ) / 2, with sqrt-3 = 2w+1
        # X^2 - (2-c) X + 1 = 0 -> X = ((2-c) +- sqrt((2-c)^2-4))/2 and
        # (2-c)^2 - 4 = -disc/3, so the root is r . sqrt-3 / 3.
        s3 = kadd(kmul(krat(2), W), K1)
        rho = kmul(krat(Fr(1, 2)),
                   kadd(krat(2 - c), kmul(krat(Fr(1, 3)), kmul(r, s3))))
        if knorm2(rho) == K1 and knorm2(ksub(K1, rho)) == krat(c):
            rots.append((c, rho))
print(f"\n{len(rots)} rotations in K: "
      + ", ".join(str(c) for c, _ in rots), flush=True)


# Sa is the closure of a 39-point seed under a 12-element group of rotations
# and reflections about the origin.  The same move here, with the six
# rotations the field provides plus complex conjugation, starting from the
# Moser spindle.
def zof(x):
    a0, a1, c0, c1 = x
    roots = None
    import cmath as _c
    w = _c.exp(2j * _c.pi / 3)
    s = _c.sqrt(-11)
    # m numerically: the largest real root of the cubic
    import math as _m
    c2, c1r, c0r = CUB[2], CUB[1], CUB[0]
    p = c1r - c2 * c2 / 3.0
    q = 2 * c2 ** 3 / 27.0 - c2 * c1r / 3.0 + c0r
    arg = 3 * q / (2 * p) * _m.sqrt(-3.0 / p)
    mm = 2 * _m.sqrt(-p / 3.0) * _m.cos(_m.acos(arg) / 3.0) - c2 / 3.0

    def lv(t):
        return float(t[0]) + float(t[1]) * mm + float(t[2]) * mm * mm

    return (lv(a0) + lv(a1) * w) + (lv(c0) + lv(c1) * w) * s


def close(seed, gens, cap=4000):
    pts, frontier = list(seed), list(seed)
    seen = set(seed)
    for _ in range(12):
        nxt = []
        for p in frontier:
            for g in gens:
                for q in (kmul(g, p), kconj(kmul(g, p))):
                    if q not in seen:
                        seen.add(q)
                        pts.append(q)
                        nxt.append(q)
        frontier = nxt
        if not nxt or len(pts) > cap:
            break
    return pts


def measure(pts, name):
    zs = [zof(p) for p in pts]
    n = len(pts)
    cells = {}
    for i, z in enumerate(zs):
        cells.setdefault((int(z.real // 1), int(z.imag // 1)), []).append(i)
    E = []
    for i, z in enumerate(zs):
        cx, cy = int(z.real // 1), int(z.imag // 1)
        for a in (-1, 0, 1):
            for b in (-1, 0, 1):
                for j in cells.get((cx + a, cy + b), ()):
                    if j <= i:
                        continue
                    if abs(abs(z - zs[j]) - 1) > 1e-7:
                        continue
                    if knorm2(ksub(pts[j], pts[i])) == K1:
                        E.append((i, j))
    chi = None
    for k in (3, 4, 5, 6):
        cls = [[1 + v * k + c for c in range(k)] for v in range(n)]
        for a, b in E:
            for c in range(k):
                cls.append([-(1 + a * k + c), -(1 + b * k + c)])
        with Solver(name="cd19", bootstrap_with=cls) as sv:
            if sv.solve():
                chi = k
                break
    vecs = []
    for a, b in E:
        for d in (ksub(pts[b], pts[a]), ksub(pts[a], pts[b])):
            flat = []
            for part in d:
                flat.extend(part)
            vecs.append(tuple(flat))
    dn = 1
    for v in vecs:
        for q in v:
            dn = dn * q.denominator // gcd(dn, q.denominator)
    iv = set()
    for v in vecs:
        w2 = tuple(int(q * dn) for q in v)
        g2 = 0
        for t in w2:
            g2 = gcd(g2, abs(t))
        iv.add(tuple(t // g2 for t in w2) if g2 > 1 else w2)
    iv = sorted(iv)
    blocks = [n2 for n2 in (2, 3, 4, 5)
              if has_homomorphism(iv, n2)[0] is None]
    print(f"  {name}: {n} points, {len(E)} edges, chi = {chi}, "
          f"{len(iv)} directions, blocks at {blocks}", flush=True)
    return pts, E, chi


t0 = time.time()
rh = [kz(), K1, Z6, kadd(K1, Z6)]
spindle = list(rh) + [kmul(RHO, q) for q in rh[1:]]
seed, seen = [], set()
for p in spindle:
    if p not in seen:
        seen.add(p)
        seed.append(p)
measure(seed, "Moser spindle")
for m in (1, 2, 3, 4, 5, 6):
    gens = [r for _, r in rots[:m]]
    pts = close(seed, gens)
    if len(pts) > 1200:
        print(f"  closure under {m} rotations: {len(pts)} points, too big",
              flush=True)
        break
    measure(pts, f"closure under {m} rotations")
    print(f"    [{time.time()-t0:.0f}s]", flush=True)
