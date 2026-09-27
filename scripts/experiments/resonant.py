"""Which rotations of Q(zeta_21) bring a lattice into contact with its copy?

Rotating the Eisenstein patch by a blocking step gave exactly twice the points
and exactly twice the edges -- the copies never touch, so the union is two
independent lattices and stays 3-chromatic.  The Moser rotation is different
because it RESONATES: |1 - rho|^2 = 1/3 exactly, so points at sqrt 3 from the
pivot land one apart.

In general the union gains an edge when |rho p - q| = 1 for lattice points
p, q, which is

    Re(rho * p conj(q))  =  (|p|^2 + |q|^2 - 1) / 2,

a line in rho meeting the unit circle.  Rather than solve it symbolically, this
enumerates modulus-one elements of Q(zeta_21) by Hilbert 90 and simply counts
the cross edges each one produces.  Any rotation with cross edges is a
candidate Moser analogue in a field that blocks; the chromatic number then
says whether it is one.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time, itertools
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
from fractions import Fraction
from math import gcd
from hn.cyclotomic import CycloField
from pysat.solvers import Solver

F = CycloField(21)
D = F.degree
one = F.rational(1)
z3 = F.zeta(7)
z6 = F.neg(F.mul(z3, z3))


def inverse(a):
    rows = [list(F.mul(a, tuple(Fraction(1 if j == i else 0) for j in range(D))))
            for i in range(D)]
    M = [[rows[j][i] for j in range(D)] + [Fraction(1 if i == 0 else 0)]
         for i in range(D)]
    for c in range(D):
        p = next(r for r in range(c, D) if M[r][c])
        M[c], M[p] = M[p], M[c]
        inv = Fraction(1) / M[c][c]
        M[c] = [v * inv for v in M[c]]
        for r in range(D):
            if r != c and M[r][c]:
                f = M[r][c]
                M[r] = [x - f * y for x, y in zip(M[r], M[c])]
    return tuple(M[i][D] for i in range(D))


def patch(r):
    pts = set()
    for a in range(-r, r + 1):
        for b in range(-r, r + 1):
            v = F.zero()
            for _ in range(abs(a)):
                v = F.add(v, one if a > 0 else F.neg(one))
            for _ in range(abs(b)):
                v = F.add(v, z6 if b > 0 else F.neg(z6))
            pts.add(v)
    return sorted(pts)


base = patch(3)
print(f"Eisenstein patch: {len(base)} points", flush=True)

roots = set()
z, p = one, one
for _ in range(21):
    p = F.mul(p, F.zeta(1))
    roots.add(p)
    roots.add(F.neg(p))

cands, seen = [], set()
for coeffs in itertools.product(range(-2, 3), repeat=5):
    if not any(coeffs):
        continue
    a, p = F.zero(), one
    for c in coeffs:
        a = F.add(a, tuple(Fraction(c) * x for x in p))
        p = F.mul(p, F.zeta(1))
    try:
        u = F.mul(a, inverse(F.conj(a)))
    except (StopIteration, ZeroDivisionError):
        continue
    if F.norm2(u) != one or u in seen or u in roots:
        continue
    seen.add(u)
    cands.append(u)
print(f"{len(cands)} non-root-of-unity rotations to try", flush=True)

t0, hits = time.time(), 0
for u in cands:
    img = [F.mul(u, q) for q in base]
    # Count only GENUINE cross edges. Counting every base-image pair at
    # distance one double-counts the edges of the shared points -- the origin
    # is fixed by any rotation, so its six edges appeared twice and read as
    # twelve "cross edges" where the totals said 2 x 120 and there were none.
    bset, iset = set(base), set(img)
    cross = sum(1 for a in bset - iset for b in iset - bset
                if F.norm2(F.sub(a, b)) == one)
    if not cross:
        continue
    hits += 1
    P = sorted(set(base) | set(img))
    E = [(i, j) for i in range(len(P)) for j in range(i + 1, len(P))
         if F.norm2(F.sub(P[j], P[i])) == one]
    k = None
    for kk in (3, 4, 5):
        cls = [[1 + v * kk + c for c in range(kk)] for v in range(len(P))]
        for a, b in E:
            for c in range(kk):
                cls.append([-(1 + a * kk + c), -(1 + b * kk + c)])
        with Solver(name="cd19", bootstrap_with=cls) as s:
            if s.solve():
                k = kk
                break
    print(f"  rotation with {cross} cross edges: {len(P)} points, {len(E)} "
          f"edges, chi = {k}"
          + ("   *** 4-CHROMATIC OVER A BLOCKING FIELD ***" if k and k >= 4
             else "") + f"  [{time.time()-t0:.0f}s]", flush=True)
    if k and k >= 4:
        break
    if hits >= 8:
        break
print(f"{hits} resonant rotations found among {len(cands)}")
