"""Which distances can K close, and which does X realise?

A rotation rho about p carries a point q at distance d to one at distance
2d sin(theta/2) from it.  Asking for that to be exactly 1 fixes the rotation:
rho + rhobar = 2 - 1/D with D = d^2, so rho is a root of

    t^2 - (2 - 1/D) t + 1,       discriminant  -(4D - 1)/D^2,

and the rotation lives in K exactly when sqrt(1 - 4D) does.  For rational D
that is a question about squarefree parts, because a rational square root
lands in K only through one of its three quadratic subfields:

    K = Q(m, sqrt-3, sqrt-11)  ->  Q(sqrt-3), Q(sqrt-11), Q(sqrt33),

the cubic Q(m) admitting none (2 does not divide 3).  So D is closable over K
iff the squarefree part of 1 - 4D is 1, -3, -11 or 33.

This is what decides whether de Grey's architecture ports.  His chain of
closing distances is 1 -> sqrt3 -> 2 -> 4, with 1 - 4D equal to -3, -11, -15,
-63 = -9.7 -- the four radicands of his field, in order.  K has the first two
and neither of the last two, so the chain breaks after sqrt3 and K needs its
own continuation.
"""
import sys, time
from fractions import Fraction as Fr
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
src = open("/tmp/claude-0/-home-user-darwin-50/"
           "aceaa9ec-f432-5848-a506-39c59179b415/scratchpad/"
           "degreyfield.py").read()
exec(src[:src.index("t0 = time.time()\nrh = ")])

t0 = time.time()


def squarefree(n):
    if n == 0:
        return 0
    s, d = (-1 if n < 0 else 1), abs(n)
    f = 2
    while f * f <= d:
        e = 0
        while d % f == 0:
            d //= f
            e += 1
        if e % 2:
            s *= f
        f += 1
    return s * d


def closable(D):
    """Is there a rotation in K taking a pair at distance sqrt(D) to 1?"""
    r = Fr(1) - 4 * Fr(D)
    if r == 0:
        return False
    s = squarefree(r.numerator * r.denominator)
    return s in (1, -3, -11, 33)


print("de Grey's chain:", [(D, Fr(1) - 4 * Fr(D),
                            squarefree((Fr(1) - 4 * Fr(D)).numerator
                                       * (Fr(1) - 4 * Fr(D)).denominator))
                           for D in (1, 3, 4, 16)], flush=True)
print("closable over K:", [D for D in (1, 3, 4, 16) if closable(D)],
      flush=True)

loe = [D for D in range(1, 130)
       if any(a * a + a * b + b * b == D
              for a in range(-13, 14) for b in range(-13, 14))]
print(f"\nLoeschian D up to 129 that K can close: "
      f"{[D for D in loe if closable(D)]}", flush=True)
small = [Fr(p, q) for q in range(1, 13) for p in range(1, 12 * q + 1)
         if Fr(p, q).denominator == q]
print(f"small rational D that K can close: "
      f"{sorted({str(D) for D in small if closable(D)}, key=lambda s: float(Fr(s)))[:28]}",
      flush=True)

# what X actually realises
rh = [kz(), K1, Z6, kadd(K1, Z6)]
spindle = list(rh) + [kmul(RHO, q) for q in rh[1:]]
seed, seen = [], set()
for p in spindle:
    if p not in seen:
        seen.add(p)
        seed.append(p)
X = close(seed, [r for _, r in rots])
print(f"\nX: {len(X)} points  [{time.time()-t0:.0f}s]", flush=True)

zs = [zof(p) for p in X]
from collections import Counter
occ = Counter()
for i in range(len(X)):
    for j in range(i + 1, len(X)):
        v = abs(zs[i] - zs[j]) ** 2
        r = Fr(round(v * 1584), 1584)
        if abs(float(r) - v) < 1e-7:
            occ[r] += 1
print(f"  {len(occ)} candidate rational squared distances  "
      f"[{time.time()-t0:.0f}s]", flush=True)
good = sorted((D for D in occ if closable(D)), key=float)
print(f"  closable ones X realises: "
      f"{[(str(D), occ[D]) for D in good]}  [{time.time()-t0:.0f}s]",
      flush=True)
