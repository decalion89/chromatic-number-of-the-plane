"""Which distances can be spindled in this field?

de Grey reaches five colours by finding a pair forced monochromatic in every
4-colouring at distance sqrt3, then spindling it with the rotation of chord
1/3 -- which is (5 +- sqrt-11)/6, the Moser rotation.  That step needs sqrt-11
in the field.

The tower here is K = Q(m, sqrt V, sqrt-3) and V must be POSITIVE at the
identity embedding for the chain to be real, so V cannot be -11 times a
square and sqrt-11 cannot enter that way.  Distance sqrt3 is therefore not
available, and the question becomes which distances ARE.

A pair at distance d spindles iff K contains a rotation with
|1 - rho|^2 = 1/d^2, i.e. rho + rhobar = t = 2 - 1/d^2 with t^2 - 4 a square
in K.  Since K = F(sqrt-3) and t^2 - 4 < 0 for 0 < 1/d^2 < 4, that is

    (4 - t^2)/3  a square in F,

one test per candidate distance.  F = Q(m)(sqrt V), so reduction at a rational
prime where the cubic splits completely and every V(t_i) is a quadratic
residue makes F mod p a product of six copies of F_p, and squareness is tested
componentwise.
"""
import sys, time
from fractions import Fraction as Fr

CUB = (8, 14, -9)                      # T^3 - 9T^2 + 14T + 8, constant first


def cmul(x, y):
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


def V_of(m):
    mm = cmul(m, m)
    m3 = cmul(mm, m)
    m4 = cmul(m3, m)

    def sc(c, x):
        return tuple(Fr(c) * a for a in x)

    def add(*xs):
        out = (Fr(0),) * 3
        for x in xs:
            out = tuple(a + b for a, b in zip(out, x))
        return out

    return add(sc(-39, m4), sc(540, m3), sc(-666, mm), sc(-324, m),
               (Fr(297), Fr(0), Fr(0)))


M = (Fr(0), Fr(1), Fr(0))
V = V_of(M)
print(f"V(m) = {V}", flush=True)

# Primes where the cubic splits completely and every V(t_i) is a QR mod p:
# then F mod p is F_p^6 and squareness is componentwise.
good = []
for p in range(11, 4000):
    if any(p % q == 0 for q in range(2, int(p ** 0.5) + 1)):
        continue
    roots = [t for t in range(p)
             if (t * t * t + CUB[2] * t * t + CUB[1] * t + CUB[0]) % p == 0]
    if len(set(roots)) != 3:
        continue
    vs = [(int(V[0]) + int(V[1]) * t + int(V[2]) * t * t) % p for t in roots]
    if any(v == 0 for v in vs):
        continue
    if all(pow(v, (p - 1) // 2, p) == 1 for v in vs):
        good.append((p, roots, vs))
    if len(good) >= 8:
        break
print(f"{len(good)} split primes with V a residue: "
      + ", ".join(str(p) for p, _, _ in good), flush=True)


def is_square_in_F(num, den):
    """num/den in Q(m) -- as an element of F it is a square iff each of the
    six components is a quadratic residue.  The element here lies in Q(m),
    so its two F-components at each root agree and only three tests are
    needed, but both signs of sqrt V are checked anyway."""
    for p, roots, vs in good:
        for t in roots:
            n = sum(int(c) * pow(t, i, p) for i, c in enumerate(num)) % p
            d = sum(int(c) * pow(t, i, p) for i, c in enumerate(den)) % p
            if d == 0:
                return None
            x = n * pow(d, -1, p) % p
            if x == 0:
                continue
            if pow(x, (p - 1) // 2, p) != 1:
                return False
    return True


# Candidate distances: d^2 rational, and d^2 in Q(m).  Start with rationals.
print("\nrational squared distances d^2 and whether they spindle:", flush=True)
t0 = time.time()
ok = []
for num in range(1, 40):
    for den in (1, 2, 3, 4, 6, 8, 12):
        d2 = Fr(num, den)
        if d2 <= Fr(1, 4):
            continue                    # need 1/d^2 < 4 for a real rotation
        t = 2 - 1 / d2
        val = (4 - t * t) / 3
        if val <= 0:
            continue
        nn = (Fr(val.numerator * val.denominator), Fr(0), Fr(0))
        dd = (Fr(val.denominator * val.denominator), Fr(0), Fr(0))
        r = is_square_in_F(nn, dd)
        if r:
            ok.append(d2)
print(f"  {len(ok)} of the rational candidates spindle: "
      + ", ".join(str(x) for x in ok[:30]) + f"  [{time.time()-t0:.0f}s]",
      flush=True)
print(f"  is sqrt3 among them (d^2 = 3)? {Fr(3) in ok}", flush=True)
