"""Verify the field that has both: de Grey's spindle and blocking.

Cubic:  x^3 - 10x^2 + 26x - 11, totally real, no rational root.

    V(m) = -39m^4 + 540m^3 - 666m^2 - 324m + 297  totally positive
    V = 33 . s^2  for some s in Q(m)      -> sqrt33 in F = Q(m, sqrt V)
    5 inert in Q(m)                       -> residue degree 6 in F
    K = F(sqrt-3) = Q(m, sqrt-3, sqrt-11) -> the Moser rotation (5+sqrt-11)/6

The norm N(V) = 33^3 . 436^2 is only necessary; the square root is constructed
here at the three real embeddings and then checked exactly.
"""
import math
from fractions import Fraction as Fr

CUB = (-11, 26, -10)                    # constant, x, x^2 (monic cubic)


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


def add(*xs):
    out = (Fr(0),) * 3
    for x in xs:
        out = tuple(a + b for a, b in zip(out, x))
    return out


def sc(c, x):
    return tuple(Fr(c) * a for a in x)


M = (Fr(0), Fr(1), Fr(0))
mm = cmul(M, M)
m3 = cmul(mm, M)
m4 = cmul(m3, M)
V = add(sc(-39, m4), sc(540, m3), sc(-666, mm), sc(-324, M),
        (Fr(297), Fr(0), Fr(0)))
print("V(m) =", V, flush=True)

roots = []
for k in range(3):
    c2, c1, c0 = CUB[2], CUB[1], CUB[0]
    p = c1 - c2 * c2 / 3.0
    q = 2 * c2 ** 3 / 27.0 - c2 * c1 / 3.0 + c0
    arg = 3 * q / (2 * p) * math.sqrt(-3.0 / p)
    roots.append(2 * math.sqrt(-p / 3.0)
                 * math.cos(math.acos(arg) / 3.0 - 2 * math.pi * k / 3.0)
                 - c2 / 3.0)
vals = [float(V[0]) + float(V[1]) * t + float(V[2]) * t * t for t in roots]
print("roots      =", [f"{t:.9f}" for t in roots], flush=True)
print("V at them  =", [f"{v:.6f}" for v in vals], flush=True)
print("all positive?", all(v > 0 for v in vals), flush=True)

# s with s^2 = V/33: solve for its coefficients at the three embeddings.
target = [v / 33.0 for v in vals]
found = None
for signs in [(1, 1, 1), (1, 1, -1), (1, -1, 1), (1, -1, -1)]:
    ys = [sg * math.sqrt(t) for sg, t in zip(signs, target)]
    A = [[1.0, t, t * t] for t in roots]
    # solve A c = ys
    Mx = [row[:] + [y] for row, y in zip(A, ys)]
    for c in range(3):
        piv = max(range(c, 3), key=lambda r: abs(Mx[r][c]))
        Mx[c], Mx[piv] = Mx[piv], Mx[c]
        f = Mx[c][c]
        Mx[c] = [u / f for u in Mx[c]]
        for r in range(3):
            if r != c:
                g = Mx[r][c]
                Mx[r] = [u - g * w for u, w in zip(Mx[r], Mx[c])]
    co = [Mx[i][3] for i in range(3)]
    rat = [Fr(x).limit_denominator(10 ** 7) for x in co]
    if all(abs(float(a) - b) < 1e-6 for a, b in zip(rat, co)):
        s = tuple(rat)
        if cmul(s, s) == sc(Fr(1, 33), V):
            found = s
            break
print("\ns with s^2 = V/33 :", found, flush=True)
if found:
    print("  exact check  s^2 * 33 == V :", cmul(found, found) == sc(Fr(1, 33), V)
          and sc(33, cmul(found, found)) == V, flush=True)

f5 = [(t ** 3 + CUB[2] * t * t + CUB[1] * t + CUB[0]) % 5 for t in range(5)]
print(f"\ncubic mod 5 has roots at {[t for t in range(5) if f5[t] == 0]} "
      f"-> {'IRREDUCIBLE, 5 inert' if all(f5) else 'reducible'}", flush=True)
print("33 mod 5 =", 33 % 5, "; a quadratic residue?",
      pow(33, 2, 5) == 4 and pow(3, (5 - 1) // 2, 5) == 1, flush=True)
print("so 5 is inert in Q(sqrt33) too, and the residue degree in "
      "F = Q(m, sqrt33) is lcm(3, 2) = 6 >= 3 -- the field BLOCKS", flush=True)
