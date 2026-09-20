"""The first 5-chromatic unit-distance graph that is also blocked.

A coset colouring is a homomorphism phi from the edge module M to Z/5 that is
nonzero on every edge vector; colouring v by phi(v - v0) is then a proper
5-colouring, so a graph with one has chi <= 5.  Contrapositive:

    EVERY 6-CHROMATIC UNIT-DISTANCE GRAPH IS BLOCKED.

Blocking at n = 5 is the first gate on the way to chi >= 6, and the only known
5-chromatic graph fails it: de Grey's G lives over Q(sqrt3, sqrt5, sqrt7,
sqrt11) and a multiquadratic field never blocks.

Take U = G u {w.G : w a denominator-29 unit step of Q(zeta_7)}, each w.G the
rotation of G about the origin by w.  U contains G, so chi(U) >= 5.  For the
blocking, two facts:

  * blocking is monotone in the direction set, so a blocking SUBSET suffices;
  * phi -> phi . r is a bijection between the homomorphisms out of M and those
    out of rM for any unit r, so blocking is invariant under a global rotation.

G contains a rotated copy r.Y of Y, and Y keeps the six hexagonal unit edges
at the origin (de Grey deletes only (1/3, 0) and (-1/3, 0)).  So U's
directions contain r.{w . zeta_6^k}, and that blocks exactly when
{w . zeta_6^k} does -- which is the 300-direction set measured here.
"""
import sys, itertools
from fractions import Fraction
from math import gcd
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
from hn.cyclotomic import CycloField
from hn.homcol import has_homomorphism, minimum_blocking_set
from hn.degrey import build_Y
from hn.graph import build_graph
from hn.geometry import DEGREY_FIELD

Y = build_Y()
g = build_graph(Y)
E = list(g.edges())
origin = next(i for i, p in enumerate(Y)
              if p.x == DEGREY_FIELD.zero() and p.y == DEGREY_FIELD.zero())
nb = [j for i, j in E if i == origin] + [i for i, j in E if j == origin]
hexd = set()
for j in nb:
    hexd.add((str(Y[j].x - Y[origin].x), str(Y[j].y - Y[origin].y)))
print(f"Y: {len(Y)} points, {len(E)} edges; the origin has {len(nb)} "
      f"neighbours spanning {len(hexd)} directions", flush=True)
for d in sorted(hexd):
    print(f"    ({d[0]}, {d[1]})", flush=True)

K = CycloField(21)
D = K.degree


def cyc_inv(a):
    rows = [list(K.mul(a, tuple(Fraction(1 if j == i else 0) for j in range(D))))
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


def denom(t):
    d = 1
    for x in t:
        d = d * x.denominator // gcd(d, x.denominator)
    return d


Z6 = K.neg(K.mul(K.zeta(7), K.zeta(7)))
raw, gens = set(), []
for coeffs in itertools.product(range(-3, 4), repeat=4):
    if not any(coeffs):
        continue
    a, p = K.zero(), K.rational(1)
    for c in coeffs:
        a = K.add(a, tuple(Fraction(c) * x for x in p))
        p = K.mul(p, K.zeta(3))
    try:
        u = K.mul(a, cyc_inv(K.conj(a)))
    except (StopIteration, ZeroDivisionError):
        continue
    if K.norm2(u) == K.rational(1) and denom(u) == 29 and u not in raw:
        gens.append(u)
        v = u
        for _ in range(6):
            raw.add(v)
            v = K.mul(v, Z6)
steps = sorted(raw)


def as_int_vecs(elts):
    den = 1
    for v in elts:
        den = den * denom(v) // gcd(den, denom(v))
    out = []
    for v in elts:
        w = tuple(int(x * den) for x in v)
        g2 = 0
        for t in w:
            g2 = gcd(g2, abs(t))
        w = tuple(t // g2 for t in w) if g2 > 1 else w
        if w not in out:
            out.append(w)
    return out


iv = as_int_vecs(steps)
print(f"\n{len(gens)} generating rotations, {len(steps)} directions "
      f"({len(iv)} projectively distinct)", flush=True)
phi, _ = has_homomorphism(iv, 5)
print("  full set: " + ("BLOCKS" if phi is None else "admits a coset colouring"),
      flush=True)
best = minimum_blocking_set(iv, n=5)
print(f"  minimum blocking subset: {len(best) if best else None} directions",
      flush=True)
