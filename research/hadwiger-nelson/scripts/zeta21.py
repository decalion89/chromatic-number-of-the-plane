"""Triangles need zeta_6, blocking needs a split prime: one field must give both.

A triangle in a unit-distance graph is equilateral, so it is a rotation by 60
degrees, so the coordinate field must contain zeta_6.  The roots of unity in
Q(zeta_7) are exactly mu_14, and 6 does not divide 14 -- so no graph over
Q(zeta_7) has a single triangle, whatever its size.  That is why the blocked
ball came back with a small chromatic number: it is triangle-free, clique
number 2, and every high-chromatic unit-distance graph known leans on
triangles throughout.

Q(zeta_21) = Q(zeta_3, zeta_7) is the smallest field that can have both.  It
contains zeta_3 hence zeta_6 = -zeta_3^2, so the whole Eisenstein lattice and
its triangles live there; and it contains zeta_7, which is where the split
primes supplying blocking directions come from.  Its degree is phi(21) = 12,
and 5 has order 6 in (Z/21)*, so 5 splits into two primes of residue degree 6
rather than staying inert -- still plenty of room in the dual.

This checks both halves: that Q(zeta_7) really is triangle-free, and whether
Q(zeta_21) blocks.
"""
import sys, itertools
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
from fractions import Fraction
from math import gcd
from hn.cyclotomic import CycloField
from hn.homcol import has_homomorphism, minimum_blocking_set

print("order of 5 in (Z/21)*:",
      next(k for k in range(1, 13) if pow(5, k, 21) == 1))

for n in (7, 21):
    F = CycloField(n)
    D = F.degree
    one = F.rational(1)

    def inverse(a, F=F, D=D):
        rows = [list(F.mul(a, tuple(Fraction(1 if j == i else 0)
                                    for j in range(D)))) for i in range(D)]
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

    us, box = set(), 3 if n == 7 else 2
    for coeffs in itertools.product(range(-box, box + 1), repeat=4):
        a = tuple(Fraction(c) for c in coeffs) + tuple(
            Fraction(0) for _ in range(D - 4))
        if not any(a):
            continue
        try:
            u = F.mul(a, inverse(F.conj(a)))
        except (StopIteration, ZeroDivisionError):
            continue
        if F.norm2(u) == one:
            us.add(u)

    # A triangle 0, u, v needs |u| = |v| = |u - v| = 1, which forces
    # v = u * zeta_6.  Counting such pairs tests the claim directly.
    tri = sum(1 for u in us for v in us
              if u != v and F.norm2(F.sub(u, v)) == one)
    print(f"\nQ(zeta_{n}): degree {D}, {len(us)} modulus-one elements sampled, "
          f"{tri} unit-apart pairs among them (triangles through 0)")

    dirs = set()
    for u in us:
        d = 1
        for x in u:
            d = d * x.denominator // gcd(d, x.denominator)
        if d % 5 == 0:
            continue
        v = tuple(int(x * d) for x in u)
        g = 0
        for t in v:
            g = gcd(g, abs(t))
        w = tuple(t // g for t in v) if g > 1 else v
        dirs.add(min(w, tuple(-t for t in w)))
    phi, why = has_homomorphism(sorted(dirs), 5)
    print(f"  {len(dirs)} directions; coset 5-colouring: "
          + ("EXISTS -- does not block" if phi else "NONE -- BLOCKS"))
