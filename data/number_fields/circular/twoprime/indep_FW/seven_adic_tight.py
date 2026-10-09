"""Sanity check of the type-7 mechanism on the 7-adic colourings themselves (own code).
For e in A'_7 the map z -> Re(conj(e) * (z mod 7)) in Z/7 (z in Z_(7)[i]) is a (7,2)-colouring of the Cayley graph of
Z[1/65][i] on the rational rotations (values {2,3,4,5} on mu_8), with constant lifts l(x, gamma) = Re(conj(e) gamma~).
Its tight steps are the gamma with value 2.  We check, for each e: the tight rotations inside G(2,1) contain a positive
relation (so the colouring has a tight cycle, as the draft's argument predicts), and whether G(1,1) already does."""
from fractions import Fraction as Fr
from itertools import combinations
from math import lcm

def cm(a, b): return (a[0]*b[0]-a[1]*b[1], a[0]*b[1]+a[1]*b[0])
def cpow(z, e):
    r = (Fr(1), Fr(0))
    if e < 0: z = (z[0], -z[1]); e = -e
    for _ in range(e): r = cm(r, z)
    return r
RHO, SIG = (Fr(3, 5), Fr(4, 5)), (Fr(5, 13), Fr(12, 13))
UNITS = [(Fr(1), Fr(0)), (Fr(0), Fr(1)), (Fr(-1), Fr(0)), (Fr(0), Fr(-1))]
def G(K, M):
    return [cm(u, cm(cpow(RHO, j), cpow(SIG, l))) for u in UNITS for j in range(-K, K+1) for l in range(-M, M+1)]
def red7(z):
    a, b = z
    inv = lambda d: pow(d % 7, -1, 7)
    return ((a.numerator * inv(a.denominator)) % 7, (b.numerator * inv(b.denominator)) % 7)
A7 = [(2, 3), (2, 4), (3, 2), (3, 5), (4, 2), (4, 5), (5, 3), (5, 4)]
def val(e, g):
    gt = red7(g)
    return (e[0]*gt[0] + e[1]*gt[1]) % 7      # Re(conj(e) g~) = e0 g0 + e1 g1
def pos_rel(T):
    best = None
    for r in (2, 3):
        for S in combinations(T, r):
            if r == 2:
                a, b = S
                if a[0] == -b[0] and a[1] == -b[1]:
                    return (S, (1, 1))
                continue
            a, b, c = S
            det = a[0]*b[1]-a[1]*b[0]
            if det == 0: continue
            x = (-c[0]*b[1]+c[1]*b[0])/det; y = (-a[0]*c[1]+a[1]*c[0])/det
            if x > 0 and y > 0:
                L = lcm(x.denominator, y.denominator)
                m = (int(x*L), int(y*L), L)
                if best is None or sum(m) < sum(best[1]): best = (S, m)
    return best
for e in A7:
    for K, M in ((1, 1), (2, 1), (2, 0), (3, 0)):
        Gs = G(K, M)
        assert all(val(e, g) in (2, 3, 4, 5) for g in Gs)
        T = [g for g in Gs if val(e, g) == 2]
        r = pos_rel(T)
        tag = f"G({K},{M})"
        if r is None:
            print(f"e={e}: {tag}: {len(T)} tight rotations, NO positive relation (they lie in a half-plane)")
        else:
            S, m = r
            assert sum(mm*s[0] for mm, s in zip(m, S)) == 0 and sum(mm*s[1] for mm, s in zip(m, S)) == 0
            assert sum(m) % 7 == 0
            print(f"e={e}: {tag}: {len(T)} tight rotations, shortest positive triple relation {m} (walk length {sum(m)})")
