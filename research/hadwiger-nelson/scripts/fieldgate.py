"""Blocking is a property of the FIELD, not of the graph.

Reduce everything mod 5.  A coset 5-colouring of a graph whose edge module is
M is an F_5-linear functional phi on M/5M that is nonzero on every edge
vector.  When 5 is unramified and does not divide the index, M/5M is the
F_5-algebra A = O/5, and the edge vectors are unit steps: elements u of K with
u.ubar = 1, where bar is the nontrivial automorphism of K over its maximal
real subfield.  By Hilbert 90 every such u is a/abar, so the reductions of ALL
possible unit steps form

    N  =  { a . sigma(a)^{-1} : a in A* },       sigma = conjugation mod 5,

a subgroup of A* computable once and for all from the field.  A graph over K
blocks exactly when its own directions meet every hyperplane ker(phi); it can
only do so if N does.  Hence

    IF SOME HYPERPLANE OF A MISSES N, NO UNIT-DISTANCE GRAPH OVER K
    BLOCKS -- whatever its size, however it is built.

and conversely if every hyperplane meets N the field is not the obstruction.
This is one finite computation per field, and it decides in advance whether a
field is worth building in.
"""
import sys, itertools
from fractions import Fraction
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
from hn.cyclotomic import CycloField

P = 5


def analyse(name, d, mul, conj, verbose=True):
    """mul and conj act on tuples in F_5^d."""
    els = list(itertools.product(range(P), repeat=d))
    zero = tuple([0] * d)
    inv = {}
    for a in els:
        if a == zero:
            continue
        for b in els:
            if b == zero:
                continue
            if mul(a, b) == tuple([1] + [0] * (d - 1)):
                inv[a] = b
                break
    N = set()
    for a in els:
        if a not in inv:
            continue
        N.add(mul(a, inv[conj(a)]))
    hyps = []
    seen = set()
    for c in els:
        if c == zero or c in seen:
            continue
        for k in range(1, P):
            seen.add(tuple((k * x) % P for x in c))
        hyps.append(c)
    missed = []
    for c in hyps:
        if not any(sum(x * y for x, y in zip(c, u)) % P == 0 for u in N):
            missed.append(c)
    verdict = ("CAN BLOCK: every hyperplane meets N"
               if not missed else
               f"CANNOT BLOCK: {len(missed)} of {len(hyps)} hyperplanes miss N")
    print(f"{name}: degree {d}, |A*| units with inverse {len(inv)}, "
          f"|N| = {len(N)}, {len(hyps)} hyperplanes", flush=True)
    print(f"    {verdict}", flush=True)
    return not missed


def from_cyclo(K):
    d = K.degree

    def mul(a, b):
        x = K.mul(tuple(Fraction(t) for t in a), tuple(Fraction(t) for t in b))
        return tuple(int(t.numerator * pow(t.denominator, -1, P) % P) for t in x)

    def conj(a):
        x = K.conj(tuple(Fraction(t) for t in a))
        return tuple(int(t.numerator * pow(t.denominator, -1, P) % P) for t in x)

    return d, mul, conj


# Q(zeta_3) = Q(sqrt-3): the Eisenstein plane, the triangular lattice.
d, mul, conj = from_cyclo(CycloField(3))
analyse("Q(zeta_3), the triangular lattice", d, mul, conj)

# Q(zeta_3, sqrt-11): the field the Moser spindle lives in and cannot leave.
K3 = CycloField(3)


def ext_mul(a, b):
    ha, hb = a[:2], a[2:]
    hc, hd = b[:2], b[2:]

    def m(x, y):
        z = K3.mul(tuple(Fraction(t) for t in x), tuple(Fraction(t) for t in y))
        return tuple(int(t.numerator * pow(t.denominator, -1, P) % P) for t in z)

    def s(x, y):
        return tuple((p + q) % P for p, q in zip(x, y))

    ac, bd = m(ha, hc), m(hb, hd)
    return s(ac, tuple((-11 * t) % P for t in bd)) + s(m(ha, hd), m(hb, hc))


def ext_conj(a):
    ha, hb = a[:2], a[2:]

    def c(x):
        z = K3.conj(tuple(Fraction(t) for t in x))
        return tuple(int(t.numerator * pow(t.denominator, -1, P) % P) for t in z)

    return c(ha) + tuple((-t) % P for t in c(hb))


analyse("Q(zeta_3, sqrt-11), the Moser spindle field", 4, ext_mul, ext_conj)

# Q(zeta_7): where the denominator-29 blocking steps live.
d, mul, conj = from_cyclo(CycloField(7))
analyse("Q(zeta_7), home of the denominator-29 steps", d, mul, conj)
