"""Blocking is decided by the field, before any graph is drawn.

Reduce mod 5.  A coset 5-colouring of a graph with edge module M is an
F_5-linear functional on M/5M nonzero on every edge vector; when M has full
rank and index prime to 5 in O, that is a functional on A = O/5.  Edge vectors
are unit steps, u.ubar = 1, and scaling the whole graph by lambda multiplies
every direction by lambda, which is a linear automorphism of A permuting the
hyperplanes -- so the question is scale-free and reduces to the norm-one set

    N = { x in A : x.sigma(x) = 1 },   sigma = conjugation mod 5,

which by Hilbert 90 is exactly the set of reductions of all unit steps of K.
A graph blocks only if its directions meet every hyperplane of A, so:

    IF ONE HYPERPLANE OF A MISSES N, NO UNIT-DISTANCE GRAPH OVER K BLOCKS,
    whatever its size, however it is built.

One finite computation per field, and it decides in advance whether a field is
worth building in at all.
"""
import sys, itertools, time
from fractions import Fraction
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
from hn.cyclotomic import CycloField

P = 5


def table(d, mulbasis, conjbasis):
    """Structure constants mod 5 and the conjugation matrix."""
    C = [[tuple(x % P for x in mulbasis(i, j)) for j in range(d)]
         for i in range(d)]
    S = [tuple(x % P for x in conjbasis(i)) for i in range(d)]
    return C, S


def analyse(name, d, C, S, show=True):
    t0 = time.time()
    one = tuple([1] + [0] * (d - 1))

    def mul(a, b):
        out = [0] * d
        for i in range(d):
            if not a[i]:
                continue
            ai = a[i]
            for j in range(d):
                if not b[j]:
                    continue
                f = ai * b[j]
                row = C[i][j]
                for k in range(d):
                    if row[k]:
                        out[k] += f * row[k]
        return tuple(x % P for x in out)

    def conj(a):
        out = [0] * d
        for i in range(d):
            if a[i]:
                row = S[i]
                for k in range(d):
                    out[k] += a[i] * row[k]
        return tuple(x % P for x in out)

    N = [a for a in itertools.product(range(P), repeat=d)
         if mul(a, conj(a)) == one]
    seen, hyps = set(), []
    for c in itertools.product(range(P), repeat=d):
        if not any(c) or c in seen:
            continue
        for k in range(1, P):
            seen.add(tuple((k * x) % P for x in c))
        hyps.append(c)
    missed = 0
    for c in hyps:
        if not any(sum(x * y for x, y in zip(c, u)) % P == 0 for u in N):
            missed += 1
    ok = missed == 0
    if show:
        print(f"{name}: degree {d}, |N| = {len(N)}, {len(hyps)} hyperplanes -> "
              + ("CAN BLOCK" if ok else
                 f"CANNOT BLOCK ({missed} hyperplanes miss N)")
              + f"  [{time.time()-t0:.0f}s]", flush=True)
    return ok


def cyclo_tables(n):
    K = CycloField(n)
    d = K.degree

    def e(i):
        return tuple(Fraction(1 if j == i else 0) for j in range(d))

    def red(x):
        return tuple(int(t.numerator * pow(t.denominator, -1, P) % P)
                     for t in x)

    C = [[red(K.mul(e(i), e(j))) for j in range(d)] for i in range(d)]
    S = [red(K.conj(e(i))) for i in range(d)]
    return d, C, S


def quad_tables(n, m):
    """K = Q(zeta_n, sqrt(m)) with m < 0, basis (cyclo basis) x (1, sqrt m)."""
    K = CycloField(n)
    h = K.degree
    d = 2 * h

    def e(i):
        return tuple(Fraction(1 if j == i else 0) for j in range(h))

    def red(x):
        return tuple(int(t.numerator * pow(t.denominator, -1, P) % P)
                     for t in x)

    C = [[None] * d for _ in range(d)]
    for i in range(d):
        for j in range(d):
            ii, si = i % h, i // h
            jj, sj = j % h, j // h
            prod = red(K.mul(e(ii), e(jj)))
            if si and sj:
                prod = tuple((m % P) * t % P for t in prod)
                C[i][j] = prod + tuple([0] * h)
            elif si or sj:
                C[i][j] = tuple([0] * h) + prod
            else:
                C[i][j] = prod + tuple([0] * h)
    S = []
    for i in range(d):
        ii, si = i % h, i // h
        c = red(K.conj(e(ii)))
        S.append(tuple([0] * h) + tuple((-t) % P for t in c) if si
                 else c + tuple([0] * h))
    return d, C, S


print("-- fields a spindle or a lattice can live in --", flush=True)
for n in (3, 4):
    d, C, S = cyclo_tables(n)
    analyse(f"Q(zeta_{n})", d, C, S)
for m in (-1, -2, -5, -6, -7, -10, -11, -13, -15):
    d, C, S = quad_tables(3, m)
    analyse(f"Q(zeta_3, sqrt{m})", d, C, S)
print("\n-- fields with a seventh root, where the blocking steps live --",
      flush=True)
for n in (5, 7, 8, 9, 12):
    d, C, S = cyclo_tables(n)
    analyse(f"Q(zeta_{n})", d, C, S)
