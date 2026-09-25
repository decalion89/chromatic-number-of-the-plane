"""Lower bounds for chi(F_q^2) from the spectrum (notes/local_colourings.md, section 14).

The unit-distance graph of F_q^2 (x^2 + y^2) is the Cayley graph of F_q^2 on the q -/+ 1 unit vectors; its
eigenvalues are lambda_xi = sum_u cos(2 pi xi.u / q).  With --inert the form is x^2 - n y^2, n the least
non-square: the anisotropic plane G_q of section 4 (q + 1 unit vectors), the same graph when q = 3 (mod 4).  Hoffman's ratio bound gives
    alpha <= n (-lambda_min) / (d - lambda_min),   chi >= chi_f = n / alpha      (vertex-transitive).
Every eigenvalue is computed here in interval arithmetic (mpmath.iv), so the bound on alpha is rigorous.
The script also prints the interval colouring of section 12 (m consecutive parallel lines per colour).

usage: python3 scripts/finite_hoffman.py [--inert] q [q ...]
"""
import sys
from mpmath import iv

iv.prec = 100


def units(q, kind='std'):
    sq = {t * t % q for t in range(q)}
    n0 = -1 if kind == 'std' else next(a for a in range(2, q) if a not in sq)
    return [(a, b) for a in range(q) for b in range(q) if (a * a - n0 * b * b) % q == 1]


def hoffman(q, kind='std'):
    """(degree, interval containing lambda_min, largest integer alpha allowed by Hoffman's bound)"""
    U = units(q, kind); d = len(U); n = q * q
    C = [iv.cos(2 * iv.pi * t / q) for t in range(q)]
    lmin = None
    for a in range(q):
        for b in range(q):
            if (a, b) == (0, 0):
                continue
            lam = iv.mpf(0)
            for x, y in U:
                lam += C[(a * x + b * y) % q]
            if lmin is None or lam.a < lmin.a:
                lmin = lam
    L = iv.mpf(lmin.a)                        # n x / (d + x) increases with x = -lambda
    bound = (n * (-L) / (d - L)).b
    amax = int(bound)
    return d, lmin, amax


def interval_colouring(q):
    """fewest colours ceil(q/m) of an interval colouring: a^2 + b^2 - r^2 a non-square for r < m"""
    sq = {t * t % q for t in range(q)}
    best = None
    for a in range(q):
        for b in range(q):
            N = (a * a + b * b) % q
            if N == 0:
                continue
            m = 0
            while m < q and (N - m * m) % q not in sq:
                m += 1
            if m and (best is None or -(-q // m) < best[0]):
                best = (-(-q // m), a, b, m)
    return best


if __name__ == '__main__':
    kind = 'inert' if '--inert' in sys.argv else 'std'
    for q in map(int, [a for a in sys.argv[1:] if a != '--inert']):
        d, lmin, amax = hoffman(q, kind)
        n = q * q
        line = (f"q = {q} ({kind}): degree {d}, lambda_min in [{float(lmin.a):.6f}, {float(lmin.b):.6f}], "
                f"alpha <= {amax}, chi >= {-(-n // amax)}")
        if kind == 'std':
            ub = interval_colouring(q)
            line += f"; interval colouring with {ub[0]} colours (a, b, m) = {ub[1:]}"
        print(line, flush=True)
