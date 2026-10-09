"""Level r of the plane over a ramified quadratic extension E = Q_p(sqrt(p m)) of Q_p, p = 3 (mod 4), m = 1 or a
non-residue: O = Z_p[pi] with pi^2 = p m, the points are (O/pi^r)^2 and the unit vectors the solutions of
x^2 + y^2 = 1 in O/pi^r (the form is anisotropic over E, so every unit vector of E^2 is integral and reduces to one).
Prints whether the Cayley graph has a proper K-colouring (CaDiCaL). notes/local_global.md, section 3.

usage: python3 ramified_levels.py p m r K"""
import sys
from pysat.solvers import Solver


def level_graph(p, m, r):
    A, B = p ** ((r + 1) // 2), p ** (r // 2)          # a + b pi with a mod A, b mod B
    pm = p * m

    def red(a, b):
        return (a % A, b % B if B > 1 else 0)

    def mul(e, f):
        return red(e[0] * f[0] + pm * e[1] * f[1], e[0] * f[1] + e[1] * f[0])

    def add(e, f):
        return red(e[0] + f[0], e[1] + f[1])

    R = [(a, b) for a in range(A) for b in range(max(B, 1))]
    one = red(1, 0)
    U = [(x, y) for x in R for y in R if add(mul(x, x), mul(y, y)) == one]
    V = [(x, y) for x in R for y in R]
    idx = {v: i for i, v in enumerate(V)}
    E = []
    for i, (x, y) in enumerate(V):
        for (u, w) in U:
            j = idx[(add(x, u), add(y, w))]
            if j > i:
                E.append((i, j))
    return V, U, E


def colourable(n, E, K):
    s = Solver(name="cadical153")
    for i in range(n):
        s.add_clause([i * K + c + 1 for c in range(K)])
    for i, j in E:
        for c in range(K):
            s.add_clause([-(i * K + c + 1), -(j * K + c + 1)])
    s.add_clause([1])
    r = s.solve()
    s.delete()
    return r


if __name__ == "__main__":
    p, m, r, K = (int(a) for a in sys.argv[1:5])
    V, U, E = level_graph(p, m, r)
    print(f"Q_{p}(sqrt({p * m})) level {r}: {len(V)} points, {len(U)} unit vectors, {len(E)} edges; "
          f"{K}-colourable: {colourable(len(V), E, K)}")
