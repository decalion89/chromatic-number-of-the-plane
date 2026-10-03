"""H_q = Cay(F_q^2, {(t, 1/t) : t in F_q^*}) for q = p^2 (p odd), the residue graph at a place where x^2 + y^2 is
isotropic with residue field F_q. Prints whether H_q has a proper K-colouring (CaDiCaL). For q = 25 there is none with
4 colours, so places with residue field F_25 are no gate at four colours (notes/local_global.md, section 5).

usage: python3 hyperbola_plane.py p K"""
import sys
from pysat.solvers import Solver


def h_q(p):
    n = next(a for a in range(2, p) if pow(a, (p - 1) // 2, p) == p - 1)     # F_q = F_p[s]/(s^2 - n)
    els = [(a, b) for a in range(p) for b in range(p)]

    def mul(x, y):
        return ((x[0] * y[0] + n * x[1] * y[1]) % p, (x[0] * y[1] + x[1] * y[0]) % p)

    def add(x, y):
        return ((x[0] + y[0]) % p, (x[1] + y[1]) % p)

    inv = {x: y for x in els for y in els if mul(x, y) == (1, 0)}
    S = [(t, inv[t]) for t in els if t != (0, 0)]
    V = [(x, y) for x in els for y in els]
    vid = {v: i for i, v in enumerate(V)}
    E = []
    for v in V:
        for a, b in S:
            w = (add(v[0], a), add(v[1], b))
            if vid[w] > vid[v]:
                E.append((vid[v], vid[w]))
    return V, S, E


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
    p, K = int(sys.argv[1]), int(sys.argv[2])
    V, S, E = h_q(p)
    print(f"H_{p * p}: {len(V)} vertices, {len(S)} unit vectors, {len(E)} edges; {K}-colourable: {colourable(len(V), E, K)}")
