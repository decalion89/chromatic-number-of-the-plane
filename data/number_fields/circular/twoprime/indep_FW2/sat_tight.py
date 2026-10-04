"""SAT checks on finite pieces of Cay(Gamma, S), Gamma = Z^r x Z/m.

hom_exists(H, p, q):            H -> K_{p/q}?
acyclic_hom_exists(H, p, q):    is there a homomorphism H -> K_{p/q} whose tight digraph D_c has no directed cycle?
                                (UNSAT  <=>  every homomorphism of H to K_{p/q} has a tight cycle  <=>  chi_c(H) = p/q,
                                 by Guichard's characterisation; Lemma 20 of the paper is one half.)
Models are decoded and re-verified independently (properness and acyclicity by DFS).
"""
from pysat.solvers import Solver
from fractions import Fraction as Fr
import itertools, sys


def piece(r, m, S_half, box):
    """Vertices: points of the box (list of (lo,hi) per free coordinate) times Z/m.  Edges v ~ v+s, s in S_half."""
    rng = [range(lo, hi + 1) for (lo, hi) in box]
    V = [tuple(x) + (y,) for x in itertools.product(*rng) for y in range(m)]
    idx = {v: i for i, v in enumerate(V)}
    E = set()
    for v in V:
        for s in S_half:
            w = tuple(a + b for a, b in zip(v[:r], s[:r])) + (((v[r] + s[r]) % m),)
            if w in idx and w != v:
                a, b = idx[v], idx[w]
                E.add((min(a, b), max(a, b)))
    return V, sorted(E)


class Enc:
    def __init__(self):
        self.n = 0
        self.cls = []

    def new(self):
        self.n += 1
        return self.n


def encode_hom(N, E, p, q, enc):
    x = [[enc.new() for c in range(p)] for v in range(N)]
    for v in range(N):
        enc.cls.append(x[v][:])
        for c1 in range(p):
            for c2 in range(c1 + 1, p):
                enc.cls.append([-x[v][c1], -x[v][c2]])
    for (u, v) in E:
        for c1 in range(p):
            for c2 in range(p):
                d = (c2 - c1) % p
                if not (q <= d <= p - q):
                    enc.cls.append([-x[u][c1], -x[v][c2]])
    return x


def encode_acyclic(N, E, p, q, x, enc, L=None):
    if L is None:
        L = N - 1
    # y[v][k] <-> level(v) >= k, k = 1..L
    y = [[enc.new() for k in range(L)] for v in range(N)]
    for v in range(N):
        for k in range(L - 1):
            enc.cls.append([-y[v][k + 1], y[v][k]])
    arcs = []
    for (u, v) in E:
        for (a, b) in ((u, v), (v, u)):
            t = enc.new()
            arcs.append((a, b, t))
            for c in range(p):
                enc.cls.append([-x[a][c], -x[b][(c + q) % p], t])
            # level(b) >= level(a) + 1
            enc.cls.append([-t, y[b][0]])
            for k in range(L - 1):
                enc.cls.append([-t, -y[a][k], y[b][k + 1]])
            enc.cls.append([-t, -y[a][L - 1]])
    return y, arcs


def decode(model_set, x, N, p):
    col = []
    for v in range(N):
        cs = [c for c in range(p) if x[v][c] in model_set]
        col.append(cs[0])
    return col


def verify_hom(col, E, p, q):
    return all(q <= (col[v] - col[u]) % p <= p - q for (u, v) in E)


def tight_digraph_has_cycle(col, E, p, q, N):
    adj = [[] for _ in range(N)]
    for (u, v) in E:
        if (col[v] - col[u]) % p == q:
            adj[u].append(v)
        if (col[u] - col[v]) % p == q:
            adj[v].append(u)
    state = [0] * N
    sys.setrecursionlimit(100000)
    for s in range(N):
        if state[s]:
            continue
        stack = [(s, iter(adj[s]))]
        state[s] = 1
        while stack:
            v, it = stack[-1]
            nxt = next(it, None)
            if nxt is None:
                state[v] = 2
                stack.pop()
            elif state[nxt] == 1:
                return True
            elif state[nxt] == 0:
                state[nxt] = 1
                stack.append((nxt, iter(adj[nxt])))
    return False


def hom_exists(N, E, p, q, solver='cadical153'):
    enc = Enc()
    x = encode_hom(N, E, p, q, enc)
    with Solver(name=solver, bootstrap_with=enc.cls) as s:
        ok = s.solve()
        if ok:
            col = decode(set(l for l in s.get_model() if l > 0), x, N, p)
            assert verify_hom(col, E, p, q)
            return col
        return None


def acyclic_hom_exists(N, E, p, q, solver='cadical153'):
    enc = Enc()
    x = encode_hom(N, E, p, q, enc)
    encode_acyclic(N, E, p, q, x, enc)
    with Solver(name=solver, bootstrap_with=enc.cls) as s:
        ok = s.solve()
        if ok:
            col = decode(set(l for l in s.get_model() if l > 0), x, N, p)
            assert verify_hom(col, E, p, q)
            assert not tight_digraph_has_cycle(col, E, p, q, N)
            return col
        return None


def farey_pred(p, q, nmax):
    """largest a/b < p/q with a <= nmax, b >= 1, a >= 2b"""
    best = None
    for a in range(2, nmax + 1):
        b = (a * q) // p + 1  # smallest b with a/b < p/q
        if a < 2 * b:
            continue
        f = Fr(a, b)
        if best is None or f > best[0]:
            best = (f, a, b)
    return best
