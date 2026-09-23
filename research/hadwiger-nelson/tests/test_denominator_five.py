"""The denominator-five facts, recomputed.

Exoo and Ismailescu's graph rebuilt from their 23 points; their rotation's
arithmetic; the module test that shows a denominator of five genuinely blocks
coset colourings while the integral 803-graph does not; the common kernel of
the 803-graph's coset colourings being exactly 5M; and the 61-point lattice
graph with three forbidden distances.
"""
import itertools
import json
from fractions import Fraction as Fr
from math import gcd

import numpy as np
from pysat.solvers import Solver

from hn.field import Field
from hn.geometry import Point
from hn.graph import build_graph
from hn.homcol import has_homomorphism

ROOT = "/home/user/darwin-50/research/hadwiger-nelson"
F = Field((3, 11))
R3, R11 = F.sqrt(3), F.sqrt(11)
S = [(0,0,0,0),(0,0,0,-4),(0,0,-6,-2),(0,0,-6,2),(-6,0,0,-2),(-4,0,0,0),(-4,0,-6,-2),
     (-4,0,-6,2),(-2,0,0,-2),(-2,0,-6,-4),(-2,0,-6,4),(0,-6,-6,0),(-5,-3,3,3),(-5,3,-3,3),
     (-2,-6,0,0),(-2,-6,0,-4),(-2,6,0,0),(-2,6,0,-4),(-6,-6,0,0),(-6,6,0,0),(-4,0,0,-4),
     (0,0,-12,0),(-8,0,0,0)]
EXTRA = [(-2,0,0,-6),(8,0,0,4),(-4,-6,-6,-4),(-4,6,6,-4),(-3,-3,-3,-5),(-4,0,-12,4),
         (-4,0,12,4),(7,-3,3,3),(7,3,-3,3)]


def _pt(a, b, c, d):
    q = lambda v: F.rational(Fr(v, 12))
    return Point(q(a) * R3 + q(b) * R11, q(c) + q(d) * R3 * R11)


def _ei_points():
    T = set()
    for a, b, c, d in S:
        T |= {(a, b, c, d), (a, b, -c, -d), (-a, -b, c, d)}
    half = F.rational(Fr(1, 2))
    V = {}
    for t in T:
        p = _pt(*t)
        for _ in range(6):
            V[p] = True
            p = Point(p.x * half - p.y * R3 * half, p.x * R3 * half + p.y * half)
    return list(V)


def _two_distance_edges(pts):
    one, four = F.rational(1), F.rational(4)
    e1, e2 = [], []
    for i, j in itertools.combinations(range(len(pts)), 2):
        dd = (pts[i].fx - pts[j].fx) ** 2 + (pts[i].fy - pts[j].fy) ** 2
        if abs(dd - 1) < 1e-9 and pts[i].dist2(pts[j]) == one:
            e1.append((i, j))
        elif abs(dd - 4) < 1e-9 and pts[i].dist2(pts[j]) == four:
            e2.append((i, j))
    return e1, e2


def test_exoo_ismailescu_counts_and_forced_pair():
    V = _ei_points()
    e1, e2 = _two_distance_edges(V)
    assert (len(V), len(e1), len(e2)) == (205, 966, 423)
    A, B = _pt(*EXTRA[0]), _pt(*EXTRA[1])
    VH = list(dict.fromkeys(V + [_pt(*e) for e in EXTRA]))
    e1, e2 = _two_distance_edges(VH)
    assert (len(VH), len(e1), len(e2)) == (214, 1004, 446)
    assert A.dist2(B) == F.rational(25)
    x = lambda v, c: 1 + v * 5 + c
    cnf = [[x(v, c) for c in range(5)] for v in range(len(VH))]
    for a, b in e1 + e2:
        cnf += [[-x(a, c), -x(b, c)] for c in range(5)]
    ia, ib = VH.index(A), VH.index(B)
    cnf += [[-x(ia, c), -x(ib, c)] for c in range(5)]
    # symmetry breaking: some triangle of the graph takes colours 0, 1, 2
    nb = {}
    for a, b in e1 + e2:
        nb.setdefault(a, set()).add(b)
        nb.setdefault(b, set()).add(a)
    tri = next((a, b, c) for a, b in e1 + e2 for c in nb[a] & nb[b])
    cnf += [[x(v, k)] for k, v in enumerate(tri)]
    assert Solver(name="cd19", bootstrap_with=cnf).solve() is False   # A = B forced


def test_lambda_arithmetic():
    ca, sa = F.rational(Fr(49, 50)), F.rational(Fr(3, 50)) * R11
    assert ca * ca + sa * sa == F.rational(1)
    # |1 - lambda|^2 = 2 - 2 cos = 1/25: it closes distance 5
    assert F.rational(2) - ca * 2 == F.rational(Fr(1, 25))
    # nu = (-1 + 3 sqrt-11)/10 squares to -lambda
    cn, sn = F.rational(Fr(-1, 10)), F.rational(Fr(3, 10)) * R11
    assert cn * cn - sn * sn == -ca and cn * sn * 2 == -sa


def _module_coordinates(g):
    out = set()
    for a, b in g.edges():
        d = (g.vertices[b].x - g.vertices[a].x, g.vertices[b].y - g.vertices[a].y)
        out.add(tuple(d[0].c) + tuple(d[1].c))
    den = 1
    for v in out:
        for q in v:
            den = den * Fr(q).denominator // gcd(den, Fr(q).denominator)
    vecs = sorted({max(w, tuple(-t for t in w))
                   for w in (tuple(int(Fr(q) * den) for q in v) for v in out)})
    # integer echelon basis of the module (unimodular row operations only)
    basis = []
    for v in vecs:
        v = list(v)
        while any(v):
            p = next(i for i, t in enumerate(v) if t)
            row = next((r for r in basis if next(i for i, t in enumerate(r) if t) == p), None)
            if row is None:
                basis.append(v if v[p] > 0 else [-t for t in v])
                break
            a, b = row[p], v[p]
            x0, x1, y0, y1, aa, bb = 1, 0, 0, 1, a, b
            while bb:
                q = aa // bb
                aa, bb = bb, aa - q * bb
                x0, x1, y0, y1 = x1, x0 - q * x1, y1, y0 - q * y1
            new = [x0 * r + y0 * w for r, w in zip(row, v)]
            v = [(a // aa) * w - (b // aa) * r for r, w in zip(row, v)]
            row[:] = new if new[p] > 0 else [-t for t in new]
    basis.sort(key=lambda r: next(i for i, t in enumerate(r) if t))

    def coords(v):
        v, c = list(v), []
        for r in basis:
            p = next(i for i, t in enumerate(r) if t)
            assert v[p] % r[p] == 0
            k = v[p] // r[p]
            c.append(k)
            v = [a - k * b for a, b in zip(v, r)]
        assert not any(v)
        return c
    return [coords(v) for v in vecs]


def _load(name):
    d = json.load(open(f"{ROOT}/data/{name}"))
    G = Field(tuple(d["field_generators"]))
    return G, [Point(G.element([Fr(a, b) for a, b in x]), G.element([Fr(a, b) for a, b in y]))
               for x, y in d["points"]]


def test_a_denominator_of_five_blocks_and_integrality_does_not():
    G, P = _load("five_247_c.json")
    g = build_graph(P)
    C = _module_coordinates(g)
    assert has_homomorphism(C, 5)[0] is not None           # integral: coset colouring
    v0 = max(range(g.n), key=lambda v: len(g.adj[v]))
    c = P[v0]
    ca, sa = G.rational(Fr(49, 50)), G.rational(Fr(3, 50)) * G.sqrt(11)
    rot = lambda p: Point(c.x + (p.x - c.x) * ca - (p.y - c.y) * sa,
                          c.y + (p.x - c.x) * sa + (p.y - c.y) * ca)
    g2 = build_graph(list(dict.fromkeys(P + [rot(p) for p in P])))
    C2 = _module_coordinates(g2)
    assert any(all(t % 5 == 0 for t in v) for v in C2)     # an edge vector in 5M
    assert has_homomorphism(C2, 5)[0] is None


def test_coset_colourings_of_the_803_span_the_dual():
    G, P = _load("five_247_c.json")
    C = np.array(_module_coordinates(build_graph(P))) % 5
    r = C.shape[1]
    psi = np.array(list(itertools.product(range(5), repeat=r)), dtype=np.int64)
    ok = np.ones(len(psi), bool)
    for c in C:
        ok &= (psi @ c) % 5 != 0
    adm = psi[ok]
    assert len(adm) == 1728
    # rank of the admissible functionals mod 5 is the full rank: common kernel 5M
    M = [[int(t) for t in row] for row in adm[:400]]
    rk = 0
    for col in range(r):
        piv = next((i for i in range(rk, len(M)) if M[i][col] % 5), None)
        if piv is None:
            continue
        M[rk], M[piv] = M[piv], M[rk]
        inv = pow(M[rk][col], 3, 5)
        M[rk] = [(t * inv) % 5 for t in M[rk]]
        for i in range(len(M)):
            if i != rk and M[i][col]:
                f = M[i][col]
                M[i] = [(a - f * b) % 5 for a, b in zip(M[i], M[rk])]
        rk += 1
    assert rk == r == 8


def test_three_distance_lattice_is_six_chromatic_on_61_points():
    norm = lambda x, y: x * x - x * y + y * y
    pts = [(x, y) for x in range(-4, 5) for y in range(-4, 5) if norm(x, y) <= 16]
    idx = {p: i for i, p in enumerate(pts)}
    steps = [(a, b) for a in range(-4, 5) for b in range(-4, 5) if norm(a, b) in (3, 4, 7)]
    E = {(min(i, idx[(x + a, y + b)]), max(i, idx[(x + a, y + b)]))
         for (x, y), i in idx.items() for a, b in steps if (x + a, y + b) in idx}
    assert len(pts) == 61
    for k, want in ((5, False), (6, True)):
        x = lambda v, c: 1 + v * k + c
        cnf = [[x(v, c) for c in range(k)] for v in range(len(pts))]
        for i, j in E:
            cnf += [[-x(i, c), -x(j, c)] for c in range(k)]
        assert Solver(name="cd19", bootstrap_with=cnf).solve() is want
