"""The ring theorems, recomputed rather than quoted.

Six facts decided this pass's direction, and each is cheap enough to re-derive
in a test: the unit square's four clashing sides, sigma's two identities, the
two-ring configuration's chromatic number, the 60-degree grading inside a real
neighbourhood, the absence of sqrt2 from an Eisenstein carrier, and the pair of
degree-4 ceiling points that priced a forced pair.
"""
import json
import math
from fractions import Fraction as Fr

import pytest
from pysat.solvers import Solver

from hn.field import Field
from hn.geometry import Point
from hn.graph import build_graph

ROOT = "/home/user/darwin-50/research/hadwiger-nelson"


def _chi(n, edges, cap=6):
    for k in range(1, cap + 1):
        x = lambda v, c: 1 + v * k + c
        cnf = [[x(v, c) for c in range(k)] for v in range(n)]
        for a, b in edges:
            for c in range(k):
                cnf.append([-x(a, c), -x(b, c)])
        s = Solver(name="cd19", bootstrap_with=cnf)
        ok = s.solve()
        s.delete()
        if ok:
            return k
    return cap + 1


# ---- the unit square: two copies close a width-two hub disjunction ----------

def test_unit_square_has_four_clashing_sides():
    """h the centre, one diagonal the partners, the other their images."""
    f = Field((3,))
    h = Point(f.rational(Fr(1, 2)), f.rational(Fr(1, 2)))
    v1 = Point(f.rational(0), f.rational(0))
    v2 = Point(f.rational(1), f.rational(1))
    rot = lambda p: Point(f.rational(1) - p.y, p.x)      # 90 degrees about h
    one = f.rational(1)
    assert (rot(h) - h).norm2() == f.rational(0)
    for a in (v1, v2):
        for b in (rot(v1), rot(v2)):
            assert (a - b).norm2() == one
    assert (h - v1).norm2() == f.rational(Fr(1, 2))
    assert (h - v2).norm2() == f.rational(Fr(1, 2))
    assert (v1 - v2).norm2() == f.rational(2)


def test_two_is_not_loeschian():
    """No Eisenstein carrier has a pair at sqrt2, so no square centre."""
    assert all(a * a + a * b + b * b != 2
               for a in range(-40, 41) for b in range(-40, 41))


# ---- sigma: the exact cross-ring map ---------------------------------------

@pytest.mark.parametrize("gens", [(3, 11), (3, 11, 247)])
def test_sigma_identities(gens):
    """|sigma u|^2 = |u|^2 / 3  and  |u - sigma u| = |u|, exactly."""
    f = Field(gens)
    r11, sixth = f.sqrt(11), f.rational(Fr(1, 6))
    sigma = lambda p: Point((p.x - r11 * p.y) * sixth, (r11 * p.x + p.y) * sixth)
    origin = Point(f.rational(0), f.rational(0))
    r3 = f.sqrt(3)
    samples = [Point(f.rational(1), f.rational(0)),
               Point(f.rational(Fr(1, 2)), r3 * f.rational(Fr(1, 2))),
               Point(f.rational(2), f.rational(Fr(-3, 5))),
               Point(r3, f.rational(Fr(7, 4)))]
    three = f.rational(3)
    for p in samples:
        n2 = (p - origin).norm2()
        assert (sigma(p) - origin).norm2() * three == n2
        assert (p - sigma(p)).norm2() == n2


def test_sigma_builds_a_cross_ring_edge():
    """Centred anywhere, sigma_c(u) is a unit from u and 1/sqrt3 from c."""
    f = Field((3, 11))
    r11, sixth = f.sqrt(11), f.rational(Fr(1, 6))
    c = Point(f.rational(Fr(2, 3)), f.sqrt(3) * f.rational(Fr(1, 2)))
    ang = f.rational(Fr(1, 2))
    u = Point(c.x + ang, c.y + f.sqrt(3) * ang)          # |u - c| = 1
    assert (u - c).norm2() == f.rational(1)
    dx, dy = u.x - c.x, u.y - c.y
    s = Point(c.x + (dx - r11 * dy) * sixth, c.y + (r11 * dx + dy) * sixth)
    assert (s - u).norm2() == f.rational(1)
    assert (s - c).norm2() == f.rational(Fr(1, 3))


# ---- the two-ring configuration is 3-chromatic, exactly --------------------

@pytest.mark.parametrize("levels", [2, 3, 6, 11, 20])
def test_two_ring_configuration_is_three_chromatic(levels):
    """Hexagons on even levels, two triangles on odd, a matching between."""
    idx = {(i, j): j * 6 + i for j in range(levels) for i in range(6)}
    E = set()
    for j in range(levels):
        step = 1 if j % 2 == 0 else 2
        for i in range(6):
            for s in (step, -step):
                E.add(tuple(sorted((idx[(i, j)], idx[((i + s) % 6, j)]))))
            if j + 1 < levels:
                E.add(tuple(sorted((idx[(i, j)], idx[(i, j + 1)]))))
    assert _chi(6 * levels, sorted(E)) == 3


def test_theta_zero_is_not_a_rational_multiple_of_pi():
    """2 cos(theta0) = 1/sqrt3 is not an algebraic integer: 3x^2 - 1 is not monic."""
    two_cos = 2 * (math.sqrt(3) / 6)
    assert abs(3 * two_cos * two_cos - 1) < 1e-12
    # were theta0 = p/q * pi, 2 cos theta0 would be an algebraic integer, so its
    # minimal polynomial would be monic over Z; 3x^2 - 1 has no monic integer
    # multiple with the same root, since x^2 = 1/3 is not integral.
    assert Fr(1, 3).denominator != 1


# ---- the 60-degree grading, inside a real neighbourhood --------------------

def _load(name):
    d = json.load(open(f"{ROOT}/data/{name}"))
    f = Field(tuple(d["field_generators"]))
    pts = [Point(f.element([Fr(a, b) for a, b in x]),
                 f.element([Fr(a, b) for a, b in y])) for x, y in d["points"]]
    return f, build_graph(pts)


def test_grading_matches_the_bipartition():
    """Inside a component of N(p), colour parity is angle parity."""
    from collections import deque

    _, g = _load("five_247_c.json")
    n = g.n
    adj = [set() for _ in range(n)]
    for x, y in g.edges():
        adj[x].add(y)
        adj[y].add(x)
    hx = [q.fx for q in g.vertices]
    hy = [q.fy for q in g.vertices]
    checked = 0
    for p in range(0, n, 37):
        nb = sorted(adj[p])
        if len(nb) < 6:
            continue
        part, root, cid = {}, {}, 0
        for u in nb:
            if u in part:
                continue
            part[u], root[u] = 0, cid
            q = deque([u])
            while q:
                w = q.popleft()
                for z in adj[w] & set(nb):
                    if z not in part:
                        part[z], root[z] = 1 - part[w], cid
                        q.append(z)
            cid += 1
        ang = {u: math.atan2(hy[u] - hy[p], hx[u] - hx[p]) for u in nb}
        for a in nb:
            for b in nb:
                if a >= b or root[a] != root[b]:
                    continue
                idx = ((ang[b] - ang[a]) % (2 * math.pi)) / (math.pi / 3)
                assert abs(idx - round(idx)) < 1e-6
                assert round(idx) % 2 == (part[a] + part[b]) % 2
                checked += 1
    assert checked > 100


# ---- the ceiling at four, and what a forced pair costs ---------------------

def test_the_two_degree_four_ceiling_points():
    """Exactly two vertices of the 803-graph have degree 4, with one edge inside."""
    from itertools import combinations

    f, g = _load("five_247_c.json")
    adj = [set() for _ in range(g.n)]
    for x, y in g.edges():
        adj[x].add(y)
        adj[y].add(x)
    four = [v for v in range(g.n) if len(adj[v]) == 4]
    assert four == [315, 316]
    one = f.rational(1)
    for p in four:
        nb = sorted(adj[p])
        inside = [(u, v) for u, v in combinations(nb, 2)
                  if (g.vertices[u] - g.vertices[v]).norm2() == one]
        assert len(inside) == 1
        d2 = sorted(float((g.vertices[u] - g.vertices[v]).norm2())
                    for u, v in combinations(nb, 2))
        assert d2[0] == pytest.approx(0.028382, abs=1e-6)
        assert any(abs(x - 1 / 3) < 1e-9 for x in d2)


def test_the_ring_of_a_hub_is_a_union_of_triangles():
    """Two points of the 1/sqrt3 ring are adjacent iff 120 degrees apart."""
    from collections import defaultdict

    f, g = _load("five_247_c.json")
    third, one = f.rational(Fr(1, 3)), f.rational(1)
    hx = [q.fx for q in g.vertices]
    hy = [q.fy for q in g.vertices]
    ring = defaultdict(list)
    for i in range(g.n):
        for j in range(i + 1, g.n):
            dd = (hx[i] - hx[j]) ** 2 + (hy[i] - hy[j]) ** 2
            if abs(dd - 1 / 3) < 1e-9 and (g.vertices[i] - g.vertices[j]).norm2() == third:
                ring[i].append(j)
                ring[j].append(i)
    assert ring, "the graph should carry 1/sqrt3 rings"
    hub = max(ring, key=lambda h: len(ring[h]))
    vs = ring[hub]
    for a in vs:
        for b in vs:
            if a >= b:
                continue
            sep = abs(math.atan2(hy[a] - hy[hub], hx[a] - hx[hub])
                      - math.atan2(hy[b] - hy[hub], hx[b] - hx[hub]))
            sep = min(sep, 2 * math.pi - sep)
            adjacent = (g.vertices[a] - g.vertices[b]).norm2() == one
            assert adjacent == (abs(sep - 2 * math.pi / 3) < 1e-7)


# ---- the two-ring theorem, upgraded from a computation to a proof ----------

def test_two_ring_period_two_colouring_is_proper_forever():
    """An explicit 3-colouring of period 2 settles the infinite configuration.

    Computing chi = 3 to thirty levels leaves the infinite graph open.  The
    colouring the solver returned has period two in j, so checking it on one
    period checks it everywhere: every edge of the configuration joins two
    vertices whose levels differ by at most one, and the pattern repeats.
    """
    even = [0, 1, 2, 0, 2, 1]        # hexagon level, i ~ i +- 1 (mod 6)
    odd = [2, 2, 1, 1, 0, 0]         # two triangles, i ~ i +- 2 (mod 6)
    for i in range(6):
        assert even[i] != even[(i + 1) % 6]          # the hexagon
        assert odd[i] != odd[(i + 2) % 6]            # the triangles
        assert even[i] != odd[i]                     # the matching between
    assert len(set(even)) == 3 and len(set(odd)) == 3
    # and three is forced from below: an odd level carries a triangle
    assert len({odd[0], odd[2], odd[4]}) == 3
