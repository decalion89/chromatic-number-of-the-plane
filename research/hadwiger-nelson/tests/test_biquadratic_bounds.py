"""Local upper bounds for the planes over two square roots (notes/local_colourings.md, section 13).

A place v of L that does not split in K = L(i) bounds chi(L^2) (Proposition A): by 2 if K_w / L_v is
ramified above 2 (the norm-one residues collapse to +1 = -1), by 3 if ramified elsewhere, and by the
chromatic number of the finite plane G_q = Cay(F_{q^2}, N_1) if unramified with residue field F_q.
For an odd prime q = 3 (mod 4), G_q is Cay(F_q^2, {a^2 + b^2 = 1}).
"""
import sys, os, itertools
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "scripts"))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "scripts", "experiments"))

from pysat.solvers import Solver

from fieldscreen import screen
from classify_biquadratic import local_bound
from test_split_places import gf


def circle_plane(q):
    """G_q for an odd prime q = 3 mod 4: vertices F_q^2, steps (a, b) with a^2 + b^2 = 1"""
    U = [(a, b) for a in range(q) for b in range(q) if (a * a + b * b) % q == 1]
    V = [(a, b) for a in range(q) for b in range(q)]
    idx = {v: i for i, v in enumerate(V)}
    E = set()
    for a, b in V:
        for u, w in U:
            i, j = idx[(a, b)], idx[((a + u) % q, (b + w) % q)]
            E.add((min(i, j), max(i, j)))
    return V, idx, sorted(E), U


def colourable(n, E, k, pin=()):
    X = lambda v, c: 1 + v * k + c
    s = Solver(name="cd19")
    for v in range(n):
        s.add_clause([X(v, c) for c in range(k)])
    for a, b in E:
        for c in range(k):
            s.add_clause([-X(a, c), -X(b, c)])
    for c, v in enumerate(pin):
        s.add_clause([X(v, c)])
    r = s.solve()
    s.delete()
    return r


def test_small_finite_planes():
    for q, chi in [(3, 3), (7, 4), (11, 5)]:
        V, _, E, _ = circle_plane(q)
        assert not colourable(len(V), E, chi - 1) and colourable(len(V), E, chi), q


def test_g19_is_5_chromatic():
    """A linear 5-colouring, (a, b) -> col[a + b mod 19], since {a + b : a^2 + b^2 = 1} avoids 0 and col
    colours that circulant; no 4-colouring (G_19 has no triangle, so an edge is pinned)."""
    V, idx, E, U = circle_plane(19)
    S = sorted({(a + b) % 19 for a, b in U})
    assert S == [1, 2, 4, 6, 9, 10, 13, 15, 17, 18]
    col = [4, 3, 2, 4, 3, 2, 1, 0, 4, 3, 0, 4, 3, 2, 1, 0, 2, 1, 0]
    c = [col[(a + b) % 19] for a, b in V]
    assert all(c[u] != c[v] for u, v in E)
    assert not colourable(len(V), E, 4, pin=E[0])


def test_the_planes_over_f2_and_f4():
    """Above 2, unramified: G_2 = Cay(F_4, F_4^*) = K_4, and G_4 = Cay(F_16, mu_5) is 4-chromatic."""
    els, add, mul, zero, one = gf(16)
    def pw(a, e):
        r = one
        for _ in range(e):
            r = mul(r, a)
        return r
    N1 = [a for a in els if a != zero and pw(a, 5) == one]
    assert len(N1) == 5
    idx = {e: i for i, e in enumerate(els)}
    E = sorted({(min(idx[a], idx[add(a, u)]), max(idx[a], idx[add(a, u)])) for a in els for u in N1})
    assert not colourable(16, E, 3) and colourable(16, E, 4)


def test_ramified_above_2_gives_a_bipartite_plane():
    """Q(sqrt2, sqrt5): the place above 2 is ramified in K / L, so chi <= 2 (Moorhouse's Lemma 8.4 is the
    quadratic case d = 1 mod 4)."""
    places = screen([1, 2, 5], 10)
    assert places[0][:2] == (2, "ramified")
    assert local_bound(2, 5)[0] == 2
    assert local_bound(5, 13)[0] == 2


def test_the_bounds_with_sqrt3_match_theorem_5():
    """Theorem 5 of notes/local_colourings.md, section 11. Q(sqrt3, sqrt q), q < 60 prime: bound 3 iff q = 1 (mod 3), bound 4 iff q = 2 or q = 11, 17 (mod 24),
    and nothing below 5 for q = 5, 23 (mod 24)."""
    for q in [2, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59]:
        bd, _ = local_bound(3, q)
        if q % 3 == 1:
            assert bd == 3, q
        elif q == 2 or q % 24 in (11, 17):
            assert bd == 4, q
        else:
            assert bd is None or bd >= 5, q


def test_q3_29_has_no_small_non_split_place():
    """Q(sqrt3, sqrt29): the first place that does not split in K lies above 23."""
    places = screen([1, 3, 29], 100)
    assert places[0][0] == 23 and places[0][1] == "unramified"
    assert local_bound(3, 29)[0] is None


def test_a_field_without_sqrt3_that_only_19_bounds():
    """Q(sqrt5, sqrt7): no triangles, and the first non-split place is above 19, where chi(G_19) = 5."""
    places = screen([1, 5, 7], 200)
    assert places[0][:3] == (19, "unramified", 19)
    assert local_bound(5, 7) == (5, (19, "unramified", 19))


def test_a_unit_five_cycle_over_q5_7():
    """1 + conj(s) - t + s - conj(t) = 0 for t = (2 + i sqrt5)/3 and s = (1 + i sqrt35)/6, since
    1 + 2 Re(s) = 2 Re(t) = 4/3: a closed walk of five unit steps, an odd cycle, so 3 <= chi(Q(sqrt5, sqrt7)^2) <= 5."""
    from fractions import Fraction as Fr
    from hn.field import Field
    from hn.geometry import Point
    F = Field((5, 7)); r = F.rational
    s5, s35 = F.sqrt(5), F.sqrt(5) * F.sqrt(7)
    t = Point(r(Fr(2, 3)), s5 * r(Fr(1, 3)))
    s = Point(r(Fr(1, 6)), s35 * r(Fr(1, 6)))
    steps = [Point(r(1), r(0)), Point(s.x, -s.y), Point(-t.x, -t.y), s, Point(-t.x, t.y)]
    for u in steps:
        assert u.x * u.x + u.y * u.y == F.one()
    p = Point(r(Fr(-4, 3)), r(0)); pts = []
    for u in steps:
        pts.append(p); p = Point(p.x + u.x, p.y + u.y)
    assert p == pts[0] and len(set(pts)) == 5
    for i in range(5):
        for j in range(i + 1, 5):
            dx, dy = pts[i].x - pts[j].x, pts[i].y - pts[j].y
            assert (dx * dx + dy * dy == F.one()) == (j - i in (1, 4))    # exactly the cycle's edges
