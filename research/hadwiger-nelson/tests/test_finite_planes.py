"""Bounds for Moorhouse's table of chi(F_q^2), q prime (notes/local_colourings.md, section 12).

The graph: vertices F_q^2, (x, y) ~ (x', y') iff (x - x')^2 + (y - y')^2 = 1.

**Upper bounds: interval colourings.** The line a x + b y = r meets the unit circle exactly when a^2 + b^2 - r^2
is a square (the discriminant of the intersection is b^2 (a^2 + b^2 - r^2)).  If a^2 + b^2 - r^2 is a non-square
for r = 0, 1, ..., m - 1, then a x + b y takes none of the values 0, +-1, ..., +-(m - 1) on a unit vector, and
(x, y) -> floor(((a x + b y) mod q) / m) is a proper colouring with ceil(q / m) colours.  For m = 2 this is
Vinh's colouring by pairs of parallel lines.

**Lower bounds.** No proper 4-colouring (SAT, with a triangle or an edge pinned); the slower cases are in
tests/test_finite_planes_slow.py.  Six and seven colours from Hoffman's ratio bound, with every eigenvalue in
interval arithmetic (section 14), here; six colours from the three-point bound in
tests/test_threepoint_certificates.py.
"""
import sys, os
sys.path.insert(0, os.path.dirname(__file__))

from test_biquadratic_bounds import colourable

# q: (a, b, m); the colouring uses ceil(q / m) colours
INTERVAL = {7: (2, 3, 2), 13: (3, 6, 3), 17: (3, 6, 3), 19: (4, 5, 4), 23: (3, 5, 3), 29: (5, 9, 5),
            31: (4, 13, 4), 37: (5, 6, 5), 41: (6, 19, 6), 43: (5, 8, 5), 47: (5, 22, 5), 53: (5, 13, 5),
            59: (6, 16, 6), 61: (6, 27, 6)}


def plane(q):
    U = [(a, b) for a in range(q) for b in range(q) if (a * a + b * b) % q == 1]
    idx = lambda x, y: (x % q) * q + y % q
    E = sorted({tuple(sorted((idx(x, y), idx(x + a, y + b)))) for x in range(q) for y in range(q) for a, b in U})
    return U, E


def test_interval_colourings():
    for q, (a, b, m) in INTERVAL.items():
        squares = {t * t % q for t in range(q)}
        assert all((a * a + b * b - r * r) % q not in squares for r in range(m)), q
        U, E = plane(q)
        assert not {(a * x + b * y) % q for x, y in U} & ({r % q for r in range(1 - m, m)}), q
        c = [((a * x + b * y) % q) // m for x in range(q) for y in range(q)]
        assert len(set(c)) == -(-q // m), q
        assert all(c[u] != c[v] for u, v in E), q


def test_a_linear_8_colouring_of_f43():
    """(x, y) -> c(x + y mod 43) beats the interval colouring (9 colours): c colours the circulant on
    {x + y : x^2 + y^2 = 1}, which avoids 0"""
    q = 43
    U, E = plane(q)
    S = {(x + y) % q for x, y in U}
    assert 0 not in S
    c = [0, 7, 6, 5, 4, 7, 6, 5, 4, 3, 7, 6, 5, 4, 7, 6, 5, 4, 3, 2, 1, 0, 3, 7, 6, 5, 4, 3, 2, 1, 0, 3, 2, 1, 0,
         7, 6, 5, 4, 0, 3, 2, 1]
    assert all(c[v] != c[(v + s) % q] for v in range(q) for s in S) and len(set(c)) == 8
    col = [c[(x + y) % q] for x in range(q) for y in range(q)]
    assert all(col[u] != col[v] for u, v in E)


def test_the_interval_colourings_are_optimal_for_7_13_19():
    """chi(F_q^2) = 4, 5, 5 for q = 7, 13, 19 (Moorhouse; tests/test_split_places.py; tests/test_biquadratic_bounds.py)"""
    assert {q: -(-q // INTERVAL[q][2]) for q in (7, 13, 19)} == {7: 4, 13: 5, 19: 5}


def pinned_triangle(q, E):
    adj = {}
    for u, v in E:
        adj.setdefault(u, set()).add(v)
        adj.setdefault(v, set()).add(u)
    b = next(v for v in sorted(adj[0]) if adj[0] & adj[v])
    return [0, b, min(adj[0] & adj[b])]


def test_the_planes_over_23_and_37_need_five_colours():
    for q in (23, 37):
        U, E = plane(q)
        assert not colourable(q * q, E, 4, pin=pinned_triangle(q, E)), q


def test_the_spectral_bound_gives_six_and_seven_colours():
    """Hoffman's ratio bound with every eigenvalue in interval arithmetic (scripts/finite_hoffman.py; section 14):
    alpha(F_q^2) < q^2/5 for q = 59, and < q^2/6 for q = 71, 97, 101."""
    sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'scripts'))
    from finite_hoffman import hoffman
    for q, k in [(59, 6), (71, 7), (97, 7), (101, 7)]:
        d, lmin, amax = hoffman(q)
        assert (k - 1) * amax < q * q, q          # so chi >= q^2 / alpha > k - 1


def test_proposition_b_where_weil_does_not_reach():
    """chi(G_q) >= 6 for the anisotropic planes G_q (x^2 - n y^2, n a non-square) with q = 53, 59, 61: Weil's
    bound gives it only for q > 62 (notes/local_colourings.md, sections 4 and 14)"""
    sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'scripts'))
    from finite_hoffman import hoffman
    for q in (53, 59, 61):
        d, lmin, amax = hoffman(q, 'inert')
        assert d == q + 1 and 5 * amax < q * q, q
