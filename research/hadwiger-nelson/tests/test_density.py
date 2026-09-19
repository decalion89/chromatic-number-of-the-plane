"""The independence-ratio certificate, and the Moser spindle's 2/7."""
from hn.density import (has_independent_set, independence_number,
                        independence_ratio, measurable_bound)
from hn.geometry import SPINDLE, eisenstein, origin
from hn.graph import build_graph


def moser():
    rh = [origin(), eisenstein(1, 0), eisenstein(0, 1), eisenstein(1, 1)]
    return build_graph(rh + [SPINDLE(p) for p in rh])


def test_moser_spindle_gives_two_sevenths():
    """The classical value, and the one every m_1 bound is measured against.

    Averaging a measurable 1-avoiding set of density d over rigid motions of a
    finite unit-distance graph G gives n d <= alpha(G), so m_1 <= alpha/n. The
    spindle gives 2/7 = 0.2857.
    """
    g = moser()
    a, S = independence_number(g)
    assert a == 2
    assert abs(independence_ratio(g) - 2 / 7) < 1e-12
    assert abs(measurable_bound(g) - 3.5) < 1e-12
    for u, v in g.edges():                     # it really is independent
        assert not (u in S and v in S)


def test_independent_set_search_is_exact_at_the_boundary():
    g = moser()
    assert has_independent_set(g, 2) is not None
    assert has_independent_set(g, 3) is None


def test_a_triangle_free_graph_has_the_ratio_it_should():
    """A single edge: alpha = 1 of 2, so the bound it gives is 2 colours."""
    from hn.geometry import origin as o

    g = build_graph([o(), eisenstein(1, 0)])
    assert independence_number(g)[0] == 1
    assert abs(measurable_bound(g) - 2.0) < 1e-12


def test_the_density_route_is_capped_below_six():
    """Croft's 1967 set has density 0.2293, so no graph beats it.

    Averaging runs both ways: a measurable 1-avoiding set of density d forces
    alpha(G)/n >= d for every finite unit-distance graph. Since 0.2293 > 1/5,
    no ratio ever clears 0.2 and chi_m >= 6 is out of reach by this argument --
    the ceiling is a construction, not a shortage of computation.
    """
    from hn.density import CROFT_DENSITY

    assert CROFT_DENSITY > 1 / 5
    assert 1 / CROFT_DENSITY < 5                  # at most chi_m >= 5, never 6
    assert independence_ratio(moser()) > CROFT_DENSITY


def test_counting_can_never_prove_six():
    """What Croft's bound says about the object being searched for.

    alpha(G)/n < 1/5 would give chi(G) >= 6 outright, with no measure theory:
    a 5-colouring splits n vertices into five independent sets and one has at
    least n/5. But averaging forces alpha(G)/n >= 0.2293 for every finite
    unit-distance graph, so no graph ever qualifies.

    Turned around, that is a fact about the target. A 6-chromatic
    unit-distance graph must have independent sets of at least 0.2293 n, so
    five of them cover at least 1.1465 n -- more than the whole graph.
    Whatever stops it from 5-colouring, it is never counting.
    """
    from hn.density import CROFT_DENSITY

    assert CROFT_DENSITY > 0.2
    assert 5 * CROFT_DENSITY > 1.0            # five classes have room to spare
