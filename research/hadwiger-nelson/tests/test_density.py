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
