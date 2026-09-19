"""Forced colour relations, and the ladder chi >= 6 reduces to."""
import pytest

from hn.coloring import is_k_colorable
from hn.degrey import build_G
from hn.forced import (ColourRelations, centre_pressure, find_rainbow,
                       min_colours_on)
from hn.graph import build_graph
from hn.mixed import joint_core_union


def test_edges_are_forced_different():
    g = build_graph(joint_core_union())
    rel = ColourRelations(g, 4)
    try:
        assert rel.colourable
        for u, v in list(g.edges())[:20]:
            assert rel.different(u, v)
            assert not rel.can_share(u, v)
    finally:
        rel.close()


def test_min_colours_recovers_the_chromatic_number():
    """On the whole vertex set the measurement is just chi."""
    g = build_graph(joint_core_union())
    rel = ColourRelations(g, 4)
    try:
        assert min_colours_on(rel, range(g.n)) == 4
        assert min_colours_on(rel, []) == 0
        assert min_colours_on(rel, [0]) == 1
    finally:
        rel.close()


def test_vacuous_queries_are_refused():
    """A graph with no k-colouring makes every query UNSAT; say so."""
    g = build_graph(joint_core_union())
    rel = ColourRelations(g, 3)
    try:
        assert not rel.colourable
    finally:
        rel.close()
    with pytest.raises(ValueError):
        find_rainbow(g, 3, 3)


def test_unit_circle_graphs_are_bipartite():
    """Two points of a unit circle are adjacent only at 60 degrees, and the
    60-degree orbit closes after six steps, so every component of a circle
    graph is a path or a 6-cycle and no circle ever needs three colours on its
    own account. Any third colour has to be forced by the ambient graph --
    which is why every neighbourhood of de Grey's G comes back at exactly two,
    however large: the two hubs of degree 60 included.
    """
    for g in (build_graph(joint_core_union()), build_G()):
        for v in range(g.n):
            circle = set(g.adj[v])
            colour = {}
            for a in sorted(circle):
                if a in colour:
                    continue
                colour[a], stack = 0, [a]
                while stack:
                    x = stack.pop()
                    for y in g.adj[x] & circle:
                        if y not in colour:
                            colour[y] = 1 - colour[x]
                            stack.append(y)
                        else:
                            assert colour[y] != colour[x], \
                                f"circle of {v} is not bipartite at {x},{y}"


@pytest.mark.slow
def test_degrey_G_has_no_local_pressure_at_five():
    """Every one of the 1581 neighbourhoods is 2-colourable at k=5.

    The ladder needs three: at three the pivot's colour is confined to two,
    which is the joint core the two-orbit block closes. This pins the measured
    fact that G supplies none of it.
    """
    g = build_G()
    rel = ColourRelations(g, 5)
    try:
        assert rel.colourable
        seen = {min_colours_on(rel, sorted(g.adj[v])) for v in range(g.n)}
        assert seen == {2}
    finally:
        rel.close()


def test_centre_pressure_reads_the_ladder():
    g = build_graph(joint_core_union())
    rel = ColourRelations(g, 4)
    try:
        hub = max(range(g.n), key=lambda v: len(g.adj[v]))
        r = centre_pressure(rel, sorted(g.adj[hub]))
        assert r["pressure"] + r["free_for_centre"] == 4
        assert r["smallest_possible_core"] == 4 - r["pressure"]
    finally:
        rel.close()


# -- the pressure theorem -------------------------------------------------

def test_pressure_bound_matches_the_construction():
    """A core of r needs pressure >= k - r; the built core of two saturates it.

    On the four points at three colours the pivot sees one neighbour, so
    pressure is 1 = k - 2, the theorem permits a core of two, and a separation
    query confirms one. The bound is tight exactly where the construction is.
    """
    from hn.forced import max_core_bound, pressure
    from hn.mixed import joint_core_configuration
    from hn.spindle import SeparationTest

    _E, pts, _rho, _sigma = joint_core_configuration()
    g = build_graph(pts)
    bp = g.vertices.index(pts[0])
    ia, ib = g.vertices.index(pts[2]), g.vertices.index(pts[3])
    rel = ColourRelations(g, 3)
    try:
        assert pressure(rel, bp) == 1
        assert max_core_bound(rel, bp) == 2
    finally:
        rel.close()
    st = SeparationTest(g, 3, bp, [ia, ib])
    try:
        assert not st.run(subset=[ia, ib])[0]
    finally:
        st.close()


def test_pressure_forbids_small_cores_where_it_is_low():
    """The contrapositive, checked directly: squeeze the neighbourhood into
    k - r - 1 colours and every set of r vertices fails to be a core."""
    from itertools import combinations

    from hn.forced import pressure
    from hn.spindle import SeparationTest

    g = build_graph(joint_core_union())
    rel = ColourRelations(g, 4)
    try:
        hub = max(range(g.n), key=lambda v: len(g.adj[v]))
        assert pressure(rel, hub) == 1          # so no core of size 1 or 2
    finally:
        rel.close()
    others = [v for v in range(g.n) if v != hub and v not in g.adj[hub]]
    st = SeparationTest(g, 4, hub, others)
    try:
        for pair in list(combinations(others, 2))[:40]:
            assert st.run(subset=list(pair))[0], f"{pair} should not be a core"
    finally:
        st.close()
