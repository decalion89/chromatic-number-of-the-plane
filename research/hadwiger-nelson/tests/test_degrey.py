"""The de Grey reconstruction, checked at every stage.

Vertex counts are a strong transcription check: 39 base points pushed through
a chain of rotations land on 397 and then 1581 only if every coordinate is
right.  A single wrong digit changes which points coincide and the counts move.
"""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

import pytest

from hn.degrey import S_POINTS, build_G, build_S, build_Sa, build_Sb, build_Y
from hn.geometry import DEGREY_FIELD


def test_S_has_the_published_size():
    assert len(S_POINTS) == 39
    assert len(build_S()) == 39


def test_S_lies_in_the_small_field():
    """S itself needs only sqrt3 and sqrt11; the later rotations are what
    force the extension to sqrt7 and sqrt15."""
    F = DEGREY_FIELD
    allowed = {0, F._prod.index(3), F._prod.index(11), F._prod.index(33)}
    for p in build_S():
        for coord in (p.x, p.y):
            for m, c in enumerate(coord.c):
                assert c == 0 or m in allowed


def test_symmetrisation_gives_397_points():
    assert len(build_Sa()) == 397


def test_Y_drops_exactly_two_vertices():
    sa, sb = build_Sa(), build_Sb()
    merged = len(set(sa) | set(sb))
    assert len(build_Y()) == merged - 2


def test_G_has_1581_vertices():
    assert len(build_G(as_graph=False)) == 1581


def test_G_is_a_unit_distance_graph_with_the_published_edge_count():
    g = build_G()
    assert g.n == 1581
    assert g.m == 7877
    assert min(g.degrees()) >= 4


@pytest.mark.slow
def test_G_is_not_four_colourable():
    """chi(R^2) >= 5, reproduced from the published recipe."""
    from hn.coloring import is_k_colorable

    assert is_k_colorable(build_G(), 4, timeout=1800)[0] is False
